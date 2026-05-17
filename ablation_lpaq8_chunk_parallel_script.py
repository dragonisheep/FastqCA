"""Ablation: chunk raw FASTQ and dispatch chunks directly to lpaq8 (no mapping/CA).

Workflow:
1) Split original FASTQ into blocks by reads_per_block (derived from BLOCK_SIZE).
2) Dispatch each chunk to worker pool.
3) Each worker directly runs lpaq8 compression on the chunk.
4) Aggregate compressed sizes and metrics into CSV.
"""

import os
import csv
import time
import tempfile
import subprocess
import multiprocessing
from pathlib import Path
from datetime import datetime

import psutil
from Bio import SeqIO

# ===== 配置区（可按需修改） =====
INPUT_DIR = '/media/compress/新加卷1/New_Test'
OUTPUT_DIR = '/media/compress/新加卷/output/New_Test_LossLess/Ablation_lpaq8_chunk_parallel'
BLOCK_SIZE = 128 * 1024 * 1024  # 与原实验保持一致
WORKERS = 4                      # 分流并行 worker 数
COMPRESSION_LEVEL = '9'
LPAQ8_PATH = str((Path(__file__).resolve().parent / 'lpaq8').resolve())
# ===========================


def get_file_size(file_path):
    return os.path.getsize(file_path)


def calc_reads_per_block(fastq_path, block_size):
    with open(fastq_path, 'r') as fh:
        try:
            first = next(SeqIO.parse(fh, 'fastq'))
        except StopIteration:
            return 1
    read_len = len(first.seq)
    bytes_per_read = max(read_len * 2, 1)
    return max(block_size // bytes_per_read, 1)


def split_fastq_to_chunks(input_fastq, chunk_dir, reads_per_block):
    chunk_paths = []
    chunk_idx = 0
    records = []

    with open(input_fastq, 'r') as handle:
        for rec in SeqIO.parse(handle, 'fastq'):
            records.append(rec)
            if len(records) >= reads_per_block:
                chunk_path = os.path.join(chunk_dir, f'chunk_{chunk_idx}.fastq')
                with open(chunk_path, 'w') as out:
                    SeqIO.write(records, out, 'fastq')
                chunk_paths.append(chunk_path)
                chunk_idx += 1
                records = []

    if records:
        chunk_path = os.path.join(chunk_dir, f'chunk_{chunk_idx}.fastq')
        with open(chunk_path, 'w') as out:
            SeqIO.write(records, out, 'fastq')
        chunk_paths.append(chunk_path)

    return chunk_paths


def compress_chunk_worker(chunk_path):
    out_path = f"{chunk_path}.lpaq8"
    cmd = [LPAQ8_PATH, COMPRESSION_LEVEL, chunk_path, out_path]
    p = subprocess.Popen(cmd)
    proc = psutil.Process(p.pid)
    peak_mb = 0.0
    while p.poll() is None:
        try:
            peak_mb = max(peak_mb, proc.memory_info().rss / 1024 / 1024)
        except Exception:
            pass
        time.sleep(0.05)
    if p.returncode != 0:
        raise RuntimeError(f"lpaq8 failed: {' '.join(cmd)}")
    return out_path, peak_mb


def run_one_file(input_file):
    start = time.time()
    input_size = get_file_size(input_file)
    file_name = os.path.splitext(os.path.basename(input_file))[0]

    with tempfile.TemporaryDirectory(prefix=f"abl_lpaq8_parallel_{file_name}_") as temp_dir:
        reads_per_block = calc_reads_per_block(input_file, BLOCK_SIZE)
        chunks = split_fastq_to_chunks(input_file, temp_dir, reads_per_block)

        with multiprocessing.Pool(processes=WORKERS) as pool:
            results = pool.map(compress_chunk_worker, chunks)

        total_out_size = 0
        peak_mem = 0.0
        for out_path, pmb in results:
            total_out_size += get_file_size(out_path)
            peak_mem = max(peak_mem, pmb)

        with open(os.path.join(OUTPUT_DIR, f"{file_name}.chunked_parallel_lpaq8.size.txt"), 'w') as f:
            f.write(str(total_out_size))

    elapsed = time.time() - start
    ratio = input_size / total_out_size if total_out_size > 0 else 0
    speed = (input_size / 1024 / 1024) / elapsed if elapsed > 0 else 0

    return {
        'file_name': file_name,
        'input_size_mb': input_size / 1024 / 1024,
        'output_size_mb': total_out_size / 1024 / 1024,
        'compression_ratio': ratio,
        'compression_time_s': elapsed,
        'compression_speed_mbs': speed,
        'max_memory_mb': peak_mem,
        'reads_per_block': reads_per_block,
        'block_size_bytes': BLOCK_SIZE,
        'workers': WORKERS,
        'compression_level': COMPRESSION_LEVEL,
    }


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    csv_path = os.path.join(OUTPUT_DIR, f"ablation_lpaq8_chunk_parallel_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv")
    header = [
        'file_name', 'input_size_mb', 'output_size_mb', 'compression_ratio',
        'compression_time_s', 'compression_speed_mbs', 'max_memory_mb',
        'reads_per_block', 'block_size_bytes', 'workers', 'compression_level'
    ]

    with open(csv_path, 'w', newline='') as f:
        csv.writer(f).writerow(header)

    for file_path in Path(INPUT_DIR).glob('*.fastq'):
        print(f'Processing {file_path} ...')
        try:
            metrics = run_one_file(str(file_path))
            with open(csv_path, 'a', newline='') as f:
                csv.writer(f).writerow([metrics[k] for k in header])
            print(f"Done: {file_path.name}, ratio={metrics['compression_ratio']:.3f}")
        except Exception as e:
            print(f'Error: {file_path} -> {e}')


if __name__ == '__main__':
    main()
