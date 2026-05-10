"""Ablation experiment 1: chunk FASTQ files then compress each chunk directly with lpaq8.

Chunking strategy follows main_new.py defaults:
- block_size = 128 * 1024 * 1024 bytes
- reads_per_block = block_size // (read_length * 2)
"""

import os
import csv
import time
import subprocess
import tempfile
from pathlib import Path
from datetime import datetime

import psutil
from Bio import SeqIO

# ===== 配置区（可按需修改） =====
INPUT_DIR = '/media/compress/新加卷1/New_Test'
OUTPUT_DIR = '/media/compress/新加卷/output/New_Test_LossLess/Ablation_lpaq8_chunk'
BLOCK_SIZE = 128 * 1024 * 1024  # 与原实验保持一致
LPAQ8_PATH = str((Path(__file__).resolve().parent / 'lpaq8').resolve())
COMPRESSION_LEVEL = '9'
# ===========================


def get_file_size(file_path):
    return os.path.getsize(file_path)


def monitor_process(process):
    cpu_percentages = []
    max_total_memory = 0

    try:
        parent = psutil.Process(process.pid)
        while process.poll() is None:
            try:
                procs = [parent] + parent.children(recursive=True)
                total_rss = 0.0
                total_cpu = 0.0
                for proc in procs:
                    try:
                        total_rss += proc.memory_info().rss / 1024 / 1024
                        total_cpu += proc.cpu_percent()
                    except Exception:
                        continue
                cpu_percentages.append(total_cpu)
                max_total_memory = max(max_total_memory, total_rss)
            except Exception:
                break
            time.sleep(0.05)
    except Exception:
        pass

    return {
        'avg_cpu': sum(cpu_percentages) / len(cpu_percentages) if cpu_percentages else 0,
        'max_memory': max_total_memory,
    }


def calc_reads_per_block(fastq_path, block_size):
    with open(fastq_path, 'r') as fh:
        try:
            first = next(SeqIO.parse(fh, 'fastq'))
        except StopIteration:
            return 1
    read_len = len(first.seq)
    bytes_per_read = max(read_len * 2, 1)
    rpb = block_size // bytes_per_read
    return max(rpb, 1)


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


def compress_chunk_with_lpaq8(chunk_path):
    out_path = f"{chunk_path}.lpaq8"
    cmd = [LPAQ8_PATH, COMPRESSION_LEVEL, chunk_path, out_path]
    process = subprocess.Popen(cmd)
    metrics = monitor_process(process)
    process.wait()
    if process.returncode != 0:
        raise RuntimeError(f"lpaq8 failed: {' '.join(cmd)}")
    return out_path, metrics


def run_one_file(input_file, output_dir):
    start = time.time()
    input_size = get_file_size(input_file)
    file_name = os.path.splitext(os.path.basename(input_file))[0]

    with tempfile.TemporaryDirectory(prefix=f"abl_lpaq8_{file_name}_") as temp_dir:
        reads_per_block = calc_reads_per_block(input_file, BLOCK_SIZE)
        chunks = split_fastq_to_chunks(input_file, temp_dir, reads_per_block)

        total_out_size = 0
        cpu_values = []
        peak_mem = 0

        for c in chunks:
            out_chunk, m = compress_chunk_with_lpaq8(c)
            total_out_size += get_file_size(out_chunk)
            cpu_values.append(m['avg_cpu'])
            peak_mem = max(peak_mem, m['max_memory'])

        # 复制一份总压缩产物大小对应的占位文件，便于审计
        summary_path = os.path.join(output_dir, f"{file_name}.chunked_lpaq8.size.txt")
        with open(summary_path, 'w') as f:
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
        'avg_cpu_percent': (sum(cpu_values) / len(cpu_values)) if cpu_values else 0,
        'max_memory_mb': peak_mem,
        'reads_per_block': reads_per_block,
        'block_size_bytes': BLOCK_SIZE,
    }


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    csv_path = os.path.join(OUTPUT_DIR, f"ablation_lpaq8_chunk_metrics_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv")
    header = [
        'file_name', 'input_size_mb', 'output_size_mb', 'compression_ratio',
        'compression_time_s', 'compression_speed_mbs', 'avg_cpu_percent',
        'max_memory_mb', 'reads_per_block', 'block_size_bytes'
    ]

    with open(csv_path, 'w', newline='') as f:
        csv.writer(f).writerow(header)

    for file_path in Path(INPUT_DIR).glob('*.fastq'):
        print(f'Processing {file_path} ...')
        try:
            metrics = run_one_file(str(file_path), OUTPUT_DIR)
            with open(csv_path, 'a', newline='') as f:
                w = csv.writer(f)
                w.writerow([metrics[k] for k in header])
            print(f"Done: {file_path.name}, ratio={metrics['compression_ratio']:.3f}")
        except Exception as e:
            print(f'Error: {file_path} -> {e}')


if __name__ == '__main__':
    main()
