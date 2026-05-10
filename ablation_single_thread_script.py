"""Ablation: single-thread compression with cross-block rule learning.

This script intentionally differs from main_new.py behavior:
- It processes blocks sequentially in one process/thread.
- A single rules table is reused across ALL blocks of one FASTQ file.

This enables cross-block learning ablation.
"""

import os
import csv
import time
from pathlib import Path
from datetime import datetime

import psutil
from Bio import SeqIO

import LossLess_thread as LL
import Lossy_thread as LY

# ===== 配置区（可按需修改） =====
COMPRESSOR = 'LossLess'  # 'LossLess' or 'Lossy'
INPUT_DIR = '/media/compress/新加卷1/New_Test'
OUTPUT_DIR = '/media/compress/新加卷/output/New_Test_LossLess/Ablation_cross_block_learning'
BLOCK_SIZE = 128 * 1024 * 1024  # 与原实验保持一致
SAVE_INTERMEDIATE = False
# ===========================


def get_file_size(file_path):
    return os.path.getsize(file_path)


def monitor_memory_start():
    return psutil.Process(os.getpid())


def monitor_memory_peak_mb(proc, current_peak):
    try:
        rss = proc.memory_info().rss / 1024 / 1024
        return max(current_peak, rss)
    except Exception:
        return current_peak


def merge_parts(output_path, total_blocks):
    missing_parts = []
    with open(output_path, "wb") as final_file:
        for i in range(1, total_blocks + 1):
            part_path = os.path.join(os.path.dirname(output_path), f"chunk_{i}.part")
            if os.path.exists(part_path):
                with open(part_path, "rb") as part_file:
                    final_file.write(part_file.read())
                os.remove(part_path)
            else:
                missing_parts.append(part_path)
        final_file.write(b"%eof")
    if missing_parts:
        raise FileNotFoundError(f"以下分块缺失，导致压缩结果不完整: {missing_parts}")


def compress_lossless_cross_block(fastq_path, output_dir, lpaq8_path, save, block_size):
    output_path = LL.get_output_path(fastq_path, output_dir)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    reads_per_block, total_reads = LL.get_reads_num_per_block(fastq_path, block_size)
    if reads_per_block <= 0:
        reads_per_block = 1

    rules_dict = LL.init_rules_dict()  # 关键：全文件仅一份规则表
    block_count = 0

    with open(fastq_path, 'r') as f:
        parser = SeqIO.parse(f, 'fastq')
        while True:
            records = []
            for _ in range(reads_per_block):
                try:
                    records.append(next(parser))
                except StopIteration:
                    break
            if not records:
                break

            block_count += 1
            read_length = len(records[0].seq)
            g_block, quality_block, id_block = LL.process_records(iter(records), len(records), read_length, rules_dict)
            if g_block is None:
                continue
            if save:
                LL.save_intermediate_files(g_block, quality_block, id_block, output_path, block_count)
            LL.back_compress_worker(g_block, quality_block, id_block, lpaq8_path, output_path, save, block_count)

    merge_parts(output_path, block_count)
    return output_path, total_reads, reads_per_block


def compress_lossy_cross_block(fastq_path, output_dir, lpaq8_path, save, block_size):
    output_path = LY.get_output_path(fastq_path, output_dir)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    reads_per_block, total_reads = LY.get_reads_num_per_block(fastq_path, block_size)
    if reads_per_block <= 0:
        reads_per_block = 1

    rules_dict = LY.init_rules_dict()      # base 规则表，全文件复用
    rules_dict_q = LY.init_rules_dict_q()  # quality 规则表，全文件复用
    block_count = 0

    with open(fastq_path, 'r') as f:
        parser = SeqIO.parse(f, 'fastq')
        while True:
            records = []
            for _ in range(reads_per_block):
                try:
                    records.append(next(parser))
                except StopIteration:
                    break
            if not records:
                break

            block_count += 1
            g_block, q_block, id_block = LY.process_records(records, rules_dict, rules_dict_q)
            if g_block is None:
                continue
            if save:
                LY.save_intermediate_files(g_block, q_block, id_block, output_path, block_count)
            LY.back_compress_worker(g_block, q_block, id_block, lpaq8_path, output_path, save, block_count)

    merge_parts(output_path, block_count)
    return output_path, total_reads, reads_per_block


def compress_and_collect_metrics(input_file, output_dir):
    start = time.time()
    input_size = get_file_size(input_file)
    file_name = os.path.splitext(os.path.basename(input_file))[0]

    proc = monitor_memory_start()
    peak_mem = monitor_memory_peak_mb(proc, 0)

    lpaq8_path = str((Path(__file__).resolve().parent / 'lpaq8').resolve())
    if COMPRESSOR.lower() == 'lossy':
        out_path, total_reads, reads_per_block = compress_lossy_cross_block(
            input_file, output_dir, lpaq8_path, SAVE_INTERMEDIATE, BLOCK_SIZE
        )
    else:
        out_path, total_reads, reads_per_block = compress_lossless_cross_block(
            input_file, output_dir, lpaq8_path, SAVE_INTERMEDIATE, BLOCK_SIZE
        )

    peak_mem = monitor_memory_peak_mb(proc, peak_mem)
    elapsed = time.time() - start
    output_size = get_file_size(out_path)
    ratio = input_size / output_size if output_size > 0 else 0
    speed = (input_size / 1024 / 1024) / elapsed if elapsed > 0 else 0

    return {
        'file_name': file_name,
        'compressor': COMPRESSOR,
        'cross_block_learning': True,
        'input_size_mb': input_size / 1024 / 1024,
        'output_size_mb': output_size / 1024 / 1024,
        'compression_ratio': ratio,
        'compression_time_s': elapsed,
        'compression_speed_mbs': speed,
        'avg_cpu_percent': 0.0,
        'max_memory_mb': peak_mem,
        'threads': 1,
        'block_size_bytes': BLOCK_SIZE,
        'total_reads': total_reads,
        'reads_per_block': reads_per_block,
    }


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    csv_path = os.path.join(OUTPUT_DIR, f"ablation_cross_block_learning_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv")
    header = [
        'file_name', 'compressor', 'cross_block_learning',
        'input_size_mb', 'output_size_mb', 'compression_ratio',
        'compression_time_s', 'compression_speed_mbs',
        'avg_cpu_percent', 'max_memory_mb',
        'threads', 'block_size_bytes', 'total_reads', 'reads_per_block'
    ]

    with open(csv_path, 'w', newline='') as f:
        csv.writer(f).writerow(header)

    for file_path in Path(INPUT_DIR).glob('*.fastq'):
        print(f'Processing {file_path} ...')
        try:
            metrics = compress_and_collect_metrics(str(file_path), OUTPUT_DIR)
            with open(csv_path, 'a', newline='') as f:
                csv.writer(f).writerow([metrics[k] for k in header])
            print(f"Done: {file_path.name}, ratio={metrics['compression_ratio']:.3f}")
        except Exception as e:
            print(f'Error: {file_path} -> {e}')


if __name__ == '__main__':
    main()
