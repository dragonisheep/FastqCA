"""Ablation experiment 2: run FastqCA with thread count = 1.

Uses the same block_size as previous experiments (default 128MB in main_new.py).
"""

import os
import csv
import time
import subprocess
from pathlib import Path
from datetime import datetime

import psutil

# ===== 配置区（可按需修改） =====
COMPRESSOR = 'LossLess'  # 可改成 'Lossy'
THREAD_COUNT = '1'
INPUT_DIR = '/media/compress/新加卷1/New_Test'
OUTPUT_DIR = '/media/compress/新加卷/output/New_Test_LossLess/Ablation_thread1'
BLOCK_SIZE = str(128 * 1024 * 1024)  # 与原实验保持一致
# ===========================


def get_file_size(file_path):
    return os.path.getsize(file_path)


def monitor_process(process):
    cpu_percentages = []
    max_total_memory = 0

    try:
        p = psutil.Process(process.pid)
        while process.poll() is None:
            try:
                all_procs = [p] + p.children(recursive=True)
                total_rss = 0
                total_cpu = 0
                for proc in all_procs:
                    try:
                        total_rss += proc.memory_info().rss / 1024 / 1024
                        total_cpu += proc.cpu_percent()
                    except Exception:
                        continue
                cpu_percentages.append(total_cpu)
                if total_rss > max_total_memory:
                    max_total_memory = total_rss
            except Exception:
                break
            time.sleep(0.05)
    except Exception:
        pass

    return {
        'avg_cpu': sum(cpu_percentages) / len(cpu_percentages) if cpu_percentages else 0,
        'max_memory': max_total_memory,
    }


def compress_and_collect_metrics(input_file, output_dir):
    start_time = time.time()
    input_size = get_file_size(input_file)
    file_name = os.path.splitext(os.path.basename(input_file))[0]

    cmd = [
        'python', 'main_new.py',
        '--compressor', COMPRESSOR,
        '--input_path', input_file,
        '--output_path', output_dir,
        '--mode', 'c',
        '--threads', THREAD_COUNT,
        '--block_size', BLOCK_SIZE,
    ]

    print(f"Executing: {' '.join(cmd)}")
    process = subprocess.Popen(cmd)
    metrics = monitor_process(process)
    process.wait()
    if process.returncode != 0:
        raise RuntimeError(f"Compression failed for {input_file}")

    end_time = time.time()
    compression_time = end_time - start_time
    time.sleep(1)

    suffix = '.lossy' if COMPRESSOR.lower() == 'lossy' else '.lossless'
    possible_outputs = [
        os.path.join(output_dir, file_name + suffix),
        os.path.join(output_dir, file_name + '.fastq' + suffix),
        os.path.join(output_dir, file_name),
    ]
    output_size = 0
    for f in possible_outputs:
        if os.path.exists(f):
            output_size = get_file_size(f)
            break

    ratio = input_size / output_size if output_size > 0 else 0
    speed = (input_size / 1024 / 1024) / compression_time if compression_time > 0 else 0

    return {
        'file_name': file_name,
        'input_size_mb': input_size / 1024 / 1024,
        'output_size_mb': output_size / 1024 / 1024,
        'compression_ratio': ratio,
        'compression_time_s': compression_time,
        'compression_speed_mbs': speed,
        'avg_cpu_percent': metrics['avg_cpu'],
        'max_memory_mb': metrics['max_memory'],
        'threads': int(THREAD_COUNT),
        'block_size_bytes': int(BLOCK_SIZE),
        'compressor': COMPRESSOR,
    }


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    csv_path = os.path.join(OUTPUT_DIR, f"ablation_thread1_metrics_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv")
    header = [
        'file_name', 'input_size_mb', 'output_size_mb', 'compression_ratio',
        'compression_time_s', 'compression_speed_mbs', 'avg_cpu_percent',
        'max_memory_mb', 'threads', 'block_size_bytes', 'compressor'
    ]

    with open(csv_path, 'w', newline='') as f:
        csv.writer(f).writerow(header)

    for file_path in Path(INPUT_DIR).glob('*.fastq'):
        print(f'Processing {file_path}...')
        try:
            metrics = compress_and_collect_metrics(str(file_path), OUTPUT_DIR)
            with open(csv_path, 'a', newline='') as f:
                csv.writer(f).writerow([metrics[k] for k in header])
            print(f'Completed: {file_path.name}')
        except Exception as e:
            print(f'Error: {e}')


if __name__ == '__main__':
    main()
