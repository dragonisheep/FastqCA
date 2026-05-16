"""Lossless control experiment runner with block size 32MB."""

import os
import subprocess
import time
import psutil
import csv
from pathlib import Path
from datetime import datetime
from collections import defaultdict

THREAD_COUNT = '4'
BLOCK_SIZE = str(32 * 1024 * 1024)
INPUT_DIR = '/media/compress/新加卷1/New_Test'
OUTPUT_DIR = '/media/compress/新加卷/output/New_Test_LossLess/FastqCA_block_32MB'


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
                mem_breakdown = defaultdict(float)

                for proc in all_procs:
                    try:
                        mem_info = proc.memory_info()
                        name = proc.name()
                        rss_mb = mem_info.rss / 1024 / 1024
                        total_rss += rss_mb
                        total_cpu += proc.cpu_percent()

                        if "python" in name.lower():
                            group = "Python"
                        elif "lpaq8" in name.lower():
                            group = "lpaq8"
                        else:
                            group = "Other"
                        mem_breakdown[group] += rss_mb
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
        '--compressor', 'LossLess',
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
        raise RuntimeError(f"Compression failed for {input_file}, returncode={process.returncode}")

    compression_time = time.time() - start_time
    time.sleep(1)

    possible_outputs = [
        os.path.join(output_dir, file_name + '.lossless'),
        os.path.join(output_dir, file_name + '.fastq.lossless'),
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
    }


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    csv_path = os.path.join(OUTPUT_DIR, f'compression_metrics_block_32MB_{datetime.now().strftime("%Y%m%d_%H%M%S")}.csv')
    header = [
        'file_name', 'input_size_mb', 'output_size_mb', 'compression_ratio',
        'compression_time_s', 'compression_speed_mbs',
        'avg_cpu_percent', 'max_memory_mb', 'threads', 'block_size_bytes'
    ]

    with open(csv_path, 'w', newline='') as f:
        csv.writer(f).writerow(header)

    for file_path in Path(INPUT_DIR).glob('*.fastq'):
        print(f'Processing {file_path}...')
        try:
            metrics = compress_and_collect_metrics(str(file_path), OUTPUT_DIR)
            with open(csv_path, 'a', newline='') as f:
                csv.writer(f).writerow([metrics[k] for k in header])
            print(f"Completed {file_path.name}, ratio={metrics['compression_ratio']:.3f}")
        except Exception as e:
            print(f'Error: {e}')


if __name__ == '__main__':
    main()
