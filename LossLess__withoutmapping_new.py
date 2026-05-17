import csv
import os
import subprocess
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import psutil

THREAD_COUNT = '4'
BLOCK_SIZE = str(128 * 1024 * 1024)


def get_file_size(file_path):
    return os.path.getsize(file_path)


def monitor_process(process):
    cpu_percentages = []
    max_total_memory = 0
    peak_snapshot = ''

    try:
        p = psutil.Process(process.pid)
        while process.poll() is None:
            all_procs = [p] + p.children(recursive=True)
            total_rss = 0
            total_cpu = 0
            proc_counts = defaultdict(int)

            for proc in all_procs:
                try:
                    rss_mb = proc.memory_info().rss / 1024 / 1024
                    total_rss += rss_mb
                    total_cpu += proc.cpu_percent()
                    proc_counts[proc.name()] += 1
                except Exception:
                    continue

            cpu_percentages.append(total_cpu)
            if total_rss > max_total_memory:
                max_total_memory = total_rss
                counts_str = ', '.join([f'{k}:{v}' for k, v in proc_counts.items()])
                peak_snapshot = f'内存: {total_rss:.2f} MB | 进程: {len(all_procs)} ({counts_str})'
            time.sleep(0.05)
    except Exception:
        pass

    print(f"\n{'=' * 50}\nDEBUG: 监控结束。\n【峰值详情】: {peak_snapshot}\n{'=' * 50}\n")
    return {
        'avg_cpu': sum(cpu_percentages) / len(cpu_percentages) if cpu_percentages else 0,
        'max_memory': max_total_memory,
    }


def compress_and_collect_metrics(input_file, output_dir):
    start_time = time.time()
    input_size = get_file_size(input_file)
    file_name = os.path.splitext(os.path.basename(input_file))[0]

    cmd = [
        'python', 'LossLess_withoutmapping_thread.py',
        '--input_path', input_file,
        '--output_path', output_dir,
        '--threads', THREAD_COUNT,
        '--block_size', BLOCK_SIZE,
    ]

    print(f"Executing: {' '.join(cmd)}")
    process = subprocess.Popen(cmd)
    metrics = monitor_process(process)
    process.wait()
    if process.returncode != 0:
        raise RuntimeError(f'compression failed: returncode={process.returncode}')

    compression_time = time.time() - start_time
    out_file = os.path.join(output_dir, file_name + '.withoutmapping.lossless')
    output_size = get_file_size(out_file) if os.path.exists(out_file) else 0

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
        'threads': THREAD_COUNT,
        'block_size': int(BLOCK_SIZE),
    }


def main():
    input_dir = '/media/compress/新加卷1/New_Test'
    output_dir = '/media/compress/新加卷/output/New_Test_LossLess/FastqCA_withoutmapping'
    os.makedirs(output_dir, exist_ok=True)

    csv_path = os.path.join(output_dir, f'compression_metrics_withoutmapping_{datetime.now().strftime("%Y%m%d_%H%M%S")}.csv')
    csv_header = [
        'file_name', 'input_size_mb', 'output_size_mb', 'compression_ratio', 'compression_time_s',
        'compression_speed_mbs', 'avg_cpu_percent', 'max_memory_mb', 'threads', 'block_size'
    ]

    with open(csv_path, 'w', newline='') as f:
        csv.writer(f).writerow(csv_header)

    for file_path in Path(input_dir).glob('*.fastq'):
        print(f'Processing {file_path}...')
        try:
            metrics = compress_and_collect_metrics(str(file_path), output_dir)
            with open(csv_path, 'a', newline='') as f:
                csv.writer(f).writerow([metrics[key] for key in csv_header])
            print(f"Completed processing {file_path}\n")
        except Exception as e:
            print(f'Error: {e}')


if __name__ == '__main__':
    main()
