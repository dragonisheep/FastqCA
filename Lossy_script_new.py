"""Lossy compression experiment runner.

Runs lossy compression with the same input/output directories as the lossless
experimental script (LossLess_script_new.py) and records CPU/memory metrics.
"""

import os
import subprocess
import time
import psutil
import csv
from pathlib import Path
from datetime import datetime
from collections import defaultdict

# 与无损实验脚本保持一致的线程数与输入/输出目录
THREAD_COUNT = '4'
INPUT_DIR = '/media/compress/新加卷/Temp_Test'
OUTPUT_DIR = '/media/compress/新加卷/output/Temp_Test/FastqCA_lossy'

def get_file_size(file_path):
    return os.path.getsize(file_path)


def monitor_process(process):
    cpu_percentages = []
    memory_usages = []
    max_total_memory = 0
    peak_snapshot = ""

    try:
        p = psutil.Process(process.pid)
        print(f"DEBUG: 主进程启动, PID: {p.pid}")

        step_count = 0
        while process.poll() is None:
            try:
                all_procs = [p] + p.children(recursive=True)
                total_rss = 0
                total_cpu = 0
                mem_breakdown = defaultdict(float)
                proc_counts = defaultdict(int)

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
                        proc_counts[group] += 1
                    except:  # noqa: E722
                        continue

                memory_usages.append(total_rss)
                cpu_percentages.append(total_cpu)

                if total_rss > max_total_memory:
                    max_total_memory = total_rss
                    counts_str = ", ".join([f"{k}:{v}" for k, v in proc_counts.items()])
                    peak_snapshot = f"内存: {total_rss:.2f} MB | 进程: {len(all_procs)} ({counts_str})"

                if step_count % 20 == 0 and len(all_procs) > 1:
                    breakdown_str = " | ".join([f"{k}: {v:.1f}MB" for k, v in mem_breakdown.items()])
                    # print(f"   [监控] {breakdown_str} (总: {total_rss:.1f}MB)")

                step_count += 1
            except:  # noqa: E722
                break
            time.sleep(0.05)
    except:  # noqa: E722
        pass

    print(f"\n{'=' * 50}\nDEBUG: 监控结束。\n【峰值详情】: {peak_snapshot}\n{'=' * 50}\n")
    return {'avg_cpu': sum(cpu_percentages) / len(cpu_percentages) if cpu_percentages else 0,
            'max_memory': max_total_memory}


def compress_and_collect_metrics(input_file, output_dir):
    start_time = time.time()
    input_size = get_file_size(input_file)
    file_name = os.path.splitext(os.path.basename(input_file))[0]

    cmd = [
        'python', 'main_new.py',
        '--compressor', 'Lossy',
        '--input_path', input_file,
        '--output_path', output_dir,
        '--mode', 'c',
        '--threads', THREAD_COUNT,
    ]

    print(f"Executing: {' '.join(cmd)}")
    process = subprocess.Popen(cmd)
    metrics = monitor_process(process)
    process.wait()

    end_time = time.time()
    compression_time = end_time - start_time
    time.sleep(1)

    possible_outputs = [
        os.path.join(output_dir, file_name + '.lossy'),
        os.path.join(output_dir, file_name + '.fastq.lossy'),
        os.path.join(output_dir, file_name)
    ]
    output_size = 0
    for f in possible_outputs:
        if os.path.exists(f):
            output_size = get_file_size(f)
            break

    ratio = input_size / output_size if output_size > 0 else 0
    speed = (input_size / 1024 / 1024) / compression_time if compression_time > 0 else 0

    return {
        'file_name': file_name, 'input_size_mb': input_size / 1024 / 1024, 'output_size_mb': output_size / 1024 / 1024,
        'compression_ratio': ratio, 'compression_time_s': compression_time, 'compression_speed_mbs': speed,
        'avg_cpu_percent': metrics['avg_cpu'], 'max_memory_mb': metrics['max_memory']
    }


def main():
    input_dir = INPUT_DIR
    output_dir = OUTPUT_DIR
    os.makedirs(output_dir, exist_ok=True)

    csv_path = os.path.join(output_dir, f'compression_metrics_{datetime.now().strftime("%Y%m%d_%H%M%S")}.csv')
    csv_header = ['file_name', 'input_size_mb', 'output_size_mb', 'compression_ratio', 'compression_time_s',
                  'compression_speed_mbs', 'avg_cpu_percent', 'max_memory_mb']

    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(csv_header)

    for file_path in Path(input_dir).glob('*.fastq'):
        print(f'Processing {file_path}...')
        try:
            metrics = compress_and_collect_metrics(str(file_path), output_dir)
            print(f"Result -> Peak Mem: {metrics['max_memory_mb']:.2f} MB")
            with open(csv_path, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([metrics[key] for key in csv_header])
            print(f'Completed processing {file_path}\n')
        except Exception as e:
            print(f"Error: {e}")


if __name__ == '__main__':
    main()
