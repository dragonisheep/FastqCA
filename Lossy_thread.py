import argparse
import os
import shutil
import sys
import time
import re
import gc
import mmap
import struct
from tqdm import tqdm
from PIL import Image, UnidentifiedImageError
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from collections import defaultdict
from itertools import product
from Bio import SeqIO
import multiprocessing

from lpaq8 import compress_file, decompress_file

Image.MAX_IMAGE_PIXELS = None

# 映射表
base_to_gray = {'A': 32, 'T': 64, 'G': 192, 'C': 224, 'N': 0}
gray_to_base = {32: 'A', 64: 'T', 192: 'G', 224: 'C', 0: 'N'}


def Q4(qsc):
    if qsc <= 7:
        return 5
    if qsc <= 13:
        return 12
    if qsc <= 19:
        return 18
    return 24


def find_delimiters(identifier):
    return re.findall(r'[.:_\s=/-]', identifier)


def split_identifier(identifier, delimiters):
    return re.split(r'[.:_\s=/-]', identifier)


def generate_regex(delimiters):
    regex = ""
    count = 1
    for delimiter in delimiters:
        regex += f"{{{{T{count}}}}}"
        regex += delimiter
        count += 1
    regex += f"{{{{T{count}}}}}"
    return regex


def init_rules_dict():
    values = [0, 32, 64, 192, 224]
    combinations = list(product(values, repeat=4))
    rules_dict = defaultdict(int)
    for combination in combinations:
        rules_dict[combination] = 0
    return rules_dict


def init_rules_dict_q():
    values = [5, 12, 18, 24]
    combinations = list(product(values, repeat=4))
    rules_dict_q = defaultdict(int)
    for combination in combinations:
        rules_dict_q[combination] = 0
    return rules_dict_q


def get_reads_num_per_block(fastq_path, block_size):
    with open(fastq_path, 'r') as file:
        try:
            first_record = next(SeqIO.parse(file, "fastq"))
            read_length = len(first_record.seq)
        except StopIteration:
            return 0, 0
    bytes_per_read = read_length * 2 if read_length else 1
    reads_per_block = block_size // bytes_per_read if block_size else 1
    if reads_per_block <= 0:
        reads_per_block = 1
    total_reads = os.path.getsize(fastq_path) // bytes_per_read
    return reads_per_block, total_reads


# --- Compression Logic ---

def generate_g_prime(G, rules_dict):
    G_prime = np.zeros_like(G)
    rows, cols = G.shape
    for i in range(rows):
        for j in range(cols):
            center = G[i, j]
            up = G[i - 1, j] if i != 0 else 0
            left = G[i, j - 1] if j != 0 else 0
            left_up = G[i - 1, j - 1] if i != 0 and j != 0 else 0
            matched_rule = (up, left_up, left, center)
            candidates = [(up, left_up, left, v) for v in [32, 224, 192, 64, 0]]
            top_rule = max(candidates, key=lambda r: rules_dict[r])
            G_prime[i, j] = 1 if top_rule[3] == center else center
            rules_dict[matched_rule] += 1
    return G_prime


def generate_q_prime(Q, rules_dict_q):
    Q_prime = np.zeros_like(Q)
    rows, cols = Q.shape
    for i in range(rows):
        for j in range(cols):
            center = Q[i, j]
            up = Q[i - 1, j] if i != 0 else 0
            left = Q[i, j - 1] if j != 0 else 0
            left_up = Q[i - 1, j - 1] if i != 0 and j != 0 else 0
            matched_rule = (up, left_up, left, center)
            candidates = [(up, left_up, left, v) for v in [5, 12, 18, 24]]
            top_rule = max(candidates, key=lambda r: rules_dict_q[r])
            Q_prime[i, j] = 1 if top_rule[3] == center else center
            rules_dict_q[matched_rule] += 1
    return Q_prime


def process_records(records, rules_dict, rules_dict_q):
    id_block, base_image_block, quality_block = [], [], []
    for record in records:
        id_str = record.description
        delimiters = find_delimiters(id_str)
        tokens = split_identifier(id_str, delimiters)
        regex = generate_regex(delimiters)
        id_block.append((tokens, regex))
        base_gray_values = [base_to_gray.get(base, 0) for base in record.seq]
        base_image_block.append(base_gray_values)
        quality_gray_values = [Q4(q) for q in record.letter_annotations["phred_quality"]]
        quality_block.append(quality_gray_values)

    if not base_image_block:
        return None, None, None
    g_prime = generate_g_prime(np.array(base_image_block, dtype=np.uint8), rules_dict)
    q_prime = generate_q_prime(np.array(quality_block, dtype=np.uint8), rules_dict_q)
    g_prime_img = Image.fromarray(g_prime.astype(np.uint8))
    q_prime_img = Image.fromarray(q_prime.astype(np.uint8))
    return g_prime_img, q_prime_img, id_block


def save_intermediate_files(g_block, q_block, id_block, output_path, block_count):
    front_dir = os.path.join(os.path.dirname(output_path), "front_compressed")
    os.makedirs(front_dir, exist_ok=True)
    g_block.save(os.path.join(front_dir, f'chunk_{block_count}_base.tiff'))
    q_block.save(os.path.join(front_dir, f'chunk_{block_count}_quality.tiff'))
    with open(os.path.join(front_dir, f"chunk_{block_count}_id_tokens.txt"), 'w') as f1, \
            open(os.path.join(front_dir, f"chunk_{block_count}_id_regex.txt"), 'w') as f2:
        for tokens, regex in id_block:
            f1.write(' '.join(tokens) + '\n')
            f2.write(regex + '\n')


def compress_worker_subprocess(temp_input_path, temp_output_path, lpaq8_path):
    process = compress_file(temp_input_path, temp_output_path, lpaq8_path)
    if process:
        process.wait()


def write_safe_chunk(output_file, tag, data):
    output_file.write(tag)
    output_file.write(struct.pack('<Q', len(data)))
    output_file.write(data)


def back_compress_worker(g_block, q_block, id_block, lpaq8_path, output_path, save, block_count):
    part_output_path = os.path.join(os.path.dirname(output_path), f"chunk_{block_count}.part")
    temp_prefix = os.path.join(os.path.dirname(output_path), f"temp_proc_{block_count}")
    temp_input_path, temp_output_path = f"{temp_prefix}_input", f"{temp_prefix}_output"

    back_dir = os.path.join(os.path.dirname(output_path), "back_compressed")
    if save:
        os.makedirs(back_dir, exist_ok=True)

    try:
        with open(part_output_path, "wb") as output_file:
            # ID Regex
            with open(temp_input_path, "w") as f:
                for _, regex in id_block:
                    f.write(regex + '\n')
            compress_worker_subprocess(temp_input_path, temp_output_path, lpaq8_path)
            if not os.path.exists(temp_output_path):
                raise RuntimeError("LPAQ8 compression failed for ID Regex")
            with open(temp_output_path, "rb") as f:
                data = f.read()
            write_safe_chunk(output_file, b"%id_regex%", data)
            if save:
                with open(os.path.join(back_dir, f"chunk_{block_count}_id_regex.lpaq8"), "wb") as sf:
                    sf.write(data)

            # ID Tokens
            with open(temp_input_path, "w") as f:
                for tokens, _ in id_block:
                    f.write(' '.join(tokens) + '\n')
            compress_worker_subprocess(temp_input_path, temp_output_path, lpaq8_path)
            if not os.path.exists(temp_output_path):
                raise RuntimeError("LPAQ8 compression failed for ID Tokens")
            with open(temp_output_path, "rb") as f:
                data = f.read()
            write_safe_chunk(output_file, b"%id_tokens%", data)
            if save:
                with open(os.path.join(back_dir, f"chunk_{block_count}_id_tokens.lpaq8"), "wb") as sf:
                    sf.write(data)

            # Base image
            g_block.save(temp_input_path, format="tiff")
            compress_worker_subprocess(temp_input_path, temp_output_path, lpaq8_path)
            if not os.path.exists(temp_output_path):
                raise RuntimeError("LPAQ8 compression failed for Base")
            with open(temp_output_path, "rb") as f:
                data = f.read()
            write_safe_chunk(output_file, b"%base_g_prime%", data)
            if save:
                with open(os.path.join(back_dir, f"chunk_{block_count}_base_g_prime.lpaq8"), "wb") as sf:
                    sf.write(data)

            # Quality image
            q_block.save(temp_input_path, format="tiff")
            compress_worker_subprocess(temp_input_path, temp_output_path, lpaq8_path)
            if not os.path.exists(temp_output_path):
                raise RuntimeError("LPAQ8 compression failed for Quality")
            with open(temp_output_path, "rb") as f:
                data = f.read()
            write_safe_chunk(output_file, b"%quality%", data)
            if save:
                with open(os.path.join(back_dir, f"chunk_{block_count}_quality.lpaq8"), "wb") as sf:
                    sf.write(data)
    finally:
        if os.path.exists(temp_input_path):
            os.remove(temp_input_path)
        if os.path.exists(temp_output_path):
            os.remove(temp_output_path)
    return block_count


def process_block_task_from_file(temp_chunk_path, block_count, output_path, lpaq8_path, save):
    records = []
    try:
        gc.collect()
        with open(temp_chunk_path, 'r') as f:
            records = list(SeqIO.parse(f, "fastq"))
        if not records:
            return block_count

        rules_dict = init_rules_dict()
        rules_dict_q = init_rules_dict_q()
        g_block, q_block, id_block = process_records(records, rules_dict, rules_dict_q)
        del records, rules_dict, rules_dict_q
        gc.collect()

        if g_block is None:
            return block_count
        if save:
            save_intermediate_files(g_block, q_block, id_block, output_path, block_count)
        back_compress_worker(g_block, q_block, id_block, lpaq8_path, output_path, save, block_count)
    finally:
        if os.path.exists(temp_chunk_path):
            try:
                os.remove(temp_chunk_path)
            except Exception:
                pass
        gc.collect()
    return block_count


def merge_parts(output_path, total_blocks):
    missing_parts = []
    tqdm.write(f"info：正在合并 {total_blocks} 个数据块...")
    with open(output_path, "wb") as final_file:
        for i in range(1, total_blocks + 1):
            part_path = os.path.join(os.path.dirname(output_path), f"chunk_{i}.part")
            if os.path.exists(part_path):
                with open(part_path, "rb") as part_file:
                    shutil.copyfileobj(part_file, final_file)
                os.remove(part_path)
            else:
                missing_parts.append(part_path)
        final_file.write(b"%eof")
    if missing_parts:
        raise FileNotFoundError(f"以下分块缺失，导致压缩结果不完整: {missing_parts}")
    tqdm.write("info：合并完成")


def compress_multithread(fastq_path, output_path, lpaq8_path, save, block_size, max_workers):
    output_path = get_output_path(fastq_path, output_path)
    out_dir = os.path.dirname(output_path)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    reads_per_block, total_reads = get_reads_num_per_block(fastq_path, block_size)
    total_block_count = int(np.ceil(os.path.getsize(fastq_path) / block_size)) if block_size > 0 else 1

    tqdm.write(f"info：多进程模式 (进程数={max_workers}) | 模式: 安全长度头封装")
    temp_chunk_dir = os.path.join(out_dir, "temp_chunks")
    os.makedirs(temp_chunk_dir, exist_ok=True)

    records, read_count_per_block, block_count = [], 0, 1
    with multiprocessing.Pool(processes=max_workers, maxtasksperchild=1) as pool:
        results = []
        errors = []
        tqdm.write("info：正在读取 FASTQ 并分发任务...")

        with open(fastq_path, 'r') as file:
            for record in tqdm(SeqIO.parse(file, "fastq"), desc="Chunking", total=total_reads, unit="reads"):
                records.append(record)
                read_count_per_block += 1
                if read_count_per_block >= reads_per_block:
                    temp_chunk_path = os.path.join(temp_chunk_dir, f"chunk_src_{block_count}.fastq")
                    SeqIO.write(records, temp_chunk_path, "fastq")
                    res = pool.apply_async(process_block_task_from_file,
                                           (temp_chunk_path, block_count, output_path, lpaq8_path, save))
                    results.append(res)
                    if len(results) > max_workers * 2:
                        results = [r for r in results if not r.ready()]
                        if len(results) > max_workers * 3:
                            time.sleep(1)
                    records, read_count_per_block = [], 0
                    block_count += 1
            if records:
                temp_chunk_path = os.path.join(temp_chunk_dir, f"chunk_src_{block_count}.fastq")
                SeqIO.write(records, temp_chunk_path, "fastq")
                res = pool.apply_async(process_block_task_from_file,
                                       (temp_chunk_path, block_count, output_path, lpaq8_path, save))
                results.append(res)
            else:
                block_count -= 1

        pool.close()
        tqdm.write("info：等待所有子进程完成...")
        for res in tqdm(results, desc="Waiting Workers"):
            try:
                res.get()
            except Exception as e:
                errors.append(e)
        pool.join()

        if errors:
            raise RuntimeError(f"检测到 {len(errors)} 个压缩块处理失败: {errors}")

    merge_parts(output_path, block_count)
    try:
        shutil.rmtree(temp_chunk_dir)
    except Exception:
        pass


# --- Decompression Logic ---

def monitor(process, temp_input_path, temp_output_path):
    while True:
        if os.path.exists(temp_input_path) and os.path.exists(temp_output_path):
            if process.poll() is not None:
                break
        time.sleep(0.1)


def decompress_with_monitor(temp_input_path, temp_output_path, lpaq8_path):
    process = decompress_file(temp_input_path, temp_output_path, lpaq8_path)
    monitor(process, temp_input_path, temp_output_path)


def process_compressed_block(output_path, lpaq8_path, id_regex_data, id_tokens_data, g_prime_data, quality_data, save,
                             block_count):
    back_compress_dir = os.path.join(os.path.dirname(output_path), "back_compressed")
    front_compress_dir = os.path.join(os.path.dirname(output_path), "front_compressed")
    if save:
        os.makedirs(back_compress_dir, exist_ok=True)
        os.makedirs(front_compress_dir, exist_ok=True)

    temp_prefix = os.path.join(os.path.dirname(output_path), f"temp_dec_{block_count}")
    temp_input_path, temp_output_path = f"{temp_prefix}_in", f"{temp_prefix}_out"
    id_regex, id_tokens, g_prime, quality = None, None, None, None

    try:
        # ID regex
        with open(temp_input_path, "wb") as temp_input_file:
            temp_input_file.write(id_regex_data)
        decompress_with_monitor(temp_input_path, temp_output_path, lpaq8_path)
        with open(temp_output_path, "r") as temp_output_file:
            id_regex = [line.strip() for line in temp_output_file.readlines()]
        if save:
            shutil.copy(temp_input_path, os.path.join(back_compress_dir, f"chunk_{block_count}_id_regex.lpaq8"))
            shutil.copy(temp_output_path, os.path.join(front_compress_dir, f"chunk_{block_count}_id_regex.txt"))

        # ID tokens
        with open(temp_input_path, "wb") as temp_input_file:
            temp_input_file.write(id_tokens_data)
        decompress_with_monitor(temp_input_path, temp_output_path, lpaq8_path)
        with open(temp_output_path, "r") as temp_output_file:
            id_tokens = [line.strip() for line in temp_output_file.readlines()]
        if save:
            shutil.copy(temp_input_path, os.path.join(back_compress_dir, f"chunk_{block_count}_id_tokens.lpaq8"))
            shutil.copy(temp_output_path, os.path.join(front_compress_dir, f"chunk_{block_count}_id_tokens.txt"))

        id_block = zip(id_tokens, id_regex)

        # Quality
        with open(temp_input_path, "wb") as temp_input_file:
            temp_input_file.write(quality_data)
        decompress_with_monitor(temp_input_path, temp_output_path, lpaq8_path)
        try:
            with Image.open(temp_output_path) as img:
                quality = img.copy()
        except UnidentifiedImageError:
            tqdm.write(f"无法识别的图像文件: {temp_output_path}")
            raise
        if save:
            shutil.copy(temp_input_path, os.path.join(back_compress_dir, f'chunk_{block_count}_quality.lpaq8'))
            shutil.copy(temp_output_path, os.path.join(front_compress_dir, f'chunk_{block_count}_quality.tiff'))

        # Base
        with open(temp_input_path, "wb") as temp_input_file:
            temp_input_file.write(g_prime_data)
        decompress_with_monitor(temp_input_path, temp_output_path, lpaq8_path)
        try:
            with Image.open(temp_output_path) as img:
                g_prime = img.copy()
        except UnidentifiedImageError:
            tqdm.write(f"无法识别的图像文件: {temp_output_path}")
            raise
        if save:
            shutil.copy(temp_input_path, os.path.join(back_compress_dir, f'chunk_{block_count}_base_g_prime.lpaq8'))
            shutil.copy(temp_output_path, os.path.join(front_compress_dir, f'chunk_{block_count}_base_g_prime.tiff'))

    finally:
        if os.path.exists(temp_input_path):
            try:
                os.remove(temp_input_path)
            except Exception:
                pass
        if os.path.exists(temp_output_path):
            try:
                os.remove(temp_output_path)
            except Exception:
                pass
    return id_block, g_prime, quality


def reconstruct_id(tokens, regex):
    reconstructed_ids = []
    for t, r in zip(tokens, regex):
        id_str = r
        token_list = t.split()
        for i, token in enumerate(token_list):
            id_str = id_str.replace(f"{{{{T{i+1}}}}}", token)
        reconstructed_ids.append(id_str)
    return reconstructed_ids


def reconstruct_g_from_g_prime(g_prime_array, rules_dict):
    De_g = np.zeros_like(g_prime_array)
    for i in range(g_prime_array.shape[0]):
        for j in range(g_prime_array.shape[1]):
            center = g_prime_array[i, j]
            up = De_g[i - 1, j] if i != 0 else 0
            left = De_g[i, j - 1] if j != 0 else 0
            left_up = De_g[i - 1, j - 1] if i != 0 and j != 0 else 0
            matched_rules = [(up, left_up, left, v) for v in [32, 224, 192, 64, 0]]
            top_rule = max(matched_rules, key=lambda rule: rules_dict[rule])
            De_g[i, j] = top_rule[3] if g_prime_array[i, j] == 1 else center
            matched_rule = (up, left_up, left, De_g[i, j])
            rules_dict[matched_rule] += 1
    return De_g


def reconstruct_q_from_q_prime(q_prime_array, rules_dict_q):
    De_q = np.zeros_like(q_prime_array)
    for i in range(q_prime_array.shape[0]):
        for j in range(q_prime_array.shape[1]):
            center = q_prime_array[i, j]
            up = De_q[i - 1, j] if i != 0 else 0
            left = De_q[i, j - 1] if j != 0 else 0
            left_up = De_q[i - 1, j - 1] if i != 0 and j != 0 else 0
            matched_rules = [(up, left_up, left, v) for v in [5, 12, 18, 24]]
            top_rule = max(matched_rules, key=lambda rule: rules_dict_q[rule])
            De_q[i, j] = top_rule[3] if q_prime_array[i, j] == 1 else center
            matched_rule = (up, left_up, left, De_q[i, j])
            rules_dict_q[matched_rule] += 1
    return De_q


def reconstruct_base_and_quality(g_prime_img, q_prime_img):
    bases = []
    qualities = []
    rules_dict = defaultdict(int)
    rules_dict_q = defaultdict(int)
    g_prime_array = np.array(g_prime_img)
    bases_array = reconstruct_g_from_g_prime(g_prime_array, rules_dict)
    q_prime_array = np.array(q_prime_img)
    quality_array = reconstruct_q_from_q_prime(q_prime_array, rules_dict_q)
    for i in range(bases_array.shape[0]):
        base_str = ''.join([gray_to_base[pixel] for pixel in bases_array[i]])
        quality_scores = [q for q in quality_array[i]]
        bases.append(base_str)
        qualities.append(quality_scores)
    return bases, qualities


def reconstruct_fastq(output_path, id_block, g_prime_img, quality_img, is_first_block=False):
    records = []
    final_fastq_path = output_path if output_path.endswith('.fastq') else os.path.splitext(output_path)[0] + '.fastq'
    id_block = list(id_block)
    id_tokens = [item[0] for item in id_block]
    id_regex = [item[1] for item in id_block]
    ids = reconstruct_id(id_tokens, id_regex)
    g_prime, quality = reconstruct_base_and_quality(g_prime_img, quality_img)
    for i in range(len(ids)):
        seq = Seq(g_prime[i])
        record = SeqRecord(seq, id=ids[i], description="")
        record.letter_annotations["phred_quality"] = quality[i]
        records.append(record)
    mode = 'w' if is_first_block else 'a'
    with open(final_fastq_path, mode) as output_handle:
        SeqIO.write(records, output_handle, 'fastq')


def read_chunk_safe(mmap_obj, tag):
    start_pos = mmap_obj.tell()
    header = mmap_obj.read(len(tag))
    if header != tag:
        raise RuntimeError(f"解压时在偏移 {start_pos} 处未找到预期的标记 {tag!r}")
    size_bytes = mmap_obj.read(8)
    if len(size_bytes) < 8:
        raise RuntimeError("解压时无法读取长度字段，压缩文件可能已损坏")
    size = struct.unpack('<Q', size_bytes)[0]
    data = mmap_obj.read(size)
    if len(data) != size:
        raise RuntimeError(f"解压时数据长度不足：期望 {size} 字节，仅读取到 {len(data)} 字节")
    return data


def decompress(compressed_path, output_path, lpaq8_path, save, gr_progress, max_workers):
    output_path = get_output_path(compressed_path, output_path)
    id_regex_tag = b"%id_regex%"
    id_tokens_tag = b"%id_tokens%"
    base_tag = b"%base_g_prime%"
    quality_tag = b"%quality%"
    eof_tag = b"%eof"

    block_count = 1
    with open(compressed_path, "r+b") as input_file:
        if os.path.getsize(compressed_path) == 0:
            return
        with mmap.mmap(input_file.fileno(), 0, access=mmap.ACCESS_READ) as mm:
            tqdm.write(f"info：开始解压 (安全长度模式 | 并行={max_workers})...")
            with multiprocessing.Pool(processes=max_workers, maxtasksperchild=1) as pool:
                pending = {}
                next_to_write = 1
                errors = []

                def flush_ready_results():
                    nonlocal next_to_write
                    while next_to_write in pending and pending[next_to_write].ready():
                        try:
                            id_block, g_prime, quality = pending[next_to_write].get()
                        except Exception as exc:
                            errors.append(exc)
                        else:
                            reconstruct_fastq(
                                output_path, id_block, g_prime, quality, is_first_block=(next_to_write == 1)
                            )
                        del pending[next_to_write]
                        next_to_write += 1

                while True:
                    remaining = mm.size() - mm.tell()
                    if remaining == len(eof_tag):
                        tail = mm.read(len(eof_tag))
                        if tail != eof_tag:
                            raise RuntimeError("解压结束时未发现EOF标记，压缩文件可能损坏")
                        tqdm.write("info：检测到EOF，解压结束。")
                        break
                    id_regex_data = read_chunk_safe(mm, id_regex_tag)
                    id_tokens_data = read_chunk_safe(mm, id_tokens_tag)
                    g_prime_data = read_chunk_safe(mm, base_tag)
                    quality_data = read_chunk_safe(mm, quality_tag)

                    tqdm.write(f"正在处理块: {block_count}")
                    pending[block_count] = pool.apply_async(
                        process_compressed_block,
                        (output_path, lpaq8_path, id_regex_data, id_tokens_data, g_prime_data, quality_data, save, block_count)
                    )
                    block_count += 1

                    if len(pending) > max_workers * 3:
                        time.sleep(0.1)
                    flush_ready_results()

                pool.close()
                pool.join()

                for idx in sorted(pending.keys()):
                    try:
                        id_block, g_prime, quality = pending[idx].get()
                    except Exception as exc:
                        errors.append(exc)
                        continue
                    reconstruct_fastq(output_path, id_block, g_prime, quality, is_first_block=(idx == 1))

                if errors:
                    raise RuntimeError(f"检测到 {len(errors)} 个解压块处理失败: {errors}")


def get_output_path(input_path, output_path):
    if input_path is None or not os.path.isfile(input_path):
        exit(1)
    if os.path.isdir(output_path):
        basename = os.path.splitext(os.path.basename(input_path))[0]
        return os.path.join(output_path, basename)
    return output_path


def delete_temp_files(output_path):
    temp_dir = os.path.dirname(output_path)
    for f in os.listdir(temp_dir):
        if f.startswith("temp_proc_") or f.startswith("temp_input") or f.startswith("temp_output") or f.startswith("temp_dec_"):
            try:
                os.remove(os.path.join(temp_dir, f))
            except Exception:
                pass
    for folder in ["back_compressed", "front_compressed", "temp_chunks"]:
        p = os.path.join(temp_dir, folder)
        if os.path.exists(p):
            try:
                shutil.rmtree(p)
            except Exception:
                pass


def main():
    lpaq8_path = f"{os.getcwd()}/lpaq8"
    parser = argparse.ArgumentParser(description='fastq lossy compress optimized (multithread)')
    parser.add_argument('--input_path', type=str, required=True, help='input_path')
    parser.add_argument('--output_path', type=str, required=True, help='output_path')
    parser.add_argument('--mode', type=str, required=True, help='compress(c) or decompress(d)')
    parser.add_argument('--save', type=str, default='False', help='save intermediate files (True/False)')
    parser.add_argument('--threads', type=int, default=os.cpu_count(), help='number of worker threads')
    parser.add_argument('--block_size', type=int, default=128 * 1024 * 1024, help='block size in bytes')
    args = parser.parse_args()
    save_flag = args.save.lower() == 'true'

    if args.mode in ['compress', 'c']:
        compress_multithread(args.input_path, args.output_path, lpaq8_path, save_flag, args.block_size, args.threads)
        if not save_flag:
            delete_temp_files(args.output_path)
    elif args.mode in ['decompress', 'd']:
        decompress(args.input_path, args.output_path, lpaq8_path, save_flag, None, args.threads)
        if not save_flag:
            delete_temp_files(args.output_path)
    else:
        print("Unknown mode")


if __name__ == '__main__':
    main()
