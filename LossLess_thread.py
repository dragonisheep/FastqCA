import argparse
import os
import shutil
import sys
import time
import re
import gc
import mmap
import struct  # 【新增】用于打包二进制长度数据
from tqdm import tqdm
from PIL import Image, UnidentifiedImageError
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from collections import defaultdict
from itertools import product, chain
from Bio import SeqIO
import multiprocessing

# 确保目录下有 lpaq8.py
from lpaq8 import compress_file, decompress_file

Image.MAX_IMAGE_PIXELS = None

base_to_gray = {'A': 32, 'T': 64, 'G': 192, 'C': 224, 'N': 0}
gray_to_base = {32: 'A', 64: 'T', 192: 'G', 224: 'C', 0: 'N'}
base_gray_lut = np.zeros(256, dtype=np.uint8)
for base, value in base_to_gray.items():
    base_gray_lut[ord(base)] = value


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


def get_reads_num_per_block(fastq_path, block_size):
    with open(fastq_path, 'r') as file:
        try:
            first_record = next(SeqIO.parse(file, "fastq"))
            read_length = len(first_record.seq)
        except StopIteration:
            return 0, 0
    bytes_per_read = read_length * 2
    if bytes_per_read == 0: bytes_per_read = 1
    reads_per_block = block_size // bytes_per_read
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


def process_records(records_iter, record_count, read_length, rules_dict):
    id_block = []
    base_block = np.empty((record_count, read_length), dtype=np.uint8)
    quality_block = np.empty((record_count, read_length), dtype=np.uint8)

    row = 0
    for record in records_iter:
        id_str = record.description
        delimiters = find_delimiters(id_str)
        tokens = split_identifier(id_str, delimiters)
        regex = generate_regex(delimiters)
        id_block.append((tokens, regex))

        seq_bytes = bytes(record.seq)
        base_row = base_gray_lut[np.frombuffer(seq_bytes, dtype=np.uint8)]
        base_block[row].fill(0)
        base_block[row, :min(base_row.size, read_length)] = base_row[:read_length]

        qualities = record.letter_annotations["phred_quality"]
        quality_row = np.fromiter((q * 2 for q in qualities), dtype=np.uint8)
        quality_block[row].fill(0)
        quality_block[row, :min(quality_row.size, read_length)] = quality_row[:read_length]
        row += 1

    if row == 0:
        return None, None, None

    g_prime = generate_g_prime(base_block[:row], rules_dict)
    g_prime_img = Image.fromarray(g_prime.astype(np.uint8))
    quality_img = Image.fromarray(quality_block[:row].astype(np.uint8))
    return g_prime_img, quality_img, id_block


def save_intermediate_files(g_block, quality_block, id_block, output_path, block_count):
    front_dir = os.path.join(os.path.dirname(output_path), "front_compressed")
    os.makedirs(front_dir, exist_ok=True)
    g_block.save(os.path.join(front_dir, f'chunk_{block_count}_base.tiff'))
    quality_block.save(os.path.join(front_dir, f'chunk_{block_count}_quality.tiff'))
    with open(os.path.join(front_dir, f"chunk_{block_count}_id_tokens.txt"), 'w') as f1, \
            open(os.path.join(front_dir, f"chunk_{block_count}_id_regex.txt"), 'w') as f2:
        for tokens, regex in id_block:
            f1.write(' '.join(tokens) + '\n')
            f2.write(regex + '\n')


def compress_worker_subprocess(temp_input_path, temp_output_path, lpaq8_path):
    process = compress_file(temp_input_path, temp_output_path, lpaq8_path)
    if process: process.wait()


def write_safe_chunk(output_file, tag, data):
    """【核心修复】写入Tag + 8字节长度 + 数据，防止分隔符冲突"""
    # 写入 tag (如 %quality%)
    output_file.write(tag)
    # 写入 8字节 无符号整数表示长度 (Little Endian)
    # 这确保解压时知道精确读取多少字节，而不需要搜索下一个 tag
    output_file.write(struct.pack('<Q', len(data)))
    # 写入数据
    output_file.write(data)


def back_compress_worker(g_block, quality_block, id_block, lpaq8_path, output_path, save, block_count):
    part_output_path = os.path.join(os.path.dirname(output_path), f"chunk_{block_count}.part")
    temp_prefix = os.path.join(os.path.dirname(output_path), f"temp_proc_{block_count}")
    temp_input_path, temp_output_path = f"{temp_prefix}_input", f"{temp_prefix}_output"

    back_dir = os.path.join(os.path.dirname(output_path), "back_compressed")
    if save: os.makedirs(back_dir, exist_ok=True)

    try:
        with open(part_output_path, "wb") as output_file:
            # 1. ID Regex
            with open(temp_input_path, "w") as f:
                for _, regex in id_block: f.write(regex + '\n')
            compress_worker_subprocess(temp_input_path, temp_output_path, lpaq8_path)
            if not os.path.exists(temp_output_path): raise RuntimeError("LPAQ8 compression failed for ID Regex")
            with open(temp_output_path, "rb") as f:
                data = f.read()
            # 使用安全写入
            write_safe_chunk(output_file, b"%id_regex%", data)
            if save:
                with open(os.path.join(back_dir, f"chunk_{block_count}_id_regex.lpaq8"), "wb") as sf: sf.write(data)

            # 2. ID Tokens
            with open(temp_input_path, "w") as f:
                for tokens, _ in id_block: f.write(' '.join(tokens) + '\n')
            compress_worker_subprocess(temp_input_path, temp_output_path, lpaq8_path)
            if not os.path.exists(temp_output_path): raise RuntimeError("LPAQ8 compression failed for ID Tokens")
            with open(temp_output_path, "rb") as f:
                data = f.read()
            write_safe_chunk(output_file, b"%id_tokens%", data)
            if save:
                with open(os.path.join(back_dir, f"chunk_{block_count}_id_tokens.lpaq8"), "wb") as sf: sf.write(data)

            # 3. Base
            g_block.save(temp_input_path, format="tiff")
            compress_worker_subprocess(temp_input_path, temp_output_path, lpaq8_path)
            if not os.path.exists(temp_output_path): raise RuntimeError("LPAQ8 compression failed for Base")
            with open(temp_output_path, "rb") as f:
                data = f.read()
            write_safe_chunk(output_file, b"%base_g_prime%", data)
            if save:
                with open(os.path.join(back_dir, f"chunk_{block_count}_base_g_prime.lpaq8"), "wb") as sf: sf.write(data)

            # 4. Quality
            quality_block.save(temp_input_path, format="tiff")
            compress_worker_subprocess(temp_input_path, temp_output_path, lpaq8_path)
            if not os.path.exists(temp_output_path): raise RuntimeError("LPAQ8 compression failed for Quality")
            with open(temp_output_path, "rb") as f:
                data = f.read()
            write_safe_chunk(output_file, b"%quality%", data)
            if save:
                with open(os.path.join(back_dir, f"chunk_{block_count}_quality.lpaq8"), "wb") as sf: sf.write(data)

    finally:
        if os.path.exists(temp_input_path): os.remove(temp_input_path)
        if os.path.exists(temp_output_path): os.remove(temp_output_path)
    return block_count


def process_block_task_from_file(temp_chunk_path, block_count, output_path, lpaq8_path, save, record_count):
    try:
        gc.collect()
        with open(temp_chunk_path, 'r') as f:
            parser = SeqIO.parse(f, "fastq")
            try:
                first_record = next(parser)
            except StopIteration:
                return block_count

            read_length = len(first_record.seq)
            rules_dict = init_rules_dict()

            # 预分配数组并逐行填充，避免 list -> numpy 的中间副本
            g_block, quality_block, id_block = process_records(
                chain([first_record], parser), record_count, read_length, rules_dict
            )
            del rules_dict
            gc.collect()

            if g_block is None:
                return block_count
            if save:
                save_intermediate_files(g_block, quality_block, id_block, output_path, block_count)
            back_compress_worker(g_block, quality_block, id_block, lpaq8_path, output_path, save, block_count)
    finally:
        if os.path.exists(temp_chunk_path):
            try:
                os.remove(temp_chunk_path)
            except:
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
    if out_dir and not os.path.exists(out_dir): os.makedirs(out_dir)

    reads_per_block, total_reads = get_reads_num_per_block(fastq_path, block_size)
    total_block_count = int(np.ceil(os.path.getsize(fastq_path) / block_size)) if block_size > 0 else 1

    tqdm.write(f"info：多进程模式 (进程数={max_workers}) | 模式: 安全长度头封装")
    temp_chunk_dir = os.path.join(out_dir, "temp_chunks")
    os.makedirs(temp_chunk_dir, exist_ok=True)

    read_count_per_block, block_count = 0, 1

    # 保持 maxtasksperchild=1 防止内存泄漏
    with multiprocessing.Pool(processes=max_workers, maxtasksperchild=1) as pool:
        results = []
        errors = []
        tqdm.write(f"info：正在读取 FASTQ 并分发任务...")

        temp_chunk_path = os.path.join(temp_chunk_dir, f"chunk_src_{block_count}.fastq")
        temp_handle = open(temp_chunk_path, 'w')

        def dispatch_current_chunk(path, count, record_total):
            temp_handle.close()
            res = pool.apply_async(
                process_block_task_from_file,
                (path, count, output_path, lpaq8_path, save, record_total)
            )
            results.append(res)

        try:
            with open(fastq_path, 'r') as file:
                for record in tqdm(SeqIO.parse(file, "fastq"), desc="Chunking", total=total_reads, unit="reads"):
                    SeqIO.write([record], temp_handle, "fastq")
                    read_count_per_block += 1
                    if read_count_per_block >= reads_per_block:
                        dispatch_current_chunk(temp_chunk_path, block_count, read_count_per_block)
                        if len(results) > max_workers * 2:
                            results = [r for r in results if not r.ready()]
                            if len(results) > max_workers * 3: time.sleep(1)
                        block_count += 1
                        temp_chunk_path = os.path.join(temp_chunk_dir, f"chunk_src_{block_count}.fastq")
                        temp_handle = open(temp_chunk_path, 'w')
                        read_count_per_block = 0
                if read_count_per_block > 0:
                    dispatch_current_chunk(temp_chunk_path, block_count, read_count_per_block)
                else:
                    temp_handle.close()
                    block_count -= 1
        finally:
            try:
                temp_handle.close()
            except Exception:
                pass

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
    except:
        pass


# --- Decompression Logic ---

def monitor(process, temp_input_path, temp_output_path):
    while True:
        if os.path.exists(temp_input_path) and os.path.exists(temp_output_path):
            if process.poll() is not None: break
        time.sleep(0.1)


def decompress_with_monitor(temp_input_path, temp_output_path, lpaq8_path):
    process = decompress_file(temp_input_path, temp_output_path, lpaq8_path)
    # 直接等待子进程结束，避免 monitor 轮询带来的额外阻塞
    process.wait()


def process_compressed_block(output_path, lpaq8_path, id_regex_path, id_tokens_path, g_prime_path, quality_path, save,
                             block_count):
    back_compress_dir = os.path.join(os.path.dirname(output_path), "back_compressed")
    front_compress_dir = os.path.join(os.path.dirname(output_path), "front_compressed")
    if save:
        os.makedirs(back_compress_dir, exist_ok=True)
        os.makedirs(front_compress_dir, exist_ok=True)

    temp_prefix = os.path.join(os.path.dirname(output_path), f"temp_dec_{block_count}")
    temp_input_path, temp_output_path = f"{temp_prefix}_in", f"{temp_prefix}_out"
    id_regex, id_tokens, g_prime, quality = None, None, None, None
    id_block = None
    temp_fastq_path = os.path.join(os.path.dirname(output_path), f"temp_dec_fastq_{block_count}.fastq")

    try:
        # 1. ID Regex
        shutil.copyfile(id_regex_path, temp_input_path)
        with open(temp_output_path, "w+") as f:
            decompress_with_monitor(temp_input_path, temp_output_path, lpaq8_path)
            f.seek(0)
            id_regex = [line.strip() for line in f.readlines()]
        if save:
            shutil.copy(temp_input_path, os.path.join(back_compress_dir, f"chunk_{block_count}_id_regex.lpaq8"))
            shutil.copy(temp_output_path, os.path.join(front_compress_dir, f"chunk_{block_count}_id_regex.txt"))

        # 2. ID Tokens
        shutil.copyfile(id_tokens_path, temp_input_path)
        with open(temp_output_path, "w+") as f:
            decompress_with_monitor(temp_input_path, temp_output_path, lpaq8_path)
            f.seek(0)
            id_tokens = [line.strip() for line in f.readlines()]
        if save:
            shutil.copy(temp_input_path, os.path.join(back_compress_dir, f"chunk_{block_count}_id_tokens.lpaq8"))
            shutil.copy(temp_output_path, os.path.join(front_compress_dir, f"chunk_{block_count}_id_tokens.txt"))

        id_block = zip(id_tokens, id_regex)

        # 3. Quality
        shutil.copyfile(quality_path, temp_input_path)
        decompress_with_monitor(temp_input_path, temp_output_path, lpaq8_path)
        with Image.open(temp_output_path) as img:
            quality = img.copy()
        if save:
            shutil.copy(temp_input_path, os.path.join(back_compress_dir, f'chunk_{block_count}_quality.lpaq8'))
            shutil.copy(temp_output_path, os.path.join(front_compress_dir, f'chunk_{block_count}_quality.tiff'))

        # 4. G Prime
        shutil.copyfile(g_prime_path, temp_input_path)
        decompress_with_monitor(temp_input_path, temp_output_path, lpaq8_path)
        with Image.open(temp_output_path) as img:
            g_prime = img.copy()
        if save:
            shutil.copy(temp_input_path, os.path.join(back_compress_dir, f'chunk_{block_count}_base_g_prime.lpaq8'))
            shutil.copy(temp_output_path, os.path.join(front_compress_dir, f'chunk_{block_count}_base_g_prime.tiff'))

    finally:
        if os.path.exists(temp_input_path): os.remove(temp_input_path)
        if os.path.exists(temp_output_path): os.remove(temp_output_path)
        # 清理分块缓存文件，避免磁盘和内存堆积
        for chunk_path in (id_regex_path, id_tokens_path, g_prime_path, quality_path):
            try:
                if os.path.exists(chunk_path):
                    os.remove(chunk_path)
            except Exception:
                pass

        if id_block is not None and g_prime is not None and quality is not None:
            reconstruct_fastq(output_path, id_block, g_prime, quality, is_first_block=True, custom_path=temp_fastq_path)
        # 及时释放块内大对象，避免在主进程堆积
        id_regex = id_tokens = id_block = g_prime = quality = None
        gc.collect()

    return temp_fastq_path


def reconstruct_id(tokens, regex):
    reconstructed_ids = []
    for t, r in zip(tokens, regex):
        id_str = r
        token_list = t.split()
        for i, token in enumerate(token_list):
            id_str = id_str.replace(f"{{{{T{i + 1}}}}}", token)
        reconstructed_ids.append(id_str)
    return reconstructed_ids


def reconstruct_g_from_g_prime(g_prime_array, rules_dict):
    De_g = np.zeros_like(g_prime_array)
    rows, cols = g_prime_array.shape
    for i in range(rows):
        for j in range(cols):
            center = g_prime_array[i, j]
            up = De_g[i - 1, j] if i != 0 else 0
            left = De_g[i, j - 1] if j != 0 else 0
            left_up = De_g[i - 1, j - 1] if i != 0 and j != 0 else 0
            candidates = [(up, left_up, left, v) for v in [32, 224, 192, 64, 0]]
            top_rule = max(candidates, key=lambda r: rules_dict[r])
            De_g[i, j] = top_rule[3] if g_prime_array[i, j] == 1 else center
            matched_rule = (up, left_up, left, De_g[i, j])
            rules_dict[matched_rule] += 1
    return De_g


def reconstruct_base_and_quality(g_prime_img, quality_img):
    bases, qualities = [], []
    rules_dict = defaultdict(int)
    g_prime_array = np.array(g_prime_img)
    bases_array = reconstruct_g_from_g_prime(g_prime_array, rules_dict)
    quality_array = np.array(quality_img)
    for i in range(bases_array.shape[0]):
        base_str = ''.join([gray_to_base.get(pixel, 'N') for pixel in bases_array[i]])
        quality_scores = [q / 2 for q in quality_array[i]]
        bases.append(base_str)
        qualities.append(quality_scores)
    return bases, qualities


def reconstruct_fastq(output_path, id_block, g_prime_img, quality_img, is_first_block=False, custom_path=None):
    """流式写入 FASTQ，避免构造大列表/SeqRecord 造成内存峰值。"""
    target_fastq_path = custom_path if custom_path else (
        output_path if output_path.endswith('.fastq') else os.path.splitext(output_path)[0] + '.fastq'
    )

    # 预计算 LUT，减少循环开销
    gray_to_char = {v: k for k, v in base_to_gray.items()}
    base_lut = np.frombuffer(
        bytes([ord(gray_to_char.get(i, 'N')) for i in range(256)]),
        dtype=np.uint8
    )

    id_block = list(id_block)
    id_tokens = [item[0] for item in id_block]
    id_regex = [item[1] for item in id_block]
    ids = reconstruct_id(id_tokens, id_regex)

    g_prime_array = np.array(g_prime_img, dtype=np.uint8)
    quality_array = np.array(quality_img, dtype=np.uint8)
    rows, cols = g_prime_array.shape

    # 行递归重建，保持常数行缓冲，避免整块字符串常驻
    rules_dict = defaultdict(int)
    prev_row = np.zeros(cols, dtype=np.uint8)
    mode = 'w' if custom_path or is_first_block else 'a'
    with open(target_fastq_path, mode, buffering=1024 * 1024) as output_handle:
        for i in range(rows):
            de_row = np.empty(cols, dtype=np.uint8)
            gp_row = g_prime_array[i]
            q_row = quality_array[i]
            for j in range(cols):
                center = gp_row[j]
                up = prev_row[j] if i != 0 else 0
                left = de_row[j - 1] if j != 0 else 0
                left_up = prev_row[j - 1] if i != 0 and j != 0 else 0
                candidates = [(up, left_up, left, v) for v in [32, 224, 192, 64, 0]]
                top_rule = max(candidates, key=lambda r: rules_dict[r])
                de_row[j] = top_rule[3] if gp_row[j] == 1 else center
                matched_rule = (up, left_up, left, de_row[j])
                rules_dict[matched_rule] += 1
            seq_bytes = base_lut[de_row].tobytes()
            seq_str = seq_bytes.decode('ascii')
            qual_bytes = (q_row // 2 + 33).astype(np.uint8).tobytes()
            qual_str = qual_bytes.decode('ascii')
            output_handle.write(f"@{ids[i]}\n{seq_str}\n+\n{qual_str}\n")
            prev_row = de_row

    # 显式释放大数组
    del g_prime_array, quality_array, prev_row, base_lut, ids
    gc.collect()
    return target_fastq_path


def read_chunk_safe(mmap_obj, tag):
    """严格按当前位置读取 tag + 长度 + 数据，避免搜索误匹配导致偏移"""
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
        raise RuntimeError(
            f"解压时数据长度不足：期望 {size} 字节，仅读取到 {len(data)} 字节")
    return data


def decompress(compressed_path, output_path, lpaq8_path, save, gr_progress, max_workers):
    output_path = get_output_path(compressed_path, output_path)
    # 标记符定义
    id_regex_tag = b"%id_regex%"
    id_tokens_tag = b"%id_tokens%"
    base_tag = b"%base_g_prime%"
    quality_tag = b"%quality%"
    eof_tag = b"%eof"

    block_count = 1
    final_fastq_path = output_path if output_path.endswith('.fastq') else os.path.splitext(output_path)[0] + '.fastq'
    with open(compressed_path, "r+b") as input_file:
        if os.path.getsize(compressed_path) == 0: return

        with mmap.mmap(input_file.fileno(), 0, access=mmap.ACCESS_READ) as mm:
            tqdm.write(f"info：开始解压 (安全长度模式 | 并行={max_workers})...")
            # 使用进程池保证每个块在独立进程中处理，内存可在子进程退出时被OS回收；参数传递仅用磁盘路径避免巨型对象序列化
            # 大块（700MB+）解压会瞬时产生 GB 级中间对象，硬限制并发为最多 2，防止主机内存暴涨
            pool_workers = max(1, min(max_workers, 2))
            with multiprocessing.Pool(processes=pool_workers, maxtasksperchild=1) as pool:
                pending = {}
                next_to_write = 1
                errors = []
                chunk_dir = os.path.join(os.path.dirname(output_path), "temp_dec_chunks")
                os.makedirs(chunk_dir, exist_ok=True)

                def flush_ready_results():
                    nonlocal next_to_write
                    while next_to_write in pending and pending[next_to_write].ready():
                        try:
                            temp_fastq_path = pending[next_to_write].get()
                        except Exception as exc:
                            errors.append(exc)
                        else:
                            mode = 'w' if next_to_write == 1 else 'a'
                            with open(temp_fastq_path, 'r') as src, open(final_fastq_path, mode) as dst:
                                shutil.copyfileobj(src, dst)
                            os.remove(temp_fastq_path)
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

                    # 将大块数据落盘，仅传递路径给子进程，避免 pickling 大对象
                    def dump_chunk(tag, data):
                        path = os.path.join(chunk_dir, f"chunk_{block_count}_{tag}.bin")
                        with open(path, "wb") as fh:
                            fh.write(data)
                        return path

                    id_regex_path = dump_chunk("id_regex", id_regex_data)
                    id_tokens_path = dump_chunk("id_tokens", id_tokens_data)
                    g_prime_path = dump_chunk("g_prime", g_prime_data)
                    quality_path = dump_chunk("quality", quality_data)

                    tqdm.write(f"正在处理块: {block_count}")

                    pending[block_count] = pool.apply_async(
                        process_compressed_block,
                        (output_path, lpaq8_path, id_regex_path, id_tokens_path, g_prime_path, quality_path, save,
                         block_count)
                    )
                    block_count += 1

                    if len(pending) > max_workers * 3:
                        time.sleep(0.1)
                    flush_ready_results()
                    # 防止主进程一次性将所有分块排队，导致输入队列累积大块二进制数据占用内存。
                    # 将排队上限收紧到“等于并行度”，避免 4 线程下出现 8~12 个大块同时驻留内存。
                    while len(pending) >= max_workers:
                        flush_ready_results()
                        if len(pending) >= max_workers:
                            time.sleep(0.05)
                    # 主进程不再持有块数据，立刻触发GC降低占用
                    del id_regex_data, id_tokens_data, g_prime_data, quality_data
                    gc.collect()

                pool.close()
                pool.join()

                for idx in sorted(pending.keys()):
                    try:
                        temp_fastq_path = pending[idx].get()
                    except Exception as exc:
                        errors.append(exc)
                        continue
                    mode = 'w' if idx == 1 else 'a'
                    with open(temp_fastq_path, 'r') as src, open(final_fastq_path, mode) as dst:
                        shutil.copyfileobj(src, dst)
                    os.remove(temp_fastq_path)

                if errors:
                    raise RuntimeError(f"检测到 {len(errors)} 个解压块处理失败: {errors}")
                try:
                    shutil.rmtree(chunk_dir)
                except Exception:
                    pass


def get_output_path(input_path, output_path):
    if input_path is None or not os.path.isfile(input_path): exit(1)
    if os.path.isdir(output_path):
        basename = os.path.splitext(os.path.basename(input_path))[0]
        return os.path.join(output_path, basename)
    return output_path


def delete_temp_files(output_path):
    temp_dir = os.path.dirname(output_path)
    for f in os.listdir(temp_dir):
        if f.startswith("temp_proc_") or f.startswith("temp_input") or f.startswith("temp_output") or f.startswith(
                "temp_dec_"):
            try:
                os.remove(os.path.join(temp_dir, f))
            except:
                pass
    for folder in ["back_compressed", "front_compressed", "temp_chunks", "temp_dec_chunks"]:
        p = os.path.join(temp_dir, folder)
        if os.path.exists(p):
            try:
                shutil.rmtree(p)
            except:
                pass


def main():
    lpaq8_path = f"{os.getcwd()}/lpaq8"
    parser = argparse.ArgumentParser(description='fastq compress optimized')
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
        if not save_flag: delete_temp_files(args.output_path)
    elif args.mode in ['decompress', 'd']:
        decompress(args.input_path, args.output_path, lpaq8_path, save_flag, None, args.threads)
        if not save_flag: delete_temp_files(args.output_path)
    else:
        print("Unknown mode")


if __name__ == '__main__':
    main()
