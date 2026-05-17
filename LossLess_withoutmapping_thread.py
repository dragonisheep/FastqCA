import argparse
import gc
import mmap
import os
import shutil
import struct
import tempfile
import time
from collections import defaultdict
from pathlib import Path

import multiprocessing
from Bio import SeqIO


MARKER = 1
SEQ_CANDIDATES = [ord('A'), ord('G'), ord('C'), ord('T'), ord('N')]
QUAL_CANDIDATES = list(range(33, 127))


def get_output_path(input_path, output_path):
    file_name = Path(input_path).stem
    return os.path.join(output_path, f"{file_name}.withoutmapping.lossless")


def get_reads_num_per_block(fastq_path, block_size):
    with open(fastq_path, 'r') as file:
        try:
            first_record = next(SeqIO.parse(file, 'fastq'))
            read_length = len(first_record.seq)
        except StopIteration:
            return 1, 0
    bytes_per_read = max(read_length * 2, 1)
    reads_per_block = max(block_size // bytes_per_read, 1)
    total_reads = max(os.path.getsize(fastq_path) // bytes_per_read, 1)
    return reads_per_block, total_reads


def ca_predict_encode(rows, candidates):
    rules = defaultdict(int)
    out = [bytearray(r) for r in rows]
    n = len(rows)
    m = len(rows[0]) if n else 0
    for i in range(n):
        for j in range(m):
            center = rows[i][j]
            up = rows[i - 1][j] if i else 0
            left = rows[i][j - 1] if j else 0
            left_up = rows[i - 1][j - 1] if i and j else 0
            cand = [(up, left_up, left, v) for v in candidates]
            top = max(cand, key=lambda r: rules[r])[3]
            out[i][j] = MARKER if top == center else center
            rules[(up, left_up, left, center)] += 1
    return [bytes(x) for x in out]


def encode_records(records):
    seq_rows = [str(r.seq).encode('ascii') for r in records]
    qual_rows = [r.letter_annotations['phred_quality'] for r in records]
    qual_rows = [bytes(q + 33 for q in row) for row in qual_rows]

    width = len(seq_rows[0]) if seq_rows else 0
    if any(len(r) != width for r in seq_rows + qual_rows):
        raise RuntimeError('Block contains variable read lengths; fixed length required inside each block')

    seq_prime = ca_predict_encode(seq_rows, SEQ_CANDIDATES)
    qual_prime = ca_predict_encode(qual_rows, QUAL_CANDIDATES)

    out = bytearray()
    out += struct.pack('<I', len(records))
    out += struct.pack('<I', width)
    for rec, sp, qp in zip(records, seq_prime, qual_prime):
        header = rec.description.encode('utf-8')
        out += struct.pack('<I', len(header)) + header
        out += sp
        out += qp
    return bytes(out)


def lpaq8_compress(inp, outp, lpaq8_path):
    import subprocess

    p = subprocess.Popen([lpaq8_path, '9', inp, outp], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    p.wait()
    if p.returncode != 0:
        raise RuntimeError(f'lpaq8 compression failed: {inp}')


def process_block_task_from_file(args):
    temp_chunk_path, block_count, part_dir, lpaq8_path = args
    part_output_path = os.path.join(part_dir, f'chunk_{block_count}.part')

    try:
        with open(temp_chunk_path, 'r') as f:
            records = list(SeqIO.parse(f, 'fastq'))
        if not records:
            return block_count

        payload = encode_records(records)
        with tempfile.TemporaryDirectory(prefix=f'withoutmapping_{block_count}_') as td:
            raw = os.path.join(td, 'raw.bin')
            out = os.path.join(td, 'raw.bin.lpaq8')
            with open(raw, 'wb') as w:
                w.write(payload)
            lpaq8_compress(raw, out, lpaq8_path)
            data = Path(out).read_bytes()

        with open(part_output_path, 'wb') as fw:
            fw.write(struct.pack('<Q', len(data)))
            fw.write(data)
    finally:
        if os.path.exists(temp_chunk_path):
            os.remove(temp_chunk_path)
        gc.collect()

    return block_count


def merge_parts(output_path, total_blocks, part_dir):
    with open(output_path, 'wb') as final_file:
        final_file.write(struct.pack('<I', total_blocks))
        for i in range(1, total_blocks + 1):
            part_path = os.path.join(part_dir, f'chunk_{i}.part')
            if not os.path.exists(part_path):
                raise FileNotFoundError(f'missing part: {part_path}')
            final_file.write(Path(part_path).read_bytes())
            os.remove(part_path)


def compress_multithread(fastq_path, output_path, lpaq8_path, save, block_size, max_workers):
    output_path = get_output_path(fastq_path, output_path)
    out_dir = os.path.dirname(output_path)
    os.makedirs(out_dir, exist_ok=True)

    reads_per_block, _ = get_reads_num_per_block(fastq_path, block_size)

    temp_chunk_dir = os.path.join(out_dir, 'temp_chunks_withoutmapping')
    part_dir = os.path.join(out_dir, 'temp_parts_withoutmapping')
    os.makedirs(temp_chunk_dir, exist_ok=True)
    os.makedirs(part_dir, exist_ok=True)

    read_count_per_block, block_count = 0, 1
    pending = []
    max_inflight = max(1, max_workers * 2)

    with multiprocessing.Pool(processes=max_workers, maxtasksperchild=1) as pool:
        temp_chunk_path = os.path.join(temp_chunk_dir, f'chunk_src_{block_count}.fastq')
        temp_handle = open(temp_chunk_path, 'w')

        def dispatch_current_chunk(path, count):
            temp_handle.close()
            return pool.apply_async(process_block_task_from_file, ((path, count, part_dir, lpaq8_path),))

        try:
            with open(fastq_path, 'r') as file:
                for record in SeqIO.parse(file, 'fastq'):
                    SeqIO.write([record], temp_handle, 'fastq')
                    read_count_per_block += 1
                    if read_count_per_block >= reads_per_block:
                        pending.append(dispatch_current_chunk(temp_chunk_path, block_count))
                        if len(pending) >= max_inflight:
                            pending.pop(0).get()
                        block_count += 1
                        temp_chunk_path = os.path.join(temp_chunk_dir, f'chunk_src_{block_count}.fastq')
                        temp_handle = open(temp_chunk_path, 'w')
                        read_count_per_block = 0

                if read_count_per_block > 0:
                    pending.append(dispatch_current_chunk(temp_chunk_path, block_count))
                else:
                    temp_handle.close()
                    block_count -= 1
        finally:
            try:
                temp_handle.close()
            except Exception:
                pass

        for job in pending:
            job.get()

    merge_parts(output_path, block_count, part_dir)
    if not save:
        shutil.rmtree(temp_chunk_dir, ignore_errors=True)
        shutil.rmtree(part_dir, ignore_errors=True)


def delete_temp_files(output_path):
    for folder in ('temp_chunks_withoutmapping', 'temp_parts_withoutmapping'):
        shutil.rmtree(os.path.join(output_path, folder), ignore_errors=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='LossLess without mapping backend ablation')
    parser.add_argument('--input_path', type=str, required=True)
    parser.add_argument('--output_path', type=str, required=True)
    parser.add_argument('--lpaq8_path', type=str, default=f"{os.getcwd()}/lpaq8")
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--block_size', type=int, default=128 * 1024 * 1024)
    parser.add_argument('--save', type=str, default='False')
    args = parser.parse_args()

    compress_multithread(
        args.input_path,
        args.output_path,
        args.lpaq8_path,
        args.save.lower() == 'true',
        args.block_size,
        args.threads,
    )
