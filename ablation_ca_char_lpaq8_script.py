"""Ablation: chunk + dispatch + CA prediction on raw FASTQ chars + lpaq8.

- No integer mapping stage.
- Sequence alphabet uses raw chars (A/G/C/T/N...).
- Quality uses raw FASTQ ASCII chars.
- Reversible transform with marker byte 0x01 (outside FASTQ printable range).
"""

import os
import csv
import time
import mmap
import struct
import shutil
import tempfile
import subprocess
import multiprocessing
from pathlib import Path
from datetime import datetime
from collections import defaultdict

import psutil

# ===== 配置区 =====
INPUT_DIR = '/media/compress/新加卷1/New_Test'
OUTPUT_DIR = '/media/compress/新加卷/output/New_Test_LossLess/Ablation_CA_Char'
BLOCK_SIZE = 128 * 1024 * 1024
WORKERS = 4
LPAQ8_PATH = str((Path(__file__).resolve().parent / 'lpaq8').resolve())
COMPRESSION_LEVEL = '9'
# ===============

MARKER = 1
SEQ_CANDIDATES = [ord('A'), ord('G'), ord('C'), ord('T'), ord('N')]
QUAL_CANDIDATES = list(range(33, 127))


def get_file_size(path):
    return os.path.getsize(path)


def reads_per_block_from_fastq(path, block_size):
    with open(path, 'rb') as f:
        h = f.readline()
        s = f.readline().rstrip(b'\n\r')
        if not h or not s:
            return 1
        read_len = max(len(s), 1)
    bytes_per_read = read_len * 2
    return max(block_size // bytes_per_read, 1)


def iter_fastq_records_binary(path):
    with open(path, 'rb') as f:
        while True:
            h = f.readline()
            if not h:
                break
            s = f.readline()
            p = f.readline()
            q = f.readline()
            if not (s and p and q):
                raise RuntimeError('FASTQ format error: incomplete 4-line record')
            yield h.rstrip(b'\n\r'), s.rstrip(b'\n\r'), p.rstrip(b'\n\r'), q.rstrip(b'\n\r')


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


def ca_predict_decode(rows_prime, candidates):
    rules = defaultdict(int)
    n = len(rows_prime)
    m = len(rows_prime[0]) if n else 0
    out = [bytearray(m) for _ in range(n)]
    for i in range(n):
        for j in range(m):
            up = out[i - 1][j] if i else 0
            left = out[i][j - 1] if j else 0
            left_up = out[i - 1][j - 1] if i and j else 0
            cand = [(up, left_up, left, v) for v in candidates]
            top = max(cand, key=lambda r: rules[r])[3]
            center_prime = rows_prime[i][j]
            center = top if center_prime == MARKER else center_prime
            out[i][j] = center
            rules[(up, left_up, left, center)] += 1
    return [bytes(x) for x in out]


def encode_block(records):
    seq_rows = [r[1] for r in records]
    qual_rows = [r[3] for r in records]
    width = len(seq_rows[0]) if seq_rows else 0
    if any(len(r) != width for r in seq_rows + qual_rows):
        raise RuntimeError('Block contains variable read lengths; current script expects fixed length in a block')
    seq_prime = ca_predict_encode(seq_rows, SEQ_CANDIDATES)
    qual_prime = ca_predict_encode(qual_rows, QUAL_CANDIDATES)

    out = bytearray()
    out += struct.pack('<I', len(records))
    out += struct.pack('<I', width)
    for (h, _s, p, _q), sp, qp in zip(records, seq_prime, qual_prime):
        out += struct.pack('<I', len(h)) + h
        out += struct.pack('<I', len(p)) + p
        out += sp
        out += qp
    return bytes(out)


def decode_block(data):
    mv = memoryview(data)
    pos = 0
    n = struct.unpack_from('<I', mv, pos)[0]; pos += 4
    w = struct.unpack_from('<I', mv, pos)[0]; pos += 4
    headers, pluses, seqp, qualp = [], [], [], []
    for _ in range(n):
        hl = struct.unpack_from('<I', mv, pos)[0]; pos += 4
        h = bytes(mv[pos:pos+hl]); pos += hl
        pl = struct.unpack_from('<I', mv, pos)[0]; pos += 4
        p = bytes(mv[pos:pos+pl]); pos += pl
        sp = bytes(mv[pos:pos+w]); pos += w
        qp = bytes(mv[pos:pos+w]); pos += w
        headers.append(h); pluses.append(p); seqp.append(sp); qualp.append(qp)
    seq = ca_predict_decode(seqp, SEQ_CANDIDATES)
    qual = ca_predict_decode(qualp, QUAL_CANDIDATES)

    out = bytearray()
    for h, s, p, q in zip(headers, seq, pluses, qual):
        out += h + b'\n' + s + b'\n' + p + b'\n' + q + b'\n'
    return bytes(out)


def lpaq8_compress(inp, outp):
    p = subprocess.Popen([LPAQ8_PATH, COMPRESSION_LEVEL, inp, outp])
    p.wait()
    if p.returncode != 0:
        raise RuntimeError(f'lpaq8 compress failed: {inp}')


def lpaq8_decompress(inp, outp):
    p = subprocess.Popen([LPAQ8_PATH, 'd', inp, outp])
    p.wait()
    if p.returncode != 0:
        raise RuntimeError(f'lpaq8 decompress failed: {inp}')


def worker_compress(args):
    idx, records, out_dir = args
    part_path = os.path.join(out_dir, f'chunk_{idx}.part')
    with tempfile.TemporaryDirectory(prefix=f'ca_char_{idx}_') as td:
        raw = os.path.join(td, 'raw.bin')
        lz = os.path.join(td, 'raw.bin.lpaq8')
        payload = encode_block(records)
        with open(raw, 'wb') as f:
            f.write(payload)
        lpaq8_compress(raw, lz)
        data = open(lz, 'rb').read()
    with open(part_path, 'wb') as f:
        f.write(struct.pack('<Q', len(data)))
        f.write(data)
    return part_path


def compress_file(in_fastq, out_dir):
    base = os.path.splitext(os.path.basename(in_fastq))[0]
    out_path = os.path.join(out_dir, f'{base}.ca_char.lossless')
    temp_parts = os.path.join(out_dir, f'{base}_parts')
    os.makedirs(temp_parts, exist_ok=True)

    rpb = reads_per_block_from_fastq(in_fastq, BLOCK_SIZE)
    tasks, cur, idx = [], [], 1
    for rec in iter_fastq_records_binary(in_fastq):
        cur.append(rec)
        if len(cur) >= rpb:
            tasks.append((idx, cur, temp_parts))
            idx += 1
            cur = []
    if cur:
        tasks.append((idx, cur, temp_parts))

    start = time.time()
    with multiprocessing.Pool(processes=WORKERS) as pool:
        parts = pool.map(worker_compress, tasks)

    with open(out_path, 'wb') as out:
        out.write(struct.pack('<I', len(parts)))
        for i in range(1, len(parts)+1):
            pth = os.path.join(temp_parts, f'chunk_{i}.part')
            out.write(open(pth, 'rb').read())
    shutil.rmtree(temp_parts, ignore_errors=True)
    return out_path, time.time() - start, rpb


def decompress_file(in_path, out_fastq):
    with tempfile.TemporaryDirectory(prefix='ca_char_dec_') as td, open(in_path, 'rb') as f, open(out_fastq, 'wb') as out:
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        pos = 0
        n = struct.unpack_from('<I', mm, pos)[0]; pos += 4
        for i in range(1, n + 1):
            sz = struct.unpack_from('<Q', mm, pos)[0]; pos += 8
            comp = bytes(mm[pos:pos+sz]); pos += sz
            cin = os.path.join(td, f'c{i}.lpaq8')
            rout = os.path.join(td, f'r{i}.bin')
            open(cin, 'wb').write(comp)
            lpaq8_decompress(cin, rout)
            block = open(rout, 'rb').read()
            out.write(decode_block(block))
        mm.close()


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    csv_path = os.path.join(OUTPUT_DIR, f"ablation_ca_char_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv")
    header = ['file_name', 'input_size_mb', 'output_size_mb', 'compression_ratio', 'compression_time_s', 'compression_speed_mbs', 'reads_per_block', 'workers']
    with open(csv_path, 'w', newline='') as f:
        csv.writer(f).writerow(header)

    for fp in Path(INPUT_DIR).glob('*.fastq'):
        name = fp.stem
        print('Processing', fp)
        try:
            out_path, t, rpb = compress_file(str(fp), OUTPUT_DIR)
            ins = get_file_size(str(fp)); outs = get_file_size(out_path)
            ratio = ins / outs if outs else 0
            speed = (ins / 1024 / 1024) / t if t > 0 else 0

            # 完整还原校验
            restored = os.path.join(OUTPUT_DIR, f'{name}.restored.fastq')
            decompress_file(out_path, restored)
            ok = (open(fp, 'rb').read() == open(restored, 'rb').read())
            if not ok:
                raise RuntimeError(f'Round-trip mismatch for {fp.name}')
            os.remove(restored)

            row = [name, ins/1024/1024, outs/1024/1024, ratio, t, speed, rpb, WORKERS]
            with open(csv_path, 'a', newline='') as f:
                csv.writer(f).writerow(row)
            print(f'Done: {fp.name}, ratio={ratio:.3f}')
        except Exception as e:
            print('Error:', fp, e)


if __name__ == '__main__':
    main()
