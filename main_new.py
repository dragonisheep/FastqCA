import argparse
import os
from typing import Optional

from Bio import SeqIO

from Lossy_thread import (
    compress_multithread as lossy_compress,
    decompress as lossy_decompress,
    delete_temp_files as lossy_cleanup,
    get_output_path as lossy_output_path,
)
from LossLess_thread import (
    compress_multithread as lossless_compress,
    decompress as lossless_decompress,
    delete_temp_files as lossless_cleanup,
)


LOSSY_COMMANDS = ["Lossy", "lossy"]
LOSSLESS_COMMANDS = ["LossLess", "lossless", "lossLess", "Lossless"]


def count_reads(fastq_path: str) -> int:
    return sum(1 for _ in SeqIO.parse(fastq_path, "fastq"))


def read_manifest(manifest_path: str) -> Optional[int]:
    if not os.path.exists(manifest_path):
        return None
    with open(manifest_path, "r") as handle:
        content = handle.read().strip()
        try:
            return int(content)
        except ValueError:
            return None


def write_manifest(manifest_path: str, read_count: int) -> None:
    with open(manifest_path, "w") as handle:
        handle.write(str(read_count))


def main() -> None:
    lpaq8_path = f"{os.getcwd()}/lpaq8"

    parser = argparse.ArgumentParser(description="fastq compress (multithread version)")
    parser.add_argument("--compressor", type=str, default="Lossy", help="Lossy or LossLess?")
    parser.add_argument("--input_path", type=str, required=True, help="input_path")
    parser.add_argument("--output_path", type=str, required=True, help="output_path")
    parser.add_argument("--mode", type=str, required=True, help="compress(c) or decompress(d)")
    parser.add_argument("--save", type=str, default="False", help="save intermediate files (True/False)")
    parser.add_argument("--threads", type=int, default=os.cpu_count(), help="number of worker threads")
    parser.add_argument("--block_size", type=int, default=128 * 1024 * 1024, help="block size in bytes")
    parser.add_argument(
        "--manifest",
        type=str,
        default=None,
        help="可选的 read 数量清单文件路径。默认为压缩文件后缀 .readcount",
    )
    args = parser.parse_args()

    save_flag = args.save.lower() == "true"

    if args.compressor in LOSSY_COMMANDS:
        archive_path = lossy_output_path(args.input_path, args.output_path)
        if args.mode in ["compress", "c"]:
            manifest_path = args.manifest or f"{archive_path}.readcount"
            read_count = count_reads(args.input_path)
            lossy_compress(
                args.input_path,
                args.output_path,
                lpaq8_path,
                save_flag,
                args.block_size,
                args.threads,
            )
            write_manifest(manifest_path, read_count)
            if not save_flag:
                lossy_cleanup(args.output_path)
        elif args.mode in ["decompress", "d"]:
            manifest_path = args.manifest or f"{args.input_path}.readcount"
            lossy_decompress(
                args.input_path,
                args.output_path,
                lpaq8_path,
                save_flag,
                None,
                args.threads,
            )
            restored_base = lossy_output_path(args.input_path, args.output_path)
            restored_path = (
                restored_base if restored_base.endswith(".fastq") else f"{restored_base}.fastq"
            )
            restored_reads = count_reads(restored_path)
            expected_reads = read_manifest(manifest_path)
            if expected_reads is not None and restored_reads != expected_reads:
                raise RuntimeError(
                    f"解压后的 read 数量 ({restored_reads}) 与清单记录的数量 ({expected_reads}) 不一致"
                )
            if not save_flag:
                lossy_cleanup(args.output_path)
        else:
            raise SystemExit("错误：没有指定正确的模式")

    elif args.compressor in LOSSLESS_COMMANDS:
        if args.mode in ["compress", "c"]:
            lossless_compress(
                args.input_path,
                args.output_path,
                lpaq8_path,
                save_flag,
                args.block_size,
                args.threads,
            )
            if not save_flag:
                lossless_cleanup(args.output_path)
        elif args.mode in ["decompress", "d"]:
            lossless_decompress(
                args.input_path,
                args.output_path,
                lpaq8_path,
                save_flag,
                None,
                args.threads,
            )
            if not save_flag:
                lossless_cleanup(args.output_path)
        else:
            raise SystemExit("错误：没有指定正确的模式")
    else:
        raise SystemExit("错误：没有指定正确的压缩器")


if __name__ == "__main__":
    main()
