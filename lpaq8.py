import os
import subprocess
import threading
import time


def compress_lpaq8(lpaq8_path, input_path, output_path, compression_level='9'):
    """
    使用lpaq8压缩单个文件。

    参数:
    - lpaq8_path: lpaq8的路径。
    - compression_level: 压缩级别。
    - input_path: 输入文件路径。
    - output_path: 输出文件路径。
    """
    command = [lpaq8_path, compression_level, input_path, output_path]
    try:
        process = subprocess.Popen(command)
        return process
        # print(f"文件 {input_path} 压缩成功，保存为 {output_path}")
    except subprocess.CalledProcessError as e:
        print(f"压缩过程中出错: {e}")
    except Exception as e:
        print(f"发生未知错误: {str(e)}")

def compress_lpaq8_test(lpaq8_path, input_stream, output_path, compression_level='9'):
    """
    使用lpaq8压缩单个文件。

    参数:
    - lpaq8_path: lpaq8的路径。
    - compression_level: 压缩级别。
    - input_path: 输入文件路径。
    - output_path: 输出文件路径。
    """
    command = [lpaq8_path, compression_level, input_stream, output_path]
    try:
        subprocess.run(command, check=True)
        print(f"文件 {input_stream} 压缩成功，保存为 {output_path}")
    except subprocess.CalledProcessError as e:
        print(f"压缩过程中出错: {e}")
    except Exception as e:
        print(f"发生未知错误: {str(e)}")


def decompress_lpaq8(lpaq8_path, input_path, output_path):
    """
    使用lpaq8压缩单个文件。

    参数:
    - lpaq8_path: lpaq8的路径。
    - compression_level: 压缩级别。
    - input_path: 输入文件路径。
    - output_path: 输出文件路径。
    """
    command = [lpaq8_path, 'd', input_path, output_path]
    try:
        process = subprocess.Popen(command)
        # print(f"文件 {input_path} 解压成功，保存为 {output_path}")
        return process
    except subprocess.CalledProcessError as e:
        print(f"解压过程中出错: {e}")
    except Exception as e:
        print(f"发生未知错误: {str(e)}")


def compress_file(input_file, output_file, lpaq8_path, compression_level='9'):
    """
    压缩指定的文件。

    参数:
    - input_file: 完整的输入文件路径。
    - output_file: 压缩文件的输出目录。
    - lpaq8_path: lpaq8压缩器的完整路径。
    - compression_level: 压缩级别（默认为9，范围0-9）。
    """
    # 调用lpaq8进行压缩
    return compress_lpaq8(lpaq8_path, input_file, output_file, compression_level)


def decompress_file(input_file, output_file, lpaq8_path):
    """
    压缩指定的文件。

    参数:
    - input_file: 完整的输入文件路径。
    - output_directory: 压缩文件的输出目录。
    - lpaq8_path: lpaq8压缩器的完整路径。
    """
    return decompress_lpaq8(lpaq8_path, input_file, output_file)


def compress_all_files_in_directory(input_directory, output_directory, lpaq8_path, compression_level='9'):
    """
    压缩目录中的所有文件。

    参数:
    - input_directory: 输入目录路径。
    - output_directory: 输出目录路径。
    - lpaq8_path: lpaq8压缩器的完整路径。
    - compression_level: 压缩级别（默认为9，范围0-9）。
    """
    # 记录开始时间
    start_time = time.time()

    # 确保输出目录存在
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    for root, dirs, files in os.walk(input_directory):
        for file in files:
            input_file_path = os.path.join(root, file)
            compressed_filename = f"{os.path.splitext(os.path.basename(input_file_path))[0]}.lpaq8"
            output_path = os.path.join(output_directory, compressed_filename)

            compress_file(input_file_path, output_path, lpaq8_path, compression_level)

    # 记录结束时间
    end_time = time.time()

    print(f"所有文件已压缩完成。总共耗时: {(end_time - start_time) / 60} 分钟。")


def decompress_all_files_in_directory(input_directory, output_directory, lpaq8_path):
    """
    解压目录中的所有文件。

    参数:
    - input_directory: 输入目录路径。
    - output_directory: 输出目录路径。
    - lpaq8_path: lpaq8压缩器的完整路径。
    - compression_level: 压缩级别（默认为9，范围0-9）。
    """
    # 记录开始时间
    start_time = time.time()

    # 确保输出目录存在
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    output_path_mapping = {
        "id_regex.lpaq8": "txt",
        "id_tokens.lpaq8": "txt",
        "base_g_prime.lpaq8": "tiff",
        "quality.lpaq8": "tiff"
    }

    for root, dirs, files in os.walk(input_directory):
        for file in files:
            error = True
            input_file_path = os.path.join(root, file)
            base_filename = os.path.splitext(os.path.basename(input_file_path))[0]

            for suffix in output_path_mapping.keys():
                if file.endswith(suffix):
                    output_filename = f"{base_filename}.{output_path_mapping[suffix]}"
                    output_path = os.path.join(output_directory, output_filename)
                    decompress_file(input_file_path, output_path, lpaq8_path)
                    error = False

            if error:
                print(f"未知文件类型: {file}")

    # 记录结束时间
    end_time = time.time()

    print(f"所有文件已解压完成。总共耗时: {(end_time - start_time) / 60} 分钟。")


def get_file_size(file_path):
    file_size = os.path.getsize(file_path)

    if file_size < 1024:
        return f"{file_size} bytes"
    elif file_size < 1024 * 1024:
        return f"{file_size / 1024:.2f} KB"
    elif file_size < 1024 * 1024 * 1024:
        return f"{file_size / (1024 * 1024):.2f} MB"
    else:
        return f"{file_size / (1024 * 1024 * 1024):.2f} GB"


def get_directory_size(directory_path):
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(directory_path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            total_size += os.path.getsize(fp)

    if total_size < 1024:
        return f"{total_size} bytes"
    elif total_size < 1024 * 1024:
        return f"{total_size / 1024:.2f} KB"
    elif total_size < 1024 * 1024 * 1024:
        return f"{total_size / (1024 * 1024):.2f} MB"
    else:
        return f"{total_size / (1024 * 1024 * 1024):.2f} GB"

def monitor_output_file(output_file):
    while True:
        file_size = os.path.getsize(output_file)
        print(f"Output file size: {file_size} bytes")
        time.sleep(1)  # 每隔一秒检查一次文件大小


if __name__ == '__main__':
    # 示例用法
    input_directory1 = r"D:\pythonProject\fastqtobmp\input\change_to_gray" # 定义需要压缩的文件路径
    destination_directory1 = r'D:\pythonProject\fastqtobmp\input\change_to_gray_lpaq8'  # 定义输出目录
    lpaq8_exe_path = f"{os.getcwd()}\lpaq8.exe"  # 确保这是正确的lpaq8路径


    input_destination = r"D:\pythonProject\fastqtobmp\input"
    output_destination = r"D:\pythonProject\fastqtobmp\output\1"

    output_file = os.path.join("output", "SRR554369")
    monitor_thread = threading.Thread(target=monitor_output_file, args=(output_file, ))
    monitor_thread.start()

    compress_file(os.path.join(os.getcwd(), "input", "SRR554369.fastq"), os.path.join(os.getcwd(), "output", "SRR554369"), lpaq8_exe_path)

    monitor_thread.join()

    # input_directory2 = r"D:\pythonProject\fastqtobmp\input\compressed" # 定义需要压缩的文件路径
    # destination_directory2 = r'D:\pythonProject\fastqtobmp\input\compressed_lpaq8'  # 定义输出目录

    # 压缩目录中的所有文件
    # compress_all_files_in_directory(input_directory1, destination_directory1, lpaq8_exe_path)
    # compress_all_files_in_directory(input_directory2, destination_directory2, lpaq8_exe_path)

    # 计算输出目录的大小，并转换为MB
    # size1 = get_directory_size(destination_directory1)
    # size2 = get_directory_size(destination_directory2)

    # 输出两个目录大小的比较结果
    # difference = size1 - size2
    # print(f"{destination_directory1} 比 {destination_directory2} 大了 {difference:.2f} MB。")
