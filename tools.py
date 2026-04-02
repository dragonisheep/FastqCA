import os


def check_output_path(output_path):
    if not os.path.exists(output_path):
        os.mkdir(output_path)
        return True
    elif not os.path.isdir(output_path):
        return True
    else:
        return False


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


def show_time(elapsed_time, des):
    hours, remainder = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(remainder, 60)

    if hours >= 1:
        print(f"{des}——程序运行时间：{int(hours)}小时")
    elif minutes >= 1:
        print(f"{des}——程序运行时间：{int(minutes)}分钟")
    else:
        print(f"{des}——程序运行时间：{int(seconds)}秒钟")
