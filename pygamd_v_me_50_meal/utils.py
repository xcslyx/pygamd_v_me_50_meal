import os
import shutil


def str2value(value):
    if value.lower() in ('yes', 'y', 'true', 'ture', 't', '1'):
        return True
    elif value.lower() in ('no', 'n', 'false', 'f', '0'):
        return False
    elif value == "avg":
        return "avg"
    elif value == "unset":
        return None
    else:
        return value


def create_folder(folder_path, path, overwrite=False):
    if folder_path not in os.listdir(path):
        os.mkdir(os.path.join(path, folder_path))
        print(f"✅ 创建{folder_path}文件夹...")
    elif overwrite:
        print(f"⚠️ {folder_path}文件夹已存在，正在删除...")
        shutil.rmtree(os.path.join(path, folder_path))
        os.mkdir(os.path.join(path, folder_path))
        print(f"✅ 创建{folder_path}文件夹...")
    else:
        print(f"❌ {folder_path}文件夹已存在")