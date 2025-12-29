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


def create_folder(folder_name, folder_path, overwrite=False):
    if folder_name not in os.listdir(folder_path):
        os.mkdir(os.path.join(folder_path, folder_name))
        print(f"✅ Create folder: {folder_name}...")
    elif overwrite:
        print(f"⚠️ {folder_name} exists, deleting...")
        shutil.rmtree(os.path.join(folder_path, folder_name))
        os.mkdir(os.path.join(folder_path, folder_name))
        print(f"✅ Create folder: {folder_name}...")
    else:
        print(f"❌ {folder_name} exists")


def check_xml_start_tag(xml_file):
    if xml_file.startswith("particles") or xml_file.startswith("monomer") or xml_file.startswith("simulation"):
        return True
    return None


def chain_in_box(chain_positions, box_size) -> list[bool]:
    box_length = list(map(lambda x: x / 2, box_size))
    in_box = [True, True, True]
    for pos in chain_positions:
        for i in range(3):
            if pos[i] < -box_length[i] or pos[i] > box_length[i]:
                in_box[i] = False
    return in_box


def backup_folder(backup_path, init_folder_name, backup_folder_name):
    if not os.path.exists(os.path.join(backup_path, backup_folder_name)):
        print("正在备份原始 XML 文件...")
        shutil.copytree(os.path.join(backup_path, init_folder_name), os.path.join(backup_path, backup_folder_name))
        print("备份完成。")
    else:
        if input("是否需要重新备份原始 XML 文件？(y/n)") == "y":
            print("正在备份原始 XML 文件...")
            shutil.rmtree(os.path.join(backup_path, backup_folder_name))
            shutil.copytree(os.path.join(backup_path, init_folder_name), os.path.join(backup_path, backup_folder_name))
            print("备份完成。")
