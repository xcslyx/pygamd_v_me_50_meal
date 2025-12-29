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
        print(f"✅ 创建 {folder_name} 文件夹...")
    elif overwrite:
        print(f"⚠️ {folder_name} 文件夹已存在，正在删除...")
        shutil.rmtree(os.path.join(folder_path, folder_name))
        os.mkdir(os.path.join(folder_path, folder_name))
        print(f"✅ 创建 {folder_name} 文件夹...")
    else:
        print(f"❌ {folder_name} 文件夹已存在")


def check_xml_start_tag(xml_file):
    if (
            xml_file.startswith("particles") or xml_file.startswith("monomer") or xml_file.startswith("simulation")
    ) or (xml_file.endswith("0.xml") and os.path.isfile(xml_file)):
        return True
    return None


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
