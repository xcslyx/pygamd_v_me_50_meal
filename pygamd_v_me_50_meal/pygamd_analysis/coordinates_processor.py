import os
import re
import shutil
import logging
import xml.etree.ElementTree as ET

import numpy as np

from tqdm import tqdm
from multiprocessing import Pool

import pygamd_v_me_50_meal.utils as utils
from pygamd_v_me_50_meal.Functions import Functions

from .get_sequence import GetSequence


# 处理坐标文件
class CoordinatesProcessor:
    def __init__(self, path: str, data, remove_ions_zhy: bool=False,
                 remove_enm=None, remove_condensate_pbc=False):
        self.path = path
        self.data = data

        self.length_dict = self.data.length_dict

        self.init_xml_path = os.path.join(self.path, "xml/")
        self.unwrapping_xml_path = os.path.join(self.path, "xml_unwrapping/")
        self.remove_pbc_condensate_xml_path = os.path.join(self.path, "xml_remove_pbc_condensate/")
        self.unwrapping_remove_ions_xml_path = os.path.join(self.path, "xml_unwrapping_remove_ions/")
        self.remove_pbc_condensate_remove_ions_xml_path = os.path.join(self.path, "xml_remove_pbc_condensate_remove_ions/")

        os.makedirs(self.init_xml_path, exist_ok=True)
        for file in os.listdir(self.path):
            if file.startswith("particles") and file.endswith(".xml"):
                file_path = os.path.join(self.path, file)
                destination_path = os.path.join(self.init_xml_path, file)

                # 在移动之前检查目标文件是否存在，如果存在则删除
                if os.path.exists(destination_path):
                    os.remove(destination_path)

                shutil.move(file_path, self.init_xml_path)

        self.remove_ions_zhy = remove_ions_zhy

        if remove_enm:
            self.remove_enm_bonds_request = "y"
        elif remove_enm is None:
            self.remove_enm_bonds_request = input("是否需要移除弹性键？(y/n)")
        else:
            self.remove_enm_bonds_request = "n"

        xml_path = os.path.join(self.path, "xml")
        self.xml_files = sorted(os.listdir(xml_path))
        for i in range(len(self.xml_files)):
            if self.xml_files[i].startswith("particles") and self.xml_files[i].endswith("0.xml"):
                init_xml_file = os.path.join(self.path, "xml", self.xml_files[i])
                break
        init_root = ET.parse(init_xml_file).getroot()
        self.box_size: list[float] = [float(init_root.find('.//box').attrib[i]) for i in ['lx', 'ly', 'lz']]

        self.unwrap_xml_flag: bool = True


        self.cal_xyz(remove_condensate_pbc=remove_condensate_pbc)


    def remove_enm_bonds_from_xml(self, xml_file):
        if xml_file.startswith("particles") and xml_file.endswith("0.xml"):
            pass
        else:
            return

        tree = ET.parse(os.path.join(self.init_xml_path, xml_file))
        root = tree.getroot()

        for bonds in root.findall('.//bond'):
            modified_bonds = []
            for line in bonds.text.strip().split('\n'):
                pattern = r"C\d+-C\d+"
                match = re.search(pattern, line)
                # 检查是否匹配
                if not match:
                    modified_bonds.append(line)
            # 将修改后的文本重新写回到 bond 元素中
            bonds.text = '\n' + '\n'.join(modified_bonds) + '\n'

        # 将修改后的 XML 对象写回到文件中
        tree.write(os.path.join(self.init_xml_path, xml_file), encoding='utf-8', xml_declaration=True)


    def abstract_coordinates_normal(self, xml_file):
        if not (xml_file.startswith("particles") and xml_file.endswith("0.xml")):
            return

        positions, unwrapped_positions = {}, {}

        tree = ET.parse(os.path.join(self.init_xml_path, xml_file))
        root = tree.getroot()

        # 读取 position 数据
        position_elem = root.find('.//position')
        position_elem_text = position_elem.text
        while "\n\n" in position_elem_text:
            position_elem_text = position_elem_text.replace("\n\n", "\n")
        position, unwrapped_position = [], []
        for line in position_elem_text.strip().split('\n'):
            position.append(list(map(float, line.strip().split())))


        for type_ in self.data.mol_class_dict:
            # if type_ in ["Na", "K", "Cl", "Li", ] and self.remove_ions_zhy:
            #     continue

            mol_positions, unwrapped_mol_positions = [], []
            count = self.data.mol_class_dict[type_][0]
            length = self.data.mol_class_dict[type_][1]

            for _ in range(count):
                cur_position = position[:length]
                mol_positions.append(cur_position)
                position = position[length:]

                if Functions.is_chain_wrapped(cur_position, self.box_size):
                    unwrapped_cur_position = Functions.unwrap_chain(cur_position, self.box_size)
                else:
                    unwrapped_cur_position = cur_position

                unwrapped_mol_positions.append(unwrapped_cur_position)
                unwrapped_position.extend(unwrapped_cur_position)

            positions[type_] = mol_positions
            unwrapped_positions[type_] = unwrapped_mol_positions

        with open(os.path.join(self.path, "chain_xyz/",
                               xml_file), 'w') as f_chain_xyz:
            f_chain_xyz.writelines(str(positions))

        with open(os.path.join(self.path, "chain_xyz_unwrapping/",
                               xml_file.replace("0.xml", "0.reimage.xml")), 'w') as f_chain_xyz:
            f_chain_xyz.writelines(str(unwrapped_positions))

            # 重新写入 position 数据
            position_elem.text = '\n' + '\n'.join([' '.join(map(str, pos)) for pos in unwrapped_position]) + '\n'
            position_elem.attrib['num'] = str(len(unwrapped_position))

            root.find('.//configuration').attrib["natoms"] = str(len(unwrapped_position))


            # 删除 angle dihedral mass charge body image 数据
            # for tag in ["angle", "dihedral", "mass", "charge", "body", "image", "velocity"]:
            for tag in ["image", "velocity"]:
                try:
                    for elem in root.findall(f'.//{tag}'):
                        # 删除 configuration 元素中的相关数据
                        root.find('.//configuration').remove(elem)
                except Exception as e:
                    pass

            # 写入修改后的 XML 对象
            tree.write(os.path.join(self.unwrapping_xml_path, xml_file.replace("0.xml", "0.reimage.xml")),
                       encoding='utf-8', xml_declaration=True)


    def abstract_coordinates_condensate_pbc(self, xml_file):
        remove_pbc_condensate_xml_path = self.remove_pbc_condensate_remove_ions_xml_path if self.remove_ions_zhy else self.remove_pbc_condensate_xml_path
        if os.path.exists(remove_pbc_condensate_xml_path):
            tree = ET.parse(os.path.join(remove_pbc_condensate_xml_path, xml_file.replace(".xml", ".reimage.new.xml")))
            root = tree.getroot()

            # 读取 position 数据
            position_elem = root.find('.//position')
            position = []
            for line in position_elem.text.strip().split('\n'):
                position.append(list(map(float, line.strip().split())))

            positions = {}
            for type_ in self.data.mol_class_dict:
                mol_positions = []
                count = self.data.mol_class_dict[type_][0]
                length = self.data.mol_class_dict[type_][1]

                for _ in range(count):
                    cur_position = position[:length]
                    mol_positions.append(cur_position)
                    position = position[length:]


                positions[type_] = mol_positions

            with open(os.path.join(self.path, "chain_xyz_remove_pbc_condensate/", xml_file), 'w') as f_chain_xyz:
                f_chain_xyz.writelines(str(positions))


    def remove_ions(self, xml_file):
        if not (xml_file.startswith("particles") and xml_file.endswith("0.reimage.xml")):
            return

        tree = ET.parse(os.path.join(self.unwrapping_xml_path, xml_file))
        root = tree.getroot()

        for tag in ["position", "type", "mass", "charge", "body", "image", "velocity"]:
            for elem in root.findall(f'.//{tag}'):
                elem_text_list = elem.text.strip().split('\n')
                for _type in self.data.mol_class_dict:
                    if _type in ["Na", "K", "Cl", "Li", ]:
                        # 用空字符代替原来的元素
                        start, end = self.data.mol_class_dict[_type][2][0], self.data.mol_class_dict[_type][2][1]
                        elem_text_list[start:end] = [''] * (end - start + 1)

                elem_text_list = [elem_text for elem_text in elem_text_list if elem_text]
                elem.text = '\n' + '\n'.join(elem_text_list) + '\n'
                elem.attrib['num'] = str(len(elem_text_list))

                if tag == "position":
                    root.find('.//configuration').attrib["natoms"] = str(len(elem_text_list))

        # 写入修改后的 XML 对象
        tree.write(os.path.join(self.unwrapping_remove_ions_xml_path, xml_file), encoding='utf-8', xml_declaration=True)


    def cal_xyz(self, remove_condensate_pbc=False):
        init_files = sorted([i for i in os.listdir(self.init_xml_path) if i.endswith("0.xml") and i.startswith("particles")])

        print("开始提取序列信息...")
        seq_output = f"{self.data.system_name}_sequence.txt"
        GetSequence(self.init_xml_path, sorted(init_files)[0], self.data, output_path=self.path,
                    output=seq_output).xml2sequence()
        print(f"序列信息已保存至 {seq_output}。")

        if self.remove_enm_bonds_request == "y":
            utils.backup_folder(self.path, 'xml', 'xml_init')
            # 使用 tqdm 和多进程移除弹性键
            with Pool(processes=4) as pool:
                list(tqdm(pool.imap(self.remove_enm_bonds_from_xml, init_files),
                          total=len(init_files),
                          desc="移除弹性键进度", colour="cyan", ncols=100))


        for _dir in ["xml_unwrapping", "chain_xyz", "chain_xyz_unwrapping"]:
            utils.create_folder(_dir, self.path, overwrite=True)

        with Pool(processes=4) as pool:
            print("开始去周期并提取坐标...")
            # 使用 tqdm 包装你的可迭代对象
            list(tqdm(pool.imap(self.abstract_coordinates_normal, init_files),
                      total=len(init_files),
                      desc="处理进度",
                      colour="cyan",
                      ncols=100))

        if self.remove_ions_zhy:
            unwrapping_files = sorted([i for i in os.listdir(self.unwrapping_xml_path) if i.endswith("0.reimage.xml") and i.startswith("particles")])
            # utils.backup_folder(self.path, 'xml', 'xml_init')
            utils.create_folder("xml_unwrapping_remove_ions", self.path, overwrite=True)
            # 使用 tqdm 和多进程移除离子
            with Pool(processes=4) as pool:
                list(tqdm(pool.imap(self.remove_ions, unwrapping_files),
                          total=len(unwrapping_files),
                          desc="移除离子进度", colour="cyan", ncols=100))

        if remove_condensate_pbc:
            print("开始进行针对凝聚体的 PBC 去除")
            for _dir in ["chain_xyz_remove_pbc_condensate"]:
                utils.create_folder(_dir, self.path, overwrite=True)
            self.remove_pbc_condensate_parallel()

        print("所有文件处理完成。")

    def remove_pbc_condensate(self, xml_file):
        """
        去除凝聚体的 PBC 条件，使得凝聚体的坐标不再受到 PBC 条件的限制。
        :param xml_file: 凝聚体的 XML 文件
        :return: None
        """
        if not (xml_file.startswith("particles") and xml_file.endswith("0.reimage.xml")):
            return

        if self.remove_ions_zhy:
            tree = ET.parse(os.path.join(self.unwrapping_remove_ions_xml_path, xml_file))
        else:
            tree = ET.parse(os.path.join(self.unwrapping_xml_path, xml_file))
        root = tree.getroot()

        # 读取 position 数据
        position = []
        position_elem = root.find('.//position')
        for line in position_elem.text.strip().split('\n'):
            position.append(list(map(float, line.strip().split())))

        positions = []
        for type_ in self.data.mol_class_dict:
            if type_ in ["Na", "K+", "Cl", "Br", "I-"]:
                continue
            count = self.data.mol_class_dict[type_][0]
            length = self.data.mol_class_dict[type_][1]
            for _ in range(count):
                positions.append(np.array(position[:length]))
                position = position[length:]
        del position
        # print(len(positions))

        indexed_positions = [(i, pos) for i, pos in enumerate(positions)]

        not_ordered_positions, max_group_core = Functions.pre_process(indexed_positions, self.box_size)

        for i in range(len(not_ordered_positions)):
            not_ordered_positions[i] = (not_ordered_positions[i][0], Functions.move_chain_into_box(not_ordered_positions[i][1], self.box_size, max_group_core))
        # 按照 idx 排序并提取 position
        re_ordered_positions = [position for idx, position in sorted(not_ordered_positions, key=lambda x: x[0])]

        # 更新后的蛋白质位置
        new_positions = np.vstack(re_ordered_positions)

        new_xml_file = xml_file.replace(".xml", ".new.xml")
        if self.remove_ions_zhy:
            tree = ET.parse(os.path.join(self.unwrapping_remove_ions_xml_path, xml_file))
        else:
            tree = ET.parse(os.path.join(self.unwrapping_xml_path, xml_file))
        root = tree.getroot()

        # 找到 position 标签
        position_element = root.find('.//position')
        pos_text = "\n".join(["    ".join(map(str, p)) for p in new_positions])
        position_element.text = f"\n{pos_text}\n"
        # 保存更新后的XML文件
        if self.remove_ions_zhy:
            tree.write(os.path.join(self.remove_pbc_condensate_remove_ions_xml_path, new_xml_file), encoding='utf-8', xml_declaration=True)
        else:
            tree.write(os.path.join(self.remove_pbc_condensate_xml_path, new_xml_file), encoding='utf-8', xml_declaration=True)


    def remove_pbc_condensate_parallel(self, ):
        """
        多进程去除凝聚体 PBC 条件。
        """
        if self.remove_ions_zhy:
            if not os.path.exists(self.remove_pbc_condensate_remove_ions_xml_path):
                os.makedirs(self.remove_pbc_condensate_remove_ions_xml_path, exist_ok=True)
                if not os.path.exists(self.unwrapping_remove_ions_xml_path):
                    print("未找到去周期后的 XML 文件，请先运行去周期脚本。")
                    return
            xml_files = sorted(os.listdir(self.unwrapping_remove_ions_xml_path))
        else:
            if not os.path.exists(self.remove_pbc_condensate_remove_ions_xml_path):
                os.makedirs(self.remove_pbc_condensate_remove_ions_xml_path, exist_ok=True)
                if not os.path.exists(self.unwrapping_xml_path):
                    print("未找到去周期后的 XML 文件，请先运行去周期脚本。")
                    return
            xml_files = sorted(os.listdir(self.unwrapping_xml_path))


        print(f"开始处理 {len(xml_files)} 个文件...")
        with Pool(processes=4) as pool:
            list(tqdm(pool.imap(self.remove_pbc_condensate, xml_files),
                      total=len(xml_files),
                      desc="去除 PBC 条件",
                      colour='cyan',
                      ncols=100))

        init_files = sorted(
            [i for i in os.listdir(self.init_xml_path) if i.endswith("0.xml") and i.startswith("particles")])
        with Pool(processes=4) as pool:
            print("开始处理去周期后的凝聚体文件...")
            # 使用 tqdm 包装你的可迭代对象
            list(tqdm(pool.imap(self.abstract_coordinates_condensate_pbc, init_files),
                      total=len(init_files),
                      desc="处理进度",
                      colour="cyan",
                      ncols=100))
