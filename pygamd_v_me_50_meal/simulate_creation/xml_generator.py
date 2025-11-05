#!-*- coding: utf-8 -*-
import os
import re
import math
import shutil

import xml.etree.ElementTree as ET

import numpy as np

# 定义类 XMLGenerator
class XMLGenerator:
    def __init__(self, path: str, filename, box_size: str | float, add_enm_bond=None, add_rigid_body=None, add_domain=None, dna_model=None, gen_run_file=None):
        if filename.endswith(".pdb"):
            if path is None:
                if os.path.exists(filename):
                    path = os.path.dirname(filename)
                    if not path:
                        path = os.getcwd()
                    print(f"系统目录路径未提供，已设置为 {path}")
                    filename = os.path.basename(filename)
                else:
                    raise ValueError("请提供体系路径！")
            print(f"开始转换 PDB 文件 {filename} 为 XML 文件...")
        else:
            raise ValueError("PDB 文件名格式不正确！")

        self.path = path
        self.filename = filename
        # pro_res_map 前四行为标准氨基酸，后面为全原子力场中的非标准氨基酸，如有缺少可自行添加
        self.pro_res_map = {'ALA': ['A', 5], 'ARG': ['R', 11], 'ASN': ['N', 8], 'ASP': ['D', 8], 'CYS': ['C', 6],
                        'GLN': ['Q', 9], 'GLU': ['E', 9], 'GLY': ['G', 4], 'HIS': ['H', 10], 'ILE': ['I', 8],
                        'LEU': ['L', 8], 'LYS': ['K', 9], 'MET': ['M', 8], 'PHE': ['F', 11], 'PRO': ['P', 7],
                        'SER': ['S', 6], 'THR': ['T', 7], 'TRP': ['W', 14], 'TYR': ['Y', 12], 'VAL': ['V', 7],
                        'HIE': ['H', 10]}
        self.mass_dict_atom = {"C": 12.011, "N": 14.007, "O": 15.999, "P": 30.974, 'H': 1.008, 'S': 32}
        self.box_size = box_size
        self.x_max, self.y_max, self.z_max = 0, 0, 0

        self.log_dir = ""
        self.xml_dir = ""

        self.mol_num = 0
        self.mol_class = []

        if add_enm_bond:
            self.add_enm_bond_flag = True
            self.enm_domain_list = add_enm_bond.split(',')
        else:
            self.add_enm_bond_flag = False

        if add_rigid_body:
            self.add_rigid_body = True
            self.rigid_domain_list = add_rigid_body.split(',')
        else:
            self.add_rigid_body = False
        self.rigid_body_index = 0

        if add_domain:
            self.add_domain_flag = True
        else:
            self.add_domain_flag = False


        self.protein_model = None
        self.dna_model = None
        self.rna_model = None
        if dna_model:
            if dna_model == "1":
                self.dna_model = "3SPN"
            elif dna_model == "2":
                self.dna_model = "2BeadMittal"

        self.output_file = filename.replace(".pdb", ".xml")
        self.output_tree = ET.ElementTree(ET.Element("galamost_xml", {"version": "1.0"}))
        self.pdb2xml()
        # 合并xml文件
        self.merge_xml()
        print("PDB 转换为 XML 完成！")

        if gen_run_file is None:
            run_file_request = input("是否要生成对应的 PYGAMD 运行文件？(y(es)/n(o))：")
            if run_file_request.lower() in ["y", "yes"]:
                print("正在生成对应的 PYGAMD 运行文件...")
                self.generate_pygamd_run_file()
        elif gen_run_file:
            print("正在生成对应的 PYGAMD 运行文件...")
            self.generate_pygamd_run_file()


    def pdb2xml(self):
        """
        将 PDB 文件转换为包含 <sequence> 和 <position> 的log文件。
        PDB 文件包含由 "TER" 标签分隔的多个分子。
        """
        with open(os.path.join(self.path, self.filename), 'r') as pdb:
            pdb_lines = pdb.readlines()
            atoms = []
            molecules = []
            for line in pdb_lines:
                if "ATOM" in line[:4]:
                    # 分别为ATOM、原子序号、原子名称、RES名称、链标识符、RES序号、X、Y、Z、Occupancy、TempFactor、元素名称
                    cur_atom = list(map(lambda x: x.strip(), [line[:6], line[6:11], line[12:16], line[17:20], line[21], line[22:26], line[30:38],
                                line[38:46], line[46:54], line[54:60], line[60:66], line[76:78]]))
                    # print(cur_atom)
                    if cur_atom[-1] != 'H' and cur_atom[2] not in ["OC2", "OXT"]:
                        atoms.append(cur_atom)
                if "TER" in line:
                    if atoms:
                        molecules.append(atoms)
                    atoms = []
                    continue
        self.mol_num = len(molecules)
        self.mol_class = [""] * self.mol_num

        self.log_dir = os.path.join(self.path, "log")
        self.xml_dir = os.path.join(self.path, "xml")
        dirs = os.listdir(self.path)
        if "log" not in dirs:
            os.mkdir(os.path.join(self.path, "log"))
        else:
            shutil.rmtree(os.path.join(self.path, "log"))
            os.mkdir(os.path.join(self.path, "log"))
        if "xml" not in dirs:
            os.mkdir(os.path.join(self.path, "xml"))
        else:
            shutil.rmtree(os.path.join(self.path, "xml"))
            os.mkdir(os.path.join(self.path, "xml"))

        for i in range(self.mol_num):
            root = ET.Element("molecule")
            root.append(ET.Element("sequence"))
            root.append(ET.Element("position"))

            root_tree = ET.ElementTree(root)
            seq_elem = root_tree.find("sequence")
            seq_elem.text = "\n"
            pos_elem = root_tree.find("position")
            pos_elem.text = "\n"

            for atom_num in range(len(molecules[i])):
                try:
                    if molecules[i][atom_num][5] != molecules[i][atom_num + 1][5]:
                        seq_elem.text += molecules[i][atom_num][3] + '\n'
                except IndexError:
                    pass
            seq_elem.text += molecules[i][-1][3] + '\n'
            if any(i in seq_elem.text.split('\n') for i in self.pro_res_map.keys()):
                self.mol_class[i] = "pro"
                if not self.protein_model:
                    print("请选择蛋白质模型：\n1. HPS 模型（HPS-Urry、CALVADOS系列等）\n2. Mpipi 模型")
                    dna_model_choice = input("请输入选择：")
                    if dna_model_choice == "1":
                        self.protein_model = "HPS"
                    elif dna_model_choice == "2":
                        self.protein_model = "Mpipi"
            elif any(i in seq_elem.text.split('\n') for i in ['DA', 'DT', 'DC', 'DG', 'DU']):
                self.mol_class[i] = "dna"
                if not self.dna_model:
                    print("请选择 DNA 模型：\n1. 3SPN 模型\n2. Mittal 2 Beads 模型")
                    dna_model_choice = input("请输入选择：")
                    if dna_model_choice == "1":
                        self.dna_model = "3SPN"
                    elif dna_model_choice == "2":
                        self.dna_model = "2BeadMittal"
            elif any(i in seq_elem.text.split('\n') for i in ['A', 'C', 'G', 'U']):
                self.mol_class[i] = "rna"
                if not self.rna_model:
                    print("请选择 RNA 模型：\n1. 3SPN 模型\n2. CAVADOS-RNA 模型(未提供)")
                    rna_model_choice = input("请输入选择：")
                    if rna_model_choice == "1":
                        self.rna_model = "3SPN"
                    elif rna_model_choice == "2":
                        self.rna_model = "CAVADOS-RNA"

            res = ""
            k = 1
            for atom in molecules[i]:
                if atom[5] == res or not res:
                    pass
                else:
                    k += 1
                pos_elem.text += "{:5}\t\t{:10}\t{:10}\t{:10}\t{}     {}\n".format(
                    atom[2], atom[6], atom[7], atom[8], k, atom[3])
                res = atom[5]
            with open(os.path.join(self.path, f"log/mol{i}_log.xml"), 'w') as f:
                f.write(f"<?xml version='1.0' encoding='utf-8'?>\n<molecule num_atoms='{k}'>")
                f.write(f"\n<{seq_elem.tag}>" + seq_elem.text + f"</{seq_elem.tag}>\n")
                f.write(f"<{pos_elem.tag}>" + pos_elem.text + f"</{pos_elem.tag}>\n")
                f.write(f"</molecule>")

        for i in range(len(self.mol_class)):
            if self.mol_class[i] == "pro":
                if not self.add_rigid_body:
                    rigid_request = input("是否要对结构域设置刚体？(y(es)/n(o))：")
                    if rigid_request.lower() in ["y", "yes"]:
                        self.add_rigid_body = True
                        self.rigid_domain_list = input("请输入该结构域的起始残基编号和末尾残基编号（从 1 开始），以-分隔，如 159-522，若有多个结构域，请以英文逗号分隔。\n").split(',')
                        # self.add_rigid_body = False

                self.pro_log2xml(file=f"mol{i}_log.xml")

                if self.add_rigid_body:
                    if not self.add_domain_flag:
                        domain_request = input("是否要对结构域单独设置粒子类型？(y(es)/n(o))：")
                        if domain_request.lower() in ["y", "yes"]:
                            self.add_domain_flag = True
                        else:
                            self.add_rigid_body = False

                    if self.add_domain_flag:
                        for domain_idx in range(len(self.rigid_domain_list)):
                            domain = list(map(int, self.rigid_domain_list[domain_idx].split('-')))
                            self.pro_add_domain(f"mol{i}_log.xml", domain[0], domain[1])
                        self.add_domain_flag = False
                        self.add_rigid_body = False
                else:
                    if not self.add_enm_bond_flag:
                        enm_request = input("是否要对结构域设置弹性网络？(y(es)/n(o))：")
                        if enm_request.lower() in ["y", "yes"]:
                            self.add_enm_bond_flag = True
                    if self.add_enm_bond_flag:
                        self.enm_domain_list = input("请输入该结构域的起始残基编号和末尾残基编号（从 1 开始），以-分隔，如 159-522，若有多个结构域，请以英文逗号分隔。\n").split(',')
                        for domain_idx in range(len(self.enm_domain_list)):
                            domain = list(map(int, self.enm_domain_list[domain_idx].split('-')))
                            if domain_idx == 0:
                                self.add_enm_bond(f"mol{i}_log.xml", domain[0], domain[1])
                            else:
                                self.add_enm_bond(f"mol{i}_log_enm.xml", domain[0], domain[1])
                        self.add_enm_bond_flag = False

                        if not self.add_domain_flag:
                            domain_request = input("是否要对结构域单独设置粒子类型？(y(es)/n(o))：")
                            if domain_request.lower() in ["y", "yes"]:
                                self.add_domain_flag = True

                        if self.add_domain_flag:
                            for domain_idx in range(len(self.enm_domain_list)):
                                domain = list(map(int, self.enm_domain_list[domain_idx].split('-')))
                                self.pro_add_domain(f"mol{i}_log_enm.xml", domain[0], domain[1])
                            self.add_domain_flag = False

            elif self.mol_class[i] == "dna":
                if self.dna_model == "3SPN":
                    self.dna_log2xml_3SPN(file=f"mol{i}_log.xml")
                elif self.dna_model == "2BeadMittal":
                    self.dna_log2xml_2BeadMittal(file=f"mol{i}_log.xml")

            elif self.mol_class[i] == "rna":
                if self.rna_model == "3SPN":
                    self.rna_log2xml_3SPN(file=f"mol{i}_log.xml")

        shutil.rmtree(os.path.join(self.path, "log"))


    def merge_xml(self):
        xml_files = os.listdir(self.xml_dir)
        if self.mol_num == 1:
            for file in xml_files:
                filename = os.path.join(self.xml_dir, file)
                if filename.endswith("enm.xml"):
                    self.output_file = os.path.join(self.path, self.filename.replace(".pdb", "_enm.xml"))
                    shutil.move(filename, str(self.output_file))
                else:
                    self.output_file = os.path.join(self.path, self.filename.replace(".pdb", ".xml"))
                    shutil.move(filename, str(self.output_file))
                    # TODO
        else:
            root = self.output_tree.getroot()
            root.text = '\n'
            root.append(ET.Element("configuration", {"time_step": "0", "dimensions": "3", "natoms": '0'}))
            root_configuration = root.find("configuration")
            root_configuration.text = '\n'
            root_configuration_tags = ["box"]
            # print(root_configuration.iter())
            # print(root.tag)
            for file in sorted(xml_files):
                file = os.path.join(self.xml_dir, file)
                file_tree = ET.parse(file)
                file_root = file_tree.getroot()

                for tag in file_root.iter():
                    if tag.tag == "galamost_xml":
                        continue
                    if tag.tag == "configuration":
                        cur_natoms = int(root_configuration.get("natoms"))
                        root_configuration.set("natoms", str(int(cur_natoms) + int(tag.get("natoms"))))
                        continue
                    if tag.tag == "box":
                        box_size = max(self.x_max, self.y_max, self.z_max)
                        # 将 box_size 向上变成最接近的 5 或 10 的倍数
                        box_size = box_size // 5 * 5 + 5
                        box_elem = root_configuration.find("box")
                        if box_elem is None:
                            root_configuration.append(ET.Element("box", {"lx": str(box_size), "ly": str(box_size), "lz": str(box_size)}))
                        else:
                            box_elem.set("lx", str(box_size))
                            box_elem.set("ly", str(box_size))
                            box_elem.set("lz", str(box_size))
                        print(f"Box size: {box_size}")
                    # 如果这个 tag 是第一次出现，则直接添加到 configuration 元素中
                    if tag.tag not in root_configuration_tags:
                        root_configuration.append(tag)
                        root_configuration_tags.append(str(tag.tag))
                    elif tag.tag == "bond":
                        cur_bond_list = tag.text.strip().split('\n')
                        root_bond_elem = root_configuration.find("bond")
                        root_bond_elem.set("num", str(int(root_bond_elem.get("num")) + int(tag.get("num"))))
                        for bond in cur_bond_list:
                            bond = bond.split()
                            bond[1] = str(int(bond[1]) + cur_natoms)
                            bond[2] = str(int(bond[2]) + cur_natoms)
                            bond = " ".join(bond)
                            root_bond_elem.text += bond + "\n"
                    elif tag.tag == "angle":
                        cur_angle_list = tag.text.strip().split('\n')
                        root_angle_elem = root_configuration.find("angle")
                        root_angle_elem.set("num", str(int(root_angle_elem.get("num")) + len(cur_angle_list)))
                        for angle in cur_angle_list:
                            angle = angle.split()
                            angle[1] = str(int(angle[1]) + cur_natoms)
                            angle[2] = str(int(angle[2]) + cur_natoms)
                            angle[3] = str(int(angle[3]) + cur_natoms)
                            angle = " ".join(angle)
                            root_angle_elem.text += angle + "\n"
                    elif tag.tag == "dihedral":
                        cur_dihedral_list = tag.text.strip().split('\n')
                        root_dihedral_elem = root_configuration.find("dihedral")
                        root_dihedral_elem.set("num", str(int(root_dihedral_elem.get("num")) + len(cur_dihedral_list)))
                        for dihedral in cur_dihedral_list:
                            dihedral = dihedral.split()
                            dihedral[1] = str(int(dihedral[1]) + cur_natoms)
                            dihedral[2] = str(int(dihedral[2]) + cur_natoms)
                            dihedral[3] = str(int(dihedral[3]) + cur_natoms)
                            dihedral[4] = str(int(dihedral[4]) + cur_natoms)
                            dihedral = " ".join(dihedral)
                            root_dihedral_elem.text += dihedral + "\n"
                    else:
                        cur_tag = root_configuration.find(tag.tag)
                        if cur_tag.text:
                            cur_tag.text += tag.text
                        if cur_tag.get("num") and tag.get("num"):
                            cur_tag.set("num", str(int(cur_tag.get("num")) + int(tag.get("num"))))

            self.output_tree.write(os.path.join(self.path, self.output_file), xml_declaration=True, encoding='utf-8', method='xml')
        shutil.rmtree(os.path.join(self.path, "xml"))


    def add_enm_bond(self, file, start, end):
        filename = os.path.join(self.xml_dir, file)
        if filename.endswith("enm.xml"):
            new_filename = filename
        else:
            new_filename = f'{filename[:-4]}_enm.xml'

        tree = ET.parse(filename)
        root = tree.getroot()
        root = root.find("configuration")

        pos_elem = root.find("position")
        positions = list(map(lambda x: list(map(float, x.split())), pos_elem.text.strip('\n').split('\n')))

        bond_elem = root.find("bond")
        init_bond_num = int(bond_elem.attrib["num"])

        enm_bond_file = os.path.join(self.path, f'{self.filename.replace(".pdb", "")}_enm_bond.py')

        if os.path.exists(enm_bond_file):
            read_mode = 'a'
        else:
            read_mode = 'w'

        with open(enm_bond_file, read_mode) as bond_file:
            bond_cnt = 0
            for i in range(start - 1, end):
                for j in range(i + 3, end):
                    distance = np.linalg.norm(np.array(positions[i]) - np.array(positions[j]))
                    if distance <= 0.9:
                        bond_cnt += 1
                        bond_type = f'C{i + 1}-C{j + 1}'
                        bond_elem.text += f"{bond_type} {i} {j}\n"
                        k_value = 500  # 默认 k 值为 500
                        command = f"bondforce_pro.setParams('{bond_type}', {k_value}, {distance})\n"
                        bond_file.write(command)
        # print(f"已生成弹性网络设置脚本：{os.path.join(self.path, f'{self.filename.replace(".pdb", "")}_enm_bond.py')}")
        bond_elem.attrib["num"] = str(init_bond_num + bond_cnt)
        tree.write(new_filename, xml_declaration=True)
        if not filename.endswith("enm.xml"):
            os.remove(filename)


    def pro_add_domain(self, file, start, end):
        filename = os.path.join(self.xml_dir, file)
        new_filename = f'{filename[:-4]}_add_domain.xml'

        tree = ET.parse(filename)
        root = tree.getroot()
        root = root.find("configuration")

        type_elem = root.find("type")
        type_list = type_elem.text.strip('\n').split('\n')
        type_list[start - 1:end] = [i +'D' for i in type_list[start - 1:end]]

        type_elem.text = '\n'.join(type_list) + '\n'

        tree.write(new_filename, xml_declaration=True)
        os.remove(filename)


    def pro_log2xml(self, file):
        tree = ET.parse(os.path.join(os.path.join(self.log_dir, file)))
        root = tree.getroot()

        seq_elem = root.find('sequence')
        sequence = seq_elem.text.strip('\n').split('\n')
        positions_elem = root.find('position')
        positions = positions_elem.text.strip('\n').split('\n')

        sequence_position: list[str, list[str]] = []
        for i in sequence:
            sequence_position.append([i])
        for line in positions:
            line = line.split()
            sequence_position[int(line[-2]) - 1].append(line)

        new_sequence_position = []
        for new_seq in sequence_position:
            m, x, y, z = 0, 0, 0, 0
            for i in new_seq[1:]:
                m += self.mass_dict_atom[i[0][0][0]]
                x += float(i[1]) * self.mass_dict_atom[i[0][0][0]]
                y += float(i[2]) * self.mass_dict_atom[i[0][0][0]]
                z += float(i[3]) * self.mass_dict_atom[i[0][0][0]]
            x /= m
            y /= m
            z /= m
            new_sequence_position.append([new_seq[0], x / 10, y / 10, z / 10, m])

        with open(os.path.join(self.xml_dir, file.replace(".log", ".xml")), 'w') as f:
            n_atoms = len(new_sequence_position)
            f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
            f.write('<galamost_xml version="1.3">\n')
            f.write(f'<configuration time_step="0" dimensions="3" natoms="{n_atoms}" >\n')
            f.write(f'<box lx="{self.box_size}" ly="{self.box_size}" lz=\"{self.box_size}"/>\n')

            # set position
            f.write(f'<position num="{n_atoms}">\n')
            for atom in new_sequence_position:
                f.write(f"{atom[1]:17.10f}{atom[2]:17.10f}{atom[3]:17.10f}\n")
                self.x_max = max(self.x_max, atom[1])
                self.y_max = max(self.y_max, atom[2])
                self.z_max = max(self.z_max, atom[3])
            f.write("</position>\n")

            # set image
            f.write('<image num="{}">\n'.format(n_atoms))
            for i in range(n_atoms):
                f.write("0 0 0\n")
            f.write("</image>\n")

            # set mass
            f.write('<mass num="{}">\n'.format(n_atoms))
            for atom in new_sequence_position:
                f.write("{:<17.10f}\n".format(atom[4]))
            f.write("</mass>\n")

            # set particle type
            f.write('<type num="{}">\n'.format(n_atoms))
            for atom in new_sequence_position:
                f.write("{}\n".format(self.pro_res_map[atom[0]][0]))
            f.write("</type>\n")

            # set bond
            f.write('<bond num="{}">\n'.format(n_atoms - 1))
            for i in range(0, n_atoms - 1):
                f.write("B-B {} {}\n".format(i, i + 1))
            f.write("</bond>\n")

            # set charge
            epsilon_r = 74.19
            k = 138.935  # 1 / (4 * pi * epsilon_0)
            q_unit = math.sqrt(k / epsilon_r)

            f.write('<charge num="{}">\n'.format(n_atoms))
            for atom in new_sequence_position:
                if self.pro_res_map[atom[0]][0] == 'R' or self.pro_res_map[atom[0]][0] == 'K':
                    if self.protein_model == "HPS":
                        f.write("{:<17.10f}\n".format(q_unit))
                    elif self.protein_model == "Mpipi":
                        f.write("{:<17.10f}\n".format(0.75 * q_unit))
                elif self.pro_res_map[atom[0]][0] == 'D' or self.pro_res_map[atom[0]][0] == 'E':
                    if self.protein_model == "HPS":
                        f.write("{:<17.10f}\n".format(-q_unit))
                    elif self.protein_model == "Mpipi":
                        f.write("{:<17.10f}\n".format(0.75 * -q_unit))
                elif self.pro_res_map[atom[0]][0] == 'H':
                    if self.protein_model == "HPS":
                        f.write("{:<17.10f}\n".format(1 / (1 + 10 ** (7.4 - 6)) * q_unit))
                    elif self.protein_model == "Mpipi":
                        f.write("{:<17.10f}\n".format(0.375 * q_unit))
                else:
                    f.write("{:<17.10f}\n".format(0))
            f.write("</charge>\n")

            # set molecule
            f.write('<molecule num="{}">\n'.format(n_atoms))
            for i in range(n_atoms):
                f.write("0\n")
            f.write("</molecule>\n")

            # set body
            f.write('<body num="{}">\n'.format(n_atoms))
            body_list = ['-1'] * n_atoms
            if self.add_rigid_body:
                # print("添加刚体！")
                for domain in self.rigid_domain_list:
                    domain = list(map(int, domain.split('-')))
                    body_list[domain[0] - 1:domain[1]] = [str(self.rigid_body_index)] * (domain[1] - domain[0] + 1)
                    self.rigid_body_index += 1
            f.write('\n'.join(body_list))
            f.write("</body>\n")

            f.write("</configuration>\n</galamost_xml>\n")


    def dna_log2xml_3SPN(self, file):
        tree = ET.parse(os.path.join(os.path.join(self.log_dir, file)))
        root = tree.getroot()

        seq_elem = root.find('sequence')
        sequence = seq_elem.text.split()
        positions_elem = root.find('position')
        positions = positions_elem.text.split('\n')

        sequence_position = []
        for i in sequence:
            sequence_position.append([i])
        for line in positions[1:-1]:
            line = line.split()
            sequence_position[int(line[-2]) - 1].append(line)
        # print(sequence_position)
        new_sequence = []
        for seq in sequence_position:
            # print(seq)
            new_sequence.append(["Ph"])
            new_sequence.append(["Su"])
            new_sequence.append([seq[0]])
            for i in seq[1:]:
                if "P" in i[0]:
                    new_sequence[-3].append(i)
                elif "'" in i[0]:
                    new_sequence[-2].append(i)
                else:
                    new_sequence[-1].append(i)

        new_sequence_position = []
        for new_seq in new_sequence[1:]:
            # print(len(new_seq), new_seq)
            # print(new_seq)
            m, x, y, z = 0, 0, 0, 0
            for i in new_seq[1:]:
                m += self.mass_dict_atom[i[0][0]]
                x += float(i[1]) * self.mass_dict_atom[i[0][0]]
                y += float(i[2]) * self.mass_dict_atom[i[0][0]]
                z += float(i[3]) * self.mass_dict_atom[i[0][0]]
            x /= m
            y /= m
            z /= m
            new_sequence_position.append([new_seq[0], x / 10, y / 10, z / 10, m])

        with open(os.path.join(self.xml_dir, file.replace(".log", ".xml")), 'w') as f:
            n_atoms = len(new_sequence_position)
            f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
            f.write('<galamost_xml version="1.3">\n')
            f.write(f'<configuration time_step="0" dimensions="3" natoms="{n_atoms}">\n')
            f.write(f'<box lx="{self.box_size}" ly="{self.box_size}" lz="{self.box_size}"/>\n')
            f.write('<position num="{}">\n'.format(n_atoms))
            for atom in new_sequence_position:
                f.write(f"{atom[1]:17.10f}{atom[2]:17.10f}{atom[3]:17.10f}\n")
                self.x_max = max(self.x_max, atom[1])
                self.y_max = max(self.y_max, atom[2])
                self.z_max = max(self.z_max, atom[3])
            f.write("</position>\n")

            f.write('<image num="{}">\n'.format(n_atoms))
            for i in range(n_atoms):
                f.write("0 0 0\n")
            f.write("</image>\n")

            # set particle type
            particles_map = {"DA": "Ab", "DT": "Tb", "DG": "Gb", "DC": "Cb", "Ph": "Ph", "Su": "Su"}
            f.write('<type num="{}">\n'.format(n_atoms))
            for atom in new_sequence_position:
                f.write("{}\n".format(particles_map[atom[0]]))
            f.write("</type>\n")

            # set mass
            mass_dict_particle = {"Ph": 94.97, "Su": 83.11, "DA": 134.1, "DT": 125.1, "DC": 110.1, "DG": 150.1}
            f.write('<mass num="{}">\n'.format(n_atoms))
            for atom in new_sequence_position:
                f.write("{:<17.10f}\n".format(mass_dict_particle[atom[0]]))
            f.write("</mass>\n")

            # set charge
            epsilon_r = 74.19
            k = 138.935  # 1 / (4 * pi * epsilon_0)
            q_unit = math.sqrt(k / epsilon_r)

            f.write('<charge num="{}">\n'.format(n_atoms))
            for atom in new_sequence_position:
                f.write(f"{-q_unit:<17.10f}\n") if atom[0] == 'DP' else f.write(f"{0.:<17.10f}\n")
            f.write('</charge>\n')

            # set molecule
            f.write('<molecule num="{}">\n'.format(n_atoms))
            for i in range(n_atoms):
                f.write("0\n")
            f.write("</molecule>\n")

            # set bond
            f.write('<bond num="{}">\n'.format(n_atoms - 1))
            for i in range(0, n_atoms, 3):
                # print(i)
                f.write(f"S-{particles_map[new_sequence_position[i + 1][0]]} {i} {i + 1}\n")
                if i + 1 != n_atoms - 1:
                    f.write("S5-P {} {}\n".format(i, i + 2))
                    f.write("S3-P {} {}\n".format(i + 2, i + 3))
            f.write("</bond>\n")

            # set angle
            f.write('<angle num="{}">\n'.format((4 * n_atoms - 11) // 3))
            for i in range(0, n_atoms, 3):
                # print(i)
                if i + 1 != n_atoms - 1:
                    f.write(f"P-5S-{particles_map[new_sequence_position[i + 1][0]]} {i + 1} {i} {i + 2}\n")
                    f.write(f"S5-P-3S {i} {i + 2} {i + 3}\n")
                    f.write(f"P-3S-{particles_map[new_sequence_position[i + 4][0]]} {i+2} {i+3} {i+4}\n")
                if i + 1 != n_atoms - 4 and i + 1 != n_atoms - 1:
                    f.write(f"P-5S3-P {i+2} {i+3} {i+5}\n")
            f.write("</angle>\n")

            # set dihedral
            f.write("<dihedral num=\"{}\">\n".format((4 * n_atoms - 14) // 3))
            for i in range(0, n_atoms, 3):
                if i + 1 != n_atoms - 1:
                    f.write(f"{particles_map[new_sequence_position[i + 1][0]]}-S3-P-5S {i+1} {i} {i+2} {i+3}\n")
                    f.write(f"S3-P-5S-{particles_map[new_sequence_position[i + 4][0]]} {i} {i+2} {i+3} {i+4}\n")
                if i + 1 != n_atoms - 4 and i + 1 != n_atoms - 1:
                    f.write(f"P-5S3-P-{particles_map[new_sequence_position[i + 4][0]]} {i} {i + 2} {i + 3} {i + 5}\n")
                    f.write(f"P-5S3-P-5S {i + 2} {i + 3} {i + 5} {i + 6}\n")
                if i + 1 != n_atoms - 7 and i + 1 != n_atoms - 4 and i + 1 != n_atoms - 1:
                    f.write(f"S3-P-5S3-P {i} {i + 2} {i + 3} {i + 5}\n")
            f.write("</dihedral>\n")

            f.write("</configuration>\n</galamost_xml>\n")


    def dna_log2xml_2BeadMittal(self, file):
        tree = ET.parse(os.path.join(self.log_dir, file))
        root = tree.getroot()

        seq_elem = root.find('sequence')
        sequence = seq_elem.text.split()
        positions_elem = root.find('position')
        positions = positions_elem.text.split('\n')

        sequence_position = []
        for i in sequence:
            sequence_position.append([i])
        for line in positions[1:-1]:
            line = line.split()
            sequence_position[int(line[-2]) - 1].append(line)
        # print(sequence_position)
        new_sequence = []
        for seq_idx in range(len(sequence_position)):
            # print(sequence_position[seq_idx])
            new_sequence.append(["Backbone"])
            new_sequence.append([sequence_position[seq_idx][0]])

            for i in sequence_position[seq_idx][1:]:
                if 'P' in i[0][0] or "'" in i[0][0]:  # 判断是否属于糖 - 磷酸主链原子
                    new_sequence[2 * seq_idx].append(i)
                else:
                    new_sequence[2 * seq_idx + 1].append(i)
            # print(new_sequence)
        new_sequence_position = []

        for new_seq in new_sequence:
            m, x, y, z = 0, 0, 0, 0
            for i in new_seq[1:]:
                m += self.mass_dict_atom[i[0][0][0]]
                x += float(i[1]) * self.mass_dict_atom[i[0][0][0]]
                y += float(i[2]) * self.mass_dict_atom[i[0][0][0]]
                z += float(i[3]) * self.mass_dict_atom[i[0][0][0]]
            x /= m
            y /= m
            z /= m
            new_sequence_position.append([new_seq[0], x / 10, y / 10, z / 10, m])
            # print(new_sequence_position)

        with open(os.path.join(self.xml_dir, file.replace(".log", ".xml")), 'w') as f:
            n_atoms = len(new_sequence_position)
            f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
            f.write('<galamost_xml version="1.3">\n')
            f.write(f'<configuration time_step="0" dimensions="3" natoms="{n_atoms}">\n')
            f.write(f'<box lx="{self.box_size}" ly="{self.box_size}" lz="{self.box_size}"/>\n')

            f.write('<position num="{}">\n'.format(n_atoms))
            for atom in new_sequence_position:
                f.write(f"{atom[1]:17.10f}{atom[2]:17.10f}{atom[3]:17.10f}\n")
                self.x_max = max(self.x_max, atom[1])
                self.y_max = max(self.y_max, atom[2])
                self.z_max = max(self.z_max, atom[3])
            f.write("</position>\n")

            f.write('<image num="{}">\n'.format(n_atoms))
            for i in range(n_atoms):
                f.write("0 0 0\n")
            f.write("</image>\n")

            # set mass
            mass_dict_particle = {"Backbone": 178.08, "DA": 134.1, "DT": 125.1, "DC": 110.1, "DG": 150.1}

            f.write('<mass num="{}">\n'.format(n_atoms))
            for atom in new_sequence_position:
                f.write(f"{mass_dict_particle[atom[0]]:<17.10f}\n")
            f.write("</mass>\n")

            # set particle type
            particles_map = {"DA": "Ab", "DT": "Tb", "DG": "Gb", "DC": "Cb", "Backbone": "PhSu"}
            f.write(f'<type num="{n_atoms}">\n')
            for atom in new_sequence_position:
                f.write(f"{particles_map[atom[0]]}\n")
            f.write("</type>\n")

            # set bond
            f.write('<bond num="{}">\n'.format(n_atoms - 1))
            for i in range(0, n_atoms, 2):
                f.write(f"PhSu-{particles_map[new_sequence_position[i+1][0]]} {i} {i+1}\n")
                if i + 2 < n_atoms:
                    f.write(f"PhSu-PhSu {i} {i+2}\n")
            f.write("</bond>\n")

            # set angle
            f.write('<angle num="{}">\n'.format(n_atoms // 2 - 4))
            for i in range(0, n_atoms, 2):
                if i != n_atoms - 2 and i != n_atoms - 4 and i != n_atoms // 2 - 2 and i != n_atoms // 2 - 4:
                    f.write(f"P-P-P {i} {i+2} {i+4}\n")
            f.write("</angle>\n")

            # set charge
            epsilon_r = 74.19
            k = 138.935  # 1 / (4 * pi * epsilon_0)
            q_unit = math.sqrt(k / epsilon_r)

            f.write('<charge num="{}">\n'.format(n_atoms))
            for atom in new_sequence_position:
                f.write(f"{-q_unit:<17.10f}\n") if atom[0] == 'Backbone' else f.write(f"{0.:<17.10f}\n")
            f.write('</charge>\n')

            # set molecule
            f.write('<molecule num="{}">\n'.format(n_atoms))
            for i in range(n_atoms):
                f.write("0\n")
            f.write("</molecule>\n")

            f.write("</configuration>\n</galamost_xml>\n")


    def rna_log2xml_3SPN(self, file):
        tree = ET.parse(os.path.join(os.path.join(self.log_dir, file)))
        root = tree.getroot()

        seq_elem = root.find('sequence')
        sequence = seq_elem.text.split()
        positions_elem = root.find('position')
        positions = positions_elem.text.split('\n')

        sequence_position = []
        for i in sequence:
            sequence_position.append([i])
        for line in positions[1:-1]:
            line = line.split()
            sequence_position[int(line[-2]) - 1].append(line)
        # print(sequence_position)
        new_sequence = []
        for seq in sequence_position:
            # print(seq)
            new_sequence.append(["Ph"])
            new_sequence.append(["Su"])
            new_sequence.append([seq[0]])
            for i in seq[1:]:
                if "P" in i[0]:
                    new_sequence[-3].append(i)
                elif "'" in i[0]:
                    new_sequence[-2].append(i)
                else:
                    new_sequence[-1].append(i)

        # print(new_sequence)

        new_sequence_position = []
        for new_seq in new_sequence[1:]:
            # print(len(new_seq), new_seq)
            # print(new_seq)
            m, x, y, z = 0, 0, 0, 0
            for i in new_seq[1:]:
                m += self.mass_dict_atom[i[0][0]]
                x += float(i[1]) * self.mass_dict_atom[i[0][0]]
                y += float(i[2]) * self.mass_dict_atom[i[0][0]]
                z += float(i[3]) * self.mass_dict_atom[i[0][0]]
            x /= m
            y /= m
            z /= m
            new_sequence_position.append([new_seq[0], x / 10, y / 10, z / 10, m])

        with open(os.path.join(self.xml_dir, file.replace(".log", ".xml")), 'w') as f:
            n_atoms = len(new_sequence_position)
            f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
            f.write('<galamost_xml version="1.3">\n')
            f.write(f'<configuration time_step="0" dimensions="3" natoms="{n_atoms}">\n')
            f.write(f'<box lx="{self.box_size}" ly="{self.box_size}" lz="{self.box_size}"/>\n')
            f.write('<position num="{}">\n'.format(n_atoms))
            for atom in new_sequence_position:
                f.write(f"{atom[1]:17.10f}{atom[2]:17.10f}{atom[3]:17.10f}\n")
                self.x_max = max(self.x_max, atom[1])
                self.y_max = max(self.y_max, atom[2])
                self.z_max = max(self.z_max, atom[3])
            f.write("</position>\n")

            f.write('<image num="{}">\n'.format(n_atoms))
            for i in range(n_atoms):
                f.write("0 0 0\n")
            f.write("</image>\n")

            # set particle type
            particles_map = {"Ph": "Ph", "Su": "Su",
                             "A": "Arb", "U": "Urb", "G": "Grb", "C": "Crb", }
            f.write('<type num="{}">\n'.format(n_atoms))
            for atom in new_sequence_position:
                f.write("{}\n".format(particles_map[atom[0]]))
            f.write("</type>\n")

            # set mass
            mass_dict_particle = {"Ph": 94.97, "Su": 83.11,
                                  "A": 134.1, "U": 125.1, "C": 110.1, "G": 150.1}
            f.write('<mass num="{}">\n'.format(n_atoms))
            for atom in new_sequence_position:
                f.write("{:<17.10f}\n".format(mass_dict_particle[atom[0]]))
            f.write("</mass>\n")

            # set charge
            epsilon_r = 74.19
            k = 138.935  # 1 / (4 * pi * epsilon_0)
            q_unit = math.sqrt(k / epsilon_r)

            f.write('<charge num="{}">\n'.format(n_atoms))
            for atom in new_sequence_position:
                f.write(f"{-q_unit:<17.10f}\n") if atom[0] == 'DP' else f.write(f"{0.:<17.10f}\n")
            f.write('</charge>\n')

            # set molecule
            f.write('<molecule num="{}">\n'.format(n_atoms))
            for i in range(n_atoms):
                f.write("0\n")
            f.write("</molecule>\n")

            # set bond
            f.write('<bond num="{}">\n'.format(n_atoms - 1))
            for i in range(0, n_atoms, 3):
                # print(i)
                f.write(f"S-{particles_map[new_sequence_position[i + 1][0]]} {i} {i + 1}\n")
                if i + 1 != n_atoms - 1:
                    f.write("S5-P {} {}\n".format(i, i + 2))
                    f.write("S3-P {} {}\n".format(i + 2, i + 3))
            f.write("</bond>\n")

            # set angle
            f.write('<angle num="{}">\n'.format((4 * n_atoms - 11) // 3))
            for i in range(0, n_atoms, 3):
                # print(i)
                if i + 1 != n_atoms - 1:
                    f.write(f"P-5S-{particles_map[new_sequence_position[i + 1][0]]} {i + 1} {i} {i + 2}\n")
                    f.write(f"S5-P-3S {i} {i + 2} {i + 3}\n")
                    f.write(f"P-3S-{particles_map[new_sequence_position[i + 4][0]]} {i + 2} {i + 3} {i + 4}\n")
                if i + 1 != n_atoms - 4 and i + 1 != n_atoms - 1:
                    f.write(f"P-5S3-P {i + 2} {i + 3} {i + 5}\n")
            f.write("</angle>\n")

            # set dihedral
            f.write("<dihedral num=\"{}\">\n".format((4 * n_atoms - 14) // 3))
            for i in range(0, n_atoms, 3):
                if i + 1 != n_atoms - 1:
                    f.write(f"{particles_map[new_sequence_position[i + 1][0]]}-S3-P-5S {i + 1} {i} {i + 2} {i + 3}\n")
                    f.write(f"S3-P-5S-{particles_map[new_sequence_position[i + 4][0]]} {i} {i + 2} {i + 3} {i + 4}\n")
                if i + 1 != n_atoms - 4 and i + 1 != n_atoms - 1:
                    f.write(f"P-5S3-P-{particles_map[new_sequence_position[i + 4][0]]} {i} {i + 2} {i + 3} {i + 5}\n")
                    f.write(f"P-5S3-P-5S {i + 2} {i + 3} {i + 5} {i + 6}\n")
                if i + 1 != n_atoms - 7 and i + 1 != n_atoms - 4 and i + 1 != n_atoms - 1:
                    f.write(f"S3-P-5S3-P {i} {i + 2} {i + 3} {i + 5}\n")
            f.write("</dihedral>\n")

            f.write("</configuration>\n</galamost_xml>\n")


    def generate_pygamd_run_file(self):
        file_path = os.path.join(self.path, f"run_{self.filename.split('.')[0]}.py")
        with open(file_path, 'w') as f:
            f.write("#!/usr/bin/python\n")
            f.write("import sys\n")
            f.write("import math\n\n")

            f.write("from poetry import cu_gala as gala\n")
            f.write("from poetry import force_field_gala\n")
            f.write("from poetry import _options\n\n")

            f.write(f"filename = \"{os.path.basename(self.output_file)}\"\n")
            f.write("build_method = gala.XMLReader(filename)\n")
            f.write("perform_config = gala.PerformConfig(_options.gpu)\n")
            f.write("all_info = gala.AllInfo(build_method, perform_config)\n\n")

            if self.add_enm_bond_flag:
                f.write("dt = 0.0002\n"); print("已设置弹性键，生成结构松弛模拟脚本，dt 首先设置为为 0.0002")
            else:
                f.write("dt = 0.02\n"); print("dt 默认设置为 0.02")
            f.write("app = gala.Application(all_info, dt)\n\n")


            f.write("Temperature = 310.0  # (k)\n"); print("温度设置为 310.0 K")
            Temperature = 310.0
            f.write("concent = 50.0  # salt concentration\n\n"); print("盐浓度设置为 50.0 mM")
            concent = 50.0

            f.write("kappaD = 13.603 * math.sqrt((50.0 / concent) * (Temperature / 300.0))  # related to salt concentration 13.603 -> 50mM NaCl\n")
            f.write("real_epsilon = 0.26 * 4184.0  # (0.26 kcal/mol = 0.26**4184.0 J/mol)\n")
            f.write("R = 8.314472  # gas constant\n")
            f.write("energy_reduce_unit = 1000.0\n")
            f.write("lenscale = 10.0  # angstrom\n")
            f.write("lsq = lenscale ** 2\n")
            f.write("epsilon = real_epsilon / energy_reduce_unit\n\n")
            kappaD = 13.603 * math.sqrt((50.0 / concent) * (Temperature / 300.0))
            print("kappaD = ", kappaD)
            real_epsilon = 0.26 * 4184.0
            R = 8.314472
            energy_reduce_unit = 1000.0
            lenscale = 10.0
            lsq = lenscale ** 2
            epsilon = real_epsilon / energy_reduce_unit
            print("epsilon = ", epsilon)

            f.write("rcut = 4.0\n\n"); print("短程相互作用 r_cut 默认为 4.0")

            f.write("neighbor_list = gala.NeighborList(all_info, rcut, 0.4)  # (,rcut,rbuffer)\n")
            f.write("neighbor_list.exclusion([\"bond\"])\n")
            f.write("neighbor_list.addExclusionsFromBonds()\n")
            if "dna" in self.mol_class:
                f.write("neighbor_list.addExclusionsFromAngles()\n")
                f.write("neighbor_list.addExclusionsFromDihedrals()\n")
                f.write("neighbor_list.setRCutPair(\"Ph\", \"Ph\", kappaD * 4.0 / lenscale)\n")
            f.write("neighbor_list.countExclusions()\n\n")

            f.write("ShortRangeEpsilon = 0.8368\n")
            f.write("debye_length = 0.794\n")
            f.write("ahdh = force_field_gala.AHDHForce(all_info, neighbor_list, rcut, ShortRangeEpsilon, debye_length, \"ahdh_DNA.force_field\", dna=True)\n")
            f.write("app.add(ahdh)\n\n")

            ahdh_params = {'A':[0.504, 0.0011162643859539204], 'R':[0.655999, 0.7249915947715212],
                           'N':[0.568, 0.43832729970272843], 'D':[0.558, 0.029182123776349763],
                           'C':[0.548, 0.610362354303913], 'Q':[0.602, 0.3268188050525212],
                           'E':[0.592, 0.006100281608649786], 'G':[0.45, 0.7012713677972457],
                           'H':[0.608, 0.46519480823469783], 'I':[0.618, 0.6075268330845265],
                           'L':[0.618, 0.5563020305733198], 'K':[0.636, 0.058617173158697924],
                           'M':[0.618, 0.7458993420826714], 'F':[0.636, 0.9216959832175945],
                           'P':[0.555999, 0.37296418535993475], 'S':[0.518, 0.46485701300656046],
                           'T':[0.562, 0.5379777613307019], 'W':[0.678, 0.9844235478393931],
                           'Y':[0.646, 0.9950108229594323], 'V':[0.586000, 0.41850068525598694],
                           "AD": [0.504, 0.729729729729730], "RD": [0.65599, 0.0],
                           "ND": [0.568, 0.432432432432432], "DD": [0.558, 0.378378378378378],
                           "CD": [0.548, 0.594594594594595], "QD": [0.602, 0.513513513513514],
                           "ED": [0.592, 0.459459459459459], "GD": [0.45, 0.648648648648649],
                           "HD": [0.608, 0.513513513513514], "ID": [0.618, 0.972972972972973],
                           "LD": [0.618, 0.972972972972973], "KD": [0.636, 0.513513513513514],
                           "MD": [0.618, 0.837837837837838], "FD": [0.636, 1.0],
                           "PD": [0.555999, 1.0], "SD": [0.518, 0.594594594594595],
                           "TD": [0.562, 0.675675675675676], "WD": [0.678, 0.945945945945946],
                           "YD": [0.646, 0.864864864864865], "VD": [0.586000, 0.891891891891892],
                           "Ph": [0.611, 0.459459], "Su": [0.611, 0.756757],
                           "Ab": [0.611, 0.351351], "Gb": [0.611, 0.540541],
                           "Cb": [0.611, 0.297297], "Tb": [0.611, 0.594595]}

            with open(os.path.join(self.path, "ahdh_DNA.force_field"), "w") as f_force_field:
                f_force_field.write("<ah_params>\n")
                if self.add_enm_bond_flag:
                    for key, value in list(ahdh_params.items())[:40]:
                        f_force_field.write(f"{key:5} {value[0]:<10} {value[1]:<30}\n")
                else:
                    for key, value in list(ahdh_params.items())[:20]:
                        f_force_field.write(f"{key:5} {value[0]:<10} {value[1]:<30}\n")
                if "dna" in self.mol_class:
                    for key, value in list(ahdh_params.items())[40:]:
                        f_force_field.write(f"{key:5} {value[0]:<10f} {value[1]:<30}\n")
                f_force_field.write("</ah_params>\n")

            f.write("bondforce_pro = gala.BondForceHarmonic(all_info)\n")
            f.write("bondforce_pro.setParams('B-B', 8033.28, 0.38)\n")
            if self.add_enm_bond_flag:
                enm_bond_file = f'{self.filename.replace(".pdb", "")}_enm_bond.py'
                f.write('with '); f.write(f'open("{enm_bond_file}", "r") as f:\n')
                f.write("    f_lines = f.readlines()\n")
                f.write("    for line in f_lines:\n")
                f.write("        eval(line.strip())\n")
            f.write("app.add(bondforce_pro)\n\n")

            f.write("T = Temperature * 8.3143 / 1000.0  # reduced unit\n\n")

            if self.add_rigid_body:
                f.write("group = gala.ParticleSet(all_info, \"all\")\n")
                f.write("comp_info = gala.ComputeInfo(all_info, group)\n\n")

                f.write("b_group = gala.ParticleSet(all_info, \"body\")\n")
                f.write("rigid_bd = gala.LangevinNVTRigid(all_info, b_group, T, 12345)\n\n")
                f.write("app.add(rigid_bd)\n\n")

                f.write("nb_group = gala.ParticleSet(all_info, \"non_body\")\n")
                f.write("non_rigid_bd = gala.LangevinNVT(all_info, nb_group, T, 12345)\n\n")

                gamma_dict = {"R": 1.5619, "D": 1.1509, "N": 1.141, "E": 1.2911, "K": 1.2817, "H": 1.3714,
                              "Q": 1.2813, "S": 0.8708, "C": 1.0314, "G": 0.5705, "T": 1.0111, "A": 0.7107,
                              "M": 1.312, "Y": 1.6318, "V": 0.9913, "W": 1.8622, "L": 1.1316, "I": 1.1316,
                              "P": 0.9712, "F": 1.4718}
                for res in gamma_dict.keys():
                    f.write(f"non_rigid_bd.setGamma(\"{res}\", {gamma_dict[res]})\n")

                f.write("app.add(non_rigid_bd)\n\n")
            else:
                f.write("groupall = gala.ParticleSet(all_info, \"all\")\n")
                f.write("comp_info = gala.ComputeInfo(all_info, groupall)\n\n")
                f.write("gp = gala.LangevinNVT(all_info, groupall, T, 12345)\n")
                f.write("app.add(gp)\n\n")

            f.write("sort_method = gala.Sort(all_info)\n")
            f.write("sort_method.setPeriod(500)\n")
            f.write("app.add(sort_method)\n\n")
            f.write("ZeroMomentum = gala.ZeroMomentum(all_info)\n")
            f.write("ZeroMomentum.setPeriod(2000)\n")
            f.write("app.add(ZeroMomentum)\n\n")
            f.write("DInfo = gala.DumpInfo(all_info, comp_info, 'data.log')\n")
            f.write("DInfo.setPeriod(int(1000))\n")
            f.write("app.add(DInfo)\n\n")

            f.write("xml = gala.XMLDump(all_info, 'particles')\n")
            f.write(f"xml.setPeriod(int(5e5))  # (period)\n"); print("输出周期设置为 5e5")
            f.write("xml.setOutputBond(True)\n")
            f.write("xml.setOutputAngle(True)\n")
            f.write("xml.setOutputDihedral(True)\n")
            f.write("xml.setOutputVelocity(True)\n")
            f.write("xml.setOutputMass(True)\n")
            f.write("xml.setOutputCharge(True)\n")
            f.write("xml.setOutputBody(True)\n")
            f.write("app.add(xml)\n\n")

            if self.add_enm_bond_flag:
                f.write("app.run(int(100000))  # (How many steps to run)\n")
                f.write("app.setDt(0.001)\n")
                f.write("app.run(int(100000))\n")
                f.write("app.setDt(0.02)\n")
                f.write("app.run(int(15000000))\n")
            else:
                f.write("app.run(int(1e9))\n"); print("模拟时间设置为 1e9")