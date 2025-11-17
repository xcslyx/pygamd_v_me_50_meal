import os

import xml.etree.ElementTree as ET

from pygamd_v_me_50_meal.pygamd_analysis.get_sequence import GetSequence


class XMLConverter:
    def __init__(self, data, path, xml_file):
        self.data = data
        self.path = path
        self.xml_file = os.path.join(path, xml_file)

        # 残基名称映射表
        self.residue_name_mapping = {
                                    'A': 'ALA', 'F': 'PHE', 'C': 'CYS', 'U': 'SEC', 'D': 'ASP',
                                    'N': 'ASN', 'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS',
                                    'L': 'LEU', 'I': 'ILE', 'K': 'LYS', 'O': 'PYL', 'M': 'MET',
                                    'P': 'PRO', 'R': 'ARG', 'S': 'SER', 'T': 'THR', 'V': 'VAL',
                                    'W': 'TRP', 'Y': 'TYR'
                                    }

        self.sequence = GetSequence(path, xml_file, data)


    @staticmethod
    def extract_position_and_type(xml_file_path):
        """
        解析XML文件，提取position和type元素中的信息。
        """
        tree = ET.parse(xml_file_path)
        root = tree.getroot()

        position_list = []
        for position in root.findall('.//position'):
            positions_for_current_element = [
                tuple(map(float, line.split()))
                for line in position.text.strip().split('\n')
            ]
            position_list.append(positions_for_current_element)

        type_list = []
        for type_element in root.findall('.//type'):
            type_list.extend(line.strip() for line in type_element.text.strip().split('\n'))

        return position_list[0], type_list


    def xml2pdb(self, output_file=None):
        """
        将XML文件转换为PDB文件。每 26 条链一个 PDB
        """
        chain_num = sum([i[0] for i in self.data.mol_class_dict.values()])
        chain_count = [26] * (chain_num // 26) + [chain_num % 26] if chain_num > 26 else [chain_num]

        positions, types = self.extract_position_and_type(self.xml_file)
        output_file = os.path.join(self.path, output_file) if output_file else self.xml_file.replace('.xml', '.pdb')
        """
        将信息写入PDB文件。
        """
        with open(output_file, 'w', encoding='utf-8') as pdb_file:
            atom_serial = 0
            residue_serial = 0

            # 生成每个链
            chain_names = [chr(65 + i % 26) for i in range(chain_count)]  # A, B, ..., Z, A, B, ...
            for chain_index, chain_name in enumerate(chain_names):
                chain_start_index = start_index + chain_index * chain_length
                chain_end_index = chain_start_index + chain_length

                for i in range(chain_start_index, chain_end_index):
                    x, y, z = positions[i]
                    atom_name = "CA"
                    residue_name = self.residue_name_mapping.get(types[i][0], "UNK")
                    occupancy = 1.00
                    b_factor = 1.00
                    element = "C"

                    # atom_serial 重置逻辑，限制最大值为 99999，超过则重置为 0
                    atom_serial = atom_serial % 100000  # 如果 atom_serial 达到 100000，则从 0 开始

                    # residue_serial 重置逻辑，限制最大值为 9999，超过则重置为 0
                    residue_serial = residue_serial % 10000  # 如果 atom_serial 达到 10000，则从 0 开始

                    pdb_file.write(
                        f"ATOM  {atom_serial:5d} {atom_name:^4} {residue_name:<3} {chain_name}{residue_serial:4d}    "
                        f"{x * 10:8.3f}{y * 10:8.3f}{z * 10:8.3f}{occupancy:6.2f}{b_factor:6.2f}          {element:>2}\n"
                    )
                    atom_serial += 1
                    residue_serial += 1

                residue_serial = 0  # 每条链的残基计数重置
            pdb_file.write("TER\n")  # 链结束
            pdb_file.write("END\n")  # 文件结束
