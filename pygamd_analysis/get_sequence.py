import os
import xml.etree.ElementTree as ET

from data import Data


# 定义类 GetSequence
class GetSequence:
    def __init__(self, path: str, filename: str, data: Data, output: str=None):
        self.filename = os.path.join(path, filename)
        self.data = data
        self.output = os.path.join(path, output) if output else self.filename.replace('.xml', '.txt')

    def pdb2sequence(self):
        with open(self.filename, 'r') as pdb:
            pdb_lines = pdb.readlines()
            atoms = []
            for line in pdb_lines:
                if "ATOM" in line[:4]:
                    atoms.append(line.replace('\n', '').split())

        pro_res_map = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                       'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                       'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                       'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
                       'HIE': 'H'}
        sequence = []

        for atom_num in range(len(atoms)):
            try:
                if atoms[atom_num][5] != atoms[atom_num + 1][5]:
                    sequence.append(pro_res_map[atoms[atom_num][3]])
            except IndexError:
                pass
        sequence.append(pro_res_map[atoms[-1][3]])
        with open(self.filename.replace(".pdb", ".txt"), "w") as log:
            mol_name = input("请输入分子名称：")
            log.write(f">{mol_name}\n" + ' '.join(sequence))

    def xml2sequence(self, mass: bool=True):
        mol_class_dict = self.data.mol_class_dict

        tree = ET.parse(self.filename)
        root = tree.getroot()
        for child in root.iter():
            if child.tag == "type":
                type_elem = child
                sequence = type_elem.text.strip().split('\n')
            if child.tag == "mass" and mass:
                mass_elem = child.text.strip().split('\n')
                while '' in mass_elem:
                    mass_elem.remove('')
                mass_list = list(map(float, mass_elem))


        with open(self.output, "w") as log:
            seq_dict = {}
            for mol in mol_class_dict.keys():
                if mass:
                    cur_sequence = [[i, j] for i, j in zip(sequence[:mol_class_dict[mol][1]], mass_list[:mol_class_dict[mol][1]])]
                else:
                    cur_sequence = [[i] for i in sequence[:mol_class_dict[mol][1]]]
                seq_dict[mol] = cur_sequence
                sequence = sequence[mol_class_dict[mol][0] * mol_class_dict[mol][1]:]
                if mass:
                    mass_list = mass_list[mol_class_dict[mol][0] * mol_class_dict[mol][1]:]
            log.write(str(seq_dict))
