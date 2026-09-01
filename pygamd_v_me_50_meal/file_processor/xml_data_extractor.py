import numpy as np
import xml.etree.ElementTree as ET

class XMLDataExtractor:
    def __init__(self, xml_file_path):
        self.xml_file_path = xml_file_path
        self.tree = ET.parse(self.xml_file_path)
        self.root = self.tree.getroot()
    

    def get_box_size(self):
        box_size: list[float] = [float(self.root.find('.//box').attrib[i]) for i in ['lx', 'ly', 'lz']]
        print(f"Your box size: a = {box_size[0]}, b = {box_size[1]}, c = {box_size[2]}.")
        return box_size
    

    def extract_position_data(self) -> np.ndarray:
        position_elem = self.root.findall(f".//position")[0]
        position_list = np.array([list(map(float, line.split())) for line in position_elem.text.splitlines()[1:]])
        return position_list
    
    def extract_bond_data(self) -> tuple:
        bond_elem = self.root.findall(f".//bond")[0]
        bond_list = [[line.split()[0], [line.split()[1], line.split()[2]]] for line in bond_elem.text.splitlines()[1:]]
        bond_dict = {}
        for bond in bond_list:
            if bond[0] not in bond_dict:
                bond_dict[bond[0]] = [bond[1]]
            bond_dict[bond[0]].append(bond[1])
        return bond_list, bond_dict
    