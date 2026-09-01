import os

import numpy as np
import pygamd_v_me_50_meal.utils as utils

from tqdm import tqdm
from multiprocessing import Pool
from pygamd_v_me_50_meal.file_processor.xml_data_extractor import XMLDataExtractor

class BondAnalysys:
    def __init__(self, path: str, data, lang: str="zh"):
        self.path, self.data, self.lang = path, data, lang

        self.xml_path = os.path.join(self.path, "xml_unwrapping/")

        self.sys_topol_path = os.path.join(self.path, "sys_topol/")
        os.makedirs(self.sys_topol_path, exist_ok=True)

        self.bond_path = os.path.join(self.sys_topol_path, "bond/")
        os.makedirs(self.bond_path, exist_ok=True)


    def get_bond_length(self, file_name: str):
        position_list: np.ndarray = XMLDataExtractor(self.xml_path + file_name).extract_position_data()
        _, bond_dict = XMLDataExtractor(self.xml_path + file_name).extract_bond_data()

        bond_length_dict = {}

        for bond in bond_dict:
            bond_length_dict[bond] = []
            for bond_pair in bond_dict[bond]:
                bond_length_dict[bond].append(np.linalg.norm(position_list[int(bond_pair[0])] - position_list[int(bond_pair[1])]))
        
        with open(os.path.join(self.bond_path, file_name), 'w') as f:
            f.write(str(bond_length_dict))

        
    def get_bond_length_parallel(self):
        files = os.listdir(self.xml_path)
        
        with Pool(processes=4) as pool:
                # 使用 tqdm 包装可迭代对象
                list(tqdm(pool.imap(self.get_bond_length, files),
                          total=len(files),
                          desc="Calculating bond length",
                          colour='cyan',
                          bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]',
                          ncols=100))
    
    
