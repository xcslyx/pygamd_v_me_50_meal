import os

import numpy as np
import pygamd_v_me_50_meal.utils as utils

from tqdm import tqdm
from multiprocessing import Pool

from pygamd_v_me_50_meal.Functions import Functions
from pygamd_v_me_50_meal.file_processor.xml_data_extractor import XMLDataExtractor

class DihedralAnalysys:
    def __init__(self, path: str, data, lang: str="zh"):
        self.path, self.data, self.lang = path, data, lang

        self.xml_path = os.path.join(self.path, "xml_unwrapping/")

        self.sys_topol_path = os.path.join(self.path, "sys_topol/")
        os.makedirs(self.sys_topol_path, exist_ok=True)

        self.dihedral_path = os.path.join(self.sys_topol_path, "dihedral/")
        os.makedirs(self.dihedral_path, exist_ok=True)


    def get_dihedral_degree(self, file_name: str):
        position_list: np.ndarray = XMLDataExtractor(self.xml_path + file_name).extract_position_data()
        _, dihedral_dict = XMLDataExtractor(self.xml_path + file_name).extract_dihedral_data()

        dihedral_dict = {}

        for dihedral in dihedral_dict:
            dihedral_dict[dihedral] = []
            for dihedral in dihedral_dict[dihedral]:
                dihedral_dict[dihedral].append(Functions.compute_dihedral(position_list[dihedral[0]], position_list[dihedral[1]], position_list[dihedral[2]], position_list[dihedral[3]]))
        
        with open(os.path.join(self.dihedral_path, file_name), 'w') as f:
            f.write(str(dihedral_dict))
        
    def get_dihedral_degree_parallel(self):
        files = os.listdir(self.xml_path)
        
        with Pool(processes=4) as pool:
                # 使用 tqdm 包装可迭代对象
                list(tqdm(pool.imap(self.get_dihedral_degree, files),
                          total=len(files),
                          desc="Calculating dihedral degree",
                          colour='cyan',
                          bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]',
                          ncols=100))
        
        dihedral_dict = {}
        for file_name in os.listdir(self.dihedral_path):
            with open(os.path.join(self.dihedral_path, file_name), 'r') as f:
                for dihedral, values in eval(f.read()).items():
                    if dihedral not in dihedral_dict:
                        dihedral_dict[dihedral] = []
                    dihedral_dict[dihedral].extend(values)
        
        with open(os.path.join(self.sys_topol_path, "dihedral_degree.txt"), 'w') as f:
            f.write(str(dihedral_dict))
    
    
