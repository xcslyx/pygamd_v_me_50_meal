import os

import numpy as np
import pygamd_v_me_50_meal.utils as utils

from tqdm import tqdm
from multiprocessing import Pool

from pygamd_v_me_50_meal.Functions import Functions
from pygamd_v_me_50_meal.file_processor.xml_data_extractor import XMLDataExtractor

class AngleAnalysys:
    def __init__(self, path: str, data, lang: str="zh"):
        self.path, self.data, self.lang = path, data, lang

        self.xml_path = os.path.join(self.path, "xml_unwrapping/")

        self.sys_topol_path = os.path.join(self.path, "sys_topol/")
        os.makedirs(self.sys_topol_path, exist_ok=True)

        self.angle_path = os.path.join(self.sys_topol_path, "angle/")
        os.makedirs(self.angle_path, exist_ok=True)


    def get_angle_degree(self, file_name: str):
        position_list: np.ndarray = XMLDataExtractor(self.xml_path + file_name).extract_position_data()
        _, angle_dict = XMLDataExtractor(self.xml_path + file_name).extract_angle_data()

        angle_dict = {}

        for angle in angle_dict:
            angle_dict[angle] = []
            for angle_index in angle_dict[angle]:
                angle_dict[angle].append(Functions.compute_angle(position_list[angle_index[0]], position_list[angle_index[1]], position_list[angle_index[2]]))
        
        with open(os.path.join(self.angle_path, file_name), 'w') as f:
            f.write(str(angle_dict))
        
    def get_angle_degree_parallel(self):
        files = os.listdir(self.xml_path)
        
        with Pool(processes=4) as pool:
                # 使用 tqdm 包装可迭代对象
                list(tqdm(pool.imap(self.get_angle_degree, files),
                          total=len(files),
                          desc="Calculating angle degree",
                          colour='cyan',
                          bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]',
                          ncols=100))
        
        angle_dict = {}
        for file_name in os.listdir(self.angle_path):
            with open(os.path.join(self.angle_path, file_name), 'r') as f:
                for angle, values in eval(f.read()).items():
                    if angle not in angle_dict:
                        angle_dict[angle] = []
                    angle_dict[angle].extend(values)
        
        with open(os.path.join(self.sys_topol_path, "angle_degree.txt"), 'w') as f:
            f.write(str(angle_dict))
    
    
