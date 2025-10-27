import os
import re

class Data:
    _instance = None
    mol_class_dict = {}
    length_dict = {
        "cGAS": 522, "mcGAS": 507,
        "Ga": 466, "Gb": 466, "Gaw": 466, "Gbw": 466,
        "FC": 528,
        "fuslcd": 162,
        "MED1": 671, "Na": 1, "Cl": 1
    }
    molecules = None

    system_name = ""

    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = super(Data, cls).__new__(cls)
        return cls._instance

    def __init__(self, path):
        if not hasattr(self, "_initialized"):
            self._initialized = True
            self.path = path
            path_list = self.path.split("/")
            self.system_name = path_list[-1] if path_list[-1] != "" else path_list[-2]

            matches = re.findall(r'(\d+)([a-zA-Z\d]+)(?:-(\d+))?', self.system_name)
            if matches:
                print("提取到以下分子信息：")
                start = 0
                for num, type_, length in matches:
                    self.mol_class_dict[type_] = [int(num)]
                    if not length:
                        if type_ not in self.length_dict:
                            self.length_dict[type_] = int(input(f"未知的分子类型: {type_}，请手动输入长度："))
                    else:
                        self.length_dict[type_] = int(length)
                    self.mol_class_dict[type_].append(self.length_dict[type_])
                    self.mol_class_dict[type_].append([start, start + int(num) * self.length_dict[type_]])
                    start += int(num) * self.length_dict[type_]
                    print(f"类型: {type_}，数量: {num}，长度: {self.mol_class_dict[type_][1]}，起止位置：{self.mol_class_dict[type_][2]}")

                    self.mol_class_list = list(self.mol_class_dict.keys())
                    self.molecules = "\n".join([f"{i}: {name}" for i, name in enumerate(self.mol_class_dict.keys())])
            else:
                raise ValueError("路径格式不正确")




