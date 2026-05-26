import os
import re
import json

class Data:
    _instance = None
    mol_class_dict = {}
    length_dict = {
        "cGAS": 522, "hcGAS": 522, "mcGAS": 507,
        "Ga": 466, "Gb": 466, "Gaw": 466, "Gbw": 466,
        "FC": 528,
        "fuslcd": 162,
        "MED1": 627, "Na": 1, "Cl": 1
    }
    molecules = None

    system_name = ""

    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = super(Data, cls).__new__(cls)
        return cls._instance

    def __init__(self, path, lang='zh'):
        if not hasattr(self, "_initialized"):
            self._initialized = True

            # 加载消息文件
            with open(os.path.join(os.path.dirname(__file__), 'message.json'), 'r', encoding='utf-8') as f:
                messages = json.load(f)
                msg = messages['data_message']

            self.path = path
            path_list = self.path.split("/")
            self.system_name = path_list[-1] if path_list[-1] != "" else path_list[-2]

            matches = re.findall(r'(\d+)([a-zA-Z\d]+)(?:-(\d+))?', self.system_name)
            if matches:
                start = 0
                for num, type_, length in matches:
                    self.mol_class_dict[type_] = [int(num)]
                    if not length:
                        if type_ not in self.length_dict:
                            self.length_dict[type_] = int(input(f"{msg['UnkonwnMoleculeType'][lang]}{type_}, {msg['InputLength'][lang]}: "))
                    else:
                        self.length_dict[type_] = int(length)
                    self.mol_class_dict[type_].append(self.length_dict[type_])
                    self.mol_class_dict[type_].append([start, start + int(num) * self.length_dict[type_]])
                    start += int(num) * self.length_dict[type_]

                    self.mol_class_list = list(self.mol_class_dict.keys())
                    self.molecules = "\n".join([f"{i+1}: {name}" for i, name in enumerate(self.mol_class_dict.keys())])

                self.particle_num = 0
                print("The following molecular information was extracted:")
                for mol in self.mol_class_dict:
                    print(f"Molecule: {mol}, number: {self.mol_class_dict[mol][0]}, length: {self.mol_class_dict[mol][1]}")
                    self.particle_num += self.mol_class_dict[mol][0] * self.mol_class_dict[mol][1]
                print(f"Total number of particles: {self.particle_num}")
            else:
                raise ValueError("Path name is not valid.")




