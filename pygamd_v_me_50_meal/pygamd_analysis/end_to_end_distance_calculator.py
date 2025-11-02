import os
import re
import shutil

from tqdm import tqdm

import numpy as np
import matplotlib.pyplot as plt

class EndToEndDistanceCalculator:
    def __init__(self, path, data):
        self.path = path
        self.data = data

        self.mol_class_dict = self.data.mol_class_dict
        self.length_dict = self.data.length_dict

        self.chain_path = os.path.join(self.path, "chain_xyz_unwrapping/")

        self.save_path = os.path.join(self.path, "draw_log/")
        os.makedirs(self.save_path, exist_ok=True)

        self.balance_cut = ""
        self.end_to_end_distance_path = os.path.join(self.save_path, "end_to_end_distance/")

        self.cal_end_to_end_distance_list = []
        if not self.cal_end_to_end_distance_list:
            print(f"\n您的分子类型有：\n{self.data.molecules}")
            self.cal_end_to_end_distance_list = input("请输入您想要计算末端距时需要包括的分子，以逗号分隔，如“1,2”。\n"
                                                      "如需计算全部分子，请输入 all 或直接回车：").split(',')
            if "all" in self.cal_end_to_end_distance_list or self.cal_end_to_end_distance_list == [""]:
                self.cal_end_to_end_distance_list = list(self.data.mol_class_dict.keys())
            else:
                self.cal_end_to_end_distance_list = [self.data.mol_class_list[int(i)] for i in self.cal_end_to_end_distance_list]
        print(f"即将计算末端距的分子：{self.cal_end_to_end_distance_list}")

        self.files = sorted(os.listdir(self.chain_path))
        if not self.balance_cut:
            self.balance_cut = input(
                "请输入需要截取的平衡后的文件索引，格式为‘开始,结束’，例如：1000,2000，直接回车则不截取：")
            if not self.balance_cut:
                self.files = sorted(os.listdir(self.chain_path))
            else:
                start, end = list(map(int, self.balance_cut.split(',')))
                self.files = sorted(os.listdir(self.chain_path))[start: end + 1]

        self.end_to_end_distance_dict = dict(zip(self.cal_end_to_end_distance_list,
                                                 [[] for _ in range(len(self.cal_end_to_end_distance_list))]))

    def cal_end_to_end_distance(self, filename):
        with open(os.path.join(self.chain_path, filename), 'rb') as f:
            data = eval(f.read())

        # match = re.search(r'\d+', filename)
        # if match:
        #     time_step = int(match.group())
        #
        # time_step = None
        for mol_name in self.cal_end_to_end_distance_list:
            mol_data = data[mol_name]
            end_to_end_distance = []
            for chain_idx in range(len(mol_data)):
                end_to_end_distance_chain = []
                for i in range(len(mol_data[chain_idx]) - 1):
                    d = np.linalg.norm(np.array(mol_data[chain_idx][i+1]) - np.array(mol_data[chain_idx][0]))
                    end_to_end_distance_chain.append(d)
                end_to_end_distance.append(end_to_end_distance_chain)
            self.end_to_end_distance_dict[mol_name].append(np.array(end_to_end_distance))

    def cal_end_to_end_distance_parallel(self):
        # 计算 EED
        if not os.path.exists(self.end_to_end_distance_path):
            os.makedirs(self.end_to_end_distance_path, exist_ok=True)
        else:
            shutil.rmtree(self.end_to_end_distance_path)
            os.makedirs(self.end_to_end_distance_path, exist_ok=True)

        # 计算 EED
        print("✅ 开始计算末端距")
        list(tqdm(map(self.cal_end_to_end_distance, self.files),
                       total=len(self.files),
                       desc="计算末端距",
                       colour='cyan',
                       bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]',
                       ncols=100))

        # 时间平均
        for mol_name in self.cal_end_to_end_distance_list:
            end_to_end_distance = self.end_to_end_distance_dict[mol_name]
            end_to_end_distance_avg = np.mean(end_to_end_distance, axis=0)
            self.end_to_end_distance_dict[mol_name] = end_to_end_distance_avg

        # 储存
        for mol_name in self.cal_end_to_end_distance_list:
            end_to_end_distance = self.end_to_end_distance_dict[mol_name]
            with open(os.path.join(self.end_to_end_distance_path, f"{mol_name}.txt"), 'wb') as f:
                f.write(str(end_to_end_distance).encode('utf-8'))

        # 绘图
        fig, ax = plt.subplots(figsize=(12, 9), dpi=300)
        fig_ln, ax_ln = plt.subplots(figsize=(12, 9), dpi=300)
        for mol_name in self.cal_end_to_end_distance_list:
            end_to_end_distance = self.end_to_end_distance_dict[mol_name]
            for chain_idx in range(len(end_to_end_distance)):
                end_to_end_distance_chain = end_to_end_distance[chain_idx]
                residue = np.arange(len(end_to_end_distance_chain))  # ns
                ax.plot(residue, end_to_end_distance_chain, label=f"{mol_name}_{chain_idx+1}")

                ln_residue = np.log(residue)
                ln_end_to_end_distance_chain = np.log(end_to_end_distance_chain)
                ax_ln.plot(ln_residue, ln_end_to_end_distance_chain, label=f"{mol_name}_{chain_idx+1}_ln")
            ax.set_xlabel('N')
            ax.set_ylabel('End-to-end distance (nm)')
            ax.set_title(f'End-to-End Distance of {mol_name}')
            ax.legend()
            fig.savefig(os.path.join(self.end_to_end_distance_path, f"draw_end_to_end_distance_{mol_name}.png"))

            ax_ln.set_xlabel('ln(N)')
            ax_ln.set_ylabel('ln(End-to-end distance)')
            ax_ln.set_title(f'End-to-End Distance of {mol_name}')
            ax_ln.legend()
            fig_ln.savefig(os.path.join(self.end_to_end_distance_path, f"draw_end_to_end_distance_ln_{mol_name}.png"))

        print("✅ 计算完成")