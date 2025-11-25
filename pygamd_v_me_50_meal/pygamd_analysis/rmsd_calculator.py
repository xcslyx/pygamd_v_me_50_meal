import os
import math

import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm
from scipy.ndimage import gaussian_filter

from pygamd_v_me_50_meal.data import Data
from pygamd_v_me_50_meal.Functions import Functions


# Define a class to calculate RMSD.
class RMSDCalculator:
    def __init__(self, path, data: Data, ref=None):
        """
        用于计算 RMSD.

        :param path: 路径
        :param data: Data 类实例
        :param ref: 参考结构的 XML 文件名
        """
        self.path = path
        self.data = data

        self.chain_path = os.path.join(self.path, "chain_xyz_unwrapping/")

        self.save_path = os.path.join(self.path, "draw_log/")
        os.makedirs(self.save_path, exist_ok=True)

        self.ref = input("请指定含有单个参考结构的 XML 文件名：").replace(".xml", "") if not ref else ref
        self.init_pos = self.get_init_pos(self.ref)  # 读取参考结构的初始位置


        self.cal_class_rmsd = []
        self.rmsd_results = {}
        self.cur_chain_class = ""
        self.domain = []

        self.balance_cut = None


    @staticmethod
    def get_init_pos(filename):
        init_pos = []
        with open(filename + ".xml", 'r') as f:
            f_lines = f.readlines()
            position_flag = 0
            for line in f_lines:
                if "<position" in line:
                    position_flag = 1
                    continue
                elif "</position" in line:
                    break

                if position_flag:
                    init_pos.append(list(map(float, line.strip('\n').split())))
        return init_pos


    def cal_rmsd(self, chain_xyz):
        if self.domain:
            init_chain_xyz = np.array(self.init_pos[self.domain[0]-1:self.domain[1]])
        else:
            init_chain_xyz = np.array(self.init_pos)
        chain_xyz = np.array(chain_xyz)

        aa_num = len(init_chain_xyz)
        if aa_num != len(chain_xyz):
            exit("Error! Length of chain is not equal to that of initial chain.!")

        init_chain_xyz, aligned_chain_xyz = Functions.kabsch_align(init_chain_xyz, chain_xyz)

        rmsd = 0
        for i in range(aa_num):
            rmsd += sum(sum([(aligned_chain_xyz[i][j] - init_chain_xyz[i][j]) ** 2]) for j in range(3))
        return math.sqrt(rmsd / aa_num)


    def process_chain_file(self, name):
        with open(os.path.join(self.chain_path, name), 'r') as fp:
            x_mat: dict = eval(fp.read())

        for i in x_mat.keys():
            if i not in self.cal_class_rmsd:
                continue
            else:
                if self.domain:
                    cur_chain_xyz = map(lambda x: x[self.domain[0]-1:self.domain[1]], x_mat[i])
                else:
                    cur_chain_xyz = x_mat[i]
            if i in self.cal_class_rmsd:
                self.rmsd_results[i] = [] if i not in self.rmsd_results else self.rmsd_results[i]
                cur_rmsd_results = list(map(self.cal_rmsd, cur_chain_xyz))
                self.rmsd_results[i] += cur_rmsd_results


    def calculate(self):
        self.cal_class_rmsd = list(
            map(lambda x: int(x) - 1, input(f"请输入想要计算 RMSD 的分子序号：\n{self.data.molecules}\n").split(',')))
        self.cal_class_rmsd = [self.data.mol_class_list[i] for i in self.cal_class_rmsd]
        print(f"即将计算 RMSD 的分子：{self.cal_class_rmsd}")
        self.rmsd_results = {}

        domain = input(
            "若只需要计算结构域，请输入该结构域的起始残基编号和末尾残基编号（从 1 开始），以-分隔，如 159-522，若有多个结构域，请以英文逗号分隔。\n"
            "否则请直接回车：")
        if domain:
            domain = list(map(int, domain.split('-')))
            self.domain = domain
            print(f"即将计算结构域：{domain}")
        else:
            self.domain = None

        if not self.balance_cut:
            self.balance_cut = input(
                "请输入需要截取的平衡后的文件索引，索引从 1 开始，格式为 'START-END', 例如：1000-2000, 直接回车则不截取：")
        if not self.balance_cut:
            chain_files = os.listdir(self.chain_path)
        else:
            start, end = list(map(int, self.balance_cut.split('-')))
            chain_files = os.listdir(self.chain_path)[start - 1: end]

        # 使用 tqdm 包装可迭代对象以显示进度条
        list(tqdm(map(self.process_chain_file, sorted(chain_files)),
                  total=len(chain_files),
                  desc="Calculating RMSD",
                  colour='cyan',
                  bar_format = '{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]',
                  ncols=100))

        # 保存结果
        result_file = os.path.join(self.save_path, f"draw_RMSD_ref_{os.path.basename(self.ref)}.log")
        with open(result_file, 'w') as f:
            f.write(str(self.rmsd_results))
        print(f"RMSD 计算完成！结果已保存至文件 draw_RMSD_ref_{os.path.basename(self.ref)}.log")

        self.draw_rmsd_distribution()


    def draw_rmsd_distribution(self):
        if not self.ref:
            self.ref = input("请指定含有单个参考结构的 XML 文件名：")
        self.ref = self.ref.replace(".xml", "")

        with open(os.path.join(self.save_path, f"draw_RMSD_ref_{os.path.basename(self.ref)}.log"), 'r') as f:
            rmsd_results = eval(f.read())  # 读取 RMSD 结果

        for mol in rmsd_results.keys():
            fig, ax = plt.subplots(figsize=(12, 9), dpi=300)
            init_rmsd_list = rmsd_results[mol]
            rmsd = np.array(sorted(init_rmsd_list))
            frame_num = int(len(init_rmsd_list) / self.data.mol_class_dict[mol][0])
            print(f"{mol} 的帧数：{frame_num}")

            # 划分 bin
            bins = 40
            hist, bin_edges = np.histogram(rmsd, bins=bins)
            # hist = gaussian_filter(hist, sigma=5)
            # 计算概率
            probabilities = hist / sum(hist) / (bin_edges[1] - bin_edges[0])

            bin_edges += (bin_edges[1] - bin_edges[0]) / 2  # 使 bin 居中
            ax.plot(bin_edges[:-1], probabilities, label=f"{mol} RMSD")
            ax.hist(init_rmsd_list, bins=bins, density=True, alpha=0.5, color="#99FFFF", edgecolor='black',
                    label=rf"{mol} RMSD (Histogram)")
            ax.set_xlabel(r'RMSD (nm)')
            ax.set_ylabel('Probability density')
            ax.set_title(f'Probability Density Function of {mol} RMSD\nref: {os.path.basename(self.ref)}')
            ax.legend()
            fig.savefig(os.path.join(self.save_path, f"draw_RMSD_{mol}_ref_{os.path.basename(self.ref)}.png"))
            plt.close(fig)
            print(f"RMSD 绘图完成！结果已保存至文件 draw_RMSD_{mol}_ref_{os.path.basename(self.ref)}.png")
