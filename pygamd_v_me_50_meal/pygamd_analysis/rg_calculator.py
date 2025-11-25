import os
import math

import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm
from scipy.ndimage import gaussian_filter

from pygamd_v_me_50_meal.data import Data


# 定义一个类，用于计算 Rg, RMSD, RMSF。
class RgCalculator:
    def __init__(self, path, data: Data):
        """
        用于计算 Rg, RMSD, RMSF。
        """
        self.path = path
        self.data = data

        self.chain_path = os.path.join(self.path, "chain_xyz_unwrapping/")

        self.save_path = os.path.join(self.path, "draw_log/")
        os.makedirs(self.save_path, exist_ok=True)

        self.balance_cut = ""

        self.cal_class = {}
        self.cal_class_rg = []
        self.rg_results = {}
        self.cur_chain_class = ""
        self.init_pos = []
        self.domain = []

    def cal_rg(self, chain_xyz):
        chain_xyz = np.array(chain_xyz)

        # calculate the core of mass of the chain
        with open(f"{os.path.join(self.path, self.data.system_name)}_sequence.txt") as f:
            if self.domain:
                domain_sequence = eval(f.read())[self.cur_chain_class][self.domain[0] - 1:self.domain[1]]
            else:
                domain_sequence = eval(f.read())[self.cur_chain_class][:]

        mass_list = list(map(float, [i[1] for i in domain_sequence]))
        # print(mass_list)
        if len(domain_sequence) != len(chain_xyz):
            exit(f"Error! Length of domain_sequence ({len(domain_sequence)}) is not equal to that of chain_xyz ({len(chain_xyz)}).!")
        aa_num = len(domain_sequence)
        mass = sum(mass_list)
        mass_core = [sum([chain_xyz[i][j] * mass_list[i] for i in range(aa_num)]) / mass for j in range(3)]
        # print(mass_core)
        r_g = 0
        for i in range(aa_num):
            r_g += sum(sum([(chain_xyz[i][j] - mass_core[j]) ** 2 * mass_list[i]]) for j in range(3))
        return math.sqrt(r_g / mass)

    def process_chain_file(self, name):
        with open(os.path.join(self.chain_path, name), 'r') as fp:
            x_mat: dict = eval(fp.read())

        for i in x_mat.keys():
            if i not in self.cal_class_rg:
                continue
            else:
                if self.domain:
                    cur_chain_xyz = map(lambda x: x[self.domain[0]-1:self.domain[1]], x_mat[i])
                else:
                    cur_chain_xyz = x_mat[i]
            if i in self.cal_class_rg:
                self.cur_chain_class = i
                self.rg_results[i] = [] if i not in self.rg_results else self.rg_results[i]
                cur_rg_results = list(map(self.cal_rg, cur_chain_xyz))
                self.rg_results[i] += cur_rg_results


    def calculate(self):
        if not self.balance_cut:
            self.balance_cut = input(
                "请输入需要截取的平衡后的文件索引，索引从 1 开始，格式为‘开始,结束’，例如：1000,2000，直接回车则不截取：")
        if not self.balance_cut:
            chain_files = os.listdir(self.chain_path)
        else:
            start, end = list(map(int, self.balance_cut.split(',')))
            chain_files = os.listdir(self.chain_path)[start - 1: end]

        self.cur_chain_class = None
        print(f"您当前的分子类型有：\n{self.data.molecules}\n")
        self.cal_class_rg = list(map(lambda x: int(x) - 1, input(f"请输入想要计算 Rg 的分子序号, 以 ',' 分隔：").split(',')))
        self.cal_class_rg = [self.data.mol_class_list[i] for i in self.cal_class_rg]
        print(f"即将计算 Rg 的分子：{self.cal_class_rg}")
        self.rg_results = {}

        domain = input(
            "若只需要计算结构域，请输入该结构域的起始残基编号和末尾残基编号（从 1 开始），以-分隔，如 159-522，若有多个结构域，请以英文逗号分隔。\n"
            "否则请直接回车：")
        if domain:
            domain = list(map(int, domain.split('-')))
            self.domain = domain
            print(f"即将计算结构域：{domain}")
        else:
            self.domain = None

        # 使用 tqdm 包装可迭代对象以显示进度条
        list(tqdm(map(self.process_chain_file, sorted(chain_files)),
                  total=len(chain_files),
                  desc="Calculating",
                  colour='cyan',
                  bar_format = '{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]',
                  ncols=100))

        # 保存结果
        result_file = os.path.join(self.save_path, f"draw_Rg.log")
        with open(result_file, 'w') as f:
            f.write(str(self.rg_results))
        print(f"Rg 计算完成！结果已保存至文件 draw_Rg.log")

        self.draw_rg_distribution()


    def draw_rg_distribution(self):
        with open(os.path.join(self.save_path, f"draw_Rg.log"), 'r') as f:
            rg_results = eval(f.read())  # 读取 RMSD 结果

        for mol in rg_results.keys():
            fig, ax = plt.subplots(figsize=(12, 9), dpi=300)
            init_rg_list = rg_results[mol]
            rg = np.array(sorted(init_rg_list))
            frame_num = int(len(init_rg_list) / self.data.mol_class_dict[mol][0])
            print(f"{mol} 的帧数：{frame_num}")

            # 计算概率密度函数
            bins = 40  # 划分 bin 的数量
            # 计算直方图
            hist, bin_edges = np.histogram(rg, bins=bins)
            # hist = gaussian_filter(hist, sigma=1)
            # 计算概率, 除以总链数和 bin 的宽度
            probabilities = hist / sum(hist) / (bin_edges[1] - bin_edges[0])

            bin_edges += (bin_edges[1] - bin_edges[0]) / 2  # 使 bin 居中
            ax.plot(bin_edges[:-1], probabilities, label=rf"{mol} $R_{{\mathrm{{g}}}}$")
            ax.hist(init_rg_list, bins=bins, density=True, alpha=0.5, color="#99FFFF", edgecolor='black',
                    label=rf"{mol} $R_{{\mathrm{{g}}}}$ (Histogram)")
            ax.set_xlabel(r'$R_{\mathrm{g}}$ (nm)')
            ax.set_ylabel('Probability density')
            ax.set_title(rf'Probability Density Function of {mol} $R_{{\mathrm{{g}}}}$')
            ax.legend()
            fig.savefig(os.path.join(self.save_path, f"draw_Rg_{mol}.png"))
            plt.close(fig)

            print(f"Rg 绘图完成！结果已保存至文件 {os.path.join(self.save_path, f'draw_Rg_{mol}.png')}")
