import os
import re

import numpy as np
import pandas as pd
import torch as torch
import multiprocessing as mp
import matplotlib.pyplot as plt

from tqdm import tqdm
from multiprocessing import Pool
from scipy.ndimage import gaussian_filter
from matplotlib.ticker import ScalarFormatter, MaxNLocator, FormatStrFormatter

from Functions import Functions

# 定义一个类用于计算和绘制 contact map
class ContactMapCalculator:
    """
    用于计算和绘制 contact map 的类。
    """
    def __init__(self, path, data, cm_choice, r_cut: float):
        """
        初始化 ContactMapCalculator 类
        :param path:
        """
        self.path = path
        self.data = data

        self.mol_class_dict = self.data.mol_class_dict
        self.mol_class_list = self.data.mol_class_list
        self.length_dict = self.data.length_dict

        self.chain_path = os.path.join(self.path, f"chain_xyz/")

        if mp.get_start_method(allow_none=True) is None:
            mp.set_start_method('spawn')

        self.device = None

        self.cm_class = []
        self.cm_class_list = []
        self.balance_cut = ""
        self.domain_cut = []
        if cm_choice != '/':
            cm_choice = cm_choice.split('/')
            self.cm_class_list = cm_choice[0].split(',')
            self.balance_cut = cm_choice[1]

        self.domain = None
        self.domains = None

        self.avg_sigma_mat = None

        self.sequence = {}
        with open(f"{os.path.join(self.path, self.data.system_name)}_sequence.txt") as f:
            sequence = eval(f.read())
            for cal_mol in self.mol_class_list:
                self.sequence[cal_mol] = [i[0] for i in sequence[cal_mol]]

        self.r_cut = r_cut

        self.draw_path = os.path.join(self.path, f"draw_log/")
        if not os.path.exists(self.draw_path):
            os.makedirs(self.draw_path, exist_ok=True)
        self.cm_path = os.path.join(self.path, f"draw_log/cm/")
        self.cur_cm_path = ""


    def calculate_contact_map(self, name):
        x_mat: list = eval(open(self.chain_path + name, 'r').read())
        if self.domain is not None:
            cm_matrix = torch.zeros((self.domain[1] - self.domain[0] + 1, self.domain[1] - self.domain[0] + 1),
                                    device=self.device)
        elif self.domains is not None:
            length = 0
            for domain in self.domains:
                length += domain[1] - domain[0] + 1
            cm_matrix = torch.zeros(length, length, device=self.device)
        else:
            cm_matrix = torch.zeros(
                (self.mol_class_dict[self.cm_class[0]][1], self.mol_class_dict[self.cm_class[1]][1]),
                                    device=self.device)

        for ii in range(len(x_mat[self.cm_class[0]])):
            range_j = range(len(x_mat[self.cm_class[1]])) if self.cm_class[0] != self.cm_class[1] else range(ii+1, len(x_mat[self.cm_class[1]]))
            for jj in range_j:
                # if ii < jj or (len(x_mat[self.cm_class[0]]) == 1 and len(x_mat[self.cm_class[1]]) == 1):
                if self.domain is not None:
                    x_a = torch.tensor(x_mat[self.cm_class[0]][ii][self.domain[0] - 1:self.domain[1]],
                                       device=self.device)
                    x_b = torch.tensor(x_mat[self.cm_class[1]][jj][self.domain[0] - 1:self.domain[1]],
                                       device=self.device)
                elif self.domains is not None:
                    x_a_list = []
                    x_b_list = []
                    for domain in self.domains:
                        x_a_list.extend(x_mat[self.cm_class[0]][ii][domain[0] - 1:domain[1]])
                        x_b_list.extend(x_mat[self.cm_class[1]][jj][domain[0] - 1:domain[1]])
                    x_a = torch.tensor(x_a_list, device=self.device)
                    x_b = torch.tensor(x_b_list, device=self.device)
                else:
                    x_a = torch.tensor(x_mat[self.cm_class[0]][ii], device=self.device)
                    x_b = torch.tensor(x_mat[self.cm_class[1]][jj], device=self.device)
                # 计算欧氏距离
                d = Functions.euclidean_distances(x_a, x_b)
                c = d < self.avg_sigma_mat  # 创建布尔数组

                if self.cm_class[0] != self.cm_class[1]:
                    cm_matrix += c
                else:
                    cm_matrix += c.transpose(0, 1)
                    cm_matrix += c

                    # cm_matrix += c.transpose(0, 1) + c

        cm_matrix = cm_matrix.cpu().numpy()
        # 保存 contact map
        with open(os.path.join(self.cur_cm_path, name), 'w') as m:
            for i in range(len(cm_matrix)):
                for j in range(len(cm_matrix[0])):
                    m.write(str(cm_matrix[i][j]))
                    m.write(' ')
                m.write('\n')
        return True


    def calculate_contact_map_parallel(self):
        if torch.cuda.is_available():
            gpu_choice = input("即将使用 GPU 加速计算 Contact Map，请指定 GPU 编号，或直接回车使用默认 0 号 GPU。若想使用 CPU，请输入 CPU：")
            if gpu_choice.strip().upper() == "CPU":
                self.device = torch.device("cpu")
            elif gpu_choice.strip():
                self.device = torch.device(f"cuda:{gpu_choice}")
            else:
                self.device = torch.device("cuda:0")
        else:
            self.device = torch.device("cpu")
        print(f"使用设备：{self.device}")
        if mp.get_start_method(allow_none=True) is None:
            mp.set_start_method('spawn')

        if not self.cm_class_list:
            print(f"\n您的分子类型有：\n{self.data.molecules}")
            self.cm_class_list = input("请输入您想要计算的 contact map 的分子对，以逗号分隔，如“1-1,2-2,1-2”。\n"
                                        "如需计算全部分子组合，请输入 all 或直接回车：").split(',')
        if "all" in self.cm_class_list or self.cm_class_list == [""]:
            self.cm_class_list = [[i, j] for i in range(len(self.data.mol_class_list)) for j in range(i, len(self.data.mol_class_list))]
        else:
            self.cm_class_list = list(map(lambda x: list(map(int, x.split('-'))), self.cm_class_list))

        domain = input(
            "若只需要计算结构域，请输入该结构域的起始残基编号和末尾残基编号（从 1 开始），以-分隔，如 159-522。\n"
            "注：仅限计算的分子类型相同的结构域，如0-0、1-1。如需计算不同分子类型的结构域，请另外计算。\n"
            "若有多个结构域，请以英文逗号分隔，否则请直接回车：")
        if domain:
            if ',' in domain:
                domains = domain.split(',')
                domains = list(map(lambda x: list(map(int, x.split('-'))), domains))
                self.domains = domains
                print(f"即将计算结构域：{domains}")
            else:
                domain = list(map(int, domain.split('-')))
                self.domain = domain
                print(f"即将计算结构域：{domain}")

        if not self.balance_cut:
            self.balance_cut = input("请输入需要截取的平衡后的文件索引，格式为‘开始,结束’，例如：1000,2000，直接回车则不截取：")

        if not self.balance_cut:
            files = os.listdir(self.chain_path)
        else:
            start, end = list(map(int, self.balance_cut.split(',')))
            files = os.listdir(self.chain_path)[start: end+1]


        for cm_class in self.cm_class_list:
            self.cm_class = [self.data.mol_class_list[cm_class[0]], self.data.mol_class_list[cm_class[1]]]
            self.avg_sigma_mat = torch.tensor(self.r_cut + Functions.cal_sigma_mat(self.sequence[self.cm_class[0]], self.sequence[self.cm_class[1]]),
                                               device=self.device)
            print(f"正在计算 {' 和 '.join(self.cm_class)} 的 contact map")
            self.cur_cm_path = os.path.join(self.cm_path, f"{self.cm_class[0]}_{self.cm_class[1]}_r_cut_{self.r_cut:.2f}")
            if os.path.exists(self.cur_cm_path):
                print(f"✅ 已存在 {self.cur_cm_path} 文件夹，正在删除...")
                os.system(f"rm -rf {self.cur_cm_path}")
            os.makedirs(self.cur_cm_path, exist_ok=True)
            with Pool(processes=4) as pool:
                # 使用 tqdm 包装可迭代对象
                list(tqdm(pool.imap(self.calculate_contact_map, files),
                          total=len(files),
                          desc="计算中",
                          colour='cyan',
                          bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]',
                          ncols=100))

        self.average_contact_map()
        print(f"所有 contact map 已经计算完成，保存在目录 {os.path.join(self.path, 'draw_log')} 中。")
        self.draw_contact_map()


    def avg_cm_file(self, cm_file):
        with open(os.path.join(self.cur_cm_path, cm_file), 'r') as f:
            data_matrix = []
            for line in f.readlines():
                float_line = list(map(float, line.strip().split(" ")))
                data_matrix.append(float_line)
        return np.array(data_matrix)


    def average_contact_map(self):
        for cm_class in self.cm_class_list:
            self.cm_class = [self.data.mol_class_list[cm_class[0]], self.data.mol_class_list[cm_class[1]]]
            self.cur_cm_path = os.path.join(self.cm_path, f"{self.cm_class[0]}_{self.cm_class[1]}_r_cut_{self.r_cut:.2f}")
            if not os.path.exists(self.cur_cm_path):
                print(f"未找到 {self.cur_cm_path} 文件夹，请先进行计算。")
                return
            cm_files = os.listdir(self.cur_cm_path)

            # 使用多进程读取和处理 CM 文件
            with Pool(processes=4) as pool:
                results = list(tqdm(pool.imap(self.avg_cm_file, cm_files),
                                    total=len(cm_files),
                                    desc="平均中",
                                    colour='cyan',
                                    bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]',
                                    ncols=100))

            # 计算 cm_mat 的维度
            cm_mat = np.zeros_like(results[0])  # 假设所有的结果有相同的形状
            for dataMat in results:
                cm_mat += dataMat

            cn_list = np.zeros(len(results))
            for i, dataMat in enumerate(results):
                cn_list[i] = np.sum(dataMat)

            # 保存 cn_list
            with open(os.path.join(self.cm_path, f"draw_cm_{self.cm_class[0]}_{self.cm_class[1]}_r_cut_{self.r_cut:.2f}_cn_list.log"),
                      'w') as save_file:
                for i in range(len(cn_list)):
                    save_file.write(f"{i+1} {cn_list[i]}\n")

            avg_cm_mat = (cm_mat / len(cm_files))

            # 保存平均后的 contact map
            with open(os.path.join(self.cm_path, f"draw_cm_{self.cm_class[0]}_{self.cm_class[1]}_r_cut_{self.r_cut:.2f}_avg_matrix.log"),
                          'w') as save_file:
                    for i in avg_cm_mat:
                        save_file.write(" ".join(map(str, i)) + '\n')

            with open(os.path.join(self.cm_path, f"draw_cm_{self.cm_class[0]}_{self.cm_class[1]}_r_cut_{self.r_cut:.2f}_avg.log"),
                      'w') as save_file:
                for i in range(len(avg_cm_mat)):
                    for j in range(len(avg_cm_mat[i])):
                        save_file.write(f"{i + 1} {j + 1} {avg_cm_mat[i][j]}\n")


    def draw_contact_map(self):
        """
        绘制 contact map
        :return:
        """
        cm_files = sorted(os.listdir(self.cm_path))
        for cm_file in cm_files:
            match = re.match(r"draw_cm_(\w+)_(\w+)_r_cut_(\d+\.\d+)_avg_matrix.log", cm_file)
            if not match:
                continue

            cm_class = [match.group(1), match.group(2)]
            r_cut = float(match.group(3))
            if len(match.groups()) != 3:
                print(f"❌ 无法解析文件名 {cm_file} {match.groups()}")
                continue

            print(f"✅ Drawing contact map of {cm_class[0]}-{cm_class[1]}, r_cut 为 {r_cut}")

            with open(os.path.join(self.cm_path, cm_file), 'r') as f:
                data_matrix = []
                if self.domain is not None:
                    # f_lines = f.readlines()[self.domain[0]-1:self.domain[1]]
                    f_lines = [' '.join(['0'] * (self.domain[1] - self.domain[0] + 1))] * (self.domain[0] - 1)
                    # f_lines = np.zeros((self.domain[0] - 1, self.domain[1] - self.domain[0] + 1))
                    # print(len(f_lines))
                    f_lines.extend(f.readlines())
                    f_lines.extend([' '.join(['0'] * (self.domain[1] - self.domain[0] + 1))] * (self.data.length_dict[cm_class[1]] - self.domain[1]))
                    # print(len(f_lines))
                    for line in f_lines:
                        cur_line = ['0'] * (self.domain[0] - 1)
                        cur_line.extend(line.strip().split())
                        cur_line.extend(['0'] * (self.data.length_dict[cm_class[1]] - self.domain[1]))
                        # print(len(cur_line))
                        data_matrix.append(cur_line)
                else:
                    f_lines = f.readlines()[:self.data.length_dict[cm_class[0]]]
                    for line in f_lines:
                        data_matrix.append(line.strip().split()[:self.data.length_dict[cm_class[1]]])


            data_mat = pd.DataFrame(data_matrix, dtype=np.float64)
            # data_mat.to_csv(os.path.join(self.path, f"draw_log/draw_cm_{cm_class[0]}_{cm_class[1]}_r_cut_{r_cut}_avg_by_{avg_class}.csv"), index=False, header=False)
            # print(data_mat.shape)
            # data_mat = data_mat / data_mat.sum().sum()  # 归一化
            data_mat = gaussian_filter(data_mat, sigma=1.5)

            # 计算数据矩阵的最大值，并根据最大值设置颜色范围
            max_value = data_mat.max().max()
            # data_mat = data_mat / max_value
            if max_value < 1e-50:
                max_value = 0.0
            if np.isnan(max_value):
                max_value = 0.0
            max_value_sci = "{:.2e}".format(max_value)
            max_value_num = float(max_value_sci.split('e')[0])
            exponent = int(max_value_sci.split('e')[1])

            # vmax = 10 * 10 ** exponent  # 最大值设为10的整数倍

            # if max_value_num < 3:
            #     vmax = 3 * 10 ** exponent
            # elif 3 <= max_value_num < 5:
            #     vmax = 5 * 10 ** exponent
            # elif 5 <= max_value_num < 8:
            #     vmax = 8 * 10 ** exponent
            # else:
            #     vmax = 10 * 10 ** exponent

            if max_value_num < 5:
                vmax = 5 * 10 ** exponent
            else:
                vmax = 10 * 10 ** exponent

            # flights = data_mat.pivot("residues", "residues", "contact number")
            fig, ax = plt.subplots(figsize=(12, 9), dpi=300)
            # print(data_mat)
            # im = ax.imshow(data_mat, cmap=plt.get_cmap('Reds'), aspect='auto',)
            im = ax.imshow(data_mat, cmap=plt.get_cmap('jet'), aspect='auto', vmin=0., vmax=vmax)
            ax.invert_yaxis()

            ax.set_title(f"{cm_class[0]}-{cm_class[1]} contact map")
            ax.set_ylabel(f"{cm_class[0]} residues")
            ax.set_xlabel(f"{cm_class[1]} residues")

            # plt.xticks(rotation=-45)  # 设置x轴表明文字的方向
            if self.domain is not None:
                # 若有结构域，则应该记录为最大值而非实际长度
                length_y = self.domain[1]
                length_x = self.domain[1]

                # 设置刻度位置和标签
                xticks = np.arange(self.domain[0], self.domain[1] + 1, 50) - 1  # 刻度位置，包括最大值
                yticks = np.arange(self.domain[0], self.domain[1] + 1, 50) - 1  # 刻度位置，包括最大值
            else:
                length_y = self.data.length_dict[cm_class[0]]
                length_x = self.data.length_dict[cm_class[1]]

                # 设置刻度位置和标签
                xticks = np.arange(1, length_x + 1, 50) - 1  # 刻度位置，包括最大值
                yticks = np.arange(1, length_y + 1, 50) - 1  # 刻度位置，包括最大值

            if length_x - xticks[-1] > 25:
                # 在最大值处添加一个刻度
                xticks_pos = np.append(xticks, length_x - 1)
                xticks_label = np.append(xticks, length_x)
            else:
                xticks[-1] = length_x - 1
                xticks_pos = xticks
                xticks[-1] = length_x
                xticks_label = xticks

            if length_y - yticks[-1] > 25:
                # 在最大值处添加一个刻度
                yticks_pos = np.append(yticks, length_y - 1)
                yticks_label = np.append(yticks, length_y)
            else:
                yticks[-1] = length_y - 1
                yticks_pos = yticks
                yticks[-1] = length_y
                yticks_label = yticks

            # 设置刻度位置和标签
            ax.set_xticks(xticks_pos)
            # print("xticks_pos:", xticks_pos)
            ax.set_yticks(yticks_pos)
            # print("yticks_pos:", yticks_pos)
            ax.set_xticklabels(xticks_label)
            ax.set_yticklabels(yticks_label)

            ax.set_xlim(xticks_pos[0] - 0.5, xticks_pos[-1] + 0.5)
            ax.set_ylim(yticks_pos[0] - 0.5, yticks_pos[-1] + 0.5)

            # 画颜色条
            cbar = ax.figure.colorbar(im, ax=ax)
            # colorbar的设置
            # colorbar标签为‘contact number’，纵向放置
            cbar.ax.set_ylabel('contact number', rotation=-90, va="bottom")
            # 设置颜色条上的刻度为每 0.5 一个刻度
            locator = MaxNLocator(steps=[1, 2, 5])  # 默认的steps可能会导致不以0.5为间隔
            cbar.ax.yaxis.set_major_locator(locator)
            # 设置颜色条上的刻度为科学计数法
            if exponent >= 2:
                formatter = ScalarFormatter()
                formatter.set_scientific(True)
                formatter.set_powerlimits((-1, 1))  # 设置科学计数法的显示范围
                formatter.set_useOffset(False)  # 确保不使用偏移量
                cbar.ax.yaxis.set_major_formatter(formatter)
            # plt.tight_layout()
            plt.savefig(os.path.join(self.cm_path, f"draw_cm_{cm_class[0]}_{cm_class[1]}_r_cut_{r_cut}_avg.png"), dpi=300)
            plt.close(fig)
