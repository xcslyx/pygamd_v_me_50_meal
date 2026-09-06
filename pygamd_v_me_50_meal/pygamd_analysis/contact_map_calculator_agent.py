import os
import re
import json
import subprocess

import numpy as np
import pandas as pd
import torch as torch
import multiprocessing as mp
import matplotlib.pyplot as plt

from tqdm import tqdm
from multiprocessing import Pool
# from scipy.ndimage import gaussian_filter
from matplotlib.ticker import ScalarFormatter, MaxNLocator

from pygamd_v_me_50_meal.data import Data
from pygamd_v_me_50_meal.Functions import Functions

# 加载消息文件
with open(os.path.join(os.path.dirname(__file__), 'message.json'), 'r', encoding='utf-8') as f:
    messages = json.load(f)
    msg = messages['contact_map_calculator_message']

# 定义一个类用于计算和绘制 contact map
class ContactMapCalculatorAgent:
    """
    用于计算和绘制 contact map 的类。
    """
    def __init__(self, 
                    path: str, data: Data, gpu_choice: str | int=0, r_cut: float=4.0,
                    cm_class_list: str="", balance_cut: str | None=None, domain: str='',
                    draw_limit: bool=False, lang: str='zh'
                    ):
        """
        初始化 ContactMapCalculatorAgent 类
        :param path: 体系路径
        :param data: 体系数据
        :param gpu_choice: GPU 选择
        :param r_cut: 截止距离
        :param cm_class_list: contact map 类型选择
        :param balance_cut: 选择一部分轨迹进行计算
        :param domain: 结构域选择
        :param draw_limit: 是否绘制 contact map 限制
        :return:
        """
        self.path, self.data, self.lang = path, data, lang

        self.mol_class_dict = self.data.mol_class_dict
        self.mol_class_list = self.data.mol_class_list
        self.length_dict = self.data.length_dict

        self.cm_class = []
        
        self.draw_limit = draw_limit

        self.chain_path = os.path.join(self.path, f"chain_xyz/")

        if torch.cuda.is_available():
            if gpu_choice.strip().upper() == "CPU":
                self.device = torch.device("cpu")
            elif gpu_choice.strip():
                self.device = torch.device(f"cuda:{gpu_choice}")
            else:
                self.device = torch.device("cuda:0")
        else:
            self.device = torch.device("cpu")
        print(f"Using device：{self.device}")
        if mp.get_start_method(allow_none=True) is None:
            mp.set_start_method('spawn')

        self.cm_class_list = cm_class_list.split(',')
        if "all" in self.cm_class_list or self.cm_class_list == [""]:
            self.cm_class_list = [[i, j] for i in range(len(self.mol_class_list)) for j in range(i, len(self.mol_class_list))]
        else:
            self.cm_class_list = list(map(lambda x: list(map(lambda y: int(y) - 1, x.split('-'))), self.cm_class_list))

        self.balance_cut = balance_cut

        self.domains, self.domain = None, None
        if domain:
            if ',' in domain:
                domains = domain.split(',')
                self.domains = list(map(lambda x: list(map(int, x.split('-'))), domains))
                print(f"Domains to be calculated: {domains}")
            else:
                self.domain = list(map(int, domain.split('-')))
                print(f"Domain to be calculated: {domain}")

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
        cm_class_0, cm_class_1 = self.cm_class

        if self.domain is not None:
            cm_matrix = torch.zeros((self.domain[1] - self.domain[0] + 1, self.domain[1] - self.domain[0] + 1),
                                    device=self.device)
        elif self.domains is not None:
            length = 0
            for domain in self.domains:
                length += domain[1] - domain[0] + 1
            cm_matrix = torch.zeros(length, length,
                                    device=self.device)
        else:
            cm_matrix = torch.zeros((self.mol_class_dict[cm_class_0][1], self.mol_class_dict[cm_class_1][1]),
                                    device=self.device)

        contact_matrix = torch.zeros((self.mol_class_dict[cm_class_0][0], self.mol_class_dict[cm_class_1][0]),
                                    device=self.device)

        for ii in range(len(x_mat[cm_class_0])):
            if cm_class_0 != cm_class_1 or (len(x_mat[cm_class_0]) == len(x_mat[cm_class_1]) == 1):
                range_j = range(len(x_mat[cm_class_1]))
            else:
                range_j = range(ii+1, len(x_mat[cm_class_1]))

            for jj in range_j:
                if self.domain is not None:
                    x_a = torch.tensor(x_mat[cm_class_0][ii][self.domain[0] - 1:self.domain[1]],
                                       device=self.device)
                    x_b = torch.tensor(x_mat[cm_class_1][jj][self.domain[0] - 1:self.domain[1]],
                                       device=self.device)
                elif self.domains is not None:
                    x_a_list = []
                    x_b_list = []
                    for domain in self.domains:
                        x_a_list.extend(x_mat[cm_class_0][ii][domain[0] - 1:domain[1]])
                        x_b_list.extend(x_mat[cm_class_1][jj][domain[0] - 1:domain[1]])
                    x_a = torch.tensor(x_a_list, device=self.device)
                    x_b = torch.tensor(x_b_list, device=self.device)
                else:
                    x_a = torch.tensor(x_mat[cm_class_0][ii], device=self.device)
                    x_b = torch.tensor(x_mat[cm_class_1][jj], device=self.device)
                # 计算欧氏距离
                d = Functions.euclidean_distances(x_a, x_b)
                c = d < self.avg_sigma_mat  # 创建布尔数组

                # 若有接触，contact_matrix[ii][jj] = 1，否则为 0
                contact_matrix[ii][jj] = 1 if c.any().item() else 0

                if cm_class_0 != cm_class_1:
                    cm_matrix += c
                else:
                    cm_matrix += c
                    if not len(x_mat[cm_class_0]) == len(x_mat[cm_class_1]) == 1:
                        cm_matrix += c.transpose(0, 1)

                    # cm_matrix += c.transpose(0, 1) + c

        cm_matrix = cm_matrix.cpu().numpy()
        contact_matrix = contact_matrix.cpu().numpy()
        # 保存 contact map
        with open(os.path.join(self.cur_cm_path, name), 'w') as m:
            for i in range(cm_matrix.shape[0]):
                for j in range(cm_matrix.shape[1]):
                    m.write(str(cm_matrix[i][j]))
                    m.write(' ')
                m.write('\n')
        
        with open(os.path.join(self.cur_cm_path, name.replace(".xml", "_contact.log")),  'w') as m:
            for i in range(contact_matrix.shape[0]):
                for j in range(contact_matrix.shape[1]):
                    m.write(str(contact_matrix[i][j]))
                    m.write(' ')
                m.write('\n')
        return True


    def calculate_contact_map_parallel(self):
        files = sorted(os.listdir(self.chain_path))
        if self.balance_cut is not None:
            start, end = list(map(int, self.balance_cut.split('-')))
            files = files[start: end+1]

        for cm_class in self.cm_class_list:
            self.cm_class = [self.data.mol_class_list[cm_class[0]], self.data.mol_class_list[cm_class[1]]]
            self.avg_sigma_mat = torch.tensor(self.r_cut + Functions.cal_sigma_mat(self.sequence[self.cm_class[0]], self.sequence[self.cm_class[1]]),
                                               device=self.device)
            print(f"Calculating contact map of {' 和 '.join(self.cm_class)}")
            self.cur_cm_path = os.path.join(self.cm_path, f"{self.cm_class[0]}_{self.cm_class[1]}_r_cut_{self.r_cut:.2f}")
            if os.path.exists(self.cur_cm_path):
                print(f"✅ {self.cur_cm_path} has already existed, removing...")
                subprocess.run(f"rm -rf {self.cur_cm_path}", shell=True)
            os.makedirs(self.cur_cm_path, exist_ok=True)
            with Pool(processes=4) as pool:
                # 使用 tqdm 包装可迭代对象
                list(tqdm(pool.imap(self.calculate_contact_map, files),
                          total=len(files),
                          desc="Calculating contact map",
                          colour='cyan',
                          bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]',
                          ncols=100))

        self.average_contact_map()
        print(f"所有 contact map 已经计算完成，保存在目录 {os.path.join(self.path, 'draw_log')} 中。")
        print(f"Finished calculating all contact maps and saved in the directory {os.path.join (self.path,'draw_log')}.")
        self.draw_contact_map()


    def avg_cm_file(self, cm_file):
        with open(os.path.join(self.cur_cm_path, cm_file), 'r') as f:
            data_matrix = []
            for line in f.readlines():
                float_line = list(map(float, line.strip().split(" ")))
                data_matrix.append(float_line)
            cur_contact_map = np.array(data_matrix)

        # 计算接触数
        contact_number = np.sum(cur_contact_map)
        with open(os.path.join(self.cur_cm_path, cm_file.replace('.xml', '_cn.log')), 'w') as save_file:
            save_file.write(f"{contact_number}\n")

        with open(os.path.join(self.cur_cm_path, cm_file.replace('.xml', '_contact.log')), 'r') as f:
            contact_matrix = []
            for line in f.readlines():
                float_line = list(map(float, line.strip().split(" ")))
                contact_matrix.append(float_line)
            cur_contact_matrix = np.array(contact_matrix)
        
        return cur_contact_map, cur_contact_matrix


    def average_contact_map(self, get_cotact_mol: bool = False):
        for cm_class in self.cm_class_list:
            self.cm_class = [self.data.mol_class_list[cm_class[0]], self.data.mol_class_list[cm_class[1]]]
            self.cur_cm_path = os.path.join(self.cm_path, f"{self.cm_class[0]}_{self.cm_class[1]}_r_cut_{self.r_cut:.2f}")
            if not os.path.exists(self.cur_cm_path):
                print(f"未找到 {self.cur_cm_path} 文件夹，请先进行计算。")
                return
            cm_files = sorted([f for f in os.listdir(self.cur_cm_path) if f.endswith('.xml')])

            # 使用多进程读取和处理 CM 文件
            with Pool(processes=4) as pool:
                results = list(tqdm(pool.imap(self.avg_cm_file, cm_files),
                                    total=len(cm_files),
                                    desc="Averaging contact map",
                                    colour='cyan',
                                    bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]',
                                    ncols=100))

            file_prefix = f"draw_cm_{self.cm_class[0]}_{self.cm_class[1]}_r_cut_{self.r_cut:.2f}"

            cm_mat = np.zeros_like(results[0][0])
            contact_mat = np.zeros_like(results[0][1])
            for dataMat in results:
                cm_mat += dataMat[0]
                contact_mat += dataMat[1]

            avg_cm_mat = cm_mat / len(cm_files)
            avg_contact_mat = contact_mat / len(cm_files)
            # 保存平均后的 contact map
            with open(os.path.join(self.cm_path, f"{file_prefix}_avg_matrix.log"), 'w') as save_file:
                    for i in avg_cm_mat:
                        save_file.write(" ".join(map(str, i)) + '\n')

            with open(os.path.join(self.cm_path, f"{file_prefix}_avg.log"), 'w') as save_file:
                for i in range(len(avg_cm_mat)):
                    for j in range(len(avg_cm_mat[i])):
                        save_file.write(f"{i + 1} {j + 1} {avg_cm_mat[i][j]}\n")
            
            with open(os.path.join(self.cm_path, f"{file_prefix}_avg_contact_matrix.log"), 'w') as save_file:
                    for i in avg_contact_mat:
                        save_file.write(" ".join(map(str, i)) + '\n')

            with open(os.path.join(self.cm_path, f"{file_prefix}_avg_contact.log"), 'w') as save_file:
                for i in range(len(avg_contact_mat)):
                    for j in range(len(avg_contact_mat[i])):
                        save_file.write(f"{i + 1} {j + 1} {avg_contact_mat[i][j]}\n")

            # 保存 cn_list
            with open(os.path.join(self.cm_path, f"{file_prefix}_cn_list.log"), 'w') as cn_save_file:
                for cm_idx, cm_file in enumerate(cm_files):
                    cn_file = cm_file.replace('.xml', '_cn.log')
                    with open(os.path.join(self.cur_cm_path, cn_file), 'r') as f:
                        cn = float(f.readline().strip())
                        cn_save_file.write(f"{cm_idx+1} {cn}\n")


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

            print(f"✅ Drawing contact map of {cm_class[0]}-{cm_class[1]}, r_cut is {r_cut}")

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
            # data_mat = gaussian_filter(data_mat, sigma=1.5)

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

            # flights = data_mat.pivot("residues", "residues", "contact number")
            fig, ax = plt.subplots(figsize=(12, 9), dpi=300)

            if self.draw_limit:
                if max_value_num <= 2:
                    vmax = 2 * 10 ** exponent
                elif max_value_num <= 3:
                    vmax = 3 * 10 ** exponent
                elif max_value_num <= 5:
                    vmax = 5 * 10 ** exponent
                elif max_value_num <= 9:
                    vmax = 9 * 10 ** exponent
                else:
                    vmax = 10 * 10 ** exponent
                im = ax.imshow(data_mat, cmap=plt.get_cmap('jet'), aspect='auto', vmin=0., vmax=vmax)
            else:
                im = ax.imshow(data_mat, cmap=plt.get_cmap('jet'), aspect='auto',)
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
            if abs(exponent) >= 2:
                formatter = ScalarFormatter()
                formatter.set_scientific(True)
                formatter.set_powerlimits((-1, 1))  # 设置科学计数法的显示范围
                formatter.set_useOffset(False)  # 确保不使用偏移量
                cbar.ax.yaxis.set_major_formatter(formatter)
            # plt.tight_layout()
            plt.savefig(os.path.join(self.cm_path, f"draw_cm_{cm_class[0]}_{cm_class[1]}_r_cut_{r_cut}_avg.png"), dpi=300)
            plt.close(fig)
