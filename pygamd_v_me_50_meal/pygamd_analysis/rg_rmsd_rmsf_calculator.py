import os
import math

import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm
from scipy.ndimage import gaussian_filter

from pygamd_v_me_50_meal.data import Data


# 定义一个类，用于计算 Rg, RMSD, RMSF。
class RgRMSDRMSFCalculator:
    def __init__(self, path, data: Data, ref):
        """
        用于计算 Rg, RMSD, RMSF。
        """
        self.path = path
        self.data = data

        self.mol_class_dict = self.data.mol_class_dict
        self.length_dict = self.data.length_dict

        self.chain_path = os.path.join(self.path, "chain_xyz_unwrapping/")

        self.save_path = os.path.join(self.path, "draw_log/")
        os.makedirs(self.save_path, exist_ok=True)

        self.ref = ref
        self.cal_class = {}
        self.cal_class_rg, self.cal_class_rmsd, self.cal_class_rmsf = [], [], []
        self.rg_results, self.rmsd_results, self.rmsf_results = {}, {}, {}
        self.cur_chain_class = ""
        self.init_pos = []
        self.domain = []


    @staticmethod
    def kabsch_align(p, q) -> tuple[np.ndarray, np.ndarray]:
        """
        Kabsch算法，用于对齐两个结构。
        :param p: 参考结构
        :param q: 待对齐的结构
        :return: 将参考结构居中的 P_centered: np.ndarray, 对齐后的 Q_aligned: np.ndarray
        """
        # 使用Kabsch算法将结构Q对齐到结构P。
        p_centroid = np.mean(p, axis=0)
        q_centroid = np.mean(q, axis=0)

        p_centered = p - p_centroid
        q_centered = q - q_centroid

        covariance_matrix = q_centered.T @ p_centered

        U, _, Vt = np.linalg.svd(covariance_matrix)
        rotation_matrix = U @ Vt

        q_aligned = np.dot(q_centered, rotation_matrix)

        return p_centered, q_aligned


    def get_init_pos(self):
        init_pos = []
        with open(self.ref + ".xml", 'r') as f:
            f_lines = f.readlines()
            position_flag = 0
            cnt = 0
            for line in f_lines:
                if "<position" in line:
                    position_flag = 1
                    continue
                elif "</position" in line:
                    break

                if position_flag:
                    cnt += 1
                    init_pos.append(list(map(float, line.strip('\n').split())))
        return init_pos


    def cal_rg(self, chain_xyz):
        chain_xyz = np.array(chain_xyz)
        # mass_dict = {'H': 137.14, 'D': 115.09, 'R': 156.19, 'F': 147.18, 'A': 71.07,
        #              'C': 103.14, 'G': 57.05, 'Q': 128.13, 'E': 129.11, 'K': 128.17,
        #              'L': 113.16, 'M': 131.20, 'N': 114.10, 'S': 87.08, 'Y': 163.18,
        #              'T': 101.11, 'I': 113.16, 'W': 186.22, 'P': 97.12, 'V': 99.13,
        #              'HD': 137.14, 'DD': 115.09, 'RD': 156.19, 'FD': 147.18, 'AD': 71.07,
        #              'CD': 103.14, 'GD': 57.05, 'QD': 128.13, 'ED': 129.11, 'KD': 128.17,
        #              'LD': 113.16, 'MD': 131.20, 'ND': 114.10, 'SD': 87.08, 'YD': 163.18,
        #              'TD': 101.11, 'ID': 113.16, 'WD': 186.22, 'PD': 97.12, 'VD': 99.13,
        #              }
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


    def cal_rmsd(self, chain_xyz):
        if self.domain:
            init_chain_xyz = np.array(self.init_pos[self.domain[0]-1:self.domain[1]])
        else:
            init_chain_xyz = np.array(self.init_pos)
        chain_xyz = np.array(chain_xyz)

        aa_num = len(init_chain_xyz)
        if aa_num != len(chain_xyz):
            exit("Error! Length of chain is not equal to that of initial chain.!")

        init_chain_xyz, aligned_chain_xyz = self.kabsch_align(init_chain_xyz, chain_xyz)

        rmsd = 0
        for i in range(aa_num):
            rmsd += sum(sum([(aligned_chain_xyz[i][j] - init_chain_xyz[i][j]) ** 2]) for j in range(3))
        return math.sqrt(rmsd / aa_num)


    def cal_rmsf(self, chain_xyz):
        if self.domain:
            init_chain_xyz = np.array(self.init_pos[self.domain[0] - 1:self.domain[1]])
        else:
            init_chain_xyz = np.array(self.init_pos)
        chain_xyz = np.array(chain_xyz)

        aa_num = len(init_chain_xyz)
        if aa_num != len(chain_xyz):
            exit("Error! Length of chain is not equal to that of initial chain.!")

        init_chain_xyz, aligned_chain_xyz = self.kabsch_align(init_chain_xyz, chain_xyz)

        # 计算 aligned_chain_xyz 和 init_chain_xyz 的差异
        difference = aligned_chain_xyz - init_chain_xyz
        # 计算均方根波动（RMSF）
        rmsf = np.sum(difference ** 2, axis=1)  # 直接计算每个残基的平方和
        return rmsf


    def process_chain_file(self, name):
        with open(os.path.join(self.chain_path, name), 'r') as fp:
            x_mat: dict = eval(fp.read())

        for i in x_mat.keys():
            if i not in self.cal_class_rg and i not in self.cal_class_rmsd and i not in self.cal_class_rmsf:
                continue
            else:
                if self.domain:
                    cur_chain_xyz = map(lambda x: x[self.domain[0]-1:self.domain[1]], x_mat[i])
                else:
                    cur_chain_xyz = x_mat[i]
            if i in self.cal_class_rg and self.cal_class["Rg"]:
                self.cur_chain_class = i
                self.rg_results[i] = [] if i not in self.rg_results else self.rg_results[i]
                cur_rg_results = list(map(self.cal_rg, cur_chain_xyz))
                self.rg_results[i] += cur_rg_results
            if i in self.cal_class_rmsd and self.cal_class["RMSD"]:
                self.rmsd_results[i] = [] if i not in self.rmsd_results else self.rmsd_results[i]
                cur_rmsd_results = list(map(self.cal_rmsd, cur_chain_xyz))
                self.rmsd_results[i] += cur_rmsd_results
            if i in self.cal_class_rmsf and self.cal_class["RMSF"]:
                cur_rmsf_results = np.array([self.cal_rmsf(x) for x in cur_chain_xyz])
                cur_rmsf_results = np.sqrt(np.mean(cur_rmsf_results, axis=0))
                self.rmsf_results[i] += cur_rmsf_results
        # print(self.rmsf_results)


    def calculate(self, cal_class):
        self.cal_class = cal_class
        chain_files = os.listdir(self.chain_path)

        if self.cal_class["Rg"]:
            self.cur_chain_class = None
            self.cal_class_rg = list(map(int, input(f"请输入想要计算 Rg 的分子序号：\n{self.data.molecules}\n").split(',')))
            self.cal_class_rg = [self.data.mol_class_list[i] for i in self.cal_class_rg]
            print(f"即将计算 Rg 的分子：{self.cal_class_rg}")
            self.rg_results = {}

        if self.cal_class["RMSD"]:
            if not self.ref:
                self.ref = input("请指定含有单个参考结构的 XML 文件名：")
            self.ref = self.ref.replace(".xml", "")

            self.init_pos = self.get_init_pos()  # 读取参考结构的初始位置
            self.cal_class_rmsd = list(
                map(int, input(f"请输入想要计算 RMSD 的分子序号：\n{self.data.molecules}\n").split(',')))
            self.cal_class_rmsd = [self.data.mol_class_list[i] for i in self.cal_class_rmsd]
            print(f"即将计算 RMSD 的分子：{self.cal_class_rmsd}")
            self.rmsd_results = {}

        if self.cal_class["RMSF"]:
            if not self.ref:
                self.ref = input("请指定含有单个参考结构的 XML 文件名：")
            self.ref = self.ref.replace(".xml", "")

            self.init_pos = self.get_init_pos()  # 读取参考结构的初始位置
            self.cal_class_rmsf = list(
                map(int, input(f"请输入想要计算 RMSF 的分子序号：\n{self.data.molecules}\n").split(',')))
            self.cal_class_rmsf = [self.data.mol_class_list[i] for i in self.cal_class_rmsf]
            print(f"即将计算 RMSF 的分子：{self.cal_class_rmsf}")

        domain = input(
            "若只需要计算结构域，请输入该结构域的起始残基编号和末尾残基编号（从 1 开始），以-分隔，如 159-522，若有多个结构域，请以英文逗号分隔。\n"
            "否则请直接回车：")
        if domain:
            domain = list(map(int, domain.split('-')))
            self.domain = domain
            print(f"即将计算结构域：{domain}")
            self.rmsf_results = {i: np.zeros(domain[1] - domain[0] + 1) for i in self.cal_class_rmsf}
        else:
            self.domain = None

        # 使用 tqdm 包装可迭代对象以显示进度条
        list(tqdm(map(self.process_chain_file, sorted(chain_files)),
                  total=len(chain_files),
                  desc="计算中",
                  colour='cyan',
                  bar_format = '{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]',
                  ncols=100))

        # 保存结果
        for cal_type in cal_class.keys():
            if cal_class[cal_type]:
                if cal_type == "Rg":
                    result_file = os.path.join(self.save_path, f"draw_{cal_type}.log")
                elif cal_type == "RMSF":
                    self.rmsf_results = {i: (self.rmsf_results[i] / len(chain_files)).tolist() for i in self.cal_class_rmsf}
                    result_file = os.path.join(self.save_path, f"draw_{cal_type}_ref_{os.path.basename(self.ref)}.log")
                else:
                    result_file = os.path.join(self.save_path, f"draw_{cal_type}_ref_{os.path.basename(self.ref)}.log")
                with open(result_file, 'w') as f:
                    exec(f"f.write(str(self.{cal_type.lower()}_results))")
                print(f"{cal_type} 计算完成！结果已保存至文件 {result_file}")

                self.draw_rg_distribution() if cal_type == "Rg" else None
                self.draw_rmsd_distribution() if cal_type == "RMSD" else None
                self.draw_rmsf() if cal_type == "RMSF" else None
                print(f"{cal_type} 绘图完成！结果已保存至文件 {os.path.join(self.save_path, f'draw_{cal_type}_ref_{os.path.basename(self.ref)}.png')}")


    def draw_rg_distribution(self):
        with open(os.path.join(self.save_path, f"draw_Rg.log"), 'r') as f:
            rg_results = eval(f.read())  # 读取 RMSD 结果

        for mol in rg_results.keys():
            fig, ax = plt.subplots(figsize=(12, 9), dpi=300)
            init_rg_list = rg_results[mol]
            rg = np.array(sorted(init_rg_list))
            frame_num = int(len(init_rg_list) / self.data.mol_class_dict[mol][0])
            print(f"{mol} 的帧数：{frame_num}")

            # 划分 bin
            bins = 400
            hist, bin_edges = np.histogram(rg, bins=bins)
            hist = gaussian_filter(hist, sigma=1)
            # 计算概率
            probabilities = hist / sum(hist)

            plt.plot(bin_edges[:-1], probabilities, label=f"{mol} RMSD")
            plt.xlabel(r'Rg ($\AA$)')
            plt.ylabel('probability')
            plt.title(f'probability distribution')
            plt.legend()
            plt.savefig(os.path.join(self.save_path, f"draw_Rg_{mol}.png"))
            plt.close(fig)


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
            bins = 400
            hist, bin_edges = np.histogram(rmsd, bins=bins)
            hist = gaussian_filter(hist, sigma=5)
            # 计算概率
            probabilities = hist / sum(hist)

            plt.plot(bin_edges[:-1], probabilities, label=f"{mol} RMSD")
            plt.xlabel(r'RMSD ($\AA$)')
            plt.ylabel('probability')
            plt.title(f'probability distribution\nref: {os.path.basename(self.ref)}')
            plt.legend()
            plt.savefig(os.path.join(self.save_path, f"draw_RMSD_{mol}_ref_{os.path.basename(self.ref)}.png"))
            plt.close(fig)


    def draw_rmsf(self):
        if not self.ref:
            self.ref = input("请指定含有单个参考结构的 XML 文件名：")
        self.ref = self.ref.replace(".xml", "")

        with open(os.path.join(self.save_path, f"draw_RMSF_ref_{os.path.basename(self.ref)}.log"), 'r') as f:
            rmsf_results = eval(f.read())  # 读取 RMSF 结果
            for mol in rmsf_results.keys():
                fig, ax = plt.subplots(figsize=(12, 9), dpi=300)
                ax.plot(rmsf_results[mol], label=f"{mol} RMSF")
                ax.set_xlabel('residue number')
                ax.set_ylabel('RMSF')
                ax.set_title(f'RMSF of {mol}\nref: {os.path.basename(self.ref)}')
                ax.legend()
                plt.savefig(os.path.join(self.save_path, f"draw_RMSF_{mol}_ref_{os.path.basename(self.ref)}.png"))
                plt.close(fig)
