import os
import math

import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm
from scipy.ndimage import gaussian_filter

from pygamd_v_me_50_meal.data import Data
from pygamd_v_me_50_meal.Functions import Functions


# 定义一个类，用于计算 Rg, RMSD, RMSF。
class RMSFCalculator:
    def __init__(self, path, data: Data, ref, balance_cut=None):
        """
        用于计算 Rg, RMSD, RMSF。
        """
        self.path = path
        self.data = data

        self.chain_path = os.path.join(self.path, "chain_xyz_unwrapping/")

        self.save_path = os.path.join(self.path, "draw_log/")
        os.makedirs(self.save_path, exist_ok=True)

        self.ref = input("请指定含有单个参考结构的 XML 文件名：").replace(".xml", "") if not ref else ref
        self.init_pos = self.get_init_pos(self.ref)  # 读取参考结构的初始位置

        self.cal_class_rmsf = []
        self.rmsf_results = {}
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


    @staticmethod
    def cal_rmsf(chain_xyz, init_pos, domain=None):
        if domain:
            init_chain_xyz = np.array(init_pos[domain[0] - 1: domain[1]])
        else:
            init_chain_xyz = np.array(init_pos)
        chain_xyz = np.array(chain_xyz)

        aa_num = len(init_chain_xyz)
        assert aa_num == len(chain_xyz), "The number of residues in the reference structure and the chain is not the same."

        init_chain_xyz, aligned_chain_xyz = Functions.kabsch_align(init_chain_xyz, chain_xyz)

        # 计算 aligned_chain_xyz 和 init_chain_xyz 的差异
        difference = aligned_chain_xyz - init_chain_xyz
        # 计算均方根波动（RMSF）
        rmsf = np.sum(difference ** 2, axis=1)  # 直接计算每个残基的平方和
        return rmsf


    def process_chain_file(self, name):
        with open(os.path.join(self.chain_path, name), 'r') as fp:
            x_mat: dict = eval(fp.read())

        for i in x_mat.keys():
            if i not in self.cal_class_rmsf:
                continue
            else:
                if self.domain:
                    cur_chain_xyz = map(lambda x: x[self.domain[0]-1:self.domain[1]], x_mat[i])
                else:
                    cur_chain_xyz = x_mat[i]

            if i in self.cal_class_rmsf:
                cur_rmsf_results = np.array([self.cal_rmsf(x, self.init_pos, domain=self.domain) for x in cur_chain_xyz])
                cur_rmsf_results = np.sqrt(np.mean(cur_rmsf_results, axis=0))
                self.rmsf_results[i] += cur_rmsf_results
        # print(self.rmsf_results)


    def calculate(self):
        print(f"您当前的分子类型有：\n{self.data.molecules}")
        self.cal_class_rmsf = list(
            map(lambda x: int(x) - 1, input(f"请输入想要计算 RMSF 的分子序号, 以 ',' 分隔: ").split(',')))
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
            self.rmsf_results = {i: np.zeros(self.data.length_dict[i]) for i in self.cal_class_rmsf}
            self.domain = None

        if not self.balance_cut:
            self.balance_cut = input(
                "请输入需要截取的平衡后的文件索引，索引从 1 开始，格式为 'START-END', 例如：1000-2000，直接回车则不截取：")
        if not self.balance_cut:
            chain_files = os.listdir(self.chain_path)
        else:
            start, end = list(map(int, self.balance_cut.split('-')))
            chain_files = os.listdir(self.chain_path)[start - 1: end]

        # 使用 tqdm 包装可迭代对象以显示进度条
        list(tqdm(map(self.process_chain_file, sorted(chain_files)),
                  total=len(chain_files),
                  desc="calculating RMSF",
                  colour='cyan',
                  bar_format = '{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]',
                  ncols=100))

        # 保存结果
        result_file = os.path.join(self.save_path, f"draw_RMSF_ref_{os.path.basename(self.ref)}.log")
        with open(result_file, 'w') as f:
            f.write(str({key: list(self.rmsf_results[key]) for key in self.rmsf_results.keys()}))
        print(f"RMSF 计算完成！结果已保存至文件 draw_RMSF_ref_{os.path.basename(self.ref)}.log")

        self.draw_rmsf()


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
                print(
                    f"RMSF 绘图完成！结果已保存至文件 draw_RMSF_{mol}_ref_{os.path.basename(self.ref)}.png")
