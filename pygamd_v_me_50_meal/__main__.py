#!-*- coding: utf-8 -*-
import os
import re
import shutil
import logging
import argparse

import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm

from pygamd_v_me_50_meal.data import Data
from pygamd_v_me_50_meal.utils import str2value

from pygamd_v_me_50_meal.simulate_creation.xml_generator import XMLGenerator
from pygamd_v_me_50_meal.pygamd_analysis.get_sequence import GetSequence
from pygamd_v_me_50_meal.pygamd_analysis.msd_calculator import MSDCalculator
from pygamd_v_me_50_meal.pygamd_analysis.coordinates_processor import CoordinatesProcessor
from pygamd_v_me_50_meal.pygamd_analysis.contact_map_calculator import ContactMapCalculator
from pygamd_v_me_50_meal.pygamd_analysis.rg_rmsd_rmsf_calculator import RgRMSDRMSFCalculator
from pygamd_v_me_50_meal.pygamd_analysis.mass_density_distribution_calculator import MassDensityDistributionCalculator


plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["axes.linewidth"] = .8
plt.rcParams["axes.labelsize"] = 22
plt.rcParams["axes.titlesize"] = 26
plt.rcParams["legend.fontsize"] = 22
plt.rcParams["xtick.minor.visible"] = True
plt.rcParams["ytick.minor.visible"] = True
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["xtick.labelsize"] = 22
plt.rcParams["ytick.labelsize"] = 22
plt.rcParams["xtick.top"] = True
plt.rcParams["ytick.right"] = True
plt.rcParams["axes.formatter.use_mathtext"] = True


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

        match = re.search(r'\d+', filename)
        if match:
            time_step = int(match.group())

        time_step = None
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

def main():
    parser = argparse.ArgumentParser(
        prog=f'{os.path.basename(__file__)} v0.0.20 增强版',
        description='Do something you want to do in your system.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)} v0.0.14')

    parser.add_argument('-p', '--path', metavar="/path/to/system",
                        type=str, default=None, help='体系目录路径.')

    parser.add_argument('-o', "--output",
                        type=str, help="输出文件名称")

    parser.add_argument('-box_size', metavar="box_size",
                        type=float, default=100.0, help="用于 pdb-xml 转换的盒子大小, 默认为 100.0 nm.")

    parser.add_argument("-pdb2xml", metavar="/path/to/pdb_file or filename.pdb with -p /path/to/system.",
                        type=str, default=None, help="将 PDB 文件转换供 GALAMOST 模拟的 XML 文件。")

    parser.add_argument("-add_enm_bond", metavar="enm_domain_list",
                        type=str2value, default=None,
                        help="若要设置弹性网络，请输入结构域的起始残基编号和末尾残基编号（从 1 开始），以-分隔，如 159-522，若有多个结构域，请以英文逗号分隔。")

    parser.add_argument("-add_rigid_body", metavar="rigid_body_domain_list",
                        type=str2value, default=None,
                        help="若要设置刚体，请输入结构域的起始残基编号和末尾残基编号（从 1 开始），以-分隔，如 159-522，若有多个结构域，请以英文逗号分隔。")

    parser.add_argument("-add_domain", metavar="domain_list",
                        type=str2value, default=None,
                        help="若要将结构域残基单独设置粒子类型（例如甘氨酸A→AD），请设置为 True，将以rigid body 或 enm bond 的结构域列表来设置单独的粒子类型。")

    parser.add_argument("-dna_model", metavar="DNA model",
                        type=str, default="unset", help="设置 DNA 蛋白质模型，可选：略.")

    parser.add_argument("-gen_run_file", metavar="T(rue)/F(alse)",
                        type=str2value, default="unset", help="是否生成 PYGAMD 动力学模拟脚本.")

    parser.add_argument('-xyz', metavar="T(rue)/F(alse)",
                        type=str2value, default="unset", help="是否提取坐标文件.")

    parser.add_argument('-remove_enm', metavar="T(rue)/F(alse)",
                        type=str2value, default="unset", help="是否移除弹性键。")

    parser.add_argument('-remove_condensate_pbc', metavar="T(rue)/F(alse)",
                        type=str2value, default="unset", help="是否移除凝聚体的 PBC。")

    parser.add_argument('-remove_ions_zhy', metavar="T(rue)/F(alse)",
                        type=str2value, default=False, help="是否移除 xml 文件中的离子。")

    parser.add_argument('-cm', metavar="T(rue)/F(alse)",
                        type=str2value, default="unset", help="是否计算接触图（contact map），计算后会自动绘图。")

    parser.add_argument('-r_cut', metavar="r_cut of contact map",
                        type=float, default=4.0, help="计算接触图的截断半径，默认值为 4.0 Å。")

    parser.add_argument('-draw', metavar="想要绘制图像的类型，如cm,rmsd,rmsf。",
                        type=str, default=None, help="用于无需计算的情况下绘图。")

    parser.add_argument('-cm_choice', metavar="分子组合/轨迹切片",
                        type=str, default="/", help="计算接触图的选择，如“0-0,1-1/1000,2000”。"
                                                    "不提供某一项代表全选，如“0-0,1-1/”代表只计算0-0和1-1之间的接触图，选取所有轨迹。"
                                                    "若不提供，则会在运行中进行提示，此选项供 nohup 使用。")

    # TODO: 增加对接触图平均的功能
    parser.add_argument('-avg', metavar="计算类型",
                        type=str, default="unset", help="用于无需计算的情况系进行平均（还不好用）。")

    parser.add_argument('-rg', metavar="T(rue)/F(alse)",
                        type=str2value, default="unset", help="是否计算 Rg.")

    parser.add_argument('-rmsd', metavar="T(rue)/F(alse)",
                        type=str2value, default="unset", help="是否计算 RMSD.")

    parser.add_argument('-rmsf', metavar="T(rue)/F(alse)",
                        type=str2value, default="unset", help="是否计算 RMSF.")

    parser.add_argument('-ref', metavar="/path/to/reference_structure.xml",
                        type=str, default=None, help="RMSD/RMSF 计算的参考结构文件路径。")

    parser.add_argument('-get_seq', metavar="/path/to/file or filename with -p /path/to/system.",
                        type=str, default=None, help="获取 XML 文件或 PDB 文件的序列。")

    parser.add_argument('-mass_density', metavar="T(rue)/F(alse)",
                        type=str2value, default="unset", help="是否计算质量密度分布。")

    parser.add_argument('-msd', metavar="T(rue)/F(alse)",
                        type=str2value, default="unset", help="是否计算 MSD。")

    parser.add_argument('-eed', metavar="T(rue)/F(alse)",
                        type=str2value, default="unset", help="是否计算末端距。")

    file_args = parser.parse_args()

    # current_dir_path = os.getcwd()
    path = file_args.path

    if path is None:
        if file_args.get_seq is not None:
            path = os.path.dirname(file_args.get_seq)

    if path is None:
        raise ValueError("请提供体系目录路径。")

    # path = os.path.join(current_dir_path, path)
    data = Data(path)

    ref = file_args.ref

    if file_args.pdb2xml:
        XMLGenerator(path, file_args.pdb2xml, file_args.box_size,
                     add_enm_bond=file_args.add_enm_bond, add_rigid_body=file_args.add_rigid_body,
                     add_domain=file_args.add_domain,
                     dna_model=file_args.dna_model)
        exit()

    if file_args.get_seq:
        GetSequence(path, file_args.get_seq, data=data)
        exit()

    if file_args.xyz:
        print("开始进行坐标提取...")
        CoordinatesProcessor(path, data, file_args.remove_ions_zhy, remove_condensate_pbc=file_args.remove_condensate_pbc)

    if file_args.cm:
        print("开始计算 contact map 文件...")
        ContactMapCalculator(path,
                             data=data,
                             cm_choice=file_args.cm_choice,
                             r_cut=file_args.r_cut,
                             ).calculate_contact_map_parallel()
        exit()

    cal_class_dict = {"Rg": file_args.rg, "RMSD": file_args.rmsd, "RMSF": file_args.rmsf}
    if True in cal_class_dict.values():
        RgRMSDRMSFCalculator(path, data, ref).calculate(cal_class_dict)
        exit()

    if file_args.mass_density:
        MassDensityDistributionCalculator(path, data).cal_mass_density_distribution_parallel()
        exit()

    if file_args.msd:
        MSDCalculator(path, data).cal_msd_parallel()
        # MSDCalculator().cal_displacement_probability_distribution_parallel()
        exit()

    if file_args.eed:
        EndToEndDistanceCalculator(path, data).cal_end_to_end_distance_parallel()
        exit()

if __name__ == '__main__':
    print(">>> Running pygamd_v_me_50_meal package...")
    main()





