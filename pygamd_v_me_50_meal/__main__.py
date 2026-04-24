#!-*- coding: utf-8 -*-
import os
import re
import json
import shutil
import logging
import argparse

import matplotlib.pyplot as plt

import pygamd_v_me_50_meal as p50
from pygamd_v_me_50_meal.utils import str2value


# 设置 matplotlib 参数
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

print(os.path.dirname(__file__))
# 加载消息文件
with open(os.path.join(os.path.dirname(__file__), 'massage.json'), 'r', encoding='utf-8') as f:
    messages = json.load(f)['main_massage']


def main():
    print(f"New version notification {p50.__version__}.")
    print("Now you can use command v50_en to run the package, which is in English.")
    run_main('zh')

def main_en():
    run_main('en')

def run_main(lang):
    msg = messages[lang]

    print(msg['checking_updates'])
    try:
        from pygamd_v_me_50_meal.version_check import check_update
        check_update()
    except Exception as e:
        print(e)
    print(msg['running_package'])

    parser = argparse.ArgumentParser(
        # prog=f'{os.path.basename(__file__)} v0.0.20 增强版',
        # description='Do something you want to do in your system.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)} v{p50.__version__}')

    parser.add_argument('-p', '--path', metavar="/path/to/system",
                        type=str, default=None, help='体系目录路径.' if lang == 'zh' else 'System directory path.')

    parser.add_argument('-o', "--output",
                        type=str, help="输出文件名称" if lang == 'zh' else "Output file name")

    parser.add_argument('-box_size', metavar="box_size",
                        type=float, default=100.0, help="用于 pdb-xml 转换的盒子大小, 默认为 100.0 nm." if lang == 'zh' else "Box size for pdb-xml conversion, default 100.0 nm.")

    parser.add_argument("-pdb2xml", metavar="/path/to/pdb_file or filename.pdb with -p /path/to/system.",
                        type=str, default=None, help="将 PDB 文件转换供 PYGAMD 模拟的 XML 文件。" if lang == 'zh' else "Convert PDB file to XML file for PYGAMD simulation.")

    parser.add_argument("-add_enm_bond", metavar="enm_domain_list",
                        type=str2value, default=None,
                        help="若要设置弹性网络，请输入结构域的起始残基编号和末尾残基编号（从 1 开始），以-分隔，如 159-522，若有多个结构域，请以英文逗号分隔。" if lang == 'zh' else "To set up elastic network, enter the starting and ending residue numbers of domains (starting from 1), separated by -, e.g., 159-522. For multiple domains, separate with commas.")

    parser.add_argument("-add_rigid_body", metavar="rigid_body_domain_list",
                        type=str2value, default=None,
                        help="若要设置刚体，请输入结构域的起始残基编号和末尾残基编号（从 1 开始），以-分隔，如 159-522，若有多个结构域，请以英文逗号分隔。" if lang == 'zh' else "To set up rigid body, enter the starting and ending residue numbers of domains (starting from 1), separated by -, e.g., 159-522. For multiple domains, separate with commas.")

    parser.add_argument('-add_domain', action='store_true', help="若要将结构域残基单独设置粒子类型（例如甘氨酸A→AD），请设置为 True，将以rigid body 或 enm bond 的结构域列表来设置单独的粒子类型。" if lang == 'zh' else "To set separate particle types for domain residues (e.g., glycine A→AD), set to True. It will use rigid body or enm bond domain lists to set separate particle types.")

    parser.add_argument("-dna_model", metavar="DNA model",
                        type=str, default="unset", help="设置 DNA 蛋白质模型，可选：略." if lang == 'zh' else "Set DNA protein model, options: see details.")

    parser.add_argument("-gen_run_file", action='store_true', help="是否生成 PYGAMD 动力学模拟脚本." if lang == 'zh' else "Whether to generate PYGAMD dynamics simulation script.")

    parser.add_argument('-xyz', action='store_true', help="是否提取坐标文件." if lang == 'zh' else "Whether to extract coordinate files.")

    parser.add_argument('-remove_enm', action='store_true', help="是否移除弹性键。" if lang == 'zh' else "Whether to remove elastic bonds.")

    parser.add_argument('-remove_condensate_pbc', action='store_true', help="是否移除凝聚体的 PBC。" if lang == 'zh' else "Whether to remove PBC and move the largest condensate to the center of the box.")

    parser.add_argument('-remove_ions_zhy', action='store_true', help="是否移除 xml 文件中的离子。" if lang == 'zh' else "Whether to remove ions from xml files.")

    parser.add_argument('-cm', action='store_true', help="是否计算接触图, 计算后会自动绘图。" if lang == 'zh' else "Whether to calculate contact map, will automatically plot after calculation.")

    parser.add_argument('-em', action='store_true', help="是否计算 Energy Matrix, 计算后会自动绘图。" if lang == 'zh' else "Whether to calculate Energy Matrix, will automatically plot after calculation.")

    parser.add_argument('-r_cut', metavar="r_cut of contact map",
                        type=float, default=4.0, help="计算接触图的截断半径，默认值为 4.0 Å。" if lang == 'zh' else "Cutoff radius for contact map calculation, default value is 4.0 Å.")

    parser.add_argument('-draw', metavar="想要绘制图像的类型, 如cm,rmsd,rmsf。",
                        type=str, default=None, help="用于无需计算的情况下绘图。" if lang == 'zh' else "For drawing without calculation.")

    parser.add_argument('-cm_choice', metavar="分子组合/轨迹切片",
                        type=str, default="/", help="计算接触图的选择，如“0-0,1-1/1000,2000”。" if lang == 'zh' else "Contact map calculation options, e.g., '0-0,1-1/1000,2000'." +
                                                    "不提供某一项代表全选，如“0-0,1-1/”代表只计算0-0和1-1之间的接触图，选取所有轨迹。" if lang == 'zh' else "Omitting an item means selecting all, e.g., '0-0,1-1/' means only calculating contact maps between 0-0 and 1-1, selecting all trajectories." +
                                                    "若不提供，则会在运行中进行提示，此选项供 nohup 使用。" if lang == 'zh' else "If not provided, it will prompt during runtime, this option is for nohup use.")

    # TODO: 增加对接触图平均的功能
    parser.add_argument('-avg', metavar="计算类型",
                        type=str, default="unset", help="用于无需计算的情况系进行平均（还不好用）。" if lang == 'zh' else "For averaging without calculation (not fully functional yet).")

    parser.add_argument('-rg', action='store_true', help="是否计算 Rg." if lang == 'zh' else "Whether to calculate Rg.")

    parser.add_argument('-rmsd', action='store_true', help="是否计算 RMSD." if lang == 'zh' else "Whether to calculate RMSD.")

    parser.add_argument('-rmsf', action='store_true', help="是否计算 RMSF." if lang == 'zh' else "Whether to calculate RMSF.")

    parser.add_argument('-ref', metavar="/path/to/reference_structure.xml",
                        type=str, default=None, help="RMSD/RMSF 计算的参考结构文件路径。" if lang == 'zh' else "Reference structure file path for RMSD/RMSF calculation.")

    parser.add_argument('-get_seq', metavar="/path/to/file or filename with -p /path/to/system.",
                        type=str, default=None, help="获取 XML 文件或 PDB 文件的序列。" if lang == 'zh' else "Get sequence from XML file or PDB file.")

    parser.add_argument('-mass_density', action='store_true', help="是否计算质量密度分布。" if lang == 'zh' else "Whether to calculate mass density distribution.")

    parser.add_argument('-amino_acid', action='store_true', help="按照氨基酸类型计算质量密度分布。" if lang == 'zh' else "Calculate mass density distribution by amino acid type.")

    parser.add_argument('-msd', action='store_true', help="是否计算 MSD。" if lang == 'zh' else "Whether to calculate MSD.")

    parser.add_argument('-eed', action='store_true', help="是否计算末端距。" if lang == 'zh' else "Whether to calculate end-to-end distance.")

    parser.add_argument('-seq_analysis', action='store_true', help="是否进行序列分析。" if lang == 'zh' else "Whether to perform sequence analysis.")

    parser.add_argument('-seq', metavar="path/to/sequence_file or sequence_string",
                        type=str, default=None, help="序列文件路径或直接传入的序列字符串。" if lang == 'zh' else "Sequence file path or directly passed sequence string.")

    parser.add_argument('-seq_window', metavar="window_size",
                        type=int, default=15, help="序列分析的滑动窗口大小，默认 15。" if lang == 'zh' else "Sliding window size for sequence analysis, default 15.")

    parser.add_argument('-seq_output', metavar="output_file",
                        type=str, default=None, help="序列分析的输出文件路径。" if lang == 'zh' else "Output file path for sequence analysis.")

    parser.add_argument('-nc', action='store_true', help="是否计算网络聚类。" if lang == 'zh' else "Whether to calculate network clustering.")

    parser.add_argument('-nc_node', metavar="node_molecule_type",
                        type=str, default=None, help="网络聚类的节点分子类型，如 'cGAS'。" if lang == 'zh' else "Node molecule type for network clustering, e.g., 'cGAS'.")

    parser.add_argument('-nc_edge', metavar="edge_molecule_type",
                        type=str, default=None, help="网络聚类的边分子类型，如 'MED1'。" if lang == 'zh' else "Edge molecule type for network clustering, e.g., 'MED1'.")

    parser.add_argument('-nc_threshold', metavar="distance_threshold",
                        type=float, default=10.0, help="网络构建的距离阈值，默认 10.0 Å。" if lang == 'zh' else "Distance threshold for network construction, default 10.0 Å.")

    file_args = parser.parse_args()

    # 序列分析作为独立功能，不需要 -p 参数
    if file_args.seq_analysis:
        print(msg['start_seq_analysis'])
        if file_args.seq:
            # 动态导入 SequenceAnalyzer，避免其他模块的依赖问题
            from pygamd_v_me_50_meal.pygamd_analysis.sequence_analysis.sequence_analyzer import SequenceAnalyzer
            
            seq_input = file_args.seq
            
            # 检查是否是文件路径
            if os.path.exists(seq_input):
                # 是文件路径，读取文件内容
                import ast
                with open(seq_input, 'r') as f:
                    seq_content = f.read()
                    try:
                        seq = list(ast.literal_eval(seq_content))
                    except:
                        seq = list(seq_content.strip())
                base_name = os.path.splitext(os.path.basename(seq_input))[0]
            else:
                # 不是文件路径，直接作为序列字符串处理
                seq = list(seq_input.strip())
                base_name = "sequence"
            
            analyzer = SequenceAnalyzer(window=file_args.seq_window)
            
            # 交互式选择分析类型
            print(msg['select_analysis_type'])
            print(msg['seq_analysis'])
            choice = input(msg['enter_option'])
            
            if choice == "1":
                analysis_type = "ncpr"
                output_path = file_args.seq_output
                if not output_path:
                    output_path = f"{base_name}_NCPR.png"
                fig, ax, values = analyzer.plot_ncpr(seq, save_path=output_path)
                print(msg['ncpr_analysis_done'].format(output_path=output_path))
                print(msg['ncpr_values_shape'].format(shape=values.shape))
            elif choice == "2":
                analysis_type = "aromatic"
                output_path = file_args.seq_output
                if not output_path:
                    output_path = f"{base_name}_aromaticity.png"
                fig, ax, values = analyzer.plot_aromaticity(seq, save_path=output_path)
                print(msg['aromatic_analysis_done'].format(output_path=output_path))
                print(msg['aromatic_values_shape'].format(shape=values.shape))
            elif choice == "3":
                analysis_type = "hydrophobicity"
                output_path = file_args.seq_output
                if not output_path:
                    output_path = f"{base_name}_hydrophobicity.png"
                fig, ax, values = analyzer.plot_hydrophobicity(seq, save_path=output_path)
                print(msg['hydrophobicity_analysis_done'].format(output_path=output_path))
                print(msg['hydrophobicity_values_shape'].format(shape=values.shape))
            else:
                print(msg['error_invalid_option'])
                exit()
            exit()
        else:
            print(msg['error_no_seq'])
            exit()

    # 其他分析功能需要 -p 参数
    # current_dir_path = os.getcwd()
    path = str(file_args.path)

    if path is None:
        if file_args.get_seq is not None:
            path = str(os.path.dirname(file_args.get_seq))
        elif file_args.pdb2xml is not None:
            path = str(os.path.dirname(file_args.pdb2xml))
            if not path:
                path = os.getcwd()
        else:
            raise ValueError(msg['please_provide_path'])

    # path = os.path.join(current_dir_path, path)

    ref = file_args.ref

    if file_args.pdb2xml:
        from pygamd_v_me_50_meal.simulate_creation.xml_generator import XMLGenerator
        XMLGenerator(path, file_args.pdb2xml, file_args.box_size,
                     add_enm_bond=file_args.add_enm_bond, add_rigid_body=file_args.add_rigid_body,
                     add_domain=file_args.add_domain,
                     dna_model=file_args.dna_model, gen_run_file=file_args.gen_run_file)
        exit()

    from pygamd_v_me_50_meal.data import Data
    data = Data(path)

    if file_args.get_seq:
        from pygamd_v_me_50_meal.pygamd_analysis.sequence_extractor import GetSequence
        GetSequence(path, file_args.get_seq, data=data)
        exit()

    if file_args.xyz:
        from pygamd_v_me_50_meal.pygamd_analysis.coordinates_processor import CoordinatesProcessor
        CoordinatesProcessor(path, data, file_args.remove_ions_zhy).cal_xyz(remove_condensate_pbc=file_args.remove_condensate_pbc)

    if file_args.cm:
        from pygamd_v_me_50_meal.pygamd_analysis.contact_map_calculator import ContactMapCalculator
        print(msg['start_calculating_cm'])
        ContactMapCalculator(path,
                            data=data,
                            cm_choice=file_args.cm_choice,
                            r_cut=file_args.r_cut,
                            ).calculate_contact_map_parallel()
        exit()

    if file_args.em:
        from pygamd_v_me_50_meal.pygamd_analysis.energy_map_calculator import EnergyMapCalculator
        print(msg['start_calculating_em'])
        EnergyMapCalculator(path,
                             data=data,
                             em_choice=file_args.cm_choice,
                             r_cut=file_args.r_cut,
                             ).calculate_energy_map_parallel()


    if file_args.rg:
        from pygamd_v_me_50_meal.pygamd_analysis.rg_calculator import RgCalculator
        RgCalculator(path, data).calculate()
    if file_args.rmsd:
        from pygamd_v_me_50_meal.pygamd_analysis.rmsd_calculator import RMSDCalculator
        RMSDCalculator(path, data, ref).calculate()
    if file_args.rmsf:
        from pygamd_v_me_50_meal.pygamd_analysis.rmsf_calculator import RMSFCalculator
        RMSFCalculator(path, data, ref).calculate()

    if file_args.mass_density:
        from pygamd_v_me_50_meal.pygamd_analysis.mass_density_distribution_calculator import MassDensityDistributionCalculator
        MassDensityDistributionCalculator(path, data, lang).cal_mass_density_distribution_parallel(amino_acid=file_args.amino_acid)

    if file_args.msd:
        from pygamd_v_me_50_meal.pygamd_analysis.msd_calculator import MSDCalculator
        MSDCalculator(path, data).cal_msd_parallel()
        # MSDCalculator().cal_displacement_probability_distribution_parallel()

    if file_args.eed:
        from pygamd_v_me_50_meal.pygamd_analysis.end_to_end_distance_calculator import EndToEndDistanceCalculator
        EndToEndDistanceCalculator(path, data).cal_end_to_end_distance_parallel()

    if file_args.nc:
        from pygamd_v_me_50_meal.pygamd_analysis.network_cluster_calculator import NetworkClusterCalculator
        print(msg['start_calculating_nc'])
        if file_args.nc_node is None or file_args.nc_edge is None:
            print(msg['your_molecule_types'].format(data=data.molecules))
            node_idx = int(input(msg['enter_node_mol_idx'])) - 1
            edge_idx = int(input(msg['enter_edge_mol_idx'])) - 1
            node_mol = data.mol_class_list[node_idx-1]
            edge_mol = data.mol_class_list[edge_idx-1]
        else:
            node_mol = file_args.nc_node
            edge_mol = file_args.nc_edge
        
        NetworkClusterCalculator(
            path=path,
            data=data,
            node_mol_type=node_mol,
            edge_mol_type=edge_mol,
            distance_threshold=file_args.nc_threshold
        ).calculate_parallel()

    if file_args.draw:
        draw_class = file_args.draw.lower().split(",")
        for draw_type in draw_class:
            if draw_type == "rg":
                from pygamd_v_me_50_meal.pygamd_analysis.rg_calculator import RgCalculator
                RgCalculator(path, data).draw_rg_distribution()
            if draw_type == "rmsd":
                from pygamd_v_me_50_meal.pygamd_analysis.rmsd_calculator import RMSDCalculator
                RMSDCalculator(path, data, ref).draw_rmsd_distribution()


if __name__ == '__main__':
    main()