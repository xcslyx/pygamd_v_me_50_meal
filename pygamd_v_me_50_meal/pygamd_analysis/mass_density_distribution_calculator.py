import os
import shutil

import numpy as np
import scipy as sp
import torch as torch
import multiprocessing as mp
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET

from pygamd_v_me_50_meal.Functions import Functions

# 计算径向质量数密度分布的类
class MassDensityDistributionCalculator:
    def __init__(self, path, data):
        self.path = path
        self.data = data

        self.mol_class_dict = self.data.mol_class_dict
        self.length_dict = self.data.length_dict

        self.init_xml_path = os.path.join(self.path, "xml/")

        remove_pbc_choice = input("是否需要使用移除 PBC后的文件？(y/n)")
        if remove_pbc_choice == 'y':
            self.chain_path = os.path.join(self.path, "chain_xyz_remove_pbc_condensate/")
        else:
            self.chain_path = os.path.join(self.path, "chain_xyz/")

        self.save_path = os.path.join(self.path, "draw_log/")
        os.makedirs(self.save_path, exist_ok=True)

        self.free_chain_save_path = os.path.join(self.path, "free_chain/")
        self.free_chain_choice = input("是否需要排除游离链? (y/n)")
        if self.free_chain_choice == "y":
            if not os.path.exists(self.free_chain_save_path):
                os.makedirs(self.free_chain_save_path, exist_ok=True)
                print(f"✅ 已创建 {self.free_chain_save_path} 目录")

            else:
                self.cover_free_chain_save_path_choice = input(f"已存在 {self.free_chain_save_path} 目录，是否覆盖？(y/n)")
                if self.cover_free_chain_save_path_choice == 'y':
                    shutil.rmtree(self.free_chain_save_path)
                    os.makedirs(self.free_chain_save_path, exist_ok=True)
                    print(f"✅ 已创建 {self.free_chain_save_path} 目录")

        self.balance_cut = ""
        self.mass_density_path = os.path.join(self.save_path, "mass_density/")
        if not os.path.exists(self.mass_density_path):
            os.makedirs(self.mass_density_path, exist_ok=True)
            print(f"✅ 已创建 {self.mass_density_path} 目录")

        self.cal_mass_density_distribution_list = []
        if not self.cal_mass_density_distribution_list:
            print(f"\n您的分子类型有：\n{self.data.molecules}")
            self.cal_mass_density_distribution_list = input("请输入您想要计算质量密度分布时需要包括的分子，以逗号分隔，如“1,2”。\n"
                                                            "如需计算全部分子，请输入 all 或直接回车：").split(',')
        if "all" in self.cal_mass_density_distribution_list or self.cal_mass_density_distribution_list == [""]:
            self.cal_mass_density_distribution_list = list(self.mol_class_dict.keys())
        else:
            self.cal_mass_density_distribution_list = [self.data.mol_class_list[int(i)] for i in self.cal_mass_density_distribution_list]
        print(f"即将计算质量密度分布的分子：{self.cal_mass_density_distribution_list}")


        self.sequence = {}
        with open(f"{os.path.join(self.path, self.data.system_name)}_sequence.txt") as f:
            sequence = eval(f.read())
            for cal_mol in self.cal_mass_density_distribution_list:
                self.sequence[cal_mol] = [i[1] for i in sequence[cal_mol]]

        self.domain = None
        self.new_name = None
        domain = input(
            "若只需要计算结构域，请输入该结构域的起始残基编号和末尾残基编号（从 1 开始），以-分隔，如 159-522，否则请直接回车。"
            "注：只支持选择单个分子计算。"
            "请输入：")
        if domain:
            domain = list(map(int, domain.split('-')))
            self.domain = domain
            print(f"即将计算结构域：{domain}")
            self.new_name = input("是否给此结构域命名？请输入名称，否则直接回车：")
            if self.new_name is not None:
                self.mass_density_save_path = os.path.join(self.mass_density_path, f"{self.new_name}/")
            else:
                self.mass_density_save_path = os.path.join(self.mass_density_path, cal_mol)

            if not os.path.exists(self.mass_density_save_path):
                os.makedirs(self.mass_density_save_path, exist_ok=True)
                print(f"✅ 已创建 {self.mass_density_save_path} 目录")
            elif file_args.avg != "mass_density":
                shutil.rmtree(self.mass_density_save_path)
                os.makedirs(self.mass_density_save_path, exist_ok=True)

            self.res_file = os.path.join(self.mass_density_path, f"mass_density_of_{self.new_name}.dat")
        else:
            for cal_mol in self.cal_mass_density_distribution_list:
                mass_density_save_path = os.path.join(self.mass_density_path, cal_mol)
                if not os.path.exists(mass_density_save_path):
                    os.makedirs(mass_density_save_path, exist_ok=True)
                    print(f"✅ 已创建 {mass_density_save_path} 目录")
                elif file_args.avg != "mass_density":
                    shutil.rmtree(mass_density_save_path)
                    os.makedirs(mass_density_save_path, exist_ok=True)
                    print(f"✅ 已重置 {mass_density_save_path} 目录")

        if not self.balance_cut:
            self.balance_cut = input(
                "请输入需要截取的平衡后的文件索引，格式为‘开始,结束’，例如：1000,2000，直接回车则不截取：")
        if not self.balance_cut:
            self.files = sorted(os.listdir(self.chain_path))
        else:
            start, end = list(map(int, self.balance_cut.split(',')))
            self.files = sorted(os.listdir(self.chain_path))[start: end + 1]

        # 选择计算哪一种 MSD（sphere/axis/slab）
        while True:
            print("可计算的 MSD 类型有：sphere/axis/slab")
            self.mass_density_choice = input("请输入需要计算的 MSD 类型，直接回车则使用默认值：")
            if not self.mass_density_choice or self.mass_density_choice == "sphere":
                self.mass_density_class = "sphere"
                break
            elif self.mass_density_choice == "axis":
                self.mass_density_class = "axis"
                break
            elif self.mass_density_choice == "slab":
                self.mass_density_class = "slab"
                break
            else:
                print("输入错误，请重新输入！")

        xml_path = os.path.join(self.path, "xml")
        self.xml_files = sorted(os.listdir(xml_path))
        init_xml_file = os.path.join(self.path, "xml", self.xml_files[0])
        init_root = ET.parse(init_xml_file).getroot()
        self.box_size = [float(init_root.find('.//box').attrib[i]) for i in ['lx', 'ly', 'lz']]

        self.dr = 2.5  # 计算质量密度分布时使用的 bin 大小

        if self.mass_density_class == "slab":
            self.r_max = self.box_size[2]
        else:
            self.r_max = self.box_size[2] / 2
        print(f"计算质量密度分布时使用的 bin 大小为 {self.dr} nm")
        # print(f"计算质量密度分布的边界为 {self.r_max} nm")

        # 初始化壳层数量
        self.num_bins = int(self.r_max / self.dr)

        if mp.get_start_method(allow_none=True) is None:
            mp.set_start_method('spawn')
        self.device = None


    def abstract_free_chain(self, name):
        x_mat = eval(open(os.path.join(self.chain_path, name), 'r').read())
        if self.free_chain_choice == "y":
            free_chain_dict = dict(zip(self.mol_class_dict.keys(), [[] for _ in range(len(self.mol_class_dict))]))
        with open(self.free_chain_save_path + name, 'w') as f:
            all_mol = []
            for mol in x_mat:
                all_mol.extend(x_mat[mol])

            mol_num = len(all_mol)
            free_chain_list = [[] for _ in range(mol_num)]
            for chain_idx_i in range(mol_num):
                for chain_idx_j in range(chain_idx_i + 1, mol_num):
                    a = torch.tensor(all_mol[chain_idx_i], device=self.device)
                    b = torch.tensor(all_mol[chain_idx_j], device=self.device)
                    c = Functions.euclidean_distances(a, b)
                    # c = Functions.euclidean_distances(all_mol[chain_idx_i], all_mol[chain_idx_j])
                    if c.min() < 0.8:
                        free_chain_list[chain_idx_i].append(chain_idx_j)
                        free_chain_list[chain_idx_j].append(chain_idx_i)

            contact_list = list(map(set, free_chain_list))
            length = len(contact_list)
            for i in range(1, length):
                for j in range(i):
                    if contact_list[i] == {0} or contact_list[j] == {0}:
                        continue
                    x = contact_list[i].union(contact_list[j])
                    y = len(contact_list[i]) + len(contact_list[j])
                    if len(x) < y:
                        contact_list[i] = x
                        contact_list[j] = {0}
            contact_list = [i for i in contact_list if i != {0}]

            contact_chain = list(max(contact_list, key=len))
            single_chain = list(set(list(range(mol_num))) - set(contact_chain))

            for i in range(len(single_chain)):
                cur_chain_idx = single_chain[i]
                for mol in self.mol_class_dict:
                    if cur_chain_idx >= 0 > cur_chain_idx - self.mol_class_dict[mol][0]:
                        break
                    else:
                        cur_chain_idx -= self.mol_class_dict[mol][0]

                free_chain_dict[mol].append(cur_chain_idx)

            f.write(str(free_chain_dict))


    def cal_mass_density_distribution(self, name):
        x_mat = eval(open(os.path.join(self.chain_path, name), 'r').read())
        if self.free_chain_choice:
            free_chain_dict = eval(open(os.path.join(self.free_chain_save_path, name), 'r').read())
            center, _ = Functions.abstract_centroid(os.path.join(self.chain_path, name),
                                                    self.cal_mass_density_distribution_list,
                                                    free_chain_dict=free_chain_dict,)
        else:
            center, _ = Functions.abstract_centroid(os.path.join(self.chain_path, name),
                                                    self.cal_mass_density_distribution_list,)

        # 存储中心坐标
        with open(os.path.join(self.mass_density_path, "center.txt"), 'w') as f:
            f.write(str(list(center)))

        for cal_mol in self.cal_mass_density_distribution_list:
            # if self.free_chain_choice == "y":
            #     free_chains = free_chain_dict[cal_mol]
            molecules = []
            for chain_idx in range(len(x_mat[cal_mol])):
                # if self.free_chain_choice == "y":
                    # if chain_idx in free_chains:
                    #     continue
                chain = x_mat[cal_mol][chain_idx]
                if self.domain:
                    chain = chain[self.domain[0] - 1:self.domain[1]]
                else:
                    chain = chain
                molecules.append(chain)
            # print(f"计算质量密度分布的分子数：{len(molecules)}")

            # 计算质量密度分布，molecules 列表中储存了链的坐标，每个链的坐标是一个二维列表，每个元素为一个三维坐标
            # 初始化质量和体积
            shell_count = np.zeros(self.num_bins)
            shell_volume = np.zeros(self.num_bins)

            # 遍历每条链
            for chain in molecules:
                for coord_idx in range(len(chain)):
                    # 计算到中心的径向距离
                    r = np.linalg.norm(np.array(chain[coord_idx]) - center)

                    # 确定壳层索引
                    shell_index = int(r / self.dr)

                    # 累加粒子质量，单位是 g/mol
                    if shell_index < self.num_bins:
                        shell_count[shell_index] += self.sequence[cal_mol][coord_idx // self.length_dict[cal_mol]]

            # 计算壳层体积，单位是 mL，即 cm^3
            for i in range(self.num_bins):
                r_inner = i * self.dr
                r_outer = (i + 1) * self.dr
                shell_volume[i] = 4 / 3 * np.pi * (r_outer ** 3 - r_inner ** 3) * 1e-21   # 单位为 mL，即 cm^3

            # 计算密度（质量 / 体积），换算单位为 mg/mL
            density = shell_count / shell_volume / sp.constants.N_A * 1000  # 单位为 mg/mL

            # 计算每个壳层的半径中心
            radii = np.linspace(self.dr / 2, self.r_max - self.dr / 2, self.num_bins)

            # 保存结果
            if self.new_name:
                with open(os.path.join(self.mass_density_save_path, name), 'w') as f:
                    for r, d in zip(radii, density):
                        f.write(f"{r:20.4f}{d:20.6f}\n")
            else:
                with open(os.path.join(self.mass_density_path, cal_mol, name), 'w') as f:
                    for r, d in zip(radii, density):
                        f.write(f"{r:20.4f}{d:20.6f}\n")


    def cal_mass_density_distribution_axis(self, name):
        x_mat = eval(open(os.path.join(self.chain_path, name), 'r').read())

        # 选择一个参考点作为原点
        try:
            free_chain_dict = eval(open(os.path.join(self.free_chain_save_path, name), 'r').read())
            center, _, height = Functions.abstract_centroid(os.path.join(self.chain_path, name),
                                                            self.cal_mass_density_distribution_list,
                                                            free_chain_dict=free_chain_dict,
                                                            get_height=True)
        except:
            center, _, height = Functions.abstract_centroid(os.path.join(self.chain_path, name),
                                                            self.cal_mass_density_distribution_list,
                                                            get_height=True)

        # 存储质心坐标
        with open(os.path.join(self.mass_density_path, "center.txt"), 'w') as f:
            f.write(str(list(center)))

        for cal_mol in self.cal_mass_density_distribution_list:
            # free_chains = free_chain_dict[cal_mol]
            molecules = []
            for chain_idx in range(len(x_mat[cal_mol])):
                # if chain_idx in free_chains:
                #     continue
                chain = x_mat[cal_mol][chain_idx]
                if self.domain:
                    chain = chain[self.domain[0]-1:self.domain[1]]
                else:
                    chain = chain
                molecules.append(chain)
            # print(f"计算质量密度分布的分子数：{len(molecules)}")

            # 展平所有坐标为一个二维数组
            all_coords = np.concatenate(molecules, axis=0)
            relative_positions = all_coords - center  # N x 3 array

            pca = PCA(n_components=1)
            pca.fit(all_coords)
            axis_vector = pca.components_[0]

            # 储存 axis_vector
            # if self.new_name:
            #     with open(os.path.join(self.mass_density_save_path, name.replace(".xml", "_axis_vector.txt")), 'w') as f:
            #         f.write(str(list(axis_vector)))
            # else:
            #     with open(os.path.join(self.mass_density_path, cal_mol, name.replace(".xml", "_axis_vector.txt")), 'w') as f:
            #         f.write(str(list(axis_vector)))

            # 计算把主轴旋转到 z 轴所需的旋转矩阵，利用旋转矩阵将所有坐标旋转到 z 轴
            R = Functions.rodrigues_rotation(axis_vector)  # 3 x 3 array
            rotated_positions = relative_positions @ R.T  # N x 3 array

            # 储存旋转后的坐标
            # if self.new_name:
            #     Functions.replace_position(os.path.join(self.init_xml_path, name), rotated_positions,
            #                                new_xml_file=os.path.join(self.mass_density_save_path, name.replace(".xml", "_rotated.xml")))
            # else:
            #     Functions.replace_position(os.path.join(self.init_xml_path, name), rotated_positions,
            #                                new_xml_file=os.path.join(self.mass_density_path, cal_mol, name.replace(".xml", "_rotated.xml")))

            # 计算质量密度分布
            # 初始化质量和体积
            shell_count = np.zeros(self.num_bins)
            shell_volume = np.zeros(self.num_bins)

            # 遍历每条链
            for coord_idx in range(len(rotated_positions)):
                # 计算到 z 轴的径向距离，即 x^2 + y^2
                r = np.linalg.norm(rotated_positions[coord_idx][:2])
                # print(r)

                # 确定壳层索引
                shell_index = int(r / self.dr)

                # 累加粒子质量，单位是 g/mol
                if shell_index < self.num_bins:
                    shell_count[shell_index] += self.sequence[cal_mol][coord_idx // self.length_dict[cal_mol]]
                    # print(shell_count[shell_index])

            # 统计质量
            # 计算壳层体积，单位是 mL
            for i in range(self.num_bins):
                r_inner = i * self.dr
                r_outer = (i + 1) * self.dr

                shell_volume[i] = height * np.pi * (r_outer ** 2 - r_inner ** 2) * 1e-21  # 单位为 mL

            # 在 shell_volume 中，如果数值为 0，则说明该壳层没有粒子，则将其体积设置为 1，以避免除 0 错误
            shell_volume[shell_volume == 0] = 1
            # 计算密度（质量 / 体积），换算单位为 mg/mL
            # print(f"shell_count: {shell_count}")
            # print(f"shell_volume: {shell_volume}")
            density = shell_count / shell_volume / sp.constants.N_A * 1000  # 单位为 mg/mL
            # density = shell_count / sp.constants.N_A / 1e-21 * 1000  # 单位为 mg/mL

            # 计算每个壳层的半径中心
            radii = np.linspace(self.dr / 2, self.r_max - self.dr / 2, self.num_bins)

            if self.new_name:
                with open(os.path.join(self.mass_density_save_path, name.replace(".xml", "_axis.xml")), 'w') as f:
                    for r, d in zip(radii, density):
                        f.write(f"{r:20.4f}{d:20.6f}\n")
            else:
                with open(os.path.join(self.mass_density_path, cal_mol, name.replace(".xml", "_axis.xml")), 'w') as f:
                    for r, d in zip(radii, density):
                        f.write(f"{r:20.4f}{d:20.6f}\n")


    def cal_mass_density_distribution_slab(self, name):
        x_mat = eval(open(os.path.join(self.chain_path, name), 'r').read())
        # free_chain_dict = eval(open(os.path.join(self.free_chain_save_path, name), 'r').read())

        # 从 z 轴最小值开始，每次向上移动 dr，计算每一层的质量密度
        for cal_mol in self.cal_mass_density_distribution_list:
            # free_chains = free_chain_dict[cal_mol]
            molecules = []
            for chain_idx in range(len(x_mat[cal_mol])):
                # if chain_idx in free_chains:
                #     continue
                chain = x_mat[cal_mol][chain_idx]
                if self.domain:
                    chain = chain[self.domain[0]-1:self.domain[1]]
                else:
                    chain = chain
                molecules.append(chain)
            # print(f"计算质量密度分布的分子数：{len(molecules)}")

            # 展平所有坐标为一个二维数组
            all_coords = np.concatenate(molecules, axis=0)
            # 减去 z 坐标均值，使粒子中心在 0 处
            all_coords[:, 2] -= all_coords[:, 2].mean()


            # 计算质量密度分布
            # 初始化质量和体积
            shell_count = np.zeros(self.num_bins)
            shell_volume = np.zeros(self.num_bins)

            # 遍历每个粒子
            for coord_idx in range(len(all_coords)):
                # 确定壳层索引
                z = all_coords[coord_idx][2] + self.box_size[2] / 2
                shell_index = int(z / self.dr)

                # 累加粒子质量，单位是 g/mol
                if shell_index < self.num_bins:
                    shell_count[shell_index] += self.sequence[cal_mol][coord_idx // self.length_dict[cal_mol]]

            # 统计质量
            # 计算壳层体积，单位是 mL
            shell_volume = self.box_size[0] * self.box_size[1] * self.dr * 1e-21  # 单位为 mL

            # 计算密度（质量 / 体积），换算单位为 mg/mL
            density = shell_count / shell_volume / sp.constants.N_A * 1000  # 单位为 mg/mL

            # 计算每个壳层的半径中心
            radii = np.linspace(self.dr / 2, self.box_size[2] - self.dr / 2, self.num_bins)

            if self.new_name:
                with open(os.path.join(self.mass_density_save_path, name.replace(".xml", "_slab.xml")), 'w') as f:
                    for r, d in zip(radii, density):
                        f.write(f"{r:20.4f}{d:20.6f}\n")
            else:
                with open(os.path.join(self.mass_density_path, cal_mol, name.replace(".xml", "_slab.xml")), 'w') as f:
                    for r, d in zip(radii, density):
                        f.write(f"{r:20.4f}{d:20.6f}\n")



    def average_mass_density_distribution(self):
        if self.mass_density_class == "sphere":
            for cal_mol in self.cal_mass_density_distribution_list:
                cur_density_dict = {}
                for file in self.files:
                    if self.new_name:
                        cur_result = os.path.join(self.mass_density_save_path, file)
                    else:
                        cur_result = os.path.join(self.mass_density_path, cal_mol, file)

                    with open(cur_result, 'r') as f:
                        f_lines = f.readlines()
                        for line in f_lines:
                            if line.strip():
                                r, d = line.strip().split()
                                r, d = float(r), float(d)
                                if r not in cur_density_dict:
                                    cur_density_dict[r] = 0
                                cur_density_dict[r] += d

                # 计算平均密度
                for r in cur_density_dict.keys():
                    cur_density_dict[r] /= len(self.files)
                cur_density_dict = sorted(cur_density_dict.items())

                if self.new_name:
                    self.res_file = os.path.join(self.mass_density_path, f"mass_density_of_{self.new_name}.dat")
                else:
                    self.res_file = os.path.join(self.mass_density_path, f"mass_density_of_{cal_mol}.dat")
                with open(self.res_file, 'w') as f:
                    for r, d in cur_density_dict:
                        f.write(f"{r:20.4f}{d:20.6f}\n")
        elif self.mass_density_class == "axis":
            for cal_mol in self.cal_mass_density_distribution_list:
                cur_density_dict = {}
                for file in self.files:
                    if self.new_name:
                        cur_result = os.path.join(self.mass_density_save_path, file.replace(".xml", "_axis.xml"))
                    else:
                        cur_result = os.path.join(self.mass_density_path, cal_mol, file.replace(".xml", "_axis.xml"))

                    with open(cur_result, 'r') as f:
                        f_lines = f.readlines()
                        for line in f_lines:
                            if line.strip():
                                r, d = line.strip().split()
                                r, d = float(r), float(d)
                                if r not in cur_density_dict:
                                    cur_density_dict[r] = 0
                                cur_density_dict[r] += d

                # 计算平均密度
                for r in cur_density_dict.keys():
                    cur_density_dict[r] /= len(self.files)
                cur_density_dict = sorted(cur_density_dict.items())

                if self.new_name:
                    self.res_file = os.path.join(self.mass_density_path, f"mass_density_of_{self.new_name}_axis.dat")
                else:
                    self.res_file = os.path.join(self.mass_density_path, f"mass_density_of_{cal_mol}_axis.dat")
                with open(self.res_file, 'w') as f:
                    for r, d in cur_density_dict:
                        f.write(f"{r:20.4f}{d:20.6f}\n")
        elif self.mass_density_class == "slab":
            for cal_mol in self.cal_mass_density_distribution_list:
                cur_density_dict = {}
                for file in self.files:
                    if self.new_name:
                        cur_result = os.path.join(self.mass_density_save_path, file.replace(".xml", "_slab.xml"))
                    else:
                        cur_result = os.path.join(self.mass_density_path, cal_mol, file.replace(".xml", "_slab.xml"))

                    with open(cur_result, 'r') as f:
                        f_lines = f.readlines()
                        for line in f_lines:
                            if line.strip():
                                r, d = line.strip().split()
                                r, d = float(r), float(d)
                                if r not in cur_density_dict:
                                    cur_density_dict[r] = 0
                                cur_density_dict[r] += d

                # 计算平均密度
                for r in cur_density_dict.keys():
                    cur_density_dict[r] /= len(self.files)
                cur_density_dict = sorted(cur_density_dict.items())

                if self.new_name:
                    self.res_file = os.path.join(self.mass_density_path, f"mass_density_of_{self.new_name}_slab.dat")
                else:
                    self.res_file = os.path.join(self.mass_density_path, f"mass_density_of_{cal_mol}_slab.dat")
                with open(self.res_file, 'w') as f:
                    for r, d in cur_density_dict:
                        f.write(f"{r:20.4f}{d:20.6f}\n")


    def cal_mass_density_distribution_parallel(self, ):
        if self.free_chain_choice == "y":
            if self.cover_free_chain_save_path_choice == "y":
                if torch.cuda.is_available():
                    gpu_choice = input("即将使用 GPU 加速提取游离链，请指定 GPU 编号，或直接回车使用默认 0 号 GPU。若想使用 CPU，请输入 CPU：")
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
                with Pool(processes=4) as pool:
                    list(tqdm(pool.imap(self.abstract_free_chain, self.files),
                            total=len(self.files),
                            desc="计算中",
                            colour='cyan',
                            bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]',
                            ncols=100))



        with Pool(processes=4) as pool:
            if self.mass_density_class == "sphere":
                list(tqdm(pool.imap(self.cal_mass_density_distribution, self.files),
                               total=len(self.files),
                               desc="计算中",
                               colour='cyan',
                               bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]',
                               ncols=100))
            elif self.mass_density_class == "axis":
                list(tqdm(pool.imap(self.cal_mass_density_distribution_axis, self.files),
                               total=len(self.files),
                               desc="计算中",
                               colour='cyan',
                               bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]',
                               ncols=100))
            elif self.mass_density_class == "slab":
                list(tqdm(pool.imap(self.cal_mass_density_distribution_slab, self.files),
                               total=len(self.files),
                               desc="计算中",
                               colour='cyan',
                               bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]',
                               ncols=100))

        self.average_mass_density_distribution()
        self.draw_mass_density_distribution()


    def draw_mass_density_distribution(self):
        for cal_mol in self.cal_mass_density_distribution_list:
            cur_density_dict = []
            mol_name = self.new_name if self.new_name else cal_mol
            if self.mass_density_class == "sphere":
                cur_result = os.path.join(self.mass_density_path, f"mass_density_of_{mol_name}.dat")
                save_name = os.path.join(self.mass_density_path, f"mass_density_of_{mol_name}.png")
            elif self.mass_density_class == "axis":
                cur_result = os.path.join(self.mass_density_path, f"mass_density_of_{mol_name}_axis.dat")
                save_name = os.path.join(self.mass_density_path, f"mass_density_of_{mol_name}_axis.png")
            elif self.mass_density_class == "slab":
                cur_result = os.path.join(self.mass_density_path, f"mass_density_of_{mol_name}_slab.dat")
                save_name = os.path.join(self.mass_density_path, f"mass_density_of_{mol_name}_slab.png")
            else:
                raise ValueError("未知的质量密度分布类型！")

            with open(cur_result, 'r') as f:
                f_lines = f.readlines()
                for line in f_lines:
                    if line.strip():
                        r, d = line.strip().split()
                        r, d = float(r), float(d)
                        cur_density_dict.append([r, d])

            fig, ax = plt.subplots(figsize=(12, 9), dpi=300)
            ax.plot(np.array(cur_density_dict)[:, 0], np.array(cur_density_dict)[:, 1], label=f"{mol_name}")
            ax.set_xlabel('radius (nm)')
            ax.set_ylabel(r"$\rho$ (mg/mL)")
            ax.set_title(f'Mass density distribution of {mol_name}')
            ax.legend()

            plt.savefig(save_name)
