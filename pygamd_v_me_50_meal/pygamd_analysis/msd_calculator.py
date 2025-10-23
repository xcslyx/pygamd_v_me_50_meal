import os
import shutil


class MSDCalculator:
    def __init__(self, path, data):
        self.path = path
        self.data = data
        
        self.mol_class_dict = self.data.mol_class_dict
        self.length_dict = self.data.length_dict

        remove_pbc_choice = input("是否需要使用移除 PBC后的文件？(y/n)")
        if remove_pbc_choice == 'y':
            self.chain_path = os.path.join(self.path, "chain_xyz_remove_pbc_condensate/")
        else:
            self.chain_path = os.path.join(self.path, "chain_xyz_unwrapping/")

        self.save_path = os.path.join(self.path, "draw_log/")
        os.makedirs(self.save_path, exist_ok=True)

        self.msd_path = os.path.join(self.save_path, "msd/")
        if not os.path.exists(self.msd_path):
            os.makedirs(self.msd_path, exist_ok=True)
        else:
            shutil.rmtree(self.msd_path)
            os.makedirs(self.msd_path, exist_ok=True)

        self.displacement_probability_distribution_path = os.path.join(self.save_path, "DisplacementProbabilityDistribution/")

        self.displacement_probability_distribution_bin_num = 100  # DPD 计算时使用的 bin 数量

        self.cal_msd_list = []
        if not self.cal_msd_list:
            print(f"\n您的分子类型有：\n{self.data.molecules}")
            self.cal_msd_list = input("请输入您想要计算 MSD 时需要包括的分子，以逗号分隔，如“1,2”。\n"
                                        "如需计算全部分子，请输入 all 或直接回车：").split(',')
        if "all" in self.cal_msd_list or self.cal_msd_list == [""]:
            self.cal_msd_list = list(self.data.mol_class_dict.keys())
        else:
            self.cal_msd_list = [self.data.mol_class_list[int(i)] for i in self.cal_msd_list]
        print(f"即将计算 MSD 的分子：{self.cal_msd_list}")

        self.sequence = []
        with open(f"{os.path.join(self.path, self.data.system_name)}_sequence.txt") as f:
            sequence = eval(f.read())
            for cal_mol in self.cal_msd_list:
                self.sequence.append(sequence[cal_mol])
        self.mass_list = [i[1]for i in self.sequence]

        self.balance_cut = ""
        if not self.balance_cut:
            self.balance_cut = input(
                "请输入需要截取的平衡后的文件索引，格式为‘开始,结束’，例如：1000,2000，直接回车则不截取：")
            if not self.balance_cut:
                self.files = sorted(os.listdir(self.chain_path))
            else:
                start, end = list(map(int, self.balance_cut.split(',')))
                self.files = sorted(os.listdir(self.chain_path))[start: end + 1]

        xml_path = os.path.join(self.path, "xml")
        self.xml_files = sorted(os.listdir(xml_path))
        init_xml_file = os.path.join(self.path, "xml", self.xml_files[0])
        root = ET.parse(init_xml_file).getroot()
        self.box_size = float(root.find('.//box').attrib['lx'])
        self.r_max = self.box_size / 2

        # self.manager = mp.Manager()
        self.centroid_dict: dict = {}
        self.condensate_centroid_dict: dict = {}

        self.abstract_centroid_parallel()

        self.msd_dict = dict(zip(self.cal_msd_list,
                                 [[[] for _ in range(self.data.mol_class_dict[cal_mol][0])] for cal_mol in
                                  self.cal_msd_list]))

        self.displacement_probability_distribution_dict = dict(zip(self.cal_msd_list,
                                 [[] for cal_mol in
                                  self.cal_msd_list]))

    def abstract_centroid(self, file_idx: int):
        file = self.files[file_idx]
        with open(os.path.join(self.chain_path, file), 'rb') as f:
            data = eval(f.read())

        match = re.search(r'\d+', file)
        if match:
            time_step = int(match.group())

        time_step = None
        for mol_name in self.cal_msd_list:
            mol_data = data[mol_name]
            mol_condensate_centroid = []
            for chain_idx in range(len(mol_data)):
                centroid = np.mean(mol_data[chain_idx], axis=0)
                self.centroid_dict[mol_name][chain_idx][file_idx] = centroid.tolist()
                mol_condensate_centroid.append(centroid)

            mol_condensate_centroid = np.mean(mol_condensate_centroid, axis=0)
            mol_condensate_centroid_variance = np.var(mol_condensate_centroid, axis=0)
            self.condensate_centroid_dict[mol_name][file_idx] = {'centroid': mol_condensate_centroid.tolist(),
                                                                 'variance': mol_condensate_centroid_variance.tolist()}

    def abstract_centroid_parallel(self):
        if not os.path.exists(os.path.join(self.msd_path, "centroid.txt")):
            self.centroid_dict: dict = dict(zip(self.cal_msd_list,
                                                [[[[] for _ in range(len(self.files))]
                                                 for _ in range(self.data.mol_class_dict[cal_mol][0])]
                                                for cal_mol in self.cal_msd_list]))
            self.condensate_centroid_dict = dict(zip(self.cal_msd_list,
                                                     [[[] for _ in range(len(self.files))]
                                                     for cal_mol in self.cal_msd_list]))

            list(tqdm(map(self.abstract_centroid, range(len(self.files))),
                      total=len(self.files),
                      desc="提取质心中",
                      colour='cyan',
                      bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]',
                      ncols=100))
            # 存储质心
            with open(os.path.join(self.msd_path, "centroid.txt"), 'wb') as f:
                f.write(str(self.centroid_dict).encode('utf-8'))
            with open(os.path.join(self.msd_path, "condensate_centroid.txt"), 'wb') as f:
                f.write(str(self.condensate_centroid_dict).encode('utf-8'))
        else:
            choice = input("已存在质心文件，是否重新计算？(y/n)")
            if choice == 'y':
                self.centroid_dict: dict = dict(zip(self.cal_msd_list,
                                                    [[[[] for _ in range(len(self.files))]
                                                     for _ in range(self.data.mol_class_dict[cal_mol][0])]
                                                    for cal_mol in self.cal_msd_list]))

                self.condensate_centroid_dict = dict(zip(self.cal_msd_list,
                                                         [[[] for _ in range(len(self.files))]
                                                          for cal_mol in self.cal_msd_list]))

                list(tqdm(map(self.abstract_centroid, range(len(self.files))),
                          total=len(self.files),
                          desc="提取质心中",
                          colour='cyan',
                          bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]',
                          ncols=100))
                # 存储质心
                with open(os.path.join(self.msd_path, "centroid.txt"), 'wb') as f:
                    f.write(str(self.centroid_dict).encode('utf-8'))
                with open(os.path.join(self.msd_path, "condensate_centroid.txt"), 'wb') as f:
                    f.write(str(self.condensate_centroid_dict).encode('utf-8'))
            else:
                with open(os.path.join(self.msd_path, "centroid.txt"), 'rb') as f:
                    self.centroid_dict = eval(f.read())
                with open(os.path.join(self.msd_path, "condensate_centroid.txt"), 'rb') as f:
                    self.condensate_centroid_dict = eval(f.read())


    def cal_msd(self, mol_name):
        traj = np.array(self.centroid_dict[mol_name])  # shape: [n_frames, 3]
        # print(traj.shape)
        mol_num, n_frames, _ = traj.shape
        for mol_idx in range(mol_num):
            msd = []
            for dt in range(1, int(n_frames / 2 + 1)):
                displacements = traj[mol_idx, dt:] - traj[mol_idx, :-dt]
                squared_displacements = np.sum(displacements ** 2, axis=1) / (n_frames - dt)
                msd_dt = np.mean(squared_displacements)
                msd.append(msd_dt)
            self.msd_dict[mol_name][mol_idx] = msd

    def cal_msd_parallel(self):
        # 计算 MSD
        print("✅ 开始计算 MSD")
        list(map(self.cal_msd, self.cal_msd_list))
        # 存储 MSD
        with open(os.path.join(self.save_path, "msd.txt"), 'wb') as f:
            f.write(str(self.msd_dict).encode('utf-8'))

        # 绘图
        fig, ax = plt.subplots(figsize=(12, 9), dpi=300)
        for mol_name in self.cal_msd_list:
            for mol_idx in range(len(self.msd_dict[mol_name])):
                msd = self.msd_dict[mol_name][mol_idx]
                time = np.arange(len(msd)) * 10  # ns
                ax.plot(time, msd, label=f"{mol_name}_{mol_idx}")
            ax.set_xlabel('time (ns)')
            ax.set_ylabel('MSD (nm^2)')
            ax.set_title(f'MSD of {mol_name}')
            # plt.legend()
            fig.savefig(os.path.join(self.save_path, f"draw_msd_{mol_name}.png"))

        print("✅ 计算完成")

    def calDisplacementProbabilityDistribution(self, mol_name):
        traj = np.array(self.centroid_dict[mol_name])  # shape: [n_mols, n_frames, 3]
        print(traj.shape)
        t0 = 0
        dt = 100

        d = traj[:, t0 + dt, :] - traj[:, t0, :]  # 计算位移
        d = np.linalg.norm(d, axis=1)  # 计算位移的模

        d_sqrt = np.sqrt(d)  # 位移的平方根

        # 储存位移
        self.displacement_probability_distribution_dict[mol_name].append(d_sqrt)

    def cal_displacement_probability_distribution_parallel(self):
        # 计算 DPD
        if not os.path.exists(self.displacement_probability_distribution_path):
            os.makedirs(self.displacement_probability_distribution_path, exist_ok=True)
        else:
            shutil.rmtree(self.displacement_probability_distribution_path)
            os.makedirs(self.displacement_probability_distribution_path, exist_ok=True)

        print(f"计算 DPD 时使用的 bin 数量为 {self.displacement_probability_distribution_bin_num} nm")
        print(f"计算 DPD 的边界为 {self.r_max} nm")

        list(tqdm(map(self.calDisplacementProbabilityDistribution, self.cal_msd_list),
                       total=len(self.cal_msd_list),
                       desc="计算 DPD",
                       colour='cyan',
                       bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]',
                       ncols=100))

        # 绘图
        fig, ax = plt.subplots(figsize=(12, 9), dpi=300)
        for mol_name in self.cal_msd_list:
            displacement_probability_distribution = self.displacement_probability_distribution_dict[mol_name]
            ax.hist(displacement_probability_distribution, bins=50, density=True, alpha=0.6, label='MD Data')
            ax.set_xlabel('r (nm)')
            ax.set_ylabel('DPD')
            ax.set_title(f'Displacement Probability Distribution of {mol_name}')
            # plt.legend()
            fig.savefig(os.path.join(self.displacement_probability_distribution_path, f"draw_displacement_probability_distribution_{mol_name}.png"))

        print("✅ 计算完成")