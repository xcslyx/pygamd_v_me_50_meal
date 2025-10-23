import os
import re
import shutil
import xml.etree.ElementTree as ET

import numpy as np

from tqdm import tqdm
from multiprocessing import Pool

from .get_sequence import GetSequence


# 处理坐标文件
class CoordinatesProcessor:
    def __init__(self, data, file_args, ):
        self.data = data
        self.file_args = file_args
        self.path = file_args.path
        self.mol_class_dict = data.mol_class_dict
        self.length_dict = data.length_dict

        self.init_xml_path = os.path.join(self.path, "xml/")
        self.unwrapping_xml_path = os.path.join(self.path, "xml_unwrapping/")
        self.remove_pbc_condensate_xml_path = os.path.join(self.path, "xml_remove_pbc_condensate/")

        self.box_size = None

        os.makedirs(self.init_xml_path, exist_ok=True)

        self.unwrap_xml_flag: bool = True

        self.max_group_start = 0
        self.max_group_core = None

        for file in os.listdir(self.path):
            if file.startswith("particles") and file.endswith(".xml"):
                file_path = os.path.join(self.path, file)
                destination_path = os.path.join(self.init_xml_path, file)

                # 在移动之前检查目标文件是否存在，如果存在则删除
                if os.path.exists(destination_path):
                    os.remove(destination_path)

                shutil.move(file_path, self.init_xml_path)


    def remove_enm_bonds_from_xml(self, xml_file):
        if xml_file.startswith("particles") and xml_file.endswith("0.xml"):
            pass
        else:
            return

        tree = ET.parse(os.path.join(self.init_xml_path, xml_file))
        root = tree.getroot()

        for bonds in root.findall('.//bond'):
            modified_bonds = []
            for line in bonds.text.strip().split('\n'):
                pattern = r"C\d+-C\d+"
                match = re.search(pattern, line)
                # 检查是否匹配
                if not match:
                    modified_bonds.append(line)
            # 将修改后的文本重新写回到 bond 元素中
            bonds.text = '\n' + '\n'.join(modified_bonds) + '\n'

        # 将修改后的 XML 对象写回到文件中
        tree.write(os.path.join(self.init_xml_path, xml_file), encoding='utf-8', xml_declaration=True)

    @staticmethod
    def is_chain_wrapped(positions, box_length, threshold=0.5):
        """
        判断链是否跨越了盒子的周期边界
        :param positions: 粒子坐标 (N, 3) numpy 数组
        :param box_length: 模拟盒子尺寸
        :param threshold: 跳跃阈值（单位：盒子长度比例），默认0.5（即跳过一半盒子视为跨盒）
        :return: True 表示链是跨盒的，False 表示连续
        """
        positions = np.array(positions)
        box_length = np.array([box_length, box_length, box_length])
        deltas = positions[1:] - positions[:-1]
        min_image_deltas = deltas - box_length * np.round(deltas / box_length)

        # 如果某个方向上距离跳跃了超过阈值 × 盒子长度，则视为“跨盒”
        jump = np.abs(deltas - min_image_deltas)
        return np.any(jump > (threshold * box_length))

    @staticmethod
    def unwrap_chain(positions: list[list[float]], box_length: float) -> list[list[float]]:
        """
        将一个周期性边界条件下的链解包 (unwrap) 为连续结构
        :param positions: 原始坐标数组，形状为 (N, 3)
        :param box_length: 模拟盒子长度
        :return: 解包后的坐标数组，形状为 (N, 3)
        """
        positions = np.array(positions)
        box_length = np.array([box_length, box_length, box_length])

        # 存储解包后的坐标
        unwrapped = np.zeros_like(positions)
        unwrapped[0] = positions[0]  # 第一个点保持不变

        for i in range(1, len(positions)):
            delta = positions[i] - positions[i - 1]
            # 最小图像法
            delta -= box_length * np.round(delta / box_length)
            unwrapped[i] = unwrapped[i - 1] + delta

        return unwrapped.tolist()

    def abstract_coordinates_1(self, xml_file):
        if xml_file.startswith("particles") and xml_file.endswith("0.xml"):
            pass
        else:
            return
        positions, unwrapped_positions = {}, {}

        tree = ET.parse(os.path.join(self.init_xml_path, xml_file))
        root = tree.getroot()

        if not self.box_size:
            self.box_size = float(root.find('.//box').attrib['lx'])

        # 读取 position 数据
        position_elem = root.find('.//position')
        position_elem_text = position_elem.text
        while "\n\n" in position_elem_text:
            position_elem_text = position_elem_text.replace("\n\n", "\n")
        position, unwrapped_position = [], []
        for line in position_elem_text.strip().split('\n'):
            position.append(list(map(float, line.strip().split())))

        for type_ in self.mol_class_dict:
            mol_positions, unwrapped_mol_positions = [], []
            count = self.mol_class_dict[type_][0]
            length = self.mol_class_dict[type_][1]


            for _ in range(count):
                cur_position = position[:length]
                mol_positions.append(cur_position)
                position = position[length:]

                if self.is_chain_wrapped(cur_position, self.box_size):
                    unwrapped_cur_position = self.unwrap_chain(cur_position, self.box_size)
                else:
                    unwrapped_cur_position = cur_position

                unwrapped_mol_positions.append(unwrapped_cur_position)
                unwrapped_position.extend(unwrapped_cur_position)

            positions[type_] = mol_positions
            unwrapped_positions[type_] = unwrapped_mol_positions

        with open(os.path.join(self.path, "chain_xyz/",
                               xml_file), 'w') as f_chain_xyz:
            f_chain_xyz.writelines(str(positions))

        with open(os.path.join(self.path, "chain_xyz_unwrapping/",
                               xml_file.replace("0.xml", "0.reimage.xml")), 'w') as f_chain_xyz:
            f_chain_xyz.writelines(str(unwrapped_positions))
            # 重新写入 position 数据
            position_elem.text = '\n' + '\n'.join([' '.join(map(str, pos)) for pos in unwrapped_position]) + '\n'

            # 删除 angle 和 dihedral 数据
            angle_elem = root.find('.//angle')
            if angle_elem is not None and angle_elem in root:
                root.remove(angle_elem)
            dihedral_elem = root.find('.//dihedral')
            if dihedral_elem is not None and dihedral_elem in root:
                root.remove(dihedral_elem)
            # 删除 charge 和 mass 数据
            charge_elem = root.find('.//charge')
            if charge_elem is not None and charge_elem in root:
                root.remove(charge_elem)
            mass_elem = root.find('.//mass')
            if mass_elem is not None and mass_elem in root:
                root.remove(mass_elem)
            # 删除 body 数据
            body_elem = root.find('.//body')
            if body_elem is not None and body_elem in root:
                root.remove(body_elem)
            # 删除 image 数据
            image_elem = root.find('.//image')
            if image_elem is not None and image_elem in root:
                root.remove(image_elem)

            # 写入修改后的 XML 对象
            tree.write(os.path.join(self.unwrapping_xml_path, xml_file.replace("0.xml", "0.reimage.xml")),
                       encoding='utf-8', xml_declaration=True)

    def abstract_coordinates_2(self, xml_file):
        if os.path.exists(self.remove_pbc_condensate_xml_path):
            tree = ET.parse(os.path.join(self.remove_pbc_condensate_xml_path, xml_file.replace(".xml", ".reimage.new.xml")))
            root = tree.getroot()

            # 读取 position 数据
            position_elem = root.find('.//position')
            position = []
            for line in position_elem.text.strip().split('\n'):
                position.append(list(map(float, line.strip().split())))

            positions = {}
            for type_ in self.mol_class_dict:
                mol_positions = []
                count = self.mol_class_dict[type_][0]
                length = self.mol_class_dict[type_][1]

                for _ in range(count):
                    cur_position = position[:length]
                    mol_positions.append(cur_position)
                    position = position[length:]


                positions[type_] = mol_positions

            with open(os.path.join(self.path, "chain_xyz_remove_pbc_condensate/", xml_file), 'w') as f_chain_xyz:
                f_chain_xyz.writelines(str(positions))

    def cal_xyz(self, remove_condensate_pbc=True):
        init_files = sorted([i for i in os.listdir(self.init_xml_path) if i.endswith("0.xml") and i.startswith("particles")])
        if self.file_args.remove_enm:
            remove_enm_bonds_request = "y"
        elif self.file_args.remove_enm is None:
            remove_enm_bonds_request = input("是否需要移除弹性键？(y/n)")
        else:
            remove_enm_bonds_request = "n"

        if remove_enm_bonds_request.lower() in ["y", "yes"]:
            # 使用 tqdm 和多进程移除弹性键
            with Pool(processes=4) as pool:
                list(tqdm(pool.imap(self.remove_enm_bonds_from_xml, init_files),
                          total=len(init_files),
                          desc="移除弹性键进度", colour="cyan", ncols=100))

        for _dir in ["chain_xyz", "chain_xyz_unwrapping"]:
            if _dir not in os.listdir(self.path):
                print(f"✅ 创建{_dir}文件夹...")
                os.mkdir(os.path.join(self.path, _dir))
            else:
                shutil.rmtree(os.path.join(self.path, _dir))
                os.mkdir(os.path.join(self.path, _dir))
        # print(init_files + unwrapping_files)

        if os.path.exists(self.unwrapping_xml_path):
            if len(os.listdir(self.init_xml_path)) == len(os.listdir(self.unwrapping_xml_path)):
                unwrap_request = input("已存在去周期后的 XML 文件，是否需要重新去周期？(y/n)")
                if unwrap_request.lower() == "n":
                    self.unwrap_xml_flag = False

        if self.unwrap_xml_flag:
            # 创建 xml_unwrapping 文件夹
            os.makedirs(self.unwrapping_xml_path, exist_ok=True)

        with Pool(processes=4) as pool:
            print("开始处理文件...")
            # 使用 tqdm 包装你的可迭代对象
            list(tqdm(pool.imap(self.abstract_coordinates_1, init_files),
                      total=len(init_files),
                      desc="处理进度",
                      colour="cyan",
                      ncols=100))

        print("所有文件处理完成。")
        print("开始提取序列信息...")
        self.file_args.output = os.path.join(self.path, f"{self.data.system_name}_sequence.txt")
        GetSequence(self.init_xml_path, sorted(init_files)[0], self.data).xml2sequence()
        print(f"序列信息已保存至 {self.file_args.output}。")

        if self.file_args.remove_condensate_pbc:
            print("开始进行针对凝聚体的 PBC 去除")
            for _dir in ["chain_xyz_remove_pbc_condensate"]:
                if _dir not in os.listdir(self.path):
                    print(f"✅ 创建{_dir}文件夹...")
                    os.mkdir(os.path.join(self.path, _dir))
                else:
                    shutil.rmtree(os.path.join(self.path, _dir))
                    os.mkdir(os.path.join(self.path, _dir))
            self.remove_pbc_condensate_parallel()

            with Pool(processes=4) as pool:
                print("开始处理去周期后的凝聚体文件...")
                # 使用 tqdm 包装你的可迭代对象
                list(tqdm(pool.imap(self.abstract_coordinates_2, init_files),
                          total=len(init_files),
                          desc="处理进度",
                          colour="cyan",
                          ncols=100))

    @staticmethod
    def move_chains_towards_com(positions, reference_com, box_size):
        """
        移动蛋白质到与参考质心尽量接近的位置。
        :param positions:
        :param reference_com:
        :param box_size:
        :return: positions: 移动后的蛋白质位置
        """
        # 定义在x, y, z轴上的位移量
        shifts = np.array([0, box_size, -box_size])

        for i in range(1, len(positions)):
            # 对x轴进行调整
            for axis in range(3):  # 分别处理x, y, z三个轴
                best_position = positions[i].copy()
                best_distance = abs(positions[i][:, axis] - reference_com[axis])

                # 依次尝试三个位移，看看是否可以减少距离
                for shift in shifts:
                    shifted_position = positions[i].copy()
                    shifted_position[:, axis] += shift  # 仅改变当前轴的坐标
                    shifted_distance = abs(shifted_position[:, axis] - reference_com[axis])

                    # 如果位移后的距离更小，则保留该位移
                    if np.all(shifted_distance < best_distance):
                        best_position = shifted_position
                        best_distance = shifted_distance

                # 更新粒子的位置
                positions[i] = best_position
        return positions


    def remove_pbc_condensate(self, xml_file):
        """
        去除凝聚体的 PBC 条件，使得凝聚体的坐标不再受到 PBC 条件的限制。
        :param xml_file: 凝聚体的 XML 文件
        :return: None
        """
        condensate_index = []
        adjusted_chains = []

        def move_chain_into_box(chain_positions):
            # 若链的质心不在盒子内，则使其位于盒子内
            # print(f"chain_positions: {chain_positions}")
            # exit()
            chain_positions -= self.max_group_core
            com = np.mean(chain_positions, axis=0)
            # print(com)
            # xyz = "xyz"
            for i in range(3):  # 3维坐标 (x, y, z)
                cur_move = np.array([0, 0, 0])
                cur_move[i] = self.box_size
                # print("cur_move: ", cur_move)
                if com[i] < -self.box_size / 2:
                    # print("cur_move: ", cur_move)
                    chain_positions = chain_positions + cur_move
                    com = np.mean(chain_positions, axis=0)
                    # print(f"new_com: {com}")
                elif com[i] > self.box_size / 2:
                    # print("cur_move: ", cur_move)
                    chain_positions = chain_positions - cur_move
                    com = np.mean(chain_positions, axis=0)
                    # print(f"new_com: {com}")
            return chain_positions


        def adjust_centroid_to_box(re_positions):
            new_re_positions = []
            # print(f"re_positions: {re_positions}")
            box_min = -self.box_size / 2
            box_max = self.box_size / 2
            positions = [pos for _, pos in re_positions]
            condensate_com = np.mean(np.array(list(map(lambda x: np.mean(x, axis=0), positions))), axis=0)  # 计算凝聚体的质心
            # print(f"condensate_com: {condensate_com}")
            for i in range(3):  # 3维坐标 (x, y, z)
                cur_move = np.array([0, 0, 0])
                cur_move[i] = self.box_size
                # print(cur_move)
                while condensate_com[i] < box_min:
                    positions = list(map(lambda x: x + cur_move, positions))
                    condensate_com = np.mean(np.array(list(map(lambda x: np.mean(x, axis=0), positions))), axis=0)  # 更新质心
                    # print(f"condensate_com: {condensate_com}")
                while condensate_com[i] > box_max:
                    positions = list(map(lambda x: x - cur_move, positions))
                    condensate_com = np.mean(np.array(list(map(lambda x: np.mean(x, axis=0), positions))), axis=0)  # 更新质心
                    # print(f"condensate_com: {condensate_com}")

            for i in range(len(re_positions)):
                new_re_positions.append((re_positions[i][0], positions[i]))
            return new_re_positions

        def pre_process(indexed_positions):
            # 计算第一个蛋白质的质心
            first_chain_com = np.mean(indexed_positions[0][1], axis=0)

            # 先移动所有链条以形成第一个凝聚体
            new_positions = self.move_chains_towards_com([pos for _, pos in indexed_positions], first_chain_com,
                                                         self.box_size)

            # 重新组合 indexed_positions
            new_indexed_positions = [(indexed_positions[i][0], new_positions[i]) for i in range(len(new_positions))]

            # 计算所有蛋白质到第一个凝聚体的距离
            distances = []
            for i, (index, pos) in enumerate(new_indexed_positions):
                current_com = np.mean(pos, axis=0)
                distance = np.linalg.norm(current_com - first_chain_com)
                distances.append((distance, index))  # 记录距离和原始索引

            # 按照距离第一个凝聚体的远近重新排序
            distances.sort(key=lambda x: x[0])

            position_dict = dict(new_indexed_positions)

            # 从 distances 中获取 re_positions
            re_positions = [(idx, position_dict[idx]) for idx in [pos[1] for pos in distances] if idx in position_dict]

            first_com = np.mean(re_positions[0][1], axis=0)  # 初始位置的质心

            group_start = None

            last_distance = 0.
            for i in range(len(re_positions)):  # 从第1个索引开始遍历
                current_com = np.mean(re_positions[i][1], axis=0)
                cur_distance = np.linalg.norm(current_com - first_com)
                distance = cur_distance - last_distance  # 与上一个质心的距离差值
                last_distance = distance
                # print(distance)
                if distance > 4:
                    group_start = i  # 如果距离差值大于 4，则认为是第二个凝聚体
                    condensate_index.append(group_start + (condensate_index[-1] if condensate_index else 0))  # 记录凝聚体的索引
                    break
                if group_start is None:
                    group_start = len(re_positions)  # 如果没有大于 4 的，所有都属于第一个凝聚体

            if group_start > self.max_group_start:
                self.max_group_start = group_start
                max_se_dict = dict(re_positions[:group_start])
                max_updated_positions = [(idx, max_se_dict[idx] if idx in max_se_dict else pos) for idx, pos in
                                         new_indexed_positions]
                max_re_ordered_positions = [position for idx, position in
                                            sorted(max_updated_positions, key=lambda x: x[0])]
                max_new_positions = np.vstack(max_re_ordered_positions)
                self.max_group_core = np.mean(max_new_positions, axis=0)

            # 递归调用：对第二个小组及第二个小组之后的蛋白质重复操作
            if group_start < len(re_positions):
                adjusted_chains.append(re_positions[:group_start])
                # 对新小组递归调用
                se_positions = pre_process(re_positions[group_start:])
                # 调整质心到盒子内
                re_positions[:group_start] = adjust_centroid_to_box(re_positions[:group_start])
                # print(f"re_positions: {re_positions[:group_start][0]}")
                # print(f"adjust_centroid_to_box(re_positions[:group_start]): {adjust_centroid_to_box(re_positions[:group_start])[0]}")
                # 将 se_positions 转换为字典，以 idx 为键，position 为值
                se_dict = dict(se_positions)

                # 替换 new_indexed_positions 中的 position
                updated_positions = [(idx, se_dict[idx] if idx in se_dict else pos) for idx, pos in
                                     new_indexed_positions]
                new_indexed_positions = updated_positions
            return new_indexed_positions

        if not (xml_file.startswith("particles") and xml_file.endswith("0.reimage.xml")):
            return

        tree = ET.parse(os.path.join(self.unwrapping_xml_path, xml_file))
        root = tree.getroot()
        self.box_size = float(root.find('.//box').attrib['lx'])

        # 读取 position 数据
        position = []
        position_elem = root.find('.//position')
        for line in position_elem.text.strip().split('\n'):
            position.append(list(map(float, line.strip().split())))

        positions = []
        for type_ in self.mol_class_dict:
            count = self.mol_class_dict[type_][0]
            length = self.mol_class_dict[type_][1]
            for _ in range(count):
                positions.append(np.array(position[:length]))
                position = position[length:]
        del position
        # print(len(positions))

        indexed_positions = [(i, pos) for i, pos in enumerate(positions)]

        not_ordered_positions = pre_process(indexed_positions)

        for i in range(len(not_ordered_positions)):
            # print(f"moving: {not_ordered_positions[i][0]}")
            not_ordered_positions[i] = (not_ordered_positions[i][0], move_chain_into_box(not_ordered_positions[i][1]))
        # 按照 idx 排序并提取 position
        re_ordered_positions = [position for idx, position in sorted(not_ordered_positions, key=lambda x: x[0])]

        # 更新后的蛋白质位置
        new_positions = np.vstack(re_ordered_positions)

        new_xml_file = xml_file.replace(".xml", ".new.xml")
        tree = ET.parse(os.path.join(self.unwrapping_xml_path, xml_file))
        root = tree.getroot()

        # 找到 position 标签
        position_element = root.find('.//position')
        pos_text = "\n".join(["    ".join(map(str, p)) for p in new_positions])
        position_element.text = f"\n{pos_text}\n"
        # 保存更新后的XML文件
        tree.write(os.path.join(self.remove_pbc_condensate_xml_path, new_xml_file), encoding='utf-8', xml_declaration=True)


    def remove_pbc_condensate_parallel(self):
        """
        多进程去除凝聚体 PBC 条件。
        """
        if not os.path.exists(self.remove_pbc_condensate_xml_path):
            os.makedirs(self.remove_pbc_condensate_xml_path, exist_ok=True)
            if not os.path.exists(self.unwrapping_xml_path):
                print("未找到去周期后的 XML 文件，请先运行去周期脚本。")
                return

        xml_files = sorted(os.listdir(self.unwrapping_xml_path))
        print(f"开始处理 {len(xml_files)} 个文件...")
        with Pool(processes=4) as pool:
            list(tqdm(pool.imap(self.remove_pbc_condensate, xml_files),
                      total=len(xml_files),
                      desc="去除 PBC 条件",
                      colour='cyan',
                      ncols=100))
