#!-*- coding: utf-8 -*-
import xml.etree.ElementTree as ET

import numpy as np
import torch as torch

from sklearn.decomposition import PCA


class Functions:
    def __init__(self):
        pass

    @staticmethod
    def replace_position(xml_file, position, new_xml_file: bool=None):
        """
        替换 xml文件中的 position 标签

        :param xml_file: 需要修改的 xml 文件路径
        :param position: 新的位置信息，格式为 [[x1, y1, z1], [x2, y2, z2],...]
        :param new_xml_file: 新 xml 文件路径，如果为 None 则直接修改原文件，否则生成新文件
        :return: None
        """
        tree = ET.parse(xml_file)
        root = tree.getroot()
        for child in root.iter():
            if child.tag == "position":
                child.text = '\n' + '\n'.join([' '.join(map(str, pos)) for pos in position]) + '\n'
        if new_xml_file is None:
            tree.write(xml_file)
        else:
            tree.write(new_xml_file)


    @staticmethod
    def euclidean_distances(a, b):
        """
        计算两个矩阵 A 和 B 的欧氏距离.
        :param a: 矩阵 A，形状为 (n, d)
        :param b: 矩阵 B，形状为 (m, d)
        :return: 矩阵 A 和 B 的距离矩阵，每个元素代表 A 中第 i 个向量和 B 中第 j 个向量之间的距离
        """
        vec_prod = torch.mm(a, b.transpose(0, 1))

        # 计算 A 和 B 的平方和
        a_sq = a ** 2
        sum_a_sq = torch.sum(a_sq, dim=1)
        sum_a_sq_extend = sum_a_sq.view(-1, 1).expand(-1, vec_prod.shape[1])
        b_sq = b ** 2
        sum_b_sq = torch.sum(b_sq, dim=1)
        sum_b_sq_extend = sum_b_sq.view(1, -1).expand(vec_prod.shape[0], -1)

        # 计算距离平方
        ed_sq = sum_b_sq_extend + sum_a_sq_extend - 2 * vec_prod

        # 确保距离平方值为非负
        ed_sq[ed_sq < 0] = 0.0

        # 计算欧氏距离
        ed = torch.sqrt(ed_sq)
        return ed


    @staticmethod
    def cal_sigma_mat(sequence_a: list[str] | str, sequence_b: list[str] | str) -> np.ndarray:
        """
        计算两个序列的平均分子体积矩阵。
        :param sequence_a: 序列 A，可以是字符串或列表
        :param sequence_b: 序列 B，可以是字符串或列表
        :return: 平均分子体积矩阵，形状为 (len(sequence_a), len(sequence_b))
        """
        sigma_dict = {'A':0.504, 'R':0.655999, 'N':0.568, 'D':0.558, 'C':0.548,
                      'Q':0.602, 'E':0.592, 'G':0.45, 'H':0.608, 'I':0.618,
                      'L':0.618, 'K':0.636, 'M':0.618, 'F':0.636, 'P':0.555999,
                      'S':0.518, 'T':0.562, 'W':0.678, 'Y':0.646, 'V':0.586000,
                      "M1": 0.618,
                      "L2": 0.618,
                      # RNA
                      'Ad': 0.844, 'Cd': 0.822, 'Gd': 0.851, 'Ud': 0.817,
                      # DNA
                      "Ph": 0.611, "Su": 0.611, "Ab": 0.611, "Gb": 0.611, "Cb": 0.611, "Tb": 0.611,
                      # Enzo np
                      "CA": 0.47, "Nda": 0.47, "COH": 0.47, "COO": 0.47, }

        # 计算平均分子体积矩阵
        seq_sigma_list_a = np.array([sigma_dict[i] for i in sequence_a])
        seq_sigma_list_b = np.array([sigma_dict[i] for i in sequence_b])
        mean_matrix = (seq_sigma_list_a[:, None] + seq_sigma_list_b[None, :]) / 2
        return mean_matrix

    @staticmethod
    def cal_lambda_mat(sequence_a: list[str] | str, sequence_b: list[str] | str) -> np.ndarray:
        """
        计算两个序列的平均分子体积矩阵。
        :param sequence_a: 序列 A，可以是字符串或列表
        :param sequence_b: 序列 B，可以是字符串或列表
        :return: 平均分子体积矩阵，形状为 (len(sequence_a), len(sequence_b))
        """
        lambda_dict = {
            # Amino acids
            'A': 0.0011162643859539204, 'R': 0.7249915947715212, 'N': 0.43832729970272843, 'D': 0.029182123776349763, 'C': 0.610362354303913,
            'Q': 0.3268188050525212, 'E': 0.006100281608649786, 'G': 0.7012713677972457, 'H': 0.46519480823469783, 'I': 0.6075268330845265,
            'L': 0.5563020305733198, 'K': 0.058617173158697924, 'M': 0.7458993420826714, 'F': 0.9216959832175945, 'P': 0.37296418535993475,
            'S': 0.46485701300656046, 'T': 0.5379777613307019, 'W': 0.9844235478393931, 'Y': 0.9950108229594323, 'V': 0.41850068525598694,

            # DNA coarse-grained beads
            "Ph": 0.459459, "Su": 0.756757, "Ab": 0.351351, "Gb": 0.540541, "Cb": 0.297297, "Tb": 0.594595,
        }

        # 计算平均分子体积矩阵
        seq_lambda_list_a = np.array([lambda_dict[i] for i in sequence_a])
        seq_lambda_list_b = np.array([lambda_dict[i] for i in sequence_b])
        lambda_mean_matrix = (seq_lambda_list_a[:, None] + seq_lambda_list_b[None, :]) / 2
        return lambda_mean_matrix


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


    @staticmethod
    def is_chain_wrapped(positions, box_length: list[float], threshold: float = 0.5):
        """
        判断链是否跨越了盒子的周期边界
        :param positions: 粒子坐标 (N, 3) numpy 数组
        :param box_length: 模拟盒子尺寸
        :param threshold: 跳跃阈值（单位：盒子长度比例），默认0.5（即跳过一半盒子视为跨盒）
        :return: True 表示链是跨盒的，False 表示连续
        """
        positions = np.array(positions)
        box_length = np.array(box_length)
        deltas = positions[1:] - positions[:-1]
        min_image_deltas = deltas - box_length * np.round(deltas / box_length)

        # 如果某个方向上距离跳跃了超过阈值 × 盒子长度，则视为“跨盒”
        jump = np.abs(deltas - min_image_deltas)
        return np.any(jump > (threshold * box_length))

    @staticmethod
    def unwrap_chain(positions: list[list[float]], box_length: list[float]) -> list[list[float]]:
        """
        将一个周期性边界条件下的链解包 (unwrap) 为连续结构
        :param positions: 原始坐标数组，形状为 (N, 3)
        :param box_length: 模拟盒子长度
        :return: 解包后的坐标数组，形状为 (N, 3)
        """
        positions = np.array(positions)
        box_length = np.array(box_length)

        # 存储解包后的坐标
        unwrapped = np.zeros_like(positions)
        unwrapped[0] = positions[0]  # 第一个点保持不变

        for i in range(1, len(positions)):
            delta = positions[i] - positions[i - 1]
            # 最小图像法
            delta -= box_length * np.round(delta / box_length)
            unwrapped[i] = unwrapped[i - 1] + delta

        return unwrapped.tolist()

    @staticmethod
    def rodrigues_rotation(v, z=np.array([0, 0, 1]), axis=None):
        """
        通过罗德里格斯旋转公式计算将一个三维向量 v 旋转到另一个三维向量 z 所需的旋转矩阵。
        :param v: 要旋转的目标向量。
        :param z: 旋转目标方向的向量，默认为 [0, 0, 1]。
        :param axis: 旋转轴，默认为 None，如果提供则直接使用。
        :return: 返回计算得到的旋转矩阵 R。
        """
        # 归一化输入向量
        v = v / np.linalg.norm(v)
        z = z / np.linalg.norm(z)

        # 计算旋转轴（叉积）
        if axis is None:
            axis = np.cross(v, z)
        axis_norm = np.linalg.norm(axis)

        # 如果已经对齐（即v和z平行），直接返回单位矩阵或反向矩阵
        if axis_norm < 1e-8:
            if np.allclose(v, z):
                return np.eye(3)  # 已经对齐
            else:
                return -np.eye(3)  # 相反方向

        axis = axis / axis_norm  # 归一化旋转轴
        angle = np.arccos(np.clip(np.dot(v, z), -1.0, 1.0))  # 计算旋转角度

        # 叉积矩阵 [u]_×
        K = np.array([[0, -axis[2], axis[1]],
                      [axis[2], 0, -axis[0]],
                      [-axis[1], axis[0], 0]])

        # 罗德里格斯旋转公式
        R = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * np.dot(K, K)

        return R

    @staticmethod
    def abstract_centroid(file_name: str, cal_mol: list[str], free_chain_dict: dict[str, list[int]]=None,
                            mass: bool=False, get_height: bool=False,):
        x_mat = eval(open(file_name, 'r').read())

        molecules = []
        for mol_name in cal_mol:
            for chain_idx in range(len(x_mat[mol_name])):
                if free_chain_dict is not None and chain_idx in free_chain_dict[mol_name]:
                    continue
                chain = x_mat[mol_name][chain_idx]
                molecules.append(chain)

        # 展平所有坐标为一个二维数组
        all_coords = np.concatenate(molecules, axis=0)

        mol_condensate_centroid = np.mean(all_coords, axis=0)
        mol_condensate_centroid_variance = np.var(all_coords, axis=0)

        if get_height:
            # 以轴向作为圆柱的轴，找到轴上离中心最远的两个点作为圆柱底面的参考
            pca = PCA(n_components=1)
            pca.fit(all_coords)
            axis_vector = pca.components_[0]

            # 计算把主轴旋转到 z 轴所需的旋转矩阵
            R = Functions.rodrigues_rotation(axis_vector)  # 3 x 3 array

            relative_positions = all_coords - mol_condensate_centroid  # N x 3 array

            # 利用旋转矩阵将所有坐标旋转到 z 轴
            rotated_positions = relative_positions @ R.T  # N x 3 array

            height = np.max(rotated_positions[:, 2]) - np.min(rotated_positions[:, 2])
            res = [mol_condensate_centroid, mol_condensate_centroid_variance, height]
        else:
            res = [mol_condensate_centroid, mol_condensate_centroid_variance]

        return res

    @staticmethod
    def move_chains_towards_com(positions, reference_com, box_size: list[float]):
        """
        移动蛋白质到与参考质心尽量接近的位置。
        :param positions:
        :param reference_com: 参考质心
        :param box_size: 模拟盒子尺寸
        :return: positions: 移动后的蛋白质位置
        """
        for i in range(1, len(positions)):
            for axis in range(3):  # 分别处理x, y, z三个轴
                # 定义在x, y, z轴上的位移量
                shifts = np.array([0, box_size[axis], -box_size[axis]])
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

    @staticmethod
    def adjust_centroid_to_box(re_positions, box_size: list[float]):
        new_re_positions = []
        # print(f"re_positions: {re_positions}")
        positions = [pos for _, pos in re_positions]
        condensate_com = np.mean(np.array(list(map(lambda x: np.mean(x, axis=0), positions))), axis=0)  # 计算凝聚体的质心
        # print(f"condensate_com: {condensate_com}")
        for i in range(3):  # 3维坐标 (x, y, z)
            box_min = -box_size[i] / 2
            box_max = box_size[i] / 2
            cur_move = np.array([0, 0, 0])
            cur_move[i] = box_size[i]
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

    @staticmethod
    def move_chain_into_box(chain_positions, box_size: list[float], max_group_core=np.array([0, 0, 0])):
        # 若链的质心不在盒子内，则使其位于盒子内
        # print(f"chain_positions: {chain_positions}")
        # exit()
        chain_positions -= max_group_core
        com = np.mean(chain_positions, axis=0)
        # print(com)
        # xyz = "xyz"
        for i in range(3):  # 3维坐标 (x, y, z)
            cur_move = np.array([0, 0, 0])
            cur_move[i] = box_size[i]
            # print("cur_move: ", cur_move)
            if com[i] < -box_size[i] / 2:
                # print("cur_move: ", cur_move)
                chain_positions = chain_positions + cur_move
                com = np.mean(chain_positions, axis=0)
                # print(f"new_com: {com}")
            elif com[i] > box_size[i] / 2:
                # print("cur_move: ", cur_move)
                chain_positions = chain_positions - cur_move
                com = np.mean(chain_positions, axis=0)
                # print(f"new_com: {com}")
        return chain_positions

    @staticmethod
    def pre_process(indexed_positions, box_size: list[float], max_group_start=0, condensate_index=None, adjusted_chains=None, max_group_core=None):
        # 计算第一个蛋白质的质心
        if condensate_index is None:
            condensate_index = []
        if adjusted_chains is None:
            adjusted_chains = []
        first_chain_com = np.mean(indexed_positions[0][1], axis=0)

        # 先移动所有链条以形成第一个凝聚体
        new_positions = Functions.move_chains_towards_com([pos for _, pos in indexed_positions], first_chain_com,
                                                          box_size)

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

        if group_start > max_group_start:
            max_group_start = group_start
            max_se_dict = dict(re_positions[:group_start])
            max_updated_positions = [(idx, max_se_dict[idx] if idx in max_se_dict else pos) for idx, pos in
                                     new_indexed_positions]
            max_re_ordered_positions = [position for idx, position in
                                        sorted(max_updated_positions, key=lambda x: x[0])]
            max_new_positions = np.vstack(max_re_ordered_positions)
            max_group_core = np.mean(max_new_positions, axis=0)

        # 递归调用：对第二个小组及第二个小组之后的蛋白质重复操作
        if group_start < len(re_positions):
            adjusted_chains.append(re_positions[:group_start])
            # 对新小组递归调用
            se_positions, _ = Functions.pre_process(re_positions[group_start:], box_size, max_group_start=max_group_start,
                                       condensate_index=condensate_index, adjusted_chains=adjusted_chains, max_group_core=max_group_core)
            # 调整质心到盒子内
            re_positions[:group_start] = Functions.adjust_centroid_to_box(re_positions[:group_start], box_size)
            # 将 se_positions 转换为字典，以 idx 为键，position 为值
            se_dict = dict(se_positions)

            # 替换 new_indexed_positions 中的 position
            updated_positions = [(idx, se_dict[idx] if idx in se_dict else pos) for idx, pos in
                                 new_indexed_positions]
            new_indexed_positions = updated_positions
        return new_indexed_positions, max_group_core
