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
                      # RNA
                      'Ad': 0.844, 'Cd': 0.822, 'Gd': 0.851, 'Ud': 0.817,
                      # DNA
                      "Ph": 0.611, "Su": 0.611, "Ab": 0.611, "Gb": 0.611, "Cb": 0.611, "Tb": 0.611}

        # 计算平均分子体积矩阵
        seq_sigma_list_a = np.array([sigma_dict[i] for i in sequence_a])
        seq_sigma_list_b = np.array([sigma_dict[i] for i in sequence_b])
        mean_matrix = (seq_sigma_list_a[:, None] + seq_sigma_list_b[None, :]) / 2
        return mean_matrix

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
    def abstract_centroid(file_name: str, cal_mol: list[str],
                          free_chain_dict: dict[str, list[int]]=None,
                          mass: bool=False,
                          get_height: bool=False,):
        x_mat = eval(open(file_name, 'r').read())

        # time_step = None
        # match = re.search(r'\d+', file_name)
        # if match:
        #     time_step = int(match.group())

        # centroid, condensate_centroid
        molecules = []
        for mol_name in cal_mol:
            if mol_name == "Na" or mol_name == "Cl":
                continue
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