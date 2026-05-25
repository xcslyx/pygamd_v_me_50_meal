import os
import json
import numpy as np

from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from pygamd_v_me_50_meal import utils

# 加载消息文件
with open(os.path.join(os.path.dirname(__file__), 'message.json'), 'r', encoding='utf-8') as f:
    messages = json.load(f)
    # msg = messages['mass_density_message']


class SpanningRatioCalculator:
    def __init__(self, path, data, lang: str="zh"):
        self.path = path
        self.data = data
        self.lang = lang
        
        self.mol_class_dict = self.data.mol_class_dict
        self.length_dict = self.data.length_dict

        xml_path = os.path.join(self.path, "xml")
        self.xml_files = sorted(os.listdir(xml_path))
        init_xml_file = ""
        for i in range(len(self.xml_files)):
            if utils.check_xml_start_tag(self.xml_files[i]) and self.xml_files[i].endswith("0.xml"):
                init_xml_file = os.path.join(self.path, "xml", self.xml_files[i])
                break
        if not init_xml_file:
            raise ValueError("Cannot find a valid initial xml file.")
        self.box_size: list[float] = XMLDataExtractor(init_xml_file).get_box_size()

        self.cutoff = 3.5
        self.tolerance = 0.95


    def calculate_slab_spanning_ratio(self, positions):
        """
        计算 Slab 模拟中凝聚体的横向贯通比 (Spanning Ratio)

        参数:
        positions : numpy.ndarray, shape (N, 3)

        返回:
        spanning_ratio : float 处于贯通网络中的分子数占总输入分子数的比例。
        """
        N = positions.shape[0]
        if N == 0:
            return 0.0


        # 1. 使用 KDTree 快速查找距离小于 cutoff 的分子对
        # 这比计算全排列距离矩阵快得多，且省内存
        tree = cKDTree(positions)
        pairs = tree.query_pairs(r=self.cutoff)


        if not pairs:
            return 0.0  # 没有任何分子相连


        # 2. 构建稀疏邻接矩阵 (Sparse Adjacency Matrix)
        row = [p[0] for p in pairs]
        col = [p[1] for p in pairs]
        data = np.ones(len(pairs))

        # 建立对称的无向图矩阵
        graph = csr_matrix((data, (row, col)), shape=(N, N))
        graph = graph + graph.T


        # 3. 寻找所有连通分量 (划分团簇)
        # n_components: 团簇总数; labels: 每个分子所属的团簇编号
        n_components, labels = connected_components(csgraph=graph, directed=False, return_labels=True)


        # 4. 统计满足横向贯通条件的团簇
        spanning_atoms_count = 0

        for i in range(n_components):
            # 找出属于第 i 个团簇的所有分子索引
            cluster_idx = np.where(labels == i)[0]

            # 优化：如果团簇内分子数极少，直接跳过，不可能贯通
            if len(cluster_idx) < 2:
                continue

            # 获取该团簇的坐标
            cluster_coords = positions[cluster_idx]

            # 计算该团簇在 x 和 y 方向的跨度
            # np.ptp (peak-to-peak) 等价于 max() - min()
            dx = np.ptp(cluster_coords[:, 0])
            dy = np.ptp(cluster_coords[:, 1])

            # 只要在 x 或 y 任意一个方向跨度超过阈值，即判定为贯通簇
            if dx >= self.box_size[0] * self.tolerance or dy >= self.box_size[1] * self.tolerance:
                spanning_atoms_count += len(cluster_idx)


        # 5. 计算贯通比
        spanning_ratio = spanning_atoms_count / N

        return spanning_ratio



    def calculate_slab_spanning_ratio_parallel(self, ):
        pass
# ==========================================
# 使用示例 (Dummy Data Test)
# ==========================================
if __name__ == "__main__":
    # 假设我们有 5000 个分子
    num_molecules = 5000

    # 随机生成一个范围在 x:[0, 100], y:[0, 100], z:[0, 50] 的未缠绕坐标体系
    # (真实情况请从你的轨迹中读取 unwrapped 坐标)
    mock_positions = np.random.rand(num_molecules, 3) * [100, 100, 50]

    Lx, Ly = 100.0, 100.0
    rc = 3.5 # 假设相连的截断距离是 3.5

    ratio = SpanningRatioCalculator().calculate_slab_spanning_ratio(mock_positions)

    print(f"当前帧的贯通比为: {ratio:.4f} ({ratio*100:.2f}%)")
