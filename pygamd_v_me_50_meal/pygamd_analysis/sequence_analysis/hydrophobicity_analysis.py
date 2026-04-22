import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from pygamd_v_me_50_meal.constants import HYDROPHOBICITY_LAMBDA


class HydrophobicityAnalyzer:
    def __init__(self, window=15):
        self.window = window
        # 从constants.py导入氨基酸疏水性值字典（lambda值）
        self.lambda_dict = HYDROPHOBICITY_LAMBDA
    
    def calculate_hydrophobicity(self, seq, window=None):
        """
        计算序列的疏水性分布
        :param seq: 氨基酸序列
        :param window: 滑动窗口大小，默认使用初始化时的窗口大小
        :return: 疏水性分布值
        """
        if window is None:
            window = self.window
        
        # 计算每个氨基酸的疏水性值
        hydrophobicity_values = []
        for aa in seq:
            if aa in self.lambda_dict:
                hydrophobicity_values.append(self.lambda_dict[aa])
            else:
                # 对于未知氨基酸，使用平均值
                hydrophobicity_values.append(np.mean(list(self.lambda_dict.values())))  # 这里使用lambda_dict的平均值
        
        # 计算滑动窗口平均
        if len(hydrophobicity_values) < window:
            return np.array(hydrophobicity_values)
        
        hydrophobicity_avg = np.convolve(hydrophobicity_values, np.ones(window)/window, mode='valid')
        return hydrophobicity_avg
    
    def plot_hydrophobicity(self, seq, save_path=None, title=None, window=None):
        """
        绘制序列的疏水性分布
        :param seq: 氨基酸序列
        :param save_path: 保存路径，默认为None
        :param title: 图表标题，默认为None
        :param window: 滑动窗口大小，默认使用初始化时的窗口大小
        :return: 图表对象
        """
        if window is None:
            window = self.window
        
        hydrophobicity_avg = self.calculate_hydrophobicity(seq, window=window)
        
        # 创建图表
        fig, ax = plt.subplots(figsize=(12, 6))
        
        # 绘制疏水性分布
        x = np.arange(len(hydrophobicity_avg)) + window // 2
        ax.plot(x, hydrophobicity_avg, 'b-', linewidth=2)
        
        # 添加颜色条
        points = np.array([x, hydrophobicity_avg]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        norm = plt.Normalize(0, 1)
        lc = LineCollection(segments, cmap='coolwarm', norm=norm)
        lc.set_array(hydrophobicity_avg)
        lc.set_linewidth(2)
        line = ax.add_collection(lc)
        cbar = fig.colorbar(line, ax=ax)
        cbar.set_label('Hydrophobicity (λ)', fontsize=12)
        
        # 设置图表属性
        ax.set_xlabel('Residue Position', fontsize=14)
        ax.set_ylabel('Hydrophobicity', fontsize=14)
        ax.set_ylim(0, 1)
        ax.set_title(title or f'Hydrophobicity Analysis (window size: {window})', fontsize=16)
        ax.grid(True, alpha=0.3)
        
        # 保存图表
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig, ax, hydrophobicity_avg
