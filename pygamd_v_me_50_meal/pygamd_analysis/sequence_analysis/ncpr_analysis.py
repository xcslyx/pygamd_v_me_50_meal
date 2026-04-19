import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.collections import LineCollection


class NCPRAnalyzer:
    def __init__(self, window=15, color_pos='#E64B35', color_neg='#4DBBD5'):
        self.window = window
        self.color_pos = color_pos
        self.color_neg = color_neg
        self._configure_matplotlib()

    def _configure_matplotlib(self):
        plt.rcParams["font.family"] = "Arial"
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

    @staticmethod
    def sliding_ncpr(seq, window=15):
        charges = {'K': 1, 'R': 1, 'D': -1, 'E': -1}
        N = len(seq)
        values = []
        for i in range(N - window + 1):
            win = seq[i:i+window]
            q = sum(charges.get(a, 0 if isinstance(a, str) else a) for a in win)
            values.append(q / window)
        return np.array(values)

    @staticmethod
    def sliding_ncpr_edge_adaptive(seq, window=15, pad_mode='edge'):
        charges = {'K': 1, 'R': 1, 'D': -1, 'E': -1}
        pad_width = window // 2
        if isinstance(seq, str) or isinstance(seq, list):
            seq_padded = ["A"] * pad_width + seq + ["A"] * pad_width
        else:
            seq_padded = np.pad(seq, pad_width, mode=pad_mode)

        N = len(seq_padded)
        values = []
        for i in range(N - window + 1):
            win = seq_padded[i:i+window]
            q = sum(charges.get(a, 0 if isinstance(a, str) else a) for a in win)
            values.append(q / window)
        return np.array(values)

    def _split_segment_by_zero(self, seg):
        (x1, y1), (x2, y2) = seg
        if y1 * y2 >= 0:
            if y1 > 0 or y2 > 0:
                color = self.color_pos
            elif y1 < 0 or y2 < 0:
                color = self.color_neg
            else:
                color = 'black'
            return [seg], [color]
        elif y1 * y2 < 0:
            x_interp = x1 + (0 - y1) * (x2 - x1) / (y2 - y1)
            interp_point = (x_interp, 0)

            seg_above = []
            seg_below = []
            if y1 > 0:
                seg_above = [[x1, y1], interp_point]
                seg_below = [interp_point, [x2, y2]]
            else:
                seg_above = [interp_point, [x2, y2]]
                seg_below = [[x1, y1], interp_point]

            sub_segs = []
            sub_colors = []
            if seg_above:
                sub_segs.append(seg_above)
                sub_colors.append(self.color_pos)
            if seg_below:
                sub_segs.append(seg_below)
                sub_colors.append(self.color_neg)
            return sub_segs, sub_colors
        else:
            color = self.color_pos if y1 > 0 or y2 > 0 else self.color_neg
            return [seg], [color]

    def _process_segments(self, ncpr_values):
        x = np.arange(len(ncpr_values))
        points = np.column_stack([x, ncpr_values])

        all_sub_segments = []
        all_colors = []

        for i in range(len(points) - 1):
            original_seg = points[i:i + 2]
            sub_segs, sub_cols = self._split_segment_by_zero(original_seg)
            all_sub_segments.extend(sub_segs)
            all_colors.extend(sub_cols)

        return np.array(all_sub_segments), all_colors, x

    def plot_ncpr(self, seq, save_path=None, title=None, window=None):
        if window is None:
            window = self.window
        ncpr_values = self.sliding_ncpr_edge_adaptive(seq, window=window)
        all_sub_segments, all_colors, x = self._process_segments(ncpr_values)

        fig, ax = plt.subplots(figsize=(12, 8), dpi=300)

        lc = LineCollection(all_sub_segments, colors=all_colors, linewidth=2)
        ax.add_collection(lc)

        ax.fill_between(x, ncpr_values, 0, where=(ncpr_values > 0), interpolate=True,
                        color=self.color_pos, alpha=0.2, zorder=2)
        ax.fill_between(x, ncpr_values, 0, where=(ncpr_values < 0), interpolate=True,
                        color=self.color_neg, alpha=0.2, zorder=2)

        ax.autoscale()
        ax.axhline(0, color='black', linestyle='--')
        ax.set_xlabel("Position")
        ax.set_ylabel("NCPR")

        seq_length = len(ncpr_values)
        if seq_length <= 100:
            tick_interval = 10
        elif seq_length <= 200:
            tick_interval = 20
        elif seq_length <= 500:
            tick_interval = 50
        elif seq_length <= 1000:
            tick_interval = 100
        else:
            tick_interval = 200

        xticks = np.arange(0, seq_length, tick_interval)
        if seq_length - xticks[-1] > 10:
            xticks = np.append(xticks, seq_length)
        ax.set_xticks(xticks)
        ax.set_xticklabels([0] + list(xticks)[1:])

        if title:
            ax.set_title(title)

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300)
        return fig, ax, ncpr_values


if __name__ == '__main__':
    analyzer = NCPRAnalyzer()

    pro_name = "hcGAS"
    hcGAS_seq = list(eval(open(f'{pro_name}_sequence.txt', 'r').read()))
    print("Sequence length:", len(hcGAS_seq))
    hcGAS_seq = hcGAS_seq[:159 if pro_name == "hcGAS" else -1]

    fig, ax, ncpr_values = analyzer.plot_ncpr(hcGAS_seq, save_path=f"{pro_name}_NCPR.png")
    print(f"NCPR values shape: {ncpr_values.shape}")

    counts_df = pd.read_csv("../cGAS跨物种/cGAS/cGAS_IDR 物种保守性.csv", index_col=0)
    fig2, ax2, conserved_ncpr = analyzer.plot_conserved_ncpr(counts_df, save_path="hcGAS_conserved_NCPR.png")
    print(f"Conserved NCPR values shape: {conserved_ncpr.shape}")