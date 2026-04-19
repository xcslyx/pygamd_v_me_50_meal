import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


class AromaticAnalyzer:
    def __init__(self, window=15, color='#8E44AD'):
        self.window = window
        self.color = color
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
    def sliding_aromaticity(seq, window=15):
        aromatic = {'F': 1, 'W': 1, 'Y': 1}
        N = len(seq)
        values = []
        for i in range(N - window + 1):
            win = seq[i:i+window]
            count = sum(aromatic.get(a, 0) for a in win)
            values.append(count / window)
        return np.array(values)

    @staticmethod
    def sliding_aromaticity_edge_adaptive(seq, window=15, pad_mode='edge'):
        aromatic = {'F': 1, 'W': 1, 'Y': 1}
        pad_width = window // 2
        if isinstance(seq, str) or isinstance(seq, list):
            seq_padded = ["A"] * pad_width + seq + ["A"] * pad_width
        else:
            seq_padded = np.pad(seq, pad_width, mode=pad_mode)

        N = len(seq_padded)
        values = []
        for i in range(N - window + 1):
            win = seq_padded[i:i+window]
            count = sum(aromatic.get(a, 0) for a in win)
            values.append(count / window)
        return np.array(values)

    def _process_segments(self, aromaticity_values):
        x = np.arange(len(aromaticity_values))
        points = np.column_stack([x, aromaticity_values])

        all_sub_segments = []
        all_colors = []

        for i in range(len(points) - 1):
            original_seg = points[i:i + 2]
            all_sub_segments.append(original_seg)
            all_colors.append(self.color)

        return np.array(all_sub_segments), all_colors, x

    def plot_aromaticity(self, seq, save_path=None, title=None, window=None):
        if window is None:
            window = self.window
        aromaticity_values = self.sliding_aromaticity_edge_adaptive(seq, window=window)
        all_sub_segments, all_colors, x = self._process_segments(aromaticity_values)

        fig, ax = plt.subplots(figsize=(12, 8), dpi=300)

        lc = LineCollection(all_sub_segments, colors=all_colors, linewidth=2)
        ax.add_collection(lc)

        ax.fill_between(x, aromaticity_values, 0, interpolate=True,
                        color=self.color, alpha=0.2, zorder=2)

        ax.autoscale()
        ax.set_xlabel("Position")
        ax.set_ylabel("Aromaticity")

        seq_length = len(aromaticity_values)
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
        return fig, ax, aromaticity_values


if __name__ == '__main__':
    analyzer = AromaticAnalyzer()

    pro_name = "hcGAS"
    try:
        hcGAS_seq = list(eval(open(f'{pro_name}_sequence.txt', 'r').read()))
        print("Sequence length:", len(hcGAS_seq))
        hcGAS_seq = hcGAS_seq[:159 if pro_name == "hcGAS" else -1]

        fig, ax, aromaticity_values = analyzer.plot_aromaticity(hcGAS_seq, save_path=f"{pro_name}_aromaticity.png")
        print(f"Aromaticity values shape: {aromaticity_values.shape}")
    except FileNotFoundError:
        # Test with a sample sequence
        sample_seq = "MKVDELVQGLLKQISAEELKKARNEIARQHLEKTHQDLKKDILTYLTDRQIKQLEDAFQKLLAEKTEENKLAQAVENSLGQLEEKLKEA"
        fig, ax, aromaticity_values = analyzer.plot_aromaticity(sample_seq, save_path="sample_aromaticity.png")
        print(f"Sample sequence aromaticity values shape: {aromaticity_values.shape}")