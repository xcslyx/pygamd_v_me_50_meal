import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from pygamd_v_me_50_meal.pygamd_analysis.sequence_analysis.ncpr_analysis import NCPRAnalyzer
from pygamd_v_me_50_meal.pygamd_analysis.sequence_analysis.aromatic_analysis import AromaticAnalyzer
from pygamd_v_me_50_meal.pygamd_analysis.sequence_analysis.hydrophobicity_analysis import HydrophobicityAnalyzer


class SequenceAnalyzer:
    def __init__(self, window=15):
        self.window = window
        self.ncpr_analyzer = NCPRAnalyzer(window=window)
        self.aromatic_analyzer = AromaticAnalyzer(window=window)
        self.hydrophobicity_analyzer = HydrophobicityAnalyzer(window=window)

    def plot_ncpr(self, seq, save_path=None, title=None, window=None):
        return self.ncpr_analyzer.plot_ncpr(seq, save_path=save_path, title=title, window=window)

    def plot_aromaticity(self, seq, save_path=None, title=None, window=None):
        return self.aromatic_analyzer.plot_aromaticity(seq, save_path=save_path, title=title, window=window)
        
    def plot_hydrophobicity(self, seq, save_path=None, title=None, window=None):
        return self.hydrophobicity_analyzer.plot_hydrophobicity(seq, save_path=save_path, title=title, window=window)