# from ..pygamd_v_me_50_meal.pygamd_analysis.coordinates_processor import CoordinatesProcessor
#
# from ..pygamd_v_me_50_meal.data import Data
#
# path = "40MED1+200polyA20"
# data = Data(path)
# CoordinatesProcessor(path, data, )

import numpy as np

import torch
#
# d = np.array([[1, 2], [0, 4]])
# print(d != 0)
# d_mask = d != 0
#
# c = np.ones_like(d)
# print(c * d_mask)

import matplotlib.pyplot as plt

epsilon = 0.8368
sig = 0.611
lam = 0.5

d = torch.tensor(np.arange(0, 2, 0.01), device="cuda:0")
print(d)
mask_minimum = d >= 0   # 生成 mask，距离非零为 True

r_safe = torch.maximum(d, torch.tensor(1e-6, device=d.device))
U_LJ = 4 * epsilon * ((sig / r_safe) ** 12 - (sig / r_safe) ** 6)
U_LJ = U_LJ * mask_minimum  # 对零距离的 pair 置 0

r_cut = 2 ** (1 / 6) * sig
mask_rep = (r_safe <= r_cut) & mask_minimum  # 排斥区
mask_att = (r_safe > r_cut) & mask_minimum  # 吸引区

U_AH = (
        (U_LJ + epsilon * (1 - lam)) * mask_rep
        + (lam * U_LJ) * mask_att
)

# 建议的平滑处理逻辑
U_max = lam * epsilon  # 设定一个你能接受的最大能量值
U_AH = torch.clamp(U_AH, max=U_max)

d = d.cpu(); U_AH = U_AH.cpu()
plt.plot(d, U_AH)
plt.savefig("U_AH.png")

# epsilon_r = 74.19
# k = 138.935 / epsilon_r  # 1 / (4 * pi * epsilon_0)
# debye_length = 0.794
# prefactor = k * self.charge_mat
# screening = torch.exp(-r_safe / debye_length)
#
# U_DH = prefactor * screening / r_safe
# U_DH = U_DH * mask_minimum  # 对零距离的 pair 置 0