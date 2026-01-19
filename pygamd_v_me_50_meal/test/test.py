# from ..pygamd_v_me_50_meal.pygamd_analysis.coordinates_processor import CoordinatesProcessor
#
# from ..pygamd_v_me_50_meal.data import Data
#
# path = "40MED1+200polyA20"
# data = Data(path)
# CoordinatesProcessor(path, data, )

import numpy as np

d = np.array([[1, 2], [0, 4]])
print(d != 0)
d_mask = d != 0

c = np.ones_like(d)
print(c * d_mask)

