from pygamd_v_me_50_meal.simulate_creation import xml_converter
from pygamd_v_me_50_meal.data import Data

path = "1hcGAS"
data = Data(path)

xml_converter.XMLConverter(data, path, "cGAS_human_init.xml").xml2pdb()
