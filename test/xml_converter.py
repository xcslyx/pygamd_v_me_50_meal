from pygamd_v_me_50_meal.simulate_creation import xml_converter
from pygamd_v_me_50_meal.data import Data

path = "2DNA-134"
data = Data(path)

xml_converter.XMLConverter(data, path, "45bpDNA.xml").xml2pdb()
