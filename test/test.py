# from pygamd_v_me_50_meal.simulation_runner.single_protein import SingleProteinSimulation
from pygamd_v_me_50_meal.pygamd_analysis.contact_map_calculator import ContactMapCalculator
from pygamd_v_me_50_meal.pygamd_analysis.coordinates_processor import CoordinatesProcessor

from pygamd_v_me_50_meal.data import Data

# SingleProteinSimulation("hcGAS", protein_pdb_file="1hcGAS/cGAS_human_init.pdb").run_simulation()
data = Data("20fus-526")
CoordinatesProcessor("20fus-526", data)
ContactMapCalculator("20fus-526", data, cm_choice='/', r_cut=4)

