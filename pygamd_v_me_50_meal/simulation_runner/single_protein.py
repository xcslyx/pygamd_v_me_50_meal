try:
    from poetry import cu_gala as gala
    from poetry import force_field_gala
except ImportError:
    print()
from poetry import _options

from pygamd_v_me_50_meal.simulate_creation.xml_generator import XMLGenerator

class SingleProteinSimulation:
    def __init__(self, protein_name, protein_sequence: str="", protein_pdb_file: str="", _gpu_id: int=0):
        self.protein_name = protein_name

        self.path = f"1{protein_name}"
        if protein_pdb_file:
            filename = XMLGenerator(self.path, protein_pdb_file, gen_run_file=False).get_output_file_name()


        if "filename" not in locals():
            filename = f"{protein_name}.xml"

        print(filename)
        build_method = gala.XMLReader(filename)
        perform_config = gala.PerformConfig(_options.gpu)
        all_info = gala.AllInfo(build_method, perform_config)

    def run_simulation(self):
        pass


