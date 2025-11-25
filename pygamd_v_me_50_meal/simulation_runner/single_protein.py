import os

try:
    from poetry import cu_gala as gala
    from poetry import force_field_gala
except ImportError:
    raise ImportError("Please run this script on Linux system with PYGAMD installed.")
from poetry import _options

from scipy.constants import R as gas_constant

from pygamd_v_me_50_meal.simulate_creation.xml_generator import XMLGenerator

class SingleProteinSimulation:
    def __init__(self, protein_name, protein_sequence: str="", protein_pdb_file: str="", _gpu_id: int=0,):
        self.protein_name = protein_name

        self.path = f"1{protein_name}"

        self.add_rigid_body = False
        if protein_pdb_file:
            xml_file = XMLGenerator(self.path, protein_pdb_file, gen_run_file=False)
            filename = xml_file.get_output_file_name()
            self.add_rigid_body = xml_file.add_rigid_body

        self.gen_force_field_file(self.path)

        if "filename" not in locals():
            filename = f"{protein_name}.xml"

        print(filename)
        self.build_method = gala.XMLReader(filename)
        self.perform_config = gala.PerformConfig(_options.gpu)


    def run_simulation(self, dt: float=0.02, temperature: float=310, r_cut = 4.0):
        os.chdir(self.path)
        all_info = gala.AllInfo(self.build_method, self.perform_config)
        app = gala.Application(all_info, dt)

        neighbor_list = gala.NeighborList(all_info, r_cut, r_cut / 10)  # (,rcut,rbuffer)
        neighbor_list.exclusion(["bond"])
        neighbor_list.addExclusionsFromBonds()
        neighbor_list.addExclusionsFromAngles()
        neighbor_list.addExclusionsFromDihedrals()

        ShortRangeEpsilon = 0.8368
        debye_length = 0.794
        ahdh = force_field_gala.AHDHForce(all_info, neighbor_list, r_cut,
                                          ShortRangeEpsilon, debye_length, "ahdh.force_field")
        app.add(ahdh)

        bondforce_pro = gala.BondForceHarmonic(all_info)
        bondforce_pro.setParams('B-B', 8033.28, 0.38)
        app.add(bondforce_pro)

        temperature = temperature * gas_constant / 1000.0  # reduced unit

        groupall = gala.ParticleSet(all_info, "all")
        comp_info = gala.ComputeInfo(all_info, groupall)

        if self.add_rigid_body:
            body_group = gala.ParticleSet(all_info, 'body')
            rigid_body = gala.LangevinNVTRigid(all_info, body_group, temperature, 12345)
            app.add(rigid_body)

            non_body_group = gala.ParticleSet(all_info, 'non_body')
            non_body = gala.LangevinNVT(all_info, non_body_group, temperature, 12345)
            app.add(non_body)
        else:
            gp = gala.LangevinNVT(all_info, groupall, temperature, 12345)
            app.add(gp)

        sort_method = gala.Sort(all_info)
        sort_method.setPeriod(500)
        app.add(sort_method)

        ZeroMomentum = gala.ZeroMomentum(all_info)
        ZeroMomentum.setPeriod(2000)  # (period)
        app.add(ZeroMomentum)

        DInfo = gala.DumpInfo(all_info, comp_info, 'data.log')
        DInfo.setPeriod(int(1000))
        app.add(DInfo)

        xml = gala.XMLDump(all_info, 'particles')
        xml.setPeriod(int(5e2))  # (period)
        xml.setOutputBond(True)
        xml.setOutputAngle(True)
        xml.setOutputDihedral(True)
        xml.setOutputVelocity(True)
        xml.setOutputMass(True)
        xml.setOutputCharge(True)
        xml.setOutputBody(True)
        app.add(xml)

        # ready ro run

        app.run(int(5e2 * 1e3))

    @staticmethod
    def gen_force_field_file(path, filename="ahdh.force_field"):
        ahdh_params = {'A': [0.504, 0.0011162643859539204], 'R': [0.655999, 0.7249915947715212],
                       'N': [0.568, 0.43832729970272843], 'D': [0.558, 0.029182123776349763],
                       'C': [0.548, 0.610362354303913], 'Q': [0.602, 0.3268188050525212],
                       'E': [0.592, 0.006100281608649786], 'G': [0.45, 0.7012713677972457],
                       'H': [0.608, 0.46519480823469783], 'I': [0.618, 0.6075268330845265],
                       'L': [0.618, 0.5563020305733198], 'K': [0.636, 0.058617173158697924],
                       'M': [0.618, 0.7458993420826714], 'F': [0.636, 0.9216959832175945],
                       'P': [0.555999, 0.37296418535993475], 'S': [0.518, 0.46485701300656046],
                       'T': [0.562, 0.5379777613307019], 'W': [0.678, 0.9844235478393931],
                       'Y': [0.646, 0.9950108229594323], 'V': [0.586000, 0.41850068525598694],
                       "AD": [0.504, 0.729729729729730], "RD": [0.65599, 0.0],
                       "ND": [0.568, 0.432432432432432], "DD": [0.558, 0.378378378378378],
                       "CD": [0.548, 0.594594594594595], "QD": [0.602, 0.513513513513514],
                       "ED": [0.592, 0.459459459459459], "GD": [0.45, 0.648648648648649],
                       "HD": [0.608, 0.513513513513514], "ID": [0.618, 0.972972972972973],
                       "LD": [0.618, 0.972972972972973], "KD": [0.636, 0.513513513513514],
                       "MD": [0.618, 0.837837837837838], "FD": [0.636, 1.0],
                       "PD": [0.555999, 1.0], "SD": [0.518, 0.594594594594595],
                       "TD": [0.562, 0.675675675675676], "WD": [0.678, 0.945945945945946],
                       "YD": [0.646, 0.864864864864865], "VD": [0.586000, 0.891891891891892],
                       "Ph": [0.611, 0.459459], "Su": [0.611, 0.756757],
                       "Ab": [0.611, 0.351351], "Gb": [0.611, 0.540541],
                       "Cb": [0.611, 0.297297], "Tb": [0.611, 0.594595]}

        with open(os.path.join(path, filename), "w") as f_force_field:
            f_force_field.write("<ah_params>\n")
            for key, value in list(ahdh_params.items())[:20]:
                f_force_field.write(f"{key:5} {value[0]:<10} {value[1]:<30}\n")

        return os.path.join(path, filename)