import os
import sys
import glob
import shutil
import subprocess


class GromacsMDRunner:
    def __init__(self, pdb_file, mdps_dir=None):
        self.pdb_file = pdb_file
        self.pdb_dir = os.path.join(os.path.dirname(self.pdb_file), self.pdb_file.split(".")[0])
        if mdps_dir is None:
            mdps_dir = os.path.join(os.path.dirname(__file__), 'all_atom_mdp')
        
        self.mdp_files = glob.glob(os.path.join(mdps_dir, '*.mdp'))

        self.charmm36_ff_dir = os.path.join(os.path.dirname(__file__), 'charmm36.ff/')

    def build_simulation_box(self, ):
        os.makedirs(self.pdb_dir, exist_ok=True)
        subprocess.run(f"rm -f {pdb_dir}/*", shell=True)

        subprocess.run(["cp", self.pdb_file, os.path.join(self.pdb_dir, "protein.pdb")])
        os.chdir(self.pdb_dir)
        shutil.copytree(self.charmm36_ff_dir, './charmm36.ff/', dirs_exist_ok=True)

        for mdp_file in self.mdp_files:
            shutil.copy(mdp_file, '.')

        pdb2gmx = ["gmx", "pdb2gmx", "-f", "protein.pdb", "-o", "protein.gro",
            "-water", "tip3p", "-ff", "charmm36", "-ter", "-ignh"]
        subprocess.run(pdb2gmx, input="2\n1\n", text=True, check=True)
        
        edit_conf = ["gmx", "editconf", "-f", "protein.gro", "-o", "newbox.gro",
            "-c", "-bt", "cubic", "-d", "1.2"]
        subprocess.run(edit_conf)
        
        solvate = ["gmx", "solvate", "-cp", "newbox.gro", "-cs", "spc216.gro",
            "-p", "topol.top", "-o", "solv.gro"]
        subprocess.run(solvate)
        
        grompp = ["gmx", "grompp", "-f", "ions.mdp",
            "-c", "solv.gro", "-p", "topol.top", "-o", "ions.tpr"]
        subprocess.run(grompp)
        
        genion = ["gmx", "genion", "-s", "ions.tpr", "-o", "solv_ions.gro", "-p", "topol.top",
            "-pname", "SOD", "-nname", "CLA", "-conc", "0.15", "-neutral"]
        subprocess.run(genion, input="SOL\n", text=True, check=True)
        
        make_ndx = ["gmx", "make_ndx", "-f", "solv_ions.gro", "-o", "index.ndx"]
        subprocess.run(make_ndx, input="q\n", text=True, check=True)


    def run_simulation(self):
        os.chdir(self.pdb_dir)

        grompp_minim = ["gmx", "grompp", "-f", "minim.mdp", "-c", "solv_ions.gro",
            "-p", "topol.top", "-o", "em.tpr", "-maxwarn", "1", "-n", "index.ndx"]
        subprocess.run(grompp_minim)
        
        mdrun_minim = ["gmx", "mdrun", "-v", "-deffnm", "em", "-nt", "16"]
        subprocess.run(mdrun_minim)
        
        grompp_nvt = ["gmx", "grompp", "-f", "nvt.mdp", "-c", "em.gro", "-r", "em.gro",
            "-p", "topol.top", "-o", "nvt.tpr", "-maxwarn", "1", "-n", "index.ndx"]
        subprocess.run(grompp_nvt)
        mdrun_nvt = ["gmx", "mdrun", "-v", "-deffnm", "nvt", "-gpu_id", "3"]
        subprocess.run(mdrun_nvt)
        
        grompp_npt = ["gmx", "grompp", "-f", "npt.mdp", "-c", "nvt.gro", "-r", "nvt.gro",
            "-t", "nvt.cpt", "-p", "topol.top", "-o", "npt.tpr", "-maxwarn", "1", "-n", "index.ndx"]
        subprocess.run(grompp_npt)
        mdrun_npt = ["gmx", "mdrun", "-v", "-deffnm", "npt", "-gpu_id", "3"]
        subprocess.run(mdrun_npt)
        
        grompp_md = ["gmx", "grompp", "-f", "mdrun.mdp", "-c", "npt.gro", "-t", "npt.cpt",
            "-p", "topol.top", "-o", "md.tpr", "-maxwarn", "1", "-n", "index.ndx"]
        subprocess.run(grompp_md)
        # mdrun_md = ["gmx", "mdrun", "-v", "-deffnm", "md", "-nt", "3"]
        # subprocess.run(mdrun_md)



