import os
import sys
import glob
import shutil
import subprocess


class GromacsMDRunner:
    def __init__(self, pdb_file, mdps_dir=None):
        self.pdb_file = pdb_file
        if mdps_dir is None:
            mdps_dir = os.path.join(os.path.dirname(__file__), 'all_atom_mdp')
        
        self.mdp_files = glob.glob(os.path.join(mdps_dir, '*.mdp'))

        self.charmm36_ff_dir = os.path.join(os.path.dirname(__file__), 'charmm36.ff/')
        

        

    def build_simulation_box(self, ):
        pdb_dir = os.path.join(os.path.dirname(self.pdb_file), self.pdb_file.split(".")[0])
        os.makedirs(pdb_dir, exist_ok=True)
        subprocess.run(f"rm -f {pdb_dir}/*", shell=True)

        subprocess.run(["cp", self.pdb_file, os.path.join(pdb_dir, "protein.pdb")])
        
        os.chdir(pdb_dir)
        
        

        shutil.copytree(self.charmm36_ff_dir, './charmm36.ff/', dirs_exist_ok=True)

        for mdp_file in self.mdp_files:
            shutil.copy(mdp_file, '.')

        pdb2gmx = [
                    "gmx", "pdb2gmx",
                    "-f", "protein.pdb", "-o", "protein.gro",
                    "-water", "tip3p",
                    "-ff", "charmm36",
                    "-ter", "-ignh"]
        subprocess.run(pdb2gmx, input="4\n6\n", text=True, check=True)
        subprocess.run(["gmx", "editconf", "-f", "protein.gro", "-o newbox.gro -c -bt cubic -d 1.2"])
        subprocess.run(["gmx", "solvate", "-cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro"])
        subprocess.run(["gmx", "grompp", " -f ions.mdp -c solv.gro -p topol.top -o ions.tpr"])
        
        subprocess.run(["echo", "SOL | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname SOD -nname CLA -conc 0.15 -neutral"])
        
        subprocess.run(["gmx make_ndx -f solv_ions.gro -o index.ndx"])


    def run(self):
        ...



