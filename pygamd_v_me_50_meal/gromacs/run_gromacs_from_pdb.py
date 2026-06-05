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

        subprocess.run(["cp", os.path.join(os.path.dirname(self.pdb_file), self.pdb_file), "protein.gro"])
        
        os.chdir(pdb_dir)
        subprocess.run(f"rm -f {pdb_dir}/*", shell=True)
        

        shutil.copy(self.charmm36_ff_dir, '.')

        for mdp_file in self.mdp_files:
            shutil.copy(mdp_file, '.')

        subprocess.run("printf 4\n6\n | gmx pdb2gmx -f protein.gro -o protein.gro -water tip3p -ff charmm36 -ter -ignh", shell=True)
        subprocess.run(["gmx", "editconf", "-f", "protein.gro", "-o newbox.gro -c -bt cubic -d 1.2"])
        subprocess.run(["gmx", "solvate", "-cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro"])
        subprocess.run(["gmx", "grompp", " -f ions.mdp -c solv.gro -p topol.top -o ions.tpr"])
        
        subprocess.run(["echo", "SOL | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname SOD -nname CLA -conc 0.15 -neutral"])
        
        subprocess.run(["gmx make_ndx -f solv_ions.gro -o index.ndx"])


    def run(self):
        ...



