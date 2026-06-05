import os
import sys
import subprocess


class GromacsMDRunner:
    def __init__(self, pdb_file, mdps_dir=None):
        self.pdb_file = pdb_file
        if mdps_dir is None:
            self.mdps_dir = os.path.join(os.path.dirname(__file__), 'all_atom_mdp')
        else:
            self.mdps_dir = mdps_dir

    def build_simulation_box(self, ):
        pdb_dir = os.path.join(os.path.dirname(self.pdb_file), self.pdb_file.split(".")[0])
        os.makedirs(pdb_dir, exist_ok=True)
        os.chdir(pdb_dir)

        # os.subprocess.run(f"rm -f {pdb_dir}/*")

        os.subprocess.run(f"cp {self.pdb_file} {os.path.join(pdb_dir, 'protein.gro')}")
        os.subprocess.run(f"rm {os.path.join(pdb_dir, 'topol.top')}")

        os.subprocess.run(f"cp -r {os.path.join(os.path.dirname(__file__), 'charmm36.ff/')} {pdb_dir}")
        os.subprocess.run(f"cp {os.path.join(self.mdps_dir, '*')} {pdb_dir}")

        os.subprocess.run(f'printf "4\n6\n" | gmx pdb2gmx -f protein.gro -o protein.gro -water tip3p -ff charmm36 -ter -ignh')
        os.subprocess.run(f"gmx editconf -f protein.gro -o newbox.gro -c -bt cubic -d 1.2")
        os.subprocess.run(f"gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro")
        os.subprocess.run(f"gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr")
        
        os.subprocess.run(f'echo "SOL" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname SOD -nname CLA -conc 0.15 -neutral')
        
        os.subprocess.run(f"gmx make_ndx -f solv_ions.gro -o index.ndx")


    def run(self):
        ...



