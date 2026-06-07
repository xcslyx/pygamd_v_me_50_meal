import os
import sys
import glob
import shutil
import subprocess
from tkinter import W

from Bio.PDB import PDBParser


class GromacsMDRunner:
    def __init__(self, pdb_file, mdps_dir=None):
        self.pdb_file = pdb_file
        self.molecule_type = self.analyze_molecule_type(self.pdb_file)
        self.pdb_dir = os.path.join(os.path.dirname(self.pdb_file), self.pdb_file.split(".")[0])
        if mdps_dir is None:
            mdps_dir = os.path.join(os.path.dirname(__file__), 'all_atom_mdp')
        
        self.mdp_files = glob.glob(os.path.join(mdps_dir, '*.mdp'))

        self.charmm36_ff_dir = os.path.join(os.path.dirname(__file__), 'charmm36.ff/')


    @staticmethod
    def analyze_molecule_type(pdb_file):
        # 标准残基集合
        protein_residues = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 
            'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 
            'SER', 'THR', 'TRP', 'TYR', 'VAL'
        }
        # RNA标准残基 (包含一些旧版或常见的修饰写法)
        rna_residues = {'A', 'U', 'C', 'G', 'ADE', 'URA', 'CYT', 'GUA'}
        # DNA标准残基 (仅作对比参考)
        dna_residues = {'DA', 'DT', 'DC', 'DG'}
        
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('molecule', pdb_file)

        # 遍历结构中的所有残基
        molecule_dict = {"protein": False, "rna": False, "dna": False}
        for residue in structure.get_residues():
            # 获取残基名称并去除可能存在的空格
            res_name = residue.get_resname().strip().upper()
            
            if res_name in protein_residues:
                molecule_dict["protein"] = True
                break
            elif res_name in rna_residues:
                molecule_dict["rna"] = True
                break
            elif res_name in dna_residues:
                molecule_dict["dna"] = True

        return next((key for key, value in molecule_dict.items()), None)


    @staticmethod
    def change_content_by_line(file, line_idx, new_content):
        with open(file, 'r') as f:
            lines = f.readlines()
        lines[line_idx] = new_content + '\n'
        with open(file, 'w') as f:
            f.writelines(lines)


    def build_simulation_box(self, ):
        os.makedirs(self.pdb_dir, exist_ok=True)
        subprocess.run(f"rm -f {self.pdb_dir}/*", shell=True)

        subprocess.run(["cp", self.pdb_file, os.path.join(self.pdb_dir, "protein.pdb")])
        os.chdir(self.pdb_dir)
        shutil.copytree(self.charmm36_ff_dir, './charmm36.ff/', dirs_exist_ok=True)

        for mdp_file in self.mdp_files:
            shutil.copy(mdp_file, '.')

        pdb2gmx = ["gmx", "pdb2gmx", "-f", "protein.pdb", "-o", "protein.gro",
            "-water", "tip3p", "-ff", "charmm36", "-ter", "-ignh"]
        if self.molecule_type == "protein":
            pdb2gmx_input = "2\n1\n"
        elif self.molecule_type == "rna":
            pdb2gmx_input = "4\n6\n"
        subprocess.run(pdb2gmx, input=pdb2gmx_input, text=True, check=True)
        
        edit_conf = ["gmx", "editconf", "-f", "protein.gro", "-o", "newbox.gro",
            "-c", "-bt", "cubic", "-d", "1.2"]
        with open("edit_conf.log", "w") as log_file:
            subprocess.run(edit_conf, stdout=log_file, stderr=subprocess.STDOUT, check=True)
        
        solvate = ["gmx", "solvate", "-cp", "newbox.gro", "-cs", "spc216.gro",
            "-p", "topol.top", "-o", "solv.gro"]
        with open("solvate.log", "w") as log_file:
            subprocess.run(solvate, stdout=log_file, stderr=subprocess.STDOUT, check=True)
        
        grompp = ["gmx", "grompp", "-f", "ions.mdp",
            "-c", "solv.gro", "-p", "topol.top", "-o", "ions.tpr"]
        with open("grompp.log", "w") as log_file:
            subprocess.run(grompp, stdout=log_file, stderr=subprocess.STDOUT, check=True)
        
        genion = ["gmx", "genion", "-s", "ions.tpr", "-o", "solv_ions.gro", "-p", "topol.top",
            "-pname", "SOD", "-nname", "CLA", "-conc", "0.15", "-neutral"]
        with open("genion.log", "w") as log_file:
            subprocess.run(genion, stdout=log_file, input="SOL\n", text=True, check=True)
        
        make_ndx = ["gmx", "make_ndx", "-f", "solv_ions.gro", "-o", "index.ndx"]
        if self.molecule_type == "protein":
            make_ndx_input = "q\n"
        # elif self.molecule_type == "rna":
        #     make_ndx_input = "4\n"
        subprocess.run(make_ndx, input=make_ndx_input, text=True, check=True)


    def run_simulation(self, temperature=300.0):
        grompp_minim = ["gmx", "grompp", "-f", "minim.mdp", "-c", "solv_ions.gro",
            "-p", "topol.top", "-o", "em.tpr", "-maxwarn", "1", "-n", "index.ndx"]
        with open("grompp_minim.log", "w") as log_file:
            subprocess.run(grompp_minim, stdout=log_file, stderr=subprocess.STDOUT, check=True)
        mdrun_minim = ["gmx", "mdrun", "-v", "-deffnm", "em", "-nt", "16", "-ntmpi", "16"]
        subprocess.run(mdrun_minim)
        
        # protein_tc_grps = "Protein Non-Protein"
        rna_tc_grps = "RNA !RNA"
        
        if self.molecule_type == "rna":
            self.change_content_by_line("nvt.mdp", 31-1, f"tc-grps                 = {rna_tc_grps}")
        if temperature != 300.0:
            self.change_content_by_line("nvt.mdp", 33-1, f"ref_t                   = {temperature}     {temperature}")
        grompp_nvt = ["gmx", "grompp", "-f", "nvt.mdp", "-c", "em.gro", "-r", "em.gro",
            "-p", "topol.top", "-o", "nvt.tpr", "-maxwarn", "1", "-n", "index.ndx"]
        with open("grompp_nvt.log", "w") as log_file:
            subprocess.run(grompp_nvt, stdout=log_file, stderr=subprocess.STDOUT, check=True)
        mdrun_nvt = ["gmx", "mdrun", "-v", "-deffnm", "nvt", "-gpu_id", "3"]
        subprocess.run(mdrun_nvt)
        
        if self.molecule_type == "rna":
            self.change_content_by_line("npt.mdp", 31-1, f"tc-grps                 = {rna_tc_grps}")
        if temperature != 300.0:
            self.change_content_by_line("npt.mdp", 33-1, f"ref_t                   = {temperature}     {temperature}")
        grompp_npt = ["gmx", "grompp", "-f", "npt.mdp", "-c", "nvt.gro", "-r", "nvt.gro",
            "-t", "nvt.cpt", "-p", "topol.top", "-o", "npt.tpr", "-maxwarn", "1", "-n", "index.ndx"]
        with open("grompp_npt.log", "w") as log_file:
            subprocess.run(grompp_npt, stdout=log_file, stderr=subprocess.STDOUT, check=True)
        mdrun_npt = ["gmx", "mdrun", "-v", "-deffnm", "npt", "-gpu_id", "3"]
        subprocess.run(mdrun_npt)
        
        grompp_md = ["gmx", "grompp", "-f", "mdrun.mdp", "-c", "npt.gro", "-t", "npt.cpt",
            "-p", "topol.top", "-o", "md.tpr", "-maxwarn", "1", "-n", "index.ndx"]
        with open("grompp_md.log", "w") as log_file:
            subprocess.run(grompp_md, stdout=log_file, stderr=subprocess.STDOUT, check=True)
        # mdrun_md = ["gmx", "mdrun", "-v", "-deffnm", "md", "-nt", "3"]
        # subprocess.run(mdrun_md)



