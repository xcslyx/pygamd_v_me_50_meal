from pdbfixer import PDBFixer
from openmm.app import PDBFile

base_name = "Gb"
fixer = PDBFixer(f'{base_name}.pdb')
# fixer.findMissingResidues()
# fixer.findMissingAtoms()
fixer.addMissingHydrogens(pH=7.0)

# 平移到几何中心
import numpy as np
positions = fixer.positions
coords = np.array([[pos.x, pos.y, pos.z] for pos in positions])
center = coords.mean(axis=0)
coords -= center

for i, pos in enumerate(fixer.positions):
    pos.x, pos.y, pos.z = coords[i]

PDBFile.writeFile(fixer.topology, fixer.positions, open(f'{base_name}_output.pdb','w'))
