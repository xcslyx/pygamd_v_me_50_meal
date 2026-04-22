# pygamd_v_me_50_meal: A script for pygamd molecular dynamics simulation software

Since PYGAMD does not provide data processing software suitable for biological macromolecule condensate simulation, I have integrated the scripts I have used into a single script with rich options to choose from.

This software consists of two parts: file generation before molecular dynamics simulation, and file processing after molecular dynamics simulation. After obtaining a PDB file, you can use this software to build the simulation system, generate running scripts, and run molecular dynamics simulations. After obtaining the results, you can continue to use this software to process the XML files obtained from the simulation.

---
## Installation
You can directly install using `pip`:

    pip install git+https://gitee.com/lyxlyxlyxxx/pygamd_v_me_50_meal.git

Or

    git clone https://gitee.com/lyxlyxlyxxx/pygamd_v_me_50_meal.git
    cd pygamd_v_me_50_meal
    pip install -e .

---
## File Generation Before Molecular Dynamics Simulation

### Generating XML Files
Generate coarse-grained XML files for simulation from PDB files, currently supporting the generation of protein, DNA, and RNA coarse-grained models.

#### Usage
    v50 -pdb2xml /path/to/pdb_file.pdb
    v50 -p /path/to/system -pdb2xml filename.pdb

For example, to convert **1kx5.pdb** in the `/home/protein` folder to an **xml** file, enter

    v50 -pdb2xml /home/protein/1kx5.pdb
    or v50 -p /home/protein -pdb2xml 1kx5.pdb

This will generate the **1kx5.xml** file.

##### The following functions are presented as interactive options during runtime
- Setting up Elastic Network Model (ENM)
- Selecting DNA/RNA coarse-grained model
- Generating simulation scripts
- Setting box size (box is a cube, default size is 100nm, use -box_size option to set other lengths)
- Setting rigid body domains
- Setting separate particle types for domains

##### Protein coarse-grained model: one amino acid coarse-grained into one particle: Dignon G L, Zheng W, Kim Y C, Best R B, Mittal J. PLOS Computational Biology, 2018, 14(1): e1005941
- **HPS model**: suitable for HPS-Urry, CALVADOS series force fields
- **Mpipi model**: another protein coarse-grained model

##### DNA coarse-grained models: currently supporting two models
1. 3SPN model
   > Knotts T A, Rathore N, Schwartz D C, De Pablo J J. A coarse grain model for DNA[J]. _The Journal of Chemical Physics_, **2007**, 126(8): 084901.

2. Mittal 2 Bead model:
   > Kapoor U, Kim Y C, Mittal J. Coarse-Grained Models to Study Protein–DNA Interactions and Liquid–Liquid Phase Separation[J]. _Journal of Chemical Theory and Computation_, **2023**: acs.jctc.3c00525.

Note: Currently, the "Mittal 2 Bead model" has not been successfully implemented in PYGAMD.

##### RNA coarse-grained models
- **3SPN model**: extension based on DNA's 3SPN model
- **CAVADOS-RNA model**: (complete implementation not provided)

#### Advanced Options

**Parameter Description:**

| Parameter | Description | Default Value |
|-----------|-------------|---------------|
| `-box_size` | Set box size (cube box, unit: nm) | `100` |
| `-add_enm_bond` | Predefine elastic network domains, format like "1-50,51-100" | `None` |
| `-add_rigid_body` | Predefine rigid body domains, format like "1-50,51-100" | `None` |
| `-gen_run_file` | Whether to automatically generate run files | `unset` |
| `-dna_model` | Specify DNA model, 1=3SPN, 2=2BeadMittal | `None` |

Usage examples:
```bash
# Set box size
v50 -pdb2xml /path/to/pdb_file.pdb -box_size 150

# Predefine elastic network
v50 -pdb2xml /path/to/pdb_file.pdb -add_enm_bond "1-50,51-100"

# Predefine rigid body domains
v50 -pdb2xml /path/to/pdb_file.pdb -add_rigid_body "1-50,51-100"

# Automatically generate run file
v50 -pdb2xml /path/to/pdb_file.pdb -gen_run_file t

# Specify DNA model
v50 -pdb2xml /path/to/pdb_file.pdb -dna_model 1  # 1=3SPN, 2=2BeadMittal
```

---
## File Processing After Molecular Dynamics Simulation

### 1. Coordinate Extraction
Note: Before any data processing, coordinates need to be extracted first.

1. Since PYGAMD does not provide topology files, you need to label your system information.
   Please ensure that the last level of the provided path is in the format "number+molecule name-molecule length", and arranged in the order of the xml files.
   For example, if the system contains 40 A molecules of length 256, 20 B molecules of length 512, and 30 C molecules of length 729,
   the system must be named "40A-256+20B-512+30C-729" (this manual will use the `40A-256+20B-512+30C-729` folder as an example).
2. Place the XML files obtained from molecular dynamics simulation in the provided directory or in the xml folder under that directory.
3. Extract coordinates: `v50 -p 40A-256+20B-512+30C-729 -xyz`
4. #### Optional Parameters

| Parameter | Description | Default Value |
|-----------|-------------|---------------|
| `-remove_enm` | Whether to remove elastic bonds | None |
| `-remove_condensate_pbc` | Whether to remove PBC and move the largest condensate to the center of the box | None |

Usage examples:
```
v50 -p 40A-256+20B-512+30C-729 -xyz -remove_enm
v50 -p 40A-256+20B-512+30C-729 -xyz -remove_condensate_pbc
```


### 2. Calculate Contact Map
Note: When calculating intermolecular Contact Map, intra-molecular interactions have been excluded.

    v50 -p 40A-256+20B-512+30C-729 -cm

**Parameter Description:**

| Parameter | Description | Default Value |
|-----------|-------------|---------------|
| `-cm` | Whether to calculate contact map, will automatically plot after calculation | None |
| `-r_cut` | Set $r_{cut}$ for contact map calculation | 4.0 |

Usage example:
```
v50 -p 40A-256+20B-512+30C-729 -cm -r_cut 4.0
```

### 3. Calculate Rg, RMSD, RMSF

**Parameter Description:**

| Parameter | Description | Default Value |
|-----------|-------------|---------------|
| `-rg` | Whether to calculate radius of gyration (Rg) | None |
| `-rmsd` | Whether to calculate root mean square deviation (RMSD) | None |
| `-rmsf` | Whether to calculate root mean square fluctuation (RMSF) | None |
| `-ref` | Reference structure file path for RMSD and RMSF calculation | `None` |

Usage examples:
```
v50 -p 40A-256+20B-512+30C-729 -rg
v50 -p 40A-256+20B-512+30C-729 -rmsd -ref reference_file.xml
v50 -p 40A-256+20B-512+30C-729 -rmsf -ref reference_file.xml
```

During runtime, you will be asked which molecule to calculate, and it will calculate the one provided by ref

### 4. Calculate System Mass Density Distribution
Note: Please remove condensate PBC before using this function.

**Parameter Description:**

| Parameter | Description | Default Value |
|-----------|-------------|---------------|
| `-mass_density` | Whether to calculate mass density distribution | None |

Usage example:
```
v50 -mass_density
```

## Other Analyses

### Sequence Analysis

The sequence analysis module provides two analysis functions: NCPR (Net Charge Per Residue) analysis and aromaticity analysis.

**Feature Highlights:**
- Unified command-line interface with interactive selection of analysis type
- Support for directly passing sequence strings or sequence file paths
- Edge-adaptive completion algorithm to keep output length equal to input length
- Segmented coloring (NCPR: red for positive charge, blue for negative charge; aromaticity: purple)
- Area filling to enhance visualization effects
- Automatic adjustment of x-axis tick intervals to adapt to sequences of different lengths

**Usage Methods:**

1. Directly pass sequence string:
```
v50 -seq_analysis -seq "MKVDELVQGLLKQISAEELKKARNEIARQHLEKTHQDLKKDILTYLTDRQIKQLEDAFQKLLAEKTEENKLAQAVENSLGQLEEKLKEA"
```

2. Pass sequence through file path:
```
v50 -seq_analysis -seq pro_sequence.txt
```

3. Specify sliding window size (default 15):
```
v50 -seq_analysis -seq pro_sequence.txt -seq_window 20
```

4. Specify output file path:
```
v50 -seq_analysis -seq pro_sequence.txt -seq_output pro_analysis.png
```

**Parameter Description:**

| Parameter | Description | Default Value |
|-----------|-------------|---------------|
| `-seq_analysis` | Whether to perform sequence analysis | None |
| `-seq` | Sequence string or sequence file path | `None` |
| `-seq_window` | Sliding window size | `15` |
| `-seq_output` | Output image path | `Auto-generated` |

**Usage Flow:**
1. After running the command, you will be prompted to select the analysis type
2. Enter `1` to select NCPR analysis, enter `2` to select aromaticity analysis
3. After analysis is complete, a corresponding PNG format image will be generated

## Troubleshooting

### Common Problems and Solutions

1. **PDB File Conversion Failure**
   - Check if the PDB file format is correct
   - Ensure the PDB file contains complete atomic information
   - For multi-chain structures, ensure "TER" tags are used to separate chains

2. **Elastic Network Setting Error**
   - Ensure domain range settings are correct (counting starts from 1)
   - Check if residue numbers in the PDB file are continuous

3. **Sequence Analysis Failure**
   - Ensure the input sequence only contains standard amino acid codes
   - For file input, ensure the file encoding is UTF-8
   - Check if the sliding window size is reasonable (recommended 10-30)

4. **Coordinate Extraction Failure**
   - Ensure the directory naming format is correct (e.g., "40A-256+20B-512")
   - Check if XML files exist in the correct location
   - For large systems, memory limits may need to be increased

## Version Information

### v1.0.0
- Initial version
- Support for protein, DNA, RNA coarse-grained model generation
- Implementation of sequence analysis functions (NCPR and aromaticity analysis)
- Provide contact map, Rg, RMSD, RMSF calculations
- Support mass density distribution analysis

## Examples

### Example 1: Generate Protein XML File and Add Elastic Network

```bash
# Generate XML file and add elastic network
v50 -pdb2xml /path/to/protein.pdb

# During runtime, select:
# 1. Protein model: HPS
# 2. Whether to set elastic network: y
# 3. Domain range: 1-100,101-200
# 4. Whether to generate run file: y
```

### Example 2: Analyze NCPR and Aromaticity of Protein Sequence

```bash
# Analyze sequence file
v50 -seq_analysis t -seq protein_sequence.txt -seq_window 20

# During runtime, select:
# 1. NCPR analysis
# 2. Aromaticity analysis

# Two analysis result images will be generated
```

### Example 3: Analyze Simulation Results

```bash
# Extract coordinates
v50 -p 40A-256+20B-512 -xyz t -remove_condensate_pbc t

# Calculate contact map
v50 -p 40A-256+20B-512 -cm t -r_cut 4.5

# Calculate Rg
v50 -p 40A-256+20B-512 -rg t
```

## Acknowledgments

This software references the following research works:

1. Dignon G L, Zheng W, Kim Y C, Best R B, Mittal J. PLOS Computational Biology, 2018, 14(1): e1005941
2. Knotts T A, Rathore N, Schwartz D C, De Pablo J J. A coarse grain model for DNA[J]. _The Journal of Chemical Physics_, **2007**, 126(8): 084901.
3. Kapoor U, Kim Y C, Mittal J. Coarse-Grained Models to Study Protein–DNA Interactions and Liquid–Liquid Phase Separation[J]. _Journal of Chemical Theory and Computation_, **2023**: acs.jctc.3c00525.

---

**Note**: This software is for research use only. If you have any questions or suggestions, please contact the developer.