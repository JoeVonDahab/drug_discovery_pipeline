### SMILES to SDF Converter using RDKit
This script reads SMILES (Simplified Molecular Input Line Entry System) strings from a text file and converts them into individual 3D-optimized SDF (Structure Data File) files using RDKit.

### Features
Parses SMILES strings from an input text file

Generates 3D molecular structures

Optimizes structures using the UFF (Universal Force Field)

Outputs individual .sdf files for each valid molecule

Skips invalid SMILES and handles embedding errors gracefully

### Requirements
Python 3.7+

RDKit (must be installed with all necessary chemistry modules)

### Installation
Use conda to install RDKit (recommended):

conda create -c rdkit -n rdkit-env rdkit python=3.8
conda activate rdkit-env
Usage
Prepare a file named inactives.txt containing one SMILES string per line.

### Run the script:


python smiles_to_sdf.py
Output .sdf files will be saved in the sdf_output_inactives/ directory.

### Example Input
inactives.txt:

CCO
CCC(=O)O
CCN(CC)C(=O)N[C@@H]1CCN2CCc3c(C)cccc3[C@@]12C
Output
sdf_output_inactives/ligand_0.sdf

sdf_output_inactives/ligand_1.sdf

...

### Notes
The script uses the ETKDG algorithm for 3D structure embedding.

Molecules that fail SMILES parsing or 3D embedding are skipped with a warning.

### License
This project is open-source and free to use under the MIT License.
