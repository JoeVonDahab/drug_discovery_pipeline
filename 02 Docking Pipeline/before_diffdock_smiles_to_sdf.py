# from rdkit import Chem
# from rdkit.Chem import AllChem

# # Input SMILES string
# smiles = "CCN(CC)C(=O)N[C@@H]1CCN2CCc3c(C)cccc3[C@@]12C"


# mol = Chem.MolFromSmiles(smiles)


# mol = Chem.AddHs(mol)
# AllChem.EmbedMolecule(mol, AllChem.ETKDG())
# AllChem.UFFOptimizeMolecule(mol)

# writer = Chem.SDWriter("ligand.sdf")
# writer.write(mol)
# writer.close()

# print("SDF written to ligand.sdf")

import os
from rdkit import Chem
from rdkit.Chem import AllChem

# Create output directory
output_dir = "sdf_output"
os.makedirs(output_dir, exist_ok=True)

# Read SMILES from file
with open("../docked_compounds_2.txt", "r") as file:
    smiles_list = [line.strip() for line in file if line.strip()]

# Generate SDF files
for idx, smiles in enumerate(smiles_list):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Skipping invalid SMILES at line {idx + 1}: {smiles}")
        continue
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
        print(f"Embedding failed for SMILES at line {idx + 1}: {smiles}")
        continue
    AllChem.UFFOptimizeMolecule(mol)
    
    sdf_path = os.path.join(output_dir, f"ligand_{idx}.sdf")
    writer = Chem.SDWriter(sdf_path)
    writer.write(mol)
    writer.close()

print(f"Finished writing {len(smiles_list)} SDF files to '{output_dir}'")

