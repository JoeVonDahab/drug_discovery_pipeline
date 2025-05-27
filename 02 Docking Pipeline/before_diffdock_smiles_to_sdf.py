from rdkit import Chem
from rdkit.Chem import AllChem

# Input SMILES string
smiles = "CCN(CC)C(=O)N[C@@H]1CCN2CCc3c(C)cccc3[C@@]12C"


mol = Chem.MolFromSmiles(smiles)


mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
AllChem.UFFOptimizeMolecule(mol)

writer = Chem.SDWriter("ligand.sdf")
writer.write(mol)
writer.close()

print("SDF written to ligand.sdf")
