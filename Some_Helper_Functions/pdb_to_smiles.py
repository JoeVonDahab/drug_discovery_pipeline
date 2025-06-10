from rdkit import Chem
from rdkit.Chem import AllChem
import os

def pdb_to_smiles(pdb_path):
    if not os.path.isfile(pdb_path):
        print(f"Warning: File {pdb_path} not found.")
        return None

    try:
        mol = Chem.MolFromPDBFile(pdb_path, removeHs=False)
        if mol is None:
            print(f"Error: RDKit could not parse {pdb_path}.")
            return None

        AllChem.Compute2DCoords(mol)
        smiles = Chem.MolToSmiles(mol)
        return smiles

    except Exception as e:
        print(f"Exception during SMILES conversion for {pdb_path}: {e}")
        return None

if __name__ == "__main__":
    input_file = r"C:\Users\zhaol\Downloads\drug_discovery_pipeline\Docking Pipeline\results_output\top_100_ligand_scores.txt"
    output_file = r"C:\Users\zhaol\Downloads\drug_discovery_pipeline\Docking Pipeline\results_output\top_100_molecule_smiles.txt"

    with open(input_file, 'r') as f:
        lines = f.readlines()[1:]  # skip header

    with open(output_file, 'w') as out:
        for line in lines:
            ligand_number = line.strip().split()[0]
            pdb_path = fr"C:\Users\zhaol\Downloads\drug_discovery_pipeline\Docking Pipeline\results_output\ligand_{ligand_number}\diffdock_actual_outcome\pose_1.pdb"
            smiles = pdb_to_smiles(pdb_path)
            if smiles:
                out.write(f"{smiles}\n")

