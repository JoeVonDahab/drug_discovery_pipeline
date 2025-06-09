input_file = r"C:\Users\zhaol\Downloads\drug_discovery_pipeline\Docking Pipeline\results_output\top_40935_ligand_smiles.txt"   # Replace with your actual file name
output_file = "DiffDock_Top_4000.txt"

with open(input_file, "r") as f:
    lines = f.readlines()

# Skip the header and extract SMILES strings from tab-separated values
smiles_lines = [line.strip().split('\t')[1] for line in lines[1:4001] if '\t' in line]

# Write SMILES strings to the output file with header
with open(output_file, "w") as f:
    f.write("SMILES\n")
    for smile in smiles_lines:
        f.write(smile + "\n")

print(f"{len(smiles_lines)} SMILES strings written to {output_file} with header 'SMILES'")

