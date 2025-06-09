# Read input file and extract SMILES strings
smiles_list = []
with open("1m.smi", "r") as infile:  # Replace with your actual input filename
    for line in infile:
        if line.strip():  # Skip empty lines
            smiles = line.strip().split()[0]  # Extract the SMILES string
            smiles_list.append(smiles)

# Write extracted SMILES to a new file
with open("docked_compounds_1_million.txt", "w") as outfile:
    for smiles in smiles_list:
        outfile.write(smiles + "\n")

print("SMILES strings have been saved to smiles_only.txt")
