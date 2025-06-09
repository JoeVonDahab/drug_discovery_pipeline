# # Define the input and output file paths
# input_file = r"C:\Users\zhaol\Downloads\drug_discovery_pipeline\Docking Pipeline\results_output_inactives\highest_confidences.txt"
# output_file = r"C:\Users\zhaol\Downloads\drug_discovery_pipeline\Docking Pipeline\results_output_inactives\highest_confidences_ranked.txt"

# # Open and read the input file
# with open(input_file, "r") as f:
#     lines = f.readlines()

# # Open the output file for writing
# with open(output_file, "w") as f:
#     f.write("Ligand\tConfidence Score\n")
#     count = 1
#     for line in lines:
#         line = line.strip()
#         # Skip header or empty lines
#         if not line or not line.replace('.', '', 1).isdigit():
#             continue
#         f.write(f"Ligand {count}\t{line}\n")
#         count += 1

# print(f"Labeled confidence scores written to '{output_file}'")
#############################
# import pandas as pd

# # File path
# file_path = r"C:\Users\zhaol\Downloads\drug_discovery_pipeline\Docking Pipeline\results_output_actives\highest_confidences_ranked.txt"

# # Read the file
# df = pd.read_csv(file_path, sep="\t")

# # Sort by confidence score in descending order
# df_sorted = df.sort_values(by="Confidence Score", ascending=False).reset_index(drop=True)

# # Add rank column if needed (optional)
# df_sorted.index += 1  # Start index from 1 for human-readable ranking
# df_sorted.insert(0, "Rank", df_sorted.index)

# # Overwrite the original file
# df_sorted.to_csv(file_path, sep="\t", index=False)

# # Calculate and print percentiles
# percentiles = [10, 25, 50, 75, 90]
# cutoffs = df["Confidence Score"].quantile([p / 100 for p in percentiles])

# print("Confidence Score Cutoffs:")
# for p, val in zip(percentiles, cutoffs):
#     print(f"{p}% cutoff: {val}")
#######################################################
# import pandas as pd

# # Read the file (tab-separated)
# df = pd.read_csv(r'C:\Users\zhaol\Downloads\drug_discovery_pipeline\Docking Pipeline\results_output_actives\highest_confidences_ranked.txt', sep='\t')

# # Define the percentiles to calculate
# percentiles = [10, 25, 50, 75, 90]

# # Calculate percentiles for the Confidence Score column
# results = df['Confidence Score'].quantile([p / 100 for p in percentiles])

# # Print results
# for p, value in zip(percentiles, results):
#     print(f"{p}th percentile: {value:.4f}")

#######################################################
#Define file paths
# confidence_file = r"C:\Users\zhaol\Downloads\drug_discovery_pipeline\Docking Pipeline\results_output\top_40935_ligand_scores.txt"
# smiles_file = r"C:\Users\zhaol\Downloads\drug_discovery_pipeline\Docking Pipeline\results_output\top_40935_ligand_smiles.txt"
# output_file = "Lipinski_before.txt"

# # Threshold
# threshold = 0.911488950252533

# qualified_ligands = set()
# with open(confidence_file, "r") as f:
#     next(f)  # Skip header
#     for line in f:
#         parts = line.strip().split()
#         if len(parts) == 2:
#             ligand, score = parts
#             if float(score) > threshold:
#                 qualified_ligands.add(ligand)

# # Step 2: Match ligand numbers to SMILES and collect only the SMILES
# smiles_list = []
# with open(smiles_file, "r") as f:
#     next(f)  # Skip header
#     for line in f:
#         parts = line.strip().split(maxsplit=1)
#         if len(parts) == 2:
#             ligand, smiles = parts
#             if ligand in qualified_ligands:
#                 smiles_list.append(smiles)

# # Step 3: Write only SMILES strings to output
# with open(output_file, "w") as f:
#     for smiles in smiles_list:
#         f.write(smiles + "\n")

# print(f"Extracted {len(smiles_list)} SMILES strings to '{output_file}'")

##############################################

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read all three files
df1 = pd.read_csv(r'C:\Users\zhaol\Downloads\drug_discovery_pipeline\Docking Pipeline\results_output_actives\highest_confidences_ranked.txt', sep='\t')
df2 = pd.read_csv(r'C:\Users\zhaol\Downloads\drug_discovery_pipeline\Docking Pipeline\results_output\top_40935_ligand_scores.txt', sep='\t')
df3 = pd.read_csv(r'C:\Users\zhaol\Downloads\drug_discovery_pipeline\Docking Pipeline\results_output_inactives\highest_confidences_ranked.txt', sep='\t')  # Add path to your third CSV file

# Extract confidence scores
scores1 = df1['Confidence Score']
scores2 = df2['Confidence Score']
scores3 = df3['Confidence Score']  # Confidence Score from the third dataset

# Compute quantiles from the first dataset
quantiles = scores1.quantile([0.10, 0.25, 0.50, 0.75, 0.90])

# Setup the figure with 3 subplots
fig, axes = plt.subplots(3, 1, figsize=(10, 15), sharex=True)  # Increase the figure height to fit 3 plots

# Plot distribution for file1
sns.histplot(scores1, bins=15, kde=True, color='skyblue', edgecolor='black', ax=axes[0])
axes[0].set_title('Distribution of Confidence Scores - Known Inhibitors (Actives)')
axes[0].set_ylabel('Frequency')

# Mark quantiles on first plot
for q, val in quantiles.items():
    axes[0].axvline(val, color='red', linestyle='--')
    axes[0].text(val, axes[0].get_ylim()[1]*0.9, f'{int(q*100)}%: {val:.3f}', rotation=90,
                 verticalalignment='center', color='red')

# Plot distribution for file2
sns.histplot(scores2, bins=50, kde=True, color='lightgreen', edgecolor='black', ax=axes[1])
axes[1].set_title('Distribution of Confidence Scores - Library')
axes[1].set_ylabel('Frequency')

# Apply same quantile lines from file1 to second plot
for q, val in quantiles.items():
    axes[1].axvline(val, color='red', linestyle='--')
    axes[1].text(val, axes[1].get_ylim()[1]*0.9, f'{int(q*100)}%', rotation=90,
                 verticalalignment='center', color='red')

# Plot distribution for file3
sns.histplot(scores3, bins=50, kde=True, color='blue', edgecolor='black', ax=axes[2])
axes[2].set_title('Distribution of Confidence Scores - Inactives')
axes[2].set_xlabel('Confidence Score')
axes[2].set_ylabel('Frequency')

# Apply same quantile lines from file1 to third plot
for q, val in quantiles.items():
    axes[2].axvline(val, color='red', linestyle='--')
    axes[2].text(val, axes[2].get_ylim()[1]*0.9, f'{int(q*100)}%', rotation=90,
                 verticalalignment='center', color='red')

plt.tight_layout()
plt.show()


#####################################
# import pandas as pd
# from rdkit import Chem
# from rdkit.Chem import Draw

# # Step 1: Load confidence scores and assign ranking
# df_conf = pd.read_csv(r'C:\Users\zhaol\Downloads\drug_discovery_pipeline\Docking Pipeline\results_output\top_40935_ligand_scores.txt', sep='\t')
# df_conf['Rank'] = df_conf['Confidence Score'].rank(method='first', ascending=False).astype(int)

# # Filter ligands with confidence > 0.9115
# filtered = df_conf[df_conf['Confidence Score'] > 0.9115]

# # Step 2: Load SMILES data
# df_smiles = pd.read_csv(r'C:\Users\zhaol\Downloads\drug_discovery_pipeline\Docking Pipeline\results_output\top_40935_ligand_smiles.txt', sep='\t')

# # Step 3: Merge to include SMILES and ranking
# merged = pd.merge(filtered, df_smiles, on='Ligand Number', how='inner')

# # Step 4: Generate and save images with rank-based filenames
# for _, row in merged.iterrows():
#     rank = row['Rank']
#     ligand_number = row['Ligand Number']
#     confidence = row['Confidence Score']
#     smiles = row['SMILES']

#     mol = Chem.MolFromSmiles(smiles)
#     if mol:
#         img = Draw.MolToImage(
#             mol, size=(300, 300),
#             legend=f'Rank: {rank}\nLigand: {ligand_number}\nScore: {confidence:.4f}'
#         )
#         img.save(f'rank_{rank}_ligand_{ligand_number}.png')

##########################################################################

# import pandas as pd
# import numpy as np

# # Load the data (assuming tab-separated)
# df = pd.read_csv(r'C:\Users\zhaol\Downloads\drug_discovery_pipeline\Docking Pipeline\results_output_actives\highest_confidences_ranked.txt', sep='\t')

# # Calculate the standard error of the mean (SEM)
# sem = np.std(df['Confidence Score'], ddof=1) / np.sqrt(len(df))

# # Print result
# print(f"Standard Error of the Mean (SEM): {sem:.6f}")
#############################################################

# import pandas as pd

# # Load the SMI file (SMILES \t ZINC ID)
# smi_df = pd.read_csv(r"C:\Users\zhaol\Downloads\drug_discovery_pipeline\actives.txt", sep='\t', names=["SMILES", "ZINC_ID"])

# # Load Ligand Number to SMILES mapping
# ligand_smiles_df = pd.read_csv(r"C:\Users\zhaol\Downloads\drug_discovery_pipeline\Docking Pipeline\results_output\top_40935_ligand_smiles.txt", sep='\t')

# # Load Confidence Scores
# confidence_df = pd.read_csv(r"C:\Users\zhaol\Downloads\drug_discovery_pipeline\Docking Pipeline\results_output\top_40935_ligand_scores.txt", sep='\t')

# # Merge ligand_smiles with confidence scores on Ligand Number
# merged_df = pd.merge(confidence_df, ligand_smiles_df, on="Ligand Number")

# # Merge the result with the SMI file on SMILES
# final_df = pd.merge(merged_df, smi_df, on="SMILES", how="left")

# # Add ranking based on descending confidence score
# final_df['Rank'] = final_df['Confidence Score'].rank(ascending=False, method='min').astype(int)

# # Sort by rank
# final_df.sort_values(by="Rank", inplace=True)

# # Reorder columns
# final_df = final_df[['Rank', 'Ligand Number', 'SMILES', 'ZINC_ID', 'Confidence Score']]

# # Output to CSV
# final_df.to_csv("ranked_ligands.csv", index=False)
############################################################

# import pandas as pd

# # Read the original CSV file into a DataFrame
# csv_file = r'C:\Users\zhaol\Downloads\drug_discovery_pipeline\ranked_ligands.csv'  # Replace with your actual CSV file
# df = pd.read_csv(csv_file)

# # Read the SMILES data from the text file
# with open(r'C:\Users\zhaol\Downloads\drug_discovery_pipeline\actives_2.txt', 'r') as f:
#     smiles_lines = f.readlines()

# # Read the confidence scores from the text file
# confidence_file = r'C:\Users\zhaol\Downloads\drug_discovery_pipeline\Docking Pipeline\results_output_actives\highest_confidences_ranked.txt'  # Replace with the actual confidence scores file
# confidence_df = pd.read_csv(confidence_file, delimiter='\t')

# # Create a dictionary to store ligand names and their corresponding confidence scores
# confidence_dict = {row['Ligand']: row['Confidence Score'] for index, row in confidence_df.iterrows()}

# # Iterate over the SMILES file and add the new entries with "known inhibitors" label
# new_entries = []

# for line in smiles_lines:
#     # Split on the colon to separate "Ligand N" from the SMILES string
#     parts = line.strip().split(":")
#     if len(parts) < 2:  # Skip any malformed lines
#         continue
    
#     ligand_name = parts[0].strip()  # Extract the ligand name, e.g., "Ligand 1"
#     smiles = parts[1].strip()  # Extract the SMILES string
    
#     try:
#         # Use the ligand number to find the confidence score
#         ligand_num = int(ligand_name.split()[1])  # Extract ligand number, e.g., 1 for "Ligand 1"
#     except (IndexError, ValueError):
#         print(f"Skipping line due to incorrect format: {line}")
#         continue
    
#     # Add confidence score from the confidence_dict
#     confidence_score = confidence_dict.get(ligand_name, None)
    
#     # Prepare the new entry
#     new_entry = {
#         'Rank': len(df) + len(new_entries) + 1,  # Increment rank by the current number of rows + 1
#         'Ligand Number': ligand_num,
#         'SMILES': smiles,
#         'ZINC_ID': f'ZINC{str(ligand_num).zfill(12)}',  # Use a placeholder for ZINC_ID (adjust if needed)
#         'Confidence Score': confidence_score,
#         'Label': 'known inhibitors' if ligand_name in confidence_dict else 'unknown',  # Label newly inserted entries as "known inhibitors" or "unknown"
#         'Known Inhibitor': 'Yes' if ligand_name in confidence_dict else 'No'  # New column for "known inhibitors" status
#     }
    
#     # Append the new entry to the list of new entries
#     new_entries.append(new_entry)

# # Convert the new entries into a DataFrame
# new_entries_df = pd.DataFrame(new_entries)

# # Append the new entries to the original DataFrame
# df_updated = pd.concat([df, new_entries_df], ignore_index=True)

# # Save the updated DataFrame to a new CSV file
# output_file = r'C:\Users\zhaol\Downloads\drug_discovery_pipeline\ranked_ligands_updated_file.csv'  # Replace with desired output file name
# df_updated.to_csv(output_file, index=False)

# print("New entries have been successfully added and saved.")
####################################################################






