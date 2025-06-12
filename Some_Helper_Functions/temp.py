

##############################################

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read all three files
df1 = pd.read_csv(r'C:\Users\zhaol\Downloads\drug_discovery_pipeline_1\Final_Results\top_37_ligand_scores_actives.txt', sep='\t')
df2 = pd.read_csv(r'C:\Users\zhaol\Downloads\drug_discovery_pipeline_1\Final_Results\top_40935_ligand_scores.txt', sep='\t')
df3 = pd.read_csv(r'C:\Users\zhaol\Downloads\drug_discovery_pipeline_1\Final_Results\top_19_ligand_scores_inactives.txt', sep='\t')  # Add path to your third CSV file

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


# input_file = r"C:\Users\zhaol\Downloads\drug_discovery_pipeline_1\DEMO\Lipinski_after.txt"
# output_smiles = r"C:\Users\zhaol\Downloads\drug_discovery_pipeline_1\Final_Results\top_19_ligand_smiles_inactives.txt"
# output_scores = r"C:\Users\zhaol\Downloads\drug_discovery_pipeline_1\Final_Results\top_19_ligand_scores_inactives.txt"

# # Open input file for reading
# with open(input_file, "r") as infile:
#     lines = infile.readlines()

# # Skip header and process remaining lines
# data_lines = lines[1:]

# # Prepare output content
# smiles_entries = []
# score_entries = []

# for line in data_lines:
#     parts = line.strip().split("\t")
#     if len(parts) < 4:
#         continue  # Skip malformed lines
#     ligand_number = parts[1]
#     smiles = parts[2]
#     score = parts[3]

#     smiles_entries.append(f"{ligand_number}\t{smiles}")
#     score_entries.append(f"{ligand_number}\t{score}")

# # Write Ligand Number and SMILES
# with open(output_smiles, "w") as f:
#     f.write("Ligand Number\tSMILES\n")
#     f.write("\n".join(smiles_entries))

# # Write Ligand Number and Confidence Score
# with open(output_scores, "w") as f:
#     f.write("Ligand Number\tConfidence Score\n")
#     f.write("\n".join(score_entries))


