# Purpose

Description of SMILE generation utilizing literature search, the Similarity Ensemble Approach, and molecular generation using Mutual Information Machine (MoIMIM) for diversification of known drugs. 

# Contents

The following table describes the contents of this project:

| Item                                         | Description                                                                                                |
| :------------------------------------------- | :--------------------------------------------------------------------------------------------------------- |
| `LibraryGeneration_MoIMIM_SingleSubmission.ipynb` | Input: 1 SMILE string Output: A table of novel SMILES with physical properties                             |
| `LibraryGeneration_MoIMIM_BatchProcessing.py`   | Input: .smi file Output: A table of novel SMILES with physical properties                                 |
| `Folder: SMILES from MolMIM`                   | Internal folder containing SMILES titled with name of starting molecule                                     |
| `molmim_combined_results_[TargetName].csv` `molmim_combined_results_[TargetName]-world.csv` | Indicates SMILES sourced from High Affinity Substances list found in SEA results Indicates SMILES sourced from Drug Substances list found in SEA results |
| `.env` file                                    | Create your own .env file containing your API Key                                                          |

# Similarity Ensemble Approach

Purpose: to use the SEA to diversify the library of small molecules that was discovered during a literature search. 

* Took molecule SMILES from our literature list and individually searched SEA 
* Of the molecules that pass the P-Value threshold of (10^-15) were considered 
* If the molecule's ZINC results page had drug substances, a .smi file was saved. 
* If the molecule's ZINC results page had no drug substances matched, molecules from the high affinity .smi files were saved from the Highest Affinity Substances list 
* Drug substances were saved because they have been noted to have biological relevance. 

# MolMIM

Purpose: To implement a controlled iterative method of molecule generation 

## Set up:

1.  Create an .env file containing your NVIDIA API key and place it in the same directory as the scripts. 

## NVIDIA API Call Parameters

The following table describes the parameters used for the NVIDIA API call: 

| API parameter  | Current options | Description                                                                                                |
| :------------- | :-------------- | :--------------------------------------------------------------------------------------------------------- |
| `smi`          | SMILE string    | Simplified molecular input line entry system, or the reference molecule                                      |
| `Algorithm`    | CMA-ES None     | Specifies the optimization algorithm used.                                                                   |
| `num_molecules` | 1-99            | Number of molecules to generate                                                                              |
| `Property_name` | QED pLogP       | The molecular property to optimize                                                                           |
| `minimize`     | True (1) False (0) | Indicates whether property\_name should be minimized or maximized. (True: minimize, False: maximize)         |
| `min_similarity` | 0-1             | Sets a minimum similarity threshold to the reference molecule                                                |
| `particles`    | 2-1000          | Number of candidate molecules evaluated in each generation                                                   |
| `iterations`   | 2-1000          | Number of steps in the optimization. You may consider monitoring convergence for each molecule input          |
| `Scaled_radius` | 0-1             | Extent of deviation from reference molecule in the latent space                                              |

## How to use:

### Single submission notebook

Expected inputs:

* SMILE Strings 
* NVIDIA API Call parameters (See table 1.) 

Expected outputs:

* Table of generated molecules with physical property scores Data frame saved as a CSV file 

To generate a molecule:

1.  In the Starting Molecule cell change the input molecule 
2.  Adjust API call 
3.  Consider monitoring convergence (Not yet implemented in this notebook) 

**(I want to insert an image here. Please leave a reminder)**

<img width="316" alt="SingleSubmissionSMILEstring" src="https://github.com/user-attachments/assets/90f48de3-159b-4f02-80f8-bc453b6df4e5" />

### Batch Processing

Purpose: To leverage an automated system for SMILE generation. 

Expected inputs:

* .smi file 
* NVIDIA API Call parameters (See table 1.) 

Expected outputs:

* Table of generated molecules with physical property scores Data frame saved as a CSV file 

To generate molecules:

1.  Change the input .smi file (at the end of the script) 
2.  Adjust API call 
3.  Consider monitoring convergence (Not yet implemented in this notebook) 

# Sources

**MolMIM**

* Nikolaus Hansen. The CMA Evolution Strategy. 2005 
* A. Auger, N. Hansen: Tutorial CMA-ES: Evolution Strategies and Covariance Matrix Adaptation. 2012 

**SEA**

* Keiser MJ, Roth BL, Armbruster BN, Ernsberger P, Irwin JJ, Shoichet BK. Relating protein pharmacology by ligand chemistry. Nat Biotech 25 (2), 197-206 (2007) 

Workshop notebook

[https://github.com/hw-ju/tutorial\_bionemo\_nim?tab=readme-ov-file](https://github.com/hw-ju/tutorial_bionemo_nim?tab=readme-ov-file)
