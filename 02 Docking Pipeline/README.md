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

# DiffDock Batch Submission Script using NVIDIA BioNeMo API

This script automates the process of submitting multiple ligand `.sdf` files and a single protein `.pdb` receptor to NVIDIA's BioNeMo DiffDock API for molecular docking and pose prediction.

---

## Features

- Uploads ligands and a receptor to the NVIDIA Cloud via the Asset API
- Sends inference requests to the BioNeMo DiffDock endpoint
- Handles API rate limiting (HTTP 429) with exponential backoff
- Saves response status and body for each ligand in a dedicated output directory
- Processes ligands concurrently using multithreading for efficiency

---

## Directory Structure

Before running, make sure your project folder contains:

project/
├── receptor_clean.pdb # Your protein receptor file
├── sdf_output/ # Folder with ligand SDF files
│ ├── ligand_0.sdf
│ ├── ligand_1.sdf
│ └── ...
├── results_output/ # Output folder (auto-created)
└── diffdock_submit.py # This script


---

## Requirements

- Python 3.7+
- `requests` library (install via `pip install requests`)

---

## Authentication

The script uses a bearer token in the header:

```python
header_auth = "Bearer <your_nvidia_api_token>"

Replace the placeholder with your NVIDIA BioNeMo API token.
```
### Usage

Place your ligand .sdf files in the sdf_output/ directory.

Ensure your receptor file is named receptor_clean.pdb and located in the root directory.

Run the script:

python diffdock_submit.py

Each ligand's results will be saved in:

results_output/ligand_<index>/
├── response_status.txt
├── request_url.txt
└── response_text.txt

### Notes
Ligands that fail upload or inference will be skipped with an error message.

API rate limits (HTTP 429) are handled with automatic retries.

By default, a maximum of 3 threads is used to reduce the risk of hitting rate limits.

### License
This script is provided "as-is" under the MIT License. You are responsible for securing and managing your NVIDIA API token properly.

### Support
For API documentation or further support, visit NVIDIA BioNeMo.

# Docking Pipeline - Ligand Confidence Extraction

This Python script processes docking results to extract ligand poses, save them into appropriate formats (PDB and SDF), and calculate confidence scores. It identifies the top 100 ligands with the highest confidence scores and provides a summary of the results.

## Features

- **Input**: JSON file containing ligand docking results (`response_text.txt`).
- **Outputs**:
  - PDB files for each ligand pose.
  - SDF files for each ligand position.
  - A summary file containing the confidence scores for each pose.
  - A ranking of the top 100 ligands based on their confidence scores.

## Prerequisites

Ensure that you have the required Python environment with the necessary packages. This script assumes that the input data follows a particular structure and the files are stored in a specific directory layout.

- Python 3.x
- Required Libraries:
  - `json`
  - `os`

## Directory Structure

The script expects the following directory structure:

<base_path>
│
├── ligand_<i>
│ ├── response_text.txt
│ └── diffdock_actual_outcome
│ ├── pose_1.pdb
│ ├── ligand_pose_1.sdf
│ ├── ...
│ └── pose_confidences.txt
│
└── top_100_ligands_by_confidence.txt


Where:
- `<base_path>`: The root folder containing the ligand subfolders (e.g., `ligand_0`, `ligand_1`, etc.).
- `response_text.txt`: JSON file containing the docking results for each ligand.
- `diffdock_actual_outcome`: Output directory for each ligand containing PDB, SDF files, and confidence scores.

## Code Explanation

### Step 1: Read Input Data
The script reads `response_text.txt` (JSON) from each ligand folder. This file contains:
- **Trajectory** (`data["trajectory"]`): PDB poses for each ligand.
- **Ligand Positions** (`data["ligand_positions"]`): SDF representations of ligand positions.
- **Confidence Scores** (`data["position_confidence"]`): Confidence values for the ligand poses.

### Step 2: Write Output Files
For each ligand:
1. PDB files (`pose_X.pdb`) are written into the `diffdock_actual_outcome` folder.
2. SDF files (`ligand_pose_X.sdf`) are written for each ligand's position.
3. A file `pose_confidences.txt` is generated with a rank and corresponding confidence for each pose.

### Step 3: Calculate Top 100 Ligands by Confidence
- For each ligand, the script extracts the highest valid confidence score.
- Ligands are ranked based on their highest confidence, and the top 100 ligands are selected.
- An average confidence score is computed for the top 100 ligands.

### Step 4: Output Summary
The results of the top 100 ligands are written into a summary file `top_100_ligands_by_confidence.txt` that includes:
- Ligand ID
- Corresponding confidence score
- Average confidence score for the top 100 ligands

### Example Output:

Top 100 ligands by highest confidence score:

Ligand ID Confidence Score
ligand_0 0.9843
ligand_1 0.9711
...
ligand_99 0.9374

Average confidence score: 0.9512


### Step 5: Summary and Printing
The script prints out the number of ligands selected and the average confidence score.

## Running the Script

1. Modify the `base_path` variable to point to the location of your ligand folders.
2. Ensure that the structure of each ligand folder matches the expected layout.
3. Execute the script in your Python environment.

## Notes

- The script automatically skips any ligand folder that does not contain the required `response_text.txt` file.
- The highest confidence for each ligand is selected, ignoring any `None` values in the confidence list.
- The script sorts the ligands by confidence score in descending order and selects the top 100.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


