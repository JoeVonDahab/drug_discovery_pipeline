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

### Usage

Place your ligand .sdf files in the sdf_output/ directory.

Ensure your receptor file is named receptor_clean.pdb and located in the root directory.

Run the script:
