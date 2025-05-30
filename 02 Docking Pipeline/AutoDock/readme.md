**Results of Docking Can Be Found on this LINK:**
**https://drive.google.com/drive/folders/1-Przbd8SiuCFw2KAwpXHUNoe6rlcoDq-?usp=drive_link**

# **Targeted Molecular Docking Workflow for HIF-2α**

This document outlines the step-by-step procedure followed to perform targeted molecular docking of ligands against specific residues of the Hypoxia-Inducible Factor 2-alpha (HIF-2α) PAS-B domain. The workflow utilizes UCSF Chimera/ChimeraX for receptor preparation and residue selection, Meeko for preparing receptor PDBQT files and targeted Grid Parameter Files (GPFs), AutoGrid4 for grid map generation, and AutoDock-GPU for the docking calculations.

---

## **1. Software and Environment Setup**

The following software and environment configuration was used:

* **Conda Environment:** A dedicated Conda environment was created to manage dependencies, specifically to resolve version compatibility issues.
    * Python version: 3.10
    * Environment name (example): `meeko_py310`
        ```bash
        # Command to create the environment (example)
        conda create -n meeko_py310 python=3.10
        conda activate meeko_py310
        ```
* **Meeko & RDKit:** Installed into the Conda environment from the `conda-forge` channel. This step was crucial to resolve an initial `ModuleNotFoundError: No module named 'rdkit.six'`.
    ```bash
    conda install -c conda-forge meeko rdkit
    ```
    *Initial attempts to upgrade Meeko using `pip install -U meeko` in a different environment did not resolve the `rdkit.six` issue, necessitating the creation of a new, clean Conda environment with specific package versions.*
* **AutoDock-GPU:** The GPU-accelerated docking engine. (Assumed to be pre-installed and accessible).
* **AutoGrid4:** Part of the AutoDock4 suite, used for generating grid maps. (Assumed to be pre-installed and accessible).
* **UCSF Chimera / ChimeraX:** Used for molecular visualization, receptor cleaning, and selection of target residues.

---

## **2. Step-by-Step Workflow**

### **2.1: Receptor Preparation (General)**

1.  **Obtain Receptor Structure:** A PDB structure of the target receptor is required. For this workflow, the receptor structure corresponding to PDB ID `5TBM` (HIF-2α PAS-B domain) was used as a base.
2.  **Clean Receptor:** The PDB file was cleaned to remove non-essential molecules (e.g., crystallographic waters, non-relevant co-solvents, or ligands). Hydrogens were added (this is also handled by Meeko later if not done explicitly).
    * **Input file example:** `receptor_ready_5tbm.pdb` (This is the cleaned receptor PDB file).

### **2.2: Identifying and Selecting Target Residues (UCSF Chimera/ChimeraX)**

To focus the docking on a specific binding site within HIF-2α, key residues known or hypothesized to interact with inhibitors were selected. Based on structural information (e.g., from literature, PDB ID 5TBM, and the provided summary table), the following residues in the HIF-2α PAS-B domain were identified as important:

* **H293:** Direct Bond (H-bond), Allosteric Modulator
* **Y281:** Direct Bond (H-bond network, n → π\*Ar); Allosteric Modulator
* **M252:** Allosteric Modulator, Hydrophobic Pocket Constituent (initially)
* **S304:** Gating Residue, Pocket Architecture
* **G323:** Resistance Mutation Site, Pocket Lining
* **Y278:** Allosteric Modulator, Hydrophobic Pocket Constituent
* **L309:** Hydrophobic Pocket Constituent
* **V321:** Hydrophobic Pocket Constituent
* **F325:** Hydrophobic Pocket Constituent
* **N288:** Gating Residue
* **L272:** Gating Residue
* **A277:** (Primarily for agonist interaction, H-bond with backbone)
* **H248:** Pocket Boundary

**Action:**
Using UCSF Chimera or ChimeraX:
1.  The `receptor_ready_5tbm.pdb` was loaded.
2.  The residues listed above (or a relevant subset defining the entire target cavity) were selected.
3.  Only the atoms of these selected residues were saved to a new PDB file.
    * **Output custom PDB file example:** `custom_residues_5tbm_cleaned.pdb`

### **2.3: Ligand Preparation (General)**

Ligands intended for docking need to be prepared in the PDBQT format. This involves generating 3D coordinates, adding hydrogens, assigning partial charges (e.g., Gasteiger), and defining rotatable bonds.
* This can be done using Meeko's `mk_prepare_ligand.py` script or AutoDockTools (ADT).
    ```bash
    # Example command (not run during our troubleshooting, but part of a full workflow)
    # mk_prepare_ligand.py -i ligand.sdf -o prepared_ligand.pdbqt
    ```

### **2.4: Generating a Focused Grid Parameter File (GPF) with Meeko**

Meeko's `mk_prepare_receptor.py` script was used to prepare the full receptor PDBQT and generate a GPF specifically targeted to the selected residues.

**Command Executed:**
```bash
mk_prepare_receptor.py --read_pdb receptor_ready_5tbm.pdb -o myreceptor_targeted -p -g --box_enveloping custom_residues_5tbm_cleaned.pdb --padding 5.0
