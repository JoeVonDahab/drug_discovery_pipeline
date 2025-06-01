![all_compared](https://github.com/user-attachments/assets/b1a6a5e9-ee84-475a-bea9-5cb6dba8037d)## **Docking Results**

This analysis ranks all docked compounds by their **predicted binding energy** (in kcal/mol), as computed from AutoDock-GPU `.dlg` output files. The compounds are also annotated with:

- Their **SMILES structure**, matched by compound name from the file `50k.smi`
- Whether their docked pose lies **within the expected pocket**, based on whether the corresponding filtered PDBQT file exists (typically generated when the pose intersects the defined docking box)

### Output Files

- **`docking_results_ranked_simple.csv`**:  
  A ranked list of compounds sorted by binding energy (lowest/best first), with the following columns:

  | Column                | Description                                                       |
  |------------------------|-------------------------------------------------------------------|
  | `Rank`               | Position in ranking (1 = best)                                     |
  | `Compound`           | Base name of the docked compound                                   |
  | `SMILES`             | SMILES string for the compound (from `50k.smi`, or "N/A" if missing) |
  | `Affinity_kcal_per_mol` | Predicted binding energy from `.dlg` file                       |
  | `DLG_file`           | Name of the `.dlg` file that was processed                         |
  | `Within_Pocket`      | `True` if the corresponding PDBQT file was found and valid         |

- **`top_100_ligands_simple/`**:  
  Contains the **top 100 PDBQT files** for ligands with the strongest predicted affinity **and** valid poses inside the pocket (as indicated by presence of filtered `_best_pose.pdbqt`).
* top compound imposed with crystal structure ligand:
![WhatsApp Image 2025-05-31 at 23 20 55_11bc8624](https://github.com/user-attachments/assets/eeaa7de0-2362-4128-9c91-4ab91ab5cabe)

* Data for top 100 compound 
  ![top_100_scores](https://github.com/user-attachments/assets/54cc5395-3fe7-4431-b315-4d65db9cd85b)


### Pocket Definition

A docking box was used to define the binding site. Only ligands whose docked pose fell within this region were considered for export. The pocket is visualized below:

![Docking Box View 1](https://github.com/user-attachments/assets/2fab8f0e-688c-4ded-8686-9e39ff91287b)
![Docking Box View 2](https://github.com/user-attachments/assets/f5b6a878-1d45-43ed-bef5-6462ad749a68)

### Notes

- Ligands **outside the pocket** were excluded from the top export list, even if their binding score was high.
- Missing SMILES entries (shown as "N/A") may indicate naming mismatches or compounds not present in `50k.smi`.

---

### Example Row from CSV

![image](https://github.com/user-attachments/assets/31396a9b-d5d4-482a-b352-d6c8a9a5b501)




