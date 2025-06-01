## **Docking Results**

This analysis ranks all docked compounds by their **predicted binding energy** (in kcal/mol), as computed from AutoDock-GPU `.dlg` output files. The compounds are also annotated with:

- Their **SMILES structure**, matched by compound name from the file `50k.smi`
- Whether their docked pose lies **within the expected pocket**, based on whether the corresponding filtered PDBQT file exists (typically generated when the pose intersects the defined docking box)

### Pocket Definition

A docking box was used to define the binding site. Only ligands whose docked pose fell within this region were considered for export. The pocket is visualized below:

![Docking Box View 1](https://github.com/user-attachments/assets/2fab8f0e-688c-4ded-8686-9e39ff91287b)
![Docking Box View 2](https://github.com/user-attachments/assets/f5b6a878-1d45-43ed-bef5-6462ad749a68)

### Output Files

1.  **`docking_results_ranked_docked.csv` (Main Docking Results):**
    * This is the primary CSV file containing the comprehensive list of **all** compounds that were successfully docked.
    * Compounds are ranked based on their predicted binding affinity (e.g., AutoDock-GPU score in kcal/mol), with the most favorable scores listed first.

  | Column                | Description                                                       |
  |------------------------|-------------------------------------------------------------------|
  | `Rank`               | Position in ranking (1 = best)                                     |
  | `Compound`           | Base name of the docked compound                                   |
  | `SMILES`             | SMILES string for the compound (from `50k.smi`, or "N/A" if missing) |
  | `Affinity_kcal_per_mol` | Predicted binding energy from `.dlg` file                       |
  | `DLG_file`           | Name of the `.dlg` file that was processed                         |
  | `Within_Pocket`      | `True` if the corresponding PDBQT file was found and valid         |

2.  **`docking_results_ranked_actives.csv`:**
    * contains the docking results of known active compounds.
    
3.  **`docking_results_ranked_inactives.csv`:**
     * contains the docking results of known inactive compounds but have a close structure.

4.  **`closest_100_to_actives.csv`:**
    **closest to the mean binding energy observed for the known active compounds** (from `docking_results_ranked_actives.csv`).
5.  **`top_100_ligands_simple.rar` (Archive File):**
    **PDBQT files for the top 100 ligands** that had the strongest predicted binding affinity *and* whose poses were found within the defined binding pocket.
    
6.  **`comparing_results.ipynb` (Jupyter Notebook):**
    it contains code for result analysis and visualization 
---
### **PT2385 superimposed:**
![WhatsApp Image 2025-05-31 at 23 20 55_11bc8624](https://github.com/user-attachments/assets/dc97889a-657a-44f0-ae64-caef58b76898)
 
 ### top compound imposed with crystal structure ligand:
![WhatsApp Image 2025-05-31 at 23 20 55_11bc8624](https://github.com/user-attachments/assets/eeaa7de0-2362-4128-9c91-4ab91ab5cabe)

#### Data for top 100 compound 
  ![top_100_scores](https://github.com/user-attachments/assets/54cc5395-3fe7-4431-b315-4d65db9cd85b)

### Notes

- Ligands **outside the pocket** were excluded from the top export list, even if their binding score was high.
- Missing SMILES entries (shown as "N/A") may indicate naming mismatches or compounds not present in `50k.smi`.

---

### Example Row from CSV

![image](https://github.com/user-attachments/assets/31396a9b-d5d4-482a-b352-d6c8a9a5b501)




