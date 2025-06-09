# Final Workflow

## File Description
Below is a description of the main code files that were used to predict which ligands were most likely to be a HIF-2α inhibitor.

| File Name | Service in the Pipeline        | Description                                                              |
|-----------|-------------------------------|--------------------------------------------------------------------------|
| 0.ipynb   | MolMIM                        | Uses CMA-ES algorithm to optimize for high QED (drug-likeness)           |
| 1.ipynb   | DiffDock and Lipinski's rule | Checks for ligands predicted to bind strongly and pass Lipinski's rule   |
| 3.ipynb   | TxGemma Toxicity             | Checks for top ligands' toxicity                                         |
| 4.ipynb   | TxGemma PK prediction        | Predicts PK properties                                                   |

