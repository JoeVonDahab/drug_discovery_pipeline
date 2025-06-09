# Drug Discovery Pipeline Project

## Objective: 
To design and implement an AI-driven pipeline for identifying novel small-molecule modulators targeting the HIF-2α dimerization interface, leveraging computational chemistry, machine learning, and large-scale screening tools. 

## Description: 

This project aims to integrate artificial intelligence techniques with computational drug discovery to prioritize and optimize candidates targeting a biologically significant protein-protein interaction. Central to this effort is the use of molecular docking to evaluate binding interactions between the HIF-2α dimerization interface and both known and newly generated compounds. 

Using existing HIF-2α inhibitors and agonists as a foundation, we will construct and screen virtual libraries. Generative AI models will be employed to design novel molecules inspired by these known binders, expanding the chemical space for exploration. 

To enhance screening and prioritization, we will incorporate large language models (LLMs) and other AI approaches capable of predicting bioactivity and optimizing compound selection. This workflow will be supported by tools including TxGemma, molecular dynamics simulations, and large-scale bioactivity datasets, which will be used to validate and refine candidate molecules. 

Rather than aiming for a singular "correct" solution, the project's goal is to develop a robust, modular pipeline that demonstrates how emerging AI technologies, particularly LLMs, can augment and extend traditional docking-based drug discovery methods. This pipeline will serve as a modular solution that has the capability to screen libraries of molecules to desired protein pockets. 

## Pipeline:

![image](https://github.com/user-attachments/assets/2dd102ab-23d1-428f-bf5a-87998a5f2f96)

The pipeline utilizes MolMIM, DiffDock, TxGemma, and various libraries such as Zinc and TLDR to produce meaningful binders to the HIF-2α protein.

## Output:

Through this pipeline, we discovered 30 compounds that had a high p value that binds to the HIF-2α protein. The final compounds we decided to recommend are located in `.txt`.
