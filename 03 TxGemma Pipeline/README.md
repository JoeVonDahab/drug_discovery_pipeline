# TxGemma Pipeline

This folder contains the files required to run TxGemma. The 2 agents that were developed are:

### 1) ```Agentic_TxGemma_Toxicity.ipynb```
Evaluates the toxicity of the SMILE strings after passing the binding pocket fit. 2 metrics were used for toxicity evaluation:

- AMES Mutagenicity: Test whether the chemical is mutagenic and therefore may act as a carcinogen. It has a reported AUROC of 0.798.
- ClinTox: Test whether the drug is toxic based on previously failed clinical trials. It has a reported AUROC of 0.831.

Combined, these 2 along with a pubmed search is fed to the agentic TxGemma for evaluation. If the drug is determined to be toxic, it will get filtered out. If the drug is determined to be not toxic, it will pass into the next agent: ```Agentic_TxGemma_PK_Prop.ipynb```.

### 2) ```Agentic_TxGemma_PK_Prop.ipynb```
Evaluates the PK properties of the drug. Due to TxGemma's generation, we were unable to reverse the normalization of the clearance for hepatocytes and microsomes. However, we are able to extract these 3 PK properties for each SMILE string provided:

- Oral Bioavailability (F)
    - Reported AUROC: 0.655
- Half Life (t<sub>1/2</sub>)
    - Reported Spearman Correlation: 0.494
- Steady State Volume of Distribution (VD<sub>ss</sub>)
    - Reported Spearman Correlation: 0.607

To run these 2 files, please create a conda environment using:
```
conda create -n myenv python=3.10.11
conda activate myenv
```
Then install all the required packages in ```requirements.txt``` using:
```
pip install -r requirements.txt
```
Our legacy files are also located in the folder ```Legacy Files``` if you wish to view what we have tried previously, but decided that it would not fit in our current workflow.