# 🧬 AI-Driven Drug Discovery Pipeline for HIF-2α Modulators

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Contributions Welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/Socks2109/drug_discovery_pipeline/issues)

## 🎯 Table of Contents

- [🎯 Objective](#-objective)
- [📋 Overview](#-overview)
- [🏗️ Project Structure](#️-project-structure)
- [🔧 Installation & Setup](#-installation--setup)
- [🚀 Quick Start](#-quick-start)
- [📊 Pipeline Components](#-pipeline-components)
- [📈 Results](#-results)
- [🔬 DEMO](#-demo)
- [📚 Documentation](#-documentation)
- [🤝 Contributing](#-contributing)
- [📄 License](#-license)
- [🙏 Acknowledgments](#-acknowledgments)

## 🎯 Objective

To design and implement an **AI-driven pipeline** for identifying novel small-molecule modulators targeting the **HIF-2α dimerization interface**, leveraging cutting-edge computational chemistry, machine learning, and large-scale screening technologies.

## 📋 Overview

This project integrates artificial intelligence with computational drug discovery to prioritize and optimize candidates targeting biologically significant protein-protein interactions. Our modular pipeline demonstrates how emerging AI technologies, particularly Large Language Models (LLMs), can augment traditional docking-based drug discovery methods.

### 🧪 Key Technologies
- **MolMIM**: Molecular generation using Mutual Information Machine
- **DiffDock**: AI-powered molecular docking
- **TxGemma**: Toxicity and pharmacokinetic prediction
- **SEA**: Similarity Ensemble Approach for library diversification
- **ZINC & TLDR**: Chemical compound databases

### 🎯 Target
**HIF-2α (Hypoxia-Inducible Factor 2-alpha)** - A critical transcription factor involved in cellular oxygen sensing and cancer progression.

## 🏗️ Project Structure

```
Drug_Discovery_Project/
├── 📁 01 Library Generation/        # Compound library generation using MolMIM
│   ├── SEA_MolMIM_SMILELibraryGeneration/
│   └── smiles/                      # SMILES validation and processing
│
├── 📁 02 Docking Pipeline/          # Molecular docking with DiffDock
│   ├── AutoDOCK Pipeline/           # Alternative docking methods
│   └── results/                     # Docking results and analysis
│
├── 📁 03 TxGemma Pipeline/          # Toxicity & PK property prediction
│   ├── Agentic_TxGemma_Toxicity.ipynb
│   └── Agentic_TxGemma_PK_Prop.ipynb
│
├── 📁 04 Final Workflow/            # Integrated pipeline execution
│
├── 📁 DEMO/                         # Complete demonstration workflow
│   ├── DEMO.ipynb                  # Step-by-step pipeline demo
│   ├── ligand_images/              # Molecular visualizations
│   └── results_output/             # Example outputs
│
├── 📁 Final_Results/                # Curated final compound recommendations
│   ├── top_37_ligand_smiles_actives.txt
│   ├── top_19_ligand_smiles_inactives.txt
│   └── top_40935_ligand_scores.txt
│
├── 📁 Some_Helper_Functions/        # Utility scripts
└── 📁 Some_Results/                 # Intermediate analysis results
```

## 🔧 Installation & Setup

### Prerequisites
- Python 3.8+
- Conda (recommended)
- NVIDIA API access (for MolMIM and DiffDock)
- RDKit

### Environment Setup

```bash
# Clone the repository
git clone https://github.com/Socks2109/drug_discovery_pipeline.git
cd drug_discovery_pipeline

# Create conda environment
conda create -n drug-discovery python=3.10.11
conda activate drug-discovery

# Install RDKit
conda install -c rdkit rdkit

# Install other dependencies
pip install -r requirements.txt
```

### API Configuration

1. **NVIDIA API Key**: Create a `.env` file with your NVIDIA BioNeMo API key:
   ```
   NVIDIA_API_KEY=your_api_key_here
   ```

2. **Authentication**: Update git credentials for collaboration:
   ```bash
   git config --global user.name "YourUsername"
   git config --global user.email "your.email@example.com"
   ```

## 🚀 Quick Start

### Option 1: Run the Complete DEMO
```bash
cd DEMO/
jupyter notebook DEMO.ipynb
```

### Option 2: Step-by-Step Execution
```bash
# 1. Generate compound library
cd "01 Library Generation/SEA_MolMIM_SMILELibraryGeneration/"
python LibraryGeneration_MolMIM_BatchProcessing.py

# 2. Convert SMILES to SDF and run docking
cd "../../02 Docking Pipeline/"
python before_diffdock_smiles_to_sdf.py
python diffdock_submit.py

# 3. Evaluate toxicity and PK properties
cd "../03 TxGemma Pipeline/"
jupyter notebook Agentic_TxGemma_Toxicity.ipynb
```

## 📊 Pipeline Components

### 🧪 Stage 1: Library Generation
- **Input**: Known HIF-2α modulators (SMILES)
- **Process**: 
  - SEA analysis for compound diversification
  - MolMIM-based molecular generation
  - Chemical space exploration
- **Output**: Diverse compound libraries

### 🎯 Stage 2: Molecular Docking
- **Input**: Generated compound libraries + HIF-2α structure
- **Process**:
  - SMILES to 3D structure conversion (RDKit)
  - DiffDock AI-powered pose prediction
  - Binding affinity scoring
- **Output**: Ranked compounds by binding affinity

### 🔍 Stage 3: ADMET Filtering
- **Input**: Top-ranked compounds
- **Process**:
  - Toxicity prediction (AMES, ClinTox)
  - PK property evaluation (F, t½, VDss)
  - Lipinski's Rule of Five filtering
- **Output**: Drug-like, non-toxic candidates

### 📋 Stage 4: Final Selection
- **Input**: Filtered candidates
- **Process**: Multi-criteria optimization
- **Output**: Final compound recommendations

## 📈 Results

### 🏆 Key Achievements
- **40,935** compounds initially screened
- **37** high-confidence active compounds identified
- **19** inactive compounds for comparison
- **Average confidence score**: 0.95+ for top candidates

### 📊 Performance Metrics
- **AMES Mutagenicity AUROC**: 0.798
- **ClinTox AUROC**: 0.831
- **Oral Bioavailability AUROC**: 0.655
- **Half-life Spearman Correlation**: 0.494
- **VDss Spearman Correlation**: 0.607

### 🧬 Top Compounds
The final recommended compounds show:
- High binding affinity to HIF-2α
- Favorable ADMET properties
- Drug-like characteristics
- Low toxicity profiles

*Detailed results available in `Final_Results/` directory*

## 🔬 DEMO

The `DEMO/` directory contains a complete walkthrough:
- **DEMO.ipynb**: Interactive Jupyter notebook
- **Sample Data**: 36 ligands with full pipeline results
- **Visualizations**: Molecular structures and binding poses
- **Performance Analysis**: Comprehensive metrics and comparisons

## 📚 Detailed Documentation

- [Library Generation Guide](./01%20Library%20Generation/SEA_MolMIM_SMILELibraryGeneration/readme.md)
- [Docking Pipeline Manual](./02%20Docking%20Pipeline/README.md)
- [TxGemma Documentation](./03%20TxGemma%20Pipeline/README.md)
- [AutoDock Documentation](https://github.com/JoeVonDahab/drug_discovery_pipeline/blob/main/02%20Docking%20Pipeline/AutoDOCK%20Pipeline/readme.md)


## 🤝 Contributing

We welcome contributions! Please see our contributing guidelines:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

### 🐛 Issues
Please report bugs and feature requests through [GitHub Issues](https://github.com/Socks2109/drug_discovery_pipeline/issues).

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

### 🔬 Scientific References
- **MolMIM**: Nikolaus Hansen. The CMA Evolution Strategy. 2005
- **SEA**: Keiser MJ, et al. Relating protein pharmacology by ligand chemistry. Nat Biotech 25 (2), 197-206 (2007)
- **DiffDock**: NVIDIA BioNeMo Platform

### 🛠️ Tools & Libraries
- **RDKit**: Open-source cheminformatics
- **NVIDIA BioNeMo**: AI-powered drug discovery platform
- **Jupyter**: Interactive development environment

### 👥 Team
- **Isaac Yeoh** (TxGemma) — [@Socks2109](https://github.com/Socks2109)  
- **Youssef Abo-Dahab** (AutoDock) — [@JoeVonDahab](https://github.com/JoeVonDahab)  
- **Shyiang** (DiffDock, Piepline connecting) — [@abcdefucsb](https://github.com/abcdefucsb)  
- **Somayeh Motevalli** (TxGemma) — [@SomayehMotevalli](https://github.com/SomayehMotevalli)  
- **Miko Mallari** (MolMIM) — [@mf-mallari](https://github.com/mf-mallari)

---

**Pipeline Workflow:**

![Pipeline Overview](https://github.com/user-attachments/assets/2dd102ab-23d1-428f-bf5a-87998a5f2f96)

*The pipeline successfully identified 37 high-confidence HIF-2α modulators through an integrated AI-driven approach combining generative chemistry, molecular docking, and ADMET prediction.*

---

**📧 Contact**: For questions or collaboration opportunities, please open an issue or contact the development team.

**🌟 Star this repository** if you find it useful for your drug discovery research!
