# nwPropagation

![nwPropagation Logo](https://raw.githubusercontent.com/Maj18/nwPropagation/refs/heads/main/doc/logo.tiff)
![Version](https://img.shields.io/badge/version-1.0.2-brightgreen) 
![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22833704.svg)(https://doi.org/10.5281/zenodo.22833704)

**nwPropagation** is a Python package that provides a flexible framework for propagating signals from differentially altered molecular features across biological interaction networks. It supports multiple network resources, Random Walk with Restart, hyperparameter selection, and visualization of propagation results.

---

## 🚀 Key Features

- Integrates five biological networks: **BioPlex3**, **HumanNet**, **PCNet**, **ProteomeHD**, and **STRING**.
- Implements **Random Walk with Restart (RWR)** for network propagation.
- Uses **power iteration** to calculate stable node-visiting probabilities.
- Provides a framework for selecting the optimal **restart probability ($\alpha$)**.
- Includes built-in functions for **visualizing propagation results**.



## Why Network Propagation

Traditional differential expression (DE) analyses evaluate proteins as independent entities, frequently failing to capture the systemic nature of complex skin pathologies. 
To overcome these limitations, we deployed network propagation via Random Walk with Restart (RWR) for three primary reasons:

1. **Prioritization of Sub-threshold Drivers:** It rescues biologically critical genes that fail to pass strict False Discovery Rate (FDR) thresholds due to low transcript abundance or statistical noise, provided they reside within a high-density disease neighborhood.
2. **Capture of Post-Translational Signaling:** It identifies key signaling hubs and structural anchors that are regulated via physical interactions or post-translational modifications (e.g., phosphorylation) rather than transcriptional changes.
3. **Functional Module Discovery:** By projecting statistical values onto a physical protein-protein interaction (PPI) network, it filters out isolated false positives and highlights highly coordinated, mechanistically linked sub-networks (e.g., cell-adhesion and extracellular matrix remodeling complexes).


---

## 📦 Installation

### Standard Installation
For researchers working locally (macOS/Linux):

```bash
# Clone the repository
git clone https://github.com/Maj18/nwPropagation.git
cd nwPropagation

# Install in editable mode
python -m pip install --no-build-isolation -v -e .
```

### Docker (Recommended)
Use our pre-built environment to avoid dependency conflicts. Requires Docker v23.0+.

```bash
docker run -it --rm -v $(pwd):/work yuanli202004/autoencoder:v.1.0.2 bash
```

---

## 🛠 Getting Started

### 1. Prepare Your Input Data
nwPropagation requires a csv file containing your differential statistics, including columns:
*   **Pvalue**
*   **NCBI_ID**
*   **Gene**
*   **Uniprot**


### 2. Basic Usage

```python
from nwPropagation import *

# Convert your features, e.g. from Uniprot to Gene symbols and NCBI_ID

# Run RWR 

# Visualize the results
    
```

## 🔬 Resources & Data
*   The package is developed based on [GWAS_NetworkPropagation](https://github.com/gvisona/GWAS_NetworkPropagation.git)

## 🤝 Contribution
Contributions are welcome! 
*   **Bugs**: Open a [GitHub Issue](https://github.com).
*   **Features**: For new utility functions or neuronal datasets, please open an issue for discussion before submitting a Pull Request.

## 👨‍💻 Developer
**Yuan Li**  
📧 [uncork-shady-next@duck.com](mailto:uncork-shady-next@duck.com), NBIS, Science for Life Laboratory (SciLifeLab), Sweden.

## 📄 License
This project is licensed under the **BSD 3-Clause License** - see the [LICENSE](LICENSE) file for details.

