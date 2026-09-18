# nwPropagation

![nwPropagation Logo](https://raw.githubusercontent.com/Maj18/nwPropagation/refs/heads/main/doc/logo.tiff)

**nwProgation** is a network propagation framework designed to increase your power in identifying candiate marker genes in the situation of small sample sizes.

---

## 🚀 Key Features

*   **Five networks**: BioPlex3, HumanNet, PCNet, ProteomeHD, STRING
*   **Random Walk with Restart network propagation**


## Why Network Propagation

- **Captures indirect molecular effects:** Identifies genes or proteins influenced through interaction networks, even when they are not directly detected as significant.

- **Improves biological signal:** Combines signals from related genes or proteins, helping distinguish pathway-level patterns from isolated statistical noise.

- **Prioritizes disease-relevant candidates:** Highlights network-connected genes that may contribute to disease mechanisms or serve as potential therapeutic targets.

- **Reveals affected pathways and modules:** Groups related molecular components to uncover biological processes associated with the phenotype.

- **Handles incomplete knowledge:** Uses information from known interactions to infer potential roles for poorly characterized genes.

- **Reduces dependence on strict significance cutoffs:** Allows moderately altered but well-connected genes to contribute to the analysis.

- **Integrates heterogeneous data:** Combines experimental results with protein–protein interaction, regulatory, pathway, or gene-association networks.

- **Identifies network proximity:** Measures whether disease-associated genes, drug targets, or differentially expressed genes are close to one another in the network.

- **Supports hypothesis generation:** Suggests candidate genes, pathways, and mechanisms for follow-up experiments.

- **Improves prioritization of omics results:** Ranks large gene or protein lists according to their network relevance rather than relying only on fold change or adjusted *p*-value.


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

