# nwPropagation

![nwPropagation Logo]()

**nwProgation** is a network propagation framework designed to increase your power in identifying candiate marker genes in the situation of small sample sizes.

---

## 🚀 Key Features

*   **Five networks**: BioPlex3, HumanNet, PCNet, ProteomeHD, STRING
*   **Random walk network propagation**: PageRank


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
docker run -it --rm -v $(pwd):/work yuanli202004/autoencoder:v.1.02 bash
```

---

## 🛠 Getting Started

### 1. Prepare Your Input Data
nwPropagation requires a csv file containing your differential statistics, including columns:
*   **Pvalue**
*   **NCBI_ID**


### 2. Basic Usage

```python
from nwPropagation import networkPropagation

    
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

