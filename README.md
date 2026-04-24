# Single-Nucleus Transcriptomics Analysis Identifies Sex-Dichotomous Pathological Axes in Alzheimer’s Disease: Female Homeostatic Failure versus Male Neuroimmune Activation

This repository contains the complete R analytical pipeline and scripts used in our study.

---

## 📋 Overview

This repository provides a comprehensive workflow for analyzing single-nucleus RNA sequencing (snRNA-seq) data to investigate sex-specific differences in Alzheimer's disease (AD). The analysis covers:

- **Quality Control**: Filtering and normalization of `Seurat` objects.
- **Cell Type Annotation**: Classification of excitatory/inhibitory neurons, astrocytes, microglia, oligodendrocytes, and OPCs.
- **Sex-Stratified DE Analysis**: AD vs. Control comparisons within male and female cohorts.
- **Pathway Enrichment**: GO/KEGG for DEGs on both sex; GSEA and GSVA focusing on estrogen signaling and metabolic axes.
- **Microglial Analysis**: Assessment of functional exhaustion and homeostatic markers.
- **Synaptic Dysfunction**: Specialized analysis of excitatory neuron DEGs.
- **Cell-Cell Communication**: Interactome analysis using `CellChat`.
- **Visualization**: Scripts for all main and supplementary figures.

---

## 💻 System Requirements

### Operating System & Environment
- **Tested OS**: macOS Ventura, Ubuntu 22.04 LTS, Windows 11.
- **R Version**: >= 4.1.0 (Developed with R 4.3.1).

### Required R Packages
| Package | Version | Purpose |
| :--- | :--- | :--- |
| `Seurat` | >= 4.3.0 | Core snRNA-seq data handling |
| `tidyverse` | >= 2.0.0 | Data manipulation and visualization |
| `clusterProfiler` | >= 4.6.0 | Gene set enrichment analysis |
| `msigdbr` | >= 7.5.1 | MSigDB gene sets retrieval |
| `GSVA` | >= 1.46.0 | Gene set variation analysis |
| `CellChat` | >= 1.6.1 | Ligand-receptor interaction analysis |
| `patchwork` | >= 1.1.2 | Figure composition |
| `pheatmap` | >= 1.0.12 | Heatmap generation |

> [!NOTE]
> A complete list of all packages and exact versions is available in [`session_info.txt`](session_info.txt).

### Hardware Recommendations
- **RAM**: Minimum 16 GB required; 32 GB+ recommended for large-scale GSVA/CellChat processing.
- **Storage**: ~5 GB for the integrated Seurat objects and outputs.

---

## 🔧 Installation

### 1. Clone the Repository
```bash
git clone https://github.com/atticuszhou288-maker/Sex-difference-analysis-of-AD-using-snRNA-seq-cohort.git
cd Sex-difference-analysis-of-AD-using-snRNA-seq-cohort
```

### 2. Install R Dependencies
Open R or RStudio in the repository directory and run the following script:

```r
# Install CRAN packages
cran_pkgs <- c("tidyverse", "patchwork", "pheatmap", "RColorBrewer", 
               "ggpubr", "ggrepel", "corrplot", "igraph", "ggraph")
install.packages(cran_pkgs)

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

bioc_pkgs <- c("Seurat", "clusterProfiler", "msigdbr", "GSVA", "CellChat")
BiocManager::install(bioc_pkgs)
```

> [!IMPORTANT]
> **CellChat** may require additional system dependencies (e.g., `NMF`, `ComplexHeatmap`). If installation fails, please refer to the [CellChat GitHub](https://github.com/sqjin/CellChat) instructions.

---

## 📂 Repository Structure

```text
.
├── README.md                # Project documentation
├── LICENSE                  # MIT License
├── .gitignore               # Files excluded from version control
├── session_info.txt         # Complete R environment snapshot
├── data/                    # Directory for input .rds files (User provided)
├── scripts/                 # Sequential R scripts (F01-09,analysis_process,supp)

```

---

## 📊 Data Availability

### Public Source Datasets
The raw single-cell RNA sequencing data used in this study are publicly available:
* **Mathys et al.**: [Synapse: syn18485175](https://www.synapse.org/#!Synapse:syn18485175)
* **Zhou et al.**: [Synapse: syn21125841](https://www.synapse.org/#!Synapse:syn21125841)
* **Lau et al.**: [GEO: GSE157827](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157827)

### Integrated Seurat Data
Integrated Seurat data is available at:
* **DOI**: [https://doi.org/10.17632/f5kx9k264n.1](https://doi.org/10.17632/f5kx9k264n.1)
* **Reference**: *Li et al. Cell 188, 5516–5534, October 2, 2025*

### Required Input Files
To reproduce the analysis, place the following in the `data/` directory:
1.  `mzl_seurat.rds`: The processed Seurat object.
2.  `metadata.csv`: Metadata table (see Supplementary Table 1 of the article).

---

## 🚀 How to Reproduce the Analysis

### Step-by-Step Execution
1.  **Environment**: Ensure all R packages are installed (see [Installation](#installation)).
2.  **Data Placement**: Move the Seurat object into the `data/` folder.
3.  **Working Directory**: Set your R working directory to the repository root.
4.  **Execution**: Run scripts in numerical order to maintain data flow:

```r
# Example Sequence
source("scripts/01_quality_control.R")
source("scripts/02_cell_type_annotation.R")
source("scripts/03_deg_analysis.R")
# ... continue
```

---

## 📉 Expected Outputs

After successful execution, the following will be generated in the `output/` directory:

| File/Folder | Description |
| :--- | :--- |
| `all_de_genes_results.rds` | Comprehensive Differential Expression (DE) results |
| `significant_degs.csv` | List of statistically significant DEGs |
| `DEG_heatmap.png` | Heatmap visualization of top DEGs |
| `GSVA_Results/` | Folder containing GSVA scores and effect size plots |
| `Microglia_Functional.png`| Functional scores for microglia |
| `synaptic_volcano_plot.png`| Volcano plot of synaptic DEGs |
| `CellChat_*` | Cell-cell communication networks and data |
...and many more figures as presented in the paper.
---

## 🔍 Reproducibility Notes

* **Session Info**: The `session_info.txt` file contains the exact R environment used to produce final results.
* **Random Seeds**: Seeds are set within each script to ensure consistent results for stochastic processes (e.g., clustering, UMAP).

---

## 📄 License & Contact

### License
This code is released under the **MIT License**. See the `LICENSE` file for full terms.

### Contact
For questions regarding the code or analysis, please contact the first author:
* **Email**: [2310319@mail.nankai.edu.cn]

*For questions about raw data, please refer to source GEO/Synapse pages or contact the original study authors.*
