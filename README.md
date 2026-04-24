# Code for: Single-Nucleus Transcriptomics Analysis Identifies Sex-Dichotomous Pathological Axes in Alzheimer’s Disease: Female Homeostatic Failure versus Male Neuroimmune Activation

This repository contains the R analysis code used in our study:

---

## 📋 Overview

This repository provides a complete pipeline for analyzing single-cell RNA sequencing (scRNA-seq) data to investigate sex-specific differences in Alzheimer's disease (AD). The analysis covers:

- Quality control and filtering of scRNA-seq data (`Seurat` object)
- Cell type annotation (excitatory/inhibitory neurons, astrocytes, microglia, oligodendrocytes, OPCs)
- Sex-stratified differential expression (DE) analysis between AD and control groups
- Pathway enrichment analysis (GSEA, GSVA) focusing on estrogen signaling and metabolic pathways
- Microglial functional exhaustion assessment
- Synaptic dysfunction analysis in excitatory neurons
- Cell-cell communication analysis using `CellChat`
- Generation of all main and supplementary figures

---

## 💻 System Requirements

### Operating System
The code has been tested on **macOS Ventura** and **Ubuntu 22.04 LTS**. It should run on any system capable of running R (Windows, macOS, Linux).

### R Version
- **R >= 4.1.0** (developed with R 4.3.1)

### Required R Packages
Below is a list of essential packages. See [Installation](#installation) for detailed instructions.

| Package | Version | Purpose |
| :------ | :------ | :------ |
| `Seurat` | >= 4.3.0 | Core scRNA-seq data handling |
| `tidyverse` (dplyr, ggplot2, tidyr) | >= 2.0.0 | Data manipulation and visualization |
| `clusterProfiler` | >= 4.6.0 | Gene set enrichment analysis |
| `msigdbr` | >= 7.5.1 | MSigDB gene sets retrieval |
| `GSVA` | >= 1.46.0 | Gene set variation analysis |
| `CellChat` | >= 1.6.1 | Ligand-receptor interaction analysis |
| `patchwork` | >= 1.1.2 | Figure composition |
| `pheatmap` | >= 1.0.12 | Heatmap generation |
| `RColorBrewer` | >= 1.1-3 | Color palettes |
| `ggpubr` | >= 0.6.0 | Publication-ready plots with statistics |
| `ggrepel` | >= 0.9.3 | Repelling text labels in plots |
| `corrplot` | >= 0.92 | Correlation matrix visualization |
| `igraph` / `ggraph` | >= 1.5.0 | Network visualization |

A complete list of all loaded packages with exact version numbers is available in [`session_info.txt`](session_info.txt).

### Hardware Requirements
- **RAM**: Minimum 16 GB recommended for loading the full Seurat object (actual size: ~2-3 GB). For large-scale processing, 32 GB or more is ideal.
- **CPU**: Any modern multi-core processor.

---

## 🔧 Installation

### 1. Clone the Repository
```bash
git clone [https://github.com/atticuszhou288-maker/Sex-difference-analysis-of-AD-using-snRNA-seq-cohort.git]
cd Sex-difference-analysis-of-AD-using-snRNA-seq-cohort

2. Install Required R Packages
Open R or RStudio in the repository directory and run:

r
# CRAN packages
cran_pkgs <- c("tidyverse", "patchwork", "pheatmap", "RColorBrewer", 
               "ggpubr", "ggrepel", "corrplot", "igraph", "ggraph")
install.packages(cran_pkgs)

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

bioc_pkgs <- c("Seurat", "clusterProfiler", "msigdbr", "GSVA", "CellChat")
BiocManager::install(bioc_pkgs)
Note: If you encounter any installation issues, please refer to the official documentation of each package. The CellChat package may require additional dependencies; follow the instructions at CellChat GitHub.

📂 Repository Structure
text
.
├── README.md                 # This file
├── LICENSE                   # MIT License
├── .gitignore                # Files/folders excluded from version control
├── session_info.txt          # Complete R session information
│
├──scripts                    # Analysis scripts(Code Folder F1-9;analysis process;supp)

scripts/: Contains all executable R scripts. They are named and numbered to reflect the logical order of analysis.

📊 Data Availability
The raw single-cell RNA sequencing data used in this study are publicly available at:

The Mathys et al. dataset is accessible at Synapse:syn18485175.
The Zhou et al. dataset is accessible at Synapse: syn21125841.
The Lau et al. dataset is available at GEO: GSE157827.

Integrated seurat data is available at:
[https://doi.org/10.17632/f5kx9k264n.1.] according to [Li et al. Cell 188, 5516–5534, October 2, 2025]

Required Data Files
To reproduce the analysis, you need the following files placed in the data/ directory:

Processed Seurat object: mzl_seurat.rds.

Metadata table: metadata.csv (included in the sumpplementry table 1 of the article for reference).

🚀 How to Reproduce the Analysis
Step-by-Step Execution
Ensure you have installed all required R packages (see Installation).

Place the Seurat object (mzl_seurat.rds) in the data/ directory.

Open R or RStudio and set your working directory to the root of this repository.

Execute the scripts in numerical order:

r
# Example
source("scripts/01_quality_control.R")
source("scripts/02_cell_type_annotation.R")
source("scripts/03_deg_analysis.R")
# ... continue with remaining scripts
Note: Some scripts may take a considerable amount of time to run (e.g., GSVA, CellChat). We recommend running them on a machine with at least 16 GB of RAM.

Expected Outputs
After successful execution, the following key files will be generated in the output/ directory:

all_de_genes_results.rds / significant_degs.csv: Full DE results.

DEG_heatmap.png: Heatmap of differentially expressed genes.

GSVA_Results/: Folder containing GSVA scores and effect size plots.

Microglia_Functional_Analysis.png: Functional scores for microglia.

synaptic_volcano_plot.png: Volcano plot of synaptic DEGs.

CellChat_*: Cell-cell communication figures and data.

...and many more figures as presented in the paper.

🔍 Reproducibility
For full reproducibility, we have included:

session_info.txt: Contains the output of sessionInfo() from the R environment used to produce the final results. This lists every package and its exact version.

Set seeds: Random seeds are set within each script to ensure consistent results for stochastic processes (e.g., clustering, UMAP).

📄 License
The code in this repository is released under the MIT License. See the LICENSE file for full terms.


📧 Contact
For questions regarding the code or the analysis, please contact the first author of this article

Email: [2310319@mail.nankai.edu.cn]

For questions about the original data, please refer to source GEO pages or contact authors of source papers.
