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
git clone https://github.com/yourusername/your-repo-name.git
cd your-repo-name
