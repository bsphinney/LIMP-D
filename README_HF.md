---
title: DE-LIMP Proteomics
emoji: 🧬
colorFrom: blue
colorTo: green
sdk: docker
sdk_version: 4.5.0
app_file: app.R
pinned: true
license: mit
tags:
  - proteomics
  - shiny
  - bioinformatics
  - differential-expression
  - mass-spectrometry
  - limma
  - dia-nn
---

# DE-LIMP: Differential Expression & Limpa Proteomics 🧬

An interactive R Shiny application for differential expression analysis of DIA-NN proteomics data. Built on **Limpa** and **Limma** for robust statistics, with **Google Gemini AI** integration.

## ✨ What's New in v2.2

**Contextual Help System** (February 2026):
- ❓ **15 info modal buttons** — click `?` on any tab for in-context guidance
- Covers QC Plots, DE Dashboard, Data Overview, Consistent DE, QC Trends, GSEA, Methodology, Data Chat

**Improved DE Dashboard**:
- 🌋 Volcano → Table filtering (select proteins in volcano to filter the table)
- 🗺️ MDS Plot legend now visible; heatmap expanded by default
- Cleaner layouts: help content in modals, not cluttering plots

**Previous highlights** (v2.1):
- 📈 XIC Chromatogram Viewer with MS2 Intensity Alignment ([local/HPC only](https://github.com/bsphinney/DE-LIMP))
- 🎯 Four-Way Selector Sync across all comparison dropdowns
- 🔬 P-value Distribution Diagnostic with automated pattern detection
- 🌋 Volcano Plot Annotations with FDR/logFC threshold legend

## 🚀 Features

### 📊 Interactive Analysis
- **Volcano Plots** - Fully interactive (Plotly). Click or box-select to highlight
- **Heatmaps** - Auto-scaled Z-score heatmaps of significant proteins
- **QC Trends** - Monitor run quality with group averages
- **Multi-Protein Violin Plots** - Compare expression distributions

### 🤖 AI-Powered Exploration
- **Chat with Your Data** - Google Gemini integration
- **Bi-Directional Sync** - Select proteins ↔ AI suggestions
- **Auto-Summary** - Generate publication-ready summaries

### 💾 Session Management
- **Save/Load Sessions** - Preserve analysis state (.rds files)
- **Reproducibility Logging** - Export complete R code
- **Example Data** - One-click demo dataset (Affinisep vs Evosep)

### 🎓 Education & Resources
- Embedded proteomics training materials
- UC Davis Proteomics video tutorials
- Methodology citations (limpa, limma, DIA-NN)

## 📖 Quick Start

1. **Load Data**: Upload DIA-NN .parquet file or use "Load Example Data"
2. **Assign Groups**: Use auto-guess or manual assignment
3. **Run Pipeline**: Click "▶ Run Pipeline" for limpa analysis
4. **Explore Results**: Interactive plots, tables, GSEA, AI chat

## 🔬 Methodology

- **Normalization**: Data Point Correspondence (DPC-CN)
- **Quantification**: maxLFQ algorithm
- **Statistics**: limma empirical Bayes moderation
- **FDR Correction**: Benjamini-Hochberg

## 📚 Resources

- **GitHub**: [github.com/bsphinney/DE-LIMP](https://github.com/bsphinney/DE-LIMP)
- **Website**: [bsphinney.github.io/DE-LIMP](https://bsphinney.github.io/DE-LIMP/)
- **YouTube**: [UC Davis Proteomics](https://www.youtube.com/channel/UCpulhf8gl-HVxACyJUEFPRw)
- **Core Facility**: [proteomics.ucdavis.edu](https://proteomics.ucdavis.edu)

## 🛠 System Requirements

- R 4.5+
- Bioconductor 3.22+
- All dependencies auto-install on first run

## 👨‍🔬 Developer

**Brett Phinney** - UC Davis Proteomics Core Facility

---

**Built with ❤️ for the proteomics community**
