# DE-LIMP Project Context for Claude

## Working Preferences
Proactively suggest updates to this file when new patterns, gotchas, or architectural decisions emerge.

## Project Overview
DE-LIMP is a Shiny proteomics data analysis pipeline using the LIMPA R package for differential expression analysis of DIA-NN data.

- **GitHub**: https://github.com/bsphinney/DE-LIMP
- **Hugging Face**: https://huggingface.co/spaces/brettsp/de-limp-proteomics
- **Local URL**: http://localhost:3838

## System Requirements
- **R**: 4.5+ (required for limpa)
- **Bioconductor**: 3.22+ (automatic with R 4.5+)
- **Download R**: https://cloud.r-project.org/

## Key Files

| File | Purpose |
|------|---------|
| `app.R` | Main Shiny app (~2800 lines) - used by HF Spaces |
| `DE-LIMP.R` | Identical copy of app.R - for GitHub releases/local users |
| `Dockerfile` | Docker container for HF Spaces and HPC deployment |
| `UI_LAYOUT_SPEC.md` | Spec for v2.1 responsive UI layout refactoring |
| `NORMALIZATION_DIAGNOSTIC_SPEC.md` | Spec for normalization diagnostic plot feature |
| `HPC_DEPLOYMENT.md` | Guide for HPC cluster deployment with Apptainer/Singularity |
| `USER_GUIDE.md` | End-user documentation |
| `README_GITHUB.md` | **SOURCE** for GitHub README (edit this!) |
| `README_HF.md` | **SOURCE** for HF README (edit this!) |
| `README.md` | Generated - content differs between remotes |
| `Citations/` | Reference PDFs (LIMMA, FDR, DIA-NN papers) |
| `docs/index.html` | GitHub Pages site |
| `.github/workflows/sync-to-hf.yml` | Auto-sync GitHub to HF |

## App Architecture

### Structure (app.R)
```
Lines 1-110:    Auto-installation & setup (R/Bioc version checks, package install)
Lines 111-168:  Server configuration (max upload size, Gemini API)
Lines 169-330:  Helper functions (QC stats, AI chat, organism detection)
Lines 331-830:  UI definition (page_sidebar with 9 main tabs)
Lines 831-2809: Server logic (~77 reactive elements)
Line 2810:      shinyApp(ui, server)
```

### Main Tab Structure
- **Data Overview** - 4 sub-tabs: Signal Distribution, Dataset Summary, Group QC Summary, Expression Grid
- **QC Trends** - 4 sub-tabs: Precursors, Proteins, MS1 Signal, Stats Table
- **QC Plots** - 4 sub-tabs: Normalization Diagnostic, DPC Fit, MDS Plot, Group Distribution
- **DE Dashboard** - Volcano plot, results table, violin plots, comparison banner
- **Consistent DE** - Cross-comparison consistency analysis
- **Reproducibility** - Code Log + Methodology sub-tabs
- **Gene Set Enrichment** - Dot Plot, Enrichment Map, Ridgeplot, Results Table
- **Data Chat** - AI-powered analysis (Google Gemini API)
- **Education** - Learning resources

### Key Reactive Values
- `values$plot_selected_proteins` - Selected proteins from table/volcano (used by all viz)
- `values$fit` - limma fit object from DE analysis
- `values$y_protein` - Protein-level quantification matrix
- `values$repro_log` - Cumulative R code for reproducibility

### LIMPA Pipeline Flow
1. `readDIANN()` - Load DIA-NN parquet file
2. `dpcCN()` - Data Point Correspondence normalization
3. `dpcQuant()` - Peptide to protein quantification
4. `dpcDE()` - Differential expression model
5. `contrasts.fit()` + `eBayes()` - Pairwise comparisons

### Selection System
- Volcano plot click/box-select and table row selection both update `values$plot_selected_proteins`
- Bidirectional sync for highlighting

### AI Chat Feature
- Uses Google Gemini API, uploads top 800 proteins via File API
- Bidirectional: user selects proteins for AI analysis, AI suggests proteins to highlight
- Selection format: `[[SELECT: P12345; P67890]]`

## Development Workflow

### Running Locally
```r
# In VS Code R terminal (recommended):
shiny::runApp('/Users/brettphinney/Documents/claude/DE-LIMP.r', port=3838, launch.browser=TRUE)

# From command line:
Rscript -e "shiny::runApp('/Users/brettphinney/Documents/claude/DE-LIMP.r', port=3838)" &
```

- **DO NOT** use `source()` - it doesn't work properly in VS Code
- Shiny apps don't hot-reload - must restart after every code change
- Stop: `pkill -f "DE-LIMP.r"` or `Ctrl+C`
- Check if running: `lsof -i :3838`

### Keeping app.R and DE-LIMP.R in Sync
These files must be identical. After editing one: `cp DE-LIMP.R app.R` (or vice versa).

## Deployment

### Three Platforms
1. **GitHub** (`origin` remote) - Source code, releases, download DE-LIMP.R
2. **Hugging Face Spaces** (`hf` remote) - Web app via Docker (app.R)
3. **HPC Clusters** - Apptainer/Singularity containers (see HPC_DEPLOYMENT.md)

### GitHub Actions Auto-Sync (Primary Workflow)
Pushing to `origin/main` automatically syncs to HF within 1-2 minutes.
- Workflow: `.github/workflows/sync-to-hf.yml`
- Syncs all files **except** README.md and .github/
- Uses `README_HF.md` as HF's README (with YAML frontmatter)
- **You only need `git push origin main`** - HF updates automatically

```bash
# Standard workflow for code changes:
git add DE-LIMP.R app.R  # be specific with files
git commit -m "Description"
git push origin main      # HF syncs automatically
```

### README Management (CRITICAL)
The two remotes have **permanently different README.md files**. GitHub Actions handles this automatically, but if you ever need to manually update READMEs:

- Edit `README_GITHUB.md` for GitHub, push to `origin` only
- Edit `README_HF.md` for HF, push to `hf` only
- **NEVER** push README.md changes to both remotes
- **NEVER** use `git add .` when README.md is modified

**Recovery if HF breaks** ("Missing configuration in README" error):
```bash
cp README.md README_GITHUB_BACKUP.md
cp README_HF.md README.md
git add README.md && git commit -m "Fix HF: Restore YAML README" && git push hf main
cp README_GITHUB_BACKUP.md README.md
git add README.md && git commit -m "Restore GitHub README" && git push origin main
rm README_GITHUB_BACKUP.md
```

### Adding New R Packages
When adding features that require new packages:
1. Add `library()` call to app.R/DE-LIMP.R
2. Update Dockerfile to install the package (CRAN in Section 2, Bioconductor in Section 3)
3. Consider system dependencies: graphics packages need `libcairo2-dev`, XML/web need `libxml2-dev`
4. Commit app.R, DE-LIMP.R, and Dockerfile together

### Minimize HF Builds
HF Docker builds take 5-10 min (cached) or 30-45 min (Dockerfile changes). Always test locally first, batch changes, and push only when fully tested.

## Key Patterns & Gotchas

### R Shiny Reactivity
- **DT table row indices** refer to the CURRENTLY DISPLAYED (possibly filtered) data
- **Avoid circular reactivity**: Don't filter data based on a reactive value that the filter updates
- **CRITICAL reactive loop**: If `renderDT` depends on a reactive that uses `values$selection`, AND selecting rows updates `values$selection`, you get an infinite loop (selection triggers re-render which resets selection). **Solution**: Make `renderDT` build data independently, not using reactives that depend on selection state.
- Use `isolate()` to break reactive dependencies when needed
- DT tables with `selection = "multiple"` support Ctrl+Click and Shift+Click

### Package Installation & Loading
- Install packages BEFORE any `library()` calls (attempting to update loaded packages causes unload errors)
- Use `update = FALSE` in `BiocManager::install()` and `upgrade = "never"` in `remotes::install_*()`
- Set CRAN mirror before installs: `options(repos = c(CRAN = "https://cloud.r-project.org"))`

### dplyr/tidyverse
- In `summarise()`, use explicit `{}` braces for multi-line if statements
- Always set `.groups = 'drop'` to avoid grouping warnings

### Reproducibility Logging
- Uses cumulative `add_to_log(action_name, code_lines)` helper function
- Logs with timestamps, appends rather than overwrites
- Include context: which button was clicked, which parameters changed

### UI Design Patterns (v2.1)
- All plot heights use viewport-relative units (`vh`, `calc()`) - no fixed pixel heights
- Multi-plot sections use `navset_card_tab` sub-tabs instead of dropdowns or modals
- DE Dashboard uses responsive grid (`.de-dashboard-grid`) that stacks at <1200px
- Fullscreen modals available as enhancement on all plot panels

## Common Issues

| Problem | Solution |
|---------|----------|
| `source()` doesn't start app in VS Code | Use `shiny::runApp()` instead |
| App doesn't reflect code changes | Restart app (no hot-reload) |
| "Please select a CRAN mirror" hangs | Already fixed: CRAN mirror set at line 8 |
| "limpa not found" / R version error | Need R 4.5+ from https://cloud.r-project.org/ |
| "Package cannot be unloaded" during install | Restart R session, packages install before loading |
| Selections disappear after clicking | Reactive loop - table must not depend on selection-derived reactives |
| HF "Missing configuration in README" | Run README recovery script above |

## Current TODO
- [ ] Volcano plot: Add FDR threshold annotation/legend
- [ ] GSEA: Add KEGG/Reactome enrichment; clarify which contrast is used
- [ ] DE Dashboard: Make comparison bar clickable to change contrasts
- [ ] Grid View: Open violin plot on protein click with bar plot toggle
- [ ] Consistent DE: Add CV histogram by condition
- [ ] Data Overview: Show which comparison is used for signal distribution
- [ ] Publication-quality plot exports (SVG/PNG/TIFF with size controls)
- [ ] Sample correlation heatmap (QC Plots tab)
- [ ] Venn diagram of significant proteins across comparisons
- [ ] P-value histogram diagnostic
- [ ] Sample CV distribution plots
- [ ] Protein numbers bar plot per sample
- [ ] Absence/presence table for on/off proteins
