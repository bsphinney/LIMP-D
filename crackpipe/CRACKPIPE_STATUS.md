# CrackPipe Build Status

## What Is This?
CrackPipe is a Streamlit proteomics data quality education app. Full spec provided by the user — see conversation history or the README.md in this folder.

## Target Repository
- **GitHub**: https://github.com/bsphinney/Crackpipe
- **Problem**: The Claude Code proxy doesn't have access to the Crackpipe repo yet. Files are temporarily stored here in `DE-LIMP/crackpipe/` on branch `claude/create-crackpipe-repo-YyAIu`.
- **Action needed**: Start a new Claude Code session with `bsphinney/Crackpipe` added to the GitHub App, then move these files there.

## Files
| File | Status |
|------|--------|
| `app.py` | ~1130 lines, runnable, modules 1-8 complete |
| `requirements.txt` | Complete |
| `README.md` | Complete |

## What's Done (Modules 1-8)

### Core Infrastructure
- Inline SVG logo (cracked Erlenmeyer flask, two sizes)
- FragPipe-inspired CSS (dark sidebar, colored cards, severity badges, status bar)
- Google Fonts (IBM Plex Mono + IBM Plex Sans)
- Helper functions: `header_bar()`, `severity_badge()`, `card()`, `styled_fig()`
- Sidebar navigation with radio buttons for all 16 modules + quiz
- Plotly chart styling with `plotly_white` template

### Completed Modules
1. **Home** — Welcome card, 2×3 category grid, status bar
2. **Chimeric Spectra** (HIGH) — Clean vs chimeric MS2 stick spectra with b/y ion labels
3. **Spray Instability** (HIGH) — Stable vs unstable TIC traces with dropouts/crash
4. **Mass Cal Drift** (MEDIUM) — Mass error scatter plots with drift trend line
5. **Batch Effects** (HIGH) — PCA plots showing biology vs batch driving PC1
6. **Missing Values** (MEDIUM) — MCAR vs MNAR heatmaps (80 proteins × 10 samples)
7. **TMT Compression** (MEDIUM) — Label-free vs TMT measured-vs-true FC scatter
8. **Contamination** (HIGH) — Sub-tabs for Keratin (bar chart), PEG (44 Da spectrum), Clean
9. **Incomplete Digestion** (MEDIUM) — Donut charts comparing MC rates (82% vs 38%)
10. **LC Carryover** (MEDIUM) — Two-panel TIC: sample run vs blank with ghost peaks

All modules have:
- Header bar with logo
- Severity badge
- Problem description card (colored border)
- Side-by-side Good ✅ vs Bad ❌ comparison plots
- Educational expander with detection/fix guidance
- Synthetic data with fixed random seeds

## What's Left (Modules 9-16 + Quiz)

### Modules Needing Full Implementation (currently stubs)
11. **DIA Pitfalls** (HIGH) — 3 sub-tabs: Window Design, Library Bias, Deconvolution Failure
12. **Missing Proteins** (HIGH) — 2 sub-tabs: FASTA Problems, Missing PTM Search (collagen/HyPro)
13. **AP-MS Normalization** (HIGH) — 3-panel: Raw counts, Median-normalized (wrong), SAINTexpress (right)
14. **Western vs MS** (MEDIUM) — 3 sub-tabs: Sensitivity Gap, Antibody Specificity, When MS Wins
15. **P-Value Shopping** (HIGH) — 2 sub-tabs: The Problem (lollipop chart), Interactive Demo (simulation)
16. **Cherry-Picking Data** (HIGH) — Two volcano plots: all data vs after removing "outlier" replicates
17. **Multiple Testing** (HIGH) — 3-panel: raw p volcano, adjusted p volcano, p-value histogram
18. **Quiz** — 15 multiple-choice questions with scoring, explanations, st.balloons() for perfect score

### Per-Module Details from Spec
Each stub module (`_stub()` function in app.py) needs to be replaced with a full implementation following the same pattern as modules 1-8. The full spec for each module (plot types, data generation formulas, educational content) was provided in the original user request. Key details:

- **DIA Pitfalls**: stepped bar chart for windows, grouped bars for library comparison, XIC Gaussian peaks
- **Missing Proteins**: Venn/grouped bars for FASTA overlap, bar chart for collagen peptides ± HyPro
- **AP-MS**: dot/bar plot of spectral counts → median-normalized → SAINT scores
- **Western vs MS**: waterfall abundance plot with LOD lines, 2×2 info cards, isoform grouped bars
- **P-Value Shopping**: horizontal lollipop of 6 test p-values, simulation bar chart (false positive rate vs # tests)
- **Cherry-Picking**: two volcano plots with strip chart inset showing removed replicates
- **Multiple Testing**: 4000 proteins, 3800 null + 200 true positives, BH adjustment, p-value histogram
- **Quiz**: 15 questions with st.radio(index=None), green/red cards, session_state score tracking

## How to Run Locally
```bash
cd crackpipe
pip install -r requirements.txt
streamlit run app.py
```

## Tech Stack
- Python 3.10+, Streamlit ≥ 1.28, Plotly ≥ 5.15, NumPy, Pandas
- Single-file app (everything in app.py)
- No database, no external APIs
- All data synthetically generated with fixed seeds
