# CrackPipe — Proteomics Data Quality Explorer

An educational Streamlit app that teaches proteomics researchers how to spot common data quality problems in mass spectrometry-based proteomics. Inspired by [FragPipe](https://github.com/Nesvilab/FragPipe), but instead of processing data, it shows what happens when things go **terribly wrong**.

All data is **synthetically generated** — no real datasets needed.

> **Not affiliated** with FragPipe or Nesvilab. The name "CrackPipe" is a humorous parody — the logo is a cracked laboratory flask.

## Modules

| # | Module | Severity | What it teaches |
|---|--------|----------|-----------------|
| 1 | Chimeric Spectra | HIGH | Mixed MS2 spectra from co-isolated peptides |
| 2 | Spray Instability | HIGH | TIC dropouts from electrospray problems |
| 3 | Mass Cal Drift | MEDIUM | Systematic mass error over LC-MS runs |
| 4 | Batch Effects | HIGH | PCA separating by batch instead of biology |
| 5 | Missing Values | MEDIUM | MCAR vs MNAR missingness patterns |
| 6 | TMT Compression | MEDIUM | Co-isolation interference compressing fold changes |
| 7 | Contamination | HIGH | Keratin and PEG polymer contamination |
| 8 | Incomplete Digestion | MEDIUM | Missed cleavage rate assessment |
| 9 | LC Carryover | MEDIUM | Ghost peaks from previous runs |
| 10 | DIA Pitfalls | HIGH | Window design, library bias, deconvolution |
| 11 | Missing Proteins | HIGH | FASTA database and PTM search issues |
| 12 | AP-MS Normalization | HIGH | Why standard normalization fails for pulldowns |
| 13 | Western vs MS | MEDIUM | Reconciling antibody and MS disagreements |
| 14 | P-Value Shopping | HIGH | Running multiple tests and cherry-picking |
| 15 | Cherry-Picking Data | HIGH | Selective sample removal to force significance |
| 16 | Multiple Testing | HIGH | Why raw p-values are meaningless at scale |

Plus a **15-question quiz** to test diagnostic skills.

## Quick Start

### Local

```bash
pip install -r requirements.txt
streamlit run app.py
```

### HuggingFace Spaces

1. Create a new Space (SDK: **Streamlit**)
2. Upload `app.py` and `requirements.txt`
3. Done — auto-deploys

### Streamlit Cloud

1. Push repo to GitHub
2. Go to [share.streamlit.io](https://share.streamlit.io)
3. Point at this repo and `app.py`
4. Deploy

## Tech Stack

- Python 3.10+
- Streamlit >= 1.28
- Plotly for interactive charts
- NumPy / Pandas for synthetic data
- No database, no external APIs, single-file app
