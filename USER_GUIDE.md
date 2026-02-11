# 📘 DE-LIMP User Guide

Welcome to **DE-LIMP** (Differential Expression & Limpa Proteomics), your interactive dashboard for analyzing DIA-NN proteomics data. This guide covers the complete workflow, from importing data to discovering biological insights with our integrated AI assistant.

---

## 1. Getting Started

### Prerequisites
* **R & RStudio:** Ensure you have R (version 4.2+) installed.
* **Gemini API Key:** Required for AI Chat features (See below).

### 🔑 1.1 How to Obtain a Free Gemini API Key
To use the "Chat with Data" features, you need a key from Google. It is free for standard use.

1.  Go to **[Google AI Studio](https://aistudio.google.com/)**.
2.  Sign in with your Google Account.
3.  In the top-left corner, click the blue button **"Get API key"**.
4.  Click **"Create API key"**.
    * If asked, select "Create API key in new project".
5.  Copy the long string of text that appears (it starts with `AIza...`).
6.  **Paste this key** into the "Gemini API Key" box in the DE-LIMP sidebar.

### 1.2 Launching the App
1.  Open `DE-LIMP.R` in RStudio.
2.  Click the **"Run App"** button at the top right of the script editor.
3.  The dashboard will launch in your default web browser.

---

## 2. Core Analysis Workflow

Follow the sidebar controls on the left to process your data.

### Step 2.1: Upload Data
* **Input File:** Click **"Browse..."** and select your DIA-NN report file.
    * *Requirement:* The file must be in **`.parquet`** format.
    * *Need Data?* Download our [Example Data](https://github.com/bsphinney/DE-LIMP/releases/tag/v1.0) from the repository.
* **Q-Value Cutoff:** Adjust the slider to set your False Discovery Rate (FDR) threshold (Default: 0.01).

### Step 2.2: Assign Groups
This is the most critical step for statistical analysis.
1.  Click **"Assign Groups"** (or it will auto-open after upload).
2.  **The "Magic" Button:** Click **"🪄 Auto-Guess"**. The app will try to detect groups (e.g., "Control", "Treatment", "WT", "KO") based on your filenames.
3.  **Manual Edit:** If the guess is wrong, click on the cells in the **Group** column and type the correct names.
4.  Click **"Save & Close"**.

### Step 2.3: Run Pipeline
* Click the blue **"Run Pipeline"** button.
* **What happens?** The app uses the `limpa` package to perform DPC normalization and the `limma` package to fit linear models for differential expression.
* Wait for the status to change to **"✅ Complete!"**.

### Step 2.4: Select Comparison
* Use the **"Comparison"** dropdown to select which contrast you want to view (e.g., `Treatment - Control`).

---

## 3. Deep Dive: The Data Overview & Grid View

### 📊 Data Overview
This is your landing page.
* **Signal Plot:** Visualizes the dynamic range of your experiment.
    * Click **"Color by DE Status"** to see which proteins are Up/Down-regulated in your current comparison.
* **Group Summary:** A quick table showing average precursor and protein counts per group.

### 🔬 The Grid View (New!)
Click the green **"Open Grid View"** button to open the deep-dive table.

#### **Key Features:**
1.  **Bi-Directional Filtering:**
    * If you select proteins in the **Volcano Plot** (DE Dashboard) or if the **AI** selects interesting proteins, the Grid View automatically filters to show *only those proteins*.
    * Click **"Show All / Clear Selection"** in the footer to reset the view.
2.  **Compact Headers:**
    * Columns are labeled with **Run Numbers** (1, 2, 3...) to save space.
    * **Hover** your mouse over a number to see the full **File Name** and **Group**.
    * Headers are **color-coded** by Experimental Group (refer to the Legend at the top).
3.  **Heatmap Coloring:** Cell values (Log2 Intensity) are colored Blue (Low) to Red (High) for identifying patterns at a glance.
4.  **UniProt Integration:** Click any **Protein ID** to open its official UniProt page in a new tab.
5.  **Click-to-Plot:** Click any row in the table to instantly open a **Violin Plot** showing that specific protein's expression across all samples.
6.  **Smart Export:** Click **"Export Full Table"** to download the data as a CSV. The export will use the **Full Filenames** in the header (not the Run Numbers) for publication use.

---

## 4. Visualizing Results

### 📉 DE Dashboard
* **Volcano Plot:** Interactive! Click points to select them. Box-select multiple points to analyze a cluster.
    * *Sync:* Selecting points here updates the Grid View and the AI context.
* **Violin Plots:** Select a protein in the volcano plot and click the **"📊 Violin Plot"** button to see its expression profile.
* **Heatmap:** Automatically scales and clusters the top 50 significant proteins (or your specific selection).

### 📐 QC Trends & Plots
* **Trends:** Check "Precursors per Run" to spot batch effects or instrument drift.
* **MDS Plot:** A multidimensional scaling plot to visualize how samples cluster. (Good samples should cluster by Group).

### 🧬 Gene Set Enrichment (GSEA)
1.  Click **"Run GSEA"**.
2.  The app automatically detects if your data is Human, Mouse, or Rat based on UniProt suffixes.
3.  View results as Dot Plots, Enrichment Maps (networks), or Ridgeplots to understand biological pathways.

---

## 5. 🤖 AI Chat (Gemini Integration)

DE-LIMP features a context-aware AI assistant.

### Setup
1.  Paste your **Gemini API Key** in the sidebar.
2.  (Optional) Change the Model Name if you want to use a specific version (Default: `gemini-3-flash-preview`).

### "Chat with Your Data"
You aren't just chatting with a bot; you are chatting with **your specific dataset**.
* **Auto-Analyze:** Click this button to generate a comprehensive report summarizing QC quality and the top biological findings.
* **Ask Questions:**
    * *"Which group has the highest variance?"*
    * *"Are there any mitochondrial proteins upregulated?"*
    * *"Generate a figure caption for the volcano plot."*

### Bi-Directional AI Sync
* **User -> AI:** Select points on the Volcano Plot. Then ask: *"What are the functions of these selected proteins?"*. The AI knows exactly which ones you clicked.
* **AI -> User:** If the AI finds interesting proteins (e.g., *"I found several glycolytic enzymes..."*), it will highlight them in your plots automatically.

---

## 6. Troubleshooting

| Issue | Solution |
| :--- | :--- |
| **App crashes on startup** | Run the installation commands in the `README.md` to ensure `limpa`, `limma`, and `shiny` are installed. |
| **GSEA fails** | Ensure you are connected to the internet (it needs to download gene ontologies). |
| **Grid View "Object not found"** | Ensure you have run the pipeline first. The Grid View requires processed data. |
| **AI says "No data"** | Click the "Run Pipeline" button first. The AI needs the statistical results to answer questions. |

---

*Happy analyzing!* 🧬
