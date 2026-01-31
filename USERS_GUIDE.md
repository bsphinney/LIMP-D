# LIMP-D User Guide

Welcome to LIMP-D, the Limpa Proteomics Dashboard. This guide will walk you through the steps to analyze your DIA-NN proteomics data, from initial upload to advanced AI-powered insights.

---

## 1. Getting Started

### Installation
Before launching the app, ensure you have installed all the required R packages. You can do this by running the installation command found in the `README.md` file in your R console.

### Launching the App
1.  Open the `LIMP-D.R` file in RStudio.
2.  Click the **"Run App"** button that appears at the top of the script editor.

This will launch the LIMP-D dashboard in a new window or in your web browser.

---

## 2. Step-by-Step Workflow

The app's sidebar on the left guides you through the main analytical steps.

### Step 2.1: Upload Data
1.  **DIA-NN Report:** Click the "Browse..." button and select your DIA-NN report file. The app requires the **`.parquet`** format.
2.  **Q-Value Cutoff:** Set the desired FDR (False Discovery Rate) cutoff for precursor-level data. The default is 0.01 (1%).

### Step 2.2: Assign Groups
1.  After your data uploads, a table will appear for assigning experimental groups. You can also open this table by clicking the **"Assign Groups"** button.
2.  For each file (`File.Name`), type the name of its corresponding experimental group in the `Group` column (e.g., "Control", "Treatment").
3.  **Auto-Guess:** Use the **"🪄 Auto-Guess"** button to automatically assign groups based on common keywords (e.g., 'control', 'treat') in the filenames.
4.  Click **"Save & Close"**. This step is crucial for performing differential expression analysis.

### Step 2.3: Run the Pipeline
*   Click the blue **"Run Pipeline"** button. This performs normalization, statistical modeling, and differential expression analysis. The status message will update to "✅ Complete!" when finished.

### Step 2.4: Explore Results
Once the pipeline is complete, you can explore your data using the main tabs.

#### Data Overview Pane
*   **Summary:** View key statistics about your dataset, including the total number of files and the calculated dynamic range of protein signals in orders of magnitude.
*   **Signal Distribution Plot:** This plot shows the signal intensity for every protein in your dataset.
    *   Click **"Color by DE Status"** to color the points based on whether they are up-regulated, down-regulated, or not significant in the currently selected comparison.
    *   If you select proteins via the AI chat, they will be automatically highlighted and labeled on this plot.
*   **Group QC Summary:** A table showing the average Precursors, MS1 Signal, and identified Proteins for each experimental group.

#### QC Panes (Trends & Plots)
*   Use these tabs to assess the quality of your data. View trends across runs, check for outliers with the MDS plot, and examine the distribution of QC metrics for each group with violin plots.

#### DE Dashboard Pane
*   **Volcano Plot:** Interactively explore differential expression results. Click on a point or draw a box to select proteins.
*   **Results Table:** This table automatically filters to show only the proteins you've selected in the Volcano Plot.
*   **Heatmap:** Visualizes the expression patterns of your selected proteins across all samples.
*   **Buttons:** Use the buttons at the top of the Results Table to generate violin plots of selected proteins or export the full DE results to a CSV file.

#### Gene Set Enrichment Analysis (GSEA) Pane
1.  Ensure you have selected a comparison in the sidebar.
2.  Click **"Run GSEA"**. The app automatically detects the organism (e.g., Human, Mouse, Rat) from your protein IDs and runs the analysis.
3.  Explore the results in the tabs:
    *   **Dot Plot:** An overview of the most significant Gene Ontology (GO) terms.
    *   **Enrichment Map:** A network showing how enriched GO terms are related to each other.
    *   **Ridgeplot:** Shows the distribution of your significant genes within the top enriched GO terms.
    *   **Results Table:** A searchable table of all enriched GO terms.

#### Data Chat Pane (Gemini AI)
This panel allows you to ask questions about your data in plain English.
1.  **API Key:** Enter your Google AI Gemini API key.
2.  **Ask a Question:** Type a question in the text box and click "Send". The AI has access to your QC data and the Top 800 DE proteins.
3.  **Auto-Analyze:** Click the **"🤖 Auto-Analyze"** button to send a pre-defined prompt that asks the AI to provide a concise summary of your dataset's technical quality and key biological findings.
4.  **Bidirectional Interaction:**
    *   When the AI mentions specific proteins in its response (e.g., `[[SELECT: P12345]]`), the Volcano Plot and the Signal Distribution plot will automatically update to highlight and label them.
    *   If you select proteins on the Volcano Plot, the AI will be notified to focus its analysis on your selection for your next question.

---

## 3. FAQ
*   **Why is the app crashing on startup?**
    *   Please ensure you have run the installation command in the `README.md` to install all required packages. An older version of R may also cause issues if it does not support certain syntax used in the app.
*   **Why did the GSEA fail?**
    *   GSEA requires an internet connection to download annotation packages. Ensure you are online. For automatic organism detection to work best, your protein IDs should be in a standard format that includes a UniProt organism suffix (e.g., `P12345_HUMAN`).
