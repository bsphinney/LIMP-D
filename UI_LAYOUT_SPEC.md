# DE-LIMP UI Layout Improvement Spec

## Problem Statement

The current UI layout uses fixed pixel heights (400-600px) for plots, which causes usability issues:
- On typical laptops (13-15"), users must scroll extensively to see plots
- Multiple stacked cards create 1500-2000px of vertical content per tab
- "View Fullscreen" becomes the default way to use the app, not an enhancement
- Poor use of available screen real estate

## Goals

1. Make plots readable without requiring fullscreen on standard monitors
2. Eliminate excessive vertical scrolling within tabs
3. Responsive design that works from 13" laptops to 27" monitors
4. Maintain ability to compare related information when screen allows

---

## Implementation Strategy

Use a **hybrid approach**:
1. **Tabbed sub-panels** for tabs with multiple plots (QC Plots, GSEA)
2. **Viewport-relative heights** (`vh` units) instead of fixed pixels
3. **Responsive grid** for side-by-side layouts on large screens
4. **Accordion panels** for secondary/optional content

---

## Tab-by-Tab Specifications

### 1. QC Plots Tab (MAJOR REFACTOR)

**Current Structure:**
```
QC Plots
├── Card: DPC Fit (400px)
├── Card: MDS Plot (400px)
├── Card: Pipeline Diagnostic (450px + guidance banners)
└── Card: Group QC Violin (400px)
Total: ~1650px vertical, always requires scrolling
```

**New Structure:**
```
QC Plots
└── navset_card_tab (fills available height)
    ├── nav_panel: "Normalization Diagnostic" (DEFAULT/FIRST)
    │   └── Full diagnostic with guidance banners
    │   └── Plot height: calc(100vh - 300px) or ~65vh
    ├── nav_panel: "DPC Fit"
    │   └── Plot height: 70vh
    ├── nav_panel: "MDS Plot"  
    │   └── Plot height: 70vh
    └── nav_panel: "Group Distribution"
        └── Metric selector + violin plot
        └── Plot height: 65vh
```

**Rationale:**
- Normalization Diagnostic is the most important QC view — make it default
- Each plot gets full attention without competing for space
- Sub-tabs allow quick switching without scrolling
- Keep fullscreen buttons as "extra large" option

**Code Pattern:**
```r
nav_panel("QC Plots", icon = icon("chart-line"),
  navset_card_tab(
    id = "qc_subtabs",
    
    nav_panel("Normalization Diagnostic",
      # Guidance banners at top (these are dynamic, keep them)
      uiOutput("norm_diag_guidance"),
      div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
        div(
          uiOutput("diann_norm_status_badge", inline = TRUE),
          radioButtons("norm_diag_type", NULL,
            choices = c("Box Plots" = "boxplot", "Density Overlay" = "density"),
            inline = TRUE
          )
        ),
        actionButton("fullscreen_norm_diag", "🔍 Fullscreen", class = "btn-outline-secondary btn-sm")
      ),
      plotlyOutput("norm_diagnostic_plot", height = "calc(100vh - 350px)"),
      # Expandable help section at bottom
      tags$details(...)
    ),
    
    nav_panel("DPC Fit",
      div(style = "text-align: right; margin-bottom: 10px;",
        actionButton("fullscreen_dpc", "🔍 Fullscreen", class = "btn-outline-secondary btn-sm")
      ),
      plotOutput("dpc_plot", height = "70vh")
    ),
    
    nav_panel("MDS Plot",
      div(style = "text-align: right; margin-bottom: 10px;",
        actionButton("fullscreen_mds", "🔍 Fullscreen", class = "btn-outline-secondary btn-sm")
      ),
      plotOutput("mds_plot", height = "70vh")
    ),
    
    nav_panel("Group Distribution",
      div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
        selectInput("qc_violin_metric", "Metric:", 
          choices = c("Precursors", "Proteins", "MS1_Signal"), 
          width = "200px"
        ),
        actionButton("fullscreen_qc_violin", "🔍 Fullscreen", class = "btn-outline-secondary btn-sm")
      ),
      plotlyOutput("qc_group_violin", height = "calc(100vh - 300px)")
    )
  )
)
```

---

### 2. QC Trends Tab (MINOR REFACTOR)

**Current Structure:**
```
QC Trends
├── Card: Trend Analysis (500px plotly)
└── Card: Stats Table (DT)
Side-by-side layout works but plot is cramped
```

**New Structure:**
```
QC Trends
└── Single card with internal layout
    ├── Header row: metric selector + sort order + fullscreen button
    ├── Plot: 60vh height
    └── Collapsible: Stats Table (accordion, collapsed by default)
```

**Code Pattern:**
```r
nav_panel("QC Trends", icon = icon("chart-bar"),
  card(
    card_header(
      div(style = "display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 10px;",
        div(style = "display: flex; gap: 15px; align-items: center;",
          selectInput("qc_metric_select", "Metric:", 
            choices = c("Precursors", "Proteins", "MS1_Signal"),
            width = "150px"
          ),
          radioButtons("qc_sort_order", "Order:", 
            choices = c("Run Order", "Group"), 
            inline = TRUE
          )
        ),
        actionButton("fullscreen_trend", "🔍 Fullscreen", class = "btn-outline-secondary btn-sm")
      )
    ),
    card_body(
      plotlyOutput("qc_trend_plot", height = "calc(100vh - 320px)")
    )
  ),
  # Stats table as collapsible accordion below
  accordion(
    id = "qc_trends_accordion",
    accordion_panel(
      "View Raw Stats Table",
      DTOutput("r_qc_table"),
      icon = icon("table")
    )
  )
)
```

---

### 3. DE Dashboard Tab (MODERATE REFACTOR)

**Current Structure:**
```
DE Dashboard
├── Comparison banner (fixed, good)
├── layout_columns(6, 6)
│   ├── Card: Results Table
│   └── Card: Volcano Plot (600px)
└── Card: Heatmap (400px)
```

**New Structure:**
```
DE Dashboard
├── Comparison banner (keep as-is)
├── Responsive two-column layout
│   ├── Left: Results Table (scrollable within card)
│   └── Right: Volcano Plot (50vh minimum, scales up)
└── Accordion: Heatmap (collapsed by default, expands inline)
```

**Key Changes:**
- Volcano plot uses viewport height, not fixed pixels
- Heatmap moves to accordion (it's secondary to volcano)
- Table gets internal scroll, doesn't push layout

**Code Pattern:**
```r
nav_panel("DE Dashboard", icon = icon("table-columns"),
  # Comparison banner (unchanged)
  div(style = "background-color: #007bff; color: white; padding: 12px 15px; border-radius: 8px; margin-bottom: 15px; display: flex; align-items: center; justify-content: center; gap: 10px;",
    icon("microscope"),
    span("Viewing Comparison:", style = "font-weight: 500;"),
    uiOutput("current_comparison_display", inline = TRUE, 
      style = "color: #ffe066; font-weight: bold;")
  ),
  
  # Main content - responsive grid
  div(class = "de-dashboard-grid",
    # Results table card
    card(
      card_header(
        div(style = "display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 8px;",
          span("Results Table"),
          div(
            actionButton("generate_ai_summary", "🤖 AI Summary", class = "btn-info btn-sm"),
            actionButton("clear_plot_selection", "Reset", class = "btn-outline-warning btn-sm"),
            actionButton("show_violin", "📊 Violin", class = "btn-outline-primary btn-sm"),
            downloadButton("download_result_csv", "💾 Export", class = "btn-outline-success btn-sm")
          )
        )
      ),
      card_body(
        style = "overflow-y: auto; max-height: calc(100vh - 350px);",
        DTOutput("de_table")
      )
    ),
    
    # Volcano plot card
    card(
      card_header(
        div(style = "display: flex; justify-content: space-between; align-items: center;",
          span("Volcano Plot"),
          actionButton("fullscreen_volcano", "🔍 Fullscreen", class = "btn-outline-secondary btn-sm")
        )
      ),
      card_body(
        plotlyOutput("volcano_plot_interactive", height = "calc(100vh - 350px)")
      )
    )
  ),
  
  # Heatmap as accordion (secondary visualization)
  accordion(
    id = "de_heatmap_accordion",
    open = FALSE,  # Collapsed by default
    accordion_panel(
      "Heatmap of Selected/Top Proteins",
      div(style = "text-align: right; margin-bottom: 10px;",
        actionButton("fullscreen_heatmap", "🔍 Fullscreen", class = "btn-outline-secondary btn-sm")
      ),
      plotOutput("heatmap_plot", height = "450px"),
      icon = icon("grip")
    )
  )
)
```

**Required CSS (add to tags$head):**
```css
.de-dashboard-grid {
  display: grid;
  grid-template-columns: 1fr 1fr;
  gap: 1rem;
  margin-bottom: 1rem;
}

@media (max-width: 1200px) {
  .de-dashboard-grid {
    grid-template-columns: 1fr;
  }
}
```

---

### 4. Data Overview Tab (MINOR REFACTOR)

**Current Structure:**
```
Data Overview
├── Button row
├── Card: Signal Distribution (500px)
└── Card: Group QC Summary Table
```

**New Structure:**
```
Data Overview
├── Compact button row
├── Card: Signal Distribution (60vh, responsive)
└── Accordion: Group QC Summary (collapsed)
```

**Code Pattern:**
```r
nav_panel("Data Overview", icon = icon("database"),
  div(style = "display: flex; gap: 10px; margin-bottom: 15px;",
    actionButton("show_summary_modal", "📋 Full Summary", class = "btn-outline-primary"),
    actionButton("show_grid_view", "📊 Grid View", class = "btn-outline-success")
  ),
  
  card(
    card_header(
      div(style = "display: flex; justify-content: space-between; align-items: center;",
        span("Signal Distribution Across All Protein Groups"),
        div(
          actionButton("color_de", "Color by DE", class = "btn-outline-info btn-sm"),
          actionButton("reset_color", "Reset", class = "btn-outline-secondary btn-sm"),
          actionButton("fullscreen_signal", "🔍 Fullscreen", class = "btn-outline-secondary btn-sm")
        )
      )
    ),
    card_body(
      plotOutput("protein_signal_plot", height = "calc(100vh - 320px)")
    )
  ),
  
  accordion(
    accordion_panel(
      "Group QC Summary",
      DTOutput("group_summary_table"),
      icon = icon("table")
    )
  )
)
```

---

### 5. Gene Set Enrichment Tab (MODERATE REFACTOR)

**Current Structure:**
```
GSEA
├── Card: Run GSEA button + status
└── Card: Results (navset_card_tab with 4 plot tabs)
    Each plot is 500-600px
```

**New Structure:**
```
GSEA
├── Compact control bar (inline button + status)
└── navset_card_tab for results (full height plots)
```

**Code Pattern:**
```r
nav_panel("Gene Set Enrichment", icon = icon("sitemap"),
  # Compact control bar
  card(
    card_body(
      div(style = "display: flex; align-items: center; gap: 15px;",
        actionButton("run_gsea", "▶ Run GSEA", class = "btn-success", icon = icon("play")),
        verbatimTextOutput("gsea_status", placeholder = TRUE) |> 
          tagAppendAttributes(style = "margin: 0; padding: 5px 10px; flex-grow: 1;")
      )
    )
  ),
  
  # Results tabs with full-height plots
  navset_card_tab(
    id = "gsea_results_tabs",
    
    nav_panel("Dot Plot",
      div(style = "text-align: right; margin-bottom: 10px;",
        actionButton("fullscreen_gsea_dot", "🔍 Fullscreen", class = "btn-outline-secondary btn-sm")
      ),
      plotOutput("gsea_dot_plot", height = "calc(100vh - 320px)")
    ),
    
    nav_panel("Enrichment Map",
      div(style = "text-align: right; margin-bottom: 10px;",
        actionButton("fullscreen_gsea_emap", "🔍 Fullscreen", class = "btn-outline-secondary btn-sm")
      ),
      plotOutput("gsea_emapplot", height = "calc(100vh - 320px)")
    ),
    
    nav_panel("Ridgeplot",
      div(style = "text-align: right; margin-bottom: 10px;",
        actionButton("fullscreen_gsea_ridge", "🔍 Fullscreen", class = "btn-outline-secondary btn-sm")
      ),
      plotOutput("gsea_ridgeplot", height = "calc(100vh - 320px)")
    ),
    
    nav_panel("Results Table",
      DTOutput("gsea_results_table")
    )
  )
)
```

---

### 6. Consistent DE Tab (MINOR REFACTOR)

**Current:** Single card with table, fixed height implied by content.

**New:** Use full available height for table, add planned CV histogram as tab.

```r
nav_panel("Consistent DE", icon = icon("check-double"),
  navset_card_tab(
    nav_panel("High-Consistency Proteins",
      p("Ranking by %CV (Coefficient of Variation) to find stable markers.", 
        class = "text-muted small mb-3"),
      DTOutput("consistent_table")
    ),
    nav_panel("CV Distribution",
      # Placeholder for future CV histogram feature
      p("CV histogram by condition (coming soon)", class = "text-muted")
    )
  )
)
```

---

## Global CSS Additions

Add these styles to the existing `tags$head(tags$style(HTML(...)))` block:

```css
/* Viewport-relative plot containers */
.plot-container-vh {
  min-height: 400px;
  max-height: 85vh;
}

/* Accordion styling - more compact */
.accordion {
  margin-top: 1rem;
}

.accordion-button {
  padding: 0.75rem 1rem;
  font-size: 0.95rem;
}

/* Card body with internal scroll for tables */
.card-body-scroll {
  overflow-y: auto;
  max-height: calc(100vh - 350px);
}

/* Responsive grid for DE Dashboard */
.de-dashboard-grid {
  display: grid;
  grid-template-columns: 1fr 1fr;
  gap: 1rem;
  margin-bottom: 1rem;
}

@media (max-width: 1200px) {
  .de-dashboard-grid {
    grid-template-columns: 1fr;
  }
}

/* Fullscreen button consistent styling */
.btn-fullscreen {
  font-size: 0.85rem;
  padding: 0.25rem 0.5rem;
}

/* Compact inline controls */
.controls-inline {
  display: flex;
  align-items: center;
  gap: 10px;
  flex-wrap: wrap;
}

/* Sub-tab navigation styling */
.nav-tabs .nav-link {
  padding: 0.5rem 1rem;
}
```

---

## Implementation Order

1. **QC Plots Tab** — Biggest impact, converts 4 stacked cards to tabbed interface
2. **DE Dashboard** — High-use tab, responsive grid + accordion heatmap
3. **GSEA Tab** — Similar pattern to QC Plots, straightforward
4. **QC Trends Tab** — Minor, accordion for stats table
5. **Data Overview Tab** — Minor, accordion for summary table
6. **Consistent DE Tab** — Minor, prep for future CV histogram

---

## Testing Checklist

After implementation, test on:
- [ ] 13" laptop (1280x800 or 1440x900)
- [ ] 15" laptop (1920x1080)
- [ ] 24" monitor (1920x1080)
- [ ] 27" monitor (2560x1440)

Verify:
- [ ] No horizontal scrolling on any screen size
- [ ] Plots readable without fullscreen on 13" laptop
- [ ] Fullscreen still works and provides extra-large view
- [ ] Sub-tabs switch without page jump
- [ ] Accordions expand/collapse smoothly
- [ ] Table internal scroll works (DE Dashboard)

---

## Notes for Claude Code

1. **Keep all existing functionality** — This is a UI refactor only, no changes to reactive logic or server functions

2. **Preserve fullscreen modals** — The fullscreen button/modal pattern stays, but becomes "extra large" rather than "usable size"

3. **renderPlot height parameters** — When using `calc()` or `vh` in UI, you may need to remove fixed heights from the `renderPlot(..., height = X)` calls in the server, OR keep them as fallbacks. Test both approaches.

4. **bslib version** — This spec assumes bslib 0.5+. Check that `accordion()`, `accordion_panel()`, and `navset_card_tab()` are available.

5. **Backward compatibility** — The changes are purely presentational. No changes to:
   - Reactive values
   - Data processing logic
   - File upload/download
   - AI chat functionality
   - Session save/load
