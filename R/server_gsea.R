server_gsea <- function(input, output, session, values, add_to_log) {

  # --- Helper: ontology label ---
  get_ont_label <- function(ont) {
    switch(ont,
      "BP" = "GO Biological Process",
      "MF" = "GO Molecular Function",
      "CC" = "GO Cellular Component",
      "KEGG" = "KEGG Pathways",
      "Enrichment"
    )
  }

  # --- GSEA Info Modal ---
  observeEvent(input$gsea_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " What is Gene Set Enrichment Analysis?"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("What is GSEA?"),
        p("Gene Set Enrichment Analysis tests whether predefined groups of genes (e.g., biological pathways) ",
          "show coordinated changes in your experiment. Instead of looking at individual proteins, ",
          "it asks: 'Are proteins in this pathway collectively going up or down?'"),
        tags$h6("Enrichment Databases"),
        tags$ul(
          tags$li(strong("GO Biological Process (BP): "),
                  "What biological programs are affected? (e.g., 'immune response', 'cell cycle')"),
          tags$li(strong("GO Molecular Function (MF): "),
                  "What molecular activities are enriched? (e.g., 'kinase activity', 'DNA binding')"),
          tags$li(strong("GO Cellular Component (CC): "),
                  "Where in the cell are changes concentrated? (e.g., 'mitochondrion', 'nucleus')"),
          tags$li(strong("KEGG Pathways: "),
                  "Which metabolic & signaling pathways are affected? (e.g., 'PI3K-Akt signaling'). ",
                  "Requires internet access.")
        ),
        tags$h6("Caching"),
        p("Results are cached per ontology. After running BP, you can switch to MF/CC/KEGG and run each separately. ",
          "Switching back to a previously computed ontology loads results instantly without re-running."),
        tags$h6("The visualization tabs"),
        tags$ul(
          tags$li(strong("Dot Plot: "), "Each dot = one enriched pathway. Size = number of genes, color = adjusted p-value."),
          tags$li(strong("Enrichment Map: "), "Network showing how enriched pathways relate to each other. Connected pathways share genes."),
          tags$li(strong("Ridgeplot: "), "Density curves showing the fold-change distribution of genes within each enriched pathway."),
          tags$li(strong("Results Table: "), "Full numeric results with statistics for each pathway.")
        ),
        tags$h6("Which comparison is used?"),
        p("GSEA uses the currently selected comparison from the DE Dashboard. ",
          "If you change the contrast after running GSEA, a warning will appear. Re-run to update results.")
      )
    ))
  })

  # --- GSEA Results Table Info Modal ---
  observeEvent(input$gsea_table_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " GSEA Results Table Columns"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("Column definitions"),
        tags$ul(
          tags$li(strong("ID: "), "Gene Ontology accession (e.g., GO:0006955) or KEGG pathway ID"),
          tags$li(strong("Description: "), "Human-readable name of the biological process or pathway"),
          tags$li(strong("setSize: "), "Number of genes in this pathway that were measured in your data"),
          tags$li(strong("enrichmentScore: "), "Running enrichment score \u2014 reflects degree of overrepresentation at the top or bottom of the ranked gene list"),
          tags$li(strong("NES: "), "Normalized Enrichment Score \u2014 accounts for pathway size. Positive = enriched in up-regulated genes, negative = enriched in down-regulated genes"),
          tags$li(strong("pvalue: "), "Raw p-value from the enrichment test"),
          tags$li(strong("p.adjust: "), "FDR-adjusted p-value (Benjamini-Hochberg). Below 0.05 = statistically significant"),
          tags$li(strong("rank: "), "Position in the ranked gene list where the enrichment score peaks")
        )
      )
    ))
  })

  # --- Contrast Indicator ---
  output$gsea_contrast_indicator <- renderUI({
    if (!is.null(input$contrast_selector) && !is.null(values$fit)) {
      stale_warning <- NULL
      if (!is.null(values$gsea_last_contrast) &&
          values$gsea_last_contrast != input$contrast_selector &&
          length(values$gsea_results_cache) > 0) {
        stale_warning <- tags$div(
          class = "alert alert-warning py-2 px-3 mb-2",
          style = "font-size: 0.85em;",
          icon("exclamation-triangle"),
          sprintf(" Cached results are from contrast '%s'. Re-run to update.",
                  values$gsea_last_contrast)
        )
      }
      tagList(
        tags$div(
          class = "alert alert-info py-2 px-3 mb-2",
          style = "font-size: 0.9em;",
          icon("exchange-alt"),
          paste(" Active contrast:", input$contrast_selector)
        ),
        stale_warning
      )
    } else {
      tags$div(
        class = "alert alert-warning py-2 px-3 mb-2",
        style = "font-size: 0.9em;",
        icon("exclamation-triangle"),
        " Run differential expression analysis first to enable enrichment."
      )
    }
  })

  # --- Load Cached Results on Ontology Switch ---
  observeEvent(input$gsea_ontology, {
    ont <- input$gsea_ontology
    cached <- values$gsea_results_cache[[ont]]
    ont_label <- get_ont_label(ont)
    if (!is.null(cached)) {
      values$gsea_results <- cached
      n_results <- nrow(as.data.frame(cached))
      output$gsea_status <- renderText(
        sprintf("Showing cached %s results: %d terms (contrast: %s)",
                ont_label, n_results, values$gsea_last_contrast)
      )
    } else {
      values$gsea_results <- NULL
      output$gsea_status <- renderText(
        sprintf("No cached results for %s. Click 'Run GSEA' to compute.", ont_label)
      )
    }
  })

  # --- Run GSEA (supports BP/MF/CC/KEGG) ---
  observeEvent(input$run_gsea, {
    req(values$fit, input$contrast_selector)
    ont <- input$gsea_ontology
    ont_label <- get_ont_label(ont)

    output$gsea_status <- renderText(
      sprintf("Running %s enrichment... This may take a few minutes.", ont_label)
    )

    withProgress(message = paste("Running", ont_label), {
      tryCatch({
        # --- Step 1: Prepare gene list ---
        incProgress(0.1, detail = "Preparing gene list...")
        de_results <- topTable(values$fit, coef = input$contrast_selector, number = Inf)
        org_db_name <- detect_organism_db(rownames(de_results))
        if (!requireNamespace(org_db_name, quietly = TRUE)) {
          BiocManager::install(org_db_name, ask = FALSE)
        }
        library(org_db_name, character.only = TRUE)

        incProgress(0.3, detail = "Converting Gene IDs...")
        protein_ids <- str_extract(rownames(de_results), "[A-Z0-9]+")
        id_map <- bitr(protein_ids, fromType = "UNIPROT", toType = "ENTREZID",
                       OrgDb = get(org_db_name))
        de_results_mapped <- de_results[rownames(de_results) %in% id_map$UNIPROT, ]
        de_results_mapped$ENTREZID <- id_map$ENTREZID[
          match(str_extract(rownames(de_results_mapped), "[A-Z0-9]+"), id_map$UNIPROT)
        ]
        gene_list <- de_results_mapped$logFC
        names(gene_list) <- de_results_mapped$ENTREZID
        gene_list <- sort(gene_list, decreasing = TRUE)
        gene_list <- gene_list[!duplicated(names(gene_list))]

        # --- Step 2: Run appropriate enrichment ---
        incProgress(0.6, detail = paste("Running", ont_label, "enrichment..."))

        if (ont %in% c("BP", "MF", "CC")) {
          gsea_res <- gseGO(
            geneList = gene_list,
            OrgDb = get(org_db_name),
            keyType = "ENTREZID",
            ont = ont,
            minGSSize = 10,
            maxGSSize = 500,
            pvalueCutoff = 0.05,
            verbose = FALSE
          )
        } else if (ont == "KEGG") {
          kegg_org <- switch(org_db_name,
            "org.Hs.eg.db" = "hsa",
            "org.Mm.eg.db" = "mmu",
            "org.Rn.eg.db" = "rno",
            "org.Dm.eg.db" = "dme",
            "org.Sc.sgd.db" = "sce",
            "org.Ce.eg.db" = "cel",
            "org.Dr.eg.db" = "dre",
            "org.Bt.eg.db" = "bta",
            "org.Cf.eg.db" = "cfa",
            "org.Gg.eg.db" = "gga",
            "org.Ss.eg.db" = "ssc",
            "hsa"  # fallback to human
          )
          if (kegg_org == "hsa" && org_db_name != "org.Hs.eg.db") {
            showNotification(
              paste("KEGG organism code not mapped for", org_db_name,
                    "\u2014 using human (hsa). Results may be incorrect."),
              type = "warning", duration = 10
            )
          }
          gsea_res <- gseKEGG(
            geneList = gene_list,
            organism = kegg_org,
            keyType = "ncbi-geneid",
            minGSSize = 10,
            maxGSSize = 500,
            pvalueCutoff = 0.05,
            verbose = FALSE
          )
        }

        # --- Step 3: Cache and display ---
        values$gsea_results_cache[[ont]] <- gsea_res
        values$gsea_results <- gsea_res
        values$gsea_last_contrast <- input$contrast_selector
        values$gsea_last_org_db <- org_db_name
        incProgress(1, detail = "Complete.")

        n_results <- nrow(as.data.frame(gsea_res))
        output$gsea_status <- renderText(
          sprintf("%s complete: %d significant terms (contrast: %s)",
                  ont_label, n_results, input$contrast_selector)
        )

        # --- Step 4: Log for reproducibility ---
        if (ont %in% c("BP", "MF", "CC")) {
          gsea_code <- c(
            "# ============================================================",
            "# GENE SET ENRICHMENT ANALYSIS",
            "# ============================================================",
            sprintf("# Contrast: %s", input$contrast_selector),
            sprintf("# Database: %s (%s)", ont, ont_label),
            sprintf("# Organism: %s", org_db_name),
            sprintf("# Date: %s", Sys.time()),
            "# Parameters: minGSSize=10, maxGSSize=500, pvalueCutoff=0.05",
            "",
            "library(clusterProfiler)",
            sprintf("library(%s)", org_db_name),
            "",
            sprintf("de_results <- topTable(fit, coef='%s', number=Inf)", input$contrast_selector),
            "protein_ids <- str_extract(rownames(de_results), '[A-Z0-9]+')",
            sprintf("id_map <- bitr(protein_ids, fromType='UNIPROT', toType='ENTREZID', OrgDb=%s)", org_db_name),
            "de_results_mapped <- de_results[rownames(de_results) %in% id_map$UNIPROT, ]",
            "de_results_mapped$ENTREZID <- id_map$ENTREZID[match(str_extract(rownames(de_results_mapped), '[A-Z0-9]+'), id_map$UNIPROT)]",
            "gene_list <- de_results_mapped$logFC",
            "names(gene_list) <- de_results_mapped$ENTREZID",
            "gene_list <- sort(gene_list, decreasing=TRUE)",
            "gene_list <- gene_list[!duplicated(names(gene_list))]",
            "",
            sprintf("gsea_res <- gseGO(geneList=gene_list, OrgDb=%s, keyType='ENTREZID',", org_db_name),
            sprintf("                  ont='%s', minGSSize=10, maxGSSize=500, pvalueCutoff=0.05)", ont)
          )
        } else if (ont == "KEGG") {
          gsea_code <- c(
            "# ============================================================",
            "# GENE SET ENRICHMENT ANALYSIS",
            "# ============================================================",
            sprintf("# Contrast: %s", input$contrast_selector),
            "# Database: KEGG Pathways",
            sprintf("# Organism: %s (%s)", org_db_name, kegg_org),
            sprintf("# Date: %s", Sys.time()),
            "# Parameters: minGSSize=10, maxGSSize=500, pvalueCutoff=0.05",
            "",
            "library(clusterProfiler)",
            "",
            sprintf("de_results <- topTable(fit, coef='%s', number=Inf)", input$contrast_selector),
            "protein_ids <- str_extract(rownames(de_results), '[A-Z0-9]+')",
            sprintf("id_map <- bitr(protein_ids, fromType='UNIPROT', toType='ENTREZID', OrgDb=%s)", org_db_name),
            "de_results_mapped <- de_results[rownames(de_results) %in% id_map$UNIPROT, ]",
            "de_results_mapped$ENTREZID <- id_map$ENTREZID[match(str_extract(rownames(de_results_mapped), '[A-Z0-9]+'), id_map$UNIPROT)]",
            "gene_list <- de_results_mapped$logFC",
            "names(gene_list) <- de_results_mapped$ENTREZID",
            "gene_list <- sort(gene_list, decreasing=TRUE)",
            "gene_list <- gene_list[!duplicated(names(gene_list))]",
            "",
            sprintf("gsea_res <- gseKEGG(geneList=gene_list, organism='%s',", kegg_org),
            "                    keyType='ncbi-geneid', minGSSize=10, maxGSSize=500, pvalueCutoff=0.05)"
          )
        }
        add_to_log(paste("GSEA:", ont_label), gsea_code)

      }, error = function(e) {
        output$gsea_status <- renderText(paste("GSEA Error:", e$message))
      })
    })
  })

  # --- Plot Renderers with Dynamic Titles ---
  output$gsea_dot_plot <- renderPlot({
    req(values$gsea_results)
    ont_label <- get_ont_label(input$gsea_ontology)
    if (nrow(values$gsea_results) > 0) {
      dotplot(values$gsea_results, showCategory = 20) + ggtitle(paste("GSEA:", ont_label))
    } else {
      plot(NULL, xlim = c(0, 1), ylim = c(0, 1),
           main = "No significant enrichment found.", xaxt = 'n', yaxt = 'n')
    }
  })

  output$gsea_emapplot <- renderPlot({
    req(values$gsea_results)
    if (nrow(values$gsea_results) > 0) emapplot(pairwise_termsim(values$gsea_results), showCategory = 20)
  })

  output$gsea_ridgeplot <- renderPlot({
    req(values$gsea_results)
    if (nrow(values$gsea_results) > 0) ridgeplot(values$gsea_results)
  })

  output$gsea_results_table <- renderDT({
    req(values$gsea_results)
    datatable(as.data.frame(values$gsea_results), options = list(pageLength = 10, scrollX = TRUE))
  })

  # --- Fullscreen Plot Renderers ---
  observeEvent(input$fullscreen_gsea_dot, {
    showModal(modalDialog(
      title = "GSEA Dot Plot - Fullscreen View",
      plotOutput("gsea_dot_plot_fs", height = "700px"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close")
    ))
  })
  output$gsea_dot_plot_fs <- renderPlot({
    req(values$gsea_results)
    ont_label <- get_ont_label(input$gsea_ontology)
    if (nrow(values$gsea_results) > 0) {
      dotplot(values$gsea_results, showCategory = 20) + ggtitle(paste("GSEA:", ont_label))
    } else {
      plot(NULL, xlim = c(0, 1), ylim = c(0, 1),
           main = "No significant enrichment found.", xaxt = 'n', yaxt = 'n')
    }
  }, height = 700)

  observeEvent(input$fullscreen_gsea_emap, {
    showModal(modalDialog(
      title = "GSEA Enrichment Map - Fullscreen View",
      plotOutput("gsea_emapplot_fs", height = "700px"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close")
    ))
  })
  output$gsea_emapplot_fs <- renderPlot({
    req(values$gsea_results)
    if (nrow(values$gsea_results) > 0) emapplot(pairwise_termsim(values$gsea_results), showCategory = 20)
  }, height = 700)

  observeEvent(input$fullscreen_gsea_ridge, {
    showModal(modalDialog(
      title = "GSEA Ridgeplot - Fullscreen View",
      plotOutput("gsea_ridgeplot_fs", height = "700px"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close")
    ))
  })
  output$gsea_ridgeplot_fs <- renderPlot({
    req(values$gsea_results)
    if (nrow(values$gsea_results) > 0) ridgeplot(values$gsea_results)
  }, height = 700)

}
