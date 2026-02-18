server_gsea <- function(input, output, session, values, add_to_log) {

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
        tags$h6("Gene Ontology (GO)"),
        p("This analysis uses GO Biological Process terms \u2014 curated annotations describing what biological functions ",
          "genes are involved in (e.g., 'cell division', 'immune response', 'protein folding')."),
        tags$h6("The visualization tabs"),
        tags$ul(
          tags$li(strong("Dot Plot: "), "Each dot = one enriched pathway. Size = number of genes, color = adjusted p-value. X-axis shows gene ratio (fraction of pathway genes affected)."),
          tags$li(strong("Enrichment Map: "), "Network showing how enriched pathways relate to each other. Connected pathways share genes."),
          tags$li(strong("Ridgeplot: "), "Density curves showing the fold-change distribution of genes within each enriched pathway."),
          tags$li(strong("Results Table: "), "Full numeric results with statistics for each pathway.")
        ),
        tags$h6("Which comparison is used?"),
        p("GSEA uses the currently selected comparison from the DE Dashboard. ",
          "Change the comparison selector on the DE Dashboard to run GSEA on a different contrast.")
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
          tags$li(strong("ID: "), "Gene Ontology accession (e.g., GO:0006955)"),
          tags$li(strong("Description: "), "Human-readable name of the biological process"),
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

  # --- 7. GSEA Dot Plot ---
  observeEvent(input$fullscreen_gsea_dot, {
    showModal(modalDialog(
      title = "GSEA Dot Plot - Fullscreen View",
      plotOutput("gsea_dot_plot_fs", height = "700px"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close")
    ))
  })
  output$gsea_dot_plot_fs <- renderPlot({
    req(values$gsea_results)
    if (nrow(values$gsea_results) > 0) dotplot(values$gsea_results, showCategory = 20) + ggtitle("GSEA GO Biological Process")
    else plot(NULL, xlim = c(0, 1), ylim = c(0, 1), main = "No significant enrichment found.", xaxt = 'n', yaxt = 'n')
  }, height = 700)

  # --- 8. GSEA Enrichment Map ---
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

  # --- 9. GSEA Ridgeplot ---
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

  observeEvent(input$run_gsea, {
    req(values$fit, input$contrast_selector)
    output$gsea_status <- renderText("Running GSEA... This may take a few minutes.")
    withProgress(message = "Running GSEA", {
      tryCatch({
        incProgress(0.1, detail = "Preparing gene list...")
        de_results <- topTable(values$fit, coef = input$contrast_selector, number = Inf)
        org_db_name <- detect_organism_db(rownames(de_results))
        if (!requireNamespace(org_db_name, quietly = TRUE)) { BiocManager::install(org_db_name, ask = FALSE) }
        library(org_db_name, character.only = TRUE)
        output$gsea_status <- renderText(paste("Detected Organism DB:", org_db_name))
        incProgress(0.3, detail = "Converting Gene IDs...")
        protein_ids <- str_extract(rownames(de_results), "[A-Z0-9]+")
        id_map <- bitr(protein_ids, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = get(org_db_name))
        de_results_mapped <- de_results[rownames(de_results) %in% id_map$UNIPROT, ]; de_results_mapped$ENTREZID <- id_map$ENTREZID[match(str_extract(rownames(de_results_mapped), "[A-Z0-9]+"), id_map$UNIPROT)]
        gene_list <- de_results_mapped$logFC; names(gene_list) <- de_results_mapped$ENTREZID
        gene_list <- sort(gene_list, decreasing = TRUE); gene_list <- gene_list[!duplicated(names(gene_list))]
        incProgress(0.6, detail = "Running GO enrichment...")
        gsea_res <- gseGO(geneList=gene_list, OrgDb=get(org_db_name), keyType="ENTREZID", ont="BP", minGSSize=10, maxGSSize=500, pvalueCutoff=0.05, verbose=FALSE)
        values$gsea_results <- gsea_res; incProgress(1, detail = "Complete.")
        output$gsea_status <- renderText("GSEA Complete.")

        # Log GSEA analysis
        gsea_code <- c(
          sprintf("# Contrast: %s", input$contrast_selector),
          sprintf("# Organism DB: %s", org_db_name),
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
          "                  ont='BP', minGSSize=10, maxGSSize=500, pvalueCutoff=0.05)"
        )
        add_to_log("GSEA Analysis", gsea_code)

      }, error = function(e) { output$gsea_status <- renderText(paste("GSEA Error:", e$message)) })
    })
  })

  output$gsea_dot_plot <- renderPlot({ req(values$gsea_results); if(nrow(values$gsea_results) > 0) dotplot(values$gsea_results, showCategory = 20) + ggtitle("GSEA GO Biological Process") else plot(NULL, xlim=c(0,1), ylim=c(0,1), main="No significant enrichment found.", xaxt='n', yaxt='n') }) # Height controlled by UI (calc(100vh - 340px))
  output$gsea_emapplot <- renderPlot({ req(values$gsea_results); if(nrow(values$gsea_results) > 0) emapplot(pairwise_termsim(values$gsea_results), showCategory=20) }) # Height controlled by UI (calc(100vh - 340px))
  output$gsea_ridgeplot <- renderPlot({ req(values$gsea_results); if(nrow(values$gsea_results) > 0) ridgeplot(values$gsea_results) }) # Height controlled by UI (calc(100vh - 340px))
  output$gsea_results_table <- renderDT({ req(values$gsea_results); datatable(as.data.frame(values$gsea_results), options=list(pageLength=10, scrollX=TRUE)) })

}
