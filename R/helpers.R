# ==============================================================================
#  HELPER FUNCTIONS — General utilities
# ==============================================================================

# --- QC Stats Calculation ---
get_diann_stats_r <- function(file_path) {
  tryCatch({
    df <- arrow::read_parquet(file_path)
    has_pg_q <- "PG.Q.Value" %in% names(df)
    has_ms1  <- "Ms1.Apex.Area" %in% names(df)
    if ("Q.Value" %in% names(df)) df <- df %>% filter(Q.Value <= 0.01)

    stats_df <- df %>%
      group_by(Run) %>%
      summarise(
        Precursors = n(),
        Proteins = if(has_pg_q) {
          n_distinct(Protein.Group[PG.Q.Value <= 0.01])
        } else {
          n_distinct(Protein.Group)
        },
        MS1_Signal = if(has_ms1) sum(Ms1.Apex.Area, na.rm = TRUE) else NA_real_,
        .groups = 'drop'
      ) %>% arrange(Run)

    return(stats_df)
  }, error = function(e) { data.frame(Run = "Error", Precursors = 0, Proteins = 0, MS1_Signal = 0) })
}

# --- Z-Score Utility ---
cal_z_score <- function(x) { (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) }

# --- Auto-detect Organism ---
detect_organism_db <- function(protein_ids) {
  ORGANISM_DB_MAP <- list(
    "_HUMAN" = "org.Hs.eg.db", "_MOUSE" = "org.Mm.eg.db", "_RAT"   = "org.Rn.eg.db",
    "_BOVIN" = "org.Bt.eg.db", "_CANLF" = "org.Cf.eg.db", "_CHICK" = "org.Gg.eg.db",
    "_DROME" = "org.Dm.eg.db", "_CAEEL" = "org.Ce.eg.db", "_DANRE" = "org.Dr.eg.db",
    "_YEAST" = "org.Sc.sgd.db", "_ARATH" = "org.At.tair.db", "_PIG"   = "org.Ss.eg.db"
  )
  for (suffix in names(ORGANISM_DB_MAP)) {
    if (any(grepl(suffix, protein_ids, ignore.case = TRUE))) {
      return(ORGANISM_DB_MAP[[suffix]])
    }
  }
  return("org.Hs.eg.db")
}
