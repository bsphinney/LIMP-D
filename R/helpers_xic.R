# ==============================================================================
#  HELPER FUNCTIONS — XIC (Extracted Ion Chromatogram) Viewer Utilities
# ==============================================================================

# Helper: Detect XIC format version from column names
detect_xic_format <- function(xic_dir) {
  xic_files <- list.files(xic_dir, pattern = "\\.xic\\.parquet$",
                          full.names = TRUE, recursive = TRUE)
  if (length(xic_files) == 0) return("unknown")
  cols <- tryCatch(names(arrow::read_parquet(xic_files[1], as_data_frame = FALSE)),
                   error = function(e) character(0))
  if (all(c("pr", "feature", "rt", "value") %in% cols)) return("v2")
  if ("Precursor.Id" %in% cols) return("v1")
  return("unknown")
}

# Helper: Load XIC data for a single protein using Arrow predicate pushdown
# Supports both DIA-NN 1.x (wide format) and 2.x (long format: pr/feature/rt/value)
load_xic_for_protein <- function(xic_dir, protein_id, report_map, xic_format = "v2") {
  # In DIA-NN 2.x, Precursor.Id in the report already matches the XIC "pr" column
  # (both use StrippedSequence+Charge format, e.g., "PEPTIDEK2")
  # In DIA-NN 1.x, Precursor.Id matches the XIC "Precursor.Id" column directly
  target_prs <- report_map %>%
    filter(Protein.Group == protein_id) %>%
    distinct(Precursor.Id) %>%
    pull(Precursor.Id)

  if (length(target_prs) == 0) return(NULL)

  # Read only .xic.parquet files (exclude mobilograms)
  xic_files <- list.files(xic_dir, pattern = "\\.xic\\.parquet$",
                          full.names = TRUE, recursive = TRUE)

  # Per-file read with source file tagging (needed for v2 which has no File.Name)
  xic_list <- lapply(xic_files, function(f) {
    tryCatch({
      df <- arrow::read_parquet(f, as_data_frame = TRUE)
      # Filter by precursor — use base R to avoid rlang issues
      if (xic_format == "v2") {
        df <- df[df$pr %in% target_prs, , drop = FALSE]
      } else {
        df <- df[df$Precursor.Id %in% target_prs, , drop = FALSE]
      }
      if (nrow(df) > 0) {
        sample_name <- sub("\\.xic\\.parquet$", "", basename(f))
        df$Source.File <- sample_name
      }
      df
    }, error = function(e2) NULL)
  })
  result <- bind_rows(Filter(function(x) !is.null(x) && nrow(x) > 0, xic_list))
  if (nrow(result) == 0) return(NULL)
  result
}

# Helper: Reshape XIC data for plotting — handles both v1 and v2 formats
reshape_xic_for_plotting <- function(xic_raw, metadata, xic_format = "v2") {

  if (xic_format == "v2") {
    # -- DIA-NN 2.x: already long format (pr, feature, info, rt, value) --
    # Rename columns using base R to avoid tidy evaluation issues with "pr"
    names(xic_raw)[names(xic_raw) == "pr"] <- "Precursor.Id"
    names(xic_raw)[names(xic_raw) == "feature"] <- "Fragment.Label"
    names(xic_raw)[names(xic_raw) == "rt"] <- "RT"
    names(xic_raw)[names(xic_raw) == "value"] <- "Intensity"

    xic_plot <- xic_raw %>%
      mutate(
        MS.Level = ifelse(Fragment.Label == "ms1", 1L, 2L),
        RT = as.numeric(RT),
        Intensity = as.numeric(Intensity)
      ) %>%
      filter(!(Intensity == 0 & RT == 0))

    # Match Source.File to metadata File.Name via fuzzy basename matching
    if ("Source.File" %in% names(xic_plot)) {
      meta_lookup <- metadata %>%
        mutate(File.Name.Base = basename(tools::file_path_sans_ext(File.Name)))
      xic_plot <- xic_plot %>%
        mutate(File.Name = Source.File,
               File.Name.Base = Source.File) %>%
        left_join(
          meta_lookup %>% dplyr::select(File.Name.Base, Group, ID),
          by = "File.Name.Base"
        )
    }

    # Keep only the columns we need — use base R to avoid arrow::select conflict
    keep_cols <- intersect(c("File.Name", "ID", "Group", "Precursor.Id",
                             "MS.Level", "Fragment.Label", "RT", "Intensity"),
                           names(xic_plot))
    xic_plot <- xic_plot[, keep_cols, drop = FALSE]

    return(xic_plot)

  } else {
    # -- DIA-NN 1.x: wide format with numbered columns --
    num_cols <- names(xic_raw)[grepl("^\\d+$", names(xic_raw))]
    if (length(num_cols) == 0) {
      warning("No numbered columns found in XIC data")
      return(NULL)
    }

    make_key <- function(df) {
      paste(df$File.Name, df$Precursor.Id, df$MS.Level,
            df$Theoretical.Mz, df$FragmentType, df$FragmentCharge,
            df$FragmentSeriesNumber, df$FragmentLossType, sep = "|")
    }

    rt_rows <- xic_raw %>% filter(Retention.Times == 1)
    int_rows <- xic_raw %>% filter(Intensities == 1)

    rt_long <- rt_rows %>%
      mutate(.key = make_key(rt_rows)) %>%
      dplyr::select(.key, all_of(num_cols)) %>%
      pivot_longer(cols = all_of(num_cols), names_to = "point_idx", values_to = "RT")

    int_long <- int_rows %>%
      mutate(.key = make_key(int_rows)) %>%
      dplyr::select(.key, File.Name, Precursor.Id, Modified.Sequence, MS.Level,
             Theoretical.Mz, Reference.Intensity, FragmentType, FragmentCharge,
             FragmentSeriesNumber, FragmentLossType, all_of(num_cols)) %>%
      pivot_longer(cols = all_of(num_cols), names_to = "point_idx", values_to = "Intensity")

    xic_plot <- inner_join(
      rt_long %>% dplyr::select(.key, point_idx, RT),
      int_long,
      by = c(".key", "point_idx")
    ) %>%
      mutate(
        RT = as.numeric(RT),
        Intensity = as.numeric(Intensity),
        Fragment.Label = case_when(
          MS.Level == 1 ~ paste0("MS1 (", round(as.numeric(Theoretical.Mz), 2), ")"),
          TRUE ~ paste0(FragmentType, FragmentSeriesNumber,
                        ifelse(as.integer(FragmentCharge) > 1,
                               paste0("+", FragmentCharge), ""),
                        ifelse(FragmentLossType != "noloss",
                               paste0("-", FragmentLossType), ""))
        )
      ) %>%
      filter(!(Intensity == 0 & RT == 0)) %>%
      mutate(File.Name.Base = basename(tools::file_path_sans_ext(
        tools::file_path_sans_ext(File.Name)))) %>%
      left_join(
        metadata %>%
          mutate(File.Name.Base = basename(tools::file_path_sans_ext(File.Name))) %>%
          dplyr::select(File.Name.Base, Group, ID),
        by = "File.Name.Base"
      ) %>%
      dplyr::select(File.Name, ID, Group, Precursor.Id, Modified.Sequence,
             MS.Level, Fragment.Label, Theoretical.Mz, Reference.Intensity,
             RT, Intensity)

    return(xic_plot)
  }
}
