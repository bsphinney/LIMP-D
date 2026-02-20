# ==============================================================================
#  helpers_mofa.R
#  Pure utility functions for MOFA2 Multi-View Integration.
#  No Shiny reactivity — called from server_mofa.R.
# ==============================================================================

# --- Default view configuration factory ---
default_view_config <- function(id, num) {
  list(
    id = id,
    num = num,
    name = paste("View", num),
    type = "proteomics_other",
    source = "upload",
    file_path = NULL,
    file_name = NULL,
    matrix = NULL,
    fit = NULL,
    n_features = 0,
    n_samples = 0,
    status = "pending"
  )
}

# --- Smart RDS Parser ---
# Detects: DE-LIMP session, limma EList, limma MArrayLM, matrix, named list, data.frame
parse_rds_for_mofa <- function(file_path) {
  obj <- readRDS(file_path)

  # Case 1: Full DE-LIMP session
  if (is.list(obj) && "raw_data" %in% names(obj)) {
    mat <- obj$raw_data
    return(list(
      matrix = mat,
      fit = obj$fit,
      metadata = obj$metadata,
      format = "delimp_session",
      source = paste("DE-LIMP", obj$session_info$app_version %||% ""),
      message = sprintf("Loaded DE-LIMP session: %d features x %d samples",
                        nrow(mat), ncol(mat))
    ))
  }

  # Case 2: Limma EList object
  if (inherits(obj, "EList")) {
    mat <- obj$E
    return(list(
      matrix = mat,
      fit = NULL,
      metadata = NULL,
      format = "limma_elist",
      source = "limma EList",
      message = sprintf("Loaded limma EList: %d features x %d samples",
                        nrow(mat), ncol(mat))
    ))
  }

  # Case 3: Limma MArrayLM (fit object)
  if (inherits(obj, "MArrayLM")) {
    mat <- obj$coefficients
    return(list(
      matrix = mat,
      fit = obj,
      metadata = NULL,
      format = "limma_fit",
      source = "limma fit object",
      message = sprintf("Loaded limma fit: %d features x %d coefficients",
                        nrow(mat), ncol(mat))
    ))
  }

  # Case 4: Simple matrix
  if (is.matrix(obj)) {
    return(list(
      matrix = obj,
      fit = NULL,
      metadata = NULL,
      format = "matrix_rds",
      source = "R matrix",
      message = sprintf("Loaded matrix: %d features x %d samples",
                        nrow(obj), ncol(obj))
    ))
  }

  # Case 5: Named list with matrix
  if (is.list(obj) && "matrix" %in% names(obj)) {
    mat <- obj$matrix
    return(list(
      matrix = mat,
      fit = obj$fit,
      metadata = obj$metadata,
      format = "list_rds",
      source = obj$source %||% "R list",
      message = sprintf("Loaded from list: %d features x %d samples",
                        nrow(mat), ncol(mat))
    ))
  }

  # Case 6: Data frame (convert to matrix)
  if (is.data.frame(obj)) {
    row_ids <- obj[[1]]
    mat <- as.matrix(obj[, -1, drop = FALSE])
    rownames(mat) <- row_ids
    return(list(
      matrix = mat,
      fit = NULL,
      metadata = NULL,
      format = "dataframe_rds",
      source = "R data.frame",
      message = sprintf("Loaded data.frame: %d features x %d samples",
                        nrow(mat), ncol(mat))
    ))
  }

  stop("Unrecognized RDS structure. Expected: matrix, data.frame, limma object, or DE-LIMP session.")
}

# --- Generic Matrix Parser (CSV, TSV, Parquet) ---
parse_matrix_file <- function(file_path, file_name) {
  ext <- tolower(tools::file_ext(file_name))

  mat <- switch(ext,
    "parquet" = {
      df <- arrow::read_parquet(file_path)
      # Check if it's a DIA-NN report
      if (all(c("Protein.Group", "Run", "Precursor.Normalised") %in% names(df))) {
        stop("This appears to be a DIA-NN report. Please select 'Upload DIA-NN report' instead.")
      }
      row_ids <- df[[1]]
      m <- as.matrix(df[, -1, drop = FALSE])
      rownames(m) <- row_ids
      m
    },
    "csv" = {
      df <- read.csv(file_path, row.names = 1, check.names = FALSE)
      as.matrix(df)
    },
    "tsv" = , "txt" = {
      df <- read.delim(file_path, row.names = 1, check.names = FALSE)
      as.matrix(df)
    },
    stop(sprintf("Unsupported file type: .%s", ext))
  )

  # Auto-detect if log transformation needed
  max_val <- max(mat, na.rm = TRUE)
  was_transformed <- FALSE

  if (max_val > 100) {
    mat <- log2(mat + 1)
    was_transformed <- TRUE
  }

  list(
    matrix = mat,
    format = paste0("matrix_", ext),
    transformed = was_transformed,
    message = sprintf("Loaded %s: %d features x %d samples%s",
                      ext, nrow(mat), ncol(mat),
                      if (was_transformed) " (log2 transformed)" else "")
  )
}

# --- Sample Overlap Computation ---
compute_sample_overlap <- function(views) {
  if (length(views) < 2) return(NULL)

  view_samples <- lapply(views, colnames)
  common_samples <- Reduce(intersect, view_samples)
  all_samples <- Reduce(union, view_samples)

  n_common <- length(common_samples)
  n_total <- length(all_samples)

  # Per-view breakdown
  per_view <- lapply(names(views), function(vname) {
    vsamp <- colnames(views[[vname]])
    list(
      name = vname,
      total = length(vsamp),
      shared = sum(vsamp %in% common_samples),
      unique = length(setdiff(vsamp, common_samples))
    )
  })

  list(
    common_samples = common_samples,
    all_samples = all_samples,
    n_common = n_common,
    n_total = n_total,
    pct_overlap = if (n_total > 0) 100 * n_common / n_total else 0,
    per_view = per_view
  )
}

# --- Reproducibility Code Generator ---
generate_mofa_code <- function(params) {
  c(
    "# ============================================================",
    "# MULTI-OMICS FACTOR ANALYSIS (MOFA2)",
    "# ============================================================",
    sprintf("# Date: %s", params$timestamp),
    sprintf("# Views: %d", params$n_views),
    sprintf("# Samples: %d", params$n_samples),
    sprintf("# Factors: %d", params$n_factors),
    "",
    "library(MOFA2)",
    "",
    "# Prepare data (list of matrices, features x samples)",
    "mofa_data <- list(",
    paste0("  # ", paste(params$view_names, collapse = ", ")),
    ")",
    "",
    "# Create MOFA object",
    "mofa_obj <- create_mofa(mofa_data)",
    "",
    "# Configure options",
    "data_opts <- get_default_data_options(mofa_obj)",
    sprintf("data_opts$scale_views <- %s", params$scale_views),
    "",
    "model_opts <- get_default_model_options(mofa_obj)",
    sprintf("model_opts$num_factors <- %d", params$n_factors),
    "",
    "train_opts <- get_default_training_options(mofa_obj)",
    sprintf("train_opts$convergence_mode <- '%s'", params$convergence),
    sprintf("train_opts$seed <- %d", params$seed),
    "",
    "mofa_obj <- prepare_mofa(mofa_obj,",
    "  data_options = data_opts,",
    "  model_options = model_opts,",
    "  training_options = train_opts",
    ")",
    "",
    "# Train model",
    "mofa_trained <- run_mofa(mofa_obj, use_basilisk = TRUE)",
    "",
    "# Extract results",
    "factors <- get_factors(mofa_trained)",
    "weights <- get_weights(mofa_trained)",
    "variance <- get_variance_explained(mofa_trained)"
  )
}
