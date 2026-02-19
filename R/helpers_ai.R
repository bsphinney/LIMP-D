# ==============================================================================
#  HELPER FUNCTIONS — Google Gemini AI Integration
# ==============================================================================

# --- CHECK AVAILABLE MODELS ---
list_google_models <- function(api_key) {
  req <- request("https://generativelanguage.googleapis.com/v1beta/models") %>%
    req_url_query(key = api_key)
  tryCatch({
    resp <- req_perform(req)
    data <- resp_body_json(resp)
    models <- sapply(data$models, function(x) x$name)
    models <- gsub("^models/", "", models)
    return(models)
  }, error = function(e) { return(paste("Error listing models:", e$message)) })
}

# --- FILE API UPLOADER ---
upload_csv_to_gemini <- function(df, api_key) {
  temp_file <- tempfile(fileext = ".csv")
  write.csv(df, temp_file, row.names = FALSE)
  file_size <- file.size(temp_file)

  req <- request("https://generativelanguage.googleapis.com/upload/v1beta/files") %>%
    req_url_query(key = api_key) %>%
    req_headers(
      "X-Goog-Upload-Protocol" = "raw",
      "X-Goog-Upload-Command" = "start, upload, finalize",
      "X-Goog-Upload-Header-Content-Length" = as.character(file_size),
      "X-Goog-Upload-Header-Content-Type" = "text/csv",
      "Content-Type" = "text/csv"
    ) %>%
    req_body_file(temp_file)

  resp <- req_perform(req)
  file_info <- resp_body_json(resp)
  return(file_info$file$uri)
}

# --- AI CHAT FUNCTION ---
ask_gemini_file_chat <- function(user_query, file_uri, qc_df, api_key, model_name, selected_ids = NULL) {

  qc_text <- "No QC Data Available"
  if(!is.null(qc_df)) {
    qc_text <- paste(capture.output(write.csv(qc_df, row.names=FALSE)), collapse="\n")
  }

  selection_context <- ""
  if (!is.null(selected_ids) && length(selected_ids) > 0) {
    selection_context <- paste0(
      "\n!!! URGENT: USER SELECTION ACTIVE !!!\n",
      "Focus analysis on these specific proteins:\n",
      paste(selected_ids, collapse=", "), "\n"
    )
  }

  system_instruction <- paste0(
    "You are a PhD-level expert in proteomics. ",
    "You have access to two data sources:\n",
    "SOURCE 1: QC METRICS (In Text Below)\n",
    "This table includes 'Group' columns. Use it to compare technical quality (Precursors, MS1) between experimental groups.\n",
    "--- START QC DATA ---\n",
    qc_text,
    "\n--- END QC DATA ---\n\n",
    "SOURCE 2: EXPRESSION DATA (In Uploaded File)\n",
    "Use this file to answer biological questions. It contains the Top 800 proteins.\n\n",
    "IMPORTANT: BI-DIRECTIONAL CONTROL.\n",
    "1. If the user asks about 'selected proteins', refer to the 'URGENT' section below.\n",
    "2. If you find interesting proteins, OUTPUT their IDs at the end like this:\n",
    "   [[SELECT: P12345; P67890]]\n"
  )

  base_url <- "https://generativelanguage.googleapis.com/v1beta/models/"
  clean_model <- gsub("^models/", "", model_name)
  full_url <- paste0(base_url, clean_model, ":generateContent")

  body <- list(contents = list(list(parts = list(
    list(text = paste0(system_instruction, selection_context, "\n\nUser Question: ", user_query)),
    list(file_data = list(file_uri = file_uri, mime_type = "text/csv"))
  ))))

  req <- request(full_url) %>%
    req_url_query(key = api_key) %>%
    req_headers("Content-Type" = "application/json") %>%
    req_body_json(body)

  tryCatch({
    resp <- req_perform(req)
    return(resp_body_json(resp)$candidates[[1]]$content$parts[[1]]$text)
  }, error = function(e) {
    err_msg <- "Unknown Error"
    if (!is.null(e$resp)) { err_msg <- tryCatch(resp_body_string(e$resp), error = function(z) e$message) } else { err_msg <- e$message }
    return(paste("API Error:", err_msg))
  })
}

# --- AI TEXT CHAT FUNCTION ---
ask_gemini_text_chat <- function(user_query, api_key, model_name) {
  base_url <- "https://generativelanguage.googleapis.com/v1beta/models/"
  clean_model <- gsub("^models/", "", model_name)
  full_url <- paste0(base_url, clean_model, ":generateContent")

  body <- list(contents = list(list(parts = list(list(text = user_query)))))

  req <- request(full_url) %>%
    req_url_query(key = api_key) %>%
    req_headers("Content-Type" = "application/json") %>%
    req_body_json(body)

  tryCatch({
    resp <- req_perform(req)
    return(resp_body_json(resp)$candidates[[1]]$content$parts[[1]]$text)
  }, error = function(e) {
    err_msg <- "Unknown Error"
    if (!is.null(e$resp)) { err_msg <- tryCatch(resp_body_string(e$resp), error = function(z) e$message) } else { err_msg <- e$message }
    return(paste("API Error:", err_msg))
  })
}
