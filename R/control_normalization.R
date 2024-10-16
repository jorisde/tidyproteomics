#####
# R/control_normalization.R
#####
#' Main function for normalizing quantitative data in a tidyproteomics data-object
#'
#' @description
#' `normalize()` Main function for normalizing quantitative data from a tidyproteomics
#' data-object. This is a _passthrough_ function as it returns the original
#' tidyproteomics data-object with an additional quantitative column labeled with the
#' normalization method(s) used.
#'
#' This function can accommodate multiple normalization methods in a single pass,
#' and it is useful for examining normalization effects on data. Often it is
#' adventitious to select a optimal normalization method based on performance.
#'
#' @param data tidyproteomics data object
#' @param ... use a subset of the data for normalization see `subset()`. This is useful when normalizing against a spike-in set of proteins
#' @param .method character vector of normalization to use
#' @param .cores number of CPU cores to use for multi-threading
#'
#' @return a tidyproteomics data-object
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyproteomics)
#' hela_proteins %>%
#'      normalize(.method = c("scaled", "median", "EigenMS")) %>%
#'      summary("sample")
#'
#' # normalize between samples according to a subset, then apply to all values
#' #   this would be recommended with a pull-down experiment wherein a conserved
#' #   protein complex acts as the majority content and individual inter-actors
#' #   are of quantitative differentiation
#' hela_proteins %>%
#'      normalize(!description %like% "Ribosome", .method = c("scaled", "median", "EigenMS")) %>%
#'      summary("sample")
#'
normalize <- function(
    data,
    ...,
    .method = c("scaled", "median", "linear", "limma", "loess", "svm", "randomforest", "EigenMS"),
    .cores = 1
){
  
  # Load required packages
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tidyr", quietly = TRUE)
  requireNamespace("cli", quietly = TRUE)
  requireNamespace("glue", quietly = TRUE)
  
  # visible bindings
  abundance_normalized <- NULL
  abundance_log2 <- NULL
  
  # Match and validate the normalization methods
  .method <- rlang::arg_match(.method, multiple = TRUE)
  check_data(data)
  
  use_quo <- tidyproteomics_quo(...)
  
  # Handle specific normalization method constraints
  if(intersect(.method, c("tmt-bridge")) %>% length() == 1){
    if(!is.null(use_quo)){
      cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
      cli::cli_abort(c("x" = "{.emph tmt-bridge} cannot be used with a subset for normalization"))
    }
  
    if(length(.method) > 1) {
      cli::cli_div(theme = list(span.emph = list(color = "#ff4500")))
      cli::cli_abort(c("x" = "{.emph tmt-bridge} can not be used along side other normalizations"))
    }
  }
  
  if(intersect(.method, c("limma")) %>% length() == 1 & !is.null(use_quo)) {
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500"), span.info = list(color = "blue")))
    cli::cli_alert_warning("{.emph limma} cannot be used with a subset for normalization")
    cli::cli_alert_info("  ignoring {.info {paste(use_quo, collapse = ' ')}} for limma")
  }
  
  # Extract and transform data
  d <- dc <- data %>%
    extract(values = 'raw') %>%
    transform_log2(values = 'abundance') %>%
    dplyr::rename(abundance = abundance_log2) %>%
    dplyr::select(!dplyr::matches("^origin$"))
  
  # Append normalization methods to operations
  data$operations <- append(data$operations, glue::glue("Data normalized via {paste(.method, collapse=', ')}."))
  
  # Define a set of factors for normalization
  if(!is.null(use_quo)){
  
    data$operations <- append(data$operations, glue::glue("Data normalized using subset {paste(use_quo, collapse = ' ')}."))
  
    pre_n <- dc$identifier %>% unique() %>% length()
    pre_range <- dc$abundance %>% range()
  
    dc <- data %>%
      subset(..., .verbose = FALSE) %>%
      extract(values = 'raw') %>%
      transform_log2(values = 'abundance') %>%
      dplyr::rename(abundance = abundance_log2) %>%
      dplyr::select(!dplyr::matches("^origin$"))
  
    pst_n <- dc$identifier %>% unique() %>% length()
    pst_range <- dc$abundance %>% range(na.rm = TRUE)
  
    cli::cli_alert_warning("  normalization based on {pst_n} of {pre_n} identifiers")
    data$operations <- append(data$operations, glue::glue(" ... based on a subset of {pst_n} out of {pre_n} identifiers"))
    if(pre_range[1] < pst_range[1] | pre_range[2] > pst_range[2]){
      cli::cli_alert_warning("  {.emph WARNING}: filter narrowed range, NAs may result")
      cli::cli_alert_warning("  {.emph WARNING}: omitting `limma` and `randomforest` - can not accomidate subsetting")
    }
  
  }
  
  cli::cli_alert_info("Normalizing quantitative data")
  cli::cli_progress_bar(type = 'tasks')
  
  # Loop over each normalization method
  for(m in .method){
    start_time <- Sys.time()
  
    # Initialize default grouping and centering methods
    g <- 'identifier'
    f <- 'regression'
    c <- 'median'
    if(m %in% c('median','scaled')) {
      f <- 'shift'
      g <- c('sample', 'replicate')
      if(m == 'scaled') {
        c <- 'sum'
      }
    }
  
    # Skip certain methods if subsetting is used
    if(!is.null(use_quo) & m %in% c('limma','randomforest')) {next}
  
    cli::cli_div(theme = list(span.emph = list(color = "#ff4500"), span.info = list(color = "blue")))
    cli::cli_progress_step(" ... using {.info {m} {f}}")
  
    # Center the data based on the method
    dc_this <- dc %>%
      center(group_by = g, values = 'abundance', method = c) %>%
      dplyr::rename(abundance = tidyselect::matches('abundance'))
  
    # Remove previous normalized variable if it exists
    data$quantitative <- data$quantitative %>% dplyr::select(!tidyselect::matches(paste0("\\_", m)))
  
    # Perform normalization based on the method
    if(m == 'median') {         
      d_norm <- d %>% normalize_median(dc_this)
    }
    if(m == 'scaled') {         
      d_norm <- d %>% normalize_scaled(dc_this)
    }
    if(m == 'linear') {         
      d_norm <- d %>% normalize_linear(dc_this)
    }
    if(m == 'loess') {          
      d_norm <- d %>% normalize_loess(dc_this)
    }
    if(m == 'limma') {          
      d_norm <- d %>% normalize_limma()
    }
    if(m == 'svm') {            
      d_norm <- d %>% normalize_svm(dc_this, .cores)
    }
    if(m == 'randomforest') {   
      d_norm <- d %>% normalize_randomforest(dc_this, .cores)
    }
    if(m == 'EigenMS') {       
      d_norm <- normalize_EigenMS(data, .cores)
    }
  
    # Invert log2 transformation and rename
    d_norm <- d_norm %>%
      dplyr::mutate(abundance_normalized = invlog2(abundance_normalized)) %>%
      dplyr::rename(abundance = abundance_normalized)
  
    # Merge the normalized data back into the main data object
    data <- data %>% merge_quantitative(d_norm, m)
  
    cli::cli_progress_done()
  }
  
  # Select normalization columns
  data <- data %>% select_normalization()
  
  return(data)
}
