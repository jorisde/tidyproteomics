#####
# R/normalization_eigenms.R
#####

#' EigenMS Normalization Function for a tidyproteomics Data-Object
#'
#' @description
#' `normalize_EigenMS()` performs normalization using the EigenMS method on a tidyproteomics data-object.
#'
#' @param data tidyproteomics data object
#' @param .cores Number of CPU cores to use for parallel processing (default is 1)
#' @param input_dir Directory containing the input TSV files (optional; if NULL, uses data from the data object)
#' @param output_dir Directory to save EigenMS results and plots (optional; defaults to a temporary directory)
#'
#' @return A tidyproteomics data-object with EigenMS normalized abundance values
#' @export
#'
#' @examples
#' \dontrun{
#' library(tidyproteomics)
#' normalized_data <- normalize_EigenMS(data, .cores = 2)
#' }
normalize_EigenMS <- function(
    data,
    .cores = 1,
    input_dir = NULL,
    output_dir = NULL
) {
  
  # Check if EigenMS script exists
  eigenms_script_path <- system.file("R", "EigenMS.R", package = "yourPackageName") # Adjust package name
  if (!file.exists(eigenms_script_path)) {
    stop("EigenMS.R script not found in the 'scripts' directory of the package.")
  }
  
  # Source the EigenMS script
  source(eigenms_script_path)
  
  # Define default output directory if not provided
  if (is.null(output_dir)) {
    output_dir <- tempdir()
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Extract necessary data from the tidyproteomics data object
  proteomics_data <- data %>%
    extract(values = 'raw') %>%
    dplyr::select(identifier, sample, replicate, abundance)
  
  # Convert data to wide format required by EigenMS
  abundance_matrix <- proteomics_data %>%
    tidyr::pivot_wider(names_from = sample, values_from = abundance) %>%
    as.data.frame()
  
  rownames(abundance_matrix) <- abundance_matrix$identifier
  abundance_matrix <- abundance_matrix[ , -1]  # Remove identifier column
  
  # Handle missing values as EigenMS does not handle NAs
  abundance_matrix[is.na(abundance_matrix)] <- 0  # Alternatively, implement a better imputation strategy
  
  # Define treatment groups (modify as per your experimental design)
  # Example: Assuming 'replicate' contains group information
  # You may need to adjust this based on your actual data structure
  treatment <- factor(data$quantitative$replicate)
  
  # Perform EigenMS normalization
  eigenms_result <- EigenMS::EigenMS(
    M = as.matrix(abundance_matrix),
    design = treatment,
    nCore = .cores
  )
  
  # Extract normalized values
  normalized_matrix <- eigenms_result$normalizedValues
  
  # Convert back to tidy format
  normalized_data <- normalized_matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "identifier") %>%
    tidyr::pivot_longer(
      cols = -identifier,
      names_to = "sample",
      values_to = "abundance_normalized"
    )
  
  # Merge normalized data back into the main data object
  data <- data %>%
    merge_quantitative(normalized_data, method = "EigenMS")
  
  # Optionally, add EigenMS-specific operations or metadata
  data$operations <- append(data$operations, "Data normalized via EigenMS.")
  
  return(data)
}
