#' Extract results from a model
#' @param x The model to extract results from
#' @param ... Additional arguments passed to the extraction method
#' @export
#' @importFrom assertthat assert_that
extract_results <- function(x, ...) {
  UseMethod("extract_results", x)
}

#' @export
extract_results.default <- function(x, ...) {
  message("No extraction method for class ", class(x))
}
