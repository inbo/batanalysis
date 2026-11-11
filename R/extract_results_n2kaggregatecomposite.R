#' @export
#' @importFrom assertthat assert_that
#' @importFrom dplyr bind_cols case_when group_by inner_join mutate select
#' @importFrom dplyr starts_with summarise
#' @importFrom git2rdata is_git2rdata update_metadata write_vc
#' @importFrom n2kanalysis get_file_fingerprint status
#' @importFrom rlang .data
#' @importFrom stats quantile
#' @importFrom tidyr pivot_longer
extract_results.n2kAggregateComposite <- function(x, root, ...) {
  assert_that(inherits(x, "n2kAggregateComposite"))
  if (status(x) != "converged") {
    return(invisible(NULL))
  }
  x@AggregatedImputed@Covariate |>
    bind_cols(x@AggregatedImputed@Imputation) |>
    pivot_longer(
      starts_with("Imputation"),
      names_to = "imputation",
      values_to = "value"
    ) -> results
  if (x@AnalysisMetadata$formula %in% c("~location_id", "~location")) {
    results |>
      group_by(.data$winter, .data$province, .data$type, .data$location_id) |>
      summarise(
        across(
          "value",
          .names = "{.fn}",
          list(
            mean = ~ mean(.x, na.rm = TRUE),
            sd = ~ sd(.x, na.rm = TRUE),
            min = ~ min(.x, na.rm = TRUE),
            p05 = ~ quantile(.x, 0.05, na.rm = TRUE),
            p20 = ~ quantile(.x, 0.20, na.rm = TRUE),
            p35 = ~ quantile(.x, 0.35, na.rm = TRUE),
            p50 = ~ quantile(.x, 0.50, na.rm = TRUE),
            p65 = ~ quantile(.x, 0.65, na.rm = TRUE),
            p80 = ~ quantile(.x, 0.80, na.rm = TRUE),
            p95 = ~ quantile(.x, 0.95, na.rm = TRUE),
            max = ~ max(.x, na.rm = TRUE)
          )
        ),
        .groups = "drop"
      ) |>
      bind_cols(
        x@AnalysisMetadata |>
          select(
            species = "species_group_id",
            "model_type",
            analysis = "file_fingerprint",
            fingerprint = "status_fingerprint"
          )
      ) |>
      write_vc(
        file = file.path("hibernation", "composite_location"),
        root = root,
        optimize = FALSE,
        append = TRUE,
        sorting = c(
          "model_type",
          "species",
          "winter",
          "province",
          "type",
          "location_id",
          "analysis"
        ),
        digits = 4
      )
    rm(x, results)
    gc(verbose = FALSE)
    return(invisible(NULL))
  }
  if (x@AnalysisMetadata$formula == "~province") {
    results |>
      group_by(.data$winter, .data$province) |>
      summarise(
        across(
          "value",
          .names = "{.fn}",
          list(
            mean = ~ mean(.x, na.rm = TRUE),
            sd = ~ sd(.x, na.rm = TRUE),
            min = ~ min(.x, na.rm = TRUE),
            p05 = ~ quantile(.x, 0.05, na.rm = TRUE),
            p20 = ~ quantile(.x, 0.20, na.rm = TRUE),
            p35 = ~ quantile(.x, 0.35, na.rm = TRUE),
            p50 = ~ quantile(.x, 0.50, na.rm = TRUE),
            p65 = ~ quantile(.x, 0.65, na.rm = TRUE),
            p80 = ~ quantile(.x, 0.80, na.rm = TRUE),
            p95 = ~ quantile(.x, 0.95, na.rm = TRUE),
            max = ~ max(.x, na.rm = TRUE)
          )
        ),
        .groups = "drop"
      ) |>
      bind_cols(
        x@AnalysisMetadata |>
          select(
            species = "species_group_id",
            "model_type",
            analysis = "file_fingerprint",
            fingerprint = "status_fingerprint"
          )
      ) |>
      write_vc(
        file = file.path("hibernation", "composite_province"), root = root,
        optimize = FALSE, append = TRUE,
        sorting = c("model_type", "species", "winter", "province", "analysis"),
        digits = 4
      )
    rm(x, results)
    gc(verbose = FALSE)
    return(invisible(NULL))
  }
  if (x@AnalysisMetadata$formula == "~type") {
    results |>
      group_by(.data$winter, .data$type) |>
      summarise(
        across(
          "value",
          .names = "{.fn}",
          list(
            mean = ~ mean(.x, na.rm = TRUE),
            sd = ~ sd(.x, na.rm = TRUE),
            min = ~ min(.x, na.rm = TRUE),
            p05 = ~ quantile(.x, 0.05, na.rm = TRUE),
            p20 = ~ quantile(.x, 0.20, na.rm = TRUE),
            p35 = ~ quantile(.x, 0.35, na.rm = TRUE),
            p50 = ~ quantile(.x, 0.50, na.rm = TRUE),
            p65 = ~ quantile(.x, 0.65, na.rm = TRUE),
            p80 = ~ quantile(.x, 0.80, na.rm = TRUE),
            p95 = ~ quantile(.x, 0.95, na.rm = TRUE),
            max = ~ max(.x, na.rm = TRUE)
          )
        ),
        .groups = "drop"
      ) |>
      bind_cols(
        x@AnalysisMetadata |>
          select(
            species = "species_group_id",
            "model_type",
            analysis = "file_fingerprint",
            fingerprint = "status_fingerprint"
          )
      ) |>
      write_vc(
        file = file.path("hibernation", "composite_type"), root = root,
        optimize = FALSE, append = TRUE,
        sorting = c("model_type", "species", "winter", "type", "analysis"),
        digits = 4
      )
    rm(x, results)
    gc(verbose = FALSE)
    return(invisible(NULL))
  }
  if (x@AnalysisMetadata$formula == "~tot") {
    results |>
      group_by(.data$winter) |>
      summarise(
        across(
          "value",
          .names = "{.fn}",
          list(
            mean = ~ mean(.x, na.rm = TRUE),
            sd = ~ sd(.x, na.rm = TRUE),
            min = ~ min(.x, na.rm = TRUE),
            p05 = ~ quantile(.x, 0.05, na.rm = TRUE),
            p20 = ~ quantile(.x, 0.20, na.rm = TRUE),
            p35 = ~ quantile(.x, 0.35, na.rm = TRUE),
            p50 = ~ quantile(.x, 0.50, na.rm = TRUE),
            p65 = ~ quantile(.x, 0.65, na.rm = TRUE),
            p80 = ~ quantile(.x, 0.80, na.rm = TRUE),
            p95 = ~ quantile(.x, 0.95, na.rm = TRUE),
            max = ~ max(.x, na.rm = TRUE)
          )
        ),
        .groups = "drop"
      ) |>
      bind_cols(
        x@AnalysisMetadata |>
          select(
            species = "species_group_id",
            "model_type",
            analysis = "file_fingerprint",
            fingerprint = "status_fingerprint"
          )
      ) |>
      write_vc(
        file = file.path("hibernation", "composite_total"),
        root = root,
        optimize = FALSE,
        append = TRUE,
        sorting = c("model_type", "species", "winter", "analysis"),
        digits = 4
      )
    rm(x, results)
    gc(verbose = FALSE)
    return(invisible(NULL))
  }
  stop("n2kAggregateComposite todo")
}
