#' @export
#' @importFrom assertthat assert_that
#' @importFrom dplyr bind_cols case_when group_by inner_join mutate select
#' @importFrom dplyr starts_with summarise
#' @importFrom git2rdata is_git2rdata write_vc
#' @importFrom n2kanalysis get_file_fingerprint status
#' @importFrom stats quantile sd
extract_results.n2kHurdleImputed <- function(x, root, ...) {
  assert_that(inherits(x, "n2kHurdleImputed"))
  if (n2kanalysis::status(x) != "converged") {
    return(invisible(NULL))
  }
  x@AnalysisMetadata |>
    select(
      species = "species_group_id",
      "model_type",
      analysis = "file_fingerprint",
      fingerprint = "status_fingerprint"
    ) |>
    bind_cols(
      x@Hurdle@Covariate |>
        transmute(
          .data$province, .data$type, .data$location_id, .data$sublocation_id,
          .data$winter
        ) |>
        bind_cols(x@Hurdle@Imputation) |>
        pivot_longer(
          starts_with("Imputation"),
          names_to = "imputation",
          values_to = "number"
        ) |>
        group_by(
          .data$province, .data$type, .data$location_id, .data$sublocation_id,
          .data$winter
        ) |>
        summarise(
          across(
            "number",
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
        mutate(delta = .data$p95 - .data$p05)
    ) |>
    write_vc(
      file = file.path("hibernation", "hurdle"),
      root = root,
      optimize = FALSE,
      append = TRUE,
      sorting = c(
        "model_type",
        "species",
        "location_id",
        "sublocation_id",
        "winter",
        "analysis"
      ),
      digits = 4
    )
  update_metadata(
    file = file.path("hibernation", "hurdle"),
    root = root,
    name = "hibernation_hurdle",
    title = "The imputed number of hibernating bats",
    description = paste(
      "The imputed total number of hibernating bats in the winter season.",
      "Missing values are imputed before calculating the total.",
      "The model is a first order random walk on the winter season with a",
      "negative binomial distribution."
    ),
    field_description = c(
      species = "The code of the species group",
      model_type = paste(
        "A short description of the model used to calculate the totals"
      ),
      analysis = "The file fingerprint of the analysis",
      fingerprint = "The status fingerprint of the analysis",
      winter = "The winter season defined by the year in which January falls",
      mean = "The average of the imputed total number of hibernating bats",
      min = "The minimum of the imputed total number of hibernating bats",
      p05 = "The 5% quantile of the imputed total number of hibernating bats",
      p20 = paste(
        "The 20% quantile of the imputed total number of hibernating bats"
      ),
      p35 = paste(
        "The 35% quantile of the imputed total number of hibernating bats"
      ),
      p50 = paste(
        "The 50% quantile of the imputed total number of hibernating bats"
      ),
      p65 = paste(
        "The 65% quantile of the imputed total number of hibernating bats"
      ),
      p80 = paste(
        "The 80% quantile of the imputed total number of hibernating bats"
      ),
      p95 = paste(
        "The 95% quantile of the imputed total number of hibernating bats"
      ),
      max = "The maximum of the imputed total number of hibernating bats"
    )
  )
  rm(x)
  gc(verbose = FALSE)
  return(invisible(NULL))
}
