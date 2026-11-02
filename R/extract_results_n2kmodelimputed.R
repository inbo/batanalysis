#' @export
#' @importFrom assertthat assert_that
#' @importFrom dplyr across bind_cols filter mutate select transmute
#' @importFrom git2rdata is_git2rdata update_metadata write_vc
#' @importFrom n2kanalysis get_file_fingerprint status
#' @importFrom stringr str_detect str_remove
#' @importFrom tidyr separate_wider_regex
extract_results.n2kModelImputed <- function(x, root, ...) {
  assert_that(inherits(x, "n2kModelImputed"))
  if (status(x) != "converged") {
    return(invisible(NULL))
  }
  x@AnalysisMetadata |>
    select(
      species = "species_group_id",
      "model_type",
      analysis = "file_fingerprint",
      fingerprint = "status_fingerprint"
    ) |>
    bind_cols(x@Results) -> results
  results |>
    filter(str_detect(.data$Parameter, "total")) |>
    transmute(
      .data$species,
      .data$model_type,
      .data$analysis,
      .data$fingerprint,
      winter = str_remove(.data$Parameter, "total: ") |>
        as.integer(),
      estimate = .data$Estimate,
      se = .data$SE
    ) |>
    write_vc(
      file = file.path("hibernation", "total_rw1"),
      root = root,
      optimize = FALSE,
      append = TRUE,
      sorting = c("model_type", "species", "winter", "analysis"),
      digits = 4
    )
  update_metadata(
    file = file.path("hibernation", "total_rw1"),
    root = root,
    name = "hibernation_total_rw1",
    title = paste(
      "Total number of hibernating bats modeled with a first order random walk"
    ),
    description = paste(
      "Model based on the total number of hibernating bats in the winter",
      "season.",
      "Missing values are imputed before calculating the total.",
      "The model is a first order random walk on the winter season with a",
      "negative binomial distribution."
    ),
    field_description = c(
      species = "The code of the species group",
      model_type = paste(
        "A short description of the model used to analyse the totals"
      ),
      analysis = "The file fingerprint of the analysis",
      fingerprint = "The status fingerprint of the analysis",
      winter = "The winter season defined by the year in which January falls",
      estimate = paste(
        "The estimated total number of hibernating bats on the log-scale"
      ),
      se = "The standard error of the estimate, also on the log-scale"
    )
  )
  results |>
    filter(str_detect(.data$Parameter, "index")) |>
    separate_wider_regex(
      "Parameter",
      patterns = c("index: ", target = "[0-9]*", "-", reference = "[0-9]*")
    ) |>
    select(
      "species",
      "model_type",
      "analysis",
      "fingerprint",
      "reference",
      "target",
      estimate = "Estimate",
      se = "SE"
    ) |>
    mutate(across(c("reference", "target"), as.integer)) |>
    write_vc(
      file = file.path("hibernation", "index"),
      root = root,
      optimize = FALSE,
      append = TRUE,
      sorting = c("model_type", "species", "reference", "target", "analysis"),
      digits = 4
    )
  update_metadata(
    file = file.path("hibernation", "index"),
    root = root,
    name = "hibernation_index",
    title = paste(
      "Relative change in total number of hibernating bats between two winters"
    ),
    description = paste(
      "Model based on the total number of hibernating bats in the winter",
      "season.",
      "Missing values are imputed before calculating the total.",
      "The model is a first order random walk on the winter season with a",
      "negative binomial distribution."
    ),
    field_description = c(
      species = "The code of the species group",
      model_type = paste(
        "A short description of the model used to analyse the totals"
      ),
      analysis = "The file fingerprint of the analysis",
      fingerprint = "The status fingerprint of the analysis",
      reference = paste(
        "The reference winter season defined by the year in which January falls"
      ),
      target = paste(
        "The target winter season defined by the year in which January falls"
      ),
      estimate = paste(
        "The log-ratio of the total number of hibernating bats in the target",
        "winter divided by those in the reference winter."
      ),
      se = "The standard error of the estimate, also on the log-scale"
    )
  )
  results |>
    filter(str_detect(.data$Parameter, "trend")) |>
    separate_wider_regex(
      "Parameter",
      patterns = c(
        "trend_",
        midpoint = "[0-9\\.]*",
        "_",
        duration = "[0-9]*"
      )
    ) |>
    mutate(
      midpoint = as.numeric(.data$midpoint),
      duration = as.integer(.data$duration)
    ) |>
    select(
      "species",
      "model_type",
      "analysis",
      "fingerprint",
      "midpoint",
      "duration",
      estimate = "Estimate",
      se = "SE"
    ) |>
    write_vc(
      file = file.path("hibernation", "trend"),
      root = root,
      optimize = FALSE,
      append = TRUE,
      sorting = c("model_type", "species", "midpoint", "duration", "analysis"),
      digits = 4
    )
  update_metadata(
    file = file.path("hibernation", "trend"),
    root = root,
    name = "hibernation_trend",
    title = "Linear trend in total number of hibernating bats",
    description = paste(
      "Model based on the total number of hibernating bats in the winter",
      "season.",
      "Missing values are imputed before calculating the total.",
      "The model is a first order random walk on the winter season with a",
      "negative binomial distribution."
    ),
    field_description = c(
      species = "The code of the species group",
      model_type = paste(
        "A short description of the model used to analyse the totals"
      ),
      analysis = "The file fingerprint of the analysis",
      fingerprint = "The status fingerprint of the analysis",
      midpoint = paste(
        "Central winter of the period over which the trend was calculated.",
        "The winter season defined by the year in which January falls"
      ),
      duration = "Number of winters over which the trend was calculated",
      estimate = "The average change per winter on the log-scale",
      se = "The standard error of the estimate, also on the log-scale"
    )
  )
  results |>
    filter(str_detect(.data$Parameter, "average")) |>
    separate_wider_regex(
      "Parameter",
      patterns = c(
        "average_",
        midpoint = "[0-9\\.]*",
        "_",
        duration = "[0-9]*"
      )
    ) |>
    mutate(
      midpoint = as.numeric(.data$midpoint),
      duration = as.integer(.data$duration)
    ) |>
    select(
      "species",
      "model_type",
      "analysis",
      "fingerprint",
      "midpoint",
      "duration",
      estimate = "Estimate",
      se = "SE"
    ) |>
    write_vc(
      file = file.path("hibernation", "average"),
      root = root,
      optimize = FALSE,
      append = TRUE,
      sorting = c("model_type", "species", "midpoint", "duration", "analysis"),
      digits = 4
    )
  update_metadata(
    file = file.path("hibernation", "average"),
    root = root,
    name = "hibernation_average",
    title = "The average of total number of hibernating bats over a period",
    description = paste(
      "Model based on the total number of hibernating bats in the winter",
      "season.",
      "Missing values are imputed before calculating the total.",
      "The model is a first order random walk on the winter season with a",
      "negative binomial distribution."
    ),
    field_description = c(
      species = "The code of the species group",
      model_type = paste(
        "A short description of the model used to analyse the totals"
      ),
      analysis = "The file fingerprint of the analysis",
      fingerprint = "The status fingerprint of the analysis",
      midpoint = paste(
        "Central winter of the period over which the average was calculated.",
        "The winter season defined by the year in which January falls"
      ),
      duration = "Number of winters over which the average was calculated",
      estimate = "The average on the log-scale",
      se = "The standard error of the estimate, also on the log-scale"
    )
  )
  results |>
    filter(str_detect(.data$Parameter, "difference")) |>
    separate_wider_regex(
      "Parameter",
      patterns = c(
        "difference_",
        reference = "[0-9\\.]*",
        "_",
        target = "[0-9\\.]*",
        "_",
        duration = "[0-9]*"
      )
    ) |>
    mutate(
      across(c("reference", "target"), as.numeric),
      duration = as.integer(.data$duration)
    ) |>
    select(
      "species",
      "model_type",
      "analysis",
      "fingerprint",
      "reference",
      "target",
      "duration",
      estimate = "Estimate",
      se = "SE"
    ) |>
    write_vc(
      file = file.path("hibernation", "difference"),
      root = root,
      optimize = FALSE,
      append = TRUE,
      sorting = c(
        "model_type",
        "species",
        "reference",
        "target",
        "duration",
        "analysis"
      ),
      digits = 4
    )
  update_metadata(
    file = file.path("hibernation", "difference"),
    root = root,
    name = "hibernation_difference",
    title = paste(
      "Difference in average number of total number of hibernating bats of",
      "two periods."
    ),
    description = paste(
      "Model based on the total number of hibernating bats in the winter",
      "season.",
      "Missing values are imputed before calculating the total.",
      "The model is a first order random walk on the winter season with a",
      "negative binomial distribution."
    ),
    field_description = c(
      species = "The code of the species group",
      model_type = paste(
        "A short description of the model used to analyse the totals"
      ),
      analysis = "The file fingerprint of the analysis",
      fingerprint = "The status fingerprint of the analysis",
      reference = paste(
        "The midpoint of reference winter season.",
        "The winter season defined by the year in which January falls"
      ),
      target = paste(
        "The midpoint of target winter season.",
        "The winter season defined by the year in which January falls"
      ),
      duration = "Number of winters over which the averages were calculated",
      estimate = paste(
        "The log-rate of the average of the target period divided by the",
        "reference period"
      ),
      se = "The standard error of the estimate, also on the log-scale"
    )
  )
  rm(x, results)
  gc(verbose = FALSE)
  return(invisible(NULL))
}
