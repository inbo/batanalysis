extract_results <- function(x, ...) {
  UseMethod("extract_results", x)
}

extract_results.default <- function(x, ...) {
  message("No extraction method for class ", class(x))
}

extract_results.character <- function(
    x, base, project = "batanalysis", root, ...
) {
  assert_that(is.string(x), noNA(x))
  manifest <- read_manifest(base = base, project = project, hash = x)
  manifest@Manifest |>
    filter(is.na(.data$parent)) |>
    transmute(.data$fingerprint, level = 0) -> level
  manifest@Manifest |>
    anti_join(level, by = "fingerprint") -> to_do
  while (nrow(to_do)) {
    to_do |>
      inner_join(level, by = c("parent" = "fingerprint")) |>
      distinct(.data$fingerprint, level = .data$level + 1) |>
      bind_rows(level) -> level
    manifest@Manifest |>
      anti_join(level, by = "fingerprint") -> to_do
  }
  rm_data(root = root, path = "hibernation")
  level |>
    filter(level == 2) |>
    slice(1) |>
    pull(fingerprint) |>
    read_model(base = base, project = project) -> x
}

extract_results.n2kModelImputed <- function(x, root, ...) {
  assert_that(inherits(x, "n2kModelImputed"))
  x@AnalysisMetadata |>
    select(
      species = "species_group_id", "model_type", analysis = "file_fingerprint",
      fingerprint = "status_fingerprint"
    ) |>
    bind_cols(x@Results) -> results
  results |>
    filter(str_detect(.data$Parameter, "total")) |>
    transmute(
      .data$species, .data$model_type, .data$analysis, .data$fingerprint,
      winter = str_remove(.data$Parameter, "total: ") |>
        as.integer(),
      estimate = .data$Estimate, se = .data$SE
    ) |>
    write_vc(
      file = "hibernation/total_rw1", root = root, optimize = FALSE,
      append = TRUE,
      sorting = c("model_type", "species", "winter", "analysis")
    )
  update_metadata(
    file = "hibernation/total_rw1", root = root, name = "hibernation_total_rw1",
    title =
      "Total number of hibernating bats modeled with a first order random walk",
    description =
"Model based on the total number of hibernating bats in the winter season.
Missing values are imputed before calculating the total. The model is a first
order random walk on the winter season with a negative binomial distribution.",
    field_description = c(
      species = "The code of the species group",
      model_type =
        "A short description of the model used to analyse the totals",
      analysis = "The file fingerprint of the analysis",
      fingerprint = "The status fingerprint of the analysis",
      winter = "The winter season defined by the year in which January falls",
      estimate =
        "The estimated total number of hibernating bats on the log-scale",
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
      "species", "model_type", "analysis", "fingerprint", "reference", "target",
      estimate = "Estimate", se = "SE"
    ) |>
    mutate(across(c("reference", "target"), as.integer)) |>
    write_vc(
      file = "hibernation/index", root = root, optimize = FALSE, append = TRUE,
      sorting = c("model_type", "species", "reference", "target", "analysis")
    )
  update_metadata(
    file = "hibernation/index", root = root, name = "hibernation_index",
    title =
      "Relative change in total number of hibernating bats between two winters",
    description =
      "Model based on the total number of hibernating bats in the winter season.
Missing values are imputed before calculating the total. The model is a first
order random walk on the winter season with a negative binomial distribution.",
    field_description = c(
      species = "The code of the species group",
      model_type =
        "A short description of the model used to analyse the totals",
      analysis = "The file fingerprint of the analysis",
      fingerprint = "The status fingerprint of the analysis",
      reference =
      "The reference winter season defined by the year in which January falls",
      target =
        "The target winter season defined by the year in which January falls",
      estimate =
"The log-ratio of the total number of hibernating bats in the target winter
divided by those in the reference winter",
      se = "The standard error of the estimate, also on the log-scale"
    )
  )
  results |>
    filter(str_detect(.data$Parameter, "trend")) |>
    separate_wider_regex(
      "Parameter",
      patterns = c(
        "trend_", midpoint = "[0-9\\.]*", "_", duration = "[0-9]*"
      )
    ) |>
    mutate(
      midpoint = as.numeric(.data$midpoint),
      duration = as.integer(.data$duration)
    ) |>
    select(
      "species", "model_type", "analysis", "fingerprint", "midpoint",
      "duration", estimate = "Estimate", se = "SE"
    ) |>
    write_vc(
      file = "hibernation/trend", root = root, optimize = FALSE,
      append = TRUE,
      sorting = c("model_type", "species", "midpoint", "duration", "analysis")
    )
  update_metadata(
    file = "hibernation/trend", root = root, name = "hibernation_trend",
    title =
      "Linear trend in total number of hibernating bats",
    description =
      "Model based on the total number of hibernating bats in the winter season.
Missing values are imputed before calculating the total. The model is a first
order random walk on the winter season with a negative binomial distribution.",
    field_description = c(
      species = "The code of the species group",
      model_type =
        "A short description of the model used to analyse the totals",
      analysis = "The file fingerprint of the analysis",
      fingerprint = "The status fingerprint of the analysis",
      midpoint =
"Central winter of the period over which the trend was calculated. The winter
season defined by the year in which January falls",
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
        "average_", midpoint = "[0-9\\.]*", "_", duration = "[0-9]*"
      )
    ) |>
    mutate(
      midpoint = as.numeric(.data$midpoint),
      duration = as.integer(.data$duration)
    ) |>
    select(
      "species", "model_type", "analysis", "fingerprint", "midpoint",
      "duration", estimate = "Estimate", se = "SE"
    ) |>
    write_vc(
      file = "hibernation/average", root = root, optimize = FALSE,
      append = TRUE,
      sorting = c("model_type", "species", "midpoint", "duration", "analysis")
    )
  update_metadata(
    file = "hibernation/average", root = root, name = "hibernation_average",
    title =
      "The average of total number of hibernating bats over a period",
    description =
      "Model based on the total number of hibernating bats in the winter season.
Missing values are imputed before calculating the total. The model is a first
order random walk on the winter season with a negative binomial distribution.",
    field_description = c(
      species = "The code of the species group",
      model_type =
        "A short description of the model used to analyse the totals",
      analysis = "The file fingerprint of the analysis",
      fingerprint = "The status fingerprint of the analysis",
      midpoint =
"Central winter of the period over which the average was calculated. The winter
season defined by the year in which January falls",
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
        "difference_", reference = "[0-9\\.]*", "_", target = "[0-9\\.]*", "_",
        duration = "[0-9]*"
      )
    ) |>
    mutate(
      across(c("reference", "target"), as.numeric),
      duration = as.integer(.data$duration)
    ) |>
    select(
      "species", "model_type", "analysis", "fingerprint", "reference", "target",
      "duration", estimate = "Estimate", se = "SE"
    ) |>
    write_vc(
      file = "hibernation/difference", root = root, optimize = FALSE,
      append = TRUE,
      sorting = c(
        "model_type", "species", "reference", "target", "duration", "analysis"
      )
    )
  update_metadata(
    file = "hibernation/difference", root = root,
    name = "hibernation_difference",
    title =
"Difference in average number of total number of hibernating bats of two
periods.",
    description =
      "Model based on the total number of hibernating bats in the winter season.
Missing values are imputed before calculating the total. The model is a first
order random walk on the winter season with a negative binomial distribution.",
    field_description = c(
      species = "The code of the species group",
      model_type =
        "A short description of the model used to analyse the totals",
      analysis = "The file fingerprint of the analysis",
      fingerprint = "The status fingerprint of the analysis",
      reference =
"The midpoint of reference winter season. The winter season defined by the year
in which January falls",
      target =
"The midpoint of target winter season. The winter season defined by the year
in which January falls",
      duration = "Number of winters over which the averages were calculated",
      estimate =
"The log-rate of the average of the target period divided by the reference
period",
      se = "The standard error of the estimate, also on the log-scale"
    )
  )
}
