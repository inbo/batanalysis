extract_results <- function(x, ...) {
  UseMethod("extract_results", x)
}

extract_results.default <- function(x, ...) {
  message("No extraction method for class ", class(x))
}

extract_results.character <- function(
  x, base, project = "batanalysis", raw_data, root, ...
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

extract_results.n2kAggregate <- function(x, root, ...) {
  assert_that(inherits(x, "n2kAggregate"))
  x@AggregatedImputed@Covariate |>
    bind_cols(x@AggregatedImputed@Imputation) |>
    pivot_longer(
      starts_with("Imputation"), names_to = "imputation", values_to = "total"
    ) -> results
  if (has_name(results, "count_location")) {
    x@AnalysisMetadata |>
      select(
        species = "species_group_id", "model_type",
        analysis = "file_fingerprint", fingerprint = "status_fingerprint"
      ) |>
      bind_cols(
        results |>
          group_by(
            winter = .data$count_winter, location = .data$count_location
          ) |>
          summarise(
            mean = mean(.data$total, na.rm = TRUE),
            min = min(.data$total, na.rm = TRUE),
            p05 = quantile(.data$total, 0.05, na.rm = TRUE),
            p20 = quantile(.data$total, 0.2, na.rm = TRUE),
            p35 = quantile(.data$total, 0.35, na.rm = TRUE),
            p50 = quantile(.data$total, 0.5, na.rm = TRUE),
            p65 = quantile(.data$total, 0.65, na.rm = TRUE),
            p80 = quantile(.data$total, 0.8, na.rm = TRUE),
            p95 = quantile(.data$total, 0.95, na.rm = TRUE),
            max = max(.data$total, na.rm = TRUE), .groups = "drop"
          ) |>
          mutate(location = as.integer(levels(.data$location))[.data$location])
      ) |>
      write_vc(
        file = "hibernation/total_location", root = root, optimize = FALSE,
        append = TRUE,
        sorting = c("model_type", "species", "winter", "location", "analysis")
      )
    update_metadata(
      file = "hibernation/total_location", root = root,
      name = "hibernation_total_location",
      title = "The imputed total number of hibernating bats per location",
      description =
"The imputed total number of hibernating bats in the winter season per locaion.
Missing values are imputed before calculating the total. The model is a first
order random walk on the winter season with a negative binomial distribution.",
      field_description = c(
        species = "The code of the species group",
        model_type =
          "A short description of the model used to calculate the totals",
        analysis = "The file fingerprint of the analysis",
        fingerprint = "The status fingerprint of the analysis",
        winter = "The winter season defined by the year in which January falls",
        mean = "The average of the imputed total number of hibernating bats",
        min = "The minimum of the imputed total number of hibernating bats",
        p05 = "The 5% quantile of the imputed total number of hibernating bats",
      p20 = "The 20% quantile of the imputed total number of hibernating bats",
      p35 = "The 35% quantile of the imputed total number of hibernating bats",
      p50 = "The 50% quantile of the imputed total number of hibernating bats",
      p65 = "The 65% quantile of the imputed total number of hibernating bats",
      p80 = "The 80% quantile of the imputed total number of hibernating bats",
      p95 = "The 95% quantile of the imputed total number of hibernating bats",
        max = "The maximum of the imputed total number of hibernating bats"
      )
    )
  } else {
    x@AnalysisMetadata |>
      select(
        species = "species_group_id", "model_type",
        analysis = "file_fingerprint", fingerprint = "status_fingerprint"
      ) |>
      bind_cols(
        results |>
          group_by(winter = .data$count_winter) |>
          summarise(
            mean = mean(.data$total, na.rm = TRUE),
            min = min(.data$total, na.rm = TRUE),
            p05 = quantile(.data$total, 0.05, na.rm = TRUE),
            p20 = quantile(.data$total, 0.2, na.rm = TRUE),
            p35 = quantile(.data$total, 0.35, na.rm = TRUE),
            p50 = quantile(.data$total, 0.5, na.rm = TRUE),
            p65 = quantile(.data$total, 0.65, na.rm = TRUE),
            p80 = quantile(.data$total, 0.8, na.rm = TRUE),
            p95 = quantile(.data$total, 0.95, na.rm = TRUE),
            max = max(.data$total, na.rm = TRUE)
          )
      ) |>
      write_vc(
        file = "hibernation/total", root = root, optimize = FALSE,
        append = TRUE,
        sorting = c("model_type", "species", "winter", "analysis")
      )
    update_metadata(
      file = "hibernation/total", root = root, name = "hibernation_total",
      title = "The imputed total number of hibernating bats",
      description =
"The imputed total number of hibernating bats in the winter season.
Missing values are imputed before calculating the total. The model is a first
order random walk on the winter season with a negative binomial distribution.",
      field_description = c(
        species = "The code of the species group",
        model_type =
          "A short description of the model used to calculate the totals",
        analysis = "The file fingerprint of the analysis",
        fingerprint = "The status fingerprint of the analysis",
        winter = "The winter season defined by the year in which January falls",
        mean = "The average of the imputed total number of hibernating bats",
        min = "The minimum of the imputed total number of hibernating bats",
      p05 = "The 5% quantile of the imputed total number of hibernating bats",
      p20 = "The 20% quantile of the imputed total number of hibernating bats",
      p35 = "The 35% quantile of the imputed total number of hibernating bats",
      p50 = "The 50% quantile of the imputed total number of hibernating bats",
      p65 = "The 65% quantile of the imputed total number of hibernating bats",
      p80 = "The 80% quantile of the imputed total number of hibernating bats",
      p95 = "The 95% quantile of the imputed total number of hibernating bats",
        max = "The maximum of the imputed total number of hibernating bats"
      )
    )
  }
}
