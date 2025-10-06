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

#' @export
#' @importFrom assertthat assert_that is.flag is.string noNA
#' @importFrom dplyr anti_join distinct
#' @importFrom git2rdata verify_vc write_vc
#' @importFrom n2kanalysis order_manifest read_manifest read_model
#' @importFrom purrr walk
extract_results.character <- function(
  x,
  base,
  project = "batanalysis",
  raw_data,
  root,
  random = FALSE,
  ...
) {
  assert_that(is.string(x), noNA(x), is.flag(random), noNA(random))
  file.path("data", "hibernation", "species") |>
    verify_vc(
      root = raw_data,
      variables = c("id", "code", "name", "scientific_name", "parent")
    ) |>
    write_vc(
      file.path("hibernation", "species"),
      root = root,
      sorting = "id",
      optimize = FALSE
    )
  update_metadata(
    file.path("hibernation", "species"),
    root = root,
    name = "hibernation_species",
    title = "Hibernating bat species",
    description = paste(
      "List of species potentially relevant for the hibernating bat",
      "monitoring in Flanders (Belgium)."
    ),
    field_description = c(
      id = "Unique identifier of the species",
      name = "Dutch vernacular name",
      scientific_name = "Scientific name of the species",
      code = "Code of the species",
      parent = "Unique identifier of the parent species"
    )
  )
  file.path("data", "hibernation", "locations") |>
    verify_vc(
      root = raw_data,
      variables = c("id", "code", "name", "parent_id")
    ) |>
    write_vc(
      file.path("hibernation", "locations"),
      root = root,
      sorting = "id",
      optimize = FALSE,
      digits = c("longitude" = 8, "latitude" = 8)
    )
  update_metadata(
    file.path("hibernation", "locations"),
    root = root,
    name = "hibernation_locations",
    title = paste(
      "Locations and sublocations surveyed during the hibernating bat",
      "monitoring"
    ),
    description = paste(
      "List of locations and sublocations surveyed during the hibernating",
      "bat monitoring in Flanders (Belgium)."
    ),
    field_description = c(
      id = "Unique identifier of the location or sublocation",
      name = "Name of the location or sublocation",
      parent_id = "Unique identifier of the parent location (-1 for locations)",
      code = "Code of the location or sublocation",
      longitude = "Longitude of the location or sublocation (EPSG:4326)",
      latitude = "Latitude of the location or sublocation (EPSG:4326)"
    )
  )
  read_manifest(base = base, project = project, hash = x) |>
    order_manifest() -> manifest
  for (type in c(
    "hurdle",
    "total",
    "total_location",
    "total_type",
    "total_rw1"
  )) {
    filename <- file.path("hibernation", type)
    if (is_git2rdata(filename, root = root)) {
      manifest <- manifest[
        !manifest %in%
          verify_vc(filename, root = root, "analysis")$analysis
      ]
    }
  }
  filename <- file.path("model_check", "hibernation", "hyperparameters")
  if (is_git2rdata(filename, root = root)) {
    manifest <- manifest[
      !manifest %in%
        verify_vc(filename, root = root, "analysis")$analysis
    ]
  }
  if (random) {
    manifest <- sample(manifest)
  } else {
    file.path(base, project) |>
      list.files(pattern = ".rds$", full.names = TRUE, recursive = TRUE) |>
      grepv(pattern = "converged") |>
      file.info() |>
      arrange(.data$size) |>
      rownames() |>
      basename() |>
      gsub(pattern = ".rds$", replacement = "") -> to_do
    to_do[to_do %in% manifest] -> manifest
  }
  for (i in manifest) {
    message(i)
    read_model(i, base = base, project = project) |>
      extract_results(root = root)
    gc(verbose = FALSE)
  }
  return(invisible(NULL))
}

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

#' @export
#' @importFrom assertthat assert_that
#' @importFrom dplyr bind_cols case_when group_by inner_join mutate select
#' starts_with summarise
#' @importFrom git2rdata is_git2rdata update_metadata write_vc
#' @importFrom n2kanalysis get_file_fingerprint status
#' @importFrom rlang .data
#' @importFrom stats quantile
#' @importFrom tidyr pivot_longer
extract_results.n2kAggregate <- function(x, root, ...) {
  assert_that(inherits(x, "n2kAggregate"))
  if (status(x) != "converged") {
    return(invisible(NULL))
  }
  x@AggregatedImputed@Covariate |>
    bind_cols(x@AggregatedImputed@Imputation) |>
    pivot_longer(
      starts_with("Imputation"),
      names_to = "imputation",
      values_to = "total"
    ) -> results
  if (has_name(results, "location_id")) {
    results |>
      group_by(.data$winter, .data$imputation) |>
      summarise(
        flanders = sum(.data$total, na.rm = TRUE),
        .groups = "drop"
      ) |>
      inner_join(results, by = c("winter", "imputation")) |>
      mutate(fraction = .data$total / .data$flanders) -> combined
    x@AnalysisMetadata |>
      select(
        species = "species_group_id",
        "model_type",
        analysis = "file_fingerprint",
        fingerprint = "status_fingerprint"
      ) |>
      bind_cols(
        combined |>
          group_by(.data$winter, .data$imputation) |>
          mutate(rank = rank(-.data$total, ties.method = "first")) |>
          group_by(.data$winter, .data$location_id) |>
          summarise(
            mean = mean(.data$rank, na.rm = TRUE),
            min = min(.data$rank, na.rm = TRUE),
            p05 = quantile(.data$rank, 0.05, na.rm = TRUE),
            p20 = quantile(.data$rank, 0.2, na.rm = TRUE),
            p35 = quantile(.data$rank, 0.35, na.rm = TRUE),
            p50 = quantile(.data$rank, 0.5, na.rm = TRUE),
            p65 = quantile(.data$rank, 0.65, na.rm = TRUE),
            p80 = quantile(.data$rank, 0.8, na.rm = TRUE),
            p95 = quantile(.data$rank, 0.95, na.rm = TRUE),
            max = max(.data$rank, na.rm = TRUE),
            .groups = "drop"
          )
      ) |>
      write_vc(
        file = file.path("hibernation", "rank_location"),
        root = root,
        optimize = FALSE,
        append = TRUE,
        sorting = c(
          "model_type",
          "species",
          "winter",
          "location_id",
          "analysis"
        ),
        digits = 6
      )
    update_metadata(
      file = file.path("hibernation", "rank_location"),
      root = root,
      name = "hibernation_rank_location",
      title = paste(
        "Rank of the location in Flanders"
      ),
      description = paste(
        "The rank of the location in Flanders",
        "Missing values are imputed before calculating the rank."
      ),
      field_description = c(
        species = "The code of the species group",
        model_type = paste(
          "A short description of the model used to calculate the totals"
        ),
        analysis = "The file fingerprint of the analysis",
        fingerprint = "The status fingerprint of the analysis",
        winter = "The winter season defined by the year in which January falls",
        location_id = "The unique identifier of the location",
        mean = "The average of rank",
        min = "The minimum rank",
        p05 = "The 5% quantile of the rank",
        p20 = "The 20% quantile of the rank",
        p35 = "The 35% quantile of the rank",
        p50 = "The 50% quantile of the rank",
        p65 = "The 65% quantile of the rank",
        p80 = "The 80% quantile of the rank",
        p95 = "The 95% quantile of the rank",
        max = "The maximum of the rank"
      )
    )
    x@AnalysisMetadata |>
      select(
        species = "species_group_id",
        "model_type",
        analysis = "file_fingerprint",
        fingerprint = "status_fingerprint"
      ) |>
      bind_cols(
        combined |>
          group_by(.data$winter, .data$location_id) |>
          summarise(
            mean = mean(.data$fraction, na.rm = TRUE),
            min = min(.data$fraction, na.rm = TRUE),
            p05 = quantile(.data$fraction, 0.05, na.rm = TRUE),
            p20 = quantile(.data$fraction, 0.2, na.rm = TRUE),
            p35 = quantile(.data$fraction, 0.35, na.rm = TRUE),
            p50 = quantile(.data$fraction, 0.5, na.rm = TRUE),
            p65 = quantile(.data$fraction, 0.65, na.rm = TRUE),
            p80 = quantile(.data$fraction, 0.8, na.rm = TRUE),
            p95 = quantile(.data$fraction, 0.95, na.rm = TRUE),
            max = max(.data$total, na.rm = TRUE),
            .groups = "drop"
          )
      ) |>
      write_vc(
        file = file.path("hibernation", "fraction_location"),
        root = root,
        optimize = FALSE,
        append = TRUE,
        sorting = c(
          "model_type",
          "species",
          "winter",
          "location_id",
          "analysis"
        ),
        digits = 6
      )
    update_metadata(
      file = file.path("hibernation", "fraction_location"),
      root = root,
      name = "hibernation_fraction_location",
      title = paste(
        "Fraction of the total number per location and the total number in",
        "Flanders"
      ),
      description = paste(
        "The imputed total number of hibernating bats in the winter season per",
        "location devided by the total number of hibernating bats in Flanders.",
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
        location_id = "The unique identifier of the location",
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

    x@AnalysisMetadata |>
      select(
        species = "species_group_id",
        "model_type",
        analysis = "file_fingerprint",
        fingerprint = "status_fingerprint"
      ) |>
      bind_cols(
        results |>
          group_by(.data$winter, .data$location_id) |>
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
            max = max(.data$total, na.rm = TRUE),
            .groups = "drop"
          )
      ) |>
      write_vc(
        file = file.path("hibernation", "total_location"),
        root = root,
        optimize = FALSE,
        append = TRUE,
        sorting = c(
          "model_type",
          "species",
          "winter",
          "location_id",
          "analysis"
        ),
        digits = 6
      )
    update_metadata(
      file = file.path("hibernation", "total_location"),
      root = root,
      name = "hibernation_total_location",
      title = "The imputed total number of hibernating bats per location",
      description = paste(
        "The imputed total number of hibernating bats in the winter season per",
        "location.",
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
        location_id = "The unique identifier of the location",
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
    rm(x, results)
    gc(verbose = FALSE)
    return(invisible(NULL))
  }
  if (
    any(
      c("fortress", "marl_quarry", "other_large", "small") %in%
        colnames(results)
    )
  ) {
    x@AnalysisMetadata |>
      select(
        species = "species_group_id",
        "model_type",
        analysis = "file_fingerprint",
        fingerprint = "status_fingerprint"
      ) |>
      bind_cols(
        results |>
          group_by(
            .data$winter,
            type = case_when(
              .data$fortress > 0 ~ "fortress",
              .data$marl_quarry > 0 ~ "marl quarry",
              .data$other_large > 0 ~ "other large",
              .data$small > 0 ~ "small"
            )
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
            max = max(.data$total, na.rm = TRUE),
            .groups = "drop"
          ) |>
          mutate(
            type = factor(
              .data$type,
              c("fortress", "marl quarry", "other large", "small")
            )
          )
      ) |>
      write_vc(
        file = file.path("hibernation", "total_type"),
        root = root,
        optimize = FALSE,
        append = TRUE,
        sorting = c("model_type", "species", "winter", "type", "analysis"),
        digits = 6
      )
    update_metadata(
      file = file.path("hibernation", "total_type"),
      root = root,
      name = "hibernation_total_location",
      title = paste(
        "The imputed total number of hibernating bats per type of location"
      ),
      description = paste(
        "The imputed total number of hibernating bats in the winter season per",
        "type of location.",
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
        type = "The type of location",
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
    rm(x, results)
    gc(verbose = FALSE)
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
      results |>
        group_by(winter = .data$winter) |>
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
      file = file.path("hibernation", "total"),
      root = root,
      optimize = FALSE,
      append = TRUE,
      sorting = c("model_type", "species", "winter", "analysis"),
      digits = 6
    )
  update_metadata(
    file = file.path("hibernation", "total"),
    root = root,
    name = "hibernation_total",
    title = "The imputed total number of hibernating bats",
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
  rm(x, results)
  gc(verbose = FALSE)
  return(invisible(NULL))
}

#' @export
#' @importFrom assertthat assert_that
#' @importFrom dplyr arrange bind_cols distinct group_by inner_join mutate
#' row_number select starts_with summarise
#' @importFrom git2rdata is_git2rdata update_metadata write_vc
#' @importFrom INLA inla.mesh.projector inla.posterior.sample inla.zmarginal
#' @importFrom inlatools prec2sd
#' @importFrom n2kanalysis get_file_fingerprint spde2mesh status
#' @importFrom purrr map map2_dfc
#' @importFrom sf st_area st_as_sf st_convex_hull st_sample st_union
#' @importFrom stringr str_detect str_remove str_replace_all
#' @importFrom stats quantile
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr everything pivot_longer pivot_wider unnest
extract_results.n2kSpde <- function(x, root, ..., n_sim = 100) {
  assert_that(inherits(x, "n2kSpde"), is.count(n_sim), noNA(n_sim))
  # skip if model is not converged
  if (status(x) != "converged") {
    return(invisible(NULL))
  }

  # generate posterior sample
  inla.posterior.sample(n = n_sim, result = x@Model) |>
    map("latent") |>
    map2_dfc(
      seq_len(n_sim),
      ~ sprintf("sim_%04i", .y) |>
        `colnames<-`(.x, value = _) |>
        as.data.frame()
    ) |>
    rownames_to_column(var = "parameter") |>
    pivot_longer(
      -"parameter",
      names_to = "sim",
      values_to = "estimate"
    ) -> post_sample
  x@Data |>
    distinct(
      .data$location_id,
      .data$fortress,
      .data$marl_quarry,
      .data$other_large,
      .data$small
    ) |>
    transmute(
      .data$location_id,
      type = case_when(
        .data$fortress > 0 ~ "fortress",
        .data$marl_quarry > 0 ~ "marl quarry",
        .data$other_large > 0 ~ "other large",
        .data$small > 0 ~ "small"
      ) |>
        factor(c("fortress", "marl quarry", "other large", "small"))
    ) -> location_type
  # extract the intercept and the winter effect
  post_sample |>
    filter(
      .data$parameter %in%
        c("fortress:1", "marl_quarry:1", "other_large:1", "small:1")
    ) |>
    transmute(
      type = str_remove(.data$parameter, ":1") |>
        str_replace_all("_", " ") |>
        factor(levels = levels(location_type$type)),
      .data$sim,
      intercept = .data$estimate
    ) -> ps_intercept
  post_sample |>
    filter(str_detect(.data$parameter, "winter_i")) |>
    transmute(
      winter_r = str_remove(.data$parameter, "winter_i:") |>
        as.integer(),
      .data$sim,
      .data$estimate
    ) |>
    inner_join(ps_intercept, by = "sim", relationship = "many-to-many") |>
    transmute(
      winter = .data$winter_r + min(x@Data$winter) - 1,
      .data$type,
      .data$sim,
      estimate = .data$intercept + .data$estimate
    ) -> ps_cwinter
  x@AnalysisMetadata |>
    select(
      species = "species_group_id",
      "model_type",
      analysis = "file_fingerprint",
      fingerprint = "status_fingerprint"
    ) |>
    bind_cols(
      ps_cwinter |>
        group_by(.data$type, .data$winter) |>
        summarise(
          across(
            "estimate",
            .names = "{.fn}",
            list(
              mean = ~ mean(.x, na.rm = TRUE) |>
                round(4),
              sd = ~ sd(.x, na.rm = TRUE) |>
                round(4)
            )
          ),
          .groups = "drop"
        )
    ) |>
    write_vc(
      file = file.path("model_check", "hibernation", "winter"),
      root = root,
      append = TRUE,
      sorting = c("model_type", "species", "winter", "type", "analysis"),
      digits = 4
    )
  update_metadata(
    file = file.path("model_check", "hibernation", "winter"),
    root = root,
    name = "hibernation_model_check_winter",
    title = "Modelled winter effect",
    description = paste(
      "Posterior distribution of the modelled winter effect at the location",
      "level.",
      "The model is a binomial spatio-temporal model with a logit link.",
      "The winter effect is modeled with a second order random walk."
    ),
    field_description = c(
      species = "The code of the species group",
      model_type = paste(
        "A short description of the model used to analyse the data"
      ),
      analysis = "The file fingerprint of the analysis",
      fingerprint = "The status fingerprint of the analysis",
      winter = "The winter season defined by the year in which January falls",
      type = "The type of location",
      mean = paste(
        "The mean of the posterior distribution on the logit scale for",
        "binomial models and the log scale for nbinomial models"
      ),
      sd = paste(
        "The standard deviation of the posterior distribution on the logit",
        "scale for binomial models and the log scale for nbinomial models"
      )
    )
  )
  # extract the location specific effects
  post_sample |>
    filter(str_detect(.data$parameter, "^matern_spde:")) |>
    mutate(
      parameter = str_remove(.data$parameter, "matern_spde:") |>
        as.integer()
    ) |>
    pivot_wider(names_from = "sim", values_from = "estimate") |>
    arrange(.data$parameter) |>
    select(-"parameter") |>
    as.matrix() -> ps_matern
  x@Model$.args$data[c("location_id", "X", "Y")] |>
    as.data.frame() |>
    distinct() |>
    filter(!is.na(.data$X)) -> loc_coordinates
  projector <- inla.mesh.projector(
    mesh = spde2mesh(x@Spde),
    loc = as.matrix(loc_coordinates[, c("X", "Y")])
  )
  x@Data |>
    distinct(
      .data$location_id,
      location_i = as.integer(.data$location_i),
      .data$winter,
      .data$winter_l,
      .data$winter_q,
      .data$winter_c
    ) -> location_winter
  post_sample |>
    filter(str_detect(.data$parameter, "^location_i:")) |>
    transmute(
      location_i = str_remove(.data$parameter, "location_i:") |>
        as.integer(),
      .data$sim,
      q0 = .data$estimate
    ) |>
    left_join(
      post_sample |>
        filter(str_detect(.data$parameter, "^location_l:")) |>
        transmute(
          location_i = str_remove(.data$parameter, "location_l:") |>
            as.integer(),
          .data$sim,
          q1 = .data$estimate
        ),
      by = c("location_i", "sim")
    ) |>
    left_join(
      post_sample |>
        filter(str_detect(.data$parameter, "^location_q:")) |>
        transmute(
          location_i = str_remove(.data$parameter, "location_q:") |>
            as.integer(),
          .data$sim,
          q2 = .data$estimate
        ),
      by = c("location_i", "sim")
    ) |>
    left_join(
      post_sample |>
        filter(str_detect(.data$parameter, "^location_c:")) |>
        transmute(
          location_i = str_remove(.data$parameter, "location_c:") |>
            as.integer(),
          .data$sim,
          q3 = .data$estimate
        ),
      by = c("location_i", "sim")
    ) |>
    inner_join(
      location_winter |>
        group_by(.data$location_i, .data$location_id) |>
        summarise(
          across(
            c("winter_l", "winter_q", "winter_c"),
            ~ !any(is.na(.x))
          ),
          .groups = "drop"
        ),
      by = "location_i"
    ) |>
    transmute(
      .data$location_id,
      .data$sim,
      .data$q0,
      q1 = ifelse(.data$winter_l, .data$q1, 0),
      q2 = ifelse(.data$winter_q, .data$q2, 0),
      q3 = ifelse(.data$winter_c, .data$q3, 0),
      across(c("q0", "q1", "q2", "q3"), ~ replace_na(.x, 0))
    ) |>
    inner_join(
      as.matrix(projector$proj$A %*% ps_matern) |>
        as.data.frame() |>
        mutate(
          location_id = loc_coordinates$location_id
        ) |>
        pivot_longer(starts_with("sim"), names_to = "sim", values_to = "mesh"),
      by = c("location_id", "sim")
    ) -> ps_location

  x@AnalysisMetadata |>
    select(
      species = "species_group_id",
      "model_type",
      analysis = "file_fingerprint",
      fingerprint = "status_fingerprint"
    ) |>
    bind_cols(
      ps_location |>
        group_by(.data$location_id) |>
        summarise(
          across(
            c("q0", "q1", "q2", "q3", "mesh"),
            list(
              mean = ~ mean(.x, na.rm = TRUE),
              sd = ~ sd(.x, na.rm = TRUE)
            )
          )
        )
    ) |>
    write_vc(
      file = file.path("model_check", "hibernation", "location_effect"),
      root = root,
      append = TRUE,
      sorting = c("model_type", "species", "location_id", "analysis"),
      digits = 4
    )
  update_metadata(
    file = file.path("model_check", "hibernation", "location_effect"),
    root = root,
    name = "hibernation_model_check_location_effect",
    title = "Location specific effects",
    description = paste(
      "Posterior distribution of the location specific effects.",
      "The effects are on the logit scale for binomial models and the log scale",
      "for nbinomial models."
    ),
    field_description = c(
      species = "The code of the species group",
      model_type = paste(
        "A short description of the model used to analyse the data"
      ),
      analysis = "The file fingerprint of the analysis",
      fingerprint = "The status fingerprint of the analysis",
      location_id = "The unique identifier of the location",
      q0_mean = paste(
        "The mean of the intercept of the location specific effect on the",
        "logit scale for binomial models and the log scale for nbinomial models"
      ),
      q0_sd = paste(
        "The standard deviation of the intercept of the location specific",
        "effect on the logit scale for binomial models and the log scale for",
        "nbinomial models"
      ),
      q1_mean = paste(
        "The mean of the linear term of the location specific effect on the",
        "logit scale for binomial models and the log scale for nbinomial models"
      ),
      q1_sd = paste(
        "The standard deviation of the linear term of the location specific",
        "effect on the logit scale for binomial models and the log scale for",
        "nbinomial models"
      ),
      q2_mean = paste(
        "The mean of the quadratic term of the location specific effect on",
        "the logit scale for binomial models and the log scale for nbinomial",
        "models"
      ),
      q2_sd = paste(
        "The standard deviation of the quadratic term of the location",
        "specific effect on the logit scale for binomial models and the log",
        "scale for nbinomial models"
      ),
      q3_mean = paste(
        "The mean of the cubic term of the location specific effect on the",
        "logit scale for binomial models and the log scale for nbinomial models"
      ),
      q3_sd = paste(
        "The standard deviation of the cubic term of the location",
        "specific effect on the logit scale for binomial models and the log",
        "scale for nbinomial models"
      )
    )
  )
  # predictions at the location level
  location_winter |>
    mutate(
      across(c("winter_l", "winter_q", "winter_c"), ~ replace_na(.x, 0))
    ) |>
    inner_join(x = location_type, by = "location_id") |>
    inner_join(
      ps_location,
      by = "location_id",
      relationship = "many-to-many"
    ) |>
    inner_join(ps_cwinter, by = c("winter", "type", "sim")) |>
    transmute(
      .data$location_id,
      .data$winter,
      .data$sim,
      relative = .data$q0 +
        .data$q1 * .data$winter_l +
        .data$q2 * .data$winter_q +
        .data$q3 * .data$winter_c +
        .data$mesh,
      absolute = .data$estimate + .data$relative
    ) -> ps_location
  x@AnalysisMetadata |>
    select(
      species = "species_group_id",
      "model_type",
      analysis = "file_fingerprint",
      fingerprint = "status_fingerprint"
    ) |>
    bind_cols(
      ps_location |>
        group_by(.data$location_id, .data$winter) |>
        summarise(
          across(
            c("relative", "absolute"),
            list(
              mean = ~ mean(.x, na.rm = TRUE),
              sd = ~ sd(.x, na.rm = TRUE)
            )
          ),
          .groups = "drop"
        )
    ) |>
    write_vc(
      file = file.path("model_check", "hibernation", "location"),
      root = root,
      append = TRUE,
      sorting = c("model_type", "species", "location_id", "winter", "analysis"),
      digits = 4
    )
  update_metadata(
    file = file.path("model_check", "hibernation", "location"),
    root = root,
    name = "hibernation_model_check_location",
    title = "Predictions at the location level",
    description = paste(
      "Predictions of the model at the location level.",
      "The predictions are at the logit scale for binomial models and the log",
      "scale for nbinomial models."
    ),
    field_description = c(
      species = "The code of the species group",
      model_type = paste(
        "A short description of the model used to analyse the data"
      ),
      analysis = "The file fingerprint of the analysis",
      fingerprint = "The status fingerprint of the analysis",
      location_id = "The unique identifier of the location",
      winter = "The winter season defined by the year in which January falls",
      relative_mean = paste(
        "The mean of the difference between the local trend and the global",
        "trend"
      ),
      absolute_mean = paste(
        "The mean of the local trend"
      ),
      relative_sd = paste(
        "The standard error of the difference between the local trend and the",
        "global trend"
      ),
      absolute_sd = paste(
        "The standard error of the local trend"
      )
    )
  )

  # sublocation effect
  x@Data |>
    filter(!is.na(.data$sublocation_id)) |>
    distinct(
      .data$location_id,
      sublocation_i = as.integer(.data$sublocation_i),
      .data$sublocation_id,
      .data$winter,
      .data$winter_l,
      .data$winter_q,
      .data$winter_c
    ) -> sublocation_winter
  post_sample |>
    filter(str_detect(.data$parameter, "^sublocation_i:")) |>
    transmute(
      sublocation_i = str_remove(.data$parameter, "sublocation_i:"),
      .data$sim,
      q0 = .data$estimate
    ) |>
    left_join(
      post_sample |>
        filter(str_detect(.data$parameter, "^sublocation_l:")) |>
        transmute(
          sublocation_i = str_remove(.data$parameter, "sublocation_l:"),
          .data$sim,
          q1 = .data$estimate
        ),
      by = c("sublocation_i", "sim")
    ) |>
    left_join(
      post_sample |>
        filter(str_detect(.data$parameter, "^sublocation_q:")) |>
        transmute(
          sublocation_i = str_remove(.data$parameter, "sublocation_q:"),
          .data$sim,
          q2 = .data$estimate
        ),
      by = c("sublocation_i", "sim")
    ) |>
    left_join(
      post_sample |>
        filter(str_detect(.data$parameter, "^sublocation_c:")) |>
        transmute(
          sublocation_i = str_remove(.data$parameter, "sublocation_c:"),
          .data$sim,
          q3 = .data$estimate
        ),
      by = c("sublocation_i", "sim")
    ) |>
    mutate(
      sublocation_i = as.integer(.data$sublocation_i)
    ) |>
    inner_join(
      sublocation_winter |>
        group_by(
          .data$location_id,
          .data$sublocation_id,
          .data$sublocation_i
        ) |>
        summarise(
          across(
            c("winter_l", "winter_q", "winter_c"),
            ~ !any(is.na(.x))
          ),
          .groups = "drop"
        ),
      by = "sublocation_i"
    ) |>
    transmute(
      .data$location_id,
      .data$sublocation_id,
      .data$sim,
      .data$q0,
      q1 = ifelse(.data$winter_l, .data$q1, 0),
      q2 = ifelse(.data$winter_q, .data$q2, 0),
      q3 = ifelse(.data$winter_c, .data$q3, 0),
      across(c("q0", "q1", "q2", "q3"), ~ replace_na(.x, 0))
    ) -> ps_sublocation
  x@AnalysisMetadata |>
    select(
      species = "species_group_id",
      "model_type",
      analysis = "file_fingerprint",
      fingerprint = "status_fingerprint"
    ) |>
    bind_cols(
      ps_sublocation |>
        group_by(.data$location_id, .data$sublocation_id) |>
        summarise(
          across(
            c("q0", "q1", "q2", "q3"),
            list(
              mean = ~ mean(.x, na.rm = TRUE),
              sd = ~ sd(.x, na.rm = TRUE)
            )
          ),
          .groups = "drop"
        )
    ) |>
    write_vc(
      file = file.path("model_check", "hibernation", "sublocation_effect"),
      root = root,
      append = TRUE,
      sorting = c(
        "model_type",
        "species",
        "location_id",
        "sublocation_id",
        "analysis"
      ),
      digits = 4
    )
  update_metadata(
    file = file.path("model_check", "hibernation", "sublocation_effect"),
    root = root,
    name = "hibernation_model_check_sublocation",
    title = "Sublocation effect",
    description = paste(
      "Posterior distribution of the sublocation effect at the sublocation",
      "level."
    ),
    field_description = c(
      species = "The code of the species group",
      model_type = paste(
        "A short description of the model used to analyse the data"
      ),
      analysis = "The file fingerprint of the analysis",
      fingerprint = "The status fingerprint of the analysis",
      location_id = "The unique identifier of the location",
      sublocation_id = paste(
        "The unique identifier of the sublocation within a location"
      ),
      q0_mean = paste(
        "The mean of the intercept of the sublocation effect on the logit",
        "scale for binomial models and the log scale for nbinomial models"
      ),
      q1_mean = paste(
        "The mean of the linear effect of the sublocation effect on the logit",
        "scale for binomial models and the log scale for nbinomial models"
      ),
      q2_mean = paste(
        "The mean of the quadratic effect of the sublocation effect on the",
        "logit scale for binomial models and the log scale for nbinomial models"
      ),
      q3_mean = paste(
        "The mean of the cubic effect of the sublocation effect on the logit",
        "scale for binomial models and the log scale for nbinomial models"
      ),
      q0_sd = paste(
        "The standard error of the intercept of the sublocation effect on the",
        "logit scale for binomial models and the log scale for nbinomial models"
      ),
      q1_sd = paste(
        "The standard error of the linear effect of the sublocation effect on",
        "the logit scale for binomial models and the log scale for nbinomial",
        "models"
      ),
      q2_sd = paste(
        "The standard error of the quadratic effect of the sublocation effect",
        "on the logit scale for binomial models and the log scale for",
        "nbinomial models"
      ),
      q3_sd = paste(
        "The standard error of the cubic effect of the sublocation effect on",
        "the logit scale for binomial models and the log scale for nbinomial",
        "models"
      )
    )
  )

  # predictions at the sublocation level
  x@AnalysisMetadata |>
    select(
      species = "species_group_id",
      "model_type",
      analysis = "file_fingerprint",
      fingerprint = "status_fingerprint"
    ) |>
    bind_cols(
      ps_sublocation |>
        inner_join(
          sublocation_winter |>
            mutate(
              across(c("winter_l", "winter_q", "winter_c"), ~ replace_na(.x, 0))
            ),
          by = c("location_id", "sublocation_id"),
          relationship = "many-to-many"
        ) |>
        inner_join(
          ps_location |>
            select("location_id", "winter", "sim", "absolute"),
          by = c("location_id", "winter", "sim")
        ) |>
        transmute(
          .data$location_id,
          .data$sublocation_id,
          .data$sim,
          .data$winter,
          relative = .data$q0 +
            .data$q1 * .data$winter_l +
            .data$q2 * .data$winter_q +
            .data$q3 * .data$winter_c,
          absolute = .data$absolute + .data$relative
        ) |>
        group_by(.data$location_id, .data$sublocation_id, .data$winter) |>
        summarise(
          across(
            c("relative", "absolute"),
            list(
              mean = ~ mean(.x, na.rm = TRUE),
              sd = ~ sd(.x, na.rm = TRUE)
            )
          ),
          .groups = "drop"
        )
    ) |>
    write_vc(
      file = file.path("model_check", "hibernation", "sublocation"),
      root = root,
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
    file = file.path("model_check", "hibernation", "sublocation"),
    root = root,
    name = "hibernation_model_check_sublocation_prediction",
    title = "Predictions at the sublocation level",
    description = paste(
      "Predictions of the model at the sublocation level.",
      "The predictions are at the logit scale for binomial models and the log",
      "scale for nbinomial models."
    ),
    field_description = c(
      species = "The code of the species group",
      model_type = paste(
        "A short description of the model used to analyse the data"
      ),
      analysis = "The file fingerprint of the analysis",
      fingerprint = "The status fingerprint of the analysis",
      location_id = "The unique identifier of the location",
      sublocation_id = "The unique identifier of the sublocation within a location",
      winter = "The winter season defined by the year in which January falls",
      relative_mean = paste(
        "The mean of the difference between the local trend and the global",
        "trend"
      ),
      absolute_mean = paste(
        "The mean of the local trend"
      ),
      relative_sd = paste(
        "The standard error of the difference between the local trend and the",
        "global trend"
      ),
      absolute_sd = paste(
        "The standard error of the local trend"
      )
    )
  )

  # hyperparameters
  x@AnalysisMetadata |>
    select(
      species = "species_group_id",
      "model_type",
      analysis = "file_fingerprint",
      fingerprint = "status_fingerprint"
    ) |>
    bind_cols(
      names(x@Model$marginals.hyperpar) |>
        lapply(
          x = x,
          function(y, x) {
            if (!grepl("^Precision for", y)) {
              inla.zmarginal(x@Model$marginals.hyperpar[[y]], silent = TRUE) |>
                as.data.frame() |>
                transmute(
                  parameter = y,
                  .data$mean,
                  lcl = .data$quant0.025,
                  ucl = .data$quant0.975
                ) -> z
              return(z)
            }
            x@Model$marginals.hyperpar[[y]] |>
              prec2sd() |>
              transmute(
                parameter = gsub("Precision", "Stdev", y),
                .data$mean,
                lcl = .data$quant0.025,
                ucl = .data$quant0.975
              )
          }
        ) |>
        bind_rows()
    ) |>
    write_vc(
      file = file.path("model_check", "hibernation", "hyperparameters"),
      root = root,
      append = TRUE,
      sorting = c("model_type", "species", "parameter", "analysis"),
      digits = 4
    )
  update_metadata(
    file = file.path("model_check", "hibernation", "hyperparameters"),
    root = root,
    name = "hibernation_model_check_hyperparameters",
    title = "Model hyperparameters",
    description = paste(
      "Posterior distribution of the hyperparameters of the model.",
      "The model is a binomial spatio-temporal model with a logit link."
    ),
    field_description = c(
      species = "The code of the species group",
      model_type = paste(
        "A short description of the model used to analyse the data"
      ),
      analysis = "The file fingerprint of the analysis",
      fingerprint = "The status fingerprint of the analysis",
      parameter = "The name of the hyperparameter",
      mean = "The mean of the posterior distribution of the hyperparameter",
      lcl = paste(
        "The lower 95% credible limit of the posterior distribution of the",
        "hyperparameter"
      ),
      ucl = paste(
        "The upper 95% credible limit of the posterior distribution of the",
        "hyperparameter"
      )
    )
  )

  # prediction for the mesh
  mesh <- spde2mesh(x@Spde)
  as.data.frame(mesh$loc) |>
    st_as_sf(coords = c("V1", "V2")) |>
    st_union() |>
    st_convex_hull() -> hull
  hull |>
    st_sample(type = "hexagonal", size = floor(st_area(hull) / 25)) |>
    st_coordinates() |>
    as.data.frame() |>
    mutate(field_id = row_number()) -> sample_field
  projector <- inla.mesh.projector(
    mesh = mesh,
    loc = as.matrix(sample_field[, c("X", "Y")])
  )
  x@AnalysisMetadata |>
    select(
      species = "species_group_id",
      "model_type",
      analysis = "file_fingerprint",
      fingerprint = "status_fingerprint"
    ) |>
    bind_cols(
      as.matrix(projector$proj$A %*% ps_matern) |>
        as.data.frame() |>
        mutate(field_id = sample_field$field_id) |>
        pivot_longer(
          starts_with("sim"),
          names_to = "sim",
          values_to = "mesh"
        ) |>
        group_by(.data$field_id) |>
        summarise(
          across(
            "mesh",
            list(
              mean = ~ mean(.x, na.rm = TRUE) |>
                round(4),
              sd = ~ sd(.x, na.rm = TRUE) |>
                round(4)
            )
          ),
          .groups = "drop"
        ) |>
        inner_join(sample_field, by = "field_id")
    ) |>
    write_vc(
      file = file.path("model_check", "hibernation", "mesh_prediction"),
      root = root,
      append = TRUE,
      sorting = c("model_type", "species", "X", "Y", "analysis"),
      digits = c(X = 6, Y = 6, mesh_mean = 4, mesh_sd = 4)
    )
  update_metadata(
    file = file.path("model_check", "hibernation", "mesh_prediction"),
    root = root,
    name = "hibernation_model_check_mesh_prediction",
    title = "Prediction for the mesh",
    description = paste(
      "Prediction of the spatial random effect at the mesh nodes.",
      "The predictions are at the logit scale for binomial models and the log",
      "scale for nbinomial models."
    ),
    field_description = c(
      species = "The code of the species group",
      model_type = paste(
        "A short description of the model used to analyse the data"
      ),
      analysis = "The file fingerprint of the analysis",
      fingerprint = "The status fingerprint of the analysis",
      field_id = "The unique identifier of the mesh node",
      X = paste(
        "The X coordinate of the prediction on the mesh in kilometres using",
        "EPSG:31370"
      ),
      Y = paste(
        "The Y coordinate of the prediction on the mesh in kilometres using",
        "EPSG:31370"
      ),
      mesh_mean = paste(
        "The mean of the spatial random effect at the mesh node on the logit",
        "scale for binomial models and the log scale for nbinomial models"
      ),
      mesh_sd = paste(
        "The standard error of the spatial random effect at the mesh node on",
        "the logit scale for binomial models and the log scale for nbinomial",
        "models"
      )
    )
  )

  rm(
    hull,
    location_type,
    location_winter,
    loc_coordinates,
    mesh,
    post_sample,
    projector,
    ps_cwinter,
    ps_intercept,
    ps_location,
    ps_matern,
    ps_sublocation,
    sample_field,
    sublocation_winter,
    x
  )
  gc(verbose = FALSE)
  return(invisible(NULL))
}

#' @export
#' @importFrom assertthat assert_that
#' @importFrom dplyr bind_cols case_when group_by inner_join mutate select
#' starts_with summarise
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
          .data$location_id,
          type = case_when(
            .data$fortress > 0 ~ "fortress",
            .data$marl_quarry > 0 ~ "marl quarry",
            .data$other_large > 0 ~ "other large",
            .data$small > 0 ~ "small"
          ) |>
            factor(c("fortress", "marl quarry", "other large", "small")),
          .data$sublocation_id,
          .data$winter
        ) |>
        bind_cols(x@Hurdle@Imputation) |>
        pivot_longer(
          starts_with("Imputation"),
          names_to = "imputation",
          values_to = "number"
        ) |>
        group_by(
          .data$location_id,
          .data$type,
          .data$sublocation_id,
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
