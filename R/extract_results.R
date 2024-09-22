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
#' @importFrom assertthat assert_that is.string noNA
#' @importFrom dplyr anti_join distinct
#' @importFrom git2rdata verify_vc write_vc
#' @importFrom n2kanalysis read_manifest read_model
#' @importFrom purrr walk
#'
extract_results.character <- function(
  x, base, project = "batanalysis", raw_data, root, ...
) {
  assert_that(is.string(x), noNA(x))
  verify_vc(
    "hibernation/species", root = raw_data,
    variables = c("id", "code", "name", "scientific_name", "parent")
  ) |>
    write_vc(
      "hibernation/species", root = root, sorting = "id", optimize = FALSE
    )
  verify_vc(
    "hibernation/locations", root = raw_data,
    variables = c("id", "code", "name", "parent_id")
  ) |>
    write_vc(
      "hibernation/locations", root = root, sorting = "id", optimize = FALSE
    )
  read_manifest(base = base, project = project, hash = x) |>
    slot("Manifest") |>
    distinct(.data$fingerprint) -> manifest
  if (is_git2rdata("model_check/sublocation", root = root)) {
    manifest |>
      anti_join(
        verify_vc("model_check/sublocation", root = root, "analysis"),
        by = c("fingerprint" = "analysis")
      ) -> manifest
  }
  if (is_git2rdata("hibernation/hurdle", root = root)) {
    manifest |>
      anti_join(
        verify_vc("hibernation/hurdle", root = root, "analysis"),
        by = c("fingerprint" = "analysis")
      ) -> manifest
  }
  if (is_git2rdata("hibernation/total", root = root)) {
    manifest |>
      anti_join(
        verify_vc("hibernation/total", root = root, "analysis"),
        by = c("fingerprint" = "analysis")
      ) -> manifest
  }
  if (is_git2rdata("hibernation/difference", root = root)) {
    manifest |>
      anti_join(
        verify_vc("hibernation/difference", root = root, "analysis"),
        by = c("fingerprint" = "analysis")
      ) -> manifest
  }
  walk(
    manifest$fingerprint,
    ~read_model(.x, base = base, project = project) |>
         extract_results(root = root)
  )
  return(invisible(NULL))
}

#' @export
#' @importFrom assertthat assert_that
#' @importFrom dplyr across bind_cols filter mutate select transmute
#' @importFrom git2rdata is_git2rdata update_metadata write_vc
#' @importFrom n2kanalysis get_file_fingerprint
#' @importFrom stringr str_detect str_remove
#' @importFrom tidyr separate_wider_regex
extract_results.n2kModelImputed <- function(x, root, ...) {
  assert_that(inherits(x, "n2kModelImputed"))
  message(get_file_fingerprint(x))
  # skip if model is already in the database
  if (is_git2rdata("hibernation/difference", root = root)) {
    if (
      get_file_fingerprint(x) %in%
      read_vc("hibernation/difference", root = root)$analysis
    ) {
      return(invisible(NULL))
    }
  }
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
  rm(x, results)
  gc(verbose = FALSE)
  return(invisible(NULL))
}

#' @export
#' @importFrom assertthat assert_that
#' @importFrom dplyr bind_cols group_by inner_join mutate select starts_with
#' summarise
#' @importFrom git2rdata is_git2rdata update_metadata write_vc
#' @importFrom n2kanalysis get_file_fingerprint
#' @importFrom rlang .data
#' @importFrom stats quantile
#' @importFrom tidyr pivot_longer
extract_results.n2kAggregate <- function(x, root, ...) {
  assert_that(inherits(x, "n2kAggregate"))
  message(get_file_fingerprint(x))
  # skip if model is already in the database
  if (is_git2rdata("hibernation/total", root = root)) {
    if (
      get_file_fingerprint(x) %in%
      read_vc("hibernation/total", root = root)$analysis
    ) {
      return(invisible(NULL))
    }
  }
  x@AggregatedImputed@Covariate |>
    bind_cols(x@AggregatedImputed@Imputation) |>
    pivot_longer(
      starts_with("Imputation"), names_to = "imputation", values_to = "total"
    ) -> results
  if (has_name(results, "location")) {
    x@AnalysisMetadata |>
      select(
        species = "species_group_id", "model_type",
        analysis = "file_fingerprint", fingerprint = "status_fingerprint"
      ) |>
      bind_cols(
        results |>
          group_by(
            winter = .data$winter, location = .data$location
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
  rm(x, results)
  gc(verbose = FALSE)
  return(invisible(NULL))
}

#' @export
#' @importFrom assertthat assert_that
#' @importFrom dplyr arrange bind_cols distinct group_by inner_join mutate
#' select starts_with summarise
#' @importFrom git2rdata is_git2rdata update_metadata write_vc
#' @importFrom INLA inla.mesh.projector inla.posterior.sample inla.zmarginal
#' @importFrom n2kanalysis get_file_fingerprint spde2mesh
#' @importFrom purrr map map2_dfc
#' @importFrom stringr str_detect
#' @importFrom stats quantile
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr everything pivot_longer pivot_wider
extract_results.n2kSpde <- function(x, root, ..., n_sim = 100) {
  assert_that(inherits(x, "n2kSpde"), is.count(n_sim), noNA(n_sim))
  message(get_file_fingerprint(x))
  # skip if model is not converged
  if (n2kanalysis::status(x) != "converged") {
    return(invisible(NULL))
  }
  # skip if model is already in the database
  if (is_git2rdata("model_check/sublocation", root = root)) {
    if (
      get_file_fingerprint(x) %in%
        read_vc("model_check/sublocation", root = root)$analysis
    ) {
      return(invisible(NULL))
    }
  }
  # generate posterior sample
  inla.posterior.sample(n = n_sim, result = x@Model) |>
    map("latent") |>
    map2_dfc(
      seq_len(n_sim),
      ~sprintf("sim_%04i", .y) |>
        `colnames<-`(.x, value = _) |>
        as.data.frame()
    ) |>
    rownames_to_column(var = "parameter") |>
    pivot_longer(
      -"parameter", names_to = "sim", values_to = "estimate"
    ) -> post_sample
  # extract the intercept and the cwinter effect
  post_sample |>
    filter(.data$parameter == "intercept:1") |>
    select("sim", intercept = "estimate") -> ps_intercept
  post_sample |>
    filter(str_detect(.data$parameter, "cwinter")) |>
    mutate(
      cwinter = str_remove(.data$parameter, "cwinter:") |>
        as.integer(),
      .data$sim, .data$estimate
    ) |>
    inner_join(ps_intercept, by = "sim") |>
    transmute(
      .data$cwinter, .data$sim, estimate = .data$intercept + .data$estimate
    ) -> ps_cwinter
  x@AnalysisMetadata |>
    select(
      species = "species_group_id", "model_type",
      analysis = "file_fingerprint", fingerprint = "status_fingerprint"
    ) |>
    bind_cols(
      ps_cwinter |>
        group_by(.data$cwinter) |>
        summarise(
          across(
            "estimate", .names = "{.fn}",
            list(
              mean = ~mean(.x, na.rm = TRUE) |>
                round(4),
              sd = ~sd(.x, na.rm = TRUE) |>
                round(4)
            )
          )
        )
    ) |>
    write_vc(
      file = "model_check/cwinter", root = root, append = TRUE,
      sorting = c("model_type", "species", "cwinter", "analysis")
    )
  ps_cwinter |>
    distinct(.data$cwinter) |>
    mutate(
      winter_p1 = (.data$cwinter - median(.data$cwinter)) / 10,
      winter_p2 = .data$winter_p1 ^ 2
    ) -> q_winter

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
  x@Model$.args$data[c("location", "X", "Y")] |>
    as.data.frame() |>
    distinct() |>
    filter(!is.na(.data$X)) -> loc_coordinates
  projector <- inla.mesh.projector(
    mesh = spde2mesh(x@Spde), loc = as.matrix(loc_coordinates[, c("X", "Y")])
  )
  post_sample |>
    filter(str_detect(.data$parameter, "^location:")) |>
    transmute(
      location = str_remove(.data$parameter, "location:") |>
        as.integer(),
      .data$sim, q0 = .data$estimate
    ) |>
    inner_join(
      post_sample |>
        filter(str_detect(.data$parameter, "^location2:")) |>
        transmute(
          location = str_remove(.data$parameter, "location2:") |>
            as.integer(),
          .data$sim, q1 = .data$estimate
        ),
      by = c("location", "sim")
    ) |>
    inner_join(
      post_sample |>
        filter(str_detect(.data$parameter, "^location3:")) |>
        transmute(
          location = str_remove(.data$parameter, "location3:") |>
            as.integer(),
          .data$sim, q2 = .data$estimate
        ),
      by = c("location", "sim")
    ) |>
    mutate(
      location_id = levels(x@Model$.args$data$location)[.data$location] |>
        as.integer()
    ) |>
    inner_join(
      as.matrix(projector$proj$A %*% ps_matern) |>
        as.data.frame() |>
        mutate(
          location_id = loc_coordinates$location,
          location_id = levels(.data$location_id)[.data$location_id] |>
            as.integer()
        ) |>
        pivot_longer(starts_with("sim"), names_to = "sim", values_to = "mesh"),
      by = c("location_id", "sim")
    ) -> ps_location
  x@AnalysisMetadata |>
    select(
      species = "species_group_id", "model_type",
      analysis = "file_fingerprint", fingerprint = "status_fingerprint"
    ) |>
    bind_cols(
      ps_location |>
        group_by(.data$location_id) |>
        summarise(
          across(
            c("q0", "q1", "q2", "mesh"),
            list(
              mean = ~mean(.x, na.rm = TRUE) |>
                round(4),
              sd = ~sd(.x, na.rm = TRUE) |>
                round(4)
            )
          )
        )
  ) |>
  write_vc(
    file = "model_check/location_effect", root = root, append = TRUE,
    sorting = c("model_type", "species", "location_id", "analysis")
  )

  # predictions at the location level
  ps_location |>
    mutate(cwinter = list(q_winter)) |>
    unnest("cwinter") |>
    inner_join(ps_cwinter, by = c("cwinter", "sim")) |>
    transmute(
      .data$location_id, .data$cwinter, .data$sim,
      relative = .data$q0 + .data$q1 * .data$winter_p1 +
        .data$q2 * .data$winter_p2 + .data$mesh,
      absolute = .data$estimate + .data$relative
    ) -> ps_location
  x@AnalysisMetadata |>
    select(
      species = "species_group_id", "model_type",
      analysis = "file_fingerprint", fingerprint = "status_fingerprint"
    ) |>
    bind_cols(
      ps_location |>
        group_by(.data$location_id, .data$cwinter) |>
        summarise(
          across(
            c("relative", "absolute"),
            list(
              mean = ~mean(.x, na.rm = TRUE) |>
                round(4),
              sd = ~sd(.x, na.rm = TRUE) |>
                round(4)
            )
          ),
          .groups = "drop"
        )
    ) |>
    write_vc(
      file = "model_check/location", root = root, append = TRUE,
      sorting = c("model_type", "species", "location_id", "cwinter", "analysis")
    )

  # sublocation effect
  post_sample |>
    filter(str_detect(.data$parameter, "^sublocation:")) |>
    transmute(
      sublocation = str_remove(.data$parameter, "sublocation:") |>
        as.integer(),
      .data$sim, q0 = .data$estimate
    ) |>
    inner_join(
      post_sample |>
        filter(str_detect(.data$parameter, "^sublocation2:")) |>
        transmute(
          sublocation = str_remove(.data$parameter, "sublocation2:") |>
            as.integer(),
          .data$sim, q1 = .data$estimate
        ),
      by = c("sublocation", "sim")
    ) |>
    left_join(
      post_sample |>
        filter(str_detect(.data$parameter, "^sublocation3:")) |>
        transmute(
          sublocation = str_remove(.data$parameter, "sublocation3:") |>
            as.integer(),
          .data$sim, q2 = .data$estimate
        ),
      by = c("sublocation", "sim")
    ) |>
    mutate(q2 = replace_na(.data$q2, 0)) -> ps_sublocation
  x@AnalysisMetadata |>
    select(
      species = "species_group_id", "model_type",
      analysis = "file_fingerprint", fingerprint = "status_fingerprint"
    ) |>
    bind_cols(
      ps_sublocation |>
        mutate(
          sublocation_id =
            levels(x@Model$.args$data$sublocation)[.data$sublocation]
        ) |>
        group_by(sublocation_id = as.integer(.data$sublocation_id)) |>
        summarise(
          across(
            c("q0", "q1", "q2"),
            list(
              mean = ~mean(.x, na.rm = TRUE) |>
                round(4),
              sd = ~sd(.x, na.rm = TRUE) |>
                round(4)
            )
          )
        )
    ) |>
    write_vc(
      file = "model_check/sublocation_effect", root = root, append = TRUE,
      sorting = c("model_type", "species", "sublocation_id", "analysis")
    )

  # predictions at the sublocation level
  x@AnalysisMetadata |>
    select(
      species = "species_group_id", "model_type",
      analysis = "file_fingerprint", fingerprint = "status_fingerprint"
    ) |>
    bind_cols(
      ps_sublocation |>
        mutate(
          sublocation_id =
            levels(x@Model$.args$data$sublocation)[.data$sublocation] |>
            as.integer(),
          cwinter = list(q_winter)
        ) |>
        unnest("cwinter") |>
        inner_join(
          read_vc("hibernation/locations", root = root) |>
            select(sublocation_id = "id", location_id = "parent_id"),
           by = "sublocation_id"
        ) |>
        inner_join(
          ps_location |>
            select("location_id", "cwinter", "sim", "absolute"),
          by = c("location_id", "cwinter", "sim")
        ) |>
        transmute(
          .data$sublocation_id, .data$sim, .data$cwinter,
          relative = .data$q0 + .data$q1 * .data$winter_p1 +
            .data$q2 * .data$winter_p2,
          absolute = .data$absolute + .data$relative
        ) |>
        group_by(.data$sublocation_id, .data$cwinter) |>
        summarise(
          across(
            c("relative", "absolute"),
            list(
              mean = ~mean(.x, na.rm = TRUE) |>
                round(4),
              sd = ~sd(.x, na.rm = TRUE) |>
                round(4)
            )
          ),
          .groups = "drop"
        )
    ) |>
    write_vc(
      file = "model_check/sublocation", root = root, append = TRUE,
      sorting = c(
        "model_type", "species", "sublocation_id", "cwinter", "analysis"
      )
    )

  # hyperparameters
  x@AnalysisMetadata |>
    select(
      species = "species_group_id", "model_type",
      analysis = "file_fingerprint", fingerprint = "status_fingerprint"
    ) |>
    bind_cols(
      names(x@Model$marginals.hyperpar) |>
        lapply(
          x = x,
          function(y, x) {
            if (!grepl("^Precision for", y)) {
              inla.zmarginal(x@Model$marginals.hyperpar[[7]], silent = TRUE) |>
                as.data.frame() |>
                transmute(
                  parameter = y, .data$mean,
                  lcl = .data$quant0.025, ucl = .data$quant0.975
                ) -> z
              return(z)
            }
            x@Model$marginals.hyperpar[[y]] |>
              prec2sd() |>
              transmute(
                parameter = gsub("Precision", "Stdev", y), .data$mean,
                lcl = .data$quant0.025, ucl = .data$quant0.975
              )
          }
        ) |>
        bind_rows()
    ) |>
    write_vc(
      file = "model_check/hyperparameters", root = root, append = TRUE,
      sorting = c(
        "model_type", "species", "sublocation_id", "cwinter", "analysis"
      )
    )

  rm(
    loc_coordinates, post_sample, projector, ps_cwinter, ps_intercept,
    ps_location, ps_matern, ps_sublocation, q_winter, x
  )
  gc(verbose = FALSE)
  return(invisible(NULL))
}

#' @export
#' @importFrom assertthat assert_that
#' @importFrom dplyr bind_cols group_by inner_join mutate select starts_with
#' summarise
#' @importFrom git2rdata is_git2rdata write_vc
#' @importFrom n2kanalysis get_file_fingerprint
#' @importFrom stats quantile sd
extract_results.n2kHurdleImputed <- function(x, root, ...) {
  assert_that(inherits(x, "n2kHurdleImputed"))
  message(get_file_fingerprint(x))
  # skip if model is already in the database
  if (is_git2rdata("hibernation/hurdle", root = root)) {
    if (
      get_file_fingerprint(x) %in%
      read_vc("hibernation/hurdle", root = root)$analysis
    ) {
      return(invisible(NULL))
    }
  }
  x@AnalysisMetadata |>
    select(
      species = "species_group_id", "model_type", analysis = "file_fingerprint",
      fingerprint = "status_fingerprint"
    ) |>
    bind_cols(
      x@Hurdle@Covariate |>
        transmute(
          sublocation_id =
            levels(.data$sublocation)[.data$sublocation] |>
              as.integer(),
          winter = .data$winter
        ) |>
        bind_cols(x@Hurdle@Imputation) |>
        pivot_longer(
          starts_with("Imputation"), names_to = "imputation",
          values_to = "number"
        ) |>
        group_by(.data$sublocation_id, .data$winter) |>
        summarise(
          across(
            "number", .names = "{.fn}",
            list(
              mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE),
              p05 = ~quantile(.x, 0.05, na.rm = TRUE),
              p20 = ~quantile(.x, 0.20, na.rm = TRUE),
              p35 = ~quantile(.x, 0.35, na.rm = TRUE),
              p50 = ~quantile(.x, 0.50, na.rm = TRUE),
              p65 = ~quantile(.x, 0.65, na.rm = TRUE),
              p80 = ~quantile(.x, 0.80, na.rm = TRUE),
              p95 = ~quantile(.x, 0.95, na.rm = TRUE)
            )
          ),
          .groups = "drop"
        ) |>
        mutate(delta = .data$p95 - .data$p05)
    ) |>
    write_vc(
      file = "hibernation/hurdle", root = root, optimize = FALSE,
      append = TRUE,
      sorting = c(
        "model_type", "species", "sublocation_id", "winter", "analysis"
      )
    )
  rm(x)
  gc(verbose = FALSE)
  return(invisible(NULL))
}
