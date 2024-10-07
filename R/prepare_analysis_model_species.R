#' Prepare the analysis models for a single species
#' @inheritParams prepare_analysis_data
#' @inheritParams n2kanalysis::store_model
#' @param max_dist The maximum distance for the range in kilometers
#' @export
#' @importFrom dplyr bind_rows distinct filter inner_join left_join mutate
#' select slice_max
#' @importFrom git2rdata recent_commit verify_vc
#' @importFrom n2kanalysis n2k_aggregate n2k_hurdle_imputed n2k_model_imputed
#' n2k_spde spde store_model
#' @importFrom rlang .data
#' @importFrom splines ns
prepare_analysis_model_species <- function(
  analysis_data, base, species, max_dist = 10, project = "batanalysis",
  overwrite = FALSE
) {
  file.path("hibernation", tolower(species), "locations") |>
    recent_commit(root = analysis_data, data = TRUE) |>
    bind_rows(
      file.path("hibernation", tolower(species), "analysis_data") |>
        recent_commit(root = analysis_data, data = TRUE)
    ) |>
    slice_max(.data$when, n = 1, with_ties = FALSE) |>
    distinct() -> rc
  file.path("hibernation", tolower(species), "locations") |>
    verify_vc(
      root = analysis_data,
      variables = c("location_id", x = "X", y = "Y")
    ) -> location
  location |>
    select(-"location_id") |>
    as.data.frame() |>
    spde(range = c(max_dist, 0.9), sigma = c(1, 0.01)) -> spde
  file.path("hibernation", tolower(species), "analysis_data") |>
    verify_vc(
      root = analysis_data,
      variables = c("location_id", "sublocation_id", "winter", "number")
    ) |>
    inner_join(location, by = "location_id") |>
    mutate(
      cwinter = .data$winter + 1 - min(.data$winter),
      observation_id = .data$sample_id
    ) -> dataset
  dataset |>
    distinct(.data$cwinter) -> splines
  n_df <- 3
  ns(splines$cwinter, df = n_df) |>
    as.data.frame() |>
    `colnames<-`(sprintf("knot_%02i", seq_len(n_df))) |>
    bind_cols(splines) -> splines

  # model for presence of the species
  dataset |>
    transmute(
      .data$location_id, .data$sublocation_id,
      location1 = .data$location_id, sublocation1 = .data$sublocation_id,
      location2 = .data$location_id, sublocation2 = .data$sublocation_id,
      location3 = .data$location_id, sublocation3 = .data$sublocation_id,
      .data$observation_id, .data$winter, .data$cwinter,
      present = as.integer(.data$number > 0), intercept = 1,
      datafield_id = "analysis_data", .data$X, .data$Y
    ) |>
    inner_join(splines, by = "cwinter") |>
    as.data.frame() -> data
  presence <- n2k_spde(
    formula = "
present ~ 0 + intercept +
    f(
      cwinter, model = \"rw1\",
      hyper = list(theta = list(prior = \"pc.prec\", param = c(0.15, 0.05)))
    ) +
    f(
      location_id, model = \"iid\",
      hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.05)))
    ) +
    f(
      location1, knot_01, model = \"iid\",
      hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.01)))
    ) +
    f(
      location2, knot_02, model = \"iid\",
      hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.01)))
    ) +
    f(
      location3, knot_03, model = \"iid\",
      hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.01)))
    ) +
    f(
      sublocation_id, model = \"iid\",
      hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.05)))
    ) +
    f(
      sublocation1, knot_01, model = \"iid\",
      hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.01)))
    ) +
    f(
      sublocation2, knot_02, model = \"iid\",
      hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.01)))
    ) +
    f(
      sublocation3, knot_03, model = \"iid\",
      hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.01)))
    )",
    model_type = "inla binomial: SPDE + Winter * (1 + Location + SubLocation)",
    data = data, result_datasource_id = "git", scheme_id = "hibernating bats",
    family = "binomial", species_group_id = species, spde = spde,
    location_group_id = "Flanders", seed = 19911204,
    spde_prior = list(range = c(max_dist, 0.9), sigma = c(1, 0.01)),
    first_imported_year = min(dataset$winter), analysis_date = rc$when,
    last_imported_year = max(dataset$winter)
  )
  store_model(presence, base = base, project = project, overwrite = overwrite)

  # model for the number of individuals conditional on the presence of the
  # species
  dataset |>
    transmute(
      .data$location_id, location1 = .data$location_id,
      location2 = .data$location_id, location3 = .data$location_id,
      .data$sublocation_id, sublocation1 = .data$sublocation_id,
      sublocation2 = .data$sublocation_id, sublocation3 = .data$sublocation_id,
      .data$observation_id, .data$winter, .data$cwinter,
      number = ifelse(.data$number > 0, .data$number, NA), intercept = 1,
      datafield_id = "analysis_data", .data$X, .data$Y
    ) |>
    inner_join(splines, by = "cwinter") |>
    as.data.frame() -> data
  file.path("hibernation", tolower(species), "rare_sublocation") |>
    verify_vc(
      root = analysis_data,
      variables = c("location_id", "sublocation_id", "winter", "number")
    ) |>
    left_join(location, by = "location_id") |>
    transmute(
      .data$location_id, location1 = .data$location_id,
      location2 = .data$location_id, location3 = .data$location_id,
      .data$sublocation_id, sublocation1 = .data$sublocation_id,
      sublocation2 = .data$sublocation_id, sublocation3 = .data$sublocation_id,
      .data$number, intercept = 1,
      .data$winter, cwinter = .data$winter + 1 - min(dataset$winter),
      observation_id = .data$sample_id, datafield_id = "analysis_data", .data$X,
      .data$Y
    ) |>
    inner_join(splines, by = "cwinter") |>
    as.data.frame() -> extra
  count <- n2k_spde(
    formula = "
number ~ 0 + intercept +
    f(
      cwinter, model = \"rw1\",
      hyper = list(theta = list(prior = \"pc.prec\", param = c(0.2, 0.05)))
    ) +
    f(
      location_id, model = \"iid\",
      hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.05)))
    ) +
    f(
      location1, knot_01, model = \"iid\",
      hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.01)))
    ) +
    f(
      location2, knot_02, model = \"iid\",
      hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.01)))
    ) +
    f(
      location3, knot_03, model = \"iid\",
      hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.01)))
    ) +
    f(
      sublocation_id, model = \"iid\",
      hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.05)))
    ) +
    f(
      sublocation1, knot_01, model = \"iid\",
      hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.01)))
    ) +
    f(
      sublocation2, knot_02, model = \"iid\",
      hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.01)))
    ) +
    f(
      sublocation3, knot_03, model = \"iid\",
      hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.01)))
    )",
    model_type =
    "inla zeroinflatednbinomial0: SPDE + Winter * (1 + Location + SubLocation)",
    data = data, result_datasource_id = "git", scheme_id = "hibernating bats",
    family = "zeroinflatednbinomial0", species_group_id = species, spde = spde,
    location_group_id = "Flanders", seed = 19911204, extra = extra,
    spde_prior = list(range = c(max_dist, 0.01), sigma = c(1, 0.01)),
    first_imported_year = min(dataset$winter), analysis_date = rc$when,
    last_imported_year = max(dataset$winter), imputation_size = 100,
    control = list(
      control.family = list(
        list(hyper = list(theta2 = list(initial = -20, fixed = TRUE)))
      )
    )
  )
  store_model(count, base = base, project = project, overwrite = overwrite)

  # combine the models into a hurdle model
  hurdle <- n2k_hurdle_imputed(presence = presence, count = count)
  store_model(hurdle, base = base, project = project, overwrite = overwrite)

  # create the aggregation by winter
  aggregated_tot <- n2k_aggregate(
    result_datasource_id = hurdle@AnalysisMetadata$result_datasource_id,
    scheme_id = hurdle@AnalysisMetadata$scheme_id,
    species_group_id = hurdle@AnalysisMetadata$species_group_id,
    location_group_id = hurdle@AnalysisMetadata$location_group_id,
    model_type = "aggregate imputed: sum ~ winter",
    formula = "~winter", fun = sum, status = "waiting",
    parent = hurdle@AnalysisMetadata$file_fingerprint,
    first_imported_year = hurdle@AnalysisMetadata$first_imported_year,
    last_imported_year = hurdle@AnalysisMetadata$last_imported_year,
    duration = hurdle@AnalysisMetadata$duration,
    last_analysed_year = hurdle@AnalysisMetadata$last_analysed_year,
    analysis_date = hurdle@AnalysisMetadata$analysis_date
  )
  store_model(
    aggregated_tot, base = base, project = project, overwrite = overwrite
  )

  extractor_fun <- function(model) {
    rbind(
      model$summary.lincomb.derived[, c("mean", "sd")],
      model$summary.random$cwinter[, c("mean", "sd")]
    )
  }

  prepare_model_args_fun <- function(model) {
    if (nrow(model@AggregatedImputed@Covariate) == 0) {
      return(NULL)
    }
    stopifnot(requireNamespace("INLA", quietly = TRUE))
    stopifnot(requireNamespace("n2kanalysis", quietly = TRUE))
    if (max(apply(model@AggregatedImputed@Imputation, 1, min)) < 5) {
      return(NULL)
    }
    apply(model@AggregatedImputed@Imputation, 1, max) |>
      aggregate(
        by = model@AggregatedImputed@Covariate["winter"], FUN = max
      ) -> mi
    mi <- mi[order(mi$winter), ]
    winters <- mi$winter[cumsum(mi$x) > 0 & rev(cumsum(rev(mi$x))) > 0]
    if (length(winters) < 5) {
      return(NULL)
    }
    length(winters) |>
      diag() |>
      list() |>
      setNames("cwinter") |>
      c("(Intercept)" = list(rep(1, length(winters)))) |>
      INLA::inla.make.lincombs() |>
      setNames(paste("total:", winters)) -> lc1
    comb <- expand.grid(
      winter1 = factor(winters), winter2 = factor(winters)
    )
    comb <- comb[as.integer(comb$winter1) < as.integer(comb$winter2), ]
    comb$label <- sprintf("index: %s-%s", comb$winter2, comb$winter1)
    nrow(comb) |>
      seq_len() |>
      rep(2) |>
      Matrix::sparseMatrix(
        j = c(as.integer(comb$winter2), as.integer(comb$winter1)),
        x = rep(c(1, -1), each = nrow(comb))
      ) |>
      list() |>
      setNames("cwinter") |>
      INLA::inla.make.lincombs() |>
      setNames(comb$label) -> lc2
    n2kanalysis::moving_trend(
      n_year = length(winters), duration = 12, first_year = min(winters)
    ) |>
      rbind(
        n2kanalysis::moving_trend(
          n_year = length(winters), duration = 10, first_year = min(winters)
        ),
        n2kanalysis::moving_trend(
          n_year = length(winters), duration = 6, first_year = min(winters)
        ),
        n2kanalysis::moving_difference(
          n_year = length(winters), duration = 6, first_year = min(winters)
        )
      ) |>
      unique() -> lc3
    INLA::inla.make.lincombs(cwinter = lc3) |>
      setNames(rownames(lc3)) -> lc3
    n2kanalysis::moving_average(
      n_year = length(winters), duration = 6, first_year = min(winters)
    ) |>
      rbind(
        n2kanalysis::moving_average(
          n_year = length(winters), duration = 12, first_year = min(winters)
        )
      ) -> ma
    list("(Intercept)" = rep(1, nrow(ma)), cwinter = ma) |>
      INLA::inla.make.lincombs() |>
      setNames(rownames(ma)) -> lc4
    return(list(lincomb = c(lc1, lc2, lc3, lc4)))
  }

  total_index <- n2k_model_imputed(
    result_datasource_id = aggregated_tot@AnalysisMetadata$result_datasource_id,
    scheme_id = aggregated_tot@AnalysisMetadata$scheme_id,
    species_group_id = aggregated_tot@AnalysisMetadata$species_group_id,
    location_group_id = aggregated_tot@AnalysisMetadata$location_group_id,
    model_type = "yearly imputed index: Total ~ Winter",
    formula = "~ 1 + f(cwinter, model = \"rw1\", scale.model = TRUE,
  hyper = list(theta = list(prior = \"pc.prec\", param = c(2, 0.01)))
)",
    model_fun = "INLA::inla",
    first_imported_year = aggregated_tot@AnalysisMetadata$first_imported_year,
    last_imported_year = aggregated_tot@AnalysisMetadata$last_imported_year,
    duration = aggregated_tot@AnalysisMetadata$duration,
    last_analysed_year = aggregated_tot@AnalysisMetadata$last_analysed_year,
    analysis_date = aggregated_tot@AnalysisMetadata$analysis_date,
    seed = aggregated_tot@AnalysisMetadata$seed,
    package = c("INLA", "dplyr"),
    extractor = extractor_fun,
    mutate = list(cwinter = "winter + 1 - min(winter)"),
    model_args = list(family = "nbinomial", safe = FALSE, silent = TRUE),
    prepare_model_args = list(prepare_model_args_fun),
    parent = aggregated_tot@AnalysisMetadata$file_fingerprint
  )
  store_model(
    total_index, base = base, project = project, overwrite = overwrite
  )

  # create the aggregate by winter and location
  aggregated_location <- n2k_aggregate(
    result_datasource_id = hurdle@AnalysisMetadata$result_datasource_id,
    scheme_id = hurdle@AnalysisMetadata$scheme_id,
    species_group_id = hurdle@AnalysisMetadata$species_group_id,
    location_group_id = hurdle@AnalysisMetadata$location_group_id,
    model_type = "aggregate imputed: sum ~ winter + location",
    formula = "~winter + location_id", fun = sum, status = "waiting",
    parent = hurdle@AnalysisMetadata$file_fingerprint,
    first_imported_year = hurdle@AnalysisMetadata$first_imported_year,
    last_imported_year = hurdle@AnalysisMetadata$last_imported_year,
    duration = hurdle@AnalysisMetadata$duration,
    last_analysed_year = hurdle@AnalysisMetadata$last_analysed_year,
    analysis_date = hurdle@AnalysisMetadata$analysis_date
  )
  store_model(
    aggregated_location, base = base, project = project, overwrite = overwrite
  )

  bind_rows(
    data.frame(
      analysis = c(
        get_file_fingerprint(presence), get_file_fingerprint(count)
      )
    ),
    hurdle@AnalysisRelation,
    aggregated_tot@AnalysisRelation, total_index@AnalysisRelation,
    aggregated_location@AnalysisRelation
  ) |>
    select(fingerprint = "analysis", parent = "parent_analysis")
}
