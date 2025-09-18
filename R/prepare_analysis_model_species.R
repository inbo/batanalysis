#' Prepare the analysis models for a single species
#' @inheritParams n2kanalysis::store_model
#' @inheritParams n2kanalysis::display
#' @param raw_data The git repository with the raw data
#' @param species the code of the species
#' @param start_winter The starting winter (year) for the analysis.
#' This is the year of first of January.
#' For example, for the winter 2000-2001, use 2001.
#' The default is 23 years before the current year.
#' This means that in 2024 the default is 2001.
#' This means that the analysis will use data from the winter
#' 2000-2001 up to the last completed winter.
#' @param n_winter Minimum number of winters in which a (sub)location is
#' monitored.
#' @param n_present Minimum number of winters in which the species is observed
#' at a (sub)location.
#' @param n_extrapolation Maximum number of winters from the nearest winter
#' where the species was observed at a (sub)location.
#' Only do imputations for a (sub)location when the difference between the
#' nearest observed winter and the missing winter is below or equal to this
#' threshold.
#' @param max_delta Only relevant in case of multiple visits within the same
#' winter to a location divided into sublocation.
#' We select for every winter per sublocation the visit that is closest to the
#' middle of the winter (15th of January).
#' All selected visits of a location for a given winter should be within
#' `max_delta` days from each other.
#' If not, we drop the visits the furthest away from the middle of the winter
#' until the condition is met.
#' This ensures that the visits to the sublocations are not too far apart in
#' time.
#' The default is 10 days.
#' This means that all visits to sublocations of a location should be
#' within 10 days from each other.
#' This allows to spread the visits over a few days.
#' @param max_dist The maximum distance for the range in kilometres
#' @export
#' @importFrom dplyr bind_rows distinct filter inner_join left_join mutate
#' select slice_max
#' @importFrom git2rdata recent_commit verify_vc
#' @importFrom n2kanalysis display n2k_aggregate n2k_hurdle_imputed
#' n2k_model_imputed n2k_spde spde store_model
#' @importFrom rlang .data
#' @importFrom sf st_as_sf st_coordinates st_transform
#' @importFrom stats poly
prepare_analysis_model_species <- function(
  raw_data,
  species,
  start_winter = as.integer(format(Sys.Date(), "%Y")) - 23,
  base,
  project = "batanalysis",
  n_winter = 4,
  n_present = 3,
  n_extrapolation = 5,
  max_delta = 10,
  max_dist = 10,
  overwrite = FALSE,
  verbose = TRUE
) {
  assert_that(is.count(start_winter), noNA(start_winter))
  sprintf("%i-10-01", start_winter - 1) |>
    as.Date() -> start

  # prepare species information
  get_child_species(target = raw_data, species = species) |>
    inner_join(
      file.path("data", "hibernation", "species") |>
        verify_vc(root = raw_data, variables = c("id", "scientific_name")) |>
        select("id", "scientific_name"),
      by = "id"
    ) -> this_species
  stopifnot("No matching species found" = nrow(this_species) > 0)
  paste("Preparing the data for", this_species$scientific_name[1]) |>
    display(verbose = verbose)

  # prepare location information
  file.path("data", "hibernation", "locations") |>
    verify_vc(
      root = raw_data,
      variables = c("id", "longitude", "latitude", "parent_id")
    ) -> all_locations
  all_locations |>
    filter(.data$parent_id < 0) |>
    select(location_id = "id", "longitude", "latitude") |>
    left_join(
      all_locations |>
        filter(.data$parent_id > 0) |>
        distinct(location_id = .data$parent_id, detailed = TRUE),
      by = "location_id"
    ) |>
    left_join(
      file.path("data", "hibernation", "location_types") |>
        verify_vc(root = raw_data, c("location_id", "type")),
      by = "location_id"
    ) |>
    mutate(
      detailed = replace_na(.data$detailed, FALSE),
      type = replace_na(.data$type, "small")
    ) |>
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
    st_transform(crs = 31370) -> locations
  locations |>
    st_coordinates() |>
    bind_cols(locations) |>
    mutate(across(c("X", "Y"), ~ . / 1000)) -> locations

  display(verbose = verbose, "  preparing timeseries without detail")
  no_detail <- select_imputation_no_detail(
    locations = locations,
    raw_data = raw_data,
    this_species = this_species,
    start = start,
    n_winter = n_winter,
    n_present = n_present,
    n_extrapolation = n_extrapolation
  )

  display(verbose = verbose, "  preparing detailed timeseries")
  detail <- select_imputation_detail(
    locations = locations,
    raw_data = raw_data,
    this_species = this_species,
    start = start,
    n_winter = n_winter,
    n_present = n_present,
    n_extrapolation = n_extrapolation,
    max_delta = max_delta
  )

  display(verbose = verbose, "  preparing analysis objects")
  # calculate polynomial coefficients for winter
  bind_rows(detail$sublocation, detail$location, no_detail) |>
    inner_join(locations, by = "location_id") -> dataset
  min(dataset$winter) |>
    seq(max(dataset$winter)) |>
    poly(degree = 3) -> poly_winter
  dataset |>
    transmute(
      .data$observation_id,
      .data$datafield_id,
      present = as.integer(.data$number > 0),
      number = ifelse(.data$number > 0, .data$number, NA_integer_),
      location = factor(.data$location_id),
      sublocation = factor(.data$sublocation_id),
      # define the type of location
      fortress = as.integer(.data$type == "fortress"),
      marl_quarry = as.integer(.data$type == "marl quarry"),
      other_large = as.integer(.data$type == "other large"),
      small = as.integer(.data$type == "small"),
      # centre winter to the starting year
      .data$winter,
      winter_r = .data$winter - start_winter + 1,
      # linear polynomial coefficient for winter
      winter_l = poly_winter[.data$winter_r, 1],
      location_l = .data$location,
      sublocation_l = .data$sublocation,
      # quadratic polynomial coefficient for winter
      winter_q = poly_winter[.data$winter_r, 2],
      location_q = .data$location,
      sublocation_q = .data$sublocation,
      # cubic polynomial coefficient for winter
      winter_c = poly_winter[.data$winter_r, 3],
      location_c = .data$location,
      sublocation_c = .data$sublocation,
      .data$X,
      .data$Y
    ) -> dataset
  dataset |>
    group_by(.data$sublocation) |>
    mutate(
      n_winter = n_distinct(.data$winter_c[!is.na(.data$present)]),
      subwinter_l = ifelse(.data$n_winter >= 6, .data$winter_l, NA),
      subwinter_q = ifelse(.data$n_winter >= 12, .data$winter_q, NA),
      subwinter_c = ifelse(.data$n_winter >= 18, .data$winter_c, NA)
    ) |>
    group_by(.data$location) |>
    mutate(
      n_winter = n_distinct(.data$winter_c[!is.na(.data$present)]),
      winter_l = ifelse(.data$n_winter >= 6, .data$winter_l, NA),
      winter_q = ifelse(.data$n_winter >= 12, .data$winter_q, NA),
      winter_c = ifelse(.data$n_winter >= 18, .data$winter_c, NA)
    ) |>
    ungroup() |>
    select(-"n_winter", -"number") -> ds_present
  dataset |>
    group_by(.data$sublocation) |>
    mutate(
      n_winter = n_distinct(.data$winter_c[!is.na(.data$number)]),
      subwinter_l = ifelse(.data$n_winter >= 6, .data$winter_l, NA),
      subwinter_q = ifelse(.data$n_winter >= 12, .data$winter_q, NA),
      subwinter_c = ifelse(.data$n_winter >= 18, .data$winter_c, NA)
    ) |>
    group_by(.data$location) |>
    mutate(
      n_winter = n_distinct(.data$winter_c[!is.na(.data$number)]),
      winter_l = ifelse(.data$n_winter >= 6, .data$winter_l, NA),
      winter_q = ifelse(.data$n_winter >= 12, .data$winter_q, NA),
      winter_c = ifelse(.data$n_winter >= 18, .data$winter_c, NA)
    ) |>
    ungroup() |>
    mutate(observation_id = as.integer(.data$observation_id)) |>
    select(-"n_winter", -"present") -> ds_number
  detail$rare_sublocations |>
    mutate(
      location = factor(NA, levels = levels(ds_number$location)),
      sublocation = factor(NA, levels = levels(ds_number$sublocation)),
      winter_r = .data$winter - start_winter + 1,
      location_l = .data$location,
      location_q = .data$location,
      location_c = .data$location,
      sublocation_l = .data$sublocation,
      sublocation_q = .data$sublocation,
      sublocation_c = .data$sublocation
    ) -> extra
  missing_cols <- colnames(ds_number)[!colnames(ds_number) %in% colnames(extra)]
  list(NA_real_) |>
    rep(length(missing_cols)) |>
    setNames(missing_cols) |>
    as.data.frame() |>
    bind_cols(extra) -> extra

  file.path("data", "hibernation", "locations") |>
    recent_commit(root = raw_data, data = TRUE) |>
    bind_rows(
      file.path("data", "hibernation", "observations") |>
        recent_commit(root = raw_data, data = TRUE),
      file.path("data", "hibernation", "samples") |>
        recent_commit(root = raw_data, data = TRUE),
      file.path("data", "hibernation", "species") |>
        recent_commit(root = raw_data, data = TRUE),
      file.path("data", "hibernation", "totals") |>
        recent_commit(root = raw_data, data = TRUE),
      file.path("data", "hibernation", "visits") |>
        recent_commit(root = raw_data, data = TRUE)
    ) |>
    slice_max(.data$when, n = 1, with_ties = FALSE) |>
    distinct() -> rc

  # prepare spde object
  locations |>
    select(c("X", "Y")) |>
    spde(range = c(max_dist, 0.9), sigma = c(1, 0.01)) -> spde
  presence <- n2k_spde(
    formula = "
present ~ 0 + fortress + marl_quarry + other_large + small +
      f(
        winter_r,
        model = \"rw1\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(0.15, 0.05)))
      ) +
      f(
        location,
        model = \"iid\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.05)))
      ) +
      f(
        location_l,
        winter_l,
        model = \"iid\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.01)))
      ) +
      f(
        location_q,
        winter_q,
        model = \"iid\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.01)))
      ) +
      f(
        location_c,
        winter_c,
        model = \"iid\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.01)))
      ) +
      f(
        sublocation,
        model = \"iid\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.05)))
      ) +
      f(
        sublocation_l,
        subwinter_l,
        model = \"iid\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(0.5, 0.01)))
      ) +
      f(
        sublocation_q,
        subwinter_q,
        model = \"iid\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(0.5, 0.01)))
      ) +
      f(
        sublocation_c,
        subwinter_c,
        model = \"iid\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(0.5, 0.01)))
      )",
    model_type = "inla binomial: SPDE + Winter * (1 + Location + SubLocation)",
    data = ds_present,
    result_datasource_id = "git",
    scheme_id = "hibernating bats",
    family = "binomial",
    species_group_id = species,
    spde = spde,
    location_group_id = "Flanders",
    seed = 19911204,
    spde_prior = list(range = c(max_dist, 0.9), sigma = c(1, 0.01)),
    first_imported_year = min(dataset$winter),
    analysis_date = rc$when,
    last_imported_year = max(dataset$winter) - 1,
    imputation_size = 100
  )
  store_model(presence, base = base, project = project, overwrite = overwrite)

  count <- n2k_spde(
    formula = "
number ~ 0 + fortress + marl_quarry + other_large + small +
      f(
        winter_r,
        model = \"rw1\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(0.15, 0.05)))
      ) +
      f(
        location,
        model = \"iid\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.05)))
      ) +
      f(
        location_l,
        winter_l,
        model = \"iid\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.01)))
      ) +
      f(
        location_q,
        winter_q,
        model = \"iid\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.01)))
      ) +
      f(
        location_c,
        winter_c,
        model = \"iid\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.01)))
      ) +
      f(
        sublocation,
        model = \"iid\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.05)))
      ) +
      f(
        sublocation_l,
        subwinter_l,
        model = \"iid\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(0.5, 0.01)))
      ) +
      f(
        sublocation_q,
        subwinter_q,
        model = \"iid\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(0.5, 0.01)))
      ) +
      f(
        sublocation_c,
        subwinter_c,
        model = \"iid\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(0.5, 0.01)))
      )",
    model_type = paste(
      "inla zeroinflatednbinomial0:",
      "SPDE + Winter * (1 + Location + SubLocation)"
    ),
    data = ds_number,
    result_datasource_id = "git",
    scheme_id = "hibernating bats",
    family = "zeroinflatednbinomial0",
    species_group_id = species,
    spde = spde,
    location_group_id = "Flanders",
    seed = 19911204,
    extra = extra,
    spde_prior = list(range = c(max_dist, 0.01), sigma = c(1, 0.01)),
    first_imported_year = min(dataset$winter),
    analysis_date = rc$when,
    last_imported_year = max(dataset$winter) - 1,
    imputation_size = 100,
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
    formula = "~winter",
    fun = sum,
    status = "waiting",
    parent = hurdle@AnalysisMetadata$file_fingerprint,
    first_imported_year = hurdle@AnalysisMetadata$first_imported_year,
    last_imported_year = hurdle@AnalysisMetadata$last_imported_year,
    duration = hurdle@AnalysisMetadata$duration,
    last_analysed_year = hurdle@AnalysisMetadata$last_analysed_year,
    analysis_date = hurdle@AnalysisMetadata$analysis_date
  )
  store_model(
    aggregated_tot,
    base = base,
    project = project,
    overwrite = overwrite
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
        by = model@AggregatedImputed@Covariate["winter"],
        FUN = max
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
      winter1 = factor(winters),
      winter2 = factor(winters)
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
      n_year = length(winters),
      duration = 12,
      first_year = min(winters)
    ) |>
      rbind(
        n2kanalysis::moving_trend(
          n_year = length(winters),
          duration = 10,
          first_year = min(winters)
        ),
        n2kanalysis::moving_trend(
          n_year = length(winters),
          duration = 6,
          first_year = min(winters)
        ),
        n2kanalysis::moving_difference(
          n_year = length(winters),
          duration = 6,
          first_year = min(winters)
        )
      ) |>
      unique() -> lc3
    INLA::inla.make.lincombs(cwinter = lc3) |>
      setNames(rownames(lc3)) -> lc3
    n2kanalysis::moving_average(
      n_year = length(winters),
      duration = 6,
      first_year = min(winters)
    ) |>
      rbind(
        n2kanalysis::moving_average(
          n_year = length(winters),
          duration = 12,
          first_year = min(winters)
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
    total_index,
    base = base,
    project = project,
    overwrite = overwrite
  )

  # aggregate by type
  aggregated_type <- n2k_aggregate(
    result_datasource_id = hurdle@AnalysisMetadata$result_datasource_id,
    scheme_id = hurdle@AnalysisMetadata$scheme_id,
    species_group_id = hurdle@AnalysisMetadata$species_group_id,
    location_group_id = hurdle@AnalysisMetadata$location_group_id,
    model_type = "aggregate imputed: sum ~ winter + type",
    formula = "~winter + fortress + marl_quarry + other_large + small",
    fun = sum,
    status = "waiting",
    parent = hurdle@AnalysisMetadata$file_fingerprint,
    first_imported_year = hurdle@AnalysisMetadata$first_imported_year,
    last_imported_year = hurdle@AnalysisMetadata$last_imported_year,
    duration = hurdle@AnalysisMetadata$duration,
    last_analysed_year = hurdle@AnalysisMetadata$last_analysed_year,
    analysis_date = hurdle@AnalysisMetadata$analysis_date
  )
  store_model(
    aggregated_type,
    base = base,
    project = project,
    overwrite = overwrite
  )

  # create the aggregate by winter and location
  aggregated_location <- n2k_aggregate(
    result_datasource_id = hurdle@AnalysisMetadata$result_datasource_id,
    scheme_id = hurdle@AnalysisMetadata$scheme_id,
    species_group_id = hurdle@AnalysisMetadata$species_group_id,
    location_group_id = hurdle@AnalysisMetadata$location_group_id,
    model_type = "aggregate imputed: sum ~ winter + location",
    formula = "~winter + location",
    fun = sum,
    status = "waiting",
    parent = hurdle@AnalysisMetadata$file_fingerprint,
    first_imported_year = hurdle@AnalysisMetadata$first_imported_year,
    last_imported_year = hurdle@AnalysisMetadata$last_imported_year,
    duration = hurdle@AnalysisMetadata$duration,
    last_analysed_year = hurdle@AnalysisMetadata$last_analysed_year,
    analysis_date = hurdle@AnalysisMetadata$analysis_date
  )
  store_model(
    aggregated_location,
    base = base,
    project = project,
    overwrite = overwrite
  )

  bind_rows(
    data.frame(
      analysis = c(
        get_file_fingerprint(presence),
        get_file_fingerprint(count)
      )
    ),
    hurdle@AnalysisRelation,
    aggregated_tot@AnalysisRelation,
    total_index@AnalysisRelation,
    aggregated_type@AnalysisRelation,
    aggregated_location@AnalysisRelation
  ) |>
    select(fingerprint = "analysis", parent = "parent_analysis")
}
