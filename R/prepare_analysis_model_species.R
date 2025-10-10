#' Prepare the analysis models for a single species
#' @inheritParams n2kanalysis::store_model
#' @inheritParams n2kanalysis::display
#' @param raw_data The git repository with the raw data
#' @param species the code of the species
#' @param n_extrapolation Impute only locations or sublocations when there is no
#' more than `n_extrapolation` winters between the observation to impute and the
#' nearest observation.
#' @param max_dist The maximum distance for the range in kilometres
#' @param visits The output of `visit_type()`
#' @export
#' @importFrom dplyr bind_rows distinct filter inner_join left_join mutate
#' n_distinct select slice_max
#' @importFrom git2rdata recent_commit verify_vc
#' @importFrom n2kanalysis display n2k_aggregate n2k_hurdle_imputed
#' n2k_model_imputed n2k_spde spde store_model
#' @importFrom rlang .data sym
#' @importFrom sf st_as_sf st_coordinates st_drop_geometry st_transform
#' @importFrom stats poly
#' @importFrom tidyr complete nesting replace_na
prepare_analysis_model_species <- function(
  raw_data,
  species,
  base,
  project = "batanalysis",
  n_extrapolation = 24,
  max_dist = 10,
  visits,
  overwrite = FALSE,
  verbose = TRUE
) {
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
  visits |>
    filter(.data$type %in% c("detail", "total")) |>
    distinct(
      .data$location_id,
      detailed = !is.na(.data$level)
    ) |>
    inner_join(
      file.path("data", "hibernation", "locations") |>
        verify_vc(
          root = raw_data,
          variables = c("id", "longitude", "latitude", "parent_id")
        ),
      by = c("location_id" = "id")
    ) |>
    left_join(
      file.path("data", "hibernation", "location_types") |>
        verify_vc(root = raw_data, c("location_id", "type")),
      by = "location_id"
    ) |>
    transmute(
      .data$location_id,
      .data$longitude,
      .data$latitude,
      .data$detailed,
      type = replace_na(.data$type, "small")
    ) |>
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
    st_transform(crs = 31370) -> locations
  locations |>
    st_coordinates() |>
    bind_cols(st_drop_geometry(locations)) |>
    mutate(across(c("X", "Y"), ~ . / 1000)) -> locations

  display(verbose = verbose, "  preparing data")
  file.path("data", "hibernation", "totals") |>
    verify_vc(
      root = raw_data,
      variables = c("visit_id", "species_id", "total")
    ) |>
    semi_join(this_species, by = c("species_id" = "id")) |>
    group_by(.data$visit_id) |>
    summarise(number = sum(.data$total), .groups = "drop") |>
    left_join(
      x = visits |>
        filter(
          .data$type == "total",
          is.na(.data$level) | .data$level != "detail"
        ),
      by = "visit_id"
    ) |>
    transmute(
      observation_id = .data$visit_id,
      datafield_id = 2L,
      .data$location_id,
      .data$winter,
      number = replace_na(.data$number, 0)
    ) |>
    group_by(.data$location_id) |>
    filter(max(.data$number, na.rm = TRUE) > 0) |>
    ungroup() -> totals
  file.path("data", "hibernation", "observations") |>
    verify_vc(
      root = raw_data,
      variables = c("sample_id", "species_id", "number")
    ) |>
    semi_join(this_species, by = c("species_id" = "id")) |>
    group_by(.data$sample_id) |>
    summarise(number = sum(.data$number), .groups = "drop") |>
    left_join(
      x = file.path("data", "hibernation", "samples") |>
        verify_vc(
          root = raw_data,
          variables = c("sample_id", "visit_id", "sublocation_id")
        ),
      by = "sample_id"
    ) |>
    left_join(
      file.path("data", "hibernation", "aggregation") |>
        verify_vc(
          root = raw_data,
          variables = c("sublocation_id", "aggregate")
        ),
      by = "sublocation_id"
    ) |>
    mutate(
      sublocation_id = ifelse(
        is.na(.data$aggregate),
        .data$sublocation_id,
        .data$aggregate
      )
    ) |>
    group_by(.data$visit_id, .data$sublocation_id) |>
    summarise(
      sample_id = min(.data$sample_id),
      number = sum(.data$number, na.rm = TRUE),
      .groups = "drop"
    ) |>
    left_join(
      x = visits |>
        filter(.data$type == "detail"),
      by = "visit_id"
    ) |>
    transmute(
      observation_id = ifelse(
        .data$level == "detail",
        .data$sample_id,
        .data$visit_id
      ),
      datafield_id = ifelse(.data$level == "detail", 1L, 2L),
      .data$location_id,
      .data$sublocation_id,
      .data$winter,
      .data$level,
      number = replace_na(.data$number, 0)
    ) |>
    group_by(.data$sublocation_id) |>
    filter(max(.data$number, na.rm = TRUE) > 0) |>
    ungroup() -> details
  details |>
    filter(.data$level == "mixed") |>
    group_by(
      .data$observation_id,
      .data$datafield_id,
      .data$location_id,
      .data$winter
    ) |>
    summarise(number = sum(.data$number), .groups = "drop") |>
    bind_rows(totals) -> relevant_totals
  if (nrow(relevant_totals) > 0) {
    relevant_totals |>
      complete(
        .data$location_id,
        winter = min(.data$winter):(max(.data$winter) + 1)
      ) -> relevant_totals
  }
  details |>
    filter(.data$level == "detail") |>
    select(-"level") |>
    complete(
      nesting(!!sym("location_id"), !!sym("sublocation_id")),
      winter = min(.data$winter):(max(.data$winter) + 1)
    ) |>
    mutate(
      observation_id = ifelse(
        is.na(.data$observation_id),
        -1000000L * .data$winter - .data$sublocation_id,
        .data$observation_id
      ),
      datafield_id = replace_na(.data$datafield_id, 3L)
    ) |>
    bind_rows(
      relevant_totals |>
        mutate(
          observation_id = ifelse(
            is.na(.data$observation_id),
            -1000000L * .data$winter - .data$location_id,
            .data$observation_id
          ),
          datafield_id = replace_na(.data$datafield_id, 4L)
        )
    ) -> full_dataset
  full_dataset |>
    filter(.data$number > 0) |>
    select("location_id", "sublocation_id", reference = "winter") |>
    inner_join(
      full_dataset |>
        filter(is.na(.data$number)),
      by = c("location_id", "sublocation_id"),
      relationship = "many-to-many"
    ) |>
    slice_min(
      abs(.data$winter - .data$reference),
      n = 1,
      with_ties = FALSE,
      by = c("location_id", "sublocation_id", "winter")
    ) |>
    filter(abs(.data$winter - .data$reference) <= n_extrapolation) |>
    bind_rows(
      full_dataset |>
        filter(!is.na(.data$number))
    ) |>
    inner_join(locations, by = "location_id") -> dataset

  display(verbose = verbose, "  preparing analysis objects")
  # calculate analysis date
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
  # prepare SPDE object
  locations |>
    select(c("X", "Y")) |>
    spde(range = c(max_dist, 0.9), sigma = c(1, 0.01)) -> spde

  # calculate polynomial coefficients for winter
  dataset |>
    distinct(.data$winter) |>
    mutate(
      winter_r = .data$winter - min(.data$winter) + 1,
      winter_l = poly(.data$winter, 3)[, 1] |>
        zapsmall(),
      winter_q = poly(.data$winter, 3)[, 2] |>
        zapsmall(),
      winter_c = poly(.data$winter, 3)[, 3] |>
        zapsmall()
    ) -> poly_winter
  dataset |>
    inner_join(poly_winter, by = "winter") |>
    transmute(
      .data$observation_id,
      .data$datafield_id,
      .data$number,
      .data$X,
      .data$Y,
      .data$location_id,
      # define the type of location
      .data$type,
      fortress = as.integer(.data$type == "fortress"),
      marl_quarry = as.integer(.data$type == "marl quarry"),
      other_large = as.integer(.data$type == "other large"),
      small = as.integer(.data$type == "small"),
      .data$sublocation_id,
      location_i = factor(.data$location_id),
      sublocation_i = factor(.data$sublocation_id),
      # centre winter to the starting year
      .data$winter,
      winter_i = .data$winter - min(.data$winter) + 1,
      .data$winter_l,
      .data$winter_q,
      .data$winter_c
    ) -> dataset
  display(verbose = verbose, "    presence")
  dataset |>
    group_by(.data$sublocation_id) |>
    mutate(
      present = as.integer(.data$number > 0),
      n_subloc = sum(.data$present == 0, na.rm = TRUE) |>
        pmin(sum(.data$present == 1, na.rm = TRUE)),
      n_subloc = ifelse(is.na(.data$sublocation_id), 0, .data$n_subloc),
      sublocation_l = ifelse(.data$n_subloc >= 3, .data$sublocation_id, NA) |>
        factor(levels = levels(.data$sublocation_i)),
      subwinter_l = ifelse(.data$n_subloc >= 3, .data$winter_l, NA_real_),
      sublocation_q = ifelse(.data$n_subloc >= 6, .data$sublocation_id, NA) |>
        factor(levels = levels(.data$sublocation_i)),
      subwinter_q = ifelse(.data$n_subloc >= 6, .data$winter_l, NA_real_),
      sublocation_c = ifelse(.data$n_subloc >= 9, .data$sublocation_id, NA) |>
        factor(levels = levels(.data$sublocation_i)),
      subwinter_c = ifelse(.data$n_subloc >= 9, .data$winter_c, NA_real_)
    ) |>
    ungroup() -> ds_present
  ds_present |>
    filter(!is.na(.data$present)) |>
    group_by(.data$location_id, .data$winter) |>
    summarise(n_loc = mean(.data$present), .groups = "drop_last") |>
    summarise(n_loc = sum(.data$n_loc), n = n()) |>
    transmute(.data$location_id, n_loc = pmin(.data$n_loc, n - .data$n_loc)) |>
    inner_join(ds_present, by = "location_id") |>
    mutate(
      location_l = ifelse(.data$n_loc >= 3, .data$location_id, NA) |>
        factor(levels = levels(.data$location_i)),
      winter_l = ifelse(.data$n_loc >= 3, .data$winter_l, NA_real_),
      location_q = ifelse(.data$n_loc >= 6, .data$location_id, NA) |>
        factor(levels = levels(.data$location_i)),
      winter_q = ifelse(.data$n_loc >= 6, .data$winter_q, NA_real_),
      location_c = ifelse(.data$n_loc >= 9, .data$location_id, NA) |>
        factor(levels = levels(.data$location_i)),
      winter_c = ifelse(.data$n_loc >= 9, .data$winter_c, NA_real_)
    ) |>
    select(-"n_loc", -"n_subloc", -"number", -"type") -> ds_present
  presence <- n2k_spde(
    formula = paste(
      "present ~ 0",
      paste(
        c(
          "fortress"[any(ds_present$fortress > 0)],
          "marl_quarry"[any(ds_present$marl_quarry > 0)],
          "other_large"[any(ds_present$other_large > 0)],
          "small"[any(ds_present$small > 0)]
        ),
        collapse = " +\n"
      ),
      "f(
        winter_i,
        model = \"rw1\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(0.15, 0.05)))
      ) +
      f(
        sublocation_i,
        model = \"iid\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.05)))
      )",
      paste(
        c(
          "      f(
        location_i,
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
      )"[length(unique(ds_present$location_id)) > 1],
          "      f(
        location_c,
        winter_c,
        model = \"iid\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.01)))
      )"[
            length(unique(ds_present$location_id)) > 1 &&
              any(!is.na(ds_present$location_c))
          ],
          "f(
        sublocation_l,
        subwinter_l,
        model = \"iid\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.01)))
      )"[any(!is.na(ds_present$sublocation_l))],
          "f(
        sublocation_q,
        subwinter_q,
        model = \"iid\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.01)))
      )"[any(!is.na(ds_present$sublocation_q))],
          "f(
        sublocation_c,
        subwinter_c,
        model = \"iid\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.01)))
      )"[any(!is.na(ds_present$sublocation_c))]
        ),
        collapse = " +\n"
      ),
      sep = " +\n "
    ),
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

  display(verbose = verbose, "    number when present")
  dataset |>
    group_by(.data$sublocation_id) |>
    mutate(
      n_subloc = ifelse(
        is.na(.data$sublocation_id),
        0,
        sum(.data$number > 0, na.rm = TRUE)
      ),
      number = ifelse(.data$number > 0, .data$number, NA_integer_),
      sublocation_l = ifelse(.data$n_subloc >= 6, .data$sublocation_id, NA) |>
        factor(levels = levels(.data$sublocation_i)),
      subwinter_l = ifelse(.data$n_subloc >= 6, .data$winter_l, 0),
      sublocation_q = ifelse(.data$n_subloc >= 12, .data$sublocation_id, NA) |>
        factor(levels = levels(.data$sublocation_i)),
      subwinter_q = ifelse(.data$n_subloc >= 12, .data$winter_l, 0),
      sublocation_c = ifelse(.data$n_subloc >= 18, .data$sublocation_id, NA) |>
        factor(levels = levels(.data$sublocation_i)),
      subwinter_c = ifelse(.data$n_subloc >= 18, .data$winter_c, 0)
    ) |>
    ungroup() -> ds_number
  ds_number |>
    filter(.data$number > 0) |>
    distinct(.data$location_id, .data$winter) |>
    count(.data$location_id, name = "n_loc") |>
    inner_join(ds_number, by = "location_id") |>
    mutate(
      location_l = ifelse(.data$n_loc >= 6, .data$location_id, NA) |>
        factor(levels = levels(.data$location_i)),
      winter_l = ifelse(.data$n_loc >= 6, .data$winter_l, 0),
      location_q = ifelse(.data$n_loc >= 12, .data$location_id, NA) |>
        factor(levels = levels(.data$location_i)),
      winter_q = ifelse(.data$n_loc >= 12, .data$winter_q, 0),
      location_c = ifelse(.data$n_loc >= 18, .data$location_id, NA) |>
        factor(levels = levels(.data$location_i)),
      winter_c = ifelse(.data$n_loc >= 18, .data$winter_c, 0)
    ) |>
    select(-"n_loc", -"n_subloc") -> ds_number
  count <- n2k_spde(
    formula = paste(
      "number ~ 0",
      paste(
        c(
          "fortress"[any(ds_number$fortress > 0)],
          "marl_quarry"[any(ds_number$marl_quarry > 0)],
          "other_large"[any(ds_number$other_large > 0)],
          "small"[any(ds_number$small > 0)]
        ),
        collapse = " +\n"
      ),

      "f(
        winter_i,
        model = \"rw1\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(0.15, 0.05)))
      )",
      "f(
        sublocation_i,
        model = \"iid\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.05)))
      )",
      paste(
        c(
          "f(
        location_i,
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
      )"[length(unique(ds_number$location_id)) > 1],
          "f(
        location_c,
        winter_c,
        model = \"iid\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(1, 0.01)))
      )"[
            length(unique(ds_number$location_id)) > 1 &&
              any(!is.na(ds_number$location_c))
          ],
          "f(
        sublocation_l,
        subwinter_l,
        model = \"iid\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(0.5, 0.01)))
      )"[any(!is.na(ds_number$sublocation_l))],
          "f(
        sublocation_q,
        subwinter_q,
        model = \"iid\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(0.5, 0.01)))
      )"[any(!is.na(ds_number$sublocation_q))],
          "f(
        sublocation_c,
        subwinter_c,
        model = \"iid\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(0.5, 0.01)))
      )"[any(!is.na(ds_number$sublocation_c))]
        ),
        collapse = " +\n"
      ),
      sep = " +\n"
    ),
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
    spde_prior = list(range = c(max_dist, 0.01), sigma = c(1, 0.01)),
    first_imported_year = min(dataset$winter),
    analysis_date = rc$when,
    last_imported_year = max(dataset$winter) - 1,
    imputation_size = 100,
    control = list(
      control.fixed = list(prec = 0.03),
      control.family = list(
        list(hyper = list(theta2 = list(initial = -20, fixed = TRUE)))
      )
    )
  )
  store_model(count, base = base, project = project, overwrite = overwrite)

  # combine the models into a hurdle model
  display(verbose = verbose, "    dependent models")
  hurdle <- n2k_hurdle_imputed(presence = presence, count = count)
  store_model(hurdle, base = base, project = project, overwrite = overwrite)

  # create the aggregate by winter and location
  aggregated_location <- n2k_aggregate(
    result_datasource_id = hurdle@AnalysisMetadata$result_datasource_id,
    scheme_id = hurdle@AnalysisMetadata$scheme_id,
    species_group_id = hurdle@AnalysisMetadata$species_group_id,
    location_group_id = hurdle@AnalysisMetadata$location_group_id,
    model_type = "aggregate imputed: sum ~ winter + type + location",
    formula = paste(
      "~winter + fortress + marl_quarry + other_large + small + location_id"
    ),
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

  # aggregate by type
  # fmt: skip
  aggregated_type <- n2k_aggregate(
    result_datasource_id =
      aggregated_location@AnalysisMetadata$result_datasource_id,
    scheme_id = aggregated_location@AnalysisMetadata$scheme_id,
    species_group_id = aggregated_location@AnalysisMetadata$species_group_id,
    location_group_id = aggregated_location@AnalysisMetadata$location_group_id,
    model_type = "aggregate imputed: sum ~ winter + type",
    formula = "~winter + fortress + marl_quarry + other_large + small",
    fun = sum,
    status = "waiting",
    parent = aggregated_location@AnalysisMetadata$file_fingerprint,
    first_imported_year =
      aggregated_location@AnalysisMetadata$first_imported_year,
    last_imported_year =
      aggregated_location@AnalysisMetadata$last_imported_year,
    duration = aggregated_location@AnalysisMetadata$duration,
    last_analysed_year =
      aggregated_location@AnalysisMetadata$last_analysed_year,
    analysis_date = aggregated_location@AnalysisMetadata$analysis_date
  )
  store_model(
    aggregated_type,
    base = base,
    project = project,
    overwrite = overwrite
  )

  # create the aggregation by winter
  # fmt: skip
  aggregated_tot <- n2k_aggregate(
    result_datasource_id =
      aggregated_type@AnalysisMetadata$result_datasource_id,
    scheme_id = aggregated_type@AnalysisMetadata$scheme_id,
    species_group_id = aggregated_type@AnalysisMetadata$species_group_id,
    location_group_id = aggregated_type@AnalysisMetadata$location_group_id,
    model_type = "aggregate imputed: sum ~ winter",
    formula = "~winter",
    fun = sum,
    status = "waiting",
    parent = aggregated_type@AnalysisMetadata$file_fingerprint,
    first_imported_year = aggregated_type@AnalysisMetadata$first_imported_year,
    last_imported_year = aggregated_type@AnalysisMetadata$last_imported_year,
    duration = aggregated_type@AnalysisMetadata$duration,
    last_analysed_year = aggregated_type@AnalysisMetadata$last_analysed_year,
    analysis_date = aggregated_type@AnalysisMetadata$analysis_date
  )
  store_model(
    aggregated_tot,
    base = base,
    project = project,
    overwrite = overwrite
  )

  display(verbose = verbose, "    index model")
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
