#' Prepare the analysis models for a single species
#' @inheritParams n2kanalysis::store_model
#' @inheritParams n2kanalysis::display
#' @inheritParams prepare_analysis_model_species
#' @param raw_data The git repository with the raw data
#' @param species the code of the species
#' @param visits The output of `visit_type()`
#' @export
#' @importFrom dplyr bind_rows distinct filter group_by inner_join left_join
#' @importFrom dplyr mutate select semi_join slice_max transmute
#' @importFrom git2rdata recent_commit verify_vc
#' @importFrom n2kanalysis display n2k_aggregate status<- store_model
#' @importFrom rlang .data
#' @importFrom stringr str_detect
#' @importFrom tidyr replace_na
prepare_analysis_no_model_species <- function(
  raw_data,
  species,
  base,
  project = "batanalysis",
  visits,
  overwrite = FALSE,
  verbose = TRUE
) {
  # prepare species information
  file.path("data", "hibernation", "species") |>
    verify_vc(root = raw_data, variables = c("id", "scientific_name")) |>
    filter(.data$code == species) |>
    select("id", "scientific_name", "code")  -> this_species
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
          variables = c("id", "longitude", "latitude", "parent_id", "province")
        ),
      by = c("location_id" = "id")
    ) |>
    left_join(
      file.path("data", "hibernation", "location_types") |>
        verify_vc(root = raw_data, variables = c("location_id", "type")),
      by = "location_id"
    ) |>
    transmute(
      .data$location_id,
      .data$detailed,
      .data$province,
      type = replace_na(.data$type, "small")
     ) -> locations
  locations |>
    distinct(.data$province) |>
    filter(
      str_detect(.data$province, "Brussel", negate = TRUE),
      str_detect(.data$province, "Waals", negate = TRUE)
    ) -> vlaanderen

  display(verbose = verbose, "  preparing data")
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
  file.path("data", "hibernation", "totals") |>
    verify_vc(
      root = raw_data,
      variables = c("visit_id", "species_id", "total")
    ) |>
    semi_join(this_species, by = c("species_id" = "id")) |>
    group_by(.data$visit_id) |>
    summarise(number = sum(.data$total), .groups = "drop") -> raw_totals
  raw_totals |>
    inner_join(
      x = visits |>
        filter(
          .data$type == "total",
          is.na(.data$level) | .data$level != "detail"
        ),
      by = "visit_id"
    ) |>
    select("location_id", "winter", "number") -> totals
  file.path("data", "hibernation", "observations") |>
    verify_vc(
      root = raw_data,
      variables = c("sample_id", "species_id", "number")
    ) |>
    semi_join(this_species, by = c("species_id" = "id")) |>
    group_by(.data$sample_id) |>
    summarise(number = sum(.data$number), .groups = "drop") |>
    inner_join(
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
    inner_join(
      x = visits |>
        filter(.data$type == "detail"),
      by = "visit_id"
    ) |>
    transmute(
      .data$location_id,
      .data$sublocation_id,
      .data$winter,
      .data$level,
      .data$number
    ) -> details
  details |>
    filter(.data$level == "mixed") |>
    group_by(.data$location_id, .data$winter) |>
    summarise(number = sum(.data$number), .groups = "drop") |>
    bind_rows(totals) -> relevant_totals

  # totals in detailed timeseries to be used as a minimum for the imputed totals
  details |>
    filter(.data$level == "detail") |>
    bind_rows(
      visits |>
        filter(.data$level == "detail") |>
        semi_join(x = visits, by = "location_id") |>
        filter(.data$type == "total") |>
        inner_join(raw_totals, by = "visit_id") |>
        transmute(
          .data$location_id,
          .data$winter,
          .data$number
        ),
      relevant_totals
    ) |>
    inner_join(locations, by = "location_id") |>
    select(-"detailed") -> dataset


  display(verbose = verbose, "  preparing analysis objects")
  # calculate analysis date

  # create the aggregate by winter and location
  aggregated_location <- n2k_aggregate(
    model_type = "aggregate imputed: sum ~ winter + location",
    formula = paste(
      "~winter + province + type + location_id",
      sep = "+"
    ),
    fun = sum,
    status = "waiting",
    parent = "",
    result_datasource_id = "git",
    scheme_id = "hibernating bats",
    species_group_id = species,
    location_group_id = "Flanders",
    seed = 19911204,
    first_imported_year = min(dataset$winter),
    analysis_date = rc$when,
    last_imported_year = max(dataset$winter) - 1
  )

  n_imp <- 100
  aggregated_location@RawImputed <- new(
    "aggregatedImputed",
     Covariate = dataset |>
            select(-"number") |>
            as.data.frame(),
          Imputation = rep(dataset$number, n_imp) |>
  matrix(ncol = n_imp) |>
  `colnames<-`(sprintf(fmt = "Imputation%04i", seq_len(n_imp)))
  )
  status(aggregated_location) <- "new"
  store_model(
    aggregated_location,
    base = base,
    project = project,
    overwrite = overwrite
  )

  # create the aggregate by winter, province and type
  aggregated_province_type <- n2k_aggregate(
    result_datasource_id = aggregated_location@AnalysisMetadata$result_datasource_id,
    scheme_id = aggregated_location@AnalysisMetadata$scheme_id,
    species_group_id = aggregated_location@AnalysisMetadata$species_group_id,
    location_group_id = aggregated_location@AnalysisMetadata$location_group_id,
    model_type = "aggregate imputed: sum ~ winter + province + type",
    formula = paste(
      "~winter + province + type"
    ),
    fun = sum,
    status = "waiting",
    parent = aggregated_location@AnalysisMetadata$file_fingerprint,
    first_imported_year = aggregated_location@AnalysisMetadata$first_imported_year,
    last_imported_year = aggregated_location@AnalysisMetadata$last_imported_year,
    duration = aggregated_location@AnalysisMetadata$duration,
    last_analysed_year = aggregated_location@AnalysisMetadata$last_analysed_year,
    analysis_date = aggregated_location@AnalysisMetadata$analysis_date
  )
  store_model(
    aggregated_province_type,
    base = base,
    project = project,
    overwrite = overwrite
  )

  # aggregate by type
  # fmt: skip
  aggregated_type <- n2k_aggregate(
    result_datasource_id =
      aggregated_province_type@AnalysisMetadata$result_datasource_id,
    scheme_id = aggregated_province_type@AnalysisMetadata$scheme_id,
    species_group_id =
      aggregated_province_type@AnalysisMetadata$species_group_id,
    location_group_id =
      aggregated_province_type@AnalysisMetadata$location_group_id,
    model_type = "aggregate imputed: sum ~ winter + type",
    formula = "~winter + type",
    fun = sum,
    join = vlaanderen,
    status = "waiting",
    parent = aggregated_province_type@AnalysisMetadata$file_fingerprint,
    first_imported_year =
      aggregated_province_type@AnalysisMetadata$first_imported_year,
    last_imported_year =
      aggregated_province_type@AnalysisMetadata$last_imported_year,
    duration = aggregated_province_type@AnalysisMetadata$duration,
    last_analysed_year =
      aggregated_province_type@AnalysisMetadata$last_analysed_year,
    analysis_date = aggregated_province_type@AnalysisMetadata$analysis_date
  )
  store_model(
    aggregated_type,
    base = base,
    project = project,
    overwrite = overwrite
  )

  # aggregate by province
  # fmt: skip
  aggregated_province <- n2k_aggregate(
    result_datasource_id =
      aggregated_province_type@AnalysisMetadata$result_datasource_id,
    scheme_id = aggregated_province_type@AnalysisMetadata$scheme_id,
    species_group_id =
      aggregated_province_type@AnalysisMetadata$species_group_id,
    location_group_id =
      aggregated_province_type@AnalysisMetadata$location_group_id,
    model_type = "aggregate imputed: sum ~ winter + province",
    formula = "~winter + province",
    fun = sum,
    status = "waiting",
    parent = aggregated_province_type@AnalysisMetadata$file_fingerprint,
    first_imported_year =
      aggregated_province_type@AnalysisMetadata$first_imported_year,
    last_imported_year =
      aggregated_province_type@AnalysisMetadata$last_imported_year,
    duration = aggregated_province_type@AnalysisMetadata$duration,
    last_analysed_year =
      aggregated_province_type@AnalysisMetadata$last_analysed_year,
    analysis_date = aggregated_province_type@AnalysisMetadata$analysis_date
  )
  store_model(
    aggregated_province,
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

  bind_rows(
    aggregated_location@AnalysisRelation |>
      mutate(type = "location"),
    aggregated_province_type@AnalysisRelation |>
      mutate(type = "province_type"),
    aggregated_province@AnalysisRelation |>
      mutate(type = "province"),
    aggregated_type@AnalysisRelation |>
      mutate(type = "type"),
    aggregated_tot@AnalysisRelation |>
      mutate(type = "tot")
  ) |>
    select(fingerprint = "analysis", parent = "parent_analysis", "type")
}
