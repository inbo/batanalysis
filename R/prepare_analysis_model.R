#' Prepare analysis models for all species
#' @inheritParams prepare_analysis_model_species
#' @inheritParams visit_type
#' @inheritParams n2kanalysis::store_model
#' @inheritParams n2kanalysis::display
#' @export
#' @importFrom n2kanalysis manifest_yaml_to_bash n2k_aggregate_composite
#' @importFrom n2kanalysis n2k_manifest
#' @importFrom n2kanalysis store_manifest_yaml
#' @importFrom purrr map_dfr pmap
#' @importFrom tidyr nest unnest
prepare_analysis_model <- function(
  raw_data,
  start_winter = as.integer(format(Sys.Date(), "%Y")) - 23,
  base,
  project = "batanalysis",
  max_delta = 10,
  n_extrapolation = 24,
  max_dist = 10,
  overwrite = FALSE,
  verbose = TRUE
) {
  c(
    "Mbec",
    "Mdas",
    "Mdau",
    "Mema",
    "Mmyo",
    "Mmysbra",
    "Mnat",
    "Pauraus",
    "Pipspec"
  ) |>
    map_dfr(
      ~ prepare_analysis_model_species(
        species = .x,
        base = base,
        n_extrapolation = n_extrapolation,
        max_dist = max_dist,
        project = project,
        overwrite = overwrite,
        raw_data = raw_data,
        visits = visit_type(
          raw_data = raw_data,
          start_winter = start_winter,
          max_delta = max_delta
        ),
        verbose = verbose
      )
    ) -> model
  c(
    "Mspec",
    "Cspec",
    "Eser",
    "Rfer"
  ) |>
    map_dfr(
      ~ prepare_analysis_no_model_species(
        species = .x,
        base = base,
        project = project,
        overwrite = overwrite,
        raw_data = raw_data,
        visits = visit_type(
          raw_data = raw_data,
          start_winter = start_winter,
          max_delta = max_delta
        ),
        verbose = verbose
      )
    ) -> no_model
  no_model |>
    mutate(
      parent = ifelse(.data$parent == "", NA_character_, .data$parent)
    ) |>
    bind_rows(model) -> all_model

  display("Sum over all species", verbose = verbose)
  tail(all_model$fingerprint, 1) |>
    read_model(base = base, project = project) |>
    slot("AnalysisMetadata") -> metadata
  all_model |>
    filter(.data$type %in% c("location", "province", "type", "tot")) |>
    distinct(
      .data$type,
      parent_analysis = .data$fingerprint,
      parentstatus_fingerprint = .data$status,
      parent_status = "waiting"
    ) |>
    nest(.by = "type") |>
    mutate(
      fingerprint = n2k_aggregate_composite(
        parent_status = .data$data[[1]],
        status = "waiting",
        result_datasource_id = metadata$result_datasource_id,
        scheme_id = metadata$scheme_id,
        species_group_id = "totaal",
        location_group_id = metadata$location_group_id,
        model_type = "totaal",
        first_imported_year = metadata$first_imported_year,
        last_imported_year = metadata$last_imported_year,
        duration = metadata$duration,
        last_analysed_year = metadata$last_analysed_year,
        analysis_date = metadata$analysis_date,
        fun = function(...) {
          sum(..., na.rm = TRUE)
        },
        seed = metadata$seed,
        formula = "count ~ 1"
      ) |>
        slot("AnalysisMetadata") |>
        select(fingerprint = "file_fingerprint", "status_fingerprint") |>
        list()
    ) |>
    unnest("fingerprint") |>
    unnest("data") |>
    select(
      "fingerprint",
      "status_fingerprint",
      parent = "parent_analysis"
    ) -> total_model

  display("Ratio of species to sum of all species", verbose = verbose)
  calc_ratio <- function(i) {
    stopifnot(length(i) == 2)
    min(i) / max(i)
  }
  all_model |>
    filter(.data$type %in% c("location", "province", "type", "tot")) |>
    select("type", parent = "fingerprint", "status") |>
    inner_join(total_model, by = "parent") |>
    transmute(
      model = pmap(
        list(
          .data$parent,
          .data$fingerprint,
          .data$status,
          .data$status_fingerprint,
          .data$type
        ),
        function(parent, fingerprint, status, status_fingerprint, type) {
          n2k_aggregate_composite(
            parent_status = data.frame(
              parent_analysis = c(parent, fingerprint),
              parentstatus_fingerprint = c(status, status_fingerprint),
              parent_status = "waiting"
            ),
            fun = calc_ratio,
            result_datasource_id = metadata$result_datasource_id,
            scheme_id = metadata$scheme_id,
            species_group_id = metadata$species_group_id,
            location_group_id = metadata$location_group_id,
            model_type = "ratio to sum of all species",
            formula = paste0("~", type),
            status = "waiting",
            first_imported_year = metadata$first_imported_year,
            last_imported_year = metadata$last_imported_year,
            duration = metadata$duration,
            last_analysed_year = metadata$last_analysed_year,
            analysis_date = metadata$analysis_date
          ) |>
            slot("AnalysisRelation")
        }
      )
    ) |>
    unnest("model") |>
    select(fingerprint = "analysis", parent = "parent_analysis") -> rel_tot

  display("Creating manifest", verbose = verbose)
  rel_tot |>
    bind_rows(
      total_model |>
        select("fingerprint", "parent"),
      all_model |>
        select("fingerprint", "parent")
    ) |>
    n2k_manifest() |>
    store_manifest_yaml(
      base = base,
      project = project,
      docker = "inbobmk/rn2k:dev-0.10",
      dependencies = c(
        "inbo/multimput@hotfix",
        "inbo/n2khelper@0.5.1",
        "inbo/n2kanalysis@0.4.1"
      )
    ) |>
    basename() |>
    manifest_yaml_to_bash(base = base, project = project, shutdown = TRUE)
}
