#' Prepare analysis models for all species
#' @inheritParams prepare_analysis_model_species
#' @inheritParams visit_type
#' @inheritParams n2kanalysis::store_model
#' @inheritParams n2kanalysis::display
#' @export
#' @importFrom n2kanalysis manifest_yaml_to_bash n2k_manifest
#' store_manifest_yaml
#' @importFrom purrr map_dfr
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
  visit_type(
    raw_data = raw_data,
    start_winter = start_winter,
    max_delta = max_delta
  ) -> visits
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
        visits = visits,
        verbose = verbose
      )
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
