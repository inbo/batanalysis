#' Prepare analysis models for all species
#' @inheritParams prepare_analysis_model_species
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
  n_winter = 4,
  n_present = 3,
  n_extrapolation = 5,
  max_delta = 10,
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
        max_dist = max_dist,
        project = project,
        overwrite = overwrite,
        raw_data = raw_data,
        start_winter = start_winter,
        n_winter = n_winter,
        n_present = n_present,
        n_extrapolation = n_extrapolation,
        max_delta = max_delta,
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
