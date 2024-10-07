#' Prepare analysis models for all species
#' @inheritParams prepare_analysis_model_species
#' @inheritParams prepare_analysis_data
#' @inheritParams n2kanalysis::store_model
#' @export
#' @importFrom assertthat assert_that
#' @importFrom n2kanalysis manifest_yaml_to_bash n2k_manifest
#' store_manifest_yaml
#' @importFrom purrr map_dfr
prepare_analysis_model <- function(
  analysis_data, base, max_dist = 10, project = "batanalysis",
  overwrite = FALSE
) {
  assert_that(inherits(analysis_data, "git_repository"))
  dirname(analysis_data$path) |>
    file.path("hibernation") |>
    list.dirs(full.names = FALSE) -> species
  species <- species[species != ""]
  map_dfr(
    species,
    ~prepare_analysis_model_species(
      species = .x, base = base, max_dist = max_dist, project = project,
      overwrite = overwrite, analysis_data = analysis_data
    )
  ) |>
    n2k_manifest() |>
    store_manifest_yaml(
      base = base, project = project, docker = "inbobmk/rn2k:dev-0.10",
      dependencies = c(
        "inbo/multimput@v0.2.14", "inbo/n2khelper@v0.5.0",
        "inbo/n2kanalysis@spde"
      )
    ) |>
    basename() |>
    manifest_yaml_to_bash(base = base, project = project, shutdown = TRUE)
}
