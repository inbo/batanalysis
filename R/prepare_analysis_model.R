#' Prepare analysis models for all species
#' @inheritParams prepare_analysis_model_species
#' @inheritParams prepare_analysis_data
#' @inheritParams n2kanalysis::store_model
#' @export
#' @importFrom assertthat assert_that
#' @importFrom n2kanalysis n2k_manifest store_manifest
#' @importFrom purrr map_dfr
prepare_analysis_model <- function(
  analysis_data, base, max_dist = 10, project = "batanalysis",
  overwrite = FALSE
) {
  assert_that(inherits(analysis_data, "git_repository"))
  dirname(analysis_data$path) |>
    file.path("hibernation") |>
    list.dirs(full.names = FALSE) -> species
  species <- species[species == "mdau"]
  map_dfr(
    species,
    ~prepare_analysis_model_species(
      species = .x, base = base, max_dist = max_dist, project = project,
      overwrite = overwrite, analysis_data = analysis_data
    )
  ) |>
    n2k_manifest() |>
    store_manifest(base = base, project = project)
}
