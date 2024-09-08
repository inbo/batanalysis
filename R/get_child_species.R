#' Returns a list of species that are children of the target species
#' @inheritParams import_raw_data
#' @inheritParams select_imputation_section
#' @export
#' @importFrom assertthat assert_that is.string noNA
#' @importFrom dplyr bind_rows filter select semi_join
#' @importFrom git2rdata read_vc
get_child_species <- function(target, species) {
  assert_that(
    inherits(target, "git_repository"), is.string(species), noNA(species)
  )
  read_vc("hibernation/species", root = target) |>
    select("code", "id", "parent") -> species_list
  species_list |>
    filter(.data$code == species) -> relevant_species
  stopifnot(
    "`species` doesn't match with a single species" =
      nrow(relevant_species) == 1
  )
  species_list |>
    semi_join(relevant_species, by = c("parent" = "id")) -> extra_species
  while (nrow(extra_species)) {
    relevant_species <- bind_rows(relevant_species, extra_species)
    species_list |>
      semi_join(extra_species, by = c("parent" = "id")) -> extra_species
  }
  return(relevant_species)
}
