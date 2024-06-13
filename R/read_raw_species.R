#' Import the species list
#' @inheritParams import_raw_data
#' @export
#' @importFrom assertthat assert_that
#' @importFrom dplyr inner_join left_join select
#' @importFrom DBI dbGetQuery
#' @importFrom rlang .data
read_raw_species <- function(origin) {
  assert_that(inherits(origin, "Microsoft SQL Server"))
  "SELECT s.id, s.name, s.scientific_name, s.code
FROM staging_Meetnetten.projects_project AS p
INNER JOIN staging_Meetnetten.projects_projectspecies AS psp
  ON psp.project_id = p.id
INNER JOIN staging_Meetnetten.species_species AS s ON s.id = psp.species_id
WHERE p.name = 'Vleermuizen - Wintertellingen' AND psp.is_primary = 'TRUE'" |>
    dbGetQuery(conn = origin) -> species
  data.frame(
    code = c(
      "Mmys", "Mmysbra", "Mbra", "Mbec", "Mnat", "Pippip", "Paur", "Paus",
      "Pauraus", "Mema", "Eser", "Mdas", "Bbar", "Pipnat", "Vmur", "Mmyo",
      "Mdau", "Mspec", "Pipspec", "Rfer", "Rhip", "Nlei", "Nlas", "Nspec",
      "Nnoc", "Malc"
    ),
    parent = c(
      "Mmysbra", "Mspec", "Mmysbra", "Mspec", "Mspec", "Pipspec", "Pauraus",
      "Pauraus", "Cspec", "Mspec", "Cspec", "Mspec", "Cspec", "Cspec", "Cspec",
      "Mspec", "Mspec", "Cspec", "Cspec", "Cspec", "Cspec", "Nspec", "Nspec",
      "Cspec", "Nspec", "Mspec"
    )
  ) |>
    inner_join(species, by = c("parent" = "code")) |>
    select("code", parent = "id") |>
    left_join(x = species, by = "code")
}
