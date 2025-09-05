#' Count the unsurveyed sublocations per location
#' @param data_repo Character; path to the root of the data repository
#' @return A tibble with columns location_id, location_name, sublocation_id,
#' sublocation, n (number of visits where the sublocation was not surveyed)
#' @importFrom dplyr anti_join count distinct filter inner_join mutate select
#' semi_join
#' @importFrom git2rdata verify_vc
#' @importFrom rlang .data
#' @export
unsurveyed_sublocations <- function(
  data_repo = "."
) {
  file.path("data", "hibernation", "locations") |>
    verify_vc(
      root = data_repo,
      variables = c("id", "name", "parent_id")
    ) |>
    select("id", "name", "parent_id") -> locations
  locations |>
    filter(.data$parent_id > 0) -> sublocations
  locations |>
    semi_join(sublocations, by = c("id" = "parent_id")) |>
    select(-"parent_id", location_name = "name") -> locations
  file.path("data", "hibernation", "samples") |>
    verify_vc(root = data_repo, variables = c("visit_id", "sublocation_id")) |>
    distinct(.data$visit_id, .data$sublocation_id) -> samples
  file.path("data", "hibernation", "visits") |>
    verify_vc(root = data_repo, variables = c("visit_id", "location_id")) |>
    select(-"date") |>
    semi_join(samples, by = "visit_id") |>
    inner_join(
      sublocations,
      by = c("location_id" = "parent_id"),
      relationship = "many-to-many"
    ) |>
    anti_join(samples, by = c("visit_id", "id" = "sublocation_id")) |>
    count(
      .data$location_id,
      sublocation_id = .data$id,
      sublocation = .data$name
    ) |>
    inner_join(locations, by = c("location_id" = "id"))
}
