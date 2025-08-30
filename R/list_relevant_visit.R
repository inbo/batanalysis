#' List relevant visits for hibernation analysis
#'
#' We apply the following rules to determine if a visit is relevant location by
#' location:
#' A relevant visit is no older than `n_winter` winters.
#' When a location is divided into sublocations, we count the number of visits
#' with observations at the sublocation level.
#' If there are at least three such winters with visits where bats are observed,
#' and the number of winters with detailed observations is larger than the
#' number with only total counts, only the visits with detailed observation are
#' relevant.
#' Otherwise, we look at all visits, and if there are at least three winters
#' with visits where bats are observed, then all visits within the last
#' `n_winter` are relevant.
#' @inheritParams import_raw_data
#' @param n_winter The number of winters to consider for relevance.
#' @return A data frame with the relevant visits and their relevance status.
#' @export
#' @importFrom assertthat assert_that is.count
#' @importFrom dplyr bind_rows count filter group_by inner_join left_join mutate
#'  select slice_min summarise transmute
#' @importFrom git2rdata verify_vc
#' @importFrom lubridate round_date year
#' @importFrom rlang .data
#' @export
list_relevant_visit <- function(data_root, n_winter = 24) {
  assert_that(is.count(n_winter))
  file.path("data", "hibernation", "locations") |>
    verify_vc(
      root = data_root,
      variables = c("id", "name", "parent_id")
    ) -> locations
  file.path("data", "hibernation", "aggregation") |>
    verify_vc(
      root = data_root,
      variables = c("sublocation_id", "aggregate")
    ) -> aggregation
  locations |>
    filter(.data$parent_id > 0) |>
    left_join(aggregation, by = c("id" = "sublocation_id")) |>
    mutate(id = ifelse(is.na(.data$aggregate), .data$id, .data$aggregate)) |>
    select("parent_id", sublocation_id = "id") -> sublocations
  locations |>
    filter(.data$parent_id < 0) |>
    left_join(
      sublocations |>
        count(.data$parent_id, name = "n_sublocations"),
      by = c("id" = "parent_id")
    ) |>
    transmute(
      location_id = .data$id,
      n_sublocations = replace_na(.data$n_sublocations, 1)
    ) -> main_locations
  file.path("data", "hibernation", "species") |>
    verify_vc(
      root = data_root,
      variables = c("id")
    ) -> species
  file.path("data", "hibernation", "visits") |>
    verify_vc(
      root = data_root,
      variables = c("visit_id", "location_id", "date")
    ) |>
    left_join(
      file.path("data", "hibernation", "totals") |>
        verify_vc(
          root = data_root,
          variables = c("visit_id", "total", "species_id")
        ) |>
        semi_join(species, by = c("species_id" = "id")) |>
        group_by(.data$visit_id) |>
        summarise(total = sum(.data$total)),
      by = "visit_id"
    ) -> visits
  file.path("data", "hibernation", "samples") |>
    verify_vc(
      root = data_root,
      variables = c("visit_id", "sample_id", "sublocation_id")
    ) -> samples
  samples |>
    left_join(aggregation, by = "sublocation_id") |>
    mutate(
      sublocation_id = ifelse(
        is.na(.data$aggregate),
        .data$sublocation_id,
        .data$aggregate
      )
    ) |>
    left_join(
      file.path("data", "hibernation", "observations") |>
        verify_vc(
          root = data_root,
          variables = c("sample_id", "number", "species_id")
        ) |>
        semi_join(species, by = c("species_id" = "id")),
      by = "sample_id"
    ) |>
    group_by(.data$visit_id, .data$sublocation_id) |>
    summarise(total = sum(.data$number, na.rm = TRUE), .groups = "drop_last") |>
    summarise(
      surveyed_sublocations = n(),
      detailed_total = sum(.data$total)
    ) |>
    left_join(x = visits, by = "visit_id") |>
    mutate(
      total = ifelse(is.na(.data$total), .data$detailed_total, .data$total) |>
        replace_na(0),
      winter = round_date(.data$date, "year") |>
        year()
    ) |>
    select(-"detailed_total") -> visits
  visits |>
    filter(max(.data$winter) - .data$winter + 1 > n_winter) |>
    mutate(relevant = FALSE) -> old_visits
  visits |>
    anti_join(old_visits, by = "visit_id") -> visits
  visits |>
    mutate(
      delta = paste0(.data$winter, "-1-14 18:0:0") |>
        as.POSIXct(),
      delta = abs(as.POSIXct(.data$date) - .data$delta) |>
        as.vector(),
      delta = .data$delta + is.na(.data$surveyed_sublocations) * 1e6
    ) |>
    group_by(.data$location_id, .data$winter) |>
    slice_min(.data$delta, n = 1, with_ties = FALSE) |>
    ungroup() -> candidate
  visits |>
    anti_join(candidate, by = "visit_id") |>
    mutate(relevant = FALSE) |>
    bind_rows(old_visits) -> old_visits
  candidate |>
    filter(.data$total > 0) |>
    group_by(.data$location_id) |>
    summarise(
      detailed = sum(!is.na(.data$surveyed_sublocations)),
      global = sum(is.na(.data$surveyed_sublocations))
    ) |>
    transmute(
      .data$location_id,
      relevant = .data$detailed >= 3 | .data$detailed + .data$global >= 3
    ) |>
    left_join(x = candidate, by = "location_id") |>
    select(-"delta") |>
    mutate(relevant = replace_na(.data$relevant, FALSE)) |>
    bind_rows(old_visits) |>
    inner_join(main_locations, by = "location_id")
}
