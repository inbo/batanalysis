#' Select the relevant data for the imputation model for a species
#' @inheritParams import_raw_data
#' @param species the ID of the species
#' @param start The oldest date to use in the analysis
#' @param n_winter Minimum number of winters in which the species is observed at
#' a location.
#' Remove locations with a number below this threshold.
#' @export
#' @importFrom assertthat assert_that is.date
#' @importFrom dplyr distinct filter group_by inner_join left_join mutate n
#' select semi_join slice_min summarise transmute
#' @importFrom git2rdata read_vc
#' @importFrom lubridate round_date year
#' @importFrom purrr map2
#' @importFrom rlang .data syms !!!
#' @importFrom tidyr complete nesting unnest
select_imputation_species <- function(
  target, species = 545, start = as.Date("2000-07-01"), n_winter = 2
) {
  assert_that(inherits(target, "git_repository"), is.date(start))
  read_vc("hibernation/visits", root = target) |>
    filter(.data$date >= start) |>
    mutate(
      winter = round_date(.data$date, unit = "year"),
      delta = abs(
        as.POSIXct(.data$winter) + (14 * 24 + 9) * 3600 - as.POSIXct(.data$date)
      ),
      winter = year(.data$winter)
    ) |>
    slice_min(.data$delta, n = 1, by = c("location_id", "winter")) |>
    select(-"delta", -"date") -> visits
  visits |>
    inner_join(
      read_vc("hibernation/samples", root = target), by = "visit_id"
    ) |>
    left_join(
      read_vc("hibernation/observations", root = target) |>
        mutate(number = ifelse(.data$species_id == species, .data$number, 0)) |>
        group_by(.data$sample_id) |>
        summarise(number = sum(.data$number)),
      by = "sample_id"
    ) |>
    complete(
      .data$winter, nesting(!!!syms(c("location_id", "sublocation_id")))
    ) -> observations
  observations |>
    filter(!is.na(.data$number)) |>
    distinct(.data$location_id, .data$winter) |>
    group_by(.data$location_id) |>
    summarise(
      winters = n(), first = min(.data$winter), last = max(.data$winter)
    ) |>
    filter(.data$winters >= n_winter) |>
    transmute(.data$location_id, winter = map2(.data$first, .data$last, seq)) |>
    unnest("winter") |>
    semi_join(x = observations, by = c("location_id", "winter")) -> observations
  observations |>
    filter(.data$number > 0) |>
    distinct(.data$location_id, .data$winter) |>
    group_by(.data$location_id) |>
    summarise(
      winters = n(), first = min(.data$winter), last = max(.data$winter)
    ) |>
    filter(.data$winters >= n_winter) |>
    transmute(.data$location_id, winter = map2(.data$first, .data$last, seq)) |>
    unnest("winter") |>
    semi_join(x = observations, by = c("location_id", "winter")) -> observations
  observations |>
    distinct(.data$location_id, .data$sublocation_id) |>
    inner_join(
      read_vc("hibernation/cluster_location", root = target),
      by = c("location_id" = "id")
    )
}
