#' Select the relevant data for the imputation model for a species by total.
#'
#' Use only the data collected using the totals based protocol.
#' @inheritParams import_raw_data
#' @inheritParams select_imputation_section
#' @export
#' @importFrom assertthat assert_that is.date is.string noNA
#' @importFrom dplyr bind_rows count filter group_by inner_join left_join
#' mutate select semi_join slice_min summarise transmute
#' @importFrom git2rdata read_vc
#' @importFrom lubridate round_date year
select_imputation_total <- function(
  target, species = "Mmysbra", start = Sys.Date() - 12 * 365,
  n_present = 3, n_extrapolation = 5
) {
  assert_that(
    inherits(target, "git_repository"), is.string(species), noNA(species),
    is.date(start), noNA(start), is.count(n_extrapolation),
    noNA(n_extrapolation), is.count(n_present), noNA(n_present)
  )
  relevant_species <- get_child_species(target = target, species = species)
  read_vc("hibernation/visits", root = target) |>
    filter(.data$date >= start) |>
    mutate(
      winter = round_date(.data$date, unit = "year"),
      delta = abs(
        as.POSIXct(.data$winter) + (14 * 24 + 9) * 3600 - as.POSIXct(.data$date)
      ),
      winter = year(.data$winter) |>
        as.integer()
    ) |>
    slice_min(.data$delta, n = 1, by = c("location_id", "winter")) |>
    select(-"delta", -"date") |>
    complete(
      .data$location_id, winter = min(.data$winter):max(.data$winter)
    ) -> visits
  visits |>
    inner_join(
      read_vc("hibernation/totals", root = target), by = "visit_id"
    ) |>
    mutate(
      total = ifelse(
        .data$species_id %in% relevant_species$id, .data$total, 0L
      )
    ) |>
    group_by(.data$location_id, .data$winter) |>
    summarise(total = sum(.data$total), .groups = "drop") |>
    complete(
      .data$location_id, winter = min(.data$winter):max(.data$winter)
    ) |>
    left_join(
      read_vc("hibernation/locations", root = target) |>
        filter(.data$parent_id > 0) |>
        count(location_id = .data$parent_id),
      by = "location_id"
    ) |>
    transmute(
      .data$location_id, .data$winter,
      minimum = ifelse(!is.na(.data$n), .data$total, NA),
      total = ifelse(is.na(.data$n), .data$total, NA)
    ) -> observations
  # remove locations with less than n_present observed winters
  observations |>
    filter(.data$minimum > 0 | .data$total > 0) |>
    count(.data$location_id) |>
    filter(.data$n >= n_present) |>
    semi_join(x = observations, by = "location_id") -> observations
  # remove location winter combinations with to much extrapolation
  observations |>
    filter(is.na(.data$total), is.na(.data$minimum)) |>
    select("winter", "location_id") |>
    inner_join(
      observations |>
        filter(.data$total > 0 | .data$minimum > 0) |>
        select(observed = "winter", "location_id"),
      by = "location_id", relationship = "many-to-many"
    ) |>
    mutate(delta = abs(.data$winter - .data$observed)) |>
    filter(.data$delta <= n_extrapolation) |>
    distinct(.data$winter, .data$location_id) |>
    semi_join(x = observations, by = c("winter", "location_id")) |>
    bind_rows(
      observations |>
        filter(!is.na(.data$total) | !is.na(.data$minimum))
    )
}
