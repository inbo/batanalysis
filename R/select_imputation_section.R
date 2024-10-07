#' Select the relevant data for the imputation model for a species by section.
#'
#' Use only the data collected using the individual based or section based
#' protocols.
#' @inheritParams import_raw_data
#' @param species the ID of the species
#' @param start The oldest date to use in the analysis
#' @param n_present Minimum number of winters in which the species is observed at
#' a location.
#' Remove locations with a number below this threshold.
#' @param n_extrapolation Maximum number of winters to extrapolate the data.
#' Only do imputations when the difference between the nearest observed winter
#' and the missing winter is below or equal to this threshold.
#' Is applied at the most detailed spatial level.
#' @export
#' @importFrom assertthat assert_that is.count is.date is.string noNA
#' @importFrom dplyr bind_rows distinct filter group_by inner_join left_join
#' mutate n_distinct select semi_join slice_min summarise transmute
#' @importFrom git2rdata read_vc
#' @importFrom lubridate round_date year
#' @importFrom purrr map2
#' @importFrom rlang .data syms !!!
#' @importFrom tidyr complete nesting replace_na unnest
select_imputation_section <- function(
  target, species = "Mmysbra", start = Sys.Date() - 12 * 365, n_present = 3,
  n_extrapolation = 5
) {
  assert_that(
    inherits(target, "git_repository"), is.string(species), noNA(species),
    is.date(start), noNA(start), is.count(n_present), noNA(n_present),
    is.count(n_extrapolation), noNA(n_extrapolation)
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
    select(-"delta", -"date") -> visits
  visits |>
    inner_join(
      read_vc("hibernation/samples", root = target), by = "visit_id"
    ) |>
    left_join(
      read_vc("hibernation/observations", root = target) |>
        mutate(
          number = ifelse(
            .data$species_id %in% relevant_species$id, .data$number, 0L
          )
        ) |>
        group_by(.data$sample_id) |>
        summarise(number = sum(.data$number)),
      by = "sample_id"
    ) |>
    mutate(number = replace_na(.data$number, 0L)) |>
    complete(
      winter = min(.data$winter):max(.data$winter),
      nesting(!!!syms(c("location_id", "sublocation_id")))
    ) -> observations
  # remove locations with less than n_present observed winters
  observations |>
    filter(.data$number > 0) |>
    group_by(.data$location_id) |>
    summarise(winters = n_distinct(.data$winter)) |>
    filter(.data$winters >= n_present) |>
    semi_join(x = observations, by = "location_id") -> observations
  # remove sublocation winter combinations with to much extrapolation
  observations |>
    filter(is.na(.data$number)) |>
    select("winter", "sublocation_id") |>
    inner_join(
      observations |>
        filter(.data$number > 0) |>
        select(observed = "winter", "sublocation_id"),
      by = "sublocation_id", relationship = "many-to-many"
    ) |>
    mutate(delta = abs(.data$winter - .data$observed)) |>
    filter(.data$delta <= n_extrapolation) |>
    distinct(.data$winter, .data$sublocation_id) |>
    semi_join(x = observations, by = c("winter", "sublocation_id")) |>
    bind_rows(
      observations |>
        filter(!is.na(.data$number))
    )
}
