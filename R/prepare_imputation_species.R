#' Prepare the imputation model for a species
#' @inheritParams import_raw_data
#' @param species the ID of the species
#' @param start The oldest date to use in the analysis
#' @export
#' @importFrom assertthat assert_that is.date
#' @importFrom dplyr filter mutate select slice_min
#' @importFrom git2rdata read_vc
#' @importFrom lubridate round_date year
#' @importFrom rlang .data
prepare_imputation_species <- function(
  target, species = 545, start = as.Date("2000-07-01")
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
    )
    read_vc("hibernation/observations", root = target) |>
      filter(.data$species_id == species)
}
