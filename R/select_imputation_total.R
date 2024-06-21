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
  target, species = "Mmysbra", start = Sys.Date() - 12 * 365
) {
  assert_that(
    inherits(target, "git_repository"), is.string(species), noNA(species),
    is.date(start), noNA(start)
  )
  read_vc("hibernation/species", root = target) |>
    select("code", "id", "parent") -> species_list
  species_list |>
    filter(.data$code == species) -> relevant_species
  species_list |>
    semi_join(relevant_species, by = c("parent" = "id")) -> extra_species
  while (nrow(extra_species)) {
    relevant_species <- bind_rows(relevant_species, extra_species)
    species_list |>
      semi_join(extra_species, by = c("parent" = "id")) -> extra_species
  }
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
        .data$species_id %in% relevant_species$id, .data$total, 0
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
    ) |>
    left_join(
      read_vc("hibernation/cluster_location", root = target),
      by = c("location_id" = "id")
    )
}
