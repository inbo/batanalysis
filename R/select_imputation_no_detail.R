#' Select imputation data from non-detailed locations
#' @inheritParams prepare_analysis_model_species
#' @inheritParams select_imputation_detail
#' @export
#' @importFrom git2rdata verify_vc
#' @importFrom dplyr filter group_by left_join inner_join mutate n select
#' semi_join slice_min ungroup
#' @importFrom lubridate round_date year
#' @importFrom tidyr complete replace_na
#' @importFrom rlang .data
select_imputation_no_detail <- function(
  locations,
  raw_data,
  this_species,
  start,
  n_winter = 4,
  n_present = 3,
  n_extrapolation = 5
) {
  # read visits from file
  file.path("data", "hibernation", "visits") |>
    verify_vc(root = raw_data, variables = c("visit_id", "location_id")) |>
    filter(.data$date >= start) |>
    # keep visits from non-detailed locations only
    semi_join(
      locations |>
        filter(!.data$detailed),
      by = "location_id"
    ) |>
    # determine closest visit to Jan 15
    mutate(
      winter = round_date(.data$date, "year") |>
        year(),
      delta = paste0(.data$winter, "-1-15 12:0:0") |>
        as.POSIXct(),
      delta = difftime(.data$date, .data$delta, units = "days") |>
        as.numeric()
    ) |>
    slice_min(
      abs(.data$delta),
      with_ties = FALSE,
      by = c("location_id", "winter")
    ) |>
    # merge numbers
    left_join(
      # read totals from file
      file.path("data", "hibernation", "totals") |>
        verify_vc(
          root = raw_data,
          variables = c("visit_id", "species_id", "total")
        ) |>
        # keep only totals for this species
        semi_join(this_species, by = c("species_id" = "id")) |>
        # aggregate in case of multiple species
        group_by(.data$visit_id) |>
        summarise(number = sum(.data$total), .groups = "drop"),
      by = "visit_id"
    ) |>
    select(-"date", -"delta", observation_id = "visit_id") |>
    # missing totals indicate no observation of this species
    # thus setting the number to zero
    mutate(number = replace_na(.data$number, 0)) |>
    # keep only locations with at least `n_winter` visits
    # and at least `n_present` visits with observed bats
    group_by(.data$location_id) |>
    filter(n() >= n_winter, sum(.data$number > 0) >= n_present) |>
    ungroup() -> observed
  observed |>
    # add unvisited winters per location
    complete(
      .data$location_id,
      winter = min(.data$winter):(max(.data$winter) + 1)
    ) |>
    inner_join(
      # find closest observed winter with bats for each unvisited winter
      observed |>
        filter(.data$number > 0) |>
        select("location_id", "present" = "winter"),
      by = "location_id",
      relationship = "many-to-many"
    ) |>
    # keep only closest observed winter with bats
    slice_min(
      abs(.data$present - .data$winter),
      n = 1,
      with_ties = FALSE,
      by = c("location_id", "winter")
    ) |>
    # remove unvisited winters with too much extrapolation
    filter(abs(.data$present - .data$winter) <= n_extrapolation) |>
    mutate(
      datafield_id = ifelse(is.na(.data$observation_id), 4L, 2L),
      observation_id = ifelse(
        is.na(.data$observation_id),
        -100000 * .data$winter - .data$location_id,
        .data$observation_id
      )
    ) |>
    select(-"present")
}
