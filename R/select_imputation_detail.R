#' Select detailed data for imputation
#' @inheritParams prepare_analysis_model_species
#' @param locations A data frame with location information.
#' Must contain the columns `location_id`, `detailed` and `type`.
#' @param this_species A data frame with one ore more relevant species ID.
#' @param start The oldest date to use in the analysis.
#' @export
#' @importFrom dplyr anti_join bind_rows count distinct filter group_by
#' inner_join left_join mutate n select semi_join slice_min summarise ungroup
#' @importFrom git2rdata verify_vc
#' @importFrom lubridate round_date year
#' @importFrom tidyr complete replace_na nesting
#' @importFrom rlang .data sym
select_imputation_detail <- function(
  locations,
  raw_data,
  this_species,
  start,
  n_winter = 4,
  n_present = 3,
  max_delta = 10,
  n_extrapolation = 5
) {
  # read samples from file
  file.path("data", "hibernation", "samples") |>
    verify_vc(
      root = raw_data,
      variables = c("visit_id", "sample_id", "sublocation_id")
    ) -> samples
  # read observations from file
  file.path("data", "hibernation", "observations") |>
    verify_vc(
      root = raw_data,
      variables = c("sample_id", "species_id", "number")
    ) |>
    # select only observations for this species
    semi_join(this_species, by = c("species_id" = "id")) |>
    # aggregate in case of multiple species
    group_by(.data$sample_id) |>
    summarise(number = sum(.data$number)) -> observations
  # read totals from file
  file.path("data", "hibernation", "totals") |>
    verify_vc(
      root = raw_data,
      variables = c("visit_id", "species_id", "total")
    ) |>
    # select only totals for this species
    semi_join(this_species, by = c("species_id" = "id")) |>
    # aggregate in case of multiple species
    group_by(.data$visit_id) |>
    summarise(number = sum(.data$total)) -> totals
  # read visits from file
  file.path("data", "hibernation", "visits") |>
    verify_vc(root = raw_data, variables = c("visit_id", "location_id")) |>
    # keep visits from start date only
    filter(.data$date >= start) |>
    # keep visits from detailed locations only
    semi_join(
      locations |>
        filter(.data$detailed),
      by = "location_id"
    ) |>
    mutate(
      winter = round_date(.data$date, "year") |>
        year(),
      delta = paste0(.data$winter, "-1-15 12:0:0") |>
        as.POSIXct(),
      delta = difftime(.data$date, .data$delta, units = "days") |>
        as.numeric()
    ) -> visits
  visits |>
    # visits with detailed counts
    semi_join(samples, by = "visit_id") -> detailed_visits
  visits |>
    # visits with only totals
    anti_join(
      detailed_visits,
      by = c("location_id", "winter")
    ) -> non_detailed_visits
  detailed_visits |>
    # ignore multiple visits per winter and location
    distinct(.data$location_id, .data$winter) |>
    # count the number of winters with detailed counts per location
    count(.data$location_id, name = "detail") |>
    # keep only locations with at least `n_winter` detailed counts
    filter(.data$detail >= n_winter) |>
    left_join(
      non_detailed_visits |>
        distinct(.data$location_id, .data$winter) |>
        # count the number of winters with only totals per location
        count(.data$location_id, name = "total"),
      by = "location_id"
    ) |>
    # keep locations with more detailed than total counts
    filter(.data$detail > .data$total | is.na(.data$total)) -> detail_locations
  samples |>
    # add visit information to samples
    inner_join(detailed_visits, by = "visit_id") |>
    # keep per sublocation and winter the closest visit to Jan 15
    slice_min(
      abs(.data$delta),
      with_ties = FALSE,
      by = c("location_id", "sublocation_id", "winter")
    ) |>
    # keep only multiple visits per winter and location if they are close enough
    group_by(.data$location_id, .data$winter) |>
    filter(diff(range(.data$delta)) <= max_delta) |>
    ungroup() |>
    # add the observations to the samples
    left_join(observations, by = "sample_id") |>
    # missing observations indicate no observation of this species
    mutate(number = replace_na(.data$number, 0)) |>
    select(
      -"date",
      -"delta",
      -"visit_id",
      observation_id = "sample_id"
    ) -> detailed_candidate
  detailed_candidate |>
    group_by(.data$sublocation_id) |>
    # keep only sublocations with at least `n_winter` visits
    # and at least `n_present` visits with observed bats
    filter(n() >= n_winter, sum(.data$number > 0) >= n_present) |>
    ungroup() -> detailed_samples
  detailed_candidate |>
    anti_join(detailed_samples, by = "sublocation_id") |>
    filter(.data$number > 0) |>
    inner_join(
      locations |>
        select("location_id", "type"),
      by = "location_id"
    ) |>
    mutate(
      fortress = as.integer(.data$type == "fortress"),
      marl_quarry = as.integer(.data$type == "marl quarry"),
      other_large = as.integer(.data$type == "other large"),
      small = 0L,
      datafield_id = 1L
    ) -> rare_sublocations

  detailed_samples |>
    complete(
      nesting(!!sym("location_id"), !!sym("sublocation_id")),
      winter = min(.data$winter):(max(.data$winter) + 1)
    ) |>
    inner_join(
      detailed_samples |>
        filter(.data$number > 1) |>
        select("sublocation_id", present = "winter"),
      by = "sublocation_id",
      relationship = "many-to-many"
    ) |>
    slice_min(
      abs(.data$present - .data$winter),
      n = 1,
      with_ties = FALSE,
      by = c("sublocation_id", "winter")
    ) |>
    filter(abs(.data$present - .data$winter) <= n_extrapolation) |>
    mutate(
      datafield_id = ifelse(is.na(.data$observation_id), 3L, 1L),
      observation_id = ifelse(
        is.na(.data$observation_id),
        -1000000 * .data$winter - .data$sublocation_id,
        .data$observation_id
      )
    ) |>
    select(-"present") -> sublocation_data
  visits |>
    anti_join(detail_locations, by = "location_id") |>
    slice_min(
      abs(.data$delta),
      with_ties = FALSE,
      by = c("location_id", "winter")
    ) |>
    select(-"date", -"delta") |>
    inner_join(samples, by = "visit_id") |>
    left_join(observations, by = "sample_id") |>
    mutate(number = replace_na(.data$number, 0)) |>
    group_by(.data$visit_id, .data$location_id, .data$winter) |>
    summarise(number = sum(.data$number), .groups = "drop") -> detail_to_total
  visits |>
    anti_join(detail_locations, by = "location_id") |>
    anti_join(detail_to_total, by = c("location_id", "winter")) |>
    slice_min(
      abs(.data$delta),
      with_ties = FALSE,
      by = c("location_id", "winter")
    ) |>
    select(-"date", -"delta") |>
    left_join(totals, by = "visit_id") |>
    mutate(number = replace_na(.data$number, 0)) |>
    bind_rows(detail_to_total) |>
    group_by(.data$location_id) |>
    filter(n() >= n_winter, sum(.data$number > 0) >= n_present) |>
    ungroup() -> only_totals
  if (nrow(only_totals) == 0) {
    location_data <- data.frame(
      location_id = integer(0),
      winter = integer(0),
      observation_id = integer(0),
      number = integer(0),
      datafield_id = integer(0)
    )
  } else {
    only_totals |>
      complete(
        .data$location_id,
        winter = min(.data$winter):(max(.data$winter) + 1)
      ) |>
      inner_join(
        only_totals |>
          filter(.data$number > 1) |>
          select("location_id", present = "winter"),
        by = "location_id",
        relationship = "many-to-many"
      ) |>
      slice_min(
        abs(.data$present - .data$winter),
        n = 1,
        with_ties = FALSE,
        by = c("location_id", "winter")
      ) |>
      filter(abs(.data$present - .data$winter) <= n_extrapolation) |>
      mutate(
        datafield_id = ifelse(is.na(.data$visit_id), 4L, 2L),
        observation_id = ifelse(
          is.na(.data$visit_id),
          -100000 * .data$winter - .data$location_id,
          .data$visit_id
        )
      ) |>
      select(-"present", -"visit_id") -> location_data
  }
  return(list(
    sublocation = sublocation_data,
    location = location_data,
    rare_sublocations = rare_sublocations
  ))
}
