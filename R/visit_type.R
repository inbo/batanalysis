#' Visit Type Classification
#' @export
#' @inheritParams prepare_analysis_model_species
#' @importFrom dplyr bind_rows filter group_by left_join inner_join mutate n
#' select semi_join transmute ungroup
#' @importFrom git2rdata verify_vc
#' @importFrom lubridate round_date year
#' @importFrom rlang .data
visit_type <- function(
  raw_data,
  start_winter = as.integer(format(Sys.Date(), "%Y")) - 23,
  n_winter = 4,
  max_delta = 10
) {
  start <- as.Date(paste0(start_winter, "-1-1"))
  file.path("data", "hibernation", "visits") |>
    verify_vc(
      root = raw_data,
      variables = c("visit_id", "location_id", "date")
    ) |>
    mutate(
      winter = round_date(.data$date, "year") |>
        year(),
      delta = paste0(.data$winter, "-1-15 9:0:0") |>
        as.POSIXct(),
      delta = difftime(.data$date, .data$delta, units = "days") |>
        as.numeric()
    ) -> visits
  file.path("data", "hibernation", "locations") |>
    verify_vc(
      root = raw_data,
      variables = c("id", "parent_id")
    ) -> all_locations
  all_locations |>
    filter(.data$parent_id > 0) -> sublocations
  all_locations |>
    filter(.data$parent_id < 0) |>
    anti_join(
      sublocations,
      by = c("id" = "parent_id")
    ) -> non_detailed_locations
  all_locations |>
    filter(.data$parent_id < 0) |>
    semi_join(sublocations, by = c("id" = "parent_id")) -> detailed_locations
  visits |>
    semi_join(
      non_detailed_locations,
      by = c("location_id" = "id")
    ) |>
    group_by(.data$location_id, .data$winter) |>
    transmute(
      .data$visit_id,
      .data$location_id,
      .data$date,
      type = ifelse(
        abs(.data$delta) == min(abs(.data$delta)),
        ifelse(.data$date < start, "old", "potential"),
        "extra"
      )
    ) |>
    group_by(.data$location_id, .data$type) |>
    mutate(
      type = ifelse(
        .data$type == "potential",
        ifelse(n() >= n_winter, "potential", "rare"),
        .data$type
      ),
      score = 1 /
        ifelse(
          .data$type == "potential",
          max(visits$winter) - .data$winter + 1,
          NA_integer_
        ) /
        sum(1 / seq_len(24))
    ) |>
    ungroup() -> non_detailed_visits
  visits |>
    anti_join(
      non_detailed_locations,
      by = c("location_id" = "id")
    ) -> remainder
  file.path("data", "hibernation", "samples") |>
    verify_vc(
      root = raw_data,
      variables = c("visit_id")
    ) |>
    semi_join(x = remainder, by = "visit_id") |>
    group_by(.data$location_id, .data$winter) |>
    transmute(
      .data$visit_id,
      .data$location_id,
      .data$date,
      type = ifelse(
        abs(.data$delta) == min(abs(.data$delta)),
        ifelse(.data$date < start, "old detail", "detail"),
        "extra detail"
      )
    ) |>
    ungroup() -> samples
  samples |>
    filter(.data$type == "extra detail") |>
    inner_join(
      samples |>
        filter(.data$type == "detail") |>
        select("location_id", "winter", best = "date"),
      by = c("location_id", "winter")
    ) |>
    mutate(
      type = ifelse(
        abs(.data$best - .data$date) <= max_delta,
        "alternative detail",
        "extra detail"
      )
    ) |>
    select(-"best") |>

    bind_rows(
      samples |>
        filter(.data$type != "extra detail"),
      remainder |>
        anti_join(samples, by = "visit_id") |>
        left_join(
          samples |>
            distinct(.data$location_id, .data$winter, type = "total_detail"),
          by = c("location_id", "winter")
        ) |>
        group_by(.data$location_id, .data$winter) |>
        mutate(
          type = ifelse(
            is.na(.data$type),
            ifelse(
              abs(.data$delta) == min(abs(.data$delta)),
              ifelse(.data$date < start, "old total", "total"),
              "extra total"
            ),
            .data$type
          )
        ) |>
        select(-"delta")
    ) -> detailed
  detailed |>
    filter(.data$type %in% c("detail", "total")) |>
    count(.data$location_id, .data$type) |>
    pivot_wider(names_from = "type", values_from = "n", values_fill = 0) |>
    transmute(
      .data$location_id,
      level = ifelse(
        .data$detail + .data$total >= n_winter,
        ifelse(.data$detail > .data$total, "detail", "mixed"),
        "rare"
      )
    ) |>
    left_join(x = detailed, by = "location_id") |>
    group_by(.data$location_id) |>
    mutate(
      score = ifelse(.data$type == "detail", 1, 0.5) /
        ifelse(
          .data$type %in% c("total", "detail") & .data$level != "rare",
          max(visits$winter) - .data$winter + 1,
          NA_integer_
        ) /
        sum(1 / seq_len(24))
    ) |>
    ungroup() |>
    bind_rows(non_detailed_visits)
}
