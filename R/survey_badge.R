#' Survey badge
#'
#' This function calculates a survey badge based on the number of visits to bat
#' monitoring locations.
#' The badge looks at the usable visits over the last 24 winters.
#' A visit is usable if it has either a total count or section counts.
#' If a location has sublocations, only visits with section counts are counted.
#' The score is calculated as a weighted sum of the visits, with more recent
#' visits having a higher weight.
#' The weights are calculated as `1 / n`, where `n` is the number of winters
#' ago the visit took place.
#' We rescale the weights so that the maximum score is 1.
#' The badge is then assigned based on the score:
#' - Gold: a score equal to or higher than the score of having a visit in
#'  the 75% most recent winters.
#'  - Silver: a score equal to or higher than the score of having a visit in
#'  the 50% most recent winters.
#'  - Bronze: a score equal to or higher than the score of having a visit in
#'  the 25% most recent winters.
#'  - None: a score lower than the score of having a visit in the 25% most
#'  recent winters.
#' @inheritParams import_raw_data
#' @param duration The number of recent winters to consider for the badge.
#' @return A data frame with the location ID, badge type, and score.
#' @export
#' @importFrom dplyr anti_join bind_rows count distinct filter inner_join
#' left_join mutate select semi_join
#' @importFrom git2rdata verify_vc
#' @importFrom lubridate round_date year
#' @importFrom rlang .data
#' @importFrom tidyr replace_na
survey_badge <- function(target, duration = 24) {
  verify_vc(
    "data/hibernation/visits",
    root = target,
    variables = c("visit_id", "location_id", "date")
  ) |>
    anti_join(
      verify_vc(
        "data/hibernation/problems",
        root = target,
        variables = "location_id"
      ),
      by = "visit_id"
    ) |>
    mutate(
      winter = round_date(.data$date, "year") |>
        year(),
      winter = max(.data$winter) - .data$winter + 1
    ) |>
    filter(.data$winter <= duration) -> relevant_visit
  verify_vc(
    "data/hibernation/locations",
    root = target,
    variables = c("id", "parent_id")
  ) -> locations
  verify_vc(
    "data/hibernation/totals",
    root = target,
    variables = "visit_id"
  ) -> totals
  verify_vc(
    "data/hibernation/samples",
    root = target,
    variables = "visit_id"
  ) -> samples
  relevant_visit |>
    semi_join(totals, by = "visit_id") |>
    mutate(type = "total") |>
    bind_rows(
      relevant_visit |>
        semi_join(samples, by = "visit_id") |>
        mutate(type = "section"),
      relevant_visit |>
        anti_join(totals, by = "visit_id") |>
        anti_join(samples, by = "visit_id") |>
        mutate(type = "zero count")
    ) |>
    inner_join(
      locations |>
        filter(.data$parent_id == -1) |>
        select("location_id" = "id", "name") |>
        left_join(
          locations |>
            filter(.data$parent_id != -1) |>
            count(location_id = .data$parent_id, name = "sublocations") |>
            mutate(sublocations = .data$sublocations > 1),
          by = "location_id"
        ) |>
        mutate(sublocations = replace_na(.data$sublocations, FALSE)),
      by = "location_id"
    ) |>
    filter(.data$type == "section" | !.data$sublocations) |>
    distinct(.data$location_id, .data$winter) -> relevant
  max_weight <- sum(1 / seq_len(duration))
  relevant |>
    group_by(.data$location_id) |>
    summarise(score = sum(1 / .data$winter) / max_weight) |>
    mutate(
      badge = ifelse(
        .data$score < sum(1 / seq_len(duration / 2)) / max_weight,
        ifelse(
          .data$score < sum(1 / seq_len(duration / 4)) / max_weight,
          "none",
          "bronze"
        ),
        ifelse(
          .data$score < sum(1 / seq_len(duration * 3 / 4)) / max_weight,
          "silver",
          "gold"
        )
      ) |>
        factor(levels = c("none", "bronze", "silver", "gold"))
    )
}
