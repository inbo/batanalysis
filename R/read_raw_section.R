#' Import the raw observations according to the section based protocol
#' @inheritParams import_raw_data
#' @export
#' @importFrom assertthat assert_that
#' @importFrom dplyr anti_join bind_rows count distinct filter inner_join select
#' transmute
#' @importFrom DBI dbGetQuery
#' @importFrom rlang .data
read_raw_section <- function(origin) {
  assert_that(inherits(origin, "Microsoft SQL Server"))
  "SELECT
  v.location_id, sa.location_id AS sublocation_id, v.id AS visit_id,
  sa.id AS sample_id, v.start_date AS date, sa.not_counted, o.species_id,
  o.number_min AS number, a.name AS activity
FROM staging_Meetnetten.projects_project AS p
INNER JOIN staging_Meetnetten.fieldwork_visit AS v ON v.project_id = p.id
INNER JOIN staging_Meetnetten.protocols_protocol AS pr ON pr.id = v.protocol_id
INNER JOIN staging_Meetnetten.fieldwork_sample AS sa ON sa.visit_id = v.id
LEFT JOIN staging_Meetnetten.fieldwork_observation AS o ON o.sample_id = sa.id
LEFT JOIN staging_Meetnetten.species_activity AS a ON a.id = o.activity_id
LEFT JOIN staging_Meetnetten.projects_projectspecies AS psp ON
  psp.species_id = o.species_id AND psp.project_id = p.id
WHERE
  p.name = 'Vleermuizen - Wintertellingen' AND
  pr.name = 'Vleermuizen - Wintertelling (totalen per kamer/sectie)' AND
  v.validation_status <> -1" |>
    dbGetQuery(conn = origin) |>
    filter(!.data$activity %in% c("awake", "dead", "flying")) -> raw_data
  raw_data |>
    distinct(.data$visit_id, .data$location_id, .data$date) -> raw_visit
  raw_data |>
    filter(!is.na(.data$species_id)) |>
    count(.data$sample_id, .data$species_id) |>
    filter(.data$n > 1) |>
    inner_join(raw_data, by = "sample_id") |>
    distinct(.data$visit_id) |>
    transmute(
      .data$visit_id,
      problem = "multiple samples per species in section based protocol"
    ) |>
    bind_rows(
      raw_visit |>
        count(.data$location_id, .data$date) |>
        filter(.data$n > 1) |>
        inner_join(raw_visit, by = c("location_id", "date")) |>
        transmute(
          .data$visit_id, problem = "duplicate visit in section based protocol"
        )
    ) -> problems
  raw_data |>
    anti_join(problems, by = "visit_id") -> raw_data
  raw_data |>
    distinct(.data$visit_id, .data$location_id, .data$date) -> visits
  raw_data |>
    filter(!.data$not_counted) |>
    distinct(.data$sample_id, .data$visit_id, .data$sublocation_id) -> samples
  raw_data |>
    filter(!is.na(.data$species_id)) |>
    select("sample_id", "species_id", "number") -> observations
  return(
    list(
      visits = visits, samples = samples, observations = observations,
      problems = problems
    )
  )
}
