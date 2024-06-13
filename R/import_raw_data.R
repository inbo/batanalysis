#' Import the raw observations
#' @param origin A `DBI` connection to the SQL Server database.
#' @export
#' @importFrom assertthat assert_that
#' @importFrom DBI dbGetQuery
#' @importFrom dplyr anti_join bind_rows count filter inner_join semi_join
#' transmute
#' @importFrom git2rdata write_vc
#' @importFrom rlang .data
import_raw_data <- function(origin, target) {
  assert_that(
    inherits(origin, "Microsoft SQL Server"), inherits(target, "git_repository")
  )
  individual <- read_raw_individual(origin = origin)
  section <- read_raw_section(origin = origin)
  total <- read_raw_total(origin = origin)

  visits <- bind_rows(individual$visits, section$visits)
  visits |>
    count(.data$location_id, .data$date) |>
    filter(.data$n > 1) |>
    inner_join(visits, by = c("location_id", "date")) |>
    transmute(
      .data$visit_id,
      problem =
        "duplicate visits between individual and section based protocols"
    ) |>
    bind_rows(
      total$visits |>
        semi_join(visits, by = c("location_id", "date")) |>
        transmute(
          .data$visit_id,
          problem = paste(
            "duplicate visits between total based and individual or section",
            "based protocols"
          )
        ),
      "SELECT
  v.id AS visit_id, 'visit without sample' AS problem
FROM staging_Meetnetten.projects_project AS p
INNER JOIN staging_Meetnetten.fieldwork_visit AS v ON v.project_id = p.id
INNER JOIN staging_Meetnetten.protocols_protocol AS pr ON pr.id = v.protocol_id
LEFT JOIN staging_Meetnetten.fieldwork_sample AS sa ON sa.visit_id = v.id
WHERE
  p.name = 'Vleermuizen - Wintertellingen' AND
  pr.name = 'Vleermuizen - Wintertelling (totalen per telobject)' AND
  v.validation_status <> -1 AND sa.id IS NULL" |>
        dbGetQuery(conn = origin)
    ) -> problems
  bind_rows(visits, total$visits) |>
    anti_join(problems, by = "visit_id") -> visits
  write_vc(
    visits, file = "hibernation/visits", sorting = "visit_id", stage = TRUE,
    force = TRUE, root = target
  )

  bind_rows(individual$samples, section$samples) |>
    semi_join(visits, by = "visit_id") -> samples
  write_vc(
    samples, file = "hibernation/samples", sorting = c("visit_id", "sample_id"),
    stage = TRUE, force = TRUE, root = target
  )

  bind_rows(individual$observations, section$observations) |>
    semi_join(samples, by = "sample_id") |>
    write_vc(
      file = "hibernation/samples", root = target,
      sorting = c("sample_id", "species_id"), stage = TRUE, force = TRUE
    )

  total$observations |>
    semi_join(visits, by = "visit_id") |>
    write_vc(
      file = "hibernation/totals", root = target,
      sorting = c("visit_id", "species_id"), stage = TRUE, force = TRUE
    )

  write_vc(
    problems, file = "hibernation/problems", sorting = c("visit_id", "problem"),
    optimize = FALSE, root = target
  )
}
