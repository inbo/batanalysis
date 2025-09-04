#' Import the raw observations
#' @param origin A `DBI` connection to the SQL Server database.
#' @param target A `git_repository` object to store the imported data.
#' @param ignore A character vector with the activities to ignore.
#' Defaults to `c("dead", "flying")`.
#' @inheritParams git2rdata::write_vc
#' @export
#' @importFrom assertthat assert_that
#' @importFrom DBI dbGetQuery
#' @importFrom dplyr anti_join bind_rows count distinct filter inner_join
#' semi_join transmute
#' @importFrom git2rdata write_vc
#' @importFrom rlang .data
import_raw_data <- function(
  origin,
  target,
  ignore = c("dead", "flying"),
  strict = TRUE
) {
  assert_that(
    inherits(origin, "Microsoft SQL Server"),
    inherits(target, "git_repository")
  )
  individual <- read_raw_individual(origin = origin, ignore = ignore)
  section <- read_raw_section(origin = origin, ignore = ignore)
  total <- read_raw_total(origin = origin, ignore = ignore)
  species <- read_raw_species(origin = origin)
  locations <- read_raw_location(origin = origin)

  visits <- bind_rows(individual$visits, section$visits)
  visits |>
    count(.data$location_id, .data$date) |>
    filter(.data$n > 1) |>
    inner_join(visits, by = c("location_id", "date")) |>
    transmute(
      .data$visit_id,
      .data$location_id,
      problem = paste(
        "duplicate visits between individual and section based",
        "protocols"
      )
    ) |>
    bind_rows(
      total$problem,
      individual$problem,
      section$problem,
      total$visits |>
        semi_join(visits, by = c("location_id", "date")) |>
        transmute(
          .data$visit_id,
          .data$location_id,
          problem = paste(
            "duplicate visits between total based and individual or section",
            "based protocols"
          )
        ),
      "SELECT
  v.id AS visit_id, v.location_id, 'visit without sample' AS problem
FROM staging_Meetnetten.projects_project AS p
INNER JOIN staging_Meetnetten.fieldwork_visit AS v ON v.project_id = p.id
INNER JOIN staging_Meetnetten.protocols_protocol AS pr ON pr.id = v.protocol_id
LEFT JOIN staging_Meetnetten.fieldwork_sample AS sa ON sa.visit_id = v.id
WHERE
  p.name = 'Vleermuizen - Wintertellingen' AND
  pr.name = 'Vleermuizen - Wintertelling (totalen per telobject)' AND
  v.validation_status <> -1 AND sa.id IS NULL" |>
        dbGetQuery(conn = origin)
    ) |>
    mutate(problem = factor(.data$problem)) -> problems

  visits <- bind_rows(visits, total$visits)
  file.path("data", "hibernation", "visits") |>
    write_vc(
      x = visits,
      sorting = "visit_id",
      stage = TRUE,
      force = TRUE,
      root = target,
      strict = strict
    )
  file.path("data", "hibernation", "visits") |>
    update_metadata(
      stage = TRUE,
      force = TRUE,
      root = target,
      name = "visits",
      title = "Available visits of the hibernating bat monitoring",
      field_description = c(
        visit_id = "Unique identifier of the visit",
        date = "Date of the visit",
        location_id = "Unique identifier of the location"
      )
    )

  samples <- bind_rows(individual$samples, section$samples)
  file.path("data", "hibernation", "samples") |>
    write_vc(
      x = samples,
      sorting = c("visit_id", "sample_id"),
      stage = TRUE,
      force = TRUE,
      root = target,
      strict = strict
    )
  file.path("data", "hibernation", "samples") |>
    update_metadata(
      stage = TRUE,
      force = TRUE,
      root = target,
      name = "samples",
      title = paste(
        "Available visits at the sublocation level of the hibernating bat",
        "monitoring"
      ),
      field_description = c(
        visit_id = "Unique identifier of the visit",
        sample_id = "Unique identifier of the sample",
        sublocation_id = "Unique identifier of the sublocation"
      )
    )

  observations <- bind_rows(individual$observations, section$observations)
  file.path("data", "hibernation", "observations") |>
    write_vc(
      x = observations,
      root = target,
      sorting = c("sample_id", "species_id"),
      stage = TRUE,
      force = TRUE,
      strict = strict
    )
  file.path("data", "hibernation", "observations") |>
    update_metadata(
      stage = TRUE,
      force = TRUE,
      root = target,
      name = "observations",
      title = paste(
        "Number of observed bats by species at the sublocation level of the",
        "hibernating bat monitoring"
      ),
      field_description = c(
        species_id = "Unique identifier of the species",
        sample_id = "Unique identifier of the sample",
        number = "Number of observed bats"
      )
    )

  file.path("data", "hibernation", "totals") |>
    write_vc(
      x = total$observations,
      root = target,
      sorting = c("visit_id", "species_id"),
      stage = TRUE,
      force = TRUE,
      strict = strict
    )
  file.path("data", "hibernation", "totals") |>
    update_metadata(
      stage = TRUE,
      force = TRUE,
      root = target,
      name = "totals",
      title = "Total number of observed bats by species at the location level.
Only given when no observations at the sublocation level are available.",
      field_description = c(
        visit_id = "Unique identifier of the visit",
        species_id = "Unique identifier of the species",
        total = "Total number of observed bats"
      )
    )

  file.path("data", "hibernation", "species") |>
    write_vc(
      x = species,
      root = target,
      sorting = "id",
      stage = TRUE,
      force = TRUE,
      strict = strict
    )
  file.path("data", "hibernation", "species") |>
    update_metadata(
      stage = TRUE,
      force = TRUE,
      root = target,
      name = "species",
      title = "Species observed during the hibernating bat monitoring",
      field_description = c(
        id = "Unique identifier of the species",
        name = "Dutch vernacular name",
        scientific_name = "Scientific name of the species",
        code = "Code of the species",
        parent = "Parent species"
      )
    )

  file.path("data", "hibernation", "locations") |>
    write_vc(
      x = locations,
      root = target,
      sorting = "id",
      stage = TRUE,
      force = TRUE,
      digits = 6,
      strict = strict
    )
  file.path("data", "hibernation", "locations") |>
    update_metadata(
      stage = TRUE,
      force = TRUE,
      root = target,
      name = "locations",
      title = paste(
        "Locations and sublocations observed during the hibernating",
        "bat monitoring"
      ),
      field_description = c(
        id = "Unique identifier of the location or sublocation",
        name = "Name of the location or sublocation",
        parent_id = "Parent location",
        code = "Code of the location or sublocation",
        longitude = "Longitude of the location or sublocation",
        latitude = "Latitude of the location or sublocation"
      )
    )

  file.path("data", "hibernation", "problems") |>
    write_vc(
      x = problems,
      sorting = c("visit_id", "problem"),
      optimize = FALSE,
      root = target,
      stage = TRUE,
      force = TRUE,
      strict = strict
    )
  file.path("data", "hibernation", "problems") |>
    update_metadata(
      stage = TRUE,
      force = TRUE,
      root = target,
      name = "problems",
      title = paste(
        "Issues found during the import of the hibernating bat",
        "monitoring data"
      ),
      field_description = c(
        visit_id = "Unique identifier of the visit",
        location_id = "Unique identifier of the location",
        problem = "Description of the issue"
      )
    )
}
