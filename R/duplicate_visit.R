#' Identify duplicate visits in monitoring data
#'
#' This function checks for duplicate visits in the individual-based,
#' section-based, and total-based monitoring protocols for hibernating bats.
#' It identifies visits that have been recorded multiple times for the same
#' location and date, and determines which visits to keep based on the total
#' counts of observed species.
#' The function returns a list containing cleaned visits and data frames of
#' problems found at the location and sublocation levels.
#' @inheritParams import_raw_data
#' @return A list with three elements:
#' - `clean`: A data frame with duplicates that have identical counts at both
#'   the location and sublocation levels, indicating which visit to keep.
#' - `location`: A data frame with duplicates having different totals at the
#'   location level.
#'   The data frame only lists the species with differing counts.
#' - `sublocation`: A data frame with duplicates having different totals at the
#'   sublocation level.
#'   The data frame only lists the species with differing counts.
#' @export
#' @importFrom assertthat assert_that
#' @importFrom dplyr bind_rows filter group_by inner_join mutate rename select
#' summarise ungroup
#' @importFrom tidyr complete nesting pivot_longer pivot_wider
#' @importFrom rlang .data sym
duplicate_visit <- function(origin, ignore = c("dead", "flying")) {
  assert_that(inherits(origin, "Microsoft SQL Server"))

  # read data from database
  individual <- read_raw_individual(origin = origin, ignore = ignore)
  section <- read_raw_section(origin = origin, ignore = ignore)
  total <- read_raw_total(origin = origin, ignore = ignore)
  read_raw_species(origin = origin) |>
    transmute(
      species_id = .data$id,
      species = c(
        "Mmysbra",
        "Mmys",
        "Mbra",
        "Mdau",
        "Mema",
        "Mnat",
        "Mdas",
        "Mbec",
        "Mmyo",
        "Malc",
        "Mspec",
        "Pipspec",
        "Pippip",
        "Pipnat",
        "Pauraus",
        "Paur",
        "Paus",
        "Eser",
        "Nspec",
        "Nnoc",
        "Nlei",
        "Nlas",
        "Bbar",
        "Vmur",
        "Rfer",
        "Rhip",
        "Cspec",
        .data$code
      ) |>
        unique() |>
        factor(x = .data$code)
    ) -> species
  read_raw_location(origin = origin) |>
    select("id", "name") -> locations

  # check individual protocol for duplicate visits
  individual$visits |>
    inner_join(individual$samples, by = "visit_id") |>
    inner_join(individual$observations, by = "sample_id") |>
    group_by(
      .data$visit_id,
      .data$location_id,
      .data$sublocation_id,
      .data$date,
      .data$species_id
    ) |>
    summarise(total = sum(.data$number), .groups = "drop") |>
    ungroup() |>
    complete(
      nesting(
        !!sym("location_id"),
        !!sym("sublocation_id"),
        !!sym("date"),
        !!sym("visit_id")
      ),
      !!sym("species_id"),
      fill = list(total = 0)
    ) -> individual_total
  individual$visits |>
    group_by(.data$location_id, .data$date) |>
    filter(n() > 1) |>
    ungroup() -> candidate_visits
  individual_total |>
    semi_join(candidate_visits, by = "visit_id") |>
    group_by(.data$sublocation_id, .data$date, .data$species_id) |>
    filter(suppressWarnings(diff(range(.data$total)) > 0)) |>
    ungroup() -> individual_subloc_problems
  candidate_visits |>
    anti_join(individual_subloc_problems, by = "visit_id") -> invididual_clean

  # check section protocol for duplicate visits
  section$visits |>
    inner_join(section$samples, by = "visit_id") |>
    inner_join(section$observations, by = "sample_id") |>
    group_by(
      .data$visit_id,
      .data$location_id,
      .data$sublocation_id,
      .data$date,
      .data$species_id
    ) |>
    summarise(total = sum(.data$number), .groups = "drop") |>
    ungroup() |>
    complete(
      nesting(
        !!sym("location_id"),
        !!sym("sublocation_id"),
        !!sym("date"),
        !!sym("visit_id")
      ),
      !!sym("species_id"),
      fill = list(total = 0)
    ) -> section_total
  section$visits |>
    group_by(.data$location_id, .data$date) |>
    filter(n() > 1) |>
    ungroup() -> candidate_visits
  section_total |>
    semi_join(candidate_visits, by = "visit_id") |>
    group_by(.data$sublocation_id, .data$date, .data$species_id) |>
    filter(suppressWarnings(diff(range(.data$total)) > 0)) |>
    ungroup() -> section_subloc_problems
  candidate_visits |>
    anti_join(section_subloc_problems, by = "visit_id") -> section_clean

  # both individual and section based protocol
  individual$visits |>
    rename(individual = "visit_id") |>
    inner_join(
      section$visits |>
        rename(section = "visit_id"),
      by = c("location_id", "date")
    ) -> candidate_visits
  individual_total |>
    semi_join(candidate_visits, by = c("visit_id" = "individual")) |>
    bind_rows(
      section_total |>
        semi_join(candidate_visits, by = c("visit_id" = "section"))
    ) |>
    group_by(.data$sublocation_id, .data$date, .data$species_id) |>
    filter(suppressWarnings(diff(range(.data$total)) > 0)) |>
    ungroup() -> ind_sec_subloc_problems
  candidate_visits |>
    anti_join(
      section_subloc_problems,
      by = c("section" = "visit_id")
    ) -> ind_sec_clean

  # check total based protocol for duplicate visits
  total$visits |>
    inner_join(total$observations, by = "visit_id") |>
    complete(
      nesting(!!sym("location_id"), !!sym("date")),
      !!sym("species_id"),
      fill = list(total = 0)
    ) -> total_total
  total$visits |>
    group_by(.data$location_id, .data$date) |>
    filter(n() > 1) |>
    ungroup() -> candidate_visits
  total_total |>
    semi_join(candidate_visits, by = "visit_id") |>
    group_by(.data$location_id, .data$date, .data$species_id) |>
    filter(suppressWarnings(diff(range(.data$total)) > 0)) |>
    ungroup() -> total_loc_problems
  candidate_visits |>
    anti_join(total_loc_problems, by = "visit_id") -> total_clean

  # individual and total based protocol combined
  individual$visits |>
    rename(individual = "visit_id") |>
    inner_join(
      total$visits |>
        rename(total = "visit_id"),
      by = c("location_id", "date")
    ) -> candidate_visits
  individual_total |>
    semi_join(candidate_visits, by = c("visit_id" = "individual")) |>
    group_by(.data$visit_id, .data$location_id, .data$date, .data$species_id) |>
    summarise(total = sum(.data$total), .groups = "drop") |>
    complete(
      nesting(!!sym("location_id"), !!sym("date"), !!sym("visit_id")),
      !!sym("species_id"),
      fill = list(total = 0)
    ) |>
    bind_rows(
      total_total |>
        semi_join(candidate_visits, by = c("visit_id" = "total"))
    ) |>
    group_by(.data$location_id, .data$date, .data$species_id) |>
    filter(suppressWarnings(diff(range(.data$total)) > 0)) |>
    ungroup() -> ind_tot_loc_problems
  candidate_visits |>
    anti_join(
      ind_tot_loc_problems,
      by = c("individual" = "visit_id")
    ) -> ind_tot_clean

  # section and total based protocol combined
  section$visits |>
    rename(section = "visit_id") |>
    inner_join(
      total$visits |>
        rename(total = "visit_id"),
      by = c("location_id", "date")
    ) -> candidate_visits
  section_total |>
    semi_join(candidate_visits, by = c("visit_id" = "section")) |>
    group_by(.data$visit_id, .data$location_id, .data$date, .data$species_id) |>
    summarise(total = sum(.data$total), .groups = "drop") |>
    complete(
      nesting(!!sym("location_id"), !!sym("date"), !!sym("visit_id")),
      !!sym("species_id"),
      fill = list(total = 0)
    ) |>
    bind_rows(
      total_total |>
        semi_join(candidate_visits, by = c("visit_id" = "total"))
    ) |>
    group_by(.data$location_id, .data$date, .data$species_id) |>
    filter(suppressWarnings(diff(range(.data$total)) > 0)) |>
    ungroup() -> sec_tot_loc_problems
  candidate_visits |>
    anti_join(
      sec_tot_loc_problems,
      by = c("section" = "visit_id")
    ) -> sec_tot_clean

  invididual_clean |>
    bind_rows(section_clean, total_clean) |>
    group_by(.data$location_id, .data$date) |>
    mutate(keep = .data$visit_id == min(.data$visit_id)) |>
    ungroup() |>
    bind_rows(
      ind_sec_clean |>
        pivot_longer(
          cols = c("individual", "section"),
          names_to = "keep",
          values_to = "visit_id"
        ) |>
        mutate(keep = .data$keep == "individual"),
      ind_tot_clean |>
        pivot_longer(
          cols = c("individual", "total"),
          names_to = "keep",
          values_to = "visit_id"
        ) |>
        mutate(keep = .data$keep == "individual"),
      sec_tot_clean |>
        pivot_longer(
          cols = c("total", "section"),
          names_to = "keep",
          values_to = "visit_id"
        ) |>
        mutate(keep = .data$keep == "section")
    ) |>
    arrange(
      .data$location_id,
      .data$date,
      .data$keep,
      .data$visit_id
    ) -> clean_duplicates

  individual_subloc_problems |>
    bind_rows(section_subloc_problems, ind_sec_subloc_problems) |>
    inner_join(
      locations |>
        select(location_id = "id", location = "name"),
      by = "location_id"
    ) |>
    inner_join(
      locations |>
        select(sublocation_id = "id", sublocation = "name"),
      by = "sublocation_id"
    ) |>
    inner_join(species, by = "species_id") |>
    arrange(
      .data$location_id,
      .data$sublocation_id,
      .data$date,
      .data$visit_id
    ) |>
    select(
      "visit_id",
      "location",
      "sublocation",
      "date",
      "species",
      "total"
    ) |>
    pivot_wider(
      names_from = "species",
      names_sort = TRUE,
      values_from = "total"
    ) -> sublocation_problems
  ind_tot_loc_problems |>
    bind_rows(sec_tot_loc_problems, total_loc_problems) |>
    inner_join(
      locations |>
        select(location_id = "id", location = "name"),
      by = "location_id"
    ) |>
    inner_join(species, by = "species_id") |>
    arrange(.data$location_id, .data$date, .data$visit_id) |>
    select(
      "visit_id",
      "location",
      "date",
      "species",
      "total"
    ) |>
    pivot_wider(
      names_from = "species",
      names_sort = TRUE,
      values_from = "total"
    ) -> location_problems
  return(list(
    clean = clean_duplicates,
    location = location_problems,
    sublocation = sublocation_problems
  ))
}
