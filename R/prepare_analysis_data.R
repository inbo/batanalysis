#' Convert raw data to analysis data for a single species
#' @param raw_data A gity repository with the raw data.
#' @param analysis_data A gity repository to store the analysis data.
#' @inheritParams select_imputation_section
#' @export
#' @importFrom assertthat assert_that is.count is.number noNA
#' @importFrom dplyr anti_join bind_cols bind_rows distinct filter group_by
#' inner_join mutate semi_join slice_min summarise ungroup %>%
#' @importFrom git2rdata read_vc update_metadata write_vc
#' @importFrom rlang .data
#' @importFrom sf st_as_sf st_coordinates st_drop_geometry st_transform
#' @importFrom tidyr complete nesting
prepare_analysis_data_species <- function(
  raw_data, analysis_data, species, start, n_winter = 4, n_present = 3
) {
  assert_that(
    is.number(start), is.count(n_winter), is.count(n_present), noNA(start),
    noNA(n_winter), noNA(n_present)
  )
  sprintf("%i-10-01", start - 1) |>
    as.Date() -> start
  get_child_species(target = raw_data, species = species) |>
    inner_join(
      read_vc("hibernation/species", root = raw_data) |>
        select("id", "scientific_name"),
      by = "id"
    ) -> this_species
  stopifnot("No matching species found" = nrow(this_species) > 0)
  message("Preparing the data for ", this_species$scientific_name[1])
  sections <- select_imputation_section(
    target = raw_data, species = species, start = as.Date(start), n_winter = 1
  )
  select_imputation_total(
    target = raw_data, species = species, start = as.Date(start)
  ) |>
    filter(!is.na(.data$total)) |>
    select(-"minimum", number = "total") |>
    mutate(sublocation_id = .data$location_id) |>
    bind_rows(sections) |>
    select(-"cluster") |>
    complete(.data$winter, nesting(location_id, sublocation_id)) |>
    inner_join(
      read_vc("hibernation/locations", root = raw_data) |>
        select(location_id = "id", "longitude", "latitude") |>
        st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
        st_transform(crs = 31370) %>%
        bind_cols(st_coordinates(.) / 1000) |>
        st_drop_geometry(),
      by = "location_id"
    ) |>
    group_by(.data$sublocation_id) |>
    mutate(
      relevant = sum(.data$number > 0, na.rm = TRUE),
      across(c("winter", "number"), ~as.integer(.x))
    ) |>
    ungroup() -> dataset

  # remove locations with too few winters
  dataset |>
    filter(!is.na(.data$number)) |>
    distinct(.data$location_id, .data$winter) |>
    count(.data$location_id) |>
    filter(.data$n >= n_winter) -> observed_winters
  observed_winters |>
    anti_join(x = dataset, by = "location_id") |>
    filter(!is.na(.data$number)) |>
    group_by(.data$location_id, .data$winter) |>
    summarise(
      total = sum(.data$number), .groups = "drop"
    ) |>
    write_vc(
      file.path("hibernation", tolower(species), "too_short"), optimize = FALSE,
      root = analysis_data, sorting = c("location_id", "winter"), stage = TRUE
    )
  file.path("hibernation", tolower(species), "too_short") |>
    update_metadata(
      root = analysis_data, name = sprintf("too_short_%s", tolower(species)),
      title = sprintf(
        "Locations with too few surveyed winters for %s",
        this_species$scientific_name[1]
      ),
      description = sprintf("
Available total number of observed invidiuals of %s per winter at locations
where the number of surveyed winters is less than %i. The dataset starts at %s.
Note that locations with multiple sublocations require raw data at the
sublocation level. We can only use surveys containing only total numbers when
the location consist of a single sublocation.",
        this_species$scientific_name[1], n_winter, start
      ),
      field_description = c(
        location_id = "The identifier of the location.",
        winter = "The year of January 1st during the winter.",
        total = "The total number of observed individuals."
      )
    )
  dataset |>
    semi_join(observed_winters, by = "location_id") -> dataset

  # remove duplicate observations
  dataset |>
    count(.data$sublocation_id, .data$winter) |>
    filter(.data$n > 1) |>
    select(-"n") -> duplicates
  write_vc(
    duplicates, file.path("hibernation", tolower(species), "duplicates"),
    optimize = FALSE, root = analysis_data,
    sorting = c( "sublocation_id", "winter"), stage = TRUE
  )
  file.path("hibernation", tolower(species), "duplicates") |>
    update_metadata(
      root = analysis_data,
      name = sprintf("duplicates_%s", tolower(species)),
      title = sprintf(
  "Combinations of sublocations and winters with duplicate observations of %s",
        this_species$scientific_name[1]
      ),
      description = sprintf("
Sublocation with duplicate observations for %s in a given per winter. The
dataset starts at %s.",
        this_species$scientific_name[1], start
      ),
      field_description = c(
        sublocation_id = "The identifier of the sublocation.",
        winter = "The year of January 1st during the winter."
      )
    )
  dataset |>
    semi_join(duplicates, by = c("sublocation_id", "winter")) |>
    group_by(.data$location_id, .data$sublocation_id, .data$winter) |>
    slice_min(.data$sample_id, n = 1) |>
    ungroup() |>
    mutate(number = NA_integer_) |>
    bind_rows(
      dataset |>
        anti_join(duplicates, by = c("sublocation_id", "winter"))
    ) -> dataset

  # remove rare sublocations
  dataset |>
    filter(
      0 < .data$relevant, .data$relevant < n_present, !is.na(.data$number)
    ) |>
    select("location_id", "sublocation_id", "winter", "sample_id", "number") |>
    write_vc(
      file.path("hibernation", tolower(species), "rare_sublocation"),
      optimize = FALSE, root = analysis_data, stage = TRUE,
      sorting = c("location_id", "sublocation_id", "winter")
    )
  file.path("hibernation", tolower(species), "rare_sublocation") |>
    update_metadata(
      root = analysis_data,
      name = sprintf("rare_sublocation_%s", tolower(species)),
      title = sprintf(
    "Time series with observations of %s at sublocations with too few winters",
        this_species$scientific_name[1]
      ),
      description = sprintf("
Time series of the number of observed invidiuals for %s per winter at
sublocations where the number of observed winters is less than %i. The dataset
starts at %s.",
        this_species$scientific_name[1], n_present, start
      ),
      field_description = c(
        location_id = "The identifier of the location.",
        sublocation_id = "The identifier of the sublocation.",
        winter = "The year of January 1st during the winter.",
        sample_id = "The identifier of the sample.",
        number = "The number of observed individuals."
      )
    )
  dataset |>
    filter(.data$relevant >= n_present) -> dataset

  # store relevant data
  dataset |>
    select("location_id", "sublocation_id", "winter", "sample_id", "number") |>
    mutate(
      sample_id = ifelse(
        !is.na(.data$sample_id), .data$sample_id,
        -.data$sublocation_id * 10000L - .data$winter
      )
    ) |>
    write_vc(
      file.path("hibernation", tolower(species), "analysis_data"),
      optimize = FALSE, root = analysis_data, stage = TRUE,
      sorting = c("location_id", "sublocation_id", "winter")
    )
  file.path("hibernation", tolower(species), "analysis_data") |>
    update_metadata(
      root = analysis_data,
      name = sprintf("analaysis_data_%s", tolower(species)),
      title = sprintf(
        "Time series with the relevant observations of %s per sublocation.",
        this_species$scientific_name[1]
      ),
      description = sprintf("
Time series of the number of observed invidiuals for %s per winter. Every
sublocation has at least %i surveyed winter. The species is present during at
least %i winters in every sublocation. The dataset starts at %s.",
        this_species$scientific_name[1], n_winter, n_present, start
      ),
      field_description = c(
        location_id = "The identifier of the location.",
        sublocation_id = "The identifier of the sublocation.",
        winter = "The year of January 1st during the winter.",
        sample_id = "The identifier of the sample in case of a positive number.
Negative numbers are used for sublocations without samples.
They are a combination of the sublocation_id and the winter.",
        number = "The number of observed individuals."
      )
    )
  dataset |>
    distinct(.data$location_id, .data$X, .data$Y) |>
    write_vc(
      file.path("hibernation", tolower(species), "locations"),
      optimize = FALSE, root = analysis_data, stage = TRUE,
      sorting = c("location_id", "X", "Y")
    )
  file.path("hibernation", tolower(species), "locations") |>
    update_metadata(
      root = analysis_data,
      name = sprintf("locations_%s", tolower(species)),
      title = sprintf(
        "Locations of the relevant observations of %s per winter.",
        this_species$scientific_name[1]
      ),
      description = sprintf("
The relevant locations for %s per winter. Every sublocation has at least %i
surveyed winter. The species is present during at least %i winters in every
sublocation. The dataset starts at %s. The coordinates are in the Belgian
Lambert 72 coordinate system expressed as kilometers.",
        this_species$scientific_name[1], n_winter, n_present, start
      ),
      field_description = c(
        location_id = "The identifier of the location.",
        X = "X coordinate of the location.",
        Y = "Y coordinate of the location."
      )
    )
}

#' Convert raw data to analysis data for all species
#' @inheritParams select_imputation_section
#' @inheritParams prepare_analysis_data_species
#' @importFrom git2rdata commit
#' @export
prepare_analysis_data <- function(
  raw_data, analysis_data, start, n_winter = 4, n_present = 3
) {
  prepare_analysis_data_species(
    raw_data = raw_data, analysis_data = analysis_data, start = start,
    n_winter = n_winter, n_present = n_present, species = "Mbec"
  )
  prepare_analysis_data_species(
    raw_data = raw_data, analysis_data = analysis_data, start = start,
    n_winter = n_winter, n_present = n_present, species = "Mdas"
  )
  prepare_analysis_data_species(
    raw_data = raw_data, analysis_data = analysis_data, start = start,
    n_winter = n_winter, n_present = n_present, species = "Mdau"
  )
  prepare_analysis_data_species(
    raw_data = raw_data, analysis_data = analysis_data, start = start,
    n_winter = n_winter, n_present = n_present, species = "Mema"
  )
  prepare_analysis_data_species(
    raw_data = raw_data, analysis_data = analysis_data, start = start,
    n_winter = n_winter, n_present = n_present, species = "Mmysbra"
  )
  prepare_analysis_data_species(
    raw_data = raw_data, analysis_data = analysis_data, start = start,
    n_winter = n_winter, n_present = n_present, species = "Mmyo"
  )
  prepare_analysis_data_species(
    raw_data = raw_data, analysis_data = analysis_data, start = start,
    n_winter = n_winter, n_present = n_present, species = "Mnat"
  )
  prepare_analysis_data_species(
    raw_data = raw_data, analysis_data = analysis_data, start = start,
    n_winter = n_winter, n_present = n_present, species = "Pauraus"
  )
  prepare_analysis_data_species(
    raw_data = raw_data, analysis_data = analysis_data, start = start,
    n_winter = n_winter, n_present = n_present, species = "Pipspec"
  )
  commit(
    message = "Automated commit from abvanalysis", repo = analysis_data,
    session = TRUE, all = TRUE
  )
}
