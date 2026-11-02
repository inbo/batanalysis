#' @export
#' @importFrom assertthat assert_that is.flag is.string noNA
#' @importFrom dplyr anti_join distinct
#' @importFrom git2rdata rm_data verify_vc write_vc
#' @importFrom n2kanalysis order_manifest read_manifest read_model
#' @importFrom purrr walk
extract_results.character <- function(
  x,
  base,
  project = "batanalysis",
  raw_data,
  root,
  random = FALSE,
  ...
) {
  assert_that(is.string(x), noNA(x), is.flag(random), noNA(random))
  rm_data(root = root, path = "hibernation", type = "all")
  rm_data(
    root = root, path = file.path("model_check", "hibernation"), type = "all"
  )
  if (!missing(raw_data)) {
    file.path("data", "hibernation", "species") |>
      verify_vc(
        root = raw_data,
        variables = c("id", "code", "name", "scientific_name", "parent")
      ) |>
      write_vc(
        file.path("hibernation", "species"),
        root = root,
        sorting = "id",
        optimize = FALSE
      )
    update_metadata(
      file.path("hibernation", "species"),
      root = root,
      name = "hibernation_species",
      title = "Hibernating bat species",
      description = paste(
        "List of species potentially relevant for the hibernating bat",
        "monitoring in Flanders (Belgium)."
      ),
      field_description = c(
        id = "Unique identifier of the species",
        name = "Dutch vernacular name",
        scientific_name = "Scientific name of the species",
        code = "Code of the species",
        parent = "Unique identifier of the parent species"
      )
    )
    file.path("data", "hibernation", "locations") |>
      verify_vc(
        root = raw_data,
        variables = c("id", "code", "name", "parent_id")
      ) |>
      write_vc(
        file.path("hibernation", "locations"),
        root = root,
        sorting = "id",
        optimize = FALSE,
        digits = c("longitude" = 8, "latitude" = 8)
      )
    update_metadata(
      file.path("hibernation", "locations"),
      root = root,
      name = "hibernation_locations",
      title = paste(
        "Locations and sublocations surveyed during the hibernating bat",
        "monitoring"
      ),
      description = paste(
        "List of locations and sublocations surveyed during the hibernating",
        "bat monitoring in Flanders (Belgium)."
      ),
      field_description = c(
        id = "Unique identifier of the location or sublocation",
        name = "Name of the location or sublocation",
        parent_id = paste(
          "Unique identifier of the parent location (-1 for locations)"
        ),
        code = "Code of the location or sublocation",
        longitude = "Longitude of the location or sublocation (EPSG:4326)",
        latitude = "Latitude of the location or sublocation (EPSG:4326)"
      )
    )
  }
  read_manifest(base = base, project = project, hash = x) |>
    order_manifest() -> manifest
  for (type in c(
    "hurdle",
    "total",
    "total_location",
    "total_type",
    "total_rw1"
  )) {
    filename <- file.path("hibernation", type)
    if (is_git2rdata(filename, root = root)) {
      manifest <- manifest[
        !manifest %in%
          verify_vc(filename, root = root, "analysis")$analysis
      ]
    }
  }
  filename <- file.path("model_check", "hibernation", "hyperparameters")
  if (is_git2rdata(filename, root = root)) {
    manifest <- manifest[
      !manifest %in%
        verify_vc(filename, root = root, "analysis")$analysis
    ]
  }
  if (random) {
    manifest <- sample(manifest)
  } else if (inherits(base, "character")) {
    file.path(base, project) |>
      list.files(pattern = ".rds$", full.names = TRUE, recursive = TRUE) |>
      grepv(pattern = "converged") |>
      file.info() |>
      arrange(.data$size) |>
      rownames() |>
      basename() |>
      gsub(pattern = ".rds$", replacement = "") -> to_do
    to_do[to_do %in% manifest] -> manifest
  }
  for (i in manifest) {
    message(i)
    read_model(i, base = base, project = project) |>
      extract_results(root = root)
    gc(verbose = FALSE)
  }
  prune_meta(root = root, path = "hibernation")
  prune_meta(root = root, path = file.path("model_check", "hibernation"))
  return(invisible(NULL))
}
