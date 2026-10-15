#' Import the location list
#' @inheritParams import_raw_data
#' @export
#' @importFrom assertthat assert_that
#' @importFrom DBI dbGetQuery
#' @importFrom dplyr transmute
#' @importFrom rlang .data
#' @importFrom sf st_as_sf st_drop_geometry st_join st_read st_transform
#' @importFrom sf st_within
read_raw_location <- function(origin) {
  assert_that(inherits(origin, "Microsoft SQL Server"))
  "SELECT
  l.id, l.code, l.name, l.parent_id, ROUND(l.geom.STX, 5) AS longitude,
  ROUND(l.geom.STY, 5) AS latitude
FROM staging_Meetnetten.locations_location AS l
INNER JOIN staging_Meetnetten.projects_projectlocation AS pl
  ON pl.location_id = l.id
INNER JOIN staging_Meetnetten.projects_project AS p ON p.id = pl.project_id
WHERE p.name = 'Vleermuizen - Wintertellingen'" |>
    dbGetQuery(conn = origin) |>
    st_as_sf(
      coords = c("longitude", "latitude"),
      crs = 4326,
      remove = FALSE
    ) -> location
  paste0(
    "WFS:https://ccff02.minfin.fgov.be/geoservices/arcgis/services/WMS/",
    "Cadastral_LayersWFS/MapServer/WFSServer?request=GetCapabilities&",
    "service=WFS&version=1.0.0"
  ) |>
    st_read(layer = "CL:Province", quiet = TRUE) |>
    st_transform(crs = 4326) |>
    transmute(province = factor(.data$NameDUT)) |>
    st_join(x = location, join = st_within) |>
    st_drop_geometry()
}
