#' Import the location list
#' @inheritParams import_raw_data
#' @export
#' @importFrom assertthat assert_that
#' @importFrom DBI dbGetQuery
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
    dbGetQuery(conn = origin)
}
