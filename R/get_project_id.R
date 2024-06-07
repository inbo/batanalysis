#' Get project id for the monitoring of hibernating bats
#' @param origin A DBI connection to the database
#' @export
#' @importFrom assertthat assert_that
#' @importFrom DBI dbGetQuery
get_project_id_hibernation <- function(origin) {
  assert_that(inherits(origin, "Microsoft SQL Server"))
  dbGetQuery(
    conn = origin,
    "SELECT id
FROM staging_Meetnetten.projects_project
WHERE name = 'Vleermuizen - Wintertellingen'"
  )
}
