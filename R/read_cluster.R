#' Import object clusters
#' @param cluster the name of the comma separated cluster file
#' @inheritParams import_raw_data
#' @export
#' @importFrom assertthat assert_that has_name noNA
#' @importFrom git2rdata read_vc
#' @importFrom stats aggregate
#' @importFrom utils file_test
read_cluster <- function(cluster, target) {
  stopifnot(requireNamespace("readr"))
  assert_that(file_test("-f", cluster))
  cluster <- readr::read_csv2(cluster)
  assert_that(
    has_name(cluster, "cluster"), has_name(cluster, "code"),
    is.numeric(cluster$cluster), noNA(cluster$cluster), noNA(cluster$code)
  )
  cluster$cluster <- as.integer(cluster$cluster)
  cluster$code <- as.character(cluster$code)
  if (!has_name(cluster, "eurobats")) {
    cluster$eurobats <- TRUE
  } else {
    assert_that(is.logical(cluster$eurobats), noNA(cluster$eurobats))
    valid <- aggregate(eurobats ~ cluster, data = cluster, FUN = range)
    valid_ok <- apply(valid$eurobats, 1, diff) == 0
    paste(valid$cluster[!valid_ok], collapse = "; ") |>
      sprintf(fmt = "clusters with different eurobats status: %s") -> msg
    assert_that(all(valid_ok), msg = msg)
  }
  location <- read_vc("hibernation/locations", root = target)
  cluster |>
    anti_join(location, by = "code") -> missing_code
  paste(missing_code$code, collapse = "; ") |>
    sprintf(
      fmt = paste(
        "codes in `cluster` without matching `code` in `hibernation/location`:",
        "%s", sep = "\n"
      )
    ) -> msg
  assert_that(nrow(missing_code) == 0, msg = msg)
  cluster |>
    select("cluster", "code") |>
    inner_join(
      location |>
        select("code", "id"), by = "code"
    ) |>
    select(-"code") |>
    write_vc(
      "hibernation/cluster_location", root = target, sorting = "id",
      stage = TRUE, force = TRUE
    )
  cluster |>
    distinct(.data$cluster, .data$eurobats) |>
    write_vc(
      "hibernation/cluster", root = target, sorting = "cluster",
      stage = TRUE, force = TRUE
    )
}
