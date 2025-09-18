library(batanalysis)
library(git2rdata)
library(quarto)
library(tidyverse)
dataroot <- keyring::key_get(service = "result_path", keyring = "batanalysis")
resultroot <- keyring::key_get(
  service = "target_path",
  keyring = "batanalysis_data"
)

relevant_visit <- list_relevant_visit(data_root)
write_vc(
  x = relevant_visit,
  file = "relevant_visit",
  root = ".",
  sorting = "visit_id",
  digits = 10
)

relevant_visit |>
  filter(.data$date >= Sys.Date() - 24 * 365) |>
  slice_max(.data$total, n = 1, with_ties = FALSE, by = "location_id") |>
  inner_join(
    file.path("data", "hibernation", "locations") |>
      verify_vc(root = dataroot, variables = c("id", "name")),
    by = c("location_id" = "id")
  ) |>
  arrange(desc(.data$total)) |>
  transmute(
    .data$location_id,
    code = ifelse(
      .data$code == "",
      str_remove(.data$name, "\\s*-.*"),
      .data$code
    ),
    name = str_remove(.data$name, "^[0-9]{4}(\\s*\\(.*?\\))?\\s*-\\s*"),
    filename = ifelse(.data$name == "", .data$code, .data$name) |>
      tolower() |>
      str_replace_all(" ", "-") |>
      str_replace_all("-+", "-") |>
      str_remove_all("[,;/'\\(\\)\\.]") |>
      sprintf(fmt = "%2$i-%1$s", .data$location_id),
    name = ifelse(
      .data$name == "",
      .data$code,
      sprintf("%s (%s)", .data$name, .data$code)
    ),
    .data$total
  ) -> to_do

dir.create("reports", showWarnings = FALSE)
template <- readLines("_quarto.yml")
for (i in seq_len(nrow(to_do))) {
  message(to_do$name[i])
  target <- file.path(
    "reports",
    gsub("-", "_", paste0(to_do$filename[i], ".pdf"))
  )
  if (file.exists(target)) {
    next
  }
  c(
    template,
    sprintf("  title: %s\n  shorttitle: %s", to_do$name[i], to_do$filename[i])
  ) |>
    writeLines("_quarto.yml")
  quarto_render(
    ".",
    quiet = FALSE,
    as_job = FALSE,
    execute_params = list(
      data = dataroot,
      result = resultroot,
      location = to_do$location_id[i],
      name = to_do$name[i]
    )
  )
  list.files("output", pattern = ".pdf$", full.names = TRUE) |>
    file.copy(target, overwrite = TRUE)
}
