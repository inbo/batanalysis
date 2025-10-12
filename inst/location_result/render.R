library(batanalysis)
library(git2rdata)
library(quarto)
library(tidyverse)
dataroot <- keyring::key_get("meetnetten", "batanalysis_data")
resultroot <- keyring::key_get("meetnetten", "batanalysis_result")

relevant_visit <- visit_type(dataroot)
write_vc(
  x = relevant_visit,
  file = "visit_type",
  root = ".",
  sorting = "visit_id",
  digits = 5,
  optimize = FALSE,
  strict = FALSE
)

relevant_visit |>
  filter(!is.na(.data$score)) |>
  group_by(.data$location_id) |>
  summarise(score = sum(.data$score), small = mean(is.na(.data$level))) |>
  inner_join(
    file.path("data", "hibernation", "locations") |>
      verify_vc(root = dataroot, variables = c("id", "name")),
    by = c("location_id" = "id")
  ) |>
  arrange(desc(.data$score)) |>
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
    .data$score,
    detail = .data$small == 0
  ) -> to_do

dir.create("reports", showWarnings = FALSE)
file.path("location_result", "_quarto.yml") |>
  system.file(package = "batanalysis") |>
  readLines() -> template
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
    paste(
      "  title: %s",
      paste(
        "  subtitle: Overzicht de overwinterende vleermuizen tijdens de",
        "periode 2002-2025"
      ),
      "  shorttitle: %s",
      sep = "\n"
    ) |>
      sprintf(to_do$name[i], to_do$filename[i]),
    "  detail: true"[to_do$detail[i]]
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
