#' @export
#' @importFrom assertthat assert_that
#' @importFrom dplyr bind_cols case_when group_by inner_join mutate select
#' @importFrom dplyr starts_with summarise
#' @importFrom git2rdata is_git2rdata update_metadata write_vc
#' @importFrom n2kanalysis get_file_fingerprint status
#' @importFrom rlang .data
#' @importFrom stats quantile
#' @importFrom tidyr pivot_longer
extract_results.n2kAggregate <- function(x, root, ...) {
  assert_that(inherits(x, "n2kAggregate"))
  if (status(x) != "converged") {
    return(invisible(NULL))
  }
  x@AggregatedImputed@Covariate |>
    bind_cols(x@AggregatedImputed@Imputation) |>
    pivot_longer(
      starts_with("Imputation"),
      names_to = "imputation",
      values_to = "total"
    ) -> results
  if (grepl("location_id", x@AnalysisMetadata$formula)) {
    x@AnalysisMetadata |>
      select(
        species = "species_group_id",
        "model_type",
        analysis = "file_fingerprint",
        fingerprint = "status_fingerprint"
      ) |>
      bind_cols(
        results |>
          group_by(.data$winter, .data$imputation) |>
          mutate(rank = rank(-.data$total, ties.method = "first")) |>
          group_by(.data$winter, .data$location_id) |>
          summarise(
            mean = mean(.data$rank, na.rm = TRUE),
            min = min(.data$rank, na.rm = TRUE),
            p05 = quantile(.data$rank, 0.05, na.rm = TRUE),
            p20 = quantile(.data$rank, 0.2, na.rm = TRUE),
            p35 = quantile(.data$rank, 0.35, na.rm = TRUE),
            p50 = quantile(.data$rank, 0.5, na.rm = TRUE),
            p65 = quantile(.data$rank, 0.65, na.rm = TRUE),
            p80 = quantile(.data$rank, 0.8, na.rm = TRUE),
            p95 = quantile(.data$rank, 0.95, na.rm = TRUE),
            max = max(.data$rank, na.rm = TRUE),
            .groups = "drop"
          )
      ) |>
      write_vc(
        file = file.path("hibernation", "rank_location"),
        root = root,
        optimize = FALSE,
        append = TRUE,
        sorting = c(
          "model_type",
          "species",
          "winter",
          "location_id",
          "analysis"
        ),
        digits = 6
      )
    update_metadata(
      file = file.path("hibernation", "rank_location"),
      root = root,
      name = "hibernation_rank_location",
      title = paste(
        "Rank of the location in Flanders"
      ),
      description = paste(
        "The rank of the location in Flanders",
        "Missing values are imputed before calculating the rank,",
        "except for Mspec, Rspec, Eser and Rfer."
      ),
      field_description = c(
        species = "The code of the species group",
        model_type = paste(
          "A short description of the model used to calculate the totals"
        ),
        analysis = "The file fingerprint of the analysis",
        fingerprint = "The status fingerprint of the analysis",
        winter = "The winter season defined by the year in which January falls",
        location_id = "The unique identifier of the location",
        mean = "The average of rank",
        min = "The minimum rank",
        p05 = "The 5% quantile of the rank",
        p20 = "The 20% quantile of the rank",
        p35 = "The 35% quantile of the rank",
        p50 = "The 50% quantile of the rank",
        p65 = "The 65% quantile of the rank",
        p80 = "The 80% quantile of the rank",
        p95 = "The 95% quantile of the rank",
        max = "The maximum of the rank"
      )
    )

    x@AnalysisMetadata |>
      select(
        species = "species_group_id",
        "model_type",
        analysis = "file_fingerprint",
        fingerprint = "status_fingerprint"
      ) |>
      bind_cols(
        results |>
          group_by(.data$winter, .data$location_id) |>
          summarise(
            mean = mean(.data$total, na.rm = TRUE),
            min = min(.data$total, na.rm = TRUE) |>
              as.numeric(),
            p05 = quantile(.data$total, 0.05, na.rm = TRUE),
            p20 = quantile(.data$total, 0.2, na.rm = TRUE),
            p35 = quantile(.data$total, 0.35, na.rm = TRUE),
            p50 = quantile(.data$total, 0.5, na.rm = TRUE),
            p65 = quantile(.data$total, 0.65, na.rm = TRUE),
            p80 = quantile(.data$total, 0.8, na.rm = TRUE),
            p95 = quantile(.data$total, 0.95, na.rm = TRUE),
            max = max(.data$total, na.rm = TRUE) |>
              as.numeric(),
            .groups = "drop"
          )
      ) |>
      write_vc(
        file = file.path("hibernation", "total_location"),
        root = root,
        optimize = FALSE,
        append = TRUE,
        sorting = c(
          "model_type",
          "species",
          "winter",
          "location_id",
          "analysis"
        ),
        digits = 6
      )
    update_metadata(
      file = file.path("hibernation", "total_location"),
      root = root,
      name = "hibernation_total_location",
      title = "The imputed total number of hibernating bats per location",
      description = paste(
        "The imputed total number of hibernating bats in the winter season per",
        "location.",
        "Missing values are imputed before calculating the total,",
        "except for Mspec, Rspec, Eser and Rfer."
      ),
      field_description = c(
        species = "The code of the species group",
        model_type = paste(
          "A short description of the model used to calculate the totals"
        ),
        analysis = "The file fingerprint of the analysis",
        fingerprint = "The status fingerprint of the analysis",
        winter = "The winter season defined by the year in which January falls",
        location_id = "The unique identifier of the location",
        mean = "The average of the imputed total number of hibernating bats",
        min = "The minimum of the imputed total number of hibernating bats",
        p05 = "The 5% quantile of the imputed total number of hibernating bats",
        p20 = paste(
          "The 20% quantile of the imputed total number of hibernating bats"
        ),
        p35 = paste(
          "The 35% quantile of the imputed total number of hibernating bats"
        ),
        p50 = paste(
          "The 50% quantile of the imputed total number of hibernating bats"
        ),
        p65 = paste(
          "The 65% quantile of the imputed total number of hibernating bats"
        ),
        p80 = paste(
          "The 80% quantile of the imputed total number of hibernating bats"
        ),
        p95 = paste(
          "The 95% quantile of the imputed total number of hibernating bats"
        ),
        max = "The maximum of the imputed total number of hibernating bats"
      )
    )
    rm(x, results)
    gc(verbose = FALSE)
    return(invisible(NULL))
  }
  if (grepl("type", x@AnalysisMetadata$formula)) {
    x@AnalysisMetadata |>
      select(
        species = "species_group_id",
        "model_type",
        analysis = "file_fingerprint",
        fingerprint = "status_fingerprint"
      ) |>
      bind_cols(
        results |>
          group_by(.data$winter, .data$type) |>
          summarise(
            mean = mean(.data$total, na.rm = TRUE),
            min = min(.data$total, na.rm = TRUE) |>
              as.numeric(),
            p05 = quantile(.data$total, 0.05, na.rm = TRUE),
            p20 = quantile(.data$total, 0.2, na.rm = TRUE),
            p35 = quantile(.data$total, 0.35, na.rm = TRUE),
            p50 = quantile(.data$total, 0.5, na.rm = TRUE),
            p65 = quantile(.data$total, 0.65, na.rm = TRUE),
            p80 = quantile(.data$total, 0.8, na.rm = TRUE),
            p95 = quantile(.data$total, 0.95, na.rm = TRUE),
            max = max(.data$total, na.rm = TRUE) |>
              as.numeric(),
            .groups = "drop"
          )
      ) |>
      write_vc(
        file = file.path("hibernation", "total_type"),
        root = root,
        optimize = FALSE,
        append = TRUE,
        sorting = c("model_type", "species", "winter", "type", "analysis"),
        digits = 6
      )
    update_metadata(
      file = file.path("hibernation", "total_type"),
      root = root,
      name = "hibernation_total_location",
      title = paste(
        "The imputed total number of hibernating bats per type of location"
      ),
      description = paste(
        "The imputed total number of hibernating bats in the winter season per",
        "type of location.",
        "Missing values are imputed before calculating the total,",
        "except for Mspec, Rspec, Eser and Rfer."
      ),
      field_description = c(
        species = "The code of the species group",
        model_type = paste(
          "A short description of the model used to calculate the totals"
        ),
        analysis = "The file fingerprint of the analysis",
        fingerprint = "The status fingerprint of the analysis",
        winter = "The winter season defined by the year in which January falls",
        type = "The type of location",
        mean = "The average of the imputed total number of hibernating bats",
        min = "The minimum of the imputed total number of hibernating bats",
        p05 = "The 5% quantile of the imputed total number of hibernating bats",
        p20 = paste(
          "The 20% quantile of the imputed total number of hibernating bats"
        ),
        p35 = paste(
          "The 35% quantile of the imputed total number of hibernating bats"
        ),
        p50 = paste(
          "The 50% quantile of the imputed total number of hibernating bats"
        ),
        p65 = paste(
          "The 65% quantile of the imputed total number of hibernating bats"
        ),
        p80 = paste(
          "The 80% quantile of the imputed total number of hibernating bats"
        ),
        p95 = paste(
          "The 95% quantile of the imputed total number of hibernating bats"
        ),
        max = "The maximum of the imputed total number of hibernating bats"
      )
    )
    rm(x, results)
    gc(verbose = FALSE)
    return(invisible(NULL))
  }
  if (grepl("province", x@AnalysisMetadata$formula)) {
    x@AnalysisMetadata |>
      select(
        species = "species_group_id",
        "model_type",
        analysis = "file_fingerprint",
        fingerprint = "status_fingerprint"
      ) |>
      bind_cols(
        results |>
          group_by(.data$winter, .data$province) |>
          summarise(
            mean = mean(.data$total, na.rm = TRUE),
            min = min(.data$total, na.rm = TRUE) |>
              as.numeric(),
            p05 = quantile(.data$total, 0.05, na.rm = TRUE),
            p20 = quantile(.data$total, 0.2, na.rm = TRUE),
            p35 = quantile(.data$total, 0.35, na.rm = TRUE),
            p50 = quantile(.data$total, 0.5, na.rm = TRUE),
            p65 = quantile(.data$total, 0.65, na.rm = TRUE),
            p80 = quantile(.data$total, 0.8, na.rm = TRUE),
            p95 = quantile(.data$total, 0.95, na.rm = TRUE),
            max = max(.data$total, na.rm = TRUE) |>
              as.numeric(),
            .groups = "drop"
          )
      ) |>
      write_vc(
        file = file.path("hibernation", "total_province"),
        root = root,
        optimize = FALSE,
        append = TRUE,
        sorting = c("model_type", "species", "winter", "province", "analysis"),
        digits = 6
      )
    update_metadata(
      file = file.path("hibernation", "total_province"),
      root = root,
      name = "hibernation_total_province",
      title = paste(
        "The imputed total number of hibernating bats per province"
      ),
      description = paste(
        "The imputed total number of hibernating bats in the winter season per",
        "province.",
        "Missing values are imputed before calculating the total,",
        "except for Mspec, Rspec, Eser and Rfer."
      ),
      field_description = c(
        species = "The code of the species group",
        model_type = paste(
          "A short description of the model used to calculate the totals"
        ),
        analysis = "The file fingerprint of the analysis",
        fingerprint = "The status fingerprint of the analysis",
        winter = "The winter season defined by the year in which January falls",
        province = "The province",
        mean = "The average of the imputed total number of hibernating bats",
        min = "The minimum of the imputed total number of hibernating bats",
        p05 = "The 5% quantile of the imputed total number of hibernating bats",
        p20 = paste(
          "The 20% quantile of the imputed total number of hibernating bats"
        ),
        p35 = paste(
          "The 35% quantile of the imputed total number of hibernating bats"
        ),
        p50 = paste(
          "The 50% quantile of the imputed total number of hibernating bats"
        ),
        p65 = paste(
          "The 65% quantile of the imputed total number of hibernating bats"
        ),
        p80 = paste(
          "The 80% quantile of the imputed total number of hibernating bats"
        ),
        p95 = paste(
          "The 95% quantile of the imputed total number of hibernating bats"
        ),
        max = "The maximum of the imputed total number of hibernating bats"
      )
    )
    rm(x, results)
    gc(verbose = FALSE)
    return(invisible(NULL))
  }
  stopifnot(x@AnalysisMetadata$formula == "~winter")
  x@AnalysisMetadata |>
    select(
      species = "species_group_id",
      "model_type",
      analysis = "file_fingerprint",
      fingerprint = "status_fingerprint"
    ) |>
    bind_cols(
      results |>
        group_by(.data$winter) |>
        summarise(
          mean = mean(.data$total, na.rm = TRUE),
          min = min(.data$total, na.rm = TRUE) |>
            as.numeric(),
          p05 = quantile(.data$total, 0.05, na.rm = TRUE),
          p20 = quantile(.data$total, 0.2, na.rm = TRUE),
          p35 = quantile(.data$total, 0.35, na.rm = TRUE),
          p50 = quantile(.data$total, 0.5, na.rm = TRUE),
          p65 = quantile(.data$total, 0.65, na.rm = TRUE),
          p80 = quantile(.data$total, 0.8, na.rm = TRUE),
          p95 = quantile(.data$total, 0.95, na.rm = TRUE),
          max = max(.data$total, na.rm = TRUE) |>
            as.numeric()
        )
    ) |>
    write_vc(
      file = file.path("hibernation", "total"),
      root = root,
      optimize = FALSE,
      append = TRUE,
      sorting = c("model_type", "species", "winter", "analysis"),
      digits = 6
    )
  update_metadata(
    file = file.path("hibernation", "total"),
    root = root,
    name = "hibernation_total",
    title = "The imputed total number of hibernating bats",
    description = paste(
      "The imputed total number of hibernating bats in the winter season.",
      "Missing values are imputed before calculating the total,",
      "except for Mspec, Rspec, Eser and Rfer."
    ),
    field_description = c(
      species = "The code of the species group",
      model_type = paste(
        "A short description of the model used to calculate the totals"
      ),
      analysis = "The file fingerprint of the analysis",
      fingerprint = "The status fingerprint of the analysis",
      winter = "The winter season defined by the year in which January falls",
      mean = "The average of the imputed total number of hibernating bats",
      min = "The minimum of the imputed total number of hibernating bats",
      p05 = "The 5% quantile of the imputed total number of hibernating bats",
      p20 = paste(
        "The 20% quantile of the imputed total number of hibernating bats"
      ),
      p35 = paste(
        "The 35% quantile of the imputed total number of hibernating bats"
      ),
      p50 = paste(
        "The 50% quantile of the imputed total number of hibernating bats"
      ),
      p65 = paste(
        "The 65% quantile of the imputed total number of hibernating bats"
      ),
      p80 = paste(
        "The 80% quantile of the imputed total number of hibernating bats"
      ),
      p95 = paste(
        "The 95% quantile of the imputed total number of hibernating bats"
      ),
      max = "The maximum of the imputed total number of hibernating bats"
    )
  )
  rm(x, results)
  gc(verbose = FALSE)
  return(invisible(NULL))
}
