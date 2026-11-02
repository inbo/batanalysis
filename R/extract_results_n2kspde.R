#' @export
#' @importFrom assertthat assert_that is.count noNA
#' @importFrom dplyr arrange bind_cols distinct group_by inner_join mutate
#' @importFrom dplyr row_number select starts_with summarise
#' @importFrom fmesher fm_evaluator
#' @importFrom git2rdata is_git2rdata update_metadata write_vc
#' @importFrom INLA inla.posterior.sample inla.zmarginal
#' @importFrom inlatools prec2sd
#' @importFrom n2kanalysis get_file_fingerprint spde2mesh status
#' @importFrom purrr map map2_dfc
#' @importFrom sf st_area st_as_sf st_convex_hull st_sample st_union
#' @importFrom stringr str_detect str_remove str_replace_all
#' @importFrom stats quantile
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr everything pivot_longer pivot_wider unnest
extract_results.n2kSpde <- function(x, root, ..., n_sim = 100) {
  assert_that(inherits(x, "n2kSpde"), is.count(n_sim), noNA(n_sim))
  # skip if model is not converged
  if (status(x) != "converged") {
    return(invisible(NULL))
  }
  # generate posterior sample
  if (Sys.getenv("POSITRON") == "") {
    num_threads <- NULL
  } else {
    num_threads <- 1
  }
  inla.posterior.sample(
    n = n_sim, result = x@Model, num.threads = num_threads
  ) |>
    map("latent") |>
    map2_dfc(
      seq_len(n_sim),
      ~ sprintf("sim_%04i", .y) |>
        `colnames<-`(.x, value = _) |>
        as.data.frame()
    ) |>
    rownames_to_column(var = "parameter") |>
    pivot_longer(
      -"parameter",
      names_to = "sim",
      values_to = "estimate"
    ) -> post_sample
  x@Data |>
    distinct(
      .data$location_id,
      .data$fortress,
      .data$marl_quarry,
      .data$other_large,
      .data$small
    ) |>
    transmute(
      .data$location_id,
      type = case_when(
        .data$fortress > 0 ~ "fortress",
        .data$marl_quarry > 0 ~ "marl quarry",
        .data$other_large > 0 ~ "other large",
        .data$small > 0 ~ "small"
      ) |>
        factor(c("fortress", "marl quarry", "other large", "small"))
    ) -> location_type
  # extract the intercept and the winter effect
  post_sample |>
    filter(
      .data$parameter %in%
        c("fortress:1", "marl_quarry:1", "other_large:1", "small:1")
    ) |>
    transmute(
      type = str_remove(.data$parameter, ":1") |>
        str_replace_all("_", " ") |>
        factor(levels = levels(location_type$type)),
      .data$sim,
      intercept = .data$estimate
    ) -> ps_intercept
  post_sample |>
    filter(str_detect(.data$parameter, "winter_i")) |>
    transmute(
      winter_r = str_remove(.data$parameter, "winter_i:") |>
        as.integer(),
      .data$sim,
      .data$estimate
    ) |>
    inner_join(ps_intercept, by = "sim", relationship = "many-to-many") |>
    transmute(
      winter = .data$winter_r + min(x@Data$winter) - 1,
      .data$type,
      .data$sim,
      estimate = .data$intercept + .data$estimate
    ) -> ps_cwinter
  x@AnalysisMetadata |>
    select(
      species = "species_group_id",
      "model_type",
      analysis = "file_fingerprint",
      fingerprint = "status_fingerprint"
    ) |>
    bind_cols(
      ps_cwinter |>
        group_by(.data$type, .data$winter) |>
        summarise(
          across(
            "estimate",
            .names = "{.fn}",
            list(
              mean = ~ mean(.x, na.rm = TRUE) |>
                round(4),
              sd = ~ sd(.x, na.rm = TRUE) |>
                round(4)
            )
          ),
          .groups = "drop"
        )
    ) |>
    write_vc(
      file = file.path("model_check", "hibernation", "winter"),
      root = root,
      append = TRUE,
      sorting = c("model_type", "species", "winter", "type", "analysis"),
      digits = 4
    )
  update_metadata(
    file = file.path("model_check", "hibernation", "winter"),
    root = root,
    name = "hibernation_model_check_winter",
    title = "Modelled winter effect",
    description = paste(
      "Posterior distribution of the modelled winter effect at the location",
      "level.",
      "The model is a binomial spatio-temporal model with a logit link.",
      "The winter effect is modeled with a second order random walk."
    ),
    field_description = c(
      species = "The code of the species group",
      model_type = paste(
        "A short description of the model used to analyse the data"
      ),
      analysis = "The file fingerprint of the analysis",
      fingerprint = "The status fingerprint of the analysis",
      winter = "The winter season defined by the year in which January falls",
      type = "The type of location",
      mean = paste(
        "The mean of the posterior distribution on the logit scale for",
        "binomial models and the log scale for nbinomial models"
      ),
      sd = paste(
        "The standard deviation of the posterior distribution on the logit",
        "scale for binomial models and the log scale for nbinomial models"
      )
    )
  )
  # extract the location specific effects
  post_sample |>
    filter(str_detect(.data$parameter, "^matern_spde:")) |>
    mutate(
      parameter = str_remove(.data$parameter, "matern_spde:") |>
        as.integer()
    ) |>
    pivot_wider(names_from = "sim", values_from = "estimate") |>
    arrange(.data$parameter) |>
    select(-"parameter") |>
    as.matrix() -> ps_matern
  x@Model$.args$data[c("location_id", "X", "Y")] |>
    as.data.frame() |>
    distinct() |>
    filter(!is.na(.data$X)) -> loc_coordinates
  projector <- fm_evaluator(
    mesh = spde2mesh(x@Spde),
    loc = as.matrix(loc_coordinates[, c("X", "Y")])
  )
  if (grepl("present", x@AnalysisMetadata$formula)) {
    x@Data |>
      distinct(
        .data$location_id,
        location_i = as.integer(.data$location_i),
        .data$winter,
        .data$winter_l,
        .data$winter_q,
        .data$winter_c
      ) -> location_winter
    location_winter |>
      mutate(
        across(c("winter_l", "winter_q", "winter_c"), ~ replace_na(.x, 0))
      ) -> lw
    post_sample |>
      filter(str_detect(.data$parameter, "^location_i:")) |>
      transmute(
        location_i = str_remove(.data$parameter, "location_i:") |>
          as.integer(),
        .data$sim,
        q0 = .data$estimate
      ) |>
      left_join(
        post_sample |>
          filter(str_detect(.data$parameter, "^location_l:")) |>
          transmute(
            location_i = str_remove(.data$parameter, "location_l:") |>
              as.integer(),
            .data$sim,
            q1 = .data$estimate
          ),
        by = c("location_i", "sim")
      ) |>
      left_join(
        post_sample |>
          filter(str_detect(.data$parameter, "^location_q:")) |>
          transmute(
            location_i = str_remove(.data$parameter, "location_q:") |>
              as.integer(),
            .data$sim,
            q2 = .data$estimate
          ),
        by = c("location_i", "sim")
      ) |>
      left_join(
        post_sample |>
          filter(str_detect(.data$parameter, "^location_c:")) |>
          transmute(
            location_i = str_remove(.data$parameter, "location_c:") |>
              as.integer(),
            .data$sim,
            q3 = .data$estimate
          ),
        by = c("location_i", "sim")
      ) |>
      inner_join(
        location_winter |>
          group_by(.data$location_i, .data$location_id) |>
          summarise(
            across(
              c("winter_l", "winter_q", "winter_c"),
              ~ !any(is.na(.x))
            ),
            .groups = "drop"
          ),
        by = "location_i"
      ) |>
      transmute(
        .data$location_id,
        .data$sim,
        .data$q0,
        q1 = ifelse(.data$winter_l, .data$q1, 0),
        q2 = ifelse(.data$winter_q, .data$q2, 0),
        q3 = ifelse(.data$winter_c, .data$q3, 0),
        across(c("q0", "q1", "q2", "q3"), ~ replace_na(.x, 0))
      ) -> ps_location
  } else {
    x@Data |>
      distinct(
        .data$location_id,
        location_i = as.integer(.data$location_i)
      ) -> location
    post_sample |>
      filter(str_detect(.data$parameter, "^location_i:")) |>
      transmute(
        location_i = str_remove(.data$parameter, "location_i:") |>
          as.integer(),
        .data$sim,
        q0 = .data$estimate
      ) |>
      left_join(
        post_sample |>
          filter(str_detect(.data$parameter, "^loc_1:")) |>
          transmute(
            location_i = str_remove(.data$parameter, "loc_1:") |>
              as.integer(),
            .data$sim,
            q1 = .data$estimate
          ),
        by = c("location_i", "sim")
      ) |>
      left_join(
        post_sample |>
          filter(str_detect(.data$parameter, "^loc_2:")) |>
          transmute(
            location_i = str_remove(.data$parameter, "loc_2:") |>
              as.integer(),
            .data$sim,
            q2 = .data$estimate
          ),
        by = c("location_i", "sim")
      ) |>
      left_join(
        post_sample |>
          filter(str_detect(.data$parameter, "^loc_3:")) |>
          transmute(
            location_i = str_remove(.data$parameter, "loc_3:") |>
              as.integer(),
            .data$sim,
            q3 = .data$estimate
          ),
        by = c("location_i", "sim")
      ) |>
      inner_join(location, by = "location_i") -> ps_location
    x@Data |>
      distinct(
        .data$location_id,
        .data$winter,
        winter_l = .data$loc_s_1,
        winter_q = .data$loc_s_2,
        winter_c = .data$loc_s_3
      ) -> lw
  }
  ps_location |>
    inner_join(
      as.matrix(projector$proj$A %*% ps_matern) |>
        as.data.frame() |>
        mutate(
          location_id = loc_coordinates$location_id
        ) |>
        pivot_longer(starts_with("sim"), names_to = "sim", values_to = "mesh"),
      by = c("location_id", "sim")
    ) -> ps_location

  x@AnalysisMetadata |>
    select(
      species = "species_group_id",
      "model_type",
      analysis = "file_fingerprint",
      fingerprint = "status_fingerprint"
    ) |>
    bind_cols(
      ps_location |>
        group_by(.data$location_id) |>
        summarise(
          across(
            c("q0", "q1", "q2", "q3", "mesh"),
            list(
              mean = ~ mean(.x, na.rm = TRUE),
              sd = ~ sd(.x, na.rm = TRUE)
            )
          )
        )
    ) |>
    write_vc(
      file = file.path("model_check", "hibernation", "location_effect"),
      root = root,
      append = TRUE,
      sorting = c("model_type", "species", "location_id", "analysis"),
      digits = 4
    )
  update_metadata(
    file = file.path("model_check", "hibernation", "location_effect"),
    root = root,
    name = "hibernation_model_check_location_effect",
    title = "Location specific effects",
    description = paste(
      "Posterior distribution of the location specific effects.",
      "The effects are on the logit scale for binomial models and the log",
      "scale for nbinomial models."
    ),
    field_description = c(
      species = "The code of the species group",
      model_type = paste(
        "A short description of the model used to analyse the data"
      ),
      analysis = "The file fingerprint of the analysis",
      fingerprint = "The status fingerprint of the analysis",
      location_id = "The unique identifier of the location",
      q0_mean = paste(
        "The mean of the intercept of the location specific effect on the",
        "logit scale for binomial models and the log scale for nbinomial models"
      ),
      q0_sd = paste(
        "The standard deviation of the intercept of the location specific",
        "effect on the logit scale for binomial models and the log scale for",
        "nbinomial models"
      ),
      q1_mean = paste(
        "The mean of the linear term of the location specific effect on the",
        "logit scale for binomial models and the log scale for nbinomial models"
      ),
      q1_sd = paste(
        "The standard deviation of the linear term of the location specific",
        "effect on the logit scale for binomial models and the log scale for",
        "nbinomial models"
      ),
      q2_mean = paste(
        "The mean of the quadratic term of the location specific effect on",
        "the logit scale for binomial models and the log scale for nbinomial",
        "models"
      ),
      q2_sd = paste(
        "The standard deviation of the quadratic term of the location",
        "specific effect on the logit scale for binomial models and the log",
        "scale for nbinomial models"
      ),
      q3_mean = paste(
        "The mean of the cubic term of the location specific effect on the",
        "logit scale for binomial models and the log scale for nbinomial models"
      ),
      q3_sd = paste(
        "The standard deviation of the cubic term of the location",
        "specific effect on the logit scale for binomial models and the log",
        "scale for nbinomial models"
      )
    )
  )
  # predictions at the location level
  lw |>
    inner_join(x = location_type, by = "location_id") |>
    inner_join(
      ps_location,
      by = "location_id",
      relationship = "many-to-many"
    ) |>
    inner_join(ps_cwinter, by = c("winter", "type", "sim")) |>
    transmute(
      .data$location_id,
      .data$winter,
      .data$sim,
      relative = .data$q0 +
        .data$q1 * .data$winter_l +
        .data$q2 * .data$winter_q +
        .data$q3 * .data$winter_c +
        .data$mesh,
      absolute = .data$estimate + .data$relative
    ) -> ps_location
  x@AnalysisMetadata |>
    select(
      species = "species_group_id",
      "model_type",
      analysis = "file_fingerprint",
      fingerprint = "status_fingerprint"
    ) |>
    bind_cols(
      ps_location |>
        group_by(.data$location_id, .data$winter) |>
        summarise(
          across(
            c("relative", "absolute"),
            list(
              mean = ~ mean(.x, na.rm = TRUE),
              sd = ~ sd(.x, na.rm = TRUE)
            )
          ),
          .groups = "drop"
        )
    ) |>
    write_vc(
      file = file.path("model_check", "hibernation", "location"),
      root = root,
      append = TRUE,
      sorting = c("model_type", "species", "location_id", "winter", "analysis"),
      digits = 4
    )
  update_metadata(
    file = file.path("model_check", "hibernation", "location"),
    root = root,
    name = "hibernation_model_check_location",
    title = "Predictions at the location level",
    description = paste(
      "Predictions of the model at the location level.",
      "The predictions are at the logit scale for binomial models and the log",
      "scale for nbinomial models."
    ),
    field_description = c(
      species = "The code of the species group",
      model_type = paste(
        "A short description of the model used to analyse the data"
      ),
      analysis = "The file fingerprint of the analysis",
      fingerprint = "The status fingerprint of the analysis",
      location_id = "The unique identifier of the location",
      winter = "The winter season defined by the year in which January falls",
      relative_mean = paste(
        "The mean of the difference between the local trend and the global",
        "trend"
      ),
      absolute_mean = paste(
        "The mean of the local trend"
      ),
      relative_sd = paste(
        "The standard error of the difference between the local trend and the",
        "global trend"
      ),
      absolute_sd = paste(
        "The standard error of the local trend"
      )
    )
  )

  # sublocation effect
  if (grepl("present", x@AnalysisMetadata$formula)) {
    x@Data |>
      filter(!is.na(.data$sublocation_id)) |>
      distinct(
        .data$location_id,
        sublocation_i = as.integer(.data$sublocation_i),
        .data$sublocation_id,
        .data$winter,
        .data$winter_l,
        .data$winter_q,
        .data$winter_c
      ) -> sublocation_winter
    sublocation_winter |>
      mutate(
        across(c("winter_l", "winter_q", "winter_c"), ~ replace_na(.x, 0))
      ) -> lw
    post_sample |>
      filter(str_detect(.data$parameter, "^sublocation_i:")) |>
      transmute(
        sublocation_i = str_remove(.data$parameter, "sublocation_i:"),
        .data$sim,
        q0 = .data$estimate
      ) |>
      left_join(
        post_sample |>
          filter(str_detect(.data$parameter, "^sublocation_l:")) |>
          transmute(
            sublocation_i = str_remove(.data$parameter, "sublocation_l:"),
            .data$sim,
            q1 = .data$estimate
          ),
        by = c("sublocation_i", "sim")
      ) |>
      left_join(
        post_sample |>
          filter(str_detect(.data$parameter, "^sublocation_q:")) |>
          transmute(
            sublocation_i = str_remove(.data$parameter, "sublocation_q:"),
            .data$sim,
            q2 = .data$estimate
          ),
        by = c("sublocation_i", "sim")
      ) |>
      left_join(
        post_sample |>
          filter(str_detect(.data$parameter, "^sublocation_c:")) |>
          transmute(
            sublocation_i = str_remove(.data$parameter, "sublocation_c:"),
            .data$sim,
            q3 = .data$estimate
          ),
        by = c("sublocation_i", "sim")
      ) |>
      mutate(
        sublocation_i = as.integer(.data$sublocation_i)
      ) |>
      inner_join(
        sublocation_winter |>
          group_by(
            .data$location_id,
            .data$sublocation_id,
            .data$sublocation_i
          ) |>
          summarise(
            across(
              c("winter_l", "winter_q", "winter_c"),
              ~ !any(is.na(.x))
            ),
            .groups = "drop"
          ),
        by = "sublocation_i"
      ) |>
      transmute(
        .data$location_id,
        .data$sublocation_id,
        .data$sim,
        .data$q0,
        q1 = ifelse(.data$winter_l, .data$q1, 0),
        q2 = ifelse(.data$winter_q, .data$q2, 0),
        q3 = ifelse(.data$winter_c, .data$q3, 0),
        across(c("q0", "q1", "q2", "q3"), ~ replace_na(.x, 0))
      ) -> ps_sublocation
  } else {
    x@Data |>
      distinct(
        .data$location_id,
        .data$sublocation_id,
        sublocation_i = as.integer(.data$sublocation_i)
      ) -> sublocation
    post_sample |>
      filter(str_detect(.data$parameter, "^sublocation_i:")) |>
      transmute(
        sublocation_i = str_remove(.data$parameter, "sublocation_i:") |>
          as.integer(),
        .data$sim,
        q0 = .data$estimate
      ) |>
      left_join(
        post_sample |>
          filter(str_detect(.data$parameter, "^subloc_1:")) |>
          transmute(
            sublocation_i = str_remove(.data$parameter, "subloc_1:") |>
              as.integer(),
            .data$sim,
            q1 = .data$estimate
          ),
        by = c("sublocation_i", "sim")
      ) |>
      left_join(
        post_sample |>
          filter(str_detect(.data$parameter, "^subloc_2:")) |>
          transmute(
            sublocation_i = str_remove(.data$parameter, "subloc_2:") |>
              as.integer(),
            .data$sim,
            q2 = .data$estimate
          ),
        by = c("sublocation_i", "sim")
      ) |>
      left_join(
        post_sample |>
          filter(str_detect(.data$parameter, "^subloc_3:")) |>
          transmute(
            sublocation_i = str_remove(.data$parameter, "subloc_3:") |>
              as.integer(),
            .data$sim,
            q3 = .data$estimate
          ),
        by = c("sublocation_i", "sim")
      ) |>
      inner_join(sublocation, by = "sublocation_i") -> ps_sublocation
    x@Data |>
      distinct(
        .data$location_id,
        .data$sublocation_id,
        .data$winter,
        winter_l = .data$loc_s_1,
        winter_q = .data$loc_s_2,
        winter_c = .data$loc_s_3
      ) -> lw
  }
  x@AnalysisMetadata |>
    select(
      species = "species_group_id",
      "model_type",
      analysis = "file_fingerprint",
      fingerprint = "status_fingerprint"
    ) |>
    bind_cols(
      ps_sublocation |>
        group_by(.data$location_id, .data$sublocation_id) |>
        summarise(
          across(
            c("q0", "q1", "q2", "q3"),
            list(
              mean = ~ mean(.x, na.rm = TRUE),
              sd = ~ sd(.x, na.rm = TRUE)
            )
          ),
          .groups = "drop"
        )
    ) |>
    write_vc(
      file = file.path("model_check", "hibernation", "sublocation_effect"),
      root = root,
      append = TRUE,
      sorting = c(
        "model_type",
        "species",
        "location_id",
        "sublocation_id",
        "analysis"
      ),
      digits = 4
    )
  update_metadata(
    file = file.path("model_check", "hibernation", "sublocation_effect"),
    root = root,
    name = "hibernation_model_check_sublocation",
    title = "Sublocation effect",
    description = paste(
      "Posterior distribution of the sublocation effect at the sublocation",
      "level."
    ),
    field_description = c(
      species = "The code of the species group",
      model_type = paste(
        "A short description of the model used to analyse the data"
      ),
      analysis = "The file fingerprint of the analysis",
      fingerprint = "The status fingerprint of the analysis",
      location_id = "The unique identifier of the location",
      sublocation_id = paste(
        "The unique identifier of the sublocation within a location"
      ),
      q0_mean = paste(
        "The mean of the intercept of the sublocation effect on the logit",
        "scale for binomial models and the log scale for nbinomial models"
      ),
      q1_mean = paste(
        "The mean of the linear effect of the sublocation effect on the logit",
        "scale for binomial models and the log scale for nbinomial models"
      ),
      q2_mean = paste(
        "The mean of the quadratic effect of the sublocation effect on the",
        "logit scale for binomial models and the log scale for nbinomial models"
      ),
      q3_mean = paste(
        "The mean of the cubic effect of the sublocation effect on the logit",
        "scale for binomial models and the log scale for nbinomial models"
      ),
      q0_sd = paste(
        "The standard error of the intercept of the sublocation effect on the",
        "logit scale for binomial models and the log scale for nbinomial models"
      ),
      q1_sd = paste(
        "The standard error of the linear effect of the sublocation effect on",
        "the logit scale for binomial models and the log scale for nbinomial",
        "models"
      ),
      q2_sd = paste(
        "The standard error of the quadratic effect of the sublocation effect",
        "on the logit scale for binomial models and the log scale for",
        "nbinomial models"
      ),
      q3_sd = paste(
        "The standard error of the cubic effect of the sublocation effect on",
        "the logit scale for binomial models and the log scale for nbinomial",
        "models"
      )
    )
  )

  # predictions at the sublocation level
  x@AnalysisMetadata |>
    select(
      species = "species_group_id",
      "model_type",
      analysis = "file_fingerprint",
      fingerprint = "status_fingerprint"
    ) |>
    bind_cols(
      ps_sublocation |>
        inner_join(
          lw,
          by = c("location_id", "sublocation_id"),
          relationship = "many-to-many"
        ) |>
        inner_join(
          ps_location |>
            select("location_id", "winter", "sim", "absolute"),
          by = c("location_id", "winter", "sim")
        ) |>
        transmute(
          .data$location_id,
          .data$sublocation_id,
          .data$sim,
          .data$winter,
          relative = .data$q0 +
            .data$q1 * .data$winter_l +
            .data$q2 * .data$winter_q +
            .data$q3 * .data$winter_c,
          absolute = .data$absolute + .data$relative
        ) |>
        group_by(.data$location_id, .data$sublocation_id, .data$winter) |>
        summarise(
          across(
            c("relative", "absolute"),
            list(
              mean = ~ mean(.x, na.rm = TRUE),
              sd = ~ sd(.x, na.rm = TRUE)
            )
          ),
          .groups = "drop"
        )
    ) |>
    write_vc(
      file = file.path("model_check", "hibernation", "sublocation"),
      root = root,
      append = TRUE,
      sorting = c(
        "model_type",
        "species",
        "location_id",
        "sublocation_id",
        "winter",
        "analysis"
      ),
      digits = 4
    )
  update_metadata(
    file = file.path("model_check", "hibernation", "sublocation"),
    root = root,
    name = "hibernation_model_check_sublocation_prediction",
    title = "Predictions at the sublocation level",
    description = paste(
      "Predictions of the model at the sublocation level.",
      "The predictions are at the logit scale for binomial models and the log",
      "scale for nbinomial models."
    ),
    field_description = c(
      species = "The code of the species group",
      model_type = paste(
        "A short description of the model used to analyse the data"
      ),
      analysis = "The file fingerprint of the analysis",
      fingerprint = "The status fingerprint of the analysis",
      location_id = "The unique identifier of the location",
      sublocation_id = paste(
        "The unique identifier of the sublocation within a location"
      ),
      winter = "The winter season defined by the year in which January falls",
      relative_mean = paste(
        "The mean of the difference between the local trend and the global",
        "trend"
      ),
      absolute_mean = paste(
        "The mean of the local trend"
      ),
      relative_sd = paste(
        "The standard error of the difference between the local trend and the",
        "global trend"
      ),
      absolute_sd = paste(
        "The standard error of the local trend"
      )
    )
  )

  # hyperparameters
  x@AnalysisMetadata |>
    select(
      species = "species_group_id",
      "model_type",
      analysis = "file_fingerprint",
      fingerprint = "status_fingerprint"
    ) |>
    bind_cols(
      names(x@Model$marginals.hyperpar) |>
        lapply(
          x = x,
          function(y, x) {
            if (!grepl("^Precision for", y)) {
              inla.zmarginal(x@Model$marginals.hyperpar[[y]], silent = TRUE) |>
                as.data.frame() |>
                transmute(
                  parameter = y,
                  .data$mean,
                  lcl = .data$quant0.025,
                  ucl = .data$quant0.975
                ) -> z
              return(z)
            }
            x@Model$marginals.hyperpar[[y]] |>
              prec2sd() |>
              transmute(
                parameter = gsub("Precision", "Stdev", y),
                .data$mean,
                lcl = .data$quant0.025,
                ucl = .data$quant0.975
              )
          }
        ) |>
        bind_rows()
    ) |>
    write_vc(
      file = file.path("model_check", "hibernation", "hyperparameters"),
      root = root,
      append = TRUE,
      sorting = c("model_type", "species", "parameter", "analysis"),
      digits = 4
    )
  update_metadata(
    file = file.path("model_check", "hibernation", "hyperparameters"),
    root = root,
    name = "hibernation_model_check_hyperparameters",
    title = "Model hyperparameters",
    description = paste(
      "Posterior distribution of the hyperparameters of the model.",
      "The model is a binomial spatio-temporal model with a logit link."
    ),
    field_description = c(
      species = "The code of the species group",
      model_type = paste(
        "A short description of the model used to analyse the data"
      ),
      analysis = "The file fingerprint of the analysis",
      fingerprint = "The status fingerprint of the analysis",
      parameter = "The name of the hyperparameter",
      mean = "The mean of the posterior distribution of the hyperparameter",
      lcl = paste(
        "The lower 95% credible limit of the posterior distribution of the",
        "hyperparameter"
      ),
      ucl = paste(
        "The upper 95% credible limit of the posterior distribution of the",
        "hyperparameter"
      )
    )
  )

  # prediction for the mesh
  mesh <- spde2mesh(x@Spde)
  as.data.frame(mesh$loc) |>
    st_as_sf(coords = c("V1", "V2")) |>
    st_union() |>
    st_convex_hull() -> hull
  hull |>
    st_sample(type = "hexagonal", size = floor(st_area(hull) / 25)) |>
    st_coordinates() |>
    as.data.frame() |>
    mutate(field_id = row_number()) -> sample_field
  projector <- fm_evaluator(
    mesh = mesh,
    loc = as.matrix(sample_field[, c("X", "Y")])
  )
  x@AnalysisMetadata |>
    select(
      species = "species_group_id",
      "model_type",
      analysis = "file_fingerprint",
      fingerprint = "status_fingerprint"
    ) |>
    bind_cols(
      as.matrix(projector$proj$A %*% ps_matern) |>
        as.data.frame() |>
        mutate(field_id = sample_field$field_id) |>
        pivot_longer(
          starts_with("sim"),
          names_to = "sim",
          values_to = "mesh"
        ) |>
        group_by(.data$field_id) |>
        summarise(
          across(
            "mesh",
            list(
              mean = ~ mean(.x, na.rm = TRUE) |>
                round(4),
              sd = ~ sd(.x, na.rm = TRUE) |>
                round(4)
            )
          ),
          .groups = "drop"
        ) |>
        inner_join(sample_field, by = "field_id")
    ) |>
    write_vc(
      file = file.path("model_check", "hibernation", "mesh_prediction"),
      root = root,
      append = TRUE,
      sorting = c("model_type", "species", "X", "Y", "analysis"),
      digits = c(X = 6, Y = 6, mesh_mean = 4, mesh_sd = 4)
    )
  update_metadata(
    file = file.path("model_check", "hibernation", "mesh_prediction"),
    root = root,
    name = "hibernation_model_check_mesh_prediction",
    title = "Prediction for the mesh",
    description = paste(
      "Prediction of the spatial random effect at the mesh nodes.",
      "The predictions are at the logit scale for binomial models and the log",
      "scale for nbinomial models."
    ),
    field_description = c(
      species = "The code of the species group",
      model_type = paste(
        "A short description of the model used to analyse the data"
      ),
      analysis = "The file fingerprint of the analysis",
      fingerprint = "The status fingerprint of the analysis",
      field_id = "The unique identifier of the mesh node",
      X = paste(
        "The X coordinate of the prediction on the mesh in kilometres using",
        "EPSG:31370"
      ),
      Y = paste(
        "The Y coordinate of the prediction on the mesh in kilometres using",
        "EPSG:31370"
      ),
      mesh_mean = paste(
        "The mean of the spatial random effect at the mesh node on the logit",
        "scale for binomial models and the log scale for nbinomial models"
      ),
      mesh_sd = paste(
        "The standard error of the spatial random effect at the mesh node on",
        "the logit scale for binomial models and the log scale for nbinomial",
        "models"
      )
    )
  )

  rm(
    hull,
    location_type,
    location_winter,
    loc_coordinates,
    mesh,
    post_sample,
    projector,
    ps_cwinter,
    ps_intercept,
    ps_location,
    ps_matern,
    ps_sublocation,
    sample_field,
    sublocation_winter,
    x
  )
  gc(verbose = FALSE)
  return(invisible(NULL))
}
