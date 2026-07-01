#######################################################
## Postprocess spatial datasets                      ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
initialize_submap_transect_set_data <- function(data_dir, out_dir, id, quiet = FALSE) {
  with_error_handling(
    {
      if (is.character(id) && all(id %in% sf_hull$submap_transect_set)) {
        hull <- filter(sf_hull, submap_transect_set %in% id)
        submap_ids <- unlist(strsplit(hull$submap_transect_ids, ",")) |> trimws()
        data <- filter(sf_submap, submap_transect_id %in% submap_ids)
      } else {
        stop("Unrecognized input for id! Use the submap_transect_set")
      }

      if (nrow(data) == 0) stop("Submap transect id: ", id, " not found in sf_submap")

      if (!quiet) cat(" .. Compiling data: '", id, "'\n", sep = "")

      df <- rowwise(data) |>
        mutate(
          hull = st_geometry(hull),
          grid = list(crop_to_bbox(sf_grid, hull)),
          ihfc2024_obs_hull = list(crop_to_bbox(sf_ihfc, hull)),
          lucazeau2019_sim_hull = list(crop_to_bbox(sf_sim, hull))
        ) |>
        ungroup()

      if (!quiet) cat(" .. Projecting HF data to transect strips: '", id, "'\n", sep = "")
      add_transect_strip_projections(df, quiet = quiet)
    },
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
compile_submap_transect_set_data <- function(cache_path,
                                             data_dir,
                                             out_dir,
                                             id,
                                             variogram_model = NULL,
                                             fitted_variogram = NULL,
                                             max_n_point_pairs = NULL,
                                             max_point_pair_distance = NULL,
                                             quiet = FALSE) {
  with_error_handling(
    {
      if (file.exists(cache_path)) {
        if (!quiet) cat(" .. Reading cache:", cache_path, "\n")
        load(cache_path)

        if (!exists("df") || is.null(df)) {
          df <- initialize_submap_transect_set_data(data_dir, out_dir, id, quiet = quiet)
          save(df, file = cache_path)
          return(df)
        }

        strip_sentinel_cols <- c("ihfc2024_obs_strip1_proj", "lucazeau2019_sim_strip1_proj")
        if (!all(strip_sentinel_cols %in% names(df))) {
          if (!quiet) cat(" .. Adding strip projections\n")
          df <- add_transect_strip_projections(df, quiet = quiet)
          save(df, file = cache_path)
        }

        if (!is.null(variogram_model)) {
          prefix <- "kerswell2025_krg_"
          model_str <- str_to_lower(variogram_model)
          krg_hull_str <- paste0(prefix, model_str, "_hull")

          if (krg_hull_str %in% names(df)) {
            return(df)
          }

          if (!is.null(fitted_variogram) && !is.null(max_n_point_pairs) && !is.null(max_point_pair_distance)) {
            if (!quiet) cat(" .. Adding krige results for variogram model '", model_str, "'\n", sep = "")

            shp_krg <- krige_submap_transect_set(
              sf_hf = df$ihfc2024_obs_hull[[1]],
              data_col = "ihfc2024_obs",
              fitted_variogram = fitted_variogram,
              sf_grid = df$grid[[1]],
              max_n_point_pairs = max_n_point_pairs,
              max_point_pair_distance = max_point_pair_distance,
              quiet = quiet
            )

            df <- df |>
              rowwise() |>
              mutate(
                !!krg_hull_str := list(crop_to_bbox(shp_krg, hull)),
                krg_sim_diff_hull = list(process_sim_krg_diff(!!sym(krg_hull_str), lucazeau2019_sim_hull))
              ) |>
              ungroup()

            df <- add_krige_strip_projections(df, krg_hull_str, quiet = quiet)

            if (!quiet) cat(" .. Updating cache with krige results: ", cache_path, "\n", sep = "")
            save(df, file = cache_path)

            return(df)
          } else {
            return(df)
          }
        } else {
          return(df)
        }
      }

      df <- initialize_submap_transect_set_data(data_dir, out_dir, id, quiet = quiet)
      if (!quiet) cat(" .. Caching data:", cache_path, "\n")
      create_path_dir(cache_path)
      save(df, file = cache_path)
      df
    },
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
project_hf_strips_to_transect <- function(transect_geometry, shp_hf, strip_widths = c(7.5e4), quiet = FALSE) {
  empty <- vector("list", length(strip_widths))

  with_error_handling(
    {
      if (is.null(shp_hf) || nrow(shp_hf) == 0) {
        return(empty)
      }

      if (!inherits(transect_geometry, "sfc")) {
        transect_geometry <- st_sfc(transect_geometry, crs = st_crs(shp_hf))
      } else if (st_crs(transect_geometry) != st_crs(shp_hf)) {
        shp_hf <- st_transform(shp_hf, st_crs(transect_geometry))
      }

      proj <- determine_projection(sf_obj = transect_geometry, proj_type = "aeqd", quiet = quiet)
      transect_proj <- st_transform(transect_geometry, proj$wkt)
      src_crs <- st_crs(shp_hf)

      map(strip_widths, ~ {
        strip_buffer <- st_buffer(transect_proj, .x, endCapStyle = "FLAT") |> st_transform(src_crs)

        shp_strip <- crop_to_bbox(shp_hf, strip_buffer)
        if (is.null(shp_strip) || nrow(shp_strip) == 0) {
          if (!quiet) cat(" !! Warning: no heat flow data to project onto transect\n", sep = "")
          return(NULL)
        }

        shp_strip_proj <- st_transform(shp_strip, proj$wkt)
        if (is.null(shp_strip_proj) || nrow(shp_strip_proj) == 0) {
          if (!quiet) cat(" !! Warning: no (projected) heat flow data to project onto transect\n", sep = "")
          return(NULL)
        }

        project_hf_to_transect(transect_geometry = transect_proj, shp_hf = shp_strip_proj, quiet = quiet)
      })
    },
    default = empty,
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
add_transect_strip_projections <- function(df, strip_widths = c(7.5e4), quiet = FALSE) {
  with_error_handling(
    {
      df <- df |>
        rowwise() |>
        mutate(
          .ihfc_strips = list(project_hf_strips_to_transect(submap_transect_geometry, ihfc2024_obs_hull, strip_widths, quiet = quiet)),
          .sim_strips = list(project_hf_strips_to_transect(submap_transect_geometry, lucazeau2019_sim_hull, strip_widths, quiet = quiet))
        ) |>
        ungroup()

      df <- reduce(seq_along(strip_widths), function(acc, i) {
        ihfc_proj <- paste0("ihfc2024_obs_strip", i, "_proj")
        sim_proj <- paste0("lucazeau2019_sim_strip", i, "_proj")
        ihfc_trend <- paste0("ihfc2024_obs_strip", i, "_trend")
        sim_trend <- paste0("lucazeau2019_sim_strip", i, "_trend")

        acc |>
          rowwise() |>
          mutate(
            !!ihfc_proj := list(.ihfc_strips[[i]]),
            !!sim_proj := list(.sim_strips[[i]]),
            !!ihfc_trend := list(fit_trends_to_projected_observations(.ihfc_strips[[i]], quiet = quiet)),
            !!sim_trend := list(fit_trends_to_projected_observations(.sim_strips[[i]], quiet = quiet))
          ) |>
          ungroup()
      }, .init = df)

      select(df, -.ihfc_strips, -.sim_strips)
    },
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
add_krige_strip_projections <- function(df, krg_hull_str, strip_widths = c(7.5e4), quiet = FALSE) {
  with_error_handling(
    {
      prefix <- "kerswell2025_krg_"
      model_str <- str_remove(str_remove(krg_hull_str, paste0("^", prefix)), "_hull$")

      df <- df |>
        rowwise() |>
        mutate(.krg_strips = list(project_hf_strips_to_transect(submap_transect_geometry, !!sym(krg_hull_str), strip_widths, quiet = quiet))) |>
        ungroup()

      df <- reduce(seq_along(strip_widths), function(acc, i) {
        krg_proj <- paste0(prefix, model_str, "_strip", i, "_proj")
        krg_trend <- paste0(prefix, model_str, "_strip", i, "_trend")

        acc |>
          rowwise() |>
          mutate(
            !!krg_proj := list(.krg_strips[[i]]),
            !!krg_trend := list(fit_trends_to_projected_observations(.krg_strips[[i]], quiet = quiet))
          ) |>
          ungroup()
      }, .init = df)

      select(df, -.krg_strips)
    },
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
project_hf_to_transect <- function(transect_geometry, shp_hf, quiet = FALSE) {
  empty <- st_sf(distance = numeric(0), hf = numeric(0), sigma = numeric(0), geometry = st_sfc())

  with_error_handling(
    {
      valid_names <- c("ihfc2024_obs", "lucazeau2019_sim_est", "kerswell2025_krg_est")
      matched_col <- intersect(names(shp_hf), valid_names)
      has_recognized_col <- length(matched_col) > 0

      if (
        is.null(shp_hf) ||
          !inherits(shp_hf, "sf") ||
          !inherits(transect_geometry, "sfc") ||
          st_crs(transect_geometry) != st_crs(shp_hf) ||
          !has_recognized_col
      ) {
        if (!quiet) cat(" !! Warning: invalid heat flow data for projection onto transect\n", sep = "")
        return(empty)
      }

      target_hf_col <- matched_col[1]
      projected_hf <- shp_hf[[target_hf_col]]

      sigma_map <- c(
        "ihfc2024_obs" = NA,
        "lucazeau2019_sim_est" = "lucazeau2019_sim_sigma",
        "kerswell2025_krg_est" = "kerswell2025_krg_sigma"
      )
      target_sigma_col <- sigma_map[target_hf_col]

      if (!is.na(target_sigma_col) && target_sigma_col %in% names(shp_hf)) {
        projected_sigma <- shp_hf[[target_sigma_col]]
      } else {
        projected_sigma <- rep(NA_real_, nrow(shp_hf))
      }

      suppressWarnings({
        suppressMessages({
          projected_distances <- st_line_project(transect_geometry, st_geometry(shp_hf), normalized = TRUE)

          st_as_sf(st_line_interpolate(transect_geometry, projected_distances)) |>
            rename(geometry = x) |>
            mutate(distance = projected_distances, hf = projected_hf, sigma = projected_sigma, .before = geometry)
        })
      })
    },
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit_trends_to_projected_observations <- function(shp, n = 1e3, minimum_observations = 10, lambda = 0.1, k = -1, quiet = FALSE) {
  empty <- tibble(distance = numeric(), hf = numeric(), method = character())

  with_error_handling(
    {
      if (is.null(shp) || nrow(shp) == 0) {
        if (!quiet) cat(" !! Warning: cannot fit trends because input data is NULL or empty\n")
        return(empty)
      }

      df <- shp |> st_set_geometry(NULL)

      if (nrow(df) < minimum_observations) {
        if (!quiet) cat(" !! Warning: cannot fit models with less than ", minimum_observations, " projected obs!\n", sep = "")
        return(empty)
      }

      df <- df |> arrange(distance)

      original_range <- range(df$distance, na.rm = TRUE)
      new_distances <- seq(original_range[1], original_range[2], length.out = n)
      pred_df <- data.frame(distance = new_distances)

      gam_preds <- tryCatch(
        {
          gam_model <- gam(hf ~ s(distance, k = k), data = df, method = "REML")
          predict(gam_model, newdata = pred_df, type = "response")
        },
        error = function(e) {
          if (!quiet) cat(" !! Warning: GAM fitting failed: ", conditionMessage(e), "\n", sep = "")
          rep(NA_real_, n)
        }
      )

      rqss_preds <- tryCatch(
        {
          rqss_model <- rqss(hf ~ qss(distance, lambda = lambda), data = df, tau = 0.5)
          as.vector(predict(rqss_model, newdata = pred_df))
        },
        error = function(e) {
          if (!quiet) cat(" !! Warning: RQSS fitting failed: ", conditionMessage(e), "\n", sep = "")
          rep(NA_real_, n)
        }
      )

      bind_rows(
        tibble(distance = new_distances, hf = as.numeric(gam_preds), method = "mean"),
        tibble(distance = new_distances, hf = as.numeric(rqss_preds), method = "median")
      )
    },
    default = empty,
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
process_sim_krg_diff <- function(shp_krg, shp_sim) {
  with_error_handling({
    if (!any(names(shp_krg) %in% c("kerswell2025_krg_est", "distance"))) stop("Unrecognized heatflow data column")
    if (!any(names(shp_sim) %in% c("lucazeau2019_sim_est", "distance"))) stop("Unrecognized heatflow data column")

    if (any(names(shp_sim) == "distance")) {
      shp_sim <-
        shp_sim |>
        rename(lucazeau2019_sim_est = hf, lucazeau2019_sim_sigma = sigma) |>
        select(lucazeau2019_sim_est, lucazeau2019_sim_sigma, geometry)
      shp_krg <-
        shp_krg |>
        rename(kerswell2025_krg_est = hf, kerswell2025_krg_sigma = sigma) |>
        select(kerswell2025_krg_est, kerswell2025_krg_sigma, geometry)
      shp_sim |>
        mutate(kerswell2025_krg_est = shp_krg$kerswell2025_krg_est, kerswell2025_krg_sigma = shp_krg$kerswell2025_krg_sigma) |>
        mutate(
          sim_krg_diff_est = lucazeau2019_sim_est - kerswell2025_krg_est,
          sim_krg_diff_sigma = lucazeau2019_sim_sigma - kerswell2025_krg_sigma,
          .before = geometry
        ) |>
        rename(projected_sim_krg_diff_geometry = geometry) |>
        select(sim_krg_diff_est, sim_krg_diff_sigma, projected_sim_krg_diff_geometry)
    } else {
      shp_sim |>
        mutate(kerswell2025_krg_est = shp_krg$kerswell2025_krg_est, kerswell2025_krg_sigma = shp_krg$kerswell2025_krg_sigma) |>
        mutate(
          sim_krg_diff_est = lucazeau2019_sim_est - kerswell2025_krg_est,
          sim_krg_diff_sigma = lucazeau2019_sim_sigma - kerswell2025_krg_sigma,
          .before = lucazeau2019_geometry
        ) |>
        rename(sim_krg_diff_geometry = lucazeau2019_geometry) |>
        select(sim_krg_diff_est, sim_krg_diff_sigma, sim_krg_diff_geometry)
    }
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_nlopt_iter_paths <- function(nlopt_iter_dir, nlopt_krige_dir, submap_transect_set, quiet = FALSE) {
  with_error_handling(
    {
      nlopt_iter_paths_run <- list.files(nlopt_iter_dir, pattern = paste0("^", submap_transect_set, ".*-run$"), full.names = TRUE)
      krige_paths <- list.files(nlopt_krige_dir, pattern = paste0("^", submap_transect_set, ".*\\.RData$"), full.names = TRUE)

      if (length(nlopt_iter_paths_run) < 1) {
        if (!quiet) cat(" !! Warning: no nlopt iteration runs found for: ", submap_transect_set, "\n", sep = "")
        return(list())
      }
      if (length(krige_paths) < 1) {
        if (!quiet) cat(" !! Warning: no krige model files found for: ", submap_transect_set, "\n", sep = "")
        return(list())
      }

      iter_basenames <- str_remove(basename(nlopt_iter_paths_run), "-run$")
      krige_basenames <- str_remove(basename(krige_paths), "\\.RData$")

      valid <- iter_basenames %in% krige_basenames
      nlopt_iter_paths_run[valid]
    },
    default = list(),
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
read_nlopt_iteration <- function(in_path) {
  empty <- tibble(
    submap_transect_set = character(),
    itr = integer(),
    n_heatflow_obs = numeric(),
    variogram_cutoff = numeric(),
    n_variogram_lags = numeric(),
    max_n_point_pairs = numeric(),
    max_point_pair_distance = numeric(),
    variogram_model = character(),
    variogram_weight = numeric(),
    variogram_rmse = numeric(),
    variogram_sigma = numeric(),
    variogram_cost = numeric(),
    interpolation_weight = numeric(),
    interpolation_rmse = numeric(),
    interpolation_sigma = numeric(),
    interpolation_cost = numeric(),
    total_cost = numeric()
  )

  if (!file.exists(in_path)) {
    return(empty)
  }

  read_file <- function(path) {
    con <- file(path)
    on.exit(close(con))
    suppressWarnings(readLines(con))
  }

  extract_value <- function(lines, pattern, as_numeric = TRUE) {
    matched_lines <- lines[grepl(pattern, lines)]

    if (length(matched_lines) == 0) {
      return(if (as_numeric) rep(NA_real_, 1) else NA_character_)
    }

    clean_values <- gsub(pattern, "", matched_lines)
    if (as_numeric) as.numeric(clean_values) else clean_values
  }

  with_error_handling(
    {
      data <- read_file(in_path)

      submap_transect_set_ext <- extract_value(data, "^Submap transect set: ", FALSE)
      n_heatflow_obs_ext <- extract_value(data, "^No. heatflow observations: ")
      variogram_cutoff_ext <- extract_value(data, "^Variogram variogram_cutoff: ")
      n_variogram_lags_ext <- extract_value(data, "^No. variogram lags: ")
      max_n_point_pairs_ext <- extract_value(data, "^Max no. point pairs: ")
      max_point_pair_distance_ext <- extract_value(data, "^Max point pair distance: ")
      variogram_model_ext <- extract_value(data, "^Variogram model: ", FALSE)
      variogram_weight_ext <- extract_value(data, "^Variogram weight: ")
      variogram_rmse_ext <- extract_value(data, "^Variogram rmse: ")
      variogram_sigma_ext <- extract_value(data, "^Variogram uncertainty: ")
      variogram_cost_ext <- extract_value(data, "^Variogram cost: ")
      interpolation_weight_ext <- extract_value(data, "^Interpolation weight: ")
      interpolation_rmse_ext <- extract_value(data, "^Interpolation rmse: ")
      interpolation_sigma_ext <- extract_value(data, "^Interpolation sd: ")
      interpolation_cost_ext <- extract_value(data, "^Interpolation cost: ")
      total_cost_ext <- extract_value(data, "^Total cost: ")

      tibble(
        submap_transect_set = submap_transect_set_ext,
        n_heatflow_obs = n_heatflow_obs_ext,
        variogram_cutoff = variogram_cutoff_ext,
        n_variogram_lags = n_variogram_lags_ext,
        max_n_point_pairs = max_n_point_pairs_ext,
        max_point_pair_distance = max_point_pair_distance_ext,
        variogram_model = variogram_model_ext,
        variogram_weight = variogram_weight_ext,
        variogram_rmse = variogram_rmse_ext,
        variogram_sigma = variogram_sigma_ext,
        variogram_cost = variogram_cost_ext,
        interpolation_weight = interpolation_weight_ext,
        interpolation_rmse = interpolation_rmse_ext,
        interpolation_sigma = interpolation_sigma_ext,
        interpolation_cost = interpolation_cost_ext,
        total_cost = total_cost_ext
      ) |>
        group_by(submap_transect_set, variogram_model) |>
        mutate(itr = row_number(), .after = submap_transect_set) |>
        ungroup()
    },
    default = empty
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_optimal_krige_model <- function(out_dir, submap_transect_set, quiet = FALSE) {
  empty <- list(
    "nlopt_iter_df" = NULL,
    "optimal_krige_model" = NULL,
    "experimental_variogram" = NULL,
    "fitted_variogram" = NULL
  )

  with_error_handling(
    {
      nlopt_dir <- file.path(out_dir, "nlopt")
      nlopt_iter_dir <- file.path(nlopt_dir, "iterations")
      nlopt_krige_dir <- file.path(nlopt_dir, "krige_models")

      nlopt_iter_paths <- get_nlopt_iter_paths(nlopt_iter_dir, nlopt_krige_dir, submap_transect_set, quiet = quiet)

      if (length(nlopt_iter_paths) < 1) {
        if (!quiet) cat(" !! Warning: no valid nlopt iter paths for: ", submap_transect_set, "\n", sep = "")
        return(list("nlopt_iter_df" = NULL, "optimal_krige_model" = NULL, "experimental_variogram" = NULL, "fitted_variogram" = NULL))
      }

      nlopt_iter_df <- map_df(nlopt_iter_paths, read_nlopt_iteration)

      if (nrow(nlopt_iter_df) < 1) {
        if (!quiet) cat(" !! Warning: empty nlopt iteration table for: ", submap_transect_set, "\n", sep = "")
        return(list("nlopt_iter_df" = NULL, "optimal_krige_model" = NULL, "experimental_variogram" = NULL, "fitted_variogram" = NULL))
      }

      optimal_krige_model <- slice_min(nlopt_iter_df, total_cost)

      if (nrow(optimal_krige_model) > 1) {
        optimal_krige_model <- slice(optimal_krige_model, nrow(optimal_krige_model))
      }

      nlopt_path <- file.path(nlopt_krige_dir, paste0(submap_transect_set, "-", optimal_krige_model$variogram_model, ".RData"))

      if (!file.exists(nlopt_path)) {
        if (!quiet) cat(" !! Warning: no nlopt RData found for: ", submap_transect_set, "\n", sep = "")
        return(list(
          "nlopt_iter_df" = nlopt_iter_df,
          "optimal_krige_model" = optimal_krige_model,
          "experimental_variogram" = NULL,
          "fitted_variogram" = NULL
        ))
      }

      load(nlopt_path)
      nlopt_out <- get(paste0(submap_transect_set, "_", optimal_krige_model$variogram_model))

      list(
        "nlopt_iter_df" = nlopt_iter_df,
        "optimal_krige_model" = optimal_krige_model,
        "experimental_variogram" = nlopt_out$experimental_variogram,
        "fitted_variogram" = nlopt_out$fitted_variogram
      )
    },
    default = empty,
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summarize_optimal_krige_models <- function(out_dir, submap_transect_sets, seed = 42, nprocs = NULL, parallel = FALSE, quiet = FALSE) {
  if (is.null(nprocs)) nprocs <- availableCores() - 2

  f <- function(submap_transect_set) {
    with_error_handling(get_optimal_krige_model(out_dir, submap_transect_set, quiet = quiet)$optimal_krige_model, default = NULL, quiet = quiet)
  }

  if (!quiet) cat(" .. Summarizing optimal krige models\n")
  with_error_handling(
    {
      if (parallel) {
        set.seed(seed)
        plan(multisession, workers = nprocs)
        as_tibble(future_map_dfr(submap_transect_sets, ~ f(.x), .options = furrr_options(seed = seed)))
      } else {
        map_df(submap_transect_sets, ~ f(.x))
      }
    },
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_closest_interp_obs_point_pairs <- function(data_dir,
                                               out_dir,
                                               submap_transect_set,
                                               threshold = 1e4,
                                               variogram_model = NULL,
                                               fitted_variogram = NULL,
                                               max_n_point_pairs = NULL,
                                               max_point_pair_distance = NULL,
                                               quiet = FALSE) {
  with_error_handling(
    {
      prefix <- "kerswell2025_krg_"
      model_str <- str_to_lower(variogram_model)
      krg_hull_str <- paste0(prefix, model_str, "_hull")

      cache_path <- file.path(out_dir, "submap", paste0(submap_transect_set, ".RData"))
      if (!is.null(variogram_model) && !is.null(fitted_variogram) && !is.null(max_n_point_pairs) && !is.null(max_point_pair_distance)) {
        data <- compile_submap_transect_set_data(
          cache_path = cache_path,
          data_dir = data_dir,
          out_dir = out_dir,
          id = submap_transect_set,
          variogram_model = variogram_model,
          fitted_variogram = fitted_variogram,
          max_n_point_pairs = max_n_point_pairs,
          max_point_pair_distance = max_point_pair_distance,
          quiet = quiet
        )
        shp_hf_interp <- data[[krg_hull_str]][[1]]
      } else {
        data <- compile_submap_transect_set_data(
          cache_path = cache_path,
          data_dir = data_dir,
          out_dir = out_dir,
          id = submap_transect_set,
          quiet = quiet
        )
        shp_hf_interp <- data$lucazeau2019_sim_hull[[1]]
      }

      shp_grid <- data$grid[[1]]
      shp_hf_obs <- data$ihfc2024_obs_hull[[1]]

      nearest_obs <- st_nearest_feature(shp_grid, shp_hf_obs)
      nearest_interp <- st_nearest_feature(shp_hf_obs, shp_grid)

      obs_within_threshold_distance <- st_distance(shp_grid, shp_hf_obs[nearest_obs, ], by_element = TRUE) < set_units(threshold, "m")
      interp_within_threshold_distance <- st_distance(shp_hf_obs, shp_grid[nearest_interp, ], by_element = TRUE) < set_units(threshold, "m")

      shp_hf_obs <- shp_hf_obs[nearest_obs, ][obs_within_threshold_distance, ]
      shp_hf_interp <- shp_hf_interp[nearest_interp, ][interp_within_threshold_distance, ]

      if (any(names(shp_hf_interp) %in% c("lucazeau2019_sim_est", "lucazeau2019_geometry"))) {
        shp_hf_interp |>
          mutate(ihfc2024_obs = shp_hf_obs[st_nearest_feature(shp_hf_interp, shp_hf_obs), ]$ihfc2024_obs, .before = lucazeau2019_geometry) |>
          select(-c(lucazeau2019_sim_sigma, lucazeau2019_obs)) |>
          rename(geometry = lucazeau2019_geometry)
      } else if (any(names(shp_hf_interp) %in% c("kerswell2025_krg_est", "kerswell2025_geometry"))) {
        shp_hf_interp |>
          mutate(ihfc2024_obs = shp_hf_obs[st_nearest_feature(shp_hf_interp, shp_hf_obs), ]$ihfc2024_obs, .before = kerswell2025_geometry) |>
          select(-c(kerswell2025_krg_sigma)) |>
          rename(geometry = kerswell2025_geometry)
      } else {
        invisible()
      }
    },
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summarize_interpolation_est_differences <- function(data_dir, out_dir, submap_transect_sets, seed = 42,
                                                    nprocs = NULL, parallel = FALSE, quiet = FALSE) {
  if (is.null(nprocs)) nprocs <- availableCores() - 2

  f <- function(submap_transect_set) {
    with_error_handling(
      {
        opt_data <- get_optimal_krige_model(out_dir, submap_transect_set, quiet = quiet)
        variogram_model <- opt_data$optimal_krige_model$variogram_model
        fitted_variogram <- opt_data$fitted_variogram
        max_n_point_pairs <- opt_data$optimal_krige_model$max_n_point_pairs
        max_point_pair_distance <- opt_data$optimal_krige_model$max_point_pair_distance

        if (is.null(variogram_model) && is.null(fitted_variogram) && is.null(max_n_point_pairs) && is.null(max_point_pair_distance)) {
          if (!quiet) cat(" !! Warning: no valid kriging model found for submap transect set '", submap_transect_set, "'\n", sep = "")
          return(tibble(
            submap_transect_set = submap_transect_set,
            n_grid = NA_integer_,
            n = NA_integer_,
            rmse = NA_real_,
            min = NA_real_,
            max = NA_real_,
            median = NA_real_,
            iqr = NA_real_,
            mean = NA_real_,
            sigma = NA_real_
          ))
        }

        cache_path <- file.path(out_dir, "submap", paste0(submap_transect_set, ".RData"))
        data <- compile_submap_transect_set_data(
          cache_path = cache_path,
          data_dir = data_dir,
          out_dir = out_dir,
          id = submap_transect_set,
          variogram_model = variogram_model,
          fitted_variogram = fitted_variogram,
          max_n_point_pairs = max_n_point_pairs,
          max_point_pair_distance = max_point_pair_distance,
          quiet = quiet
        )

        data$krg_sim_diff_hull[[1]] |>
          st_set_geometry(NULL) |>
          summarize(
            submap_transect_set = submap_transect_set,
            n_grid = n(),
            n = sum(!is.na(sim_krg_diff_est)),
            rmse = sqrt(mean(sim_krg_diff_est^2, na.rm = TRUE)),
            min = min(sim_krg_diff_est, na.rm = TRUE),
            max = max(sim_krg_diff_est, na.rm = TRUE),
            median = median(sim_krg_diff_est, na.rm = TRUE),
            iqr = IQR(sim_krg_diff_est, na.rm = TRUE),
            mean = mean(sim_krg_diff_est, na.rm = TRUE),
            sigma = sd(sim_krg_diff_est, na.rm = TRUE)
          )
      },
      default = NULL,
      quiet = quiet
    )
  }

  if (!quiet) cat(" .. Summarizing interpolation differences\n")
  with_error_handling(
    {
      if (parallel) {
        set.seed(seed)
        plan(multisession, workers = nprocs)
        as_tibble(future_map_dfr(submap_transect_sets, ~ f(.x), .options = furrr_options(seed = seed)))
      } else {
        map_df(submap_transect_sets, ~ f(.x))
      }
    },
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summarize_interpolation_sigma_differences <- function(data_dir, out_dir, submap_transect_sets, seed = 42,
                                                      nprocs = NULL, parallel = FALSE, quiet = FALSE) {
  if (is.null(nprocs)) nprocs <- availableCores() - 2

  f <- function(submap_transect_set) {
    with_error_handling(
      {
        opt_data <- get_optimal_krige_model(out_dir, submap_transect_set, quiet = quiet)
        variogram_model <- opt_data$optimal_krige_model$variogram_model
        fitted_variogram <- opt_data$fitted_variogram
        max_n_point_pairs <- opt_data$optimal_krige_model$max_n_point_pairs
        max_point_pair_distance <- opt_data$optimal_krige_model$max_point_pair_distance

        if (is.null(variogram_model) && is.null(fitted_variogram) && is.null(max_n_point_pairs) && is.null(max_point_pair_distance)) {
          if (!quiet) cat(" !! Warning: no valid kriging model found for submap transect set '", submap_transect_set, "'\n", sep = "")
          return(tibble(
            submap_transect_set = submap_transect_set,
            n_grid = NA_integer_,
            n = NA_integer_,
            rmse = NA_real_,
            min = NA_real_,
            max = NA_real_,
            median = NA_real_,
            iqr = NA_real_,
            mean = NA_real_,
            sigma = NA_real_
          ))
        }

        cache_path <- file.path(out_dir, "submap", paste0(submap_transect_set, ".RData"))
        data <- compile_submap_transect_set_data(
          cache_path = cache_path,
          data_dir = data_dir,
          out_dir = out_dir,
          id = submap_transect_set,
          variogram_model = variogram_model,
          fitted_variogram = fitted_variogram,
          max_n_point_pairs = max_n_point_pairs,
          max_point_pair_distance = max_point_pair_distance,
          quiet = quiet
        )

        data$krg_sim_diff_hull[[1]] |>
          st_set_geometry(NULL) |>
          summarize(
            submap_transect_set = submap_transect_set,
            n_grid = n(),
            n = sum(!is.na(sim_krg_diff_sigma)),
            rmse = sqrt(mean(sim_krg_diff_sigma^2, na.rm = TRUE)),
            min = min(sim_krg_diff_sigma, na.rm = TRUE),
            max = max(sim_krg_diff_sigma, na.rm = TRUE),
            median = median(sim_krg_diff_sigma, na.rm = TRUE),
            iqr = IQR(sim_krg_diff_sigma, na.rm = TRUE),
            mean = mean(sim_krg_diff_sigma, na.rm = TRUE),
            sigma = sd(sim_krg_diff_sigma, na.rm = TRUE)
          )
      },
      default = NULL,
      quiet = quiet
    )
  }

  if (!quiet) cat(" .. Summarizing uncertainty differences\n")
  with_error_handling(
    {
      if (parallel) {
        set.seed(seed)
        plan(multisession, workers = nprocs)
        as_tibble(future_map_dfr(submap_transect_sets, ~ f(.x), .options = furrr_options(seed = seed)))
      } else {
        map_df(submap_transect_sets, ~ f(.x))
      }
    },
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summarize_interpolation_accuracies <- function(data_dir, out_dir, submap_transect_sets, seed = 42, nprocs = NULL, parallel = FALSE, quiet = FALSE) {
  if (is.null(nprocs)) nprocs <- availableCores() - 2

  f <- function(submap_transect_set) {
    suppressWarnings({
      suppressMessages({
        with_error_handling(
          {
            opt_data <- get_optimal_krige_model(out_dir, submap_transect_set, quiet = quiet)
            variogram_model <- opt_data$optimal_krige_model$variogram_model
            fitted_variogram <- opt_data$fitted_variogram
            max_n_point_pairs <- opt_data$optimal_krige_model$max_n_point_pairs
            max_point_pair_distance <- opt_data$optimal_krige_model$max_point_pair_distance

            if (is.null(variogram_model) && is.null(fitted_variogram) && is.null(max_n_point_pairs) && is.null(max_point_pair_distance)) {
              if (!quiet) cat(" !! Warning: no valid kriging model found for submap transect set '", submap_transect_set, "'\n", sep = "")
              return(tibble(
                submap_transect_set = submap_transect_set,
                n_similarity = NA_integer_,
                min_similarity = NA_real_,
                max_similarity = NA_real_,
                mean_similarity = NA_real_,
                sd_similarity = NA_real_,
                med_similarity = NA_real_,
                iqr_similarity = NA_real_,
                rmse_similarity = NA_real_,
                n_krige = NA_integer_,
                min_krige = NA_real_,
                max_krige = NA_real_,
                mean_krige = NA_real_,
                sd_krige = NA_real_,
                med_krige = NA_real_,
                iqr_krige = NA_real_,
                rmse_krige = NA_real_
              ))
            }

            point_pairs_sim <- get_closest_interp_obs_point_pairs(data_dir, out_dir, submap_transect_set, quiet = quiet)
            point_pairs_krg <- get_closest_interp_obs_point_pairs(
              data_dir,
              out_dir,
              submap_transect_set,
              variogram_model = variogram_model,
              fitted_variogram = fitted_variogram,
              max_n_point_pairs = max_n_point_pairs,
              max_point_pair_distance = max_point_pair_distance,
              quiet = quiet
            )

            tibble(
              submap_transect_set = submap_transect_set,
              n_similarity = nrow(point_pairs_sim[!is.na(point_pairs_sim$lucazeau2019_sim_est), ]),
              min_similarity = min(point_pairs_sim$lucazeau2019_sim_est - point_pairs_sim$ihfc2024_obs, na.rm = TRUE),
              max_similarity = max(point_pairs_sim$lucazeau2019_sim_est - point_pairs_sim$ihfc2024_obs, na.rm = TRUE),
              mean_similarity = mean(point_pairs_sim$lucazeau2019_sim_est - point_pairs_sim$ihfc2024_obs, na.rm = TRUE),
              sd_similarity = sd(point_pairs_sim$lucazeau2019_sim_est - point_pairs_sim$ihfc2024_obs, na.rm = TRUE),
              med_similarity = median(point_pairs_sim$lucazeau2019_sim_est - point_pairs_sim$ihfc2024_obs, na.rm = TRUE),
              iqr_similarity = IQR(point_pairs_sim$lucazeau2019_sim_est - point_pairs_sim$ihfc2024_obs, na.rm = TRUE),
              rmse_similarity = sqrt(mean((point_pairs_sim$lucazeau2019_sim_est - point_pairs_sim$ihfc2024_obs)^2, na.rm = TRUE)),
              n_krige = nrow(point_pairs_krg[!is.na(point_pairs_krg$kerswell2025_krg_est), ]),
              min_krige = min(point_pairs_krg$kerswell2025_krg_est - point_pairs_krg$ihfc2024_obs, na.rm = TRUE),
              max_krige = max(point_pairs_krg$kerswell2025_krg_est - point_pairs_krg$ihfc2024_obs, na.rm = TRUE),
              mean_krige = mean(point_pairs_krg$kerswell2025_krg_est - point_pairs_krg$ihfc2024_obs, na.rm = TRUE),
              sd_krige = sd(point_pairs_krg$kerswell2025_krg_est - point_pairs_krg$ihfc2024_obs, na.rm = TRUE),
              med_krige = median(point_pairs_krg$kerswell2025_krg_est - point_pairs_krg$ihfc2024_obs, na.rm = TRUE),
              iqr_krige = IQR(point_pairs_krg$kerswell2025_krg_est - point_pairs_krg$ihfc2024_obs, na.rm = TRUE),
              rmse_krige = sqrt(mean((point_pairs_krg$kerswell2025_krg_est - point_pairs_krg$ihfc2024_obs)^2, na.rm = TRUE))
            )
          },
          default = NULL,
          quiet = quiet
        )
      })
    })
  }

  if (!quiet) cat(" .. Summarizing interpolation accuracies\n")
  with_error_handling(
    {
      if (parallel) {
        set.seed(seed)
        plan(multisession, workers = nprocs)
        as_tibble(future_map_dfr(submap_transect_sets, ~ f(.x), .options = furrr_options(seed = seed)))
      } else {
        map_df(submap_transect_sets, ~ f(.x))
      }
    },
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summarize_global_ccf <- function(data_dir, out_dir, submap_transect_sets, strip_idx = 1, n_grid = 500, max_lag_frac = 0.1, quiet = FALSE) {
  empty <- tibble(
    submap_transect_set = character(),
    region = character(),
    pair_idx = numeric(),
    pair_label = character(),
    n_total = numeric(),
    obs_ccf = numeric(),
    sim_ccf = numeric(),
    krg_ccf = numeric(),
    global_pair_uid = factor()
  )

  if (!quiet) cat(" .. Summarizing global CCF profiles across submap transect sets\n")

  precompile_set_data(data_dir = data_dir, out_dir = out_dir, submap_transect_sets = submap_transect_sets, quiet = quiet)

  global_ccf_summary <- suppressWarnings(suppressMessages({
    map_dfr(submap_transect_sets, function(set_id) {
      with_error_handling(
        {
          result <- load_set_data(data_dir = data_dir, out_dir = out_dir, submap_transect_set = set_id, quiet = quiet)
          if (is.null(result)) {
            return(invisible())
          }

          data <- result$data
          variogram_model <- result$variogram_model
          n_tr <- nrow(data)
          if (n_tr < 2) {
            return(invisible())
          }

          stp_col <- paste0("ihfc2024_obs_strip", strip_idx, "_proj")
          obs_col <- paste0("ihfc2024_obs_strip", strip_idx, "_trend")
          sim_col <- paste0("lucazeau2019_sim_strip", strip_idx, "_trend")
          krg_col <- make_krg_trend_col(variogram_model, strip_idx)

          obs_vecs <- get_method_trend_vectors(data, obs_col, n_grid)
          sim_vecs <- get_method_trend_vectors(data, sim_col, n_grid)
          krg_vecs <- get_method_trend_vectors(data, krg_col, n_grid)

          map_dfr(seq_len(n_tr - 1), function(i) {
            id_current <- str_remove_all(data$submap_transect_id[i], "^[^0-9]+")
            id_next <- str_remove_all(data$submap_transect_id[i + 1], "^[^0-9]+")

            obs_stp <- data[[stp_col]]

            n_current <- if (length(obs_stp[[i]]) == 0) 0 else nrow(obs_stp[[i]]) %||% 0
            n_next <- if (length(obs_stp[[i + 1]]) == 0) 0 else nrow(obs_stp[[i + 1]]) %||% 0

            n_current <- if (is.null(n_current) || length(n_current) == 0) 0 else n_current
            n_next <- if (is.null(n_next) || length(n_next) == 0) 0 else n_next

            n_total <- sum(n_current, n_next, na.rm = TRUE)

            tibble(
              submap_transect_set = set_id,
              region = str_extract(set_id, "^[^_]+"),
              pair_idx = i,
              pair_label = paste0(id_current, "--", id_next),
              n_total = n_total,
              obs_ccf = max_ccf_value(obs_vecs[[i]], obs_vecs[[i + 1]], max_lag_frac, n_grid),
              sim_ccf = max_ccf_value(sim_vecs[[i]], sim_vecs[[i + 1]], max_lag_frac, n_grid),
              krg_ccf = max_ccf_value(krg_vecs[[i]], krg_vecs[[i + 1]], max_lag_frac, n_grid)
            )
          })
        },
        default = NULL,
        quiet = quiet
      )
    })
  }))

  with_error_handling(
    {
      if (is.null(global_ccf_summary) || nrow(global_ccf_summary) == 0) {
        if (!quiet) cat(" !! Warning: no CCF data computed for global overview\n")
        return(empty)
      }

      global_ccf_summary |>
        mutate(set_num = as.integer(str_extract(submap_transect_set, "(?<=SET)\\d+"))) |>
        arrange(region, set_num, pair_idx) |>
        select(-set_num) |>
        group_by(region) |>
        mutate(
          global_pair_uid = paste(submap_transect_set, pair_idx, sep = "_"),
          global_pair_uid = factor(global_pair_uid, levels = unique(global_pair_uid))
        ) |>
        ungroup()
    },
    default = empty,
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summarize_method_ccf <- function(data_dir, out_dir, submap_transect_sets, strip_idx = 1, n_grid = 500, max_lag_frac = 0.1, quiet = FALSE) {
  empty <- tibble(
    submap_transect_set = character(),
    region = character(),
    transect_id = character(),
    obs_sim = numeric(),
    obs_krg = numeric(),
    sim_krg = numeric()
  )

  if (!quiet) cat(" .. Summarizing cross-method comparison CCF profiles\n")

  precompile_set_data(data_dir = data_dir, out_dir = out_dir, submap_transect_sets = submap_transect_sets, quiet = quiet)

  method_ccf_summary <- suppressWarnings(suppressMessages({
    map_dfr(submap_transect_sets, function(set_id) {
      with_error_handling(
        {
          result <- load_set_data(data_dir = data_dir, out_dir = out_dir, submap_transect_set = set_id, quiet = quiet)
          if (is.null(result)) {
            return(invisible())
          }

          data <- result$data
          variogram_model <- result$variogram_model

          stp_col <- paste0("ihfc2024_obs_strip", strip_idx, "_proj")
          obs_col <- paste0("ihfc2024_obs_strip", strip_idx, "_trend")
          sim_col <- paste0("lucazeau2019_sim_strip", strip_idx, "_trend")
          krg_col <- make_krg_trend_col(variogram_model, strip_idx)

          obs_vecs <- get_method_trend_vectors(data, obs_col, n_grid)
          sim_vecs <- get_method_trend_vectors(data, sim_col, n_grid)
          krg_vecs <- get_method_trend_vectors(data, krg_col, n_grid)

          map_dfr(seq_len(nrow(data)), function(i) {
            obs_stp <- data[[stp_col]]

            n <- if (length(obs_stp[[i]]) == 0) 0 else nrow(obs_stp[[i]]) %||% 0
            n <- if (is.null(n) || length(n) == 0) 0 else n

            tibble(
              submap_transect_set = set_id,
              region = str_extract(set_id, "^[^_]+"),
              transect_id = data$submap_transect_id[i],
              n = n,
              obs_sim = max_ccf_value(obs_vecs[[i]], sim_vecs[[i]], max_lag_frac, n_grid),
              obs_krg = max_ccf_value(obs_vecs[[i]], krg_vecs[[i]], max_lag_frac, n_grid),
              sim_krg = max_ccf_value(sim_vecs[[i]], krg_vecs[[i]], max_lag_frac, n_grid)
            )
          })
        },
        default = NULL,
        quiet = quiet
      )
    })
  }))

  with_error_handling(
    {
      if (is.null(method_ccf_summary) || nrow(method_ccf_summary) == 0) {
        if (!quiet) cat(" !! Warning: no CCF data computed for method comparison\n")
        return(empty)
      }

      method_ccf_summary |>
        mutate(
          transect_idx = as.integer(str_extract(transect_id, "\\d+$")),
          set_num = as.integer(str_extract(submap_transect_set, "(?<=SET)\\d+"))
        ) |>
        arrange(region, set_num, transect_idx) |>
        group_by(region) |>
        mutate(
          global_pair_uid = paste(submap_transect_set, transect_idx, sep = "_"),
          global_pair_uid = factor(global_pair_uid, levels = unique(global_pair_uid))
        ) |>
        ungroup() |>
        select(-transect_idx, -set_num)
    },
    default = empty,
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
add_transect_set_regime_classifications <- function(global_ccf_summary, method_ccf_summary, quiet = FALSE) {
  if (nrow(global_ccf_summary) == 0 || nrow(method_ccf_summary) == 0) {
    return(list(global = global_ccf_summary, method = method_ccf_summary))
  }

  with_error_handling(
    {
      method_metrics <- method_ccf_summary |>
        group_by(region, submap_transect_set) |>
        summarize(
          max_n = max(n, na.rm = TRUE),
          mean_n = mean(n, na.rm = TRUE),
          median_n = median(n, na.rm = TRUE),
          mean_obs_krg = mean(obs_krg, na.rm = TRUE),
          .groups = "drop"
        )

      global_metrics <- global_ccf_summary |>
        group_by(region, submap_transect_set) |>
        summarize(
          valid_obs_ccf_count = sum(!is.na(obs_ccf)),
          min_obs_ccf = min(obs_ccf, na.rm = TRUE),
          min_krg_ccf = min(krg_ccf, na.rm = TRUE),
          min_sim_ccf = min(sim_ccf, na.rm = TRUE),
          max_sim_ccf = max(sim_ccf, na.rm = TRUE),
          .groups = "drop"
        )

      regime_map <- method_metrics |>
        left_join(global_metrics, by = c("region", "submap_transect_set")) |>
        mutate(
          interpretation = case_when(
            mean_n < 25 ~ "ambiguous",
            mean_obs_krg >= 0.70 & (min_obs_ccf < 0 | min_krg_ccf < 0.30) ~ "discont",
            (mean_obs_krg < 0.30) | (min_obs_ccf < 0 & max_sim_ccf > 0.50) ~ "disagree",
            (
              mean_obs_krg >= 0.70 &
                min_obs_ccf >= 0.50 &
                min_krg_ccf >= 0.50 &
                min_sim_ccf >= 0.50 &
                !is.na(valid_obs_ccf_count) &
                valid_obs_ccf_count >= 2
            ) ~ "cont",
            TRUE ~ "ambiguous"
          )
        ) |>
        select(region, submap_transect_set, interpretation)

      global_enriched <- global_ccf_summary |> left_join(regime_map, by = c("region", "submap_transect_set"))
      method_enriched <- method_ccf_summary |> left_join(regime_map, by = c("region", "submap_transect_set"))

      return(list(global = global_enriched, method = method_enriched))
    },
    default = list(global = NULL, method = NULL),
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summarize_case_study_ccf <- function(data_dir, out_dir, submap_transect_sets, strip_idx = 1, n_grid = 500, max_lag_frac = 0.1, quiet = FALSE) {
  empty <- tibble(
    submap_transect_set = character(),
    region = character(),
    pair_idx = numeric(),
    pair_label = character(),
    obs_ccf = numeric(),
    sim_ccf = numeric(),
    krg_ccf = numeric(),
    global_pair_uid = factor()
  )

  if (!quiet) cat(" .. Summarizing case study CCF matrices across submap transect sets\n")

  precompile_set_data(data_dir = data_dir, out_dir = out_dir, submap_transect_sets = submap_transect_sets, quiet = quiet)

  case_ccf_summary <- suppressWarnings(suppressMessages({
    map_dfr(submap_transect_sets, function(set_id) {
      with_error_handling(
        {
          result <- load_set_data(data_dir = data_dir, out_dir = out_dir, submap_transect_set = set_id, quiet = quiet)
          if (is.null(result)) {
            return(invisible())
          }

          data <- result$data
          variogram_model <- result$variogram_model
          n_tr <- nrow(data)
          if (n_tr < 2) {
            return(invisible())
          }

          obs_col <- paste0("ihfc2024_obs_strip", strip_idx, "_trend")
          sim_col <- paste0("lucazeau2019_sim_strip", strip_idx, "_trend")
          krg_col <- make_krg_trend_col(variogram_model, strip_idx)

          obs_vecs <- get_method_trend_vectors(data, obs_col, n_grid)
          sim_vecs <- get_method_trend_vectors(data, sim_col, n_grid)
          krg_vecs <- get_method_trend_vectors(data, krg_col, n_grid)

          grid_pairs <- expand_grid(i = seq_len(n_tr), j = seq_len(n_tr)) |>
            filter(i < j)

          map_dfr(seq_len(nrow(grid_pairs)), function(idx) {
            i <- grid_pairs$i[idx]
            j <- grid_pairs$j[idx]

            id_i <- str_remove_all(data$submap_transect_id[i], "^[^0-9]+")
            id_j <- str_remove_all(data$submap_transect_id[j], "^[^0-9]+")

            tibble(
              submap_transect_set = set_id,
              region = str_extract(set_id, "^[^_]+"),
              pair_idx = idx,
              pair_label = paste0(id_i, "--", id_j),
              obs_ccf = max_ccf_value(obs_vecs[[i]], obs_vecs[[j]], max_lag_frac, n_grid),
              sim_ccf = max_ccf_value(sim_vecs[[i]], sim_vecs[[j]], max_lag_frac, n_grid),
              krg_ccf = max_ccf_value(krg_vecs[[i]], krg_vecs[[j]], max_lag_frac, n_grid)
            )
          })
        },
        default = NULL,
        quiet = quiet
      )
    })
  }))

  with_error_handling(
    {
      if (is.null(case_ccf_summary) || nrow(case_ccf_summary) == 0) {
        if (!quiet) cat(" !! Warning: no CCF data computed for case study summary\n")
        return(empty)
      }

      case_ccf_summary |>
        mutate(set_num = as.integer(str_extract(submap_transect_set, "(?<=SET)\\d+"))) |>
        arrange(region, set_num, pair_idx) |>
        select(-set_num) |>
        group_by(region) |>
        mutate(
          global_pair_uid = paste(submap_transect_set, pair_idx, sep = "_"),
          global_pair_uid = factor(global_pair_uid, levels = unique(global_pair_uid))
        ) |>
        ungroup()
    },
    default = empty,
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_nlopt_summary_tables <- function(out_dir, nlopt_summary, quiet = FALSE) {
  if (is.null(nlopt_summary) || nrow(nlopt_summary) == 0) {
    return(invisible())
  }

  with_error_handling(
    {
      csv_df <- nlopt_summary |>
        mutate(
          max_point_pair_distance = max_point_pair_distance / 1e3,
          region = str_extract(submap_transect_set, "^[^_]+")
        ) |>
        select(
          submap_transect_set,
          region,
          variogram_model,
          itr,
          max_point_pair_distance,
          variogram_cutoff,
          n_variogram_lags,
          n_heatflow_obs,
          max_n_point_pairs,
          total_cost
        ) |>
        mutate(set_num = as.integer(str_extract(submap_transect_set, "(?<=SET)\\d+"))) |>
        arrange(region, set_num) |>
        select(-set_num)

      csv_dir <- file.path(out_dir, "csv")
      create_dir(csv_dir)

      out_path_csv <- file.path(csv_dir, "nlopt-summary.csv")
      if (!quiet) cat(" -> ", out_path_csv, "\n", sep = "")

      write_csv(csv_df, out_path_csv)

      md_dir <- file.path(out_dir, "markdown")
      create_dir(md_dir)

      out_path_md <- file.path(md_dir, "nlopt-summary.md")
      if (!quiet) cat(" -> ", out_path_md, "\n", sep = "")

      formatted_summary <- nlopt_summary |>
        mutate(
          max_point_pair_distance = max_point_pair_distance / 1e3,
          Region = str_extract(submap_transect_set, "^[^_]+")
        ) |>
        select(
          `SubMap Set` = submap_transect_set,
          Region,
          Model = variogram_model,
          Iteration = itr,
          `Max Distance` = max_point_pair_distance,
          `Vgram Cutoff` = variogram_cutoff,
          `Vgram Lags` = n_variogram_lags,
          `Observations` = n_heatflow_obs,
          `Max Pairs` = max_n_point_pairs,
          `Total Cost` = total_cost
        ) |>
        mutate(
          `Max Distance` = formatC(`Max Distance`, format = "f", digits = 2),
          `Vgram Cutoff` = formatC(`Vgram Cutoff`, format = "f", digits = 1),
          `Vgram Lags` = formatC(`Vgram Lags`, format = "f", digits = 1),
          `Max Pairs` = formatC(`Max Pairs`, format = "f", digits = 1),
          `Total Cost` = formatC(`Total Cost`, format = "f", digits = 4)
        ) |>
        mutate(set_num = as.integer(str_extract(`SubMap Set`, "(?<=SET)\\d+"))) |>
        arrange(Region, set_num) |>
        select(-set_num)

      col_align <- c("l", rep("r", ncol(formatted_summary) - 1))
      md_table <- kable(formatted_summary, format = "pipe", escape = FALSE, align = col_align)
      writeLines(md_table, out_path_md)
    },
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_ihfc2024_summary_tables <- function(data_dir, out_dir, submap_transect_sets, quiet = FALSE) {
  f <- function(submap_transect_set) {
    with_error_handling(
      {
        cache_path <- file.path(out_dir, "submap", paste0(submap_transect_set, ".RData"))

        data <- compile_submap_transect_set_data(
          cache_path = cache_path,
          data_dir = data_dir,
          out_dir = out_dir,
          id = submap_transect_set,
          quiet = quiet
        )

        hf <- data$ihfc2024_obs_hull[[1]]
        obs <- hf$ihfc2024_obs

        row_summary <-
          tibble(
            submap_transect_set = submap_transect_set,
            region = str_extract(submap_transect_set, "^[^_]+"),
            n = length(na.omit(obs)),
            min = min(obs, na.rm = TRUE),
            max = max(obs, na.rm = TRUE),
            median = median(obs, na.rm = TRUE),
            iqr = IQR(obs, na.rm = TRUE),
            mean = mean(obs, na.rm = TRUE),
            sigma = sd(obs, na.rm = TRUE)
          )

        return(row_summary)
      },
      default = NULL,
      quiet = quiet
    )
  }

  if (!quiet) cat(" .. Summarizing IHFC 2024 observations across submap transect sets\n")

  summary_data <- with_error_handling(map(submap_transect_sets, ~ f(.x)) |> bind_rows(), default = NULL, quiet = quiet)

  if (!is.null(summary_data) && nrow(summary_data) > 0) {
    with_error_handling(
      {
        csv_dir <- file.path(out_dir, "csv")
        create_dir(csv_dir)

        out_path_csv <- file.path(csv_dir, "ihfc2024-summary.csv")
        if (!quiet) cat(" -> ", out_path_csv, "\n", sep = "")

        write_csv(summary_data, out_path_csv)

        md_dir <- file.path(out_dir, "markdown")
        create_dir(md_dir)

        out_path_md <- file.path(md_dir, "ihfc2024-summary.md")
        if (!quiet) cat(" -> ", out_path_md, "\n", sep = "")

        formatted_summary <- summary_data |>
          select(
            `SubMap Set` = submap_transect_set,
            Region = region,
            n,
            Min = min,
            Max = max,
            Median = median,
            IQR = iqr,
            Mean = mean,
            `$\\sigma$` = sigma
          ) |>
          mutate(set_num = as.integer(str_extract(`SubMap Set`, "(?<=SET)\\d+"))) |>
          arrange(Region, set_num) |>
          select(-set_num) |>
          mutate(across(where(is.numeric) & !n, ~ formatC(.x, format = "f", digits = 1)))

        col_align <- c("l", rep("r", ncol(formatted_summary) - 1))
        md_table <- kable(formatted_summary, format = "pipe", escape = FALSE, align = col_align)
        writeLines(md_table, out_path_md)
      },
      quiet = quiet
    )
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_submap_transect_ihfc2024_summary_tables <- function(data_dir, out_dir, submap_transect_sets, strip_idx = 1, quiet = FALSE) {
  f <- function(submap_transect_set) {
    with_error_handling(
      {
        cache_path <- file.path(out_dir, "submap", paste0(submap_transect_set, ".RData"))

        data <- compile_submap_transect_set_data(
          cache_path = cache_path,
          data_dir = data_dir,
          out_dir = out_dir,
          id = submap_transect_set,
          quiet = quiet
        )

        if (is.null(data) || nrow(data) == 0) {
          return(NULL)
        }

        region_name <- str_extract(submap_transect_set, "^[^_]+")
        obs_col <- paste0("ihfc2024_obs_strip", strip_idx, "_proj")

        transect_summaries <- map(seq_len(nrow(data)), function(idx) {
          transect_data <- data[idx, ]
          transect_id <- transect_data$submap_transect_id
          trench_name <- transect_data$trench_name

          empty <- tibble(
            submap_transect_set = submap_transect_set,
            submap_transect_id = transect_id,
            trench_name = trench_name,
            region = region_name,
            n = 0,
            min = NA_real_,
            max = NA_real_,
            median = NA_real_,
            iqr = NA_real_,
            mean = NA_real_,
            sigma = NA_real_
          )

          if (obs_col %in% names(transect_data) && !is.null(transect_data[[obs_col]][[1]])) {
            strip_data <- transect_data[[obs_col]][[1]]
            obs <- if ("hf" %in% names(strip_data)) strip_data$hf else NA_real_
          } else {
            obs <- NA_real_
          }

          clean_obs <- na.omit(obs)
          n_obs <- length(clean_obs)

          if (n_obs == 0) {
            return(empty)
          }

          tibble(
            submap_transect_set = submap_transect_set,
            submap_transect_id = transect_id,
            trench_name = trench_name,
            region = region_name,
            n = n_obs,
            min = min(clean_obs, na.rm = TRUE),
            max = max(clean_obs, na.rm = TRUE),
            median = median(clean_obs, na.rm = TRUE),
            iqr = IQR(clean_obs, na.rm = TRUE),
            mean = mean(clean_obs, na.rm = TRUE),
            sigma = sd(clean_obs, na.rm = TRUE)
          )
        }) |> bind_rows()

        return(transect_summaries)
      },
      default = NULL,
      quiet = quiet
    )
  }

  if (!quiet) cat(" .. Summarizing IHFC 2024 observations across individual submap transects\n")

  summary_data <- with_error_handling(map(submap_transect_sets, ~ f(.x)) |> bind_rows(), default = NULL, quiet = quiet)

  if (!is.null(summary_data) && nrow(summary_data) > 0) {
    with_error_handling(
      {
        csv_dir <- file.path(out_dir, "csv")
        create_dir(csv_dir)

        out_path_csv <- file.path(csv_dir, "submap-transect-ihfc2024-summary.csv")
        if (!quiet) cat(" -> ", out_path_csv, "\n", sep = "")

        write_csv(summary_data, out_path_csv)

        md_dir <- file.path(out_dir, "markdown")
        create_dir(md_dir)

        out_path_md <- file.path(md_dir, "submap-transect-ihfc2024-summary.md")
        if (!quiet) cat(" -> ", out_path_md, "\n", sep = "")

        formatted_summary <- summary_data |>
          select(
            `SubMap Set` = submap_transect_set,
            `Transect ID` = submap_transect_id,
            `Trench Name` = trench_name,
            Region = region,
            n,
            Min = min,
            Max = max,
            Median = median,
            IQR = iqr,
            Mean = mean,
            `$\\sigma$` = sigma
          ) |>
          mutate(
            set_num = as.integer(str_extract(`SubMap Set`, "\\d+")),
            transect_num = as.integer(str_extract(`Transect ID`, "\\d+"))
          ) |>
          arrange(Region, set_num, transect_num) |>
          select(-set_num, -transect_num) |>
          mutate(across(where(is.numeric) & !n, ~ ifelse(is.na(.x), "", formatC(.x, format = "f", digits = 1))))

        col_align <- c("l", rep("r", ncol(formatted_summary) - 1))
        md_table <- kable(formatted_summary, format = "pipe", escape = FALSE, align = col_align)
        writeLines(md_table, out_path_md)
      },
      quiet = quiet
    )
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_est_diff_summary_tables <- function(out_dir, interp_est_diff_summary, quiet = FALSE) {
  if (is.null(interp_est_diff_summary) || nrow(interp_est_diff_summary) == 0) {
    return(invisible())
  }

  with_error_handling(
    {
      csv_dir <- file.path(out_dir, "csv")
      create_dir(csv_dir)

      out_path_csv <- file.path(csv_dir, "point-by-point-est-diff-summary.csv")
      if (!quiet) cat(" -> ", out_path_csv, "\n", sep = "")

      write_csv(interp_est_diff_summary, out_path_csv)

      md_dir <- file.path(out_dir, "markdown")
      create_dir(md_dir)

      out_path_md <- file.path(md_dir, "point-by-point-est-diff-summary.md")
      if (!quiet) cat(" -> ", out_path_md, "\n", sep = "")

      formatted_summary <- interp_est_diff_summary |>
        mutate(Region = str_extract(submap_transect_set, "^[^_]+")) |>
        select(
          `SubMap Set` = submap_transect_set,
          Region,
          n,
          Min = min,
          Max = max,
          Median = median,
          IQR = iqr,
          Mean = mean,
          `$\\sigma$` = sigma
        ) |>
        mutate(set_num = as.integer(str_extract(`SubMap Set`, "(?<=SET)\\d+"))) |>
        arrange(Region, set_num) |>
        select(-set_num) |>
        mutate(across(where(is.numeric) & !n, ~ formatC(.x, format = "f", digits = 1)))

      col_align <- c("l", rep("r", ncol(formatted_summary) - 1))
      md_table <- kable(formatted_summary, format = "pipe", escape = FALSE, align = col_align)
      writeLines(md_table, out_path_md)
    },
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_sigma_diff_summary_tables <- function(out_dir, interp_sigma_diff_summary, quiet = FALSE) {
  if (is.null(interp_sigma_diff_summary) || nrow(interp_sigma_diff_summary) == 0) {
    return(invisible())
  }

  with_error_handling(
    {
      csv_dir <- file.path(out_dir, "csv")
      create_dir(csv_dir)

      out_path_csv <- file.path(csv_dir, "point-by-point-sigma-diff-summary.csv")
      if (!quiet) cat(" -> ", out_path_csv, "\n", sep = "")

      write_csv(interp_sigma_diff_summary, out_path_csv)

      md_dir <- file.path(out_dir, "markdown")
      create_dir(md_dir)

      out_path_md <- file.path(md_dir, "point-by-point-sigma-diff-summary.md")
      if (!quiet) cat(" -> ", out_path_md, "\n", sep = "")

      formatted_summary <- interp_sigma_diff_summary |>
        mutate(Region = str_extract(submap_transect_set, "^[^_]+")) |>
        select(
          `SubMap Set` = submap_transect_set,
          Region,
          n,
          Min = min,
          Max = max,
          Median = median,
          IQR = iqr,
          Mean = mean,
          `$\\sigma$` = sigma
        ) |>
        mutate(set_num = as.integer(str_extract(`SubMap Set`, "(?<=SET)\\d+"))) |>
        arrange(Region, set_num) |>
        select(-set_num) |>
        mutate(across(where(is.numeric) & !n, ~ formatC(.x, format = "f", digits = 1)))

      col_align <- c("l", rep("r", ncol(formatted_summary) - 1))
      md_table <- kable(formatted_summary, format = "pipe", escape = FALSE, align = col_align)
      writeLines(md_table, out_path_md)
    },
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_similarity_residual_summary_tables <- function(out_dir, interp_accuracy_summary, quiet = FALSE) {
  if (is.null(interp_accuracy_summary) || nrow(interp_accuracy_summary) == 0) {
    return(invisible())
  }

  with_error_handling(
    {
      csv_dir <- file.path(out_dir, "csv")
      create_dir(csv_dir)

      out_path_csv <- file.path(csv_dir, "point-by-point-sim-residual-summary.csv")
      if (!quiet) cat(" -> ", out_path_csv, "\n", sep = "")

      df_sim <- interp_accuracy_summary |>
        select(
          submap_transect_set,
          n_similarity,
          min_similarity,
          max_similarity,
          med_similarity,
          iqr_similarity,
          mean_similarity,
          sd_similarity,
          rmse_similarity
        )

      write_csv(df_sim, out_path_csv)

      md_dir <- file.path(out_dir, "markdown")
      create_dir(md_dir)

      out_path_md <- file.path(md_dir, "point-by-point-sim-residual-summary.md")
      if (!quiet) cat(" -> ", out_path_md, "\n", sep = "")

      formatted_summary <- interp_accuracy_summary |>
        mutate(Region = str_extract(submap_transect_set, "^[^_]+")) |>
        select(
          `SubMap Set` = submap_transect_set,
          Region,
          n = n_similarity,
          Min = min_similarity,
          Max = max_similarity,
          Median = med_similarity,
          IQR = iqr_similarity,
          Mean = mean_similarity,
          `$\\sigma$` = sd_similarity,
          RMSE = rmse_similarity
        ) |>
        mutate(set_num = as.integer(str_extract(`SubMap Set`, "(?<=SET)\\d+"))) |>
        arrange(Region, set_num) |>
        select(-set_num) |>
        mutate(across(where(is.numeric) & !n, ~ formatC(.x, format = "f", digits = 1)))

      col_align <- c("l", rep("r", ncol(formatted_summary) - 1))
      md_table <- kable(formatted_summary, format = "pipe", escape = FALSE, align = col_align)
      writeLines(md_table, out_path_md)
    },
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_krige_residual_summary_tables <- function(out_dir, interp_accuracy_summary, quiet = FALSE) {
  if (is.null(interp_accuracy_summary) || nrow(interp_accuracy_summary) == 0) {
    return(invisible())
  }

  with_error_handling(
    {
      csv_dir <- file.path(out_dir, "csv")
      create_dir(csv_dir)

      out_path_csv <- file.path(csv_dir, "point-by-point-krg-residual-summary.csv")
      if (!quiet) cat(" -> ", out_path_csv, "\n", sep = "")

      df_krg <- interp_accuracy_summary |>
        select(
          submap_transect_set,
          n_krige,
          min_krige,
          max_krige,
          med_krige,
          iqr_krige,
          mean_krige,
          sd_krige,
          rmse_krige
        )

      write_csv(df_krg, out_path_csv)

      md_dir <- file.path(out_dir, "markdown")
      create_dir(md_dir)

      out_path_md <- file.path(md_dir, "point-by-point-krg-residual-summary.md")
      if (!quiet) cat(" -> ", out_path_md, "\n", sep = "")

      formatted_summary <- interp_accuracy_summary |>
        mutate(Region = str_extract(submap_transect_set, "^[^_]+")) |>
        select(
          `SubMap Set` = submap_transect_set,
          Region,
          n = n_krige,
          Min = min_krige,
          Max = max_krige,
          Median = med_krige,
          IQR = iqr_krige,
          Mean = mean_krige,
          `$\\sigma$` = sd_krige,
          RMSE = rmse_krige
        ) |>
        mutate(set_num = as.integer(str_extract(`SubMap Set`, "(?<=SET)\\d+"))) |>
        arrange(Region, set_num) |>
        select(-set_num) |>
        mutate(across(where(is.numeric) & !n, ~ formatC(.x, format = "f", digits = 1)))

      col_align <- c("l", rep("r", ncol(formatted_summary) - 1))
      md_table <- kable(formatted_summary, format = "pipe", escape = FALSE, align = col_align)
      writeLines(md_table, out_path_md)
    },
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_global_ccf_summary_tables <- function(out_dir, global_ccf_summary, strip_idx = 1, quiet = FALSE) {
  if (is.null(global_ccf_summary) || nrow(global_ccf_summary) == 0) {
    return(invisible())
  }

  with_error_handling(
    {
      csv_dir <- file.path(out_dir, "csv")
      create_dir(csv_dir)

      out_path_csv <- file.path(csv_dir, paste0("global-ccf-summary-strip", strip_idx, ".csv"))
      if (!quiet) cat(" -> ", out_path_csv, "\n", sep = "")

      write_csv(global_ccf_summary, out_path_csv)

      md_dir <- file.path(out_dir, "markdown")
      create_dir(md_dir)

      out_path_md <- file.path(md_dir, paste0("global-ccf-summary-strip", strip_idx, ".md"))
      if (!quiet) cat(" -> ", out_path_md, "\n", sep = "")

      formatted_summary <- global_ccf_summary |>
        mutate(Region = str_extract(submap_transect_set, "^[^_]+")) |>
        select(
          `SubMap Set` = submap_transect_set,
          Region,
          `Transect Pair` = pair_label,
          `Obs CCF` = obs_ccf,
          `Sim CCF` = sim_ccf,
          `Krg CCF` = krg_ccf,
          `Interpretation` = interpretation,
          `UID` = global_pair_uid
        ) |>
        mutate(set_num = as.integer(str_extract(`SubMap Set`, "(?<=SET)\\d+"))) |>
        arrange(Region, set_num) |>
        select(-set_num) |>
        mutate(across(c(`Obs CCF`, `Sim CCF`, `Krg CCF`), ~ ifelse(is.na(.x), "", formatC(.x, format = "f", digits = 4))))

      col_align <- c("l", rep("r", ncol(formatted_summary) - 1))
      md_table <- kable(formatted_summary, format = "pipe", escape = FALSE, align = col_align)
      writeLines(md_table, out_path_md)
    },
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_method_ccf_summary_tables <- function(out_dir, method_ccf_summary, strip_idx = 1, quiet = FALSE) {
  if (is.null(method_ccf_summary) || nrow(method_ccf_summary) == 0) {
    return(invisible())
  }

  with_error_handling(
    {
      csv_dir <- file.path(out_dir, "csv")
      create_dir(csv_dir)

      out_path_csv <- file.path(csv_dir, paste0("method-ccf-summary-strip", strip_idx, ".csv"))
      if (!quiet) cat(" -> ", out_path_csv, "\n", sep = "")

      write_csv(method_ccf_summary, out_path_csv)

      md_dir <- file.path(out_dir, "markdown")
      create_dir(md_dir)

      out_path_md <- file.path(md_dir, paste0("method-ccf-summary-strip", strip_idx, ".md"))
      if (!quiet) cat(" -> ", out_path_md, "\n", sep = "")

      formatted_summary <- method_ccf_summary |>
        mutate(Region = str_extract(submap_transect_set, "^[^_]+")) |>
        select(
          `SubMap Set` = submap_transect_set,
          Region,
          `Transect ID` = transect_id,
          `Obs--Sim CCF` = obs_sim,
          `Obs--Krg CCF` = obs_krg,
          `Sim--Krg CCF` = sim_krg,
          `Interpretation` = interpretation
        ) |>
        mutate(set_num = as.integer(str_extract(`SubMap Set`, "(?<=SET)\\d+"))) |>
        arrange(Region, set_num) |>
        select(-set_num) |>
        mutate(across(c(`Obs--Sim CCF`, `Obs--Krg CCF`, `Sim--Krg CCF`), ~ ifelse(is.na(.x), "", formatC(.x, format = "f", digits = 4))))

      col_align <- c("l", rep("r", ncol(formatted_summary) - 1))
      md_table <- kable(formatted_summary, format = "pipe", escape = FALSE, align = col_align)
      writeLines(md_table, out_path_md)
    },
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_case_study_ccf_summary_tables <- function(out_dir, case_ccf_summary, strip_idx = 1, quiet = FALSE) {
  if (is.null(case_ccf_summary) || nrow(case_ccf_summary) == 0) {
    return(invisible())
  }

  with_error_handling(
    {
      csv_dir <- file.path(out_dir, "csv")
      create_dir(csv_dir)

      out_path_csv <- file.path(csv_dir, paste0("case-study-ccf-summary-strip", strip_idx, ".csv"))
      if (!quiet) cat(" -> ", out_path_csv, "\n", sep = "")
      write_csv(case_ccf_summary, out_path_csv)

      md_dir <- file.path(out_dir, "markdown")
      create_dir(md_dir)

      out_path_md <- file.path(md_dir, paste0("case-study-ccf-summary-strip", strip_idx, ".md"))
      if (!quiet) cat(" -> ", out_path_md, "\n", sep = "")

      formatted_summary <- case_ccf_summary |>
        select(
          `SubMap Set` = submap_transect_set,
          Region = region,
          `Transect Pair` = pair_label,
          `Obs CCF` = obs_ccf,
          `Sim CCF` = sim_ccf,
          `Krg CCF` = krg_ccf,
          `UID` = global_pair_uid
        ) |>
        mutate(set_num = as.integer(str_extract(`SubMap Set`, "(?<=SET)\\d+"))) |>
        arrange(Region, set_num) |>
        select(-set_num) |>
        mutate(across(c(`Obs CCF`, `Sim CCF`, `Krg CCF`), ~ ifelse(is.na(.x), "", formatC(.x, format = "f", digits = 4))))

      col_align <- c("l", rep("r", ncol(formatted_summary) - 1))
      md_table <- kable(formatted_summary, format = "pipe", escape = FALSE, align = col_align)
      writeLines(md_table, out_path_md)
    },
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_submap_transect_set_ihfc2024_obs_csv <- function(data_dir, out_dir, submap_transect_sets, quiet = FALSE) {
  f <- function(submap_transect_set) {
    with_error_handling(
      {
        csv_dir <- file.path(out_dir, "csv", "submap")
        create_dir(csv_dir)

        out_path <- file.path(csv_dir, paste0(submap_transect_set, "-ihfc2024-obs.csv"))

        if (file.exists(out_path)) {
          return(invisible())
        }

        cache_path <- file.path(out_dir, "submap", paste0(submap_transect_set, ".RData"))
        data <- compile_submap_transect_set_data(
          cache_path = cache_path,
          data_dir = data_dir,
          out_dir = out_dir,
          id = submap_transect_set,
          quiet = quiet
        )

        hf <- data$ihfc2024_obs_hull[[1]]
        coords <- st_coordinates(hf)

        hull <- data$hull[[1]]

        if (!quiet) cat(" -> ", out_path, "\n", sep = "")
        hf |>
          mutate(lon = coords[, 1], lat = coords[, 2], .before = ihfc2024_obs) |>
          mutate(across(where(is.numeric), ~ formatC(.x, format = "f", digits = 3))) |>
          st_set_geometry(NULL) |>
          mutate(
            submap_transect_set = data$submap_transect_set,
            submap_transect_ids = hull$submap_transect_ids,
            .before = ihfc2024_obs
          ) |>
          write_csv(out_path)
      },
      default = NULL,
      quiet = quiet
    )
  }

  if (!quiet) cat(" .. Processing IHFC 2024 observations by submap transect sets\n")
  with_error_handling(walk(submap_transect_sets, ~ f(.x)), quiet = quiet)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_submap_transect_set_lucazeau2019_sim_csv <- function(data_dir, out_dir, submap_transect_sets, quiet = FALSE) {
  f <- function(submap_transect_set) {
    with_error_handling(
      {
        csv_dir <- file.path(out_dir, "csv", "submap")
        create_dir(csv_dir)

        out_path <- file.path(csv_dir, paste0(submap_transect_set, "-lucazeau2019-sim.csv"))

        if (file.exists(out_path)) {
          return(invisible())
        }

        cache_path <- file.path(out_dir, "submap", paste0(submap_transect_set, ".RData"))
        data <- compile_submap_transect_set_data(
          cache_path = cache_path,
          data_dir = data_dir,
          out_dir = out_dir,
          id = submap_transect_set,
          quiet = quiet
        )

        hf <- data$lucazeau2019_sim_hull[[1]]
        coords <- st_coordinates(hf)

        hull <- data$hull[[1]]

        if (!quiet) cat(" -> ", out_path, "\n", sep = "")
        hf |>
          mutate(lon = coords[, 1], lat = coords[, 2], .before = lucazeau2019_sim_est) |>
          mutate(across(where(is.numeric), ~ formatC(.x, format = "f", digits = 3))) |>
          st_set_geometry(NULL) |>
          mutate(
            submap_transect_set = data$submap_transect_set,
            submap_transect_ids = hull$submap_transect_ids,
            .before = lucazeau2019_sim_est
          ) |>
          write_csv(out_path)
      },
      default = NULL,
      quiet = quiet
    )
  }

  if (!quiet) cat(" .. Processing Lucazeau 2019 similarity interpolation by submap transect sets\n")
  with_error_handling(walk(submap_transect_sets, ~ f(.x)), quiet = quiet)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_submap_transect_set_kerswell2025_krg_csv <- function(data_dir, out_dir, submap_transect_sets, quiet = FALSE) {
  f <- function(submap_transect_set) {
    with_error_handling(
      {
        csv_dir <- file.path(out_dir, "csv", "submap")
        create_dir(csv_dir)

        out_path <- file.path(csv_dir, paste0(submap_transect_set, "-kerswell2025-krg.csv"))

        if (file.exists(out_path)) {
          return(invisible())
        }

        opt_data <- get_optimal_krige_model(out_dir, submap_transect_set, quiet = quiet)
        variogram_model <- opt_data$optimal_krige_model$variogram_model
        fitted_variogram <- opt_data$fitted_variogram
        max_n_point_pairs <- opt_data$optimal_krige_model$max_n_point_pairs
        max_point_pair_distance <- opt_data$optimal_krige_model$max_point_pair_distance

        if (is.null(variogram_model) && is.null(fitted_variogram) && is.null(max_n_point_pairs) && is.null(max_point_pair_distance)) {
          if (!quiet) cat(" !! Warning: no valid kriging model found for submap transect '", submap_transect_set, "'\n", sep = "")
          return(invisible())
        }

        prefix <- "kerswell2025_krg_"
        model_str <- str_to_lower(variogram_model)
        krg_hull_str <- paste0(prefix, model_str, "_hull")

        cache_path <- file.path(out_dir, "submap", paste0(submap_transect_set, ".RData"))
        data <- compile_submap_transect_set_data(
          cache_path = cache_path,
          data_dir = data_dir,
          out_dir = out_dir,
          id = submap_transect_set,
          variogram_model = variogram_model,
          fitted_variogram = fitted_variogram,
          max_n_point_pairs = max_n_point_pairs,
          max_point_pair_distance = max_point_pair_distance,
          quiet = quiet
        )

        hf <- data[[krg_hull_str]][[1]]
        coords <- st_coordinates(hf)

        hull <- data$hull[[1]]

        if (!quiet) cat(" -> ", out_path, "\n", sep = "")
        hf |>
          mutate(lon = coords[, 1], lat = coords[, 2], .before = kerswell2025_krg_est) |>
          mutate(across(where(is.numeric), ~ formatC(.x, format = "f", digits = 3))) |>
          st_set_geometry(NULL) |>
          mutate(
            submap_transect_set = data$submap_transect_set,
            submap_transect_ids = hull$submap_transect_ids,
            .before = kerswell2025_krg_est
          ) |>
          write_csv(out_path)
      },
      default = NULL,
      quiet = quiet
    )
  }

  if (!quiet) cat(" .. Processing Kerswell 2025 krige interpolation by submap transect sets\n")
  with_error_handling(walk(submap_transect_sets, ~ f(.x)), quiet = quiet)
}
