#######################################################
## Postprocess spatial datasets                      ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
initialize_submap_transect_set_data <- function(data_dir, out_dir, id, ...) {
  with_error_handling({
    dots <- list(...)
    quiet <- dots$quiet %||% FALSE

    if (is.character(id) && all(id %in% sf_hull$submap_transect_set)) {
      hull <- filter(sf_hull, submap_transect_set %in% id)
      submap_ids <- unlist(strsplit(hull$submap_transect_ids, ",")) |> trimws()
      data <- filter(sf_submap, submap_transect_id %in% submap_ids)
    } else {
      stop("Unrecognized input for id! Use the submap_transect_set")
    }

    if (nrow(data) == 0) stop("Submap transect id: ", id, " not found in sf_submap")

    if (!quiet) cat("    Compiling data: '", id, "'\n", sep = "")
    rowwise(data) |>
      mutate(
        hull = st_geometry(hull),
        grid = list(crop_to_bbox(sf_grid, hull)),
        ihfc2024_obs_hull = list(crop_to_bbox(sf_ihfc, hull)),
        lucazeau2019_sim_hull = list(crop_to_bbox(sf_sim, hull))
      ) |>
      ungroup()
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
compile_submap_transect_set_data <- function(cache_path, data_dir, out_dir, id, ...) {
  with_error_handling({
    dots <- list(...)
    variogram_model <- dots$variogram_model %||% NULL
    fitted_variogram <- dots$fitted_variogram %||% NULL
    max_n_point_pairs <- dots$max_n_point_pairs %||% NULL
    max_point_pair_distance <- dots$max_point_pair_distance %||% NULL
    quiet <- dots$quiet %||% FALSE

    if (file.exists(cache_path)) {
      if (!quiet) cat("    Reading cache:", cache_path, "\n")
      load(cache_path)

      if (!exists("df") || is.null(df)) {
        if (!quiet) cat(" -- Found cache, but 'df' not found, will recompile\n")
        df <- initialize_submap_transect_set_data(data_dir, out_dir, id, ...)
        save(df, file = cache_path)
        return(df)
      }

      if (!is.null(variogram_model)) {
        prefix <- "kerswell2025_krg_"
        model_str <- str_to_lower(variogram_model)
        expected_cols <- paste0(prefix, model_str, c("_hull"))
        missing_cols <- expected_cols[!expected_cols %in% names(df)]

        if (!is.null(fitted_variogram) && !is.null(max_n_point_pairs) && !is.null(max_point_pair_distance)) {
          if (!quiet) cat(" -- Adding krige results for variogram model '", model_str, "'\n", sep = "")

          shp_krg <- krige_submap_transect_set(
            df$ihfc2024_obs_hull[[1]],
            "ihfc2024_obs",
            fitted_variogram,
            df$grid[[1]],
            max_n_point_pairs,
            max_point_pair_distance
          )

          krg_hull_str <- paste0(prefix, model_str, "_hull")

          df <- df |>
            rowwise() |>
            mutate(
              !!krg_hull_str := list(crop_to_bbox(shp_krg, hull)),
              krg_sim_diff_hull = list(process_sim_krg_diff(!!sym(krg_hull_str), lucazeau2019_sim_hull))
            ) |>
            ungroup()

          return(df)
        } else {
          return(df)
        }
      } else {
        return(df)
      }
    }

    df <- initialize_submap_transect_set_data(data_dir, out_dir, id, ...)
    if (!quiet) cat("    Caching data:", cache_path, "\n")
    create_path_dir(cache_path)
    save(df, file = cache_path)
    df
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
project_hf_to_transect <- function(transect_geometry, shp_hf) {
  with_error_handling({
    if (st_crs(transect_geometry) != st_crs(shp_hf)) {
      stop("Required identical crs")
    }
    if (!inherits(transect_geometry, "sfc")) {
      stop("Transect_geometry must be an sfc geometry")
    }
    if (!inherits(shp_hf, "sf")) {
      stop("shp_hf must be an sf object")
    }

    if (any(names(shp_hf) == "ihfc2024_obs")) {
      hf <- shp_hf$ihfc2024_obs
      sigma <- NA
    } else if (any(names(shp_hf) == "lucazeau2019_sim_est")) {
      hf <- shp_hf$lucazeau2019_sim_est
      sigma <- shp_hf$lucazeau2019_sim_sigma
    } else if (any(names(shp_hf) == "kerswell2025_krg_est")) {
      hf <- shp_hf$kerswell2025_krg_est
      sigma <- shp_hf$kerswell2025_krg_sigma
    } else {
      stop("Unrecognized heatflow data column")
    }

    suppressWarnings({
      suppressMessages({
        projected_distances <- st_line_project(transect_geometry, st_geometry(shp_hf), normalized = TRUE)

        st_as_sf(st_line_interpolate(transect_geometry, projected_distances)) |>
          rename(geometry = x) |>
          mutate(distance = projected_distances, hf = hf, sigma = sigma, .before = geometry)
      })
    })
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit_loess_to_projected_obervations <- function(shp, n = 1e3, minimum_observations = 10, span = 0.65, ...) {
  with_error_handling({
    dots <- list(...)
    quiet <- dots$quiet %||% FALSE

    df <- shp |> st_set_geometry(NULL)

    if (nrow(df) < minimum_observations) {
      if (!quiet) cat(" !! Warning: cannot fit loess with less than ", minimum_observations, " projected obs!", "\n", sep = "")
      return(invisible())
    } else {
      loess_model <- NULL

      while (span <= 0.9 && is.null(loess_model)) {
        loess_model <- tryCatch(
          {
            loess(hf ~ distance, data = df, span = span)
          },
          error = function(e) {
            invisible()
          }
        )

        if (is.null(loess_model)) {
          span <- span + 0.05
        }
      }

      if (is.null(loess_model) || is.null(loess_model$fitted)) {
        if (!quiet) cat(" !! Warning: loess fitting failed with ", nrow(df), " observations\n", sep = "")
        return(invisible())
      }

      new_distances <- seq(0, 1, length.out = n)
      original_range <- range(df$distance, na.rm = TRUE)

      loess_predictions <- tryCatch(
        {
          predict(loess_model, newdata = new_distances)
        },
        error = function(e) {
          if (!quiet) cat(" !! Warning: loess predictions failed with ", nrow(df), " observations\n", sep = "")
          invisible()
        }
      )

      if (is.null(loess_predictions)) {
        return(invisible())
      }

      loess_predictions[new_distances < original_range[1] | new_distances > original_range[2]] <- NA
      tibble(distance = new_distances, hf = loess_predictions)
    }
  })
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
get_nlopt_iter_paths <- function(nlopt_iter_dir, nlopt_krige_dir, submap_transect_set, quiet = TRUE) {
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
    default = list()
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
read_nlopt_iteration <- function(in_path) {
  empty_tibble <- tibble(
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
    return(empty_tibble)
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
    default = empty_tibble
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_optimal_krige_model <- function(out_dir, submap_transect_set, ...) {
  empty_list <- list(
    "nlopt_iter_df" = NULL,
    "optimal_krige_model" = NULL,
    "experimental_variogram" = NULL,
    "fitted_variogram" = NULL
  )

  with_error_handling(
    {
      dots <- list(...)
      quiet <- dots$quiet %||% FALSE

      nlopt_dir <- file.path(out_dir, "nlopt")
      nlopt_iter_dir <- file.path(nlopt_dir, "iterations")
      nlopt_krige_dir <- file.path(nlopt_dir, "krige_models")

      nlopt_iter_paths <- get_nlopt_iter_paths(nlopt_iter_dir, nlopt_krige_dir, submap_transect_set, ...)

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
    default = empty_list
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summarize_optimal_krige_models <- function(out_dir, submap_transect_sets, ...) {
  with_error_handling({
    dots <- list(...)
    seed <- dots$seed %||% 42
    nprocs <- dots$nprocs %||% availableCores() - 2
    parallel <- dots$parallel %||% FALSE

    inner_dots <- dots
    inner_dots$parallel <- NULL
    inner_dots$nprocs <- NULL
    inner_dots$seed <- NULL

    f <- function(submap_transect_set, ...) get_optimal_krige_model(out_dir, submap_transect_set, ...)$optimal_krige_model

    cat(" -> Summarizing optimal krige models\n")
    f_partial <- partial(f, !!!inner_dots)

    if (parallel) {
      set.seed(seed)
      plan(multisession, workers = nprocs)
      as_tibble(future_map_dfr(submap_transect_sets, f_partial, .options = furrr_options(seed = seed)))
    } else {
      map_df(submap_transect_sets, f_partial)
    }
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summarize_interpolation_differences <- function(data_dir, out_dir, submap_transect_sets, ...) {
  with_error_handling({
    dots <- list(...)
    seed <- dots$seed %||% 42
    nprocs <- dots$nprocs %||% availableCores() - 2
    parallel <- dots$parallel %||% FALSE

    inner_dots <- dots
    inner_dots$parallel <- NULL
    inner_dots$nprocs <- NULL
    inner_dots$seed <- NULL

    f <- function(submap_transect_set, ...) {
      with_error_handling(
        {
          dots <- list(...)
          quiet <- dots$quiet %||% FALSE

          opt_data <- get_optimal_krige_model(out_dir, submap_transect_set, ...)
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
              med = NA_real_,
              iqr = NA_real_,
              mean = NA_real_,
              sigma = NA_real_
            ))
          }

          cache_path <- file.path(out_dir, "submap", paste0(submap_transect_set, ".RData"))
          data <- compile_submap_transect_set_data(
            cache_path,
            data_dir,
            out_dir,
            submap_transect_set,
            variogram_model = variogram_model,
            fitted_variogram = fitted_variogram,
            max_n_point_pairs = max_n_point_pairs,
            max_point_pair_distance = max_point_pair_distance,
            ...
          )

          data$krg_sim_diff_hull[[1]] |>
            st_set_geometry(NULL) |>
            summarise(
              submap_transect_set = submap_transect_set,
              n_grid = n(),
              n = sum(!is.na(sim_krg_diff_est)),
              rmse = sqrt(mean(sim_krg_diff_est^2, na.rm = TRUE)),
              min = min(sim_krg_diff_est, na.rm = TRUE),
              max = max(sim_krg_diff_est, na.rm = TRUE),
              med = median(sim_krg_diff_est, na.rm = TRUE),
              iqr = IQR(sim_krg_diff_est, na.rm = TRUE),
              mean = mean(sim_krg_diff_est, na.rm = TRUE),
              sigma = sd(sim_krg_diff_est, na.rm = TRUE)
            )
        },
        default = NULL
      )
    }

    cat(" -> Summarizing interpolation differences\n")
    f_partial <- partial(f, !!!inner_dots)

    if (parallel) {
      set.seed(seed)
      plan(multisession, workers = nprocs)
      as_tibble(future_map_dfr(submap_transect_sets, f_partial, .options = furrr_options(seed = seed)))
    } else {
      map_df(submap_transect_sets, f_partial)
    }
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_closest_interp_obs_point_pairs <- function(data_dir, out_dir, submap_transect_set, threshold = 1e4, ...) {
  with_error_handling({
    dots <- list(...)
    variogram_model <- dots$variogram_model %||% NULL
    fitted_variogram <- dots$fitted_variogram %||% NULL
    max_n_point_pairs <- dots$max_n_point_pairs %||% NULL
    max_point_pair_distance <- dots$max_point_pair_distance %||% NULL

    prefix <- "kerswell2025_krg_"
    model_str <- str_to_lower(variogram_model)
    krg_hull_str <- paste0(prefix, model_str, "_hull")

    cache_path <- file.path(out_dir, "submap", paste0(submap_transect_set, ".RData"))
    if (!is.null(variogram_model) && !is.null(fitted_variogram) && !is.null(max_n_point_pairs) && !is.null(max_point_pair_distance)) {
      data <- compile_submap_transect_set_data(
        cache_path,
        data_dir,
        out_dir,
        submap_transect_set,
        variogram_model = variogram_model,
        fitted_variogram = fitted_variogram,
        max_n_point_pairs = max_n_point_pairs,
        max_point_pair_distance = max_point_pair_distance,
        ...
      )
      shp_hf_interp <- data[[krg_hull_str]][[1]]
    } else {
      data <- compile_submap_transect_set_data(cache_path, data_dir, out_dir, submap_transect_set, ...)
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
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summarize_interpolation_accuracies <- function(data_dir, out_dir, submap_transect_sets, ...) {
  with_error_handling({
    dots <- list(...)
    seed <- dots$seed %||% 42
    nprocs <- dots$nprocs %||% availableCores() - 2
    parallel <- dots$parallel %||% FALSE

    inner_dots <- dots
    inner_dots$parallel <- NULL
    inner_dots$nprocs <- NULL
    inner_dots$seed <- NULL

    f <- function(submap_transect_set, ...) {
      suppressWarnings({
        suppressMessages({
          with_error_handling(
            {
              dots <- list(...)
              quiet <- dots$quiet %||% FALSE

              opt_data <- get_optimal_krige_model(out_dir, submap_transect_set, ...)
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

              point_pairs_sim <- get_closest_interp_obs_point_pairs(data_dir, out_dir, submap_transect_set, ...)
              point_pairs_krg <- get_closest_interp_obs_point_pairs(
                data_dir,
                out_dir,
                submap_transect_set,
                variogram_model = variogram_model,
                fitted_variogram = fitted_variogram,
                max_n_point_pairs = max_n_point_pairs,
                max_point_pair_distance = max_point_pair_distance,
                ...
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
            default = NULL
          )
        })
      })
    }

    cat(" -> Summarizing interpolation accuracies\n")
    f_partial <- partial(f, !!!inner_dots)

    if (parallel) {
      set.seed(seed)
      plan(multisession, workers = nprocs)
      as_tibble(future_map_dfr(submap_transect_sets, f_partial, .options = furrr_options(seed = seed)))
    } else {
      map_df(submap_transect_sets, f_partial)
    }
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_submap_transect_set_ihfc2024_obs <- function(data_dir, out_dir, submap_transect_sets, ...) {
  with_error_handling({
    f <- function(submap_transect_set, ...) {
      with_error_handling(
        {
          dots <- list(...)
          quiet <- dots$quiet %||% FALSE

          save_dir <- file.path(out_dir, "heatflow", "submap_transects")
          create_dir(save_dir)

          out_path <- file.path(save_dir, paste0(submap_transect_set, "-ihfc2024-obs.csv"))

          if (file.exists(out_path)) {
            if (!quiet) cat(" -- Found data: ", out_path, "\n", sep = "")
            return(invisible())
          }

          cache_path <- file.path(out_dir, "submap", paste0(submap_transect_set, ".RData"))
          data <- compile_submap_transect_set_data(cache_path, data_dir, out_dir, submap_transect_set, ...)

          hf <- data$ihfc2024_obs_hull[[1]]
          coords <- st_coordinates(hf)

          hull <- data$hull[[1]]

          if (!quiet) cat("    Writing data: ", out_path, "\n", sep = "")
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
        default = NULL
      )
    }

    cat(" -> Processing IHFC 2024 observations by submap sransect sets\n")
    walk(submap_transect_sets, f, ...)
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_submap_transect_set_lucazeau2019_sim <- function(data_dir, out_dir, submap_transect_sets, ...) {
  with_error_handling({
    f <- function(submap_transect_set, ...) {
      with_error_handling(
        {
          dots <- list(...)
          quiet <- dots$quiet %||% FALSE

          save_dir <- file.path(out_dir, "heatflow", "submap_transects")
          create_dir(save_dir)

          out_path <- file.path(save_dir, paste0(submap_transect_set, "-lucazeau2019-sim.csv"))

          if (file.exists(out_path)) {
            if (!quiet) cat(" -- Found data: ", out_path, "\n", sep = "")
            return(invisible())
          }

          cache_path <- file.path(out_dir, "submap", paste0(submap_transect_set, ".RData"))
          data <- compile_submap_transect_set_data(cache_path, data_dir, out_dir, submap_transect_set, ...)

          hf <- data$lucazeau2019_sim_hull[[1]]
          coords <- st_coordinates(hf)

          hull <- data$hull[[1]]

          if (!quiet) cat("    Writing data: ", out_path, "\n", sep = "")
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
        default = NULL
      )
    }

    cat(" -> Processing Lucazeau 2019 similarity interpolation by submap transect sets\n")
    walk(submap_transect_sets, f, ...)
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_submap_transect_set_kerswell2025_krg <- function(data_dir, out_dir, submap_transect_sets, ...) {
  with_error_handling({
    f <- function(submap_transect_set, ...) {
      with_error_handling(
        {
          dots <- list(...)
          quiet <- dots$quiet %||% FALSE

          save_dir <- file.path(out_dir, "heatflow", "submap_transects")
          create_dir(save_dir)

          out_path <- file.path(save_dir, paste0(submap_transect_set, "-kerswell2025-krg.csv"))

          if (file.exists(out_path)) {
            if (!quiet) cat(" -- Found data: ", out_path, "\n", sep = "")
            return(invisible())
          }

          opt_data <- get_optimal_krige_model(out_dir, submap_transect_set, ...)
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
            cache_path,
            data_dir,
            out_dir,
            submap_transect_set,
            variogram_model = variogram_model,
            fitted_variogram = fitted_variogram,
            max_n_point_pairs = max_n_point_pairs,
            max_point_pair_distance = max_point_pair_distance,
            ...
          )

          hf <- data[[krg_hull_str]][[1]]
          coords <- st_coordinates(hf)

          hull <- data$hull[[1]]

          if (!quiet) cat("    Writing data: ", out_path, "\n", sep = "")
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
        default = NULL
      )
    }

    cat(" -> Processing Kerswell 2025 krige interpolation by submap transect sets\n")
    walk(submap_transect_sets, f, ...)
  })
}
