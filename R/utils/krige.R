#######################################################
## Kriging with non-linear optimization              ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
krige_submap_transect_set <- function(sf_hf, data_col, fitted_variogram, sf_grid, max_n_point_pairs = 10, max_point_pair_distance = 1e5) {
  with_error_handling({
    cat("    Kriging: interpolating ", nrow(sf_hf), " observations onto ", length(sf_grid), " grid points\n", sep = "")
    proj <- determine_projection(sf_grid, proj_type = "aeqd")
    sf_hf_local <- st_transform(sf_hf, proj$wkt)
    sf_grid_local <- st_transform(sf_grid, proj$wkt)

    suppressWarnings({
      suppressMessages({
        krige_result_local <- krige(
          reformulate("1", response = data_col),
          sf_hf_local,
          newdata = sf_grid_local,
          model = fitted_variogram,
          nmax = max_n_point_pairs,
          maxdist = max_point_pair_distance,
          debug.level = 0
        )

        krige_result_ll <- st_transform(krige_result_local, crs = 4326)

        krige_result_ll |>
          as_tibble() |>
          st_as_sf(crs = 4326) |>
          mutate(
            kerswell2025_krg_est = ifelse(var1.pred > 0 & var1.pred <= 250, var1.pred, NA),
            kerswell2025_krg_sigma = ifelse(sqrt(var1.var) > 0 & sqrt(var1.var) <= 250, sqrt(var1.var), NA),
            .before = geometry
          ) |>
          rename(kerswell2025_geometry = geometry) |>
          select(-c(var1.pred, var1.var))
      })
    })
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
make_experimental_vgrm <- function(sf_hf, data_col, variogram_cutoff, n_variogram_lags) {
  with_error_handling({
    bbox <- st_bbox(sf_hf)
    bbox_diagonal_distance <- sqrt((bbox$xmax - bbox$xmin)^2 + (bbox$ymax - bbox$ymin)^2)
    lag_cutoff <- as.vector(bbox_diagonal_distance / variogram_cutoff)
    bin_width <- lag_cutoff / n_variogram_lags
    with_error_handling(variogram(reformulate("1", response = data_col), locations = sf_hf, cutoff = lag_cutoff, width = bin_width), default = NULL)
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit_variogram_model <- function(experimental_variogram, variogram_model, quiet = TRUE) {
  with_error_handling({
    if (nrow(experimental_variogram) < 2) {
      return(invisible())
    }

    fitted_variogram <- with_error_handling(
      {
        fit.variogram(experimental_variogram, vgm(model = variogram_model), fit.method = 6)
      },
      default = NULL
    )

    if (!is.null(fitted_variogram) && fitted_variogram$range < 0) {
      if (!quiet) cat(" !! Warning: variogram range is negative: ", fitted_variogram$range, "\n", sep = "")
      invisible()
    } else {
      fitted_variogram
    }
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
krige_cross_validation <- function(sf_hf, data_col, fitted_variogram, max_n_point_pairs, max_point_pair_distance, fold_vector) {
  with_error_handling({
    with_error_handling(
      {
        krige.cv(
          reformulate("1", response = data_col),
          sf_hf,
          model = fitted_variogram,
          nmax = max_n_point_pairs,
          maxdist = max_point_pair_distance,
          nfold = fold_vector
        )
      },
      default = NULL
    )
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log_data <- function(path, data) {
  with_error_handling({
    create_path_dir(path)
    sink(path, append = TRUE)
    cat(paste(names(data), data, sep = ": ", collapse = "\n"), "\n\n", sep = "")
    sink()
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sample_and_shuffle_equal_nfolds <- function(n, ...) {
  with_error_handling({
    dots <- list(...)
    seed <- dots$seed %||% 42
    n_folds <- dots$n_folds %||% 10

    set.seed(seed)
    fold_size <- n %/% n_folds
    remainder <- n %% n_folds
    fold <- rep(1:n_folds, each = fold_size)
    if (remainder > 0) fold <- c(fold, sample(1:n_folds, remainder))
    sample(fold)
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cost_function <- function(
  sf_hf,
  out_dir,
  variogram_cutoff = 1,
  n_variogram_lags = 30,
  max_n_point_pairs = 10,
  max_point_pair_distance = 10,
  variogram_model = "Sph",
  interpolation_weight = 0.5,
  variogram_weight = 0.5,
  submap_transect_set = NULL,
  fold_vector = NULL,
  na_threshold_proportion = 0.5
) {
  with_error_handling(
    {
      penalty_cost <- 1e1

      experimental_variogram <- make_experimental_vgrm(sf_hf, "ihfc2024_obs", variogram_cutoff, n_variogram_lags)
      if (is.null(experimental_variogram)) {
        return(list(penalty_cost, penalty_cost, penalty_cost))
      }

      fitted_variogram <- fit_variogram_model(experimental_variogram, variogram_model)
      if (is.null(fitted_variogram)) {
        return(list(penalty_cost, penalty_cost, penalty_cost))
      }

      krige_cv <- krige_cross_validation(sf_hf, "ihfc2024_obs", fitted_variogram, max_n_point_pairs, max_point_pair_distance, fold_vector)
      if (is.null(krige_cv)) {
        return(list(penalty_cost, penalty_cost, penalty_cost))
      }

      na_sum <- sum(is.na(krige_cv$residual))
      na_thresh <- floor(nrow(krige_cv) * na_threshold_proportion)

      if (na_sum >= na_thresh) {
        return(list(penalty_cost, penalty_cost, penalty_cost))
      }

      variogram_rmse <- sqrt(sqrt(attr(fitted_variogram, "SSErr") / nrow(experimental_variogram)))
      variogram_sigma <- sqrt(sd(experimental_variogram$gamma, na.rm = TRUE))
      if (variogram_sigma == 0) {
        return(list(penalty_cost, penalty_cost, penalty_cost))
      }

      interpolation_rmse <- sqrt(sum(krige_cv$residual^2, na.rm = TRUE) / nrow(krige_cv))
      interpolation_sd <- sd(krige_cv$var1.pred, na.rm = TRUE)
      if (interpolation_sd == 0) {
        return(list(penalty_cost, penalty_cost, penalty_cost))
      }

      variogram_cost <- variogram_weight * variogram_rmse / variogram_sigma
      interpolation_cost <- interpolation_weight * interpolation_rmse / interpolation_sd

      total_cost <- variogram_cost + interpolation_cost

      if (variogram_sigma == 0 || interpolation_sd == 0) {
        return(list(penalty_cost, penalty_cost, penalty_cost))
      }

      nlopt_dir <- file.path(out_dir, "nlopt")
      iter_dir <- file.path(nlopt_dir, "iterations")
      out_path <- file.path(iter_dir, paste0(submap_transect_set, "-", variogram_model, "-run"))

      data <-
        c(
          "Submap transect set" = submap_transect_set,
          "No. heatflow observations" = nrow(sf_hf),
          "Variogram variogram_cutoff" = variogram_cutoff,
          "No. variogram lags" = n_variogram_lags,
          "Max no. point pairs" = max_n_point_pairs,
          "Max point pair distance" = max_point_pair_distance,
          "Variogram model" = variogram_model,
          "Variogram weight" = variogram_weight,
          "Variogram rmse" = variogram_rmse,
          "Variogram uncertainty" = variogram_sigma,
          "Variogram cost" = variogram_cost,
          "Interpolation weight" = interpolation_weight,
          "Interpolation rmse" = interpolation_rmse,
          "Interpolation sd" = interpolation_sd,
          "Interpolation cost" = interpolation_cost,
          "Total cost" = total_cost
        )

      log_data(out_path, data)
      list(variogram_cost, interpolation_cost, total_cost)
    },
    default = list(1e9, 1e9, 1e9)
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
translate_nlopt_results_to_variogram_model <- function(sf_hf, data_col, variogram_model, nlopt_object) {
  with_error_handling({
    suppressWarnings({
      experimental_variogram <- make_experimental_vgrm(sf_hf, data_col, nlopt_object$solution[1], nlopt_object$solution[2])

      if (is.null(experimental_variogram)) {
        return(invisible())
      }

      fitted_variogram <- fit_variogram_model(experimental_variogram, variogram_model)

      if (is.null(fitted_variogram)) {
        return(invisible())
      }
    })
    list("experimental_variogram" = as_tibble(experimental_variogram), "fitted_variogram" = fitted_variogram)
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
optimize_krige_model <- function(data_dir, out_dir, submap_transect_set, variogram_model, ...) {
  with_error_handling({
    dots <- list(...)
    n_starts <- dots$n_starts %||% 100
    n_top <- dots$n_top %||% 3
    n_folds <- dots$n_folds %||% 5
    data_col <- dots$data_col %||% "ihfc2024_obs"
    interpolation_weight <- dots$interpolation_weight %||% 0.5
    variogram_weight <- dots$variogram_weight %||% 0.5
    na_threshold_proportion <- dots$na_threshold_proportion %||% 0.5
    initial_guess <- dots$initial_guess %||% c(3, 30, 10, 20) # cutoff, lags, n pairs, pair distance (m) / 1e4
    lower_bounds <- dots$lower_bounds %||% c(3, 20, 2, 10)
    upper_bounds <- dots$upper_bounds %||% c(12, 50, 50, 50)
    nlopt_algorithm_global <- dots$nlopt_algorithm_global %||% "NLOPT_GN_MLSL_LDS"
    max_nlopt_evaluations_global <- dots$max_nlopt_evaluations_global %||% 100
    x_tolerance_global <- dots$x_tolerance_global %||% 1e-4 # parameter-space tolerance for global search
    f_tolerance_global <- dots$f_tolerance_global %||% 1e-4 # cost function tolerance for global search
    nlopt_algorithm_local <- dots$nlopt_algorithm_local %||% "NLOPT_LN_NELDERMEAD"
    max_nlopt_evaluations_local <- dots$max_nlopt_evaluations_local %||% 500
    x_tolerance_local <- dots$x_tolerance_local %||% 1e-5
    f_tolerance_local <- dots$f_tolerance_local %||% 1e-5
    nlopt_algorithm_fine_tune <- dots$nlopt_algorithm_fine_tune %||% "NLOPT_LN_NELDERMEAD"
    max_nlopt_evaluations_fine_tune <- dots$max_nlopt_evaluations_fine_tune %||% 500
    x_tolerance_fine_tune <- dots$x_tolerance_fine_tune %||% 1e-8
    f_tolerance_fine_tune <- dots$f_tolerance_fine_tune %||% 1e-8
    quiet <- dots$quiet %||% FALSE

    dirs <- list(
      iter = file.path(out_dir, "nlopt", "iterations"),
      krige = file.path(out_dir, "nlopt", "krige_models"),
      cache = file.path(out_dir, "submap")
    )

    nlopt_id <- paste0(submap_transect_set, "-", variogram_model)
    nlopt_path <- file.path(dirs$krige, paste0(nlopt_id, ".RData"))
    nlopt_iter_path_run <- file.path(dirs$iter, paste0(submap_transect_set, "-", variogram_model, "-run"))
    nlopt_iter_path_fail <- file.path(dirs$iter, paste0(submap_transect_set, "-", variogram_model, "-fail"))

    if (file.exists(nlopt_path) && file.exists(nlopt_iter_path_run)) {
      if (!quiet) {
        cat(" -- Found nlopt data:\n")
        cat("    ", nlopt_path, "\n", sep = "")
        cat("    ", nlopt_iter_path_run, "\n", sep = "")
      }
      return(invisible())
    }

    if (!file.exists(nlopt_path) && file.exists(nlopt_iter_path_run)) {
      file.remove(nlopt_iter_path_run)
    }

    if (file.exists(nlopt_iter_path_fail)) {
      if (!quiet) cat(" !! Warning: optimized ", variogram_model, " kriging model failed for: ", submap_transect_set, "\n", sep = "")
      return(invisible())
    }

    create_dir(dirs$krige)

    if (nlopt_algorithm_global %in% c("NLOPT_GN_MLSL_LDS", "NLOPT_GN_MLSL")) {
      nlopt_options <- list(
        print_level = 0,
        maxeval = max_nlopt_evaluations_global,
        algorithm = nlopt_algorithm_global,
        local_opts = list(
          print_level = 0,
          maxeval = max_nlopt_evaluations_local,
          algorithm = nlopt_algorithm_local,
          xtol_rel = x_tolerance_local,
          ftol_rel = f_tolerance_local
        )
      )
    } else {
      nlopt_options <- list(
        print_level = 0,
        maxeval = max_nlopt_evaluations_global,
        algorithm = nlopt_algorithm_global,
        xtol_rel = x_tolerance_global,
        ftol_rel = f_tolerance_global
      )
    }

    nlopt_options_fine_tune <- list(
      print_level = 0,
      maxeval = max_nlopt_evaluations_fine_tune,
      algorithm = nlopt_algorithm_fine_tune,
      xtol_rel = x_tolerance_fine_tune,
      ftol_rel = f_tolerance_fine_tune
    )

    cache_path <- file.path(dirs$cache, paste0(submap_transect_set, ".RData"))
    data <- compile_submap_transect_set_data(cache_path, data_dir, out_dir, submap_transect_set, ...)

    proj <- determine_projection(data$hull, proj_type = "aeqd")
    sf_hf <- st_transform(data$ihfc2024_obs_hull[[1]], proj$wkt)
    fold_vector <- sample_and_shuffle_equal_nfolds(nrow(sf_hf), n_folds = n_folds, ...)

    nlopt_fun <- function(x) {
      itr <<- itr + 1

      current_cost <- unlist(cost_function(
        sf_hf = sf_hf,
        out_dir = out_dir,
        variogram_cutoff = x[1],
        n_variogram_lags = x[2],
        max_n_point_pairs = x[3],
        max_point_pair_distance = x[4] * 1e4,
        variogram_model = variogram_model,
        interpolation_weight = interpolation_weight,
        variogram_weight = variogram_weight,
        submap_transect_set = submap_transect_set,
        fold_vector = fold_vector,
        na_threshold_proportion = na_threshold_proportion
      ))

      if (itr > 1) {
        delta_x <- sqrt(sum((x - prev_x)^2))
        delta_f <- abs(current_cost[3] - prev_cost[3])
        prev_x <<- x
        prev_cost <<- current_cost
      } else {
        delta_x <- NA
        delta_f <- NA
        prev_x <<- x
        prev_cost <<- current_cost
      }

      if (!quiet) {
        cat("\033[7A\033[J")
        format_str <- paste(
          "    Iteration                        : %.0f\n",
          "    Parameters                       : cutoff=%.6f, lags=%.6f, n pairs=%.6f, pair dist=%.6f\n",
          "    Variogram cost                   : %.6f\n",
          "    Interpolation cost               : %.6f\n",
          "    Total cost                       : %.6f\n",
          "    Delta X (moving to)              : %.2e\n",
          "    Delta F (change in cost)         : %.2e\n",
          sep = ""
        )
        cat(sprintf(format_str, itr, x[1], x[2], x[3], x[4], current_cost[1], current_cost[2], current_cost[3], delta_x, delta_f))
        flush.console()
      }

      current_cost[3]
    }

    itr <- 0
    prev_x <- NULL
    prev_f <- NULL
    exploration <- explore_parameter_space(n_starts, n_top, lower_bounds, upper_bounds, nlopt_fun, quiet = quiet)
    best_starts <- head(exploration, n_top)

    if (!quiet) cat(" -> Refining best ", n_top, " starting points with nloptr\n", sep = "")

    opt_results <- best_starts |>
      split(f = seq_len(nrow(best_starts))) |>
      imap(function(row, i) {
        start_vec <- as.numeric(row[1:4])
        with_error_handling(
          {
            local <- nlopt_algorithm_global %in% c("NLOPT_GN_MLSL_LDS", "NLOPT_GN_MLSL")
            format_str <- "| %16s | %16s | %17s | %23s |"
            init_str <- as.character(start_vec)
            lb_str <- as.character(lower_bounds)
            ub_str <- as.character(upper_bounds)
            init_str_formatted <- sprintf(format_str, init_str[1], init_str[2], init_str[3], init_str[4])
            lb_str_formatted <- sprintf(format_str, lb_str[1], lb_str[2], lb_str[3], lb_str[4])
            ub_str_formatted <- sprintf(format_str, ub_str[1], ub_str[2], ub_str[3], ub_str[4])

            if (!quiet) {
              cat(" ", i, ". Optimizing ", variogram_model, " kriging model for : ", submap_transect_set, "\n", sep = "")
              if (i == 1) cat("    n observations                   : ", nrow(sf_hf), "\n", sep = "")
              if (i == 1) cat("    n folds                          : ", n_folds, "\n", sep = "")
              if (i == 1) cat("    Interpolation weight             : ", interpolation_weight, "\n", sep = "")
              if (i == 1) cat("    Variogram weight                 : ", variogram_weight, "\n", sep = "")
              if (i == 1) cat("    nlopt algorithm       (global)   : ", nlopt_algorithm_global, "\n", sep = "")
              if (i == 1) cat("    Maximum evaluations   (global)   : ", max_nlopt_evaluations_global, "\n", sep = "")
              if (i == 1 && !local) cat("    X tolerance           (global)   : ", x_tolerance_global, " (parameter-space tolerance)\n", sep = "")
              if (i == 1 && !local) cat("    F tolerance           (global)   : ", f_tolerance_global, " (cost function tolerance)\n", sep = "")
              if (i == 1 && local) cat("    nlopt algorithm       (local)    : ", nlopt_algorithm_local, "\n", sep = "")
              if (i == 1 && local) cat("    Maximum evaluations   (local)    : ", max_nlopt_evaluations_local, "\n", sep = "")
              if (i == 1 && local) cat("    X tolerance           (local)    : ", x_tolerance_local, "\n", sep = "")
              if (i == 1 && local) cat("    F tolerance           (local)    : ", f_tolerance_local, "\n", sep = "")
              cat("\n")
              cat("    |------------------|------------------|-------------------|-------------------------|\n", sep = "")
              cat("    | variogram_cutoff | n_variogram_lags | max_n_point_pairs | max_point_pair_distance |\n", sep = "")
              cat("    |------------------|------------------|-------------------|-------------------------|\n", sep = "")
              cat("    ", init_str_formatted, " Initial parameters\n", sep = "")
              cat("    ", lb_str_formatted, " Lower bounds\n", sep = "")
              cat("    ", ub_str_formatted, " Upper bounds\n", sep = "")
              cat("    |------------------|------------------|-------------------|-------------------------|\n", sep = "")
              cat(rep("\n", 8))
            }

            itr <<- -2
            prev_x <<- NULL
            prev_f <<- NULL
            opt <- nloptr(start_vec, nlopt_fun, lb = lower_bounds, ub = upper_bounds, opts = nlopt_options)

            if (is.null(opt) || opt$status < 0 || opt$iterations < 10) {
              if (!quiet) {
                cat("    --------------------------------------------------\n")
                cat(" !! Warning: nlopt failed to converge:\n")
                cat("    nlopt status: ", opt$status, "\n", sep = "")
                cat("    nlopt iterations: ", opt$iterations, "\n", sep = "")
                cat("    ", opt$message, "\n", sep = "")
                cat("    --------------------------------------------------\n")
              }
              return(opt)
            } else {
              if (all(opt$objective >= 1e8)) {
                if (!quiet) {
                  cat("    --------------------------------------------------\n")
                  cat(" !! Warning: nlopt only produced penalties (no valid models)\n")
                  cat("    --------------------------------------------------\n")
                }
                return(invisible())
              }

              if (!quiet) {
                cat("    --------------------------------------------------\n")
                cat("    Nlopt converged with status '", opt$status, "' after ", opt$iterations, " iterations:\n", sep = "")
                cat("    Nlopt message: ", opt$message, "\n", sep = "")
                cat("    --------------------------------------------------\n")
              }
              return(opt)
            }
          },
          default = NULL
        )
      })

    valid_results <- compact(opt_results) |> keep(~ is.finite(.x$objective))

    if (length(valid_results) == 0) {
      if (!quiet) {
        cat("    --------------------------------------------------\n")
        cat(" !! All refinements failed to converge\n")
        cat("    --------------------------------------------------\n")
      }
      return(invisible())
    }

    best_opt <- valid_results |> reduce(function(best, current) if (current$objective < best$objective) current else best)

    if (!quiet) cat(" -> Fine tune ", variogram_model, " kriging model for  : ", submap_transect_set, rep("\n", 8), sep = "")
    shrink_factor <- 0.2
    lower_bounds_fine <- pmax(lower_bounds, best_opt$solution - shrink_factor * (upper_bounds - lower_bounds))
    upper_bounds_fine <- pmin(upper_bounds, best_opt$solution + shrink_factor * (upper_bounds - lower_bounds))

    itr <- -2
    prev_x <- NULL
    prev_f <- NULL
    fine_tune_opt <- nloptr(best_opt$solution, nlopt_fun, lb = lower_bounds_fine, ub = upper_bounds_fine, opts = nlopt_options_fine_tune)

    if (is.null(fine_tune_opt) || fine_tune_opt$status < 0 || fine_tune_opt$iterations < 10) {
      if (!quiet) {
        cat("    --------------------------------------------------\n")
        cat(" !! Warning: nlopt failed to converge:\n")
        cat("    nlopt status: ", fine_tune_opt$status, "\n", sep = "")
        cat("    nlopt iterations: ", fine_tune_opt$iterations, "\n", sep = "")
        cat("    ", fine_tune_opt$message, "\n", sep = "")
        cat("    --------------------------------------------------\n")
      }
    } else {
      if (all(fine_tune_opt$objective >= 1e8)) {
        if (!quiet) {
          cat("    --------------------------------------------------\n")
          cat(" !! Warning: nlopt only produced penalties (no valid models)\n")
          cat("    --------------------------------------------------\n")
        }
      }

      if (!quiet) {
        cat("    --------------------------------------------------\n")
        cat("    Nlopt converged with status '", fine_tune_opt$status, "' after ", fine_tune_opt$iterations, " iterations:\n", sep = "")
        cat("    Nlopt message: ", fine_tune_opt$message, "\n", sep = "")
        cat("    --------------------------------------------------\n")
      }
    }

    if (is.null(fine_tune_opt) || fine_tune_opt$status < 0 || fine_tune_opt$iterations < 10) {
      if (!is.null(fine_tune_opt)) {
        data <- list(
          nlopt_status = fine_tune_opt$status,
          nlopt_iterations = fine_tune_opt$iterations,
          message = fine_tune_opt$message,
          info = "nlopt failed to converge"
        )
        log_data(nlopt_iter_path_fail, data)
      } else {
        data <- list(
          nlopt_status = NULL,
          nlopt_iterations = fine_tune_opt$iterations,
          message = fine_tune_opt$message,
          info = "nlopt failed to converge"
        )
        log_data(nlopt_iter_path_fail, data)
      }
      return(invisible())
    } else {
      if (all(fine_tune_opt$objective >= 1e8)) {
        data <- list(
          nlopt_status = fine_tune_opt$status,
          nlopt_iterations = fine_tune_opt$iterations,
          message = fine_tune_opt$message,
          info = "nlopt only produced penalties"
        )
        log_data(nlopt_iter_path_fail, data)
        return(invisible())
      }

      nlopt_out <- translate_nlopt_results_to_variogram_model(sf_hf, data_col, variogram_model, fine_tune_opt)
      assign(str_replace_all(nlopt_id, "-", "_"), nlopt_out)
      save(list = str_replace_all(nlopt_id, "-", "_"), file = nlopt_path)
    }
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
explore_parameter_space <- function(n_starts, n_top, lower, upper, cost_fun, ...) {
  with_error_handling(
    {
      dots <- list(...)
      seed <- dots$seed %||% 42
      quiet <- dots$quiet %||% FALSE

      starts <- randomLHS(n_starts, length(lower))
      starts <- t(t(starts) * (upper - lower) + lower)
      param_names <- paste0("x", seq_along(lower))
      starts <- as_tibble(starts, .name_repair = ~param_names)

      if (!quiet) cat(" -> Exploring parameter space from ", n_starts, " random Latin Hypercube Samples (LHS)", rep("\n", 8), sep = "")

      eval_fun <- function(x) tryCatch(cost_fun(as.numeric(x)), error = function(e) Inf)
      costs <- map_dbl(seq_len(n_starts), ~ eval_fun(starts[.x, ]))
      df <- tibble(starts, cost = costs)
      df <- df[order(df$cost), ]

      if (!quiet) {
        cat(" -> Best exploratory candidates :\n\n", sep = "")
        format_str <- "    | %16s | %16s | %17s | %23s | %10s |\n"
        cat("    |------------------|------------------|-------------------|-------------------------|-------------------|\n")
        cat("    | variogram_cutoff | n_variogram_lags | max_n_point_pairs | max_point_pair_distance | total cost        |\n")
        cat("    |------------------|------------------|-------------------|-------------------------|-------------------|\n")
        pmap(head(df, n_top), function(x1, x2, x3, x4, cost) {
          cat(sprintf(format_str, x1, x2, x3, x4, cost), sep = "")
        })
        cat("    |------------------|------------------|-------------------|-------------------------|-------------------|\n\n")
      }
      df
    },
    default = NULL
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
process_submap_transect_sets <- function(data_dir, out_dir, submap_transect_sets, ...) {
  with_error_handling({
    dots <- list(...)
    variogram_models <- dots$variogram_models %||% c("Sph", "Exp", "Gau")
    seed <- dots$seed %||% 42
    nprocs <- dots$nprocs %||% availableCores() - 2
    parallel <- dots$parallel %||% TRUE
    quiet <- dots$quiet %||% FALSE

    relief_dir <- file.path(out_dir, "relief")
    data <- expand.grid(set = submap_transect_sets, model = variogram_models, stringsAsFactors = FALSE) |> arrange(set, model)
    f_partial <- partial(optimize_krige_model, data_dir = data_dir, out_dir = out_dir, !!!dots)

    if (parallel) {
      set.seed(seed)
      plan(multisession, workers = nprocs)
      future_walk2(data$set, data$model, f_partial, .options = furrr_options(seed = seed))
    } else {
      walk2(data$set, data$model, f_partial)
    }
  })
}
