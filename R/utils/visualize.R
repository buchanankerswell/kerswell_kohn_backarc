#######################################################
## Helper functions (visualizations)                 ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
viridis_scale_color <- function() {
  list(
    scale_color_viridis_c(
      option = "magma",
      name = bquote("Q" ~ (mWm^-2)),
      limits = c(0, 250),
      breaks = c(0, 125, 250),
      na.value = "transparent",
      guide = guide_colorbar(title.vjust = 1, frame.colour = "black", ticks.colour = "black")
    )
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
viridis_scale_fill <- function() {
  list(
    scale_fill_viridis_c(
      option = "magma",
      name = bquote("Q" ~ (mWm^-2)),
      limits = c(0, 250),
      breaks = c(0, 125, 250),
      na.value = "transparent",
      guide = guide_colorbar(title.vjust = 1, frame.colour = "black", ticks.colour = "black")
    )
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
elevation_scale_fill <- function(palette = "meyers", limits = NULL, breaks = NULL, alpha = 0.5) {
  list(
    scale_fill_hypso_tint_c(
      name = "Elevation (m)",
      palette = palette,
      alpha = alpha,
      labels = label_number(),
      breaks = breaks,
      limits = limits,
      oob = oob_squish,
      guide = guide_colorbar(title.vjust = 1, frame.colour = "black", ticks.colour = "black", theme = theme(legend.margin = margin(7, 0, 0, 0)))
    )
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
theme_facet <- function(base_size = 14, show_legend = FALSE) {
  theme <- list(
    theme_bw(base_size = base_size),
    theme(
      panel.grid.major = element_line(linewidth = 0.4, color = "white"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "grey90"),
      panel.spacing = unit(1, "lines"),
      plot.margin = margin(5, 15, 5, 5),
      plot.title = element_text(hjust = 0.5),
      axis.ticks = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size),
      legend.justification = "center",
      legend.position = "bottom",
      legend.box.margin = margin(2, 2, 2, 2),
      legend.margin = margin(),
      legend.title = element_text(vjust = 0.5),
      legend.title.position = "top",
      legend.background = element_blank(),
      legend.text = element_text(size = base_size * 0.8)
    )
  )

  if (!show_legend) {
    theme <- c(theme, list(theme(legend.position = "none")))
  }

  theme
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
theme_vgrm <- function(base_size = 14, show_legend = FALSE) {
  theme <- list(
    theme_bw(base_size = base_size),
    theme(
      panel.grid.major = element_line(linewidth = 0.4, color = "white"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "grey90"),
      panel.spacing = unit(1, "lines"),
      plot.margin = margin(5, 15, 5, 5),
      plot.title = element_text(size = base_size * 1.2, hjust = 0, margin = margin(10, 0, -10, 60)),
      axis.ticks = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size),
      legend.justification = "right",
      legend.position = "top",
      legend.key.height = unit(0.7, "lines"),
      legend.key.width = unit(3.0, "lines"),
      legend.position.inside = c(0.98, 0.75),
      legend.box.margin = margin(2, 2, 2, 2),
      legend.margin = margin(-20, 0, 0, 0),
      legend.title = element_text(hjust = 1, margin = margin(0, 5, 0, 0)),
      legend.title.position = "left",
      legend.background = element_blank(),
      legend.text = element_text(size = base_size * 0.8)
    )
  )

  if (!show_legend) {
    theme <- c(theme, list(theme(legend.position = "none")))
  }

  theme
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
theme_loess <- function(base_size = 14, show_legend = TRUE) {
  theme <- list(
    scale_y_continuous(limits = c(0, 250), breaks = seq(0, 250, 50)),
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)),
    theme_bw(base_size = base_size),
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "grey90"),
      plot.margin = margin(5, 5, 5, 5),
      legend.justification = "right",
      legend.position = "inside",
      legend.position.inside = c(0.92, 0.85),
      legend.direction = "horizontal",
      legend.key.height = unit(0.7, "lines"),
      legend.key.width = unit(1.8, "lines"),
      legend.box.margin = margin(2, 2, 2, 2),
      legend.margin = margin(),
      legend.title.position = "top",
      legend.title = element_text(hjust = 0.5, margin = margin()),
      legend.background = element_blank()
    )
  )

  if (!show_legend) {
    theme <- c(theme, list(theme(legend.position = "none")))
  }

  theme
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
theme_globe <- function(proj = NULL, base_size = 18, show_legend = TRUE) {
  if (!is.null(proj)) wkt <- proj$wkt else st_crs("4326")

  theme <- list(
    coord_sf(crs = wkt, expand = FALSE),
    theme_map(font_size = base_size),
    theme(
      axis.text = element_blank(),
      plot.margin = margin(5, 5, 5, 5),
      legend.position = "bottom",
      legend.justification = "left",
      legend.direction = "horizontal",
      legend.margin = margin(0, 0, 10, 0),
      legend.key.height = unit(0.7, "lines"),
      legend.key.width = unit(2, "lines"),
      legend.title = element_text(hjust = 0.5, margin = margin()),
      legend.title.position = "top"
    )
  )

  if (!show_legend) {
    theme <- c(theme, list(theme(legend.position = "none")))
  }

  theme
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
theme_transect_set <- function(sf_obj, proj = NULL, base_size = 18, show_legend = TRUE) {
  if (!is.null(proj)) wkt <- proj$wkt else st_crs("4326")
  if (!is.null(proj)) bbox <- st_bbox(st_transform(sf_obj, wkt)) else st_bbox(sf_obj)

  theme <- list(
    annotation_scale(location = "bl", width_hint = 0.33, text_cex = 1.3, style = "bar", line_width = 3.5),
    annotation_north_arrow(
      location = "bl",
      which_north = "true",
      pad_x = unit(0, "lines"),
      height = unit(4, "lines"),
      width = unit(4, "lines"),
      pad_y = unit(1, "lines"),
      style = north_arrow_fancy_orienteering
    ),
    scale_x_continuous(breaks = seq(-180, 180, by = 15)),
    # scale_y_continuous(breaks = seq(-90, 90, by = 10)),
    coord_sf(crs = wkt, xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE),
    theme_minimal(base_size = base_size),
    theme(
      panel.grid.major = element_line(linewidth = 0.2, color = "grey70"),
      plot.margin = margin(5, 5, 5, 5),
      plot.title = element_text(size = base_size * 1.2, hjust = 0.5, margin = margin(0, 0, 5, 0)),
      axis.title = element_blank()
    )
  )

  if (!show_legend) {
    theme <- c(theme, list(theme(legend.position = "none")))
  }

  theme
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
basemap_layers <- function(data_dir, out_dir, proj = NULL, major_grats = TRUE, minor_grats = TRUE) {
  with_error_handling(
    {
      wrap_opts <- c("WRAPDATELINE=YES", "DATELINEOFFSET=180")
      if (!is.null(proj)) lon0 <- proj$lon0 else 0

      etopo_dir <- file.path(data_dir, "mapping", "relief")
      geotif_dir <- file.path(out_dir, "relief")

      sr_relief <- with_error_handling(get_noaa_relief(-180, 180, 90, -90, etopo_dir, geotif_dir), default = NULL)
      sf_coast <- with_error_handling(st_break_antimeridian(st_wrap_dateline(sf_coast, options = wrap_opts), lon0), default = NULL)
      sf_ridge <- with_error_handling(st_break_antimeridian(st_wrap_dateline(sf_ridge, options = wrap_opts), lon0), default = NULL)
      sf_transform <- with_error_handling(st_break_antimeridian(st_wrap_dateline(sf_transform, options = wrap_opts), lon0), default = NULL)
      sf_trench <- with_error_handling(st_break_antimeridian(st_wrap_dateline(sf_trench, options = wrap_opts), lon0), default = NULL)
      sf_graticules_minor <- with_error_handling(
        st_break_antimeridian(st_wrap_dateline(st_graticule(lat = seq(-90, 90, by = 5), lon = seq(-180, 180, by = 5)), options = wrap_opts), lon0),
        default = NULL
      )
      sf_graticules_major <- with_error_handling(
        st_break_antimeridian(st_wrap_dateline(st_graticule(lat = seq(-90, 90, by = 15), lon = seq(-180, 180, by = 15)), options = wrap_opts), lon0),
        default = NULL
      )

      if (!minor_grats) sf_graticules_minor <- NULL
      if (!major_grats) sf_graticules_major <- NULL

      relief_range <- minmax(sr_relief)

      list(
        geom_spatraster(data = sr_relief),
        elevation_scale_fill(),
        new_scale_fill(),
        geom_sf(data = sf_graticules_minor, color = "grey50", linewidth = 0.1),
        geom_sf(data = sf_graticules_major, color = "black", linewidth = 0.1),
        geom_sf(data = sf_coast, linewidth = 0.2, color = "black"),
        geom_sf(data = sf_ridge, linewidth = 0.6, color = "white"),
        geom_sf(data = sf_transform, linewidth = 0.6, color = "white"),
        geom_sf(data = sf_trench, linewidth = 0.6, color = "white")
      )
    },
    default = NULL
  )
}

#######################################################
## Plotting functions                                ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_global_dataset_composition <- function(data_dir, out_dir, fig_dir, base_size = 14) {
  out_path <- file.path(fig_dir, paste0("global-dataset-composition.png"))

  if (check_plot_path(out_path)) {
    return(invisible())
  }

  blank_plot_list <- list(
    transect = ggplot() +
      theme_void(),
    heatflow = ggplot() +
      theme_void()
  )

  plots <- map(c(120, -180, -90), ~ {
    with_error_handling(
      {
        proj <- NULL
        proj$wkt <- paste0("+proj=ortho +lon_0=", .x, " +lat_0=0 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs")
        proj$lon0 <- .x

        p1 <- ggplot() +
          basemap_layers(data_dir, out_dir, proj, minor_grats = FALSE) +
          geom_sf(data = sf_hull, fill = "black", color = "black", linewidth = 0.3, alpha = 0.2, show.legend = FALSE) +
          geom_sf(data = sf_submap, color = "black", linewidth = 0.6) +
          theme_globe(proj, base_size)

        p2 <- ggplot() +
          basemap_layers(data_dir, out_dir, proj, minor_grats = FALSE) +
          geom_sf(data = sf_ihfc_raw, aes(color = ihfc2024_obs), size = 0.2, shape = 20) +
          viridis_scale_color() +
          theme_globe(proj, base_size)

        list(transect = p1, heatflow = p2)
      },
      default = blank_plot_list
    )
  })

  suppressWarnings({
    row1 <- wrap_plots(lapply(plots, `[[`, "transect"), nrow = 1)
    row2 <- wrap_plots(lapply(plots, `[[`, "heatflow"), nrow = 1)
    p <- row1 / row2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
    ggsave(file = out_path, plot = p, width = 6.5, height = 5.5, dpi = 300, bg = "white")
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_nlopt_summary <- function(out_dir, fig_dir, base_size = 16) {
  with_error_handling({
    nlopt_data <- file.path(out_dir, "nlopt-data.RData")
    out_path <- file.path(fig_dir, "nlopt.png")

    if (check_plot_path(out_path)) {
      return(invisible())
    }

    suppressWarnings({
      suppressMessages({
        p <- nlopt_summary |>
          mutate(max_point_pair_distance = max_point_pair_distance / 1e3) |>
          select(-c(
            variogram_weight,
            variogram_cost,
            variogram_rmse,
            variogram_sigma,
            interpolation_weight,
            interpolation_cost,
            interpolation_rmse,
            interpolation_sigma
          )) |>
          rename(
            "Max Pt. Pair Dist." = max_point_pair_distance,
            "V'gram Cutoff" = variogram_cutoff,
            "Nlopt Iter's" = itr,
            "No. V'gram Lags" = n_variogram_lags,
            "No. HF Obs." = n_heatflow_obs,
            "Max No. Pt. Pairs" = max_n_point_pairs
          ) |>
          pivot_longer(-c(submap_transect_set, variogram_model, total_cost)) |>
          ggplot(aes(value, total_cost, fill = submap_transect_set)) +
          geom_point(size = 2.5, shape = 21, color = "black") +
          scale_fill_d3(palette = "category20") +
          labs(x = NULL, y = "Cost", fill = "Submap Transect Set") +
          scale_x_continuous(breaks = pretty_breaks(n = 4), guide = guide_axis(check.overlap = TRUE), expand = expansion(mult = 0.10)) +
          guides(fill = guide_legend(ncol = 4)) +
          facet_wrap(~name, scales = "free_x") +
          theme_facet(base_size, show_legend = TRUE)

        ggsave(file = out_path, plot = p, width = 6.5, height = 7.0, dpi = 300, bg = "white")
      })
    })
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_interpolation_accuracy_summary <- function(out_dir, fig_dir, base_size = 16) {
  with_error_handling({
    map_data <- file.path(out_dir, "map-data.RData")
    nlopt_data <- file.path(out_dir, "nlopt-data.RData")
    out_path <- file.path(fig_dir, "sim-krg-accuracies.png")

    if (check_plot_path(out_path)) {
      return(invisible())
    }

    suppressWarnings({
      suppressMessages({
        df <- interp_accuracy_summary |>
          filter(!is.na(rmse_krige)) |>
          pivot_longer(-c(submap_transect_set)) |>
          mutate(method = ifelse(str_detect(name, "sim"), "sim", "krg"), name = str_split(name, "_", simplify = TRUE)[, 1]) |>
          rename(metric = name) |>
          mutate(method = ifelse(method == "sim", "Similarity", "Krige")) |>
          filter(metric == "rmse") |>
          rename(rmse = value) |>
          select(-metric) |>
          left_join(select(nlopt_summary, submap_transect_set, n_heatflow_obs), by = "submap_transect_set") |>
          left_join(select(interp_accuracy_summary, submap_transect_set, n_krige), by = "submap_transect_set") |>
          rename(n_control = n_krige) |>
          left_join(select(interp_diff_summary, submap_transect_set, n_grid), by = "submap_transect_set") |>
          left_join(select(interp_diff_summary, submap_transect_set, n), by = "submap_transect_set") |>
          mutate(n = ifelse(method == "Similarity", n_grid, n)) |>
          rename(n_itp = n) |>
          mutate(coverage = n_control / n_grid * 100)

        p <- df |>
          ggplot() +
          geom_point(aes(coverage, rmse, fill = submap_transect_set), size = 2.5, shape = 21, color = "black") +
          scale_fill_d3(palette = "category20") +
          labs(x = "Observational Coverage (%)", y = "RMSE", fill = "Submap Transect Set") +
          guides(fill = guide_legend(ncol = 4)) +
          facet_wrap(~method, scales = "fixed") +
          theme_facet(base_size, show_legend = TRUE)

        ggsave(file = out_path, plot = p, width = 6.5, height = 5.5, dpi = 300, bg = "white")
      })
    })
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_optimal_krige_model <- function(data_dir, out_dir, fig_dir, submap_transect_set, ..., base_size = 20) {
  with_error_handling({
    suppressWarnings({
      suppressMessages({
        out_path <- file.path(fig_dir, "submap", paste0(submap_transect_set, "-vgrm.png"))

        if (check_plot_path(out_path)) {
          return(invisible())
        }

        opt_data <- get_optimal_krige_model(out_dir, submap_transect_set, ...)
        opt_summary <- opt_data$optimal_krige_model
        nlopt_iter_df <- opt_data$nlopt_iter_df
        experimental_variogram <- opt_data$experimental_variogram
        fitted_variogram <- opt_data$fitted_variogram
        fitted_variogram_line <- variogramLine(fitted_variogram, maxdist = max(experimental_variogram$dist))
        opt_variogram_model <- opt_summary$variogram_model
        max_n_point_pairs <- opt_summary$max_n_point_pairs
        max_point_pair_distance <- opt_summary$max_point_pair_distance

        all_models <- unique(nlopt_iter_df$variogram_model)
        other_models <- setdiff(all_models, opt_variogram_model)
        other_lt_list <- c("dashed", "dotted", "dotdash", "longdash", "twodash")
        linetype_map <- setNames(rep(other_lt_list, length.out = length(other_models)), other_models)
        linetype_map[opt_variogram_model] <- "solid"

        p_title <- paste0("Kriging Optimization: ", submap_transect_set, "")

        p0 <- nlopt_iter_df |>
          ggplot(aes(x = itr, y = max_point_pair_distance / 1e3, group = variogram_model, linetype = variogram_model)) +
          geom_path(color = "forestgreen", linewidth = 1.5) +
          labs(x = "NLopt Iter's", y = "Max Pt. Pair Dist. (km)") +
          theme_vgrm(base_size)

        p1 <- nlopt_iter_df |>
          ggplot(aes(x = itr, y = max_n_point_pairs, group = variogram_model, linetype = variogram_model)) +
          geom_path(color = "firebrick", linewidth = 1.5) +
          labs(x = "NLopt Iter's", y = "Max No. Pt. Pairs", linetype = "Model") +
          guides(linetype = guide_legend(override.aes = list(color = "black"))) +
          theme_vgrm(base_size, show_legend = TRUE)

        p2 <- nlopt_iter_df |>
          ggplot(aes(x = itr, y = variogram_cutoff, group = variogram_model, linetype = variogram_model)) +
          geom_path(color = "orchid4", linewidth = 1.5) +
          labs(x = "NLopt Iter's", y = "V'gram Cutoff") +
          theme_vgrm(base_size)

        p3 <- nlopt_iter_df |>
          ggplot(aes(x = itr, y = n_variogram_lags, group = variogram_model, linetype = variogram_model)) +
          geom_path(color = "saddlebrown", linewidth = 1.5) +
          labs(x = "NLopt Iter's", y = "No. V'gram Lags") +
          theme_vgrm(base_size)

        p4 <- nlopt_iter_df |>
          ggplot(aes(x = itr, y = total_cost, group = variogram_model, linetype = variogram_model)) +
          geom_path(color = "black", linewidth = 1.5) +
          labs(x = "NLopt Iter's", y = "Total Cost") +
          theme_vgrm(base_size)

        p5 <- experimental_variogram |>
          ggplot(aes(x = dist / 1e3, y = sqrt(gamma))) +
          geom_line(data = fitted_variogram_line, linewidth = 1.5) +
          geom_point(shape = 20, size = 3) +
          labs(x = "Lag Distance (km)", y = bquote("Variance" ~ (mWm^-2))) +
          theme_vgrm(base_size)

        p <- (p0 + p1) / (p2 + p3) / (p4 + p5) +
          plot_annotation(
            title = p_title,
            theme = theme(plot.title = element_text(size = base_size * 1.2, hjust = 0, margin = margin(10, 0, -10, 60)))
          ) &
          scale_linetype_manual(name = "Model", values = linetype_map)

        ggsave(file = out_path, plot = p, width = 13, height = 10, dpi = 300, bg = "white")
      })
    })
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_submap_transect_set_composition <- function(data_dir, out_dir, fig_dir, submap_transect_set, ..., base_size = 14) {
  with_error_handling({
    suppressWarnings({
      suppressMessages({
        out_paths <- file.path(fig_dir, "submap", paste0(submap_transect_set, c("-obs.png", "-sim.png", "-krg.png", "-comp.png")))

        if (all(map_lgl(out_paths, check_plot_path))) {
          return(invisible())
        }

        opt_data <- get_optimal_krige_model(out_dir, submap_transect_set, ...)
        variogram_model <- opt_data$optimal_krige_model$variogram_model
        fitted_variogram <- opt_data$fitted_variogram
        max_n_point_pairs <- opt_data$optimal_krige_model$max_n_point_pairs
        max_point_pair_distance <- opt_data$optimal_krige_model$max_point_pair_distance

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

        proj <- determine_projection(data$hull, proj_type = "albers")

        sf_transect <- data$submap_transect_geometry
        sf_hull <- data$hull
        sf_ihfc_hull <- data$ihfc2024_obs_hull[[1]]

        sr_sim_hull <- with_error_handling(to_spatraster(data$lucazeau2019_sim_hull[[1]], data_col = "lucazeau2019_sim_est"), default = NULL)

        prefix <- "kerswell2025_krg_"
        model_str <- str_to_lower(variogram_model)
        krg_hull_str <- paste0(prefix, model_str, "_hull")
        sr_krige_hull <- with_error_handling(to_spatraster(data[[krg_hull_str]][[1]], data_col = "kerswell2025_krg_est"), default = NULL)

        p1 <- ggplot() +
          basemap_layers(data_dir, out_dir, proj) +
          geom_sf(data = sf_ihfc_hull, aes(color = ihfc2024_obs), shape = 20, size = 2) +
          ggtitle("Global HF Observations") +
          viridis_scale_color() +
          theme_transect_set(sf_hull, proj, base_size, show_legend = FALSE)

        ggsave(file = out_paths[1], plot = p1, width = 5.5, height = 5, dpi = 300, bg = "white")

        p2 <- ggplot() +
          basemap_layers(data_dir, out_dir, proj) +
          geom_spatraster(data = sr_sim_hull) +
          ggtitle("Similarity Interpolation") +
          viridis_scale_color() +
          viridis_scale_fill() +
          theme_transect_set(sf_hull, proj, base_size, show_legend = FALSE)

        ggsave(file = out_paths[2], plot = p2, width = 5.5, height = 5, dpi = 300, bg = "white")

        if (is.null(sr_krige_hull)) {
          p3 <- ggplot() +
            basemap_layers(data_dir, out_dir, proj) +
            geom_sf(data = sf_hull, linewidth = 1, fill = NA, color = "black") +
            ggtitle("Krige Interpolation (failed)") +
            theme_transect_set(sf_hull, proj, base_size, show_legend = FALSE)
        } else {
          p3 <- ggplot() +
            basemap_layers(data_dir, out_dir, proj) +
            geom_spatraster(data = sr_krige_hull) +
            ggtitle("Krige Interpolation") +
            viridis_scale_color() +
            viridis_scale_fill() +
            theme_transect_set(sf_hull, proj, base_size, show_legend = FALSE)
        }

        ggsave(file = out_paths[3], plot = p3, width = 5.5, height = 5, dpi = 300, bg = "white")

        no_lat <- list(theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()))
        p <- (p1 | p2 + no_lat | p3 + no_lat)

        ggsave(file = out_paths[4], plot = p, width = 13, height = 5, dpi = 300, bg = "white")
      })
    })
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_submap_transect_composition <- function(data_dir, out_dir, fig_dir, submap_transect_id, ..., base_size = 14) {
  with_error_handling({
    suppressWarnings({
      suppressMessages({
        out_paths <- file.path(fig_dir, "submap", paste0(submap_transect_id, c("-obs.png", "-sim.png", "-krg.png", "-comp.png")))

        if (all(map_lgl(out_paths, check_plot_path))) {
          return(invisible())
        }

        opt_data <- get_optimal_krige_model(out_dir, submap_transect_id, ...)
        variogram_model <- opt_data$optimal_krige_model$variogram_model
        fitted_variogram <- opt_data$fitted_variogram
        max_n_point_pairs <- opt_data$optimal_krige_model$max_n_point_pairs
        max_point_pair_distance <- opt_data$optimal_krige_model$max_point_pair_distance

        cache_path <- file.path(out_dir, "submap", paste0(submap_transect_id, ".RData"))
        data <- compile_submap_transect_data(
          cache_path,
          data_dir,
          out_dir,
          submap_transect_id,
          variogram_model = variogram_model,
          fitted_variogram = fitted_variogram,
          max_n_point_pairs = max_n_point_pairs,
          max_point_pair_distance = max_point_pair_distance,
          ...
        )

        proj <- determine_projection(data$hull, proj_type = "albers")

        sf_transect <- data$submap_transect_geometry
        sf_ihfc_hull <- data$ihfc2024_obs_hull[[1]]
        sf_ihfc_proj <- data$ihfc2024_obs_projected[[1]]
        ihfc_loess <- data$ihfc2024_obs_loess[[1]]

        sr_sim_hull <- with_error_handling(to_spatraster(data$lucazeau2019_sim_hull[[1]], data_col = "lucazeau2019_sim_est"), default = NULL)
        sf_sim_proj <- data$lucazeau2019_sim_projected[[1]]
        sim_loess <- data$lucazeau2019_sim_loess[[1]]

        prefix <- "kerswell2025_krg_"
        model_str <- str_to_lower(variogram_model)
        krg_hull_str <- paste0(prefix, model_str, "_hull")
        krg_projected_str <- paste0(prefix, model_str, "_projected")
        krg_loess_str <- paste0(prefix, model_str, "_loess")

        sr_krige_hull <- with_error_handling(to_spatraster(data[[krg_hull_str]][[1]], data_col = "kerswell2025_krg_est"), default = NULL)
        sf_krg_proj <- data[[krg_projected_str]][[1]]
        krg_loess <- data[[krg_loess_str]][[1]]

        p1 <- ggplot() +
          basemap_layers_transect(data_dir, out_dir, data, 1) +
          geom_sf(data = sf_ihfc_hull, aes(color = ihfc2024_obs), shape = 20, size = 2) +
          geom_sf(data = sf_transect, linewidth = 1, color = "black") +
          ggtitle("Global HF Observations") +
          viridis_scale_color() +
          theme_transect(proj$wkt, base_size, show_legend = FALSE)

        if (is.null(sf_ihfc_proj)) {
          p2 <- ggplot() +
            theme_void()
        } else {
          p2 <-
            ggplot() +
            geom_point(data = sf_ihfc_proj, aes(distance, hf), shape = 20, color = "grey20", size = 0.9)
          if (!is.null(ihfc_loess)) {
            p2 <- p2 +
              geom_path(data = ihfc_loess, aes(distance, hf), linewidth = 2, color = "black")
          }
        }

        p2 <- p2 +
          labs(x = "Normalized Distance", y = bquote("Q" ~ (mWm^-2))) +
          viridis_scale_color() +
          theme_loess(base_size)

        pp1 <- (p1 / p2) + plot_layout(widths = 1, heights = c(1.8, 1))
        ggsave(file = out_paths[1], plot = pp1, width = 6.5, height = 10, dpi = 300, bg = "white")

        p3 <- ggplot() +
          basemap_layers_transect(data_dir, out_dir, data, 1) +
          geom_spatraster(data = sr_sim_hull) +
          geom_sf(data = sf_transect, linewidth = 1, color = "black") +
          ggtitle("Similarity Interpolation") +
          viridis_scale_color() +
          viridis_scale_fill() +
          theme_transect(proj$wkt, base_size, show_legend = FALSE)

        if (is.null(sf_sim_proj)) {
          p4 <- ggplot() +
            theme_void()
        } else {
          p4 <- ggplot() +
            geom_point(data = sf_sim_proj, aes(distance, hf), shape = 20, color = "grey20", size = 0.9)
          if (!is.null(sim_loess)) {
            p4 <- p4 +
              geom_path(data = sim_loess, aes(distance, hf), linewidth = 2, color = "black")
          }
        }

        p4 <- p4 +
          labs(x = "Normalized Distance", y = bquote("Q" ~ (mWm^-2))) +
          viridis_scale_color() +
          theme_loess(base_size)

        pp2 <- (p3 / p4) + plot_layout(widths = 1, heights = c(1.8, 1))
        ggsave(file = out_paths[2], plot = pp2, width = 6.5, height = 10, dpi = 300, bg = "white")

        if (is.null(sr_krige_hull)) {
          p5 <- ggplot() +
            basemap_layers_transect(data_dir, out_dir, data, 1) +
            geom_sf(data = sf_transect, linewidth = 1, color = "black") +
            ggtitle("Krige Interpolation (failed)") +
            theme_transect(proj$wkt, base_size, show_legend = FALSE)
        } else {
          p5 <- ggplot() +
            basemap_layers_transect(data_dir, out_dir, data, 1) +
            geom_spatraster(data = sr_krige_hull) +
            geom_sf(data = sf_transect, linewidth = 1, color = "black") +
            ggtitle("Krige Interpolation") +
            viridis_scale_color() +
            viridis_scale_fill() +
            theme_transect(proj$wkt, base_size, show_legend = FALSE)
        }

        if (is.null(sf_krg_proj)) {
          p6 <- ggplot() +
            theme_void()
        } else {
          p6 <- ggplot() +
            geom_point(data = sf_krg_proj, aes(distance, hf), shape = 20, color = "grey20", size = 0.9)
          if (!is.null(krg_loess)) {
            p6 <- p6 +
              geom_path(data = krg_loess, aes(distance, hf), linewidth = 2, color = "black")
          }
        }

        p6 <- p6 +
          labs(x = "Normalized Distance", y = bquote("Q" ~ (mWm^-2))) +
          viridis_scale_color() +
          theme_loess(base_size)

        pp3 <- (p5 / p6) + plot_layout(widths = 1, heights = c(1.8, 1))
        ggsave(file = out_paths[3], plot = pp3, width = 6.5, height = 10, dpi = 300, bg = "white")

        p <- pp1 | pp2 | pp3 + plot_annotation(title = p_title)
        ggsave(file = out_paths[4], plot = p, width = 19.5, height = 10, dpi = 300, bg = "white")
      })
    })
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_submap_transect_sets <- function(data_dir, out_dir, fig_dir, submap_transect_sets, ...) {
  with_error_handling({
    dots <- list(...)
    seed <- dots$seed %||% 42
    nprocs <- dots$nprocs %||% availableCores() - 2
    parallel <- dots$parallel %||% FALSE

    inner_dots <- dots
    inner_dots$parallel <- NULL
    inner_dots$nprocs <- NULL
    inner_dots$seed <- NULL

    f <- function(data_dir, out_dir, fig_dir, submap_transect_set, ...) {
      with_error_handling(
        {
          suppressWarnings({
            suppressMessages({
              draw_submap_transect_set_composition(data_dir, out_dir, fig_dir, submap_transect_set, ...)
              draw_optimal_krige_model(data_dir, out_dir, fig_dir, submap_transect_set, ...)
              # future_walk(submap_transect_sets, draw_interp_accuracy_summary, .options = furrr_options(seed = seed))
            })
          })
        },
        default = NULL
      )
    }

    f_partial <- partial(f, data_dir = data_dir, out_dir = out_dir, fig_dir = fig_dir, !!!inner_dots)

    if (parallel) {
      set.seed(seed)
      plan(multisession, workers = nprocs)
      future_walk(submap_transect_sets, f_partial, .options = furrr_options(seed = seed))
    } else {
      walk(submap_transect_sets, f_partial)
    }
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_submap_transects <- function(data_dir, out_dir, fig_dir, submap_transect_ids, ...) {
  with_error_handling({
    dots <- list(...)
    seed <- dots$seed %||% 42
    nprocs <- dots$nprocs %||% availableCores() - 2
    parallel <- dots$parallel %||% FALSE

    zones <- c("NPA", "SAM", "SEA", "SWP")

    f <- function(data_dir, out_dir, fig_dir, x, ...) {
      with_error_handling(
        {
          suppressWarnings({
            suppressMessages({
              draw_submap_transect_composition(data_dir, out_dir, fig_dir, x, ...)
              # future_walk(ids, draw_optimal_krige_model, .options = furrr_options(seed = seed))
              # future_walk(zones, draw_interp_accuracy_summary, .options = furrr_options(seed = seed))
            })
          })
        },
        default = NULL
      )
    }

    if (parallel) {
      set.seed(seed)
      plan(multisession, workers = nprocs)

      future_walk(submap_transect_ids, function(x) {
        f(data_dir, out_dir, fig_dir, x, !!!dots)
      }, .options = furrr_options(seed = seed))
    } else {
      walk(submap_transect_ids, function(x) {
        f(data_dir, out_dir, fig_dir, x, ...)
      })
    }
  })
}
