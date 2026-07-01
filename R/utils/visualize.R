#######################################################
## Helper functions (visualizations)                 ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
viridis_scale_color <- function() {
  list(
    scale_color_viridis_c(
      option = "magma",
      name = bquote(mWm^-2),
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
      name = bquote(mWm^-2),
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
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "grey90"),
      panel.spacing = unit(1, "lines"),
      plot.margin = margin(5, 15, 5, 5),
      plot.title = element_text(hjust = 0.5, size = base_size * 1.1),
      axis.ticks = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size * 1.1),
      legend.justification = "center",
      legend.position = "bottom",
      legend.key.width = unit(2, "lines"),
      legend.box.margin = margin(2, 2, 2, 2),
      legend.margin = margin(),
      legend.title = element_text(vjust = 0.5),
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
theme_vgrm <- function(base_size = 24, show_legend = FALSE) {
  theme <- list(
    theme_bw(base_size = base_size),
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "grey90"),
      panel.spacing = unit(1, "lines"),
      plot.margin = margin(5, 15, 5, 5),
      plot.title = element_text(hjust = 0, margin = margin(10, 0, -10, 60)),
      axis.ticks = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size),
      legend.justification = "right",
      legend.position = "top",
      legend.key.height = unit(0.7, "lines"),
      legend.key.width = unit(2, "lines"),
      legend.position.inside = c(0.98, 0.73),
      legend.box.margin = margin(2, 2, 2, 2),
      legend.margin = margin(),
      legend.title = element_text(hjust = 1, margin = margin(0, 15, 0, 0)),
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
theme_trend <- function(base_size = 32, show_legend = TRUE) {
  theme <- list(
    scale_y_continuous(limits = c(0, 250), breaks = seq(0, 250, 50)),
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)),
    theme_bw(base_size = base_size),
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "grey90"),
      plot.margin = margin(5, 5, 5, 5),
      axis.ticks = element_blank(),
      legend.justification = "left",
      legend.position = "inside",
      legend.position.inside = c(0.05, 0.82),
      legend.direction = "horizontal",
      legend.key.height = unit(1.2, "lines"),
      legend.key.width = unit(2, "lines"),
      legend.box.margin = margin(2, 2, 2, 2),
      legend.margin = margin(),
      legend.title.position = "left",
      legend.title = element_text(vjust = 1.0, margin = margin(0, 10, 0, 0), size = base_size * 0.8),
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
      legend.margin = margin(),
      legend.key.height = unit(0.8, "lines"),
      legend.key.width = unit(2, "lines"),
      legend.frame = element_rect(linewidth = 0.6),
      legend.ticks = element_line(linewidth = 0.6),
      legend.title = element_text(vjust = 1.2, margin = margin(0, 5, 0, 0)),
      legend.title.position = "left"
    )
  )

  if (!show_legend) {
    theme <- c(theme, list(theme(legend.position = "none")))
  }

  theme
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
theme_transect_set <- function(sf_obj, proj = NULL, base_size = 32, show_legend = TRUE) {
  if (!is.null(proj)) wkt <- proj$wkt else st_crs("4326")
  if (!is.null(proj)) bbox <- st_bbox(st_transform(sf_obj, wkt)) else st_bbox(sf_obj)

  theme <- list(
    annotation_scale(location = "bl", width_hint = 0.33, height = unit(0.8, "lines"), text_cex = 2, style = "bar", line_width = 4),
    annotation_north_arrow(
      location = "bl",
      which_north = "true",
      pad_x = unit(0, "lines"),
      height = unit(6, "lines"),
      width = unit(6, "lines"),
      pad_y = unit(2, "lines"),
      style = north_arrow_fancy_orienteering
    ),
    scale_x_continuous(breaks = seq(-180, 180, by = 5)),
    scale_y_continuous(breaks = seq(-90, 90, by = 5)),
    coord_sf(crs = wkt, xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE),
    theme_minimal(base_size = base_size),
    theme(
      panel.grid.major = element_blank(),
      plot.margin = margin(5, 5, 5, 5),
      plot.title = element_text(hjust = 0.5, size = base_size * 1.1),
      axis.title = element_blank()
    )
  )

  if (!show_legend) {
    theme <- c(theme, list(theme(legend.position = "none")))
  }

  theme
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
basemap_layers <- function(data_dir, out_dir, proj = NULL, major_grats = TRUE, minor_grats = TRUE, quiet = FALSE) {
  with_error_handling(
    {
      wrap_opts <- c("WRAPDATELINE=YES", "DATELINEOFFSET=180")
      if (!is.null(proj)) lon0 <- proj$lon0 else 0

      etopo_dir <- file.path(data_dir, "mapping", "relief")
      geotif_dir <- file.path(out_dir, "relief")

      sr_relief <- with_error_handling(get_noaa_relief(-180, 180, 90, -90, etopo_dir, geotif_dir), default = NULL, quiet = quiet)
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
basemap_layers_transect <- function(data_dir, out_dir, data, row = 1, quiet = FALSE) {
  with_error_handling(
    {
      hull_geom <- data$hull[[row]]

      etopo_dir <- file.path(data_dir, "mapping", "relief")
      geotif_dir <- file.path(out_dir, "relief")
      sr_relief <- with_error_handling(get_noaa_relief(-180, 180, 90, -90, etopo_dir, geotif_dir), default = NULL, quiet = quiet)

      crop_sf <- function(shp) {
        if (is.null(shp)) {
          return(NULL)
        }
        suppressWarnings(st_crop(shp, st_bbox(hull_geom)))
      }

      sf_coast_c <- crop_sf(sf_coast)
      sf_ridge_c <- crop_sf(sf_ridge)
      sf_transform_c <- crop_sf(sf_transform)
      sf_trench_c <- crop_sf(sf_trench)

      list(
        if (!is.null(sr_relief)) geom_spatraster(data = sr_relief) else NULL,
        if (!is.null(sr_relief)) elevation_scale_fill() else NULL,
        if (!is.null(sr_relief)) new_scale_fill() else NULL,
        geom_sf(data = sf_coast_c, linewidth = 0.2, color = "black"),
        geom_sf(data = sf_ridge_c, linewidth = 0.6, color = "white"),
        geom_sf(data = sf_transform_c, linewidth = 0.6, color = "white"),
        geom_sf(data = sf_trench_c, linewidth = 0.6, color = "white")
      )
    },
    default = NULL
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
build_strip_buffer_outlines <- function(transect_geometry, proj, strip_widths = c(7.5e4)) {
  with_error_handling(
    {
      transect_proj <- st_transform(transect_geometry, proj$wkt)
      map(strip_widths, ~ {
        st_buffer(transect_proj, .x, endCapStyle = "FLAT") |> st_transform(4326)
      })
    },
    default = vector("list", length(strip_widths))
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
build_hf_profile_plot <- function(proj_list, trend_list, y_lab = bquote("Q" ~ (mWm^-2)), base_size = 32, show_legend = TRUE) {
  if (length(proj_list) == 0 || is.null(proj_list[[1]]) || nrow(proj_list[[1]]) == 0) {
    return(
      ggplot() +
        labs(x = "Normalized Distance", y = y_lab) +
        theme_trend(base_size, show_legend = FALSE)
    )
  }

  combined_data <- map2_dfr(proj_list, trend_list, function(proj, trend) {
    if (is.null(proj) || nrow(proj) == 0) {
      return(NULL)
    }
    bind_rows(proj, trend)
  })

  pts_data <- proj_list |>
    keep(~ !is.null(.x) && nrow(.x) > 0) |>
    bind_rows()

  valid_indices <- which(map_lgl(proj_list, ~ !is.null(.x) && nrow(.x) > 0))
  rug_data <- proj_list[[max(valid_indices)]]

  p <- ggplot() +
    geom_point(data = pts_data, aes(x = distance, y = hf, color = hf), shape = 16, size = 3.5, alpha = 0.3) +
    viridis_scale_color() +
    labs(x = "Normalized Distance", y = y_lab) +
    theme_trend(base_size, show_legend = show_legend)

  if (!is.null(trend_list) && length(trend_list) > 0) {
    trends_bound <- trend_list |>
      keep(~ !is.null(.x) && nrow(.x) > 0) |>
      bind_rows()

    if (nrow(trends_bound) > 0) {
      p + geom_line(data = filter(trends_bound, method == "mean"), aes(x = distance, y = hf, group = submap_transect_id), na.rm = TRUE)
    }
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
extract_trend_vector <- function(trend_entry, n_grid = 500) {
  if (is.null(trend_entry) || !is.data.frame(trend_entry) || nrow(trend_entry) == 0) {
    return(rep(NA_real_, n_grid))
  }
  tr <- filter(trend_entry, method == "mean")
  if (nrow(tr) < 5 || all(is.na(tr$hf))) {
    return(rep(NA_real_, n_grid))
  }
  approx(tr$distance, tr$hf, xout = seq(0, 1, length.out = n_grid), rule = 2)$y
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
max_ccf_value <- function(v1, v2, max_lag_frac = 0.1, n_grid = 500) {
  if (all(is.na(v1)) || all(is.na(v2))) {
    return(NA_real_)
  }
  max_lags <- max(1L, floor(n_grid * max_lag_frac))
  tryCatch(
    {
      res <- ccf(v1, v2, lag.max = max_lags, na.action = na.pass, plot = FALSE)
      max(res$acf, na.rm = TRUE)
    },
    error = function(e) NA_real_
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
make_krg_trend_col <- function(variogram_model, strip_idx) {
  if (is.null(variogram_model) || is.na(variogram_model)) {
    return(NULL)
  }
  paste0("kerswell2025_krg_", str_to_lower(variogram_model), "_strip", strip_idx, "_trend")
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
make_krg_proj_col <- function(variogram_model, strip_idx) {
  if (is.null(variogram_model) || is.na(variogram_model)) {
    return(NULL)
  }
  paste0("kerswell2025_krg_", str_to_lower(variogram_model), "_strip", strip_idx, "_proj")
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
precompile_set_data <- function(data_dir, out_dir, submap_transect_sets, quiet = FALSE) {
  if (!quiet) cat(" .. Pre-compiling Kriging data cache files\n")

  walk(submap_transect_sets, function(set_id) {
    with_error_handling(
      load_set_data(data_dir = data_dir, out_dir = out_dir, submap_transect_set = set_id, quiet = quiet),
      default = NULL,
      quiet = quiet
    )
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load_set_data <- function(data_dir, out_dir, submap_transect_set, quiet = FALSE) {
  with_error_handling(
    {
      opt <- get_optimal_krige_model(out_dir = out_dir, submap_transect_set = submap_transect_set, quiet = quiet)
      variogram_model <- opt$optimal_krige_model$variogram_model
      fitted_variogram <- opt$fitted_variogram
      max_n_point_pairs <- opt$optimal_krige_model$max_n_point_pairs
      max_point_pair_distance <- opt$optimal_krige_model$max_point_pair_distance

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
      list(data = data, variogram_model = variogram_model)
    },
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_method_trend_vectors <- function(data, trend_col, n_grid = 500) {
  empty <- vector("list", nrow(data)) |> map(~ rep(NA_real_, n_grid))
  if (is.null(trend_col) || !trend_col %in% names(data)) {
    return(empty)
  }
  map(data[[trend_col]], extract_trend_vector, n_grid = n_grid)
}

#######################################################
## Plotting functions                                ##
#######################################################
draw_global_ccf_summary <- function(fig_dir, global_ccf_summary, method_ccf_summary, base_size = 22, quiet = FALSE) {
  if (is.null(global_ccf_summary) || nrow(global_ccf_summary) == 0) {
    return(invisible())
  }

  if (is.null(method_ccf_summary) || nrow(method_ccf_summary) == 0) {
    return(invisible())
  }

  with_error_handling(
    {
      out_path <- file.path(fig_dir, "summary", "global-ccf-summary.png")
      if (check_plot_path(out_path)) {
        return(invisible())
      }

      ccf_long <- global_ccf_summary |>
        pivot_longer(c(obs_ccf, sim_ccf, krg_ccf), names_to = "method", values_to = "max_ccf") |>
        mutate(
          method = recode(method, obs_ccf = "Observation Profiles", sim_ccf = "Similarity Profiles", krg_ccf = "Krige Profiles"),
          method = factor(method, levels = c("Observation Profiles", "Krige Profiles", "Similarity Profiles"))
        )

      ccf_long <- ccf_long |>
        mutate(set_num = as.integer(str_extract(submap_transect_set, "(?<=SET)\\d+"))) |>
        arrange(region, set_num, global_pair_uid)

      uid_order <- ccf_long |>
        distinct(region, submap_transect_set, set_num, global_pair_uid) |>
        arrange(region, set_num, global_pair_uid) |>
        pull(global_pair_uid)

      ccf_long <- ccf_long |> mutate(global_pair_uid = factor(global_pair_uid, levels = uid_order))
      region_order <- ccf_long |>
        distinct(region) |>
        pull(region) |>
        sort() |>
        as.character()

      region_palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
      cols <- setNames(region_palette[seq_along(region_order)], region_order)

      region_shapes <- c(21, 22, 23, 24)
      shapes <- setNames(region_shapes[seq_along(region_order)], region_order)

      ccf_long <- ccf_long |> mutate(region = factor(region, levels = region_order))

      label_map <- ccf_long |>
        distinct(global_pair_uid, pair_label) |>
        deframe()
      last_facet_method <- levels(ccf_long$method)[length(levels(ccf_long$method))]

      text_data <- ccf_long |>
        filter(method == last_facet_method) |>
        arrange(region, set_num, global_pair_uid) |>
        group_by(region, submap_transect_set) |>
        slice(1) |>
        ungroup() |>
        mutate(fontface = if_else(interpretation %in% c("cont", "discont", "disagree"), "bold", "plain"))

      ccf_indexed <- ccf_long |> mutate(x_index = match(global_pair_uid, uid_order))

      rect_data <- ccf_indexed |>
        group_by(submap_transect_set, interpretation) |>
        summarise(xmin = min(x_index) - 0.5, xmax = max(x_index) + 0.5, .groups = "drop") |>
        mutate(ymin = -1, ymax = 1) |>
        mutate(interpretation = recode(interpretation,
          "cont" = "Inferred Continuous",
          "discont" = "Inferred Discontinuous",
          "disagree" = "Method Disagreement",
          "ambiguous" = "Ambiguous"
        ))

      rect_data <- rect_data |>
        mutate(
          interpretation = factor(interpretation, levels = c("Inferred Continuous", "Inferred Discontinuous", "Method Disagreement", "Ambiguous")),
          pattern_type = case_when(
            interpretation == "Inferred Continuous" ~ "crosshatch",
            interpretation == "Inferred Discontinuous" ~ "stripe",
            interpretation == "Method Disagreement" ~ "stripe",
            interpretation == "Ambiguous" ~ "none"
          ),
          pattern_angle = case_when(
            interpretation == "Inferred Discontinuous" ~ 45,
            interpretation == "Method Disagreement" ~ 135,
            TRUE ~ 45
          )
        )

      interpretation_colors <- c(
        "Inferred Continuous" = "#BFE5D9",
        "Inferred Discontinuous" = "#FFD1B3",
        "Method Disagreement" = "#E1D5E7",
        "Ambiguous" = "#E6E6E6"
      )

      p <- ggplot(
        ccf_long,
        aes(
          x = global_pair_uid,
          y = max_ccf,
          group = interaction(region, submap_transect_set),
          color = region,
          fill = region,
          size = n_total
        )
      ) +
        geom_rect_pattern(
          data = rect_data,
          aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, pattern = interpretation, pattern_angle = interpretation),
          inherit.aes = FALSE,
          fill = NA,
          color = NA,
          pattern_color = "black",
          pattern_fill = NA,
          pattern_density = 0.08,
          pattern_spacing = 0.025,
          pattern_alpha = 0.5,
          pattern_size = 0.3
        ) +
        scale_pattern_manual(
          name = "Classification",
          values = c(
            "Inferred Continuous" = "crosshatch",
            "Inferred Discontinuous" = "stripe",
            "Method Disagreement" = "stripe",
            "Ambiguous" = "none"
          ),
          guide = guide_legend(override.aes = list(fill = "transparent"))
        ) +
        scale_pattern_angle_manual(
          name = "Classification",
          values = c(
            "Inferred Continuous" = 45,
            "Inferred Discontinuous" = 45,
            "Method Disagreement" = 135,
            "Ambiguous" = 0
          )
        ) +
        geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
        geom_line(linewidth = 1.2, na.rm = TRUE) +
        geom_point(aes(shape = region), na.rm = TRUE, color = "black", stroke = 0.8) +
        geom_text(
          data = text_data,
          aes(x = global_pair_uid, y = -1.0, label = submap_transect_set, fontface = fontface),
          inherit.aes = FALSE,
          angle = 90,
          hjust = -0.05,
          vjust = 0.5,
          size = base_size * 0.22
        ) +
        scale_size_continuous(name = "Observations", range = c(0, 7), breaks = c(0, 50, 100, 1000)) +
        scale_fill_manual(name = "SubMap Region", values = cols, limits = region_order, breaks = region_order) +
        scale_color_manual(name = "SubMap Region", values = cols, limits = region_order, breaks = region_order) +
        scale_shape_manual(name = "SubMap Region", values = shapes, limits = region_order, breaks = region_order) +
        scale_y_continuous(name = "Cross Correlation Function (CCF)", limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
        scale_x_discrete(name = "Adjacent SubMap Transect Pairs", labels = label_map, expand = expansion(add = c(2, 2))) +
        coord_cartesian(clip = "off") +
        facet_wrap(~method, scales = "free_x", ncol = 1) +
        theme_facet(base_size, show_legend = TRUE) +
        theme(
          axis.text.x = element_blank(),
          axis.title.x = element_text(margin = margin(20, 0, 0, 0)),
          axis.ticks.x = element_line(linewidth = 0.3),
          axis.ticks.length.x = unit(-0.3, "lines"),
          legend.title.position = "top",
          legend.direction = "vertical",
          legend.box = "horizontal"
        )

      ggsave(file = out_path, plot = p, width = 14, height = 12, dpi = 300, bg = "white")
    },
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_method_ccf_summary <- function(fig_dir, global_ccf_summary, method_ccf_summary, base_size = 22, quiet = FALSE) {
  if (is.null(global_ccf_summary) || nrow(global_ccf_summary) == 0) {
    return(invisible())
  }

  if (is.null(method_ccf_summary) || nrow(method_ccf_summary) == 0) {
    return(invisible())
  }

  with_error_handling(
    {
      out_path <- file.path(fig_dir, "summary", "method-ccf-summary.png")
      if (check_plot_path(out_path)) {
        return(invisible())
      }

      ccf_long <- method_ccf_summary |>
        pivot_longer(c(obs_sim, obs_krg, sim_krg), names_to = "pair", values_to = "max_ccf") |>
        mutate(
          method_pair = recode(pair, obs_sim = "Obs–Sim Profile Pairs", obs_krg = "Obs–Krg Profile Pairs", sim_krg = "Sim–Krg Profile Pairs"),
          method_pair = factor(method_pair, levels = c("Obs–Sim Profile Pairs", "Obs–Krg Profile Pairs", "Sim–Krg Profile Pairs")),
          set_num = as.integer(str_extract(submap_transect_set, "(?<=SET)\\d+"))
        )

      x_levels <- ccf_long |>
        arrange(region, set_num, global_pair_uid) |>
        pull(global_pair_uid) |>
        unique()
      ccf_long <- ccf_long |> mutate(global_pair_uid = factor(global_pair_uid, levels = x_levels))

      region_order <- ccf_long |>
        distinct(region) |>
        pull(region) |>
        sort() |>
        as.character()

      region_palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
      cols <- setNames(region_palette[seq_along(region_order)], region_order)

      region_shapes <- c(21, 22, 23, 24)
      shapes <- setNames(region_shapes[seq_along(region_order)], region_order)

      ccf_long <- ccf_long |> mutate(region = factor(region, levels = region_order))

      label_map <- ccf_long |>
        distinct(global_pair_uid, transect_id) |>
        mutate(short_label = str_remove_all(transect_id, "^[^0-9]+")) |>
        select(global_pair_uid, short_label) |>
        deframe()

      last_facet_pair <- levels(ccf_long$method_pair)[length(levels(ccf_long$method_pair))]

      text_data <- ccf_long |>
        filter(method_pair == last_facet_pair) |>
        arrange(region, set_num, global_pair_uid) |>
        group_by(region, submap_transect_set) |>
        slice(1) |>
        ungroup() |>
        mutate(fontface = if_else(interpretation %in% c("cont", "discont", "disagree"), "bold", "plain"))

      df_indexed <- ccf_long |> mutate(x_index = match(global_pair_uid, x_levels))

      rect_data <- df_indexed |>
        group_by(submap_transect_set, interpretation) |>
        summarise(xmin = min(x_index) - 0.5, xmax = max(x_index) + 0.5, .groups = "drop") |>
        mutate(ymin = -1, ymax = 1) |>
        mutate(interpretation = recode(interpretation,
          "cont" = "Inferred Continuous",
          "discont" = "Inferred Discontinuous",
          "disagree" = "Method Disagreement",
          "ambiguous" = "Ambiguous"
        ))

      rect_data <- rect_data |>
        mutate(
          interpretation = factor(interpretation, levels = c("Inferred Continuous", "Inferred Discontinuous", "Method Disagreement", "Ambiguous")),
          pattern_type = case_when(
            interpretation == "Inferred Continuous" ~ "crosshatch",
            interpretation == "Inferred Discontinuous" ~ "stripe",
            interpretation == "Method Disagreement" ~ "stripe",
            interpretation == "Ambiguous" ~ "none"
          ),
          pattern_angle = case_when(
            interpretation == "Inferred Discontinuous" ~ 45,
            interpretation == "Method Disagreement" ~ 135,
            TRUE ~ 45
          )
        )

      interpretation_colors <- c(
        "Inferred Continuous" = "#BFE5D9",
        "Inferred Discontinuous" = "#FFD1B3",
        "Method Disagreement" = "#E1D5E7",
        "Ambiguous" = "#E6E6E6"
      )

      p <- ggplot(
        ccf_long,
        aes(
          x = global_pair_uid,
          y = max_ccf,
          group = interaction(region, submap_transect_set),
          color = region,
          fill = region,
          size = n
        )
      ) +
        geom_rect_pattern(
          data = rect_data,
          aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, pattern = interpretation, pattern_angle = interpretation),
          inherit.aes = FALSE,
          fill = NA,
          color = NA,
          pattern_color = "black",
          pattern_fill = NA,
          pattern_density = 0.08,
          pattern_spacing = 0.025,
          pattern_alpha = 0.5,
          pattern_size = 0.3
        ) +
        scale_pattern_manual(
          name = "Classification",
          values = c(
            "Inferred Continuous" = "crosshatch",
            "Inferred Discontinuous" = "stripe",
            "Method Disagreement" = "stripe",
            "Ambiguous" = "none"
          ),
          guide = guide_legend(override.aes = list(fill = "transparent"))
        ) +
        scale_pattern_angle_manual(
          name = "Classification",
          values = c(
            "Inferred Continuous" = 45,
            "Inferred Discontinuous" = 45,
            "Method Disagreement" = 135,
            "Ambiguous" = 0
          )
        ) +
        geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
        geom_line(linewidth = 1.2, na.rm = TRUE) +
        geom_point(aes(shape = region), na.rm = TRUE, color = "black", stroke = 0.8) +
        geom_text(
          data = text_data,
          aes(x = global_pair_uid, y = -1.0, label = submap_transect_set, fontface = fontface),
          inherit.aes = FALSE,
          angle = 90,
          hjust = -0.05,
          vjust = 0.5,
          size = base_size * 0.22
        ) +
        scale_size_continuous(name = "Observations", range = c(0, 7), breaks = c(0, 50, 100, 1000)) +
        scale_fill_manual(name = "SubMap Region", values = cols, limits = region_order, breaks = region_order) +
        scale_color_manual(name = "SubMap Region", values = cols, limits = region_order, breaks = region_order) +
        scale_shape_manual(name = "SubMap Region", values = shapes, limits = region_order, breaks = region_order) +
        scale_y_continuous(name = "Cross Correlation Function (CCF)", limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
        scale_x_discrete(name = "SubMap Transects", labels = label_map, expand = expansion(add = c(2, 2))) +
        coord_cartesian(clip = "off") +
        facet_wrap(~method_pair, scales = "free_x", ncol = 1) +
        theme_facet(base_size, show_legend = TRUE) +
        theme(
          axis.text.x = element_blank(),
          axis.title.x = element_text(margin = margin(20, 0, 0, 0)),
          axis.ticks.x = element_line(linewidth = 0.3),
          axis.ticks.length.x = unit(-0.3, "lines"),
          legend.title.position = "top",
          legend.direction = "vertical",
          legend.box = "horizontal"
        )

      ggsave(file = out_path, plot = p, width = 14, height = 12, dpi = 300, bg = "white")
    },
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_case_study_composition <- function(data_dir,
                                        out_dir,
                                        fig_dir,
                                        submap_transect_set,
                                        strip_idx = 1,
                                        n_grid = 500,
                                        max_lag_frac = 0.1,
                                        max_matrix_size = 5,
                                        base_size = 16,
                                        quiet = FALSE) {
  with_error_handling(
    {
      suppressWarnings(suppressMessages({
        result <- load_set_data(data_dir = data_dir, out_dir = out_dir, submap_transect_set = submap_transect_set, quiet = quiet)
        if (is.null(result)) {
          return(invisible())
        }

        if (!"submap_transect_id" %in% names(result$data)) {
          if (!quiet) cat(" !! Error: 'submap_transect_id' column missing in dataset '", submap_transect_set, "'\n", sep = "")
          return(invisible())
        }

        all_case_ids <- sort(unique(result$data$submap_transect_id))
        n_ids <- length(all_case_ids)

        if (n_ids <= max_matrix_size) {
          id_chunks <- list(all_case_ids)
        } else {
          id_chunks <- tibble(case_id = all_case_ids) |>
            mutate(chunk_idx = ceiling(row_number() / max_matrix_size)) |>
            nest(data = case_id) |>
            pull(data) |>
            map(\(df) df$case_id)

          if (length(id_chunks) > 1 && length(id_chunks[[length(id_chunks)]]) == 1) {
            last_idx <- length(id_chunks)
            id_chunks[[last_idx - 1]] <- c(id_chunks[[last_idx - 1]], id_chunks[[last_idx]])
            id_chunks[[last_idx]] <- NULL # Delete the broken 1-item chunk
          }
        }

        iwalk(id_chunks, function(case_ids, chunk_idx) {
          set_label <- paste(case_ids, collapse = "-")

          suffix <- if (length(id_chunks) > 1) paste0("-part", chunk_idx) else ""
          out_path <- file.path(fig_dir, "submap", "case_study", paste0(set_label, "-composition", suffix, ".png"))

          if (check_plot_path(out_path)) {
            return(invisible())
          }

          data <- result$data[result$data$submap_transect_id %in% case_ids, ]
          variogram_model <- result$variogram_model

          if (nrow(data) < 2) {
            if (!quiet) cat(" !! Warning: fewer than 2 matching transects for case study chunk '", set_label, "'\n", sep = "")
            return(invisible())
          }

          ids <- intersect(case_ids, data$submap_transect_id)
          data <- data[match(ids, data$submap_transect_id), ]

          obs_trend_col <- paste0("ihfc2024_obs_strip", strip_idx, "_trend")
          sim_trend_col <- paste0("lucazeau2019_sim_strip", strip_idx, "_trend")
          krg_trend_col <- make_krg_trend_col(variogram_model, strip_idx)

          obs_proj_col <- paste0("ihfc2024_obs_strip", strip_idx, "_proj")
          sim_proj_col <- paste0("lucazeau2019_sim_strip", strip_idx, "_proj")
          krg_proj_col <- make_krg_proj_col(variogram_model, strip_idx)

          build_ccf_matrix <- function(vecs) {
            n <- length(vecs)

            grid <- expand_grid(i = seq_len(n), j = seq_len(n)) |>
              filter(i <= j) |>
              mutate(val = map2_dbl(i, j, \(idx_i, idx_j) {
                if (idx_i == idx_j) {
                  return(1.0)
                }
                max_ccf_value(vecs[[idx_i]], vecs[[idx_j]], max_lag_frac, n_grid)
              }))

            mat <- matrix(NA_real_, n, n, dimnames = list(ids, ids))
            mat[as.matrix(grid[, c("i", "j")])] <- grid$val

            mat
          }

          mat_obs <- build_ccf_matrix(get_method_trend_vectors(data, obs_trend_col, n_grid))
          mat_sim <- build_ccf_matrix(get_method_trend_vectors(data, sim_trend_col, n_grid))
          mat_krg <- build_ccf_matrix(get_method_trend_vectors(data, krg_trend_col, n_grid))

          mat_to_long <- function(mat, method_label) {
            as_tibble(mat, rownames = "id1") |>
              pivot_longer(-id1, names_to = "id2", values_to = "max_ccf") |>
              mutate(method = method_label, id1 = factor(id1, levels = ids), id2 = factor(id2, levels = rev(ids)))
          }

          heatmap_df <- bind_rows(
            mat_to_long(mat_obs, "Observations"),
            mat_to_long(mat_krg, "Krige"),
            mat_to_long(mat_sim, "Similarity")
          ) |>
            mutate(method = factor(method, levels = c("Observations", "Krige", "Similarity")))

          heatmap_df <- heatmap_df |>
            mutate(cell_label = ifelse(as.character(id1) != as.character(id2) & !is.na(max_ccf), sprintf("%.2f", max_ccf), ""))

          p_heat <- ggplot(heatmap_df, aes(id1, id2, fill = max_ccf)) +
            geom_tile() +
            geom_text(aes(label = cell_label), size = base_size * 0.28, color = "black") +
            scale_fill_distiller(palette = "RdYlGn", direction = 1, limits = c(-1, 1), na.value = "grey90") +
            facet_wrap(~method, nrow = 1) +
            labs(x = NULL, y = NULL) +
            coord_cartesian(expand = c(0, 0)) +
            theme_facet(base_size) +
            theme(axis.text = element_text(angle = 45, hjust = 1))

          collect_proj_pts <- function(proj_col, method_label) {
            if (is.null(proj_col) || !proj_col %in% names(data)) {
              return(NULL)
            }
            map2_dfr(data[[proj_col]], ids, function(shp, tid) {
              if (is.null(shp) || nrow(shp) == 0) {
                return(NULL)
              }
              shp |>
                st_set_geometry(NULL) |>
                mutate(transect_id = tid, method = method_label)
            })
          }

          collect_trend_lines <- function(trend_col, method_label) {
            if (is.null(trend_col) || !trend_col %in% names(data)) {
              return(NULL)
            }
            map2_dfr(data[[trend_col]], ids, function(tr, tid) {
              if (is.null(tr) || nrow(tr) == 0) {
                return(NULL)
              }
              tr |>
                filter(method == "mean") |>
                mutate(transect_id = tid, method = method_label)
            })
          }

          pts_df <- bind_rows(
            collect_proj_pts(obs_proj_col, "Observations"),
            collect_proj_pts(krg_proj_col, "Krige"),
            collect_proj_pts(sim_proj_col, "Similarity")
          ) |> mutate(method = factor(method, levels = c("Observations", "Krige", "Similarity")))

          trends_df <- bind_rows(
            collect_trend_lines(obs_trend_col, "Observations"),
            collect_trend_lines(krg_trend_col, "Krige"),
            collect_trend_lines(sim_trend_col, "Similarity")
          ) |> mutate(method = factor(method, levels = c("Observations", "Krige", "Similarity")))

          p_prof <- ggplot() +
            geom_point(data = pts_df, aes(distance, hf, color = transect_id), shape = 16, size = 1.8, alpha = 0.20, na.rm = TRUE) +
            geom_line(data = trends_df, aes(distance, hf, color = transect_id, group = transect_id), linewidth = 1.1, na.rm = TRUE) +
            scale_color_d3(palette = "category10") +
            guides(color = guide_legend(nrow = 1)) +
            scale_y_continuous(limits = c(0, 250), breaks = seq(0, 250, 50)) +
            scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
            labs(x = "Normalized Distance", y = bquote("Q" ~ (mWm^-2)), color = NULL) +
            facet_wrap(~method, nrow = 1) +
            theme_facet(base_size, show_legend = TRUE) +
            theme(strip.text = element_blank(), axis.title.y = element_text(margin = margin(0, -15, 0, 0)))

          p_comp <- p_heat / p_prof + plot_layout(heights = c(1.8, 1.0))

          ggsave(file = out_path, plot = p_comp, width = 10.0, height = 6.0, dpi = 300, bg = "white")
        })
      }))
    },
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_global_dataset_summary <- function(data_dir, out_dir, fig_dir, base_size = 12, quiet = FALSE) {
  out_path <- file.path(fig_dir, "summary", paste0("global-dataset.png"))

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
          basemap_layers(data_dir, out_dir, proj, minor_grats = FALSE, quiet = quiet) +
          geom_sf(data = sf_hull, fill = "black", color = "black", linewidth = 0.3, alpha = 0.1, show.legend = FALSE) +
          geom_sf(data = sf_submap, color = "black", linewidth = 0.6) +
          theme_globe(proj, base_size)

        p2 <- ggplot() +
          basemap_layers(data_dir, out_dir, proj, minor_grats = FALSE, quiet = quiet) +
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
    p <- row1 / row2 + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, -50))
    ggsave(file = out_path, plot = p, width = 6.5, height = 5.5, dpi = 300, bg = "white")
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_nlopt_summary <- function(fig_dir, nlopt_summary, base_size = 14, quiet = FALSE) {
  if (is.null(nlopt_summary) || nrow(nlopt_summary) == 0) {
    return(invisible())
  }

  with_error_handling(
    {
      out_path <- file.path(fig_dir, "summary", "nlopt.png")
      if (check_plot_path(out_path)) {
        return(invisible())
      }

      suppressWarnings({
        suppressMessages({
          df <- nlopt_summary |>
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
              "Max Size" = max_point_pair_distance,
              "Vgram Cutoff" = variogram_cutoff,
              "Iteration" = itr,
              "Vgram Lags" = n_variogram_lags,
              "Observations" = n_heatflow_obs,
              "Max Obs" = max_n_point_pairs
            ) |>
            pivot_longer(-c(submap_transect_set, variogram_model, total_cost)) |>
            mutate(region = str_extract(submap_transect_set, "^[^_]+"))

          region_order <- df |>
            distinct(region) |>
            pull(region) |>
            sort() |>
            as.character()

          region_palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
          region_shapes <- c(21, 22, 23, 24)

          cols <- setNames(region_palette[seq_along(region_order)], region_order)
          shapes <- setNames(region_shapes[seq_along(region_order)], region_order)

          p <- df |>
            mutate(region = factor(region, levels = region_order)) |>
            ggplot(aes(value, total_cost, fill = region, shape = region)) +
            geom_point(size = 2.5, color = "black") +
            scale_fill_manual(name = "Region", values = cols, limits = region_order, breaks = region_order) +
            scale_shape_manual(name = "Region", values = shapes, limits = region_order, breaks = region_order) +
            labs(x = NULL, y = "Cost", fill = NULL, shape = NULL) +
            scale_x_continuous(breaks = pretty_breaks(n = 4), guide = guide_axis(check.overlap = TRUE), expand = expansion(mult = 0.10)) +
            scale_y_continuous(expand = expansion(mult = 0.10)) +
            guides(fill = guide_legend(ncol = 4, override.aes = list(size = 3.5)), shape = guide_legend(ncol = 4)) +
            facet_wrap(~name, scales = "free_x") +
            theme_facet(base_size, show_legend = TRUE) +
            theme(strip.text = element_text(size = base_size * 1.0))

          ggsave(file = out_path, plot = p, width = 6.5, height = 5.0, dpi = 300, bg = "white")
        })
      })
    },
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_point_by_point_summary <- function(data_dir, out_dir, fig_dir, submap_transect_sets, base_size = 26, quiet = FALSE) {
  with_error_handling(
    {
      compile_set_distributions <- function(submap_transect_set) {
        with_error_handling(
          {
            cache_path <- file.path(out_dir, "submap", paste0(submap_transect_set, ".RData"))
            if (!file.exists(cache_path)) {
              return(NULL)
            }

            env <- new.env()
            load(cache_path, envir = env)

            if (!"df" %in% names(env) || !"krg_sim_diff_hull" %in% names(env$df)) {
              return(NULL)
            }

            df_hull <- env$df$krg_sim_diff_hull[[1]]
            grid_data <- tibble(submap_transect_set = submap_transect_set)

            est_list <- if ("sim_krg_diff_est" %in% names(df_hull)) {
              df_hull$sim_krg_diff_est[!is.na(df_hull$sim_krg_diff_est)]
            } else {
              numeric(0)
            }

            sigma_list <- if ("sim_krg_diff_sigma" %in% names(df_hull)) {
              df_hull$sim_krg_diff_sigma[!is.na(df_hull$sim_krg_diff_sigma)]
            } else {
              numeric(0)
            }

            df_grid_long <- bind_rows(
              tibble(submap_transect_set = submap_transect_set, metric = "Sim - Krg (est)", value = est_list),
              tibble(submap_transect_set = submap_transect_set, metric = "Sim - Krg (std)", value = sigma_list)
            )

            opt_data <- get_optimal_krige_model(out_dir, submap_transect_set, quiet = quiet)
            variogram_model <- opt_data$optimal_krige_model$variogram_model
            fitted_variogram <- opt_data$fitted_variogram
            max_n_point_pairs <- opt_data$optimal_krige_model$max_n_point_pairs
            max_point_pair_distance <- opt_data$optimal_krige_model$max_point_pair_distance

            point_pairs_sim <- get_closest_interp_obs_point_pairs(data_dir, out_dir, submap_transect_set, quiet = TRUE)
            point_pairs_krg <- get_closest_interp_obs_point_pairs(
              data_dir, out_dir, submap_transect_set,
              variogram_model = variogram_model, fitted_variogram = fitted_variogram,
              max_n_point_pairs = max_n_point_pairs, max_point_pair_distance = max_point_pair_distance,
              quiet = quiet
            )

            sim_res <- if (!is.null(point_pairs_sim) && all(c("lucazeau2019_sim_est", "ihfc2024_obs") %in% names(point_pairs_sim))) {
              (point_pairs_sim$lucazeau2019_sim_est - point_pairs_sim$ihfc2024_obs) |>
                na.omit() |>
                as.numeric()
            } else {
              numeric(0)
            }

            krg_res <- if (!is.null(point_pairs_krg) && all(c("kerswell2025_krg_est", "ihfc2024_obs") %in% names(point_pairs_krg))) {
              (point_pairs_krg$kerswell2025_krg_est - point_pairs_krg$ihfc2024_obs) |>
                na.omit() |>
                as.numeric()
            } else {
              numeric(0)
            }

            df_acc_long <- bind_rows(
              tibble(submap_transect_set = submap_transect_set, metric = "Sim - Obs (est)", value = sim_res),
              tibble(submap_transect_set = submap_transect_set, metric = "Krg - Obs (est)", value = krg_res)
            )

            return(bind_rows(df_grid_long, df_acc_long))
          },
          default = NULL,
          quiet = quiet
        )
      }

      combined_data <- map_df(submap_transect_sets, ~ compile_set_distributions(.x))

      if (is.null(combined_data) || nrow(combined_data) == 0) {
        return(invisible())
      }

      summary_dir <- file.path(fig_dir, "summary")
      if (!dir.exists(summary_dir)) dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

      processed_data <- combined_data |>
        mutate(
          region_group_sort = str_extract(submap_transect_set, "^[^_]+"),
          set_label = str_extract(submap_transect_set, "SET\\d+"),
          set_label = factor(set_label, levels = paste0("SET", 12:1))
        ) |>
        add_count(submap_transect_set, metric, name = "n_count") |>
        mutate(n_label = paste0("n = ", n_count))

      metric_pipeline <- tibble(
        metric_name = c("Sim - Krg (est)", "Sim - Krg (std)", "Krg - Obs (est)", "Sim - Obs (est)"),
        file_suffix = c("est.png", "sigma.png", "krg-residuals.png", "sim-residuals.png"),
        axis_title = c(
          "'Sim - Krg Prediction Differences' ~ (mWm^-2)",
          "'Sim - Krg Uncertainty Differences' ~ (mWm^-2)",
          "'Krg - Obs Residuals' ~ (mWm^-2)",
          "'Sim - Obs Residuals' ~ (mWm^-2)"
        )
      )

      max_submaps_in_facet <- processed_data |>
        distinct(region_group_sort, set_label) |>
        count(region_group_sort) |>
        pull(n) |>
        max()

      pwalk(metric_pipeline, function(metric_name, file_suffix, axis_title) {
        out_path <- file.path(summary_dir, paste0("point-by-point-summary-", file_suffix))

        if (check_plot_path(out_path)) {
          return()
        }

        if (!quiet) cat(" -> ", out_path, "\n", sep = "")

        plot_df <- processed_data |> filter(metric == metric_name)

        text_df <- plot_df |>
          distinct(set_label, region_group_sort, n_label) |>
          mutate(x_pos = -120)

        p <- ggplot(plot_df, aes(x = value, y = set_label, fill = factor(after_stat(quantile)))) +
          geom_density_ridges_gradient(calc_ecdf = TRUE, quantiles = 4, quantile_lines = TRUE, scale = 1.2, show.legend = FALSE, linewidth = 0.4) +
          geom_text(data = text_df, aes(x = x_pos, y = set_label, label = n_label), inherit.aes = FALSE, hjust = 0, vjust = 1.5, size = 7.0) +
          scale_fill_viridis_d(option = "magma") +
          scale_x_continuous(limits = c(-120, 120)) +
          facet_wrap(~region_group_sort, scales = "fixed", ncol = 4) +
          labs(x = parse(text = axis_title), y = "") +
          theme_facet(base_size = base_size)

        ggsave(file = out_path, plot = p, width = 16.0, height = 8.5, dpi = 300, bg = "white")
      })
    },
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_optimal_variogram <- function(data_dir, out_dir, fig_dir, submap_transect_set, base_size = 20, quiet = FALSE) {
  with_error_handling({
    suppressWarnings({
      suppressMessages({
        out_path <- file.path(fig_dir, "submap", "optimal", paste0(submap_transect_set, "-variogram.png"))

        if (check_plot_path(out_path)) {
          return(invisible())
        }

        opt_data <- get_optimal_krige_model(out_dir = out_dir, submap_transect_set = submap_transect_set, quiet = quiet)
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

        p0 <- nlopt_iter_df |>
          ggplot(aes(x = itr, y = max_point_pair_distance / 1e3, group = variogram_model, linetype = variogram_model)) +
          geom_path(color = "forestgreen", linewidth = 1.5) +
          scale_linetype_manual(values = linetype_map) +
          labs(x = "Iteration", y = "Max Distance (km)") +
          theme_vgrm(base_size)

        p1 <- nlopt_iter_df |>
          ggplot(aes(x = itr, y = max_n_point_pairs, group = variogram_model, linetype = variogram_model)) +
          geom_path(color = "firebrick", linewidth = 1.5) +
          scale_linetype_manual(values = linetype_map) +
          labs(x = "Iteration", y = "Max Pairs", linetype = "Variogram Model") +
          guides(linetype = guide_legend(override.aes = list(color = "black"))) +
          theme_vgrm(base_size, show_legend = TRUE)

        p2 <- nlopt_iter_df |>
          ggplot(aes(x = itr, y = variogram_cutoff, group = variogram_model, linetype = variogram_model)) +
          geom_path(color = "orchid4", linewidth = 1.5) +
          scale_linetype_manual(values = linetype_map) +
          labs(x = "Iteration", y = "Vgram Cutoff") +
          theme_vgrm(base_size)

        p3 <- nlopt_iter_df |>
          ggplot(aes(x = itr, y = n_variogram_lags, group = variogram_model, linetype = variogram_model)) +
          geom_path(color = "saddlebrown", linewidth = 1.5) +
          scale_linetype_manual(values = linetype_map) +
          labs(x = "Iteration", y = "Vgram Lags") +
          theme_vgrm(base_size)

        p4 <- nlopt_iter_df |>
          ggplot(aes(x = itr, y = total_cost, group = variogram_model, linetype = variogram_model)) +
          geom_path(color = "black", linewidth = 1.5) +
          scale_linetype_manual(values = linetype_map) +
          labs(x = "Iteration", y = "Total Cost") +
          theme_vgrm(base_size)

        p5 <- experimental_variogram |>
          ggplot(aes(x = dist / 1e3, y = sqrt(gamma))) +
          geom_line(data = fitted_variogram_line, linewidth = 1.5) +
          geom_point(shape = 20, size = 3) +
          labs(x = "Lag Distance (km)", y = bquote("Variance" ~ (mWm^-2))) +
          theme_vgrm(base_size)

        p <- (p0 + p1) / (p2 + p3) / (p4 + p5)

        ggsave(file = out_path, plot = p, width = 13, height = 10, dpi = 300, bg = "white")
      })
    })
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_submap_transect_set_composition <- function(data_dir,
                                                 out_dir,
                                                 fig_dir,
                                                 submap_transect_set,
                                                 strip_widths = c(7.5e4),
                                                 base_size = 30,
                                                 quiet = quiet) {
  with_error_handling(
    {
      suppressWarnings({
        suppressMessages({
          out_paths <- file.path(fig_dir, "submap", "composition", paste0(submap_transect_set, c("-obs.png", "-sim.png", "-krg.png", "-comp.png")))

          if (all(map_lgl(out_paths, check_plot_path))) {
            return(invisible())
          }

          result <- load_set_data(data_dir = data_dir, out_dir = out_dir, submap_transect_set = submap_transect_set, quiet = quiet)
          if (is.null(result)) {
            return(invisible())
          }
          data <- result$data
          variogram_model <- result$variogram_model

          proj <- determine_projection(sf_obj = data$hull, proj_type = "albers", quiet = quiet)
          sf_hull <- data$hull
          sf_ihfc_hull <- data$ihfc2024_obs_hull[[1]]

          sr_sim_hull <- with_error_handling(to_spatraster(data$lucazeau2019_sim_hull[[1]], data_col = "lucazeau2019_sim_est"), default = NULL)

          prefix <- "kerswell2025_krg_"
          model_str <- str_to_lower(variogram_model)
          krg_hull_str <- paste0(prefix, model_str, "_hull")
          sr_krige_hull <- with_error_handling(to_spatraster(data[[krg_hull_str]][[1]], data_col = "kerswell2025_krg_est"), default = NULL)

          strip_buffers <- build_strip_buffer_outlines(st_geometry(data), proj)

          collect_strip_pts <- function(proj_col) {
            rows <- data[[proj_col]]
            valid <- Filter(Negate(is.null), rows)
            if (length(valid) == 0) {
              return(NULL)
            }
            bind_rows(map(valid, st_set_geometry, NULL))
          }

          collect_strip_trends <- function(trend_col) {
            rows <- data[[trend_col]]
            ids <- data$submap_transect_id

            valid_indices <- which(map_lgl(rows, ~ !is.null(.x) && nrow(.x) > 0 && !all(is.na(.x$hf))))

            if (length(valid_indices) == 0) {
              return(NULL)
            }

            map(valid_indices, function(i) {
              rows[[i]] |> mutate(submap_transect_id = ids[i])
            }) |>
              bind_rows()
          }

          ihfc_proj_list <- map(seq_along(strip_widths), ~ collect_strip_pts(paste0("ihfc2024_obs_strip", .x, "_proj")))
          ihfc_trend_list <- map(seq_along(strip_widths), ~ collect_strip_trends(paste0("ihfc2024_obs_strip", .x, "_trend")))

          sim_proj_list <- map(seq_along(strip_widths), ~ collect_strip_pts(paste0("lucazeau2019_sim_strip", .x, "_proj")))
          sim_trend_list <- map(seq_along(strip_widths), ~ collect_strip_trends(paste0("lucazeau2019_sim_strip", .x, "_trend")))

          krg_proj_list <- map(seq_along(strip_widths), ~ {
            col <- paste0(prefix, model_str, "_strip", .x, "_proj")
            if (col %in% names(data)) collect_strip_pts(col) else NULL
          })
          krg_trend_list <- map(seq_along(strip_widths), ~ {
            col <- paste0(prefix, model_str, "_strip", .x, "_trend")
            if (col %in% names(data)) collect_strip_trends(col) else NULL
          })

          strip_overlay_black <- map(strip_buffers, ~ geom_sf(data = .x, fill = NA, linewidth = 0.4, color = "black"))
          strip_overlay_white <- map(strip_buffers, ~ geom_sf(data = .x, fill = NA, linewidth = 0.4, color = "white", alpha = 0.2))

          p1 <- ggplot() +
            basemap_layers(data_dir, out_dir, proj) +
            geom_sf(data = sf_ihfc_hull, aes(color = ihfc2024_obs), shape = 20, size = 4) +
            strip_overlay_black +
            ggtitle("Observations") +
            viridis_scale_color() +
            theme_transect_set(sf_hull, proj, base_size, show_legend = FALSE)

          p3 <- ggplot() +
            basemap_layers(data_dir, out_dir, proj) +
            geom_spatraster(data = sr_sim_hull) +
            strip_overlay_white +
            ggtitle("Similarity") +
            viridis_scale_color() +
            viridis_scale_fill() +
            theme_transect_set(sf_hull, proj, base_size, show_legend = FALSE)

          if (is.null(sr_krige_hull)) {
            p5 <- ggplot() +
              basemap_layers(data_dir, out_dir, proj) +
              geom_sf(data = sf_hull, linewidth = 1, fill = NA, color = "black") +
              strip_overlay_white +
              ggtitle("Krige (failed)") +
              theme_transect_set(sf_hull, proj, base_size, show_legend = FALSE)
          } else {
            p5 <- ggplot() +
              basemap_layers(data_dir, out_dir, proj) +
              geom_spatraster(data = sr_krige_hull) +
              strip_overlay_white +
              ggtitle("Krige") +
              viridis_scale_color() +
              viridis_scale_fill() +
              theme_transect_set(sf_hull, proj, base_size, show_legend = FALSE)
          }

          no_y_lab <- list(labs(y = NULL))
          no_y_axis <- theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

          p2 <- build_hf_profile_plot(ihfc_proj_list, ihfc_trend_list, base_size = base_size, show_legend = TRUE)
          p4 <- build_hf_profile_plot(sim_proj_list, sim_trend_list, base_size = base_size, show_legend = FALSE) + no_y_lab
          p6 <- build_hf_profile_plot(krg_proj_list, krg_trend_list, base_size = base_size, show_legend = FALSE) + no_y_lab

          pp1 <- (p1 / p2) + plot_layout(widths = 1, heights = c(1.8, 1))
          ggsave(file = out_paths[1], plot = pp1, width = 6.5, height = 10, dpi = 300, bg = "white")

          pp2 <- (p3 / p4) + plot_layout(widths = 1, heights = c(1.8, 1))
          ggsave(file = out_paths[2], plot = pp2, width = 6.5, height = 10, dpi = 300, bg = "white")

          pp3 <- (p5 / p6) + plot_layout(widths = 1, heights = c(1.8, 1))
          ggsave(file = out_paths[3], plot = pp3, width = 6.5, height = 10, dpi = 300, bg = "white")

          comp_design <- "
            ABC
            DEF
          "

          row_spacer <- theme(plot.margin = margin(25, 5, 5, 5))

          p_comp <-
            p1 + (p5 + no_y_axis) + (p3 + no_y_axis) +
            p2 + (p6 + no_y_axis + row_spacer) + (p4 + no_y_axis + row_spacer) +
            plot_layout(design = comp_design, heights = c(1.8, 1), widths = c(1, 1, 1))

          ggsave(file = out_paths[4], plot = p_comp, width = 19.5, height = 11, dpi = 300, bg = "white")
        })
      })
    },
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_submap_transect_neighbors_composition <- function(data_dir,
                                                       out_dir,
                                                       fig_dir,
                                                       neighbor_sets,
                                                       strip_widths = c(7.5e4),
                                                       base_size = 32,
                                                       quiet = quiet) {
  with_error_handling(
    {
      walk(neighbor_sets, function(ids) {
        with_error_handling(
          {
            suppressWarnings({
              suppressMessages({
                if (length(ids) != 3) {
                  stop("Each neighbor set must contain exactly 3 submap_transect_set IDs")
                }

                set_label <- paste0(ids[1], "-", ids[2], "-", ids[3])
                out_paths <- file.path(
                  fig_dir, "submap", "neighbor",
                  paste0(set_label, c("-neighbor-ghf.png", "-neighbor-sim.png", "-neighbor-krg.png"))
                )

                if (all(map_lgl(out_paths, check_plot_path))) {
                  return(invisible())
                }

                datasets <- map(ids, function(id) {
                  with_error_handling(
                    {
                      result <- load_set_data(data_dir = data_dir, out_dir = out_dir, submap_transect_set = id, quiet = quiet)
                      if (is.null(result)) {
                        return(invisible())
                      }
                      data <- result$data
                      variogram_model <- result$variogram_model

                      list(id = id, data = data, variogram_model = variogram_model)
                    },
                    default = NULL,
                    quiet = quiet
                  )
                })

                build_col <- function(ds, method, show_legend = FALSE) {
                  if (is.null(ds)) {
                    blank <- ggplot() +
                      theme_void()
                    return(list(map = blank, profile = blank))
                  }

                  data <- ds$data
                  variogram_model <- ds$variogram_model
                  proj <- determine_projection(sf_obj = data$hull, proj_type = "albers", quiet = quiet)
                  sf_hull <- data$hull

                  strip_buffers <- build_strip_buffer_outlines(st_geometry(data), proj)

                  strip_overlay_black <- map(strip_buffers, ~ geom_sf(data = .x, fill = NA, linewidth = 0.4, color = "black"))
                  strip_overlay_white <- map(strip_buffers, ~ geom_sf(data = .x, fill = NA, linewidth = 0.4, color = "white", alpha = 0.4))

                  collect_strip_pts <- function(proj_col) {
                    rows <- data[[proj_col]]
                    valid <- Filter(Negate(is.null), rows)
                    if (length(valid) == 0) {
                      return(NULL)
                    }
                    bind_rows(map(valid, st_set_geometry, NULL))
                  }

                  collect_strip_trends <- function(trend_col) {
                    rows <- data[[trend_col]]
                    ids <- data$submap_transect_id

                    valid_indices <- which(map_lgl(rows, ~ !is.null(.x) && nrow(.x) > 0 && !all(is.na(.x$hf))))

                    if (length(valid_indices) == 0) {
                      return(NULL)
                    }

                    map(valid_indices, function(i) {
                      rows[[i]] |> mutate(submap_transect_id = ids[i])
                    }) |>
                      bind_rows()
                  }

                  if (method == "obs") {
                    sf_ihfc_hull <- data$ihfc2024_obs_hull[[1]]
                    proj_list <- map(seq_along(strip_widths), ~ collect_strip_pts(paste0("ihfc2024_obs_strip", .x, "_proj")))
                    trend_list <- map(seq_along(strip_widths), ~ collect_strip_trends(paste0("ihfc2024_obs_strip", .x, "_trend")))
                    p_map <- ggplot() +
                      basemap_layers(data_dir, out_dir, proj) +
                      geom_sf(data = sf_ihfc_hull, aes(color = ihfc2024_obs), shape = 20, size = 4) +
                      strip_overlay_black +
                      ggtitle(paste0(ds$id, "\nObservations")) +
                      viridis_scale_color() +
                      theme_transect_set(sf_hull, proj, base_size, show_legend = FALSE)
                  } else if (method == "sim") {
                    sr_sim <- with_error_handling(to_spatraster(data$lucazeau2019_sim_hull[[1]], data_col = "lucazeau2019_sim_est"), default = NULL)
                    proj_list <- map(seq_along(strip_widths), ~ collect_strip_pts(paste0("lucazeau2019_sim_strip", .x, "_proj")))
                    trend_list <- map(seq_along(strip_widths), ~ collect_strip_trends(paste0("lucazeau2019_sim_strip", .x, "_trend")))
                    p_map <- ggplot() +
                      basemap_layers(data_dir, out_dir, proj) +
                      (if (!is.null(sr_sim)) geom_spatraster(data = sr_sim) else NULL) +
                      strip_overlay_white +
                      ggtitle(paste0(ds$id, "\nSimilarity")) +
                      viridis_scale_color() +
                      viridis_scale_fill() +
                      theme_transect_set(sf_hull, proj, base_size, show_legend = FALSE)
                  } else {
                    prefix <- "kerswell2025_krg_"
                    model_str <- str_to_lower(variogram_model %||% "")
                    krg_hull_str <- paste0(prefix, model_str, "_hull")
                    sr_krg <- with_error_handling(to_spatraster(data[[krg_hull_str]][[1]], data_col = "kerswell2025_krg_est"), default = NULL)
                    krg_s_proj <- paste0(prefix, model_str, "_strip", seq_along(strip_widths), "_proj")
                    krg_s_trend <- paste0(prefix, model_str, "_strip", seq_along(strip_widths), "_trend")
                    proj_list <- map(krg_s_proj, ~ if (.x %in% names(data)) collect_strip_pts(.x) else NULL)
                    trend_list <- map(krg_s_trend, ~ if (.x %in% names(data)) collect_strip_trends(.x) else NULL)
                    title_str <- if (is.null(sr_krg)) {
                      "Krige (failed)"
                    } else {
                      "Krige"
                    }
                    p_map <- ggplot() +
                      basemap_layers(data_dir, out_dir, proj) +
                      (if (!is.null(sr_krg)) {
                        geom_spatraster(data = sr_krg)
                      } else {
                        geom_sf(data = sf_hull, linewidth = 1, fill = NA, color = "black")
                      }) +
                      strip_overlay_white +
                      ggtitle(paste0(ds$id, "\n", title_str)) +
                      viridis_scale_color() +
                      viridis_scale_fill() +
                      theme_transect_set(sf_hull, proj, base_size, show_legend = FALSE)
                  }

                  p_profile <- build_hf_profile_plot(proj_list, trend_list, base_size = base_size, show_legend = show_legend)

                  list(map = p_map, profile = p_profile)
                }

                walk2(
                  c("obs", "sim", "krg"),
                  out_paths,
                  function(method, out_path) {
                    cols <- map2(datasets, c(TRUE, FALSE, FALSE), ~ build_col(.x, method, show_legend = .y))

                    p_comp <- (cols[[1]]$map + cols[[2]]$map + cols[[3]]$map) /
                      (cols[[1]]$profile + cols[[2]]$profile + cols[[3]]$profile) +
                      plot_layout(widths = 1, heights = c(1.8, 1))

                    ggsave(file = out_path, plot = p_comp, width = 19.5, height = 11, dpi = 300, bg = "white")
                  }
                )
              })
            })
          },
          quiet = quiet
        )
      })
    },
    quiet = quiet
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_submap_transect_sets <- function(data_dir,
                                      out_dir,
                                      fig_dir,
                                      submap_transect_sets,
                                      strip_widths = c(7.5e4),
                                      seed = 42,
                                      nprocs = NULL,
                                      parallel = FALSE,
                                      quiet = FALSE) {
  if (is.null(nprocs)) nprocs <- availableCores() - 2

  f <- function(data_dir, out_dir, fig_dir, submap_transect_set) {
    with_error_handling(
      {
        suppressWarnings({
          suppressMessages({
            draw_submap_transect_set_composition(data_dir, out_dir, fig_dir, submap_transect_set, strip_widths, quiet = quiet)
            draw_optimal_variogram(data_dir, out_dir, fig_dir, submap_transect_set, quiet = quiet)
          })
        })
      },
      default = NULL,
      quiet = quiet
    )
  }

  with_error_handling(
    {
      precompile_set_data(data_dir = data_dir, out_dir = out_dir, submap_transect_sets = submap_transect_sets, quiet = quiet)

      if (parallel) {
        set.seed(seed)
        plan(multisession, workers = nprocs)
        future_walk(submap_transect_sets,
          ~ f(data_dir = data_dir, out_dir = out_dir, fig_dir = fig_dir, submap_transect_set = .x),
          .options = furrr_options(seed = seed)
        )
      } else {
        walk(submap_transect_sets, ~ f(data_dir = data_dir, out_dir = out_dir, fig_dir = fig_dir, submap_transect_set = .x))
      }

      n_sets <- length(submap_transect_sets)
      if (n_sets >= 3) {
        neighbor_sets <- map(seq_len(n_sets - 2), ~ submap_transect_sets[.x:(.x + 2)])
        draw_submap_transect_neighbors_composition(
          data_dir = data_dir, out_dir = out_dir, fig_dir = fig_dir, neighbor_sets = neighbor_sets, strip_widths = strip_widths, quiet = quiet
        )
      }
    },
    quiet = quiet
  )
}
