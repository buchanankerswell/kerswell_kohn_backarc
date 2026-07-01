#######################################################
## Load libraries                                    ##
#######################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(marmap)
  library(cowplot)
  library(ggnewscale)
  library(patchwork)
})

#######################################################
## Helper functions                                  ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load_map_data <- function(in_path) {
  load(in_path, envir = parent.frame())

  if (!exists("shp_submap", envir = parent.frame())) {
    stop(" !! ERROR: missing map data!")
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
check_plot_exists <- function(out_path) {
  if (file.exists(out_path)) {
    cat(" -- Found plot: ", out_path, "\n", sep = "")
    TRUE
  } else {
    cat(" -> Plotting: ", out_path, "\n", sep = "")
    FALSE
  }
}

#######################################################
## Plotting functions                                ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_ghf <- function(in_path, fig_dir) {
  load_map_data(in_path)

  out_path <- file.path(fig_dir, "ghf.png")

  if (check_plot_exists(out_path)) {
    return(invisible())
  }

  if (!dir.exists(fig_dir)) {
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  }

  suppressWarnings({
    suppressMessages({
      p1 <-
        ggplot() +
        geom_sf(data = shp_relief_world, aes(color = elev), shape = 15, size = 0.01) +
        scale_color_etopo(guide = "none") +
        geom_sf(data = shp_ridge, linewidth = 0.5, color = "white") +
        geom_sf(data = shp_transform, linewidth = 0.5, color = "white") +
        geom_sf(data = shp_trench, linewidth = 0.5, color = "white") +
        geom_sf(data = shp_submap, color = "black") +
        ggtitle("Submap Transects") +
        coord_sf(expand = FALSE, lims_method = "geometry_bbox") +
        theme_map(font_size = 14)

      p2 <-
        ggplot() +
        geom_sf(data = shp_relief_world, aes(color = elev), shape = 15, size = 0.01) +
        scale_color_etopo(guide = "none") +
        new_scale_color() +
        geom_sf(data = shp_ridge, linewidth = 0.5, color = "white") +
        geom_sf(data = shp_transform, linewidth = 0.5, color = "white") +
        geom_sf(data = shp_trench, linewidth = 0.5, color = "white") +
        geom_sf(data = shp_ghf_raw_wgs, aes(color = obs), size = 0.5, shape = 20) +
        scale_color_viridis_c(
          option = "magma",
          name = bquote("Q" ~ (mWm^-2)),
          limits = c(0, 250),
          breaks = c(0, 125, 250),
          na.value = "transparent",
          guide = guide_colorbar(
            title.vjust = 1, show.limits = TRUE,
            frame.colour = "black",
            ticks.colour = "black"
          )
        ) +
        ggtitle("Global Heatflow Observations") +
        coord_sf(expand = FALSE, lims_method = "geometry_bbox") +
        theme_map(font_size = 14)

      p <- p1 / p2 + plot_annotation(tag_levels = "a") &
        theme(
          axis.text = element_blank(),
          panel.grid = element_line(linewidth = 0.05, color = "grey20"),
          plot.title = element_text(vjust = 0, hjust = 0.5, margin = margin(10, 10, 10, 10)),
          plot.margin = margin(),
          plot.tag = element_text(face = "bold", size = 20, margin = margin(15, 0, 0, 10)),
          legend.position = "bottom",
          legend.justification = "center",
          legend.direction = "horizontal",
          legend.margin = margin(),
          legend.box.margin = margin(5, 5, 5, 5),
          legend.key.height = unit(0.5, "cm"),
          legend.key.width = unit(1, "cm"),
          legend.title = element_text(vjust = 0, color = "black", size = 14)
        )

      ggsave(file = out_path, plot = p, width = 6.5, height = 8.0, dpi = 300, bg = "white")
    })
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot transect buff comp !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_transect_buff_comp <- function(trans_id, base_size = 20) {
  fig_dir <- "figs/transect_buff/"
  fig_path1 <- paste0(fig_dir, trans_id, "-transect-buff-ghf.png")
  fig_path2 <- paste0(fig_dir, trans_id, "-transect-buff-sim.png")
  fig_path3 <- paste0(fig_dir, trans_id, "-transect-buff-krg.png")
  fig_path_comp <- paste0(fig_dir, trans_id, "-transect-buff-comp.png")
  if (!dir.exists(fig_dir)) {
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (file.exists(fig_path_comp) && file.exists(fig_path1) && file.exists(fig_path2) && file.exists(fig_path3)) {
    cat("\n   ", fig_path1, " already exists ...", sep = "")
    cat("\n   ", fig_path2, " already exists ...", sep = "")
    cat("\n   ", fig_path3, " already exists ...", sep = "")
    cat("\n   ", fig_path_comp, " already exists ...", sep = "")
    return(invisible())
  }
  tryCatch(
    {
      cat("\n   Plotting: ", fig_path1, " ...", sep = "")
      cat("\n   Plotting: ", fig_path2, " ...", sep = "")
      cat("\n   Plotting: ", fig_path3, " ...", sep = "")
      cat("\n   Plotting: ", fig_path_comp, " ...", sep = "")
      km <- get_optimal_krige_model(trans_id)
      fv <- km$fitted_vgrm
      np <- km$opt_krige_mod_summary$n_pairs
      x <- compile_trans_data(trans_id, fv = fv, np = np)
      comps_sim <- get_closest_interp_obs(trans_id)
      comps_krg <- get_closest_interp_obs(trans_id, fv = fv, np = np)
      hf_pt_size <- 3
      pt_stroke <- 0.8
      p_title <- paste0("Submap Transect: ", trans_id, " ", x$trench_name)
      map_theme <- list(
        theme_bw(base_size = base_size),
        theme(
          plot.margin = margin(5, 5, 5, 5),
          plot.title = element_text(vjust = 0, hjust = 0.5, margin = margin(0, 5, 10, 5))
        ),
        annotation_scale(
          location = "bl", width_hint = 0.33, text_cex = 1.6, style = "ticks",
          line_width = 4, text_face = "bold"
        ),
        annotation_north_arrow(
          location = "bl", which_north = "true", pad_x = unit(0.0, "cm"),
          height = unit(2, "cm"), width = unit(2, "cm"),
          pad_y = unit(0.5, "cm"), style = north_arrow_fancy_orienteering
        )
      )
      profile_theme <-
        list(
          theme_bw(base_size = base_size),
          theme(
            panel.grid = element_blank(), panel.background = element_rect(fill = "grey90"),
            plot.margin = margin(5, 5, 5, 5),
            legend.justification = "right", legend.position = "inside",
            legend.position.inside = c(0.92, 0.85), legend.direction = "horizontal",
            legend.key.height = unit(0.5, "cm"), legend.key.width = unit(1, "cm"),
            legend.box.margin = margin(2, 2, 2, 2), legend.margin = margin(),
            legend.title = element_text(vjust = 0, size = base_size),
            legend.background = element_blank()
          )
        )
      suppressWarnings({
        suppressMessages({
          p1 <-
            ggplot(x) +
            geom_sf(data = x$bathy[[1]], aes(color = elev), size = 0.5, shape = 15) +
            scale_color_etopo(guide = "none") +
            new_scale_color() +
            geom_sf(aes(geometry = large_buffer), fill = NA, linewidth = 0.5) +
            geom_sf(aes(geometry = ridge), color = "white") +
            geom_sf(aes(geometry = transform), color = "white") +
            geom_sf(aes(geometry = trench), color = "white", linewidth = 1.5) +
            geom_sf(aes(geometry = transect), color = "black", linewidth = 1.5) +
            geom_sf(aes(geometry = volcano), color = "black", fill = "white", shape = 24) +
            geom_sf(data = x$ghf_large_buff[[1]], aes(color = obs), shape = 20, size = hf_pt_size) +
            geom_sf(aes(geometry = small_buffer1), fill = NA, linewidth = 0.5, color = "black") +
            geom_sf(aes(geometry = small_buffer2), fill = NA, linewidth = 0.5, color = "black") +
            geom_sf(aes(geometry = small_buffer3), fill = NA, linewidth = 0.5, color = "black") +
            ggtitle("Global HF Observations") +
            scale_color_viridis_c(
              option = "magma", name = bquote("Q" ~ (mWm^-2)), limits = c(0, 250),
              breaks = c(0, 125, 250), na.value = "transparent", guide = "none"
            ) +
            coord_sf(expand = FALSE, lims_method = "geometry_bbox") +
            map_theme
          if (is.null(x$ghf_projected1[[1]])) {
            p2 <- ggplot(data = data.frame(), aes())
          } else {
            p2 <-
              ggplot() +
              geom_point(
                data = x$ghf_projected1[[1]], aes(projected_distances, obs), shape = 20,
                color = "grey20", size = 2
              ) +
              geom_point(
                data = x$ghf_projected2[[1]], aes(projected_distances, obs), shape = 20,
                color = "grey20", size = 2
              ) +
              geom_point(
                data = x$ghf_projected3[[1]], aes(projected_distances, obs), shape = 20,
                color = "grey20", size = 2
              ) +
              geom_path(data = x$ghf_loess1[[1]], aes(projected_distances, obs), color = "black") +
              geom_path(data = x$ghf_loess2[[1]], aes(projected_distances, obs), color = "black") +
              geom_path(data = x$ghf_loess3[[1]], aes(projected_distances, obs), color = "black") +
              geom_rug(
                data = x$ghf_loess3[[1]], aes(projected_distances, obs, color = obs),
                sides = "b", length = unit(0.06, "npc")
              )
          }
          p2 <- p2 +
            labs(x = "Normalized Distance", y = NULL) +
            scale_color_viridis_c(
              option = "magma", name = bquote("Q" ~ (mWm^-2)), limits = c(0, 250),
              breaks = c(0, 125, 250), na.value = "transparent",
              guide = guide_colorbar(
                title.vjust = 1, show.limits = TRUE,
                frame.colour = "black",
                ticks.colour = "black"
              )
            ) +
            scale_y_continuous(limits = c(0, 250), breaks = seq(0, 250, 50)) +
            scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) + profile_theme
          p <- (p1 / p2) + plot_layout(widths = 1, heights = c(1.5, 1)) &
            theme(plot.title = element_text(size = base_size * 1.5, margin = margin()))
          ggsave(file = fig_path1, plot = p, width = 6.5, height = 10, dpi = 300, bg = "white")
          p3 <-
            ggplot(x) +
            geom_sf(data = x$bathy[[1]], aes(color = elev), size = 0.5, shape = 15) +
            scale_color_etopo(guide = "none") +
            new_scale_color() +
            geom_sf(aes(geometry = large_buffer), fill = NA, linewidth = 0.5) +
            geom_sf(aes(geometry = ridge), color = "white") +
            geom_sf(aes(geometry = transform), color = "white") +
            geom_sf(aes(geometry = trench), color = "white", linewidth = 1.5) +
            geom_sf(aes(geometry = transect), color = "black", linewidth = 1.5) +
            geom_sf(aes(geometry = volcano), color = "black", fill = "white", shape = 24) +
            geom_sf(data = x$sim_large_buff[[1]], aes(color = est_sim), shape = 20, size = hf_pt_size) +
            geom_sf(
              data = comps_sim[!is.na(comps_sim$est_sim), ], aes(fill = est_sim), shape = 21,
              size = hf_pt_size, stroke = pt_stroke, color = "white"
            ) +
            geom_sf(aes(geometry = small_buffer1), fill = NA, linewidth = 0.5, color = "black") +
            geom_sf(aes(geometry = small_buffer2), fill = NA, linewidth = 0.5, color = "black") +
            geom_sf(aes(geometry = small_buffer3), fill = NA, linewidth = 0.5, color = "black") +
            ggtitle("Similarity Interpolation") +
            scale_color_viridis_c(
              option = "magma", name = bquote("Q" ~ (mWm^-2)), limits = c(0, 250),
              breaks = c(0, 125, 250), na.value = "transparent", guide = "none"
            ) +
            scale_fill_viridis_c(
              option = "magma", name = bquote("Q" ~ (mWm^-2)), limits = c(0, 250),
              breaks = c(0, 125, 250), na.value = "transparent", guide = "none"
            ) +
            coord_sf(expand = FALSE, lims_method = "geometry_bbox") +
            map_theme
          if (is.null(x$sim_projected1[[1]])) {
            p4 <- ggplot(data = data.frame(), aes())
          } else {
            p4 <-
              ggplot() +
              geom_point(
                data = x$sim_projected1[[1]], aes(projected_distances, obs), shape = 20,
                color = "grey20", size = 2
              ) +
              geom_point(
                data = x$sim_projected2[[1]], aes(projected_distances, obs), shape = 20,
                color = "grey20", size = 2
              ) +
              geom_point(
                data = x$sim_projected3[[1]], aes(projected_distances, obs), shape = 20,
                color = "grey20", size = 2
              ) +
              geom_path(data = x$sim_loess1[[1]], aes(projected_distances, obs), color = "black") +
              geom_path(data = x$sim_loess2[[1]], aes(projected_distances, obs), color = "black") +
              geom_path(data = x$sim_loess3[[1]], aes(projected_distances, obs), color = "black") +
              geom_rug(
                data = x$sim_loess3[[1]], aes(projected_distances, obs, color = obs),
                sides = "b", length = unit(0.06, "npc")
              )
          }
          p4 <- p4 +
            labs(x = "Normalized Distance", y = NULL) +
            scale_color_viridis_c(
              option = "magma", name = bquote("Q" ~ (mWm^-2)), limits = c(0, 250),
              breaks = c(0, 125, 250), na.value = "transparent",
              guide = guide_colorbar(
                title.vjust = 1, show.limits = TRUE,
                frame.colour = "black",
                ticks.colour = "black"
              )
            ) +
            scale_y_continuous(limits = c(0, 250), breaks = seq(0, 250, 50)) +
            scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) + profile_theme
          p <- (p3 / p4) + plot_layout(widths = 1, heights = c(1.5, 1)) &
            theme(plot.title = element_text(size = base_size * 1.5, margin = margin()))
          ggsave(file = fig_path2, plot = p, width = 6.5, height = 10, dpi = 300, bg = "white")
          p5 <-
            ggplot(x) +
            geom_sf(data = x$bathy[[1]], aes(color = elev), size = 0.5, shape = 15) +
            scale_color_etopo(guide = "none") +
            new_scale_color() +
            geom_sf(aes(geometry = large_buffer), fill = NA, linewidth = 0.5) +
            geom_sf(aes(geometry = ridge), color = "white") +
            geom_sf(aes(geometry = transform), color = "white") +
            geom_sf(aes(geometry = trench), color = "white", linewidth = 1.5) +
            geom_sf(aes(geometry = transect), color = "black", linewidth = 1.5) +
            geom_sf(aes(geometry = volcano), color = "black", fill = "white", shape = 24) +
            geom_sf(data = x$krg_large_buff[[1]], aes(color = est_krg), shape = 20, size = hf_pt_size) +
            geom_sf(
              data = comps_krg[!is.na(comps_krg$est_krg), ], aes(fill = est_krg), shape = 21,
              size = hf_pt_size, stroke = pt_stroke, color = "white"
            ) +
            geom_sf(aes(geometry = small_buffer1), fill = NA, linewidth = 0.5, color = "black") +
            geom_sf(aes(geometry = small_buffer2), fill = NA, linewidth = 0.5, color = "black") +
            geom_sf(aes(geometry = small_buffer3), fill = NA, linewidth = 0.5, color = "black") +
            ggtitle("Krige Interpolation") +
            scale_color_viridis_c(
              option = "magma", name = bquote("Q" ~ (mWm^-2)), limits = c(0, 250),
              breaks = c(0, 125, 250), na.value = "transparent", guide = "none"
            ) +
            scale_fill_viridis_c(
              option = "magma", name = bquote("Q" ~ (mWm^-2)), limits = c(0, 250),
              breaks = c(0, 125, 250), na.value = "transparent", guide = "none"
            ) +
            coord_sf(expand = FALSE, lims_method = "geometry_bbox") +
            map_theme
          if (is.null(x$krg_projected1[[1]])) {
            p6 <- ggplot(data = data.frame(), aes())
          } else {
            p6 <-
              ggplot() +
              geom_point(
                data = x$krg_projected1[[1]], aes(projected_distances, obs), shape = 20,
                color = "grey20", size = 2
              ) +
              geom_point(
                data = x$krg_projected2[[1]], aes(projected_distances, obs), shape = 20,
                color = "grey20", size = 2
              ) +
              geom_point(
                data = x$krg_projected3[[1]], aes(projected_distances, obs), shape = 20,
                color = "grey20", size = 2
              ) +
              geom_path(data = x$krg_loess1[[1]], aes(projected_distances, obs), color = "black") +
              geom_path(data = x$krg_loess2[[1]], aes(projected_distances, obs), color = "black") +
              geom_path(data = x$krg_loess3[[1]], aes(projected_distances, obs), color = "black") +
              geom_rug(
                data = x$krg_loess3[[1]], aes(projected_distances, obs, color = obs),
                sides = "b", length = unit(0.06, "npc")
              )
          }
          p6 <- p6 +
            labs(x = "Normalized Distance", y = NULL) +
            scale_color_viridis_c(
              option = "magma", name = bquote("Q" ~ (mWm^-2)), limits = c(0, 250),
              breaks = c(0, 125, 250), na.value = "transparent",
              guide = guide_colorbar(
                title.vjust = 1, show.limits = TRUE,
                frame.colour = "black",
                ticks.colour = "black"
              )
            ) +
            scale_y_continuous(limits = c(0, 250), breaks = seq(0, 250, 50)) +
            scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) + profile_theme
          p <- (p5 / p6) + plot_layout(widths = 1, heights = c(1.5, 1)) &
            theme(plot.title = element_text(size = base_size * 1.5, margin = margin()))
          ggsave(file = fig_path3, plot = p, width = 6.5, height = 10, dpi = 300, bg = "white")
          p7 <-
            (p1 + p5 + p3) / (p2 + p6 + p4) +
              plot_layout(widths = 1, heights = c(1.5, 1)) +
              plot_annotation(
                title = p_title,
                tag_levels = list(c("a)", "b)", "c)", "  ", "  ", "  ")),
                theme = theme(plot.title = element_text(
                  size = base_size * 1.5,
                  margin = margin()
                ))
              ) &
              theme(plot.tag = element_text(size = base_size * 1.5, margin = margin(0, 0, -10, 0)))
          ggsave(file = fig_path_comp, plot = p7, width = 19.5, height = 11, dpi = 300, bg = "white")
        })
      })
    },
    error = function(e) {
      cat("\n!! ERROR occurred in plot_transect:\n!!", conditionMessage(e))
    }
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot transect neighbors comp !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_transect_neighbors_comp <- function(trans_ids, base_size = 20) {
  if (length(trans_ids) < 3) {
    stop("\nNeed at least 3 transect ids!")
  }
  if (length(trans_ids) > 3) {
    trans_ids <- trans_ids[1:3]
  }
  fig_dir <- "figs/transect_neighbors/"
  trans_id_lab <- paste0(trans_ids[1], "-", trans_ids[2], "-", trans_ids[3])
  fig_paths <- list(
    paste0(fig_dir, trans_id_lab, "-transect-neighbors-ghf.png"),
    paste0(fig_dir, trans_id_lab, "-transect-neighbors-sim.png"),
    paste0(fig_dir, trans_id_lab, "-transect-neighbors-krg.png")
  )
  if (!dir.exists(fig_dir)) {
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (file.exists(fig_paths[[1]]) && file.exists(fig_paths[[2]]) && file.exists(fig_paths[[3]])) {
    cat("\n   ", fig_paths[[1]], " already exists ...", sep = "")
    cat("\n   ", fig_paths[[2]], " already exists ...", sep = "")
    cat("\n   ", fig_paths[[3]], " already exists ...", sep = "")
    return(invisible())
  }
  f <- function(trans_id) {
    tryCatch(
      {
        km <- get_optimal_krige_model(trans_id)
        fv <- km$fitted_vgrm
        np <- km$opt_krige_mod_summary$n_pairs
        x <- compile_trans_data(trans_id, fv = fv, np = np)
        comps_sim <- get_closest_interp_obs(trans_id)
        comps_krg <- get_closest_interp_obs(trans_id, fv = fv, np = np)
        hf_pt_size <- 3
        pt_stroke <- 0.8
        map_theme <- list(
          theme_bw(base_size = base_size),
          theme(
            plot.margin = margin(5, 5, 5, 5),
            plot.title = element_text(vjust = 0, hjust = 0.5, margin = margin(0, 5, 10, 5))
          ),
          annotation_scale(
            location = "bl", width_hint = 0.33, text_cex = 1.6, style = "ticks",
            line_width = 4, text_face = "bold"
          ),
          annotation_north_arrow(
            location = "bl", which_north = "true", pad_x = unit(0.0, "cm"),
            height = unit(2, "cm"), width = unit(2, "cm"),
            pad_y = unit(0.5, "cm"), style = north_arrow_fancy_orienteering
          )
        )
        profile_theme <-
          list(
            theme_bw(base_size = base_size),
            theme(
              panel.grid = element_blank(), panel.background = element_rect(fill = "grey90"),
              plot.margin = margin(5, 5, 5, 5),
              legend.justification = "right", legend.position = "inside",
              legend.position.inside = c(0.92, 0.85), legend.direction = "horizontal",
              legend.key.height = unit(0.5, "cm"), legend.key.width = unit(1, "cm"),
              legend.box.margin = margin(2, 2, 2, 2), legend.margin = margin(),
              legend.title = element_text(vjust = 0, size = base_size),
              legend.background = element_blank()
            )
          )
        suppressWarnings({
          suppressMessages({
            p1 <-
              ggplot(x) +
              geom_sf(data = x$bathy[[1]], aes(color = elev), size = 0.5, shape = 15) +
              scale_color_etopo(guide = "none") +
              new_scale_color() +
              geom_sf(aes(geometry = large_buffer), fill = NA, linewidth = 0.5) +
              geom_sf(aes(geometry = ridge), color = "white") +
              geom_sf(aes(geometry = transform), color = "white") +
              geom_sf(aes(geometry = trench), color = "white", linewidth = 1.5) +
              geom_sf(aes(geometry = transect), color = "black", linewidth = 1.5) +
              geom_sf(aes(geometry = volcano), color = "black", fill = "white", shape = 24) +
              geom_sf(data = x$ghf_large_buff[[1]], aes(color = obs), shape = 20, size = hf_pt_size) +
              geom_sf(aes(geometry = small_buffer1), fill = NA, linewidth = 0.5, color = "black") +
              geom_sf(aes(geometry = small_buffer2), fill = NA, linewidth = 0.5, color = "black") +
              geom_sf(aes(geometry = small_buffer3), fill = NA, linewidth = 0.5, color = "black") +
              ggtitle("Global HF Observations") +
              scale_color_viridis_c(
                option = "magma", name = bquote("Q" ~ (mWm^-2)), limits = c(0, 250),
                breaks = c(0, 125, 250), na.value = "transparent",
                guide = "none"
              ) +
              coord_sf(expand = FALSE, lims_method = "geometry_bbox") +
              map_theme
            if (is.null(x$ghf_projected1[[1]])) {
              p2 <- ggplot(data = data.frame(), aes())
            } else {
              p2 <-
                ggplot() +
                geom_point(
                  data = x$ghf_projected1[[1]], aes(projected_distances, obs), shape = 20,
                  color = "grey20", size = 2
                ) +
                geom_point(
                  data = x$ghf_projected2[[1]], aes(projected_distances, obs), shape = 20,
                  color = "grey20", size = 2
                ) +
                geom_point(
                  data = x$ghf_projected3[[1]], aes(projected_distances, obs), shape = 20,
                  color = "grey20", size = 2
                ) +
                geom_path(
                  data = x$ghf_loess1[[1]], aes(projected_distances, obs),
                  color = "black"
                ) +
                geom_path(
                  data = x$ghf_loess2[[1]], aes(projected_distances, obs),
                  color = "black"
                ) +
                geom_path(
                  data = x$ghf_loess3[[1]], aes(projected_distances, obs),
                  color = "black"
                ) +
                geom_rug(
                  data = x$ghf_loess3[[1]], aes(projected_distances, obs, color = obs),
                  sides = "b", length = unit(0.06, "npc")
                )
            }
            p2 <- p2 +
              labs(x = "Normalized Distance", y = NULL) +
              scale_color_viridis_c(
                option = "magma", name = bquote("Q" ~ (mWm^-2)), limits = c(0, 250),
                breaks = c(0, 125, 250), na.value = "transparent",
                guide = guide_colorbar(
                  title.vjust = 1, show.limits = TRUE,
                  frame.colour = "black",
                  ticks.colour = "black"
                )
              ) +
              scale_y_continuous(limits = c(0, 250), breaks = seq(0, 250, 50)) +
              scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) + profile_theme
            p3 <-
              ggplot(x) +
              geom_sf(data = x$bathy[[1]], aes(color = elev), size = 0.5, shape = 15) +
              scale_color_etopo(guide = "none") +
              new_scale_color() +
              geom_sf(aes(geometry = large_buffer), fill = NA, linewidth = 0.5) +
              geom_sf(aes(geometry = ridge), color = "white") +
              geom_sf(aes(geometry = transform), color = "white") +
              geom_sf(aes(geometry = trench), color = "white", linewidth = 1.5) +
              geom_sf(aes(geometry = transect), color = "black", linewidth = 1.5) +
              geom_sf(aes(geometry = volcano), color = "black", fill = "white", shape = 24) +
              geom_sf(
                data = x$sim_large_buff[[1]], aes(color = est_sim), shape = 20,
                size = hf_pt_size
              ) +
              geom_sf(
                data = comps_sim[!is.na(comps_sim$est_sim), ], aes(fill = est_sim),
                shape = 21, size = hf_pt_size, stroke = pt_stroke, color = "white"
              ) +
              geom_sf(aes(geometry = small_buffer1), fill = NA, linewidth = 0.5, color = "black") +
              geom_sf(aes(geometry = small_buffer2), fill = NA, linewidth = 0.5, color = "black") +
              geom_sf(aes(geometry = small_buffer3), fill = NA, linewidth = 0.5, color = "black") +
              ggtitle("Similarity Interpolation") +
              scale_color_viridis_c(
                option = "magma", name = bquote("Q" ~ (mWm^-2)), limits = c(0, 250),
                breaks = c(0, 125, 250), na.value = "transparent",
                guide = "none"
              ) +
              scale_fill_viridis_c(
                option = "magma", name = bquote("Q" ~ (mWm^-2)), limits = c(0, 250),
                breaks = c(0, 125, 250), na.value = "transparent",
                guide = "none"
              ) +
              coord_sf(expand = FALSE, lims_method = "geometry_bbox") +
              map_theme
            if (is.null(x$sim_projected1[[1]])) {
              p4 <- ggplot(data = data.frame(), aes())
            } else {
              p4 <-
                ggplot() +
                geom_point(
                  data = x$sim_projected1[[1]], aes(projected_distances, obs), shape = 20,
                  color = "grey20", size = 2
                ) +
                geom_point(
                  data = x$sim_projected2[[1]], aes(projected_distances, obs), shape = 20,
                  color = "grey20", size = 2
                ) +
                geom_point(
                  data = x$sim_projected3[[1]], aes(projected_distances, obs), shape = 20,
                  color = "grey20", size = 2
                ) +
                geom_path(data = x$sim_loess1[[1]], aes(projected_distances, obs), color = "black") +
                geom_path(data = x$sim_loess2[[1]], aes(projected_distances, obs), color = "black") +
                geom_path(data = x$sim_loess3[[1]], aes(projected_distances, obs), color = "black") +
                geom_rug(
                  data = x$sim_loess3[[1]], aes(projected_distances, obs, color = obs),
                  sides = "b", length = unit(0.06, "npc")
                )
            }
            p4 <- p4 +
              labs(x = "Normalized Distance", y = NULL) +
              scale_color_viridis_c(
                option = "magma", name = bquote("Q" ~ (mWm^-2)), limits = c(0, 250),
                breaks = c(0, 125, 250), na.value = "transparent",
                guide = guide_colorbar(
                  title.vjust = 1, show.limits = TRUE,
                  frame.colour = "black",
                  ticks.colour = "black"
                )
              ) +
              scale_y_continuous(limits = c(0, 250), breaks = seq(0, 250, 50)) +
              scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) + profile_theme
            p5 <-
              ggplot(x) +
              geom_sf(data = x$bathy[[1]], aes(color = elev), size = 0.5, shape = 15) +
              scale_color_etopo(guide = "none") +
              new_scale_color() +
              geom_sf(aes(geometry = large_buffer), fill = NA, linewidth = 0.5) +
              geom_sf(aes(geometry = ridge), color = "white") +
              geom_sf(aes(geometry = transform), color = "white") +
              geom_sf(aes(geometry = trench), color = "white", linewidth = 1.5) +
              geom_sf(aes(geometry = transect), color = "black", linewidth = 1.5) +
              geom_sf(aes(geometry = volcano), color = "black", fill = "white", shape = 24) +
              geom_sf(
                data = x$krg_large_buff[[1]], aes(color = est_krg), shape = 20,
                size = hf_pt_size
              ) +
              geom_sf(
                data = comps_krg[!is.na(comps_krg$est_krg), ], aes(fill = est_krg),
                shape = 21, size = hf_pt_size, stroke = pt_stroke, color = "white"
              ) +
              geom_sf(aes(geometry = small_buffer1), fill = NA, linewidth = 0.5, color = "black") +
              geom_sf(aes(geometry = small_buffer2), fill = NA, linewidth = 0.5, color = "black") +
              geom_sf(aes(geometry = small_buffer3), fill = NA, linewidth = 0.5, color = "black") +
              ggtitle("Krige Interpolation") +
              scale_color_viridis_c(
                option = "magma", name = bquote("Q" ~ (mWm^-2)), limits = c(0, 250),
                breaks = c(0, 125, 250), na.value = "transparent",
                guide = "none"
              ) +
              scale_fill_viridis_c(
                option = "magma", name = bquote("Q" ~ (mWm^-2)), limits = c(0, 250),
                breaks = c(0, 125, 250), na.value = "transparent",
                guide = "none"
              ) +
              coord_sf(expand = FALSE, lims_method = "geometry_bbox") +
              map_theme
            if (is.null(x$krg_projected1[[1]])) {
              p6 <- ggplot(data = data.frame(), aes())
            } else {
              p6 <-
                ggplot() +
                geom_point(
                  data = x$krg_projected1[[1]], aes(projected_distances, obs), shape = 20,
                  color = "grey20", size = 2
                ) +
                geom_point(
                  data = x$krg_projected2[[1]], aes(projected_distances, obs), shape = 20,
                  color = "grey20", size = 2
                ) +
                geom_point(
                  data = x$krg_projected3[[1]], aes(projected_distances, obs), shape = 20,
                  color = "grey20", size = 2
                ) +
                geom_path(data = x$krg_loess1[[1]], aes(projected_distances, obs), color = "black") +
                geom_path(data = x$krg_loess2[[1]], aes(projected_distances, obs), color = "black") +
                geom_path(data = x$krg_loess3[[1]], aes(projected_distances, obs), color = "black") +
                geom_rug(
                  data = x$krg_loess3[[1]], aes(projected_distances, obs, color = obs),
                  sides = "b", length = unit(0.06, "npc")
                )
            }
            p6 <- p6 +
              labs(x = "Normalized Distance", y = NULL) +
              scale_color_viridis_c(
                option = "magma", name = bquote("Q" ~ (mWm^-2)), limits = c(0, 250),
                breaks = c(0, 125, 250), na.value = "transparent",
                guide = guide_colorbar(
                  title.vjust = 1, show.limits = TRUE,
                  frame.colour = "black",
                  ticks.colour = "black"
                )
              ) +
              scale_y_continuous(limits = c(0, 250), breaks = seq(0, 250, 50)) +
              scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) + profile_theme
            return(list(p1, p2, p3, p4, p5, p6))
          })
        })
      },
      error = function(e) {
        cat("\n!! ERROR occurred in plot_transect:\n!!", conditionMessage(e))
      }
    )
  }
  plan(multisession, workers = availableCores() - 2)
  plots <- future_map(trans_ids, f, .options = furrr_options(seed = seed)) |> reduce(c)
  suppressWarnings({
    suppressMessages({
      p_title <- paste0("Submap Transect: ", trans_id_lab, " ")
      p1 <-
        (plots[[1]] + plots[[7]] + plots[[13]]) / (plots[[2]] + plots[[8]] + plots[[14]]) +
          plot_layout(widths = 1, heights = c(1.5, 1)) +
          plot_annotation(
            title = p_title, tag_levels = list(c("a)", "b)", "c)", "  ", "  ", "  ")),
            theme = theme(plot.title = element_text(
              size = base_size * 1.5,
              margin = margin()
            ))
          ) &
          theme(plot.tag = element_text(size = base_size * 1.5, margin = margin(0, 0, -10, 0)))
      cat("\n   Plotting: ", fig_paths[[1]], " ...", sep = "")
      ggsave(file = fig_paths[[1]], plot = p1, width = 19.5, height = 11, dpi = 300, bg = "white")
      p2 <-
        (plots[[3]] + plots[[9]] + plots[[15]]) / (plots[[4]] + plots[[10]] + plots[[16]]) +
          plot_layout(widths = 1, heights = c(1.5, 1)) +
          plot_annotation(
            title = p_title, tag_levels = list(c("a)", "b)", "c)", "  ", "  ", "  ")),
            theme = theme(plot.title = element_text(
              size = base_size * 1.5,
              margin = margin()
            ))
          ) &
          theme(plot.tag = element_text(size = base_size * 1.5, margin = margin(0, 0, -10, 0)))
      cat("\n   Plotting: ", fig_paths[[2]], " ...", sep = "")
      ggsave(file = fig_paths[[2]], plot = p2, width = 19.5, height = 11, dpi = 300, bg = "white")
      p3 <-
        (plots[[5]] + plots[[11]] + plots[[17]]) / (plots[[6]] + plots[[12]] + plots[[18]]) +
          plot_layout(widths = 1, heights = c(1.5, 1)) +
          plot_annotation(
            title = p_title, tag_levels = list(c("a)", "b)", "c)", "  ", "  ", "  ")),
            theme = theme(plot.title = element_text(
              size = base_size * 1.5,
              margin = margin()
            ))
          ) &
          theme(plot.tag = element_text(size = base_size * 1.5, margin = margin(0, 0, -10, 0)))
      cat("\n   Plotting: ", fig_paths[[3]], " ...", sep = "")
      ggsave(file = fig_paths[[3]], plot = p3, width = 19.5, height = 11, dpi = 300, bg = "white")
    })
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot transect strip comp !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_transect_strip_comp <- function(trans_ids = NULL, fname = NULL, buff = 1, base_size = 14) {
  if (is.null(trans_ids)) {
    stop("\nMissing submap transect ids!")
  }
  if (is.null(fname)) {
    fname <- "test"
  }
  if (buff == 1) {
    select_cols <- c("short_name", "ghf_loess1", "krg_loess1", "sim_loess1")
  } else if (buff == 2) {
    select_cols <- c("short_name", "ghf_loess2", "krg_loess2", "sim_loess2")
  } else if (buff == 3) {
    select_cols <- c("short_name", "ghf_loess3", "krg_loess3", "sim_loess3")
  } else {
    select_cols <- c("short_name", "ghf_loess1", "krg_loess1", "sim_loess1")
  }
  fig_dir <- "figs/summary/"
  if (!dir.exists(fig_dir)) {
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  }
  fig_path <- paste0(fig_dir, "strip-", fname, ".png")
  if (file.exists(fig_path)) {
    cat("\n   ", fig_path, " already exists ...", sep = "")
    return(invisible())
  }
  tryCatch(
    {
      cat("\n   Plotting: ", fig_path, " ...", sep = "")
      f <- function(x) {
        km <- get_optimal_krige_model(x)
        fv <- km$fitted_vgrm
        np <- km$opt_krige_mod_summary$n_pairs
        compile_trans_data(x, fv = fv, np = np)
      }
      plan(multisession, workers = availableCores() - 2)
      df <-
        future_map(trans_ids, f, .options = furrr_options(seed = seed)) |>
        reduce(rbind) |>
        st_set_geometry(NULL) |>
        select(all_of(select_cols)) |>
        pivot_longer(-c(short_name)) |>
        unnest(value)
      p <-
        ggplot(df) +
        geom_vline(aes(xintercept = projected_distances, color = obs)) +
        scale_x_continuous(breaks = c(0, 1)) +
        labs(x = "Normalized Distance", y = NULL) +
        facet_grid(vars(short_name), vars(name)) +
        scale_color_viridis_c(
          option = "magma", name = bquote("Q" ~ (mWm^-2)), limits = c(0, 250),
          breaks = c(0, 125, 250), na.value = "transparent", guide = "none"
        ) +
        theme_bw(base_size = base_size) +
        theme(
          panel.grid = element_blank(), panel.background = element_rect(fill = "grey90"),
          plot.margin = margin(5, 5, 5, 5),
          legend.justification = "right", legend.position = "inside",
          legend.position.inside = c(0.92, 0.85), legend.direction = "horizontal",
          legend.key.height = unit(0.5, "cm"), legend.key.width = unit(1, "cm"),
          legend.box.margin = margin(2, 2, 2, 2), legend.margin = margin(),
          legend.title = element_text(vjust = 0, size = base_size),
          legend.background = element_blank()
        )
      ggsave(
        file = fig_path, plot = p, width = 6.5, height = length(trans_ids) * 0.75, dpi = 300,
        bg = "white"
      )
    },
    error = function(e) {
      cat("\n!! ERROR occurred in plot_transect_strip_comp:\n!!", conditionMessage(e))
    }
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot cross correlation !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_cross_correlation <- function(trans_id1, trans_id2, ccf_thresh = 0.8, lag_thresh = 0.2,
                                   max_lags = 1e3, base_size = 14) {
  fig_dir <- "figs/transect_xcorr/"
  fig_path <- paste0(fig_dir, trans_id1, "-", trans_id2, "-transect-xcorr.png")
  if (!dir.exists(fig_dir)) {
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (file.exists(fig_path)) {
    cat("\n   ", fig_path, " already exists ...", sep = "")
    return(invisible())
  }
  tryCatch(
    {
      cat("\n   Plotting: ", fig_path, " ...", sep = "")
      kmx <- get_optimal_krige_model(trans_id1)
      kmy <- get_optimal_krige_model(trans_id2)
      fvx <- kmx$fitted_vgrm
      fvy <- kmy$fitted_vgrm
      npx <- kmx$opt_krige_mod_summary$n_pairs
      npy <- kmy$opt_krige_mod_summary$n_pairs
      x <- compile_trans_data(trans_id1, fv = fvx, np = npx)
      y <- compile_trans_data(trans_id2, fv = fvy, np = npy)
      color_mapping <- setNames(c("darkorange", "navy"), c(trans_id1, trans_id2))
      p_theme <- list(
        theme_bw(base_size = base_size),
        theme(
          panel.grid = element_blank(), panel.background = element_rect(fill = "grey90"),
          plot.margin = margin(5, 5, 5, 5), plot.title = element_text(hjust = 0.5),
          legend.justification = "right", legend.position = "inside",
          legend.position.inside = c(0.92, 0.85), legend.direction = "horizontal",
          legend.key.height = unit(0.5, "cm"), legend.key.width = unit(1, "cm"),
          legend.box.margin = margin(2, 2, 2, 2), legend.margin = margin(),
          legend.title = element_text(vjust = 0, size = base_size),
          legend.background = element_blank()
        )
      )
      smooth1_ghf <- x$ghf_loess3[[1]]
      points1_ghf <- x$ghf_projected3[[1]]
      smooth2_ghf <- y$ghf_loess3[[1]]
      points2_ghf <- y$ghf_projected3[[1]]
      ccf_result_ghf <-
        ccf(smooth1_ghf$obs, smooth2_ghf$obs, na.action = na.pass, lag.max = max_lags, plot = FALSE)
      ccf_result_ghf$lag <- ccf_result_ghf$lag * (1 / (nrow(smooth1_ghf) - 1))
      max_ccf_ghf <- max(ccf_result_ghf$acf, na.rm = TRUE)
      max_lag_ghf <- max(ccf_result_ghf$lag[which(ccf_result_ghf$acf == max_ccf_ghf)])
      smooth1_sim <- x$sim_loess3[[1]]
      points1_sim <- x$sim_projected3[[1]]
      smooth2_sim <- y$sim_loess3[[1]]
      points2_sim <- y$sim_projected3[[1]]
      ccf_result_sim <-
        ccf(smooth1_sim$obs, smooth2_sim$obs, na.action = na.pass, lag.max = max_lags, plot = FALSE)
      ccf_result_sim$lag <- ccf_result_sim$lag * (1 / (nrow(smooth1_sim) - 1))
      max_ccf_sim <- max(ccf_result_sim$acf, na.rm = TRUE)
      max_lag_sim <- max(ccf_result_sim$lag[which(ccf_result_sim$acf == max_ccf_sim)])
      smooth1_krg <- x$krg_loess3[[1]]
      points1_krg <- x$krg_projected3[[1]]
      smooth2_krg <- y$krg_loess3[[1]]
      points2_krg <- y$krg_projected3[[1]]
      ccf_result_krg <-
        ccf(smooth1_krg$obs, smooth2_krg$obs, na.action = na.pass, lag.max = max_lags, plot = FALSE)
      ccf_result_krg$lag <- ccf_result_krg$lag * (1 / (nrow(smooth1_krg) - 1))
      max_ccf_krg <- max(ccf_result_krg$acf, na.rm = TRUE)
      max_lag_krg <- max(ccf_result_krg$lag[which(ccf_result_krg$acf == max_ccf_krg)])
      #    if ((max_lag_ghf > -lag_thresh & max_lag_ghf < lag_thresh & max_ccf_ghf > ccf_thresh) |
      #        (max_lag_sim > -lag_thresh & max_lag_sim < lag_thresh & max_ccf_sim > ccf_thresh) |
      #        (max_lag_krg > -lag_thresh & max_lag_krg < lag_thresh & max_ccf_krg > ccf_thresh)) {
      p1 <-
        ggplot() +
        geom_point(
          data = points1_ghf, aes(projected_distances, obs, color = trans_id1),
          na.rm = TRUE, shape = 20, size = 0.7
        ) +
        geom_point(
          data = points2_ghf, aes(projected_distances, obs, color = trans_id2),
          na.rm = TRUE, shape = 20, size = 0.7
        ) +
        geom_path(
          data = smooth1_ghf, aes(projected_distances, obs, color = trans_id1),
          linewidth = 1, na.rm = TRUE
        ) +
        geom_path(
          data = smooth2_ghf, aes(projected_distances, obs, color = trans_id2),
          linewidth = 1, na.rm = TRUE
        ) +
        geom_path(
          data = smooth2_ghf,
          aes(projected_distances + max_lag_ghf, obs, color = trans_id2),
          linewidth = 1, linetype = 2, na.rm = TRUE
        ) +
        labs(x = "Normalized Distance", y = bquote("Q" ~ (mWm^-2))) +
        ggtitle("Global Heat Flow") +
        scale_color_manual(name = NULL, values = color_mapping) +
        p_theme
      p2 <-
        ggplot(tibble(acf = ccf_result_ghf$acf, lag = ccf_result_ghf$lag)) +
        geom_hline(yintercept = ccf_thresh, linetype = 2) +
        geom_hline(yintercept = -ccf_thresh, linetype = 2) +
        geom_hline(yintercept = 0) +
        geom_path(aes(lag, acf), linewidth = 1, na.rm = TRUE) +
        geom_point(data = tibble(acf = max_ccf_ghf, lag = max_lag_ghf), aes(lag, acf)) +
        labs(x = "Normalized Lag Distance", y = "Correlation") +
        scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
        scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
        p_theme
      p3 <-
        ggplot() +
        geom_point(
          data = points1_krg, aes(projected_distances, obs, color = trans_id1),
          na.rm = TRUE, shape = 20, size = 0.7
        ) +
        geom_point(
          data = points2_krg, aes(projected_distances, obs, color = trans_id2),
          na.rm = TRUE, shape = 20, size = 0.7
        ) +
        geom_path(
          data = smooth1_krg, aes(projected_distances, obs, color = trans_id1),
          linewidth = 1, na.rm = TRUE
        ) +
        geom_path(
          data = smooth2_krg, aes(projected_distances, obs, color = trans_id2),
          linewidth = 1, na.rm = TRUE
        ) +
        geom_path(
          data = smooth2_krg,
          aes(projected_distances + max_lag_krg, obs, color = trans_id2),
          linewidth = 1, linetype = 2, na.rm = TRUE
        ) +
        labs(x = "Normalized Distance", y = bquote("Q" ~ (mWm^-2))) +
        ggtitle("Krige") +
        scale_color_manual(name = NULL, values = color_mapping) +
        p_theme
      p4 <-
        ggplot(tibble(acf = ccf_result_krg$acf, lag = ccf_result_krg$lag)) +
        geom_hline(yintercept = ccf_thresh, linetype = 2) +
        geom_hline(yintercept = -ccf_thresh, linetype = 2) +
        geom_hline(yintercept = 0) +
        geom_path(aes(lag, acf), linewidth = 1, na.rm = TRUE) +
        geom_point(data = tibble(acf = max_ccf_krg, lag = max_lag_krg), aes(lag, acf)) +
        labs(x = "Normalized Lag Distance", y = "Correlation") +
        scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
        scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
        p_theme
      p5 <-
        ggplot() +
        geom_point(
          data = points1_sim, aes(projected_distances, obs, color = trans_id1),
          na.rm = TRUE, shape = 20, size = 0.7
        ) +
        geom_point(
          data = points2_sim, aes(projected_distances, obs, color = trans_id2),
          na.rm = TRUE, shape = 20, size = 0.7
        ) +
        geom_path(
          data = smooth1_sim, aes(projected_distances, obs, color = trans_id1),
          linewidth = 1, na.rm = TRUE
        ) +
        geom_path(
          data = smooth2_sim, aes(projected_distances, obs, color = trans_id2),
          linewidth = 1, na.rm = TRUE
        ) +
        geom_path(
          data = smooth2_sim,
          aes(projected_distances + max_lag_sim, obs, color = trans_id2),
          linewidth = 1, linetype = 2, na.rm = TRUE
        ) +
        labs(x = "Normalized Distance", y = bquote("Q" ~ (mWm^-2))) +
        ggtitle("Similarity") +
        scale_color_manual(name = NULL, values = color_mapping) +
        p_theme
      p6 <-
        ggplot(tibble(acf = ccf_result_sim$acf, lag = ccf_result_sim$lag)) +
        geom_hline(yintercept = ccf_thresh, linetype = 2) +
        geom_hline(yintercept = -ccf_thresh, linetype = 2) +
        geom_hline(yintercept = 0) +
        geom_path(aes(lag, acf), linewidth = 1, na.rm = TRUE) +
        geom_point(data = tibble(acf = max_ccf_sim, lag = max_lag_sim), aes(lag, acf)) +
        labs(x = "Normalized Lag Distance", y = "Correlation") +
        scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
        scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
        p_theme
      p_title <- paste0("Cross-correlation: ", trans_id1, " ", trans_id2)
      p_rax <- theme(
        axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
      p <-
        (p1 + (p3 + p_rax) + (p5 + p_rax)) / (p2 + (p4 + p_rax) + (p6 + p_rax)) +
          plot_annotation(
            title = p_title,
            tag_levels = list(c("a)", "b)", "c)", "  ", "  ", "  ")),
            theme = theme(plot.title = element_text(
              size = base_size * 1.5,
              margin = margin()
            ))
          ) &
          theme(plot.tag = element_text(size = base_size * 1.5, margin = margin(0, 0, -10, 0)))
      ggsave(file = fig_path, plot = p, width = 13, height = 5, dpi = 300, bg = "white")
      #    } else {
      #      cat('\n   No significant correlation found between lag threshold!')
      #    }
    },
    error = function(e) {
      cat("\n!! ERROR occurred in plot_cross_correlation:\n!!", conditionMessage(e))
    }
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot optimal krige model !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_optimal_krige_model <- function(trans_id = NULL, base_size = 22) {
  if (is.null(trans_id)) {
    stop("\nMissing submap transect ids!")
  }
  fig_dir <- "figs/nlopt/"
  fig_path <- paste0(fig_dir, trans_id, "-opt-krige.png")
  if (!dir.exists(fig_dir)) {
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (file.exists(fig_path)) {
    cat("\n   ", fig_path, " already exists ...", sep = "")
    return(invisible())
  }
  tryCatch(
    {
      cat("\n   Plotting:", fig_path)
      opt_krige_mod <- get_optimal_krige_model(trans_id)
      opt_kmod <- opt_krige_mod$opt_krige_mod_summary
      nlopt_itr <- opt_krige_mod$nlopt_itr
      ev <- opt_krige_mod$experimental_vgrm
      fv <- opt_krige_mod$fitted_vgrm
      fv_line <- variogramLine(fv, maxdist = max(ev$dist))
      p_theme <-
        list(
          theme_bw(base_size = base_size),
          theme(
            panel.grid = element_blank(), panel.background = element_rect(fill = "grey90"),
            plot.margin = margin(5, 5, 5, 5),
            plot.title = element_text(vjust = 0, hjust = 0.5, margin = margin(10, 10, 10, 10)),
            legend.position = "none"
          )
        )
      suppressWarnings({
        suppressMessages({
          p0 <-
            ggplot(filter(nlopt_itr, v_mod == opt_kmod$v_mod)) +
            geom_path(aes(itr, max_dist / 1e3), color = "forestgreen", linewidth = 1) +
            geom_point(data = opt_kmod, aes(itr, max_dist / 1e3), color = "forestgreen", size = 3) +
            geom_label_repel(
              data = opt_kmod, size = base_size * 0.28,
              aes(itr, max_dist / 1e3, label = round(max_dist / 1e3))
            ) +
            annotate("text",
              x = Inf, y = Inf, label = "Variogram Max Distance", hjust = 1.05, vjust = 1.5,
              size = base_size * 0.35
            ) +
            labs(x = "NLopt Iteration", y = "Distance (km)", color = NULL) +
            p_theme
          p1 <-
            ggplot(filter(nlopt_itr, v_mod == opt_kmod$v_mod)) +
            geom_path(aes(itr, n_pairs), color = "firebrick", linewidth = 1) +
            geom_point(data = opt_kmod, aes(itr, n_pairs), color = "firebrick", size = 3) +
            geom_label_repel(
              data = opt_kmod, aes(itr, n_pairs, label = round(n_pairs)),
              size = base_size * 0.28
            ) +
            annotate("text",
              x = Inf, y = Inf, label = "Variogram Pairs", hjust = 1.05, vjust = 1.5,
              size = base_size * 0.35
            ) +
            labs(x = "NLopt Iteration", y = "Pairs", color = NULL) +
            p_theme
          p2 <-
            ggplot(filter(nlopt_itr, v_mod == opt_kmod$v_mod)) +
            geom_path(aes(itr, cutoff), color = "orchid4", linewidth = 1) +
            geom_point(data = opt_kmod, aes(itr, cutoff), color = "orchid4", size = 3) +
            geom_label_repel(
              data = opt_kmod, aes(itr, cutoff, label = round(cutoff)),
              size = base_size * 0.28
            ) +
            annotate("text",
              x = Inf, y = Inf, label = "Variogram Cutoff", hjust = 1.05, vjust = 1.5,
              size = base_size * 0.35
            ) +
            labs(x = "NLopt Iteration", y = "Cutoff", color = NULL) +
            p_theme
          p3 <-
            ggplot(filter(nlopt_itr, v_mod == opt_kmod$v_mod)) +
            geom_path(aes(itr, n_lags), color = "saddlebrown", linewidth = 1) +
            geom_point(data = opt_kmod, aes(itr, n_lags), color = "saddlebrown", size = 3) +
            geom_label_repel(
              data = opt_kmod, aes(itr, n_lags, label = round(n_lags)),
              size = base_size * 0.28
            ) +
            annotate("text",
              x = Inf, y = Inf, label = "Variogram Lags", hjust = 1.05, vjust = 1.5,
              size = base_size * 0.35
            ) +
            labs(x = "NLopt Iteration", y = "Lags", color = NULL) +
            p_theme
          p4 <-
            ggplot(filter(nlopt_itr, v_mod == opt_kmod$v_mod)) +
            geom_path(aes(itr, cost), color = "black", linewidth = 1) +
            geom_point(data = opt_kmod, aes(itr, cost), color = "black", shape = 20, size = 5) +
            geom_label_repel(
              data = opt_kmod, aes(itr, cost, label = paste0(round(cost, 3))),
              size = base_size * 0.28
            ) +
            annotate("text",
              x = Inf, y = Inf, label = "Cost Function", hjust = 1.05, vjust = 1.5,
              size = base_size * 0.35
            ) +
            labs(x = "NLopt Iteration", y = "Cost", color = NULL) +
            p_theme
          p5 <-
            ggplot(ev) +
            geom_point(aes(x = dist / 1e3, y = sqrt(gamma)), shape = 20) +
            geom_line(data = fv_line, aes(x = dist / 1e3, y = sqrt(gamma)), linewidth = 1) +
            annotate("text",
              x = Inf, y = -Inf, label = paste0("Optimal model (", opt_kmod$v_mod, ")"),
              hjust = 1.05, vjust = -0.5, size = base_size * 0.35
            ) +
            labs(x = "Lag Distance (km)", y = bquote("Variance" ~ (mWm^-2))) +
            p_theme
          p6 <- (p0 + p1) / (p2 + p3) / (p4 + p5) +
            plot_annotation(
              title = paste0("Kriging Optimization: ", unique(nlopt_itr$short_name)),
              theme = theme(plot.title = element_text(size = base_size * 1.2, hjust = 0.5))
            ) +
            plot_layout(widths = 1, heights = 1) &
            theme(plot.tag = element_text(size = base_size * 1.5))
          ggsave(file = fig_path, plot = p6, width = 13, height = 10, dpi = 300, bg = "white")
        })
      })
    },
    error = function(e) {
      cat("\n!! ERROR occurred in plot_optimal_variogram:\n!!", conditionMessage(e))
    }
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot interp accuracy summary !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_interp_accuracy_summary <- function(submap_zone = NULL, base_size = 22) {
  load_nlopt_interpolation_data("assets/nlopt_data/interpolation-summary.RData")
  if (is.null(submap_zone)) {
    stop("\nMissing submap zone!")
  }
  fig_dir <- "figs/summary/"
  if (!dir.exists(fig_dir)) {
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  }
  f <- function(x, y) {
    if (file.exists(y)) {
      cat("\n   ", y, " already exists ...", sep = "")
      return(invisible())
    }
    tryCatch(
      {
        cat("\n   Plotting:", y)
        p_theme <- list(
          theme_bw(base_size = 14),
          theme(
            panel.grid = element_blank(), panel.background = element_rect(fill = "grey90"),
            axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "bottom",
            plot.title = element_text(hjust = 0.5)
          )
        )
        outliers_diff <-
          interp_diff_summary$short_name[which(interp_diff_summary$mean + 2 * interp_diff_summary$sigma > 300)]
        outliers_acc <-
          interp_accuracy_summary$short_name[which(interp_accuracy_summary$rmse_obs_sim > 300 | interp_accuracy_summary$rmse_obs_krg > 300)]
        outliers <- unique(c(outliers_diff, outliers_acc))
        df_acc_rmse <-
          interp_accuracy_summary |>
          filter(!is.na(rmse_obs_krg)) |>
          pivot_longer(-c(short_name)) |>
          mutate(
            method = ifelse(str_detect(name, "sim"), "sim", "krg"),
            name = str_split(name, "_", simplify = TRUE)[, 1],
            zone = str_sub(short_name, 1, 3)
          ) |>
          rename(metric = name) |>
          filter(zone == submap_zone) |>
          filter(short_name %in% x & !(short_name %in% outliers)) |>
          mutate(method = ifelse(method == "sim", "Similarity", "Krige")) |>
          filter(metric == "rmse") |>
          rename(rmse = value) |>
          select(-metric) |>
          left_join(select(nlopt_summary, short_name, n_obs), by = "short_name") |>
          left_join(select(interp_accuracy_summary, short_name, n_krg), by = "short_name") |>
          rename(n_control = n_krg) |>
          left_join(select(interp_diff_summary, short_name, n_grid), by = "short_name") |>
          left_join(select(interp_diff_summary, short_name, n), by = "short_name") |>
          mutate(n = ifelse(method == "Similarity", n_grid, n)) |>
          rename(n_itp = n) |>
          mutate(coverage = n_itp / n_grid * 100)
        suppressWarnings({
          suppressMessages({
            p0 <-
              nlopt_summary |>
              filter(short_name %in% x & !(short_name %in% outliers)) |>
              ggplot() +
              geom_col(aes(short_name, n_obs, fill = "Observation"), color = "black") +
              geom_col(
                data = filter(df_acc_rmse, method == "Krige"),
                aes(short_name, n_control, fill = "Control point"), color = "black"
              ) +
              scale_fill_manual(values = c(
                "Observation" = "forestgreen",
                "Control point" = "darkred"
              )) +
              labs(x = NULL, y = "N", color = NULL, fill = "Observation type") +
              ggtitle("Heat flow observations") +
              p_theme +
              theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
            p1 <-
              interp_diff_summary |>
              filter(short_name %in% x & !(short_name %in% outliers)) |>
              ggplot() +
              geom_hline(yintercept = 0) +
              geom_crossbar(aes(short_name, mean,
                ymin = mean - 2 * sigma,
                ymax = mean + 2 * sigma
              )) +
              labs(x = NULL, y = bquote("Difference" ~ (mWm^-2)), color = NULL) +
              ggtitle("Point-by-Point Interpolation Differences") +
              p_theme +
              theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
            p2 <-
              interp_diff_summary |>
              filter(short_name %in% x & !(short_name %in% outliers)) |>
              ggplot() +
              geom_col(aes(short_name, 100, fill = "Similarity"), color = "black") +
              geom_col(aes(short_name, n / n_grid * 100, fill = "Krige"), color = "black") +
              scale_fill_manual(
                values = c("Krige" = "darkorange", "Similarity" = "navy"),
                guide = "none"
              ) +
              labs(x = NULL, y = "Coverage (%)", fill = "Method") +
              ggtitle("Spatial Coverage") +
              p_theme +
              theme(
                axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                legend.position = "none"
              )
            p3 <-
              df_acc_rmse |>
              ggplot() +
              geom_hline(yintercept = 0) +
              geom_col(aes(short_name, rmse, fill = method), color = "black", position = "identity") +
              scale_fill_manual(values = c("Krige" = "darkorange", "Similarity" = "navy")) +
              labs(x = NULL, y = bquote("RMSE" ~ (mWm^-2)), fill = "Method") +
              ggtitle("Interpolation Accuracies") +
              p_theme
            p4 <- p0 / p1 / p2 / p3 + plot_layout(guides = "collect") &
              theme(legend.position = "bottom")
            ggsave(file = y, plot = p4, width = 13, height = 10, dpi = 300, bg = "white")
          })
        })
      },
      error = function(e) {
        cat("\n!! ERROR occurred in plot_interp_accuracy_summary:\n!!", conditionMessage(e))
      }
    )
  }
  x <- nlopt_summary$short_name[str_detect(nlopt_summary$short_name, submap_zone)]
  if (submap_zone %in% c("SAM", "SEA")) {
    midpoint <- length(x) %/% 2
    walk2(list(x[1:midpoint], x[(midpoint + 1):length(x)]), c(1, 2), ~ {
      f(.x, paste0(fig_dir, submap_zone, "-interpolation-accuracy-summary", .y, ".png"))
    })
  } else {
    f(x, paste0(fig_dir, submap_zone, "-interpolation-accuracy-summary.png"))
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot nlopt summary !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_nlopt_summary <- function(base_size = 14) {
  load_nlopt_interpolation_data("assets/nlopt_data/interpolation-summary.RData")
  fig_dir <- "figs/summary/"
  fig_path <- paste0(fig_dir, "nlopt-summary.png")
  if (!dir.exists(fig_dir)) {
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (file.exists(fig_path)) {
    cat("\n   ", fig_path, " already exists ...", sep = "")
    return(invisible())
  }
  tryCatch(
    {
      cat("\n   Plotting:", fig_path)
      outliers_diff <-
        interp_diff_summary$short_name[which(interp_diff_summary$mean + 2 * interp_diff_summary$sigma > 300)]
      outliers_acc <-
        interp_accuracy_summary$short_name[which(interp_accuracy_summary$rmse_obs_sim > 300 | interp_accuracy_summary$rmse_obs_krg > 300)]
      outliers <- unique(c(outliers_diff, outliers_acc))
      suppressWarnings({
        suppressMessages({
          p1 <-
            nlopt_summary |>
            filter(!(short_name %in% outliers)) |>
            mutate(zone = str_sub(short_name, 1, 3)) |>
            select(-c(vgrm_wt, vgrm_cost, vgrm_rmse, cv_wt, cv_cost, cv_rmse)) |>
            pivot_longer(-c(short_name, zone, v_mod, cost)) |>
            ggplot() +
            geom_point(aes(value, cost, fill = zone), shape = 21, color = "black") +
            labs(x = NULL, y = "Cost", color = "Zone") +
            scale_fill_manual(values = c(
              "NPA" = "darkorange", "SAM" = "navy", "SEA" = "darkred",
              "SWP" = "forestgreen"
            )) +
            facet_wrap(~name, scales = "free_x") +
            theme_bw(base_size = 14) +
            theme(
              panel.grid = element_blank(), panel.background = element_rect(fill = "grey90"),
              strip.background = element_blank()
            )
          ggsave(file = fig_path, plot = p1, width = 6.5, height = 4, dpi = 300, bg = "white")
        })
      })
    },
    error = function(e) {
      cat("\n!! ERROR occurred in plot_nlopt_summary:\n!!", conditionMessage(e))
    }
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot control point summary !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_control_point_summary <- function(base_size = 14) {
  load_nlopt_interpolation_data("assets/nlopt_data/interpolation-summary.RData")
  fig_dir <- "figs/summary/"
  fig_path <- paste0(fig_dir, "control-point-summary.png")
  if (!dir.exists(fig_dir)) {
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (file.exists(fig_path)) {
    cat("\n   ", fig_path, " already exists ...", sep = "")
    return(invisible())
  }
  p_theme <- list(
    theme_bw(base_size = 14),
    theme(
      panel.grid = element_blank(), panel.background = element_rect(fill = "grey90"),
      plot.title = element_text(hjust = 0.5)
    )
  )
  df <-
    interp_accuracy_summary |>
    filter(!is.na(rmse_obs_krg)) |>
    pivot_longer(-c(short_name)) |>
    mutate(
      method = ifelse(str_detect(name, "sim"), "sim", "krg"),
      name = str_split(name, "_", simplify = TRUE)[, 1],
      zone = str_sub(short_name, 1, 3)
    ) |>
    rename(metric = name) |>
    mutate(method = ifelse(method == "sim", "Similarity", "Krige")) |>
    filter(metric == "rmse") |>
    rename(rmse = value) |>
    select(-metric) |>
    left_join(select(nlopt_summary, short_name, n_obs), by = "short_name") |>
    left_join(select(interp_accuracy_summary, short_name, n_krg), by = "short_name") |>
    rename(n_control = n_krg) |>
    left_join(select(interp_diff_summary, short_name, n_grid), by = "short_name") |>
    left_join(select(interp_diff_summary, short_name, n), by = "short_name") |>
    mutate(n = ifelse(method == "Similarity", n_grid, n)) |>
    rename(n_itp = n) |>
    mutate(coverage = n_control / n_grid * 100)
  tryCatch(
    {
      cat("\n   Plotting:", fig_path)
      suppressWarnings({
        suppressMessages({
          p1 <-
            ggplot(df) +
            geom_point(aes(coverage, rmse, fill = zone), shape = 21, color = "black") +
            scale_fill_manual(values = c(
              "NPA" = "darkorange", "SAM" = "navy", "SEA" = "darkred",
              "SWP" = "forestgreen"
            )) +
            facet_wrap(~method, scales = "free_x") +
            labs(x = "Control coverage (%)", y = "RMSE", fill = "Zone") +
            p_theme +
            theme(axis.title.x = element_blank())
          p2 <-
            ggplot(filter(df, coverage <= 20)) +
            geom_point(aes(coverage, rmse, fill = zone), shape = 21, color = "black") +
            scale_fill_manual(values = c(
              "NPA" = "darkorange", "SAM" = "navy", "SEA" = "darkred",
              "SWP" = "forestgreen"
            ), guide = "none") +
            facet_wrap(~method, scales = "free_x") +
            labs(x = "Control coverage (%)", y = "RMSE", fill = "Zone") +
            p_theme +
            theme(axis.title.x = element_blank())
          p3 <-
            ggplot(filter(df, coverage <= 5)) +
            geom_point(aes(coverage, rmse, fill = zone), shape = 21, color = "black") +
            scale_fill_manual(values = c(
              "NPA" = "darkorange", "SAM" = "navy", "SEA" = "darkred",
              "SWP" = "forestgreen"
            ), guide = "none") +
            facet_wrap(~method, scales = "free_x") +
            labs(x = "Control coverage (%)", y = "RMSE", fill = "Zone") +
            p_theme
          p4 <- p1 / p2 / p3 + plot_layout(guides = "collect")
          ggsave(file = fig_path, plot = p4, width = 6.5, height = 8, dpi = 300, bg = "white")
        })
      })
    },
    error = function(e) {
      cat("\n!! ERROR occurred in plot_control_point_summary:\n!!", conditionMessage(e))
    }
  )
}
