#!/usr/bin/env Rscript

main <- function() {
  # Load packages and data
  source('R/functions.R')
  load_map_data('assets/map_data/map-data.RData')

  # Setup
  seed <- 42
  set.seed(seed)
  plan(multisession, workers=availableCores() - 2)
  fpath <- 'assets/nlopt_data/interpolation-summary.RData'

  # Submap transects/zones
  ids <- shp_submap$short_name
  zones <- c('NPA', 'SAM', 'SEA', 'SWP')

  # Optimize krige models
  if (!file.exists(fpath)) {
    nlopt_transects(ids)
    nlopt_summary <- summarize_optimal_krige_models(ids)
    interp_diff_summary <- summarize_interp_differences(ids)
    interp_accuracy_summary <- summarize_interp_accuracy(ids)
    save(nlopt_summary, interp_diff_summary, interp_accuracy_summary, file=fpath)
  } else {
    cat('\nOptimal krige model summary found at:', fpath)
  }

  # Visualize interpolations
  if (file.exists(fpath)) {
    cat('\nVisualizing interpolations ...')
    ids <- c(ids[1:3], ids[12:14], ids[43:45], ids[120:122], ids[179:181])
    id_comb <- combn(ids, 2, simplify=F)
    plot_nlopt_summary()
    plot_control_point_summary()
    future_walk(ids, plot_transect_buff_comp, .options=furrr_options(seed=seed))
    future_walk(ids, plot_optimal_krige_model, .options=furrr_options(seed=seed))
    future_walk(zones, plot_interp_accuracy_summary, .options=furrr_options(seed=seed))
    plot_transect_neighbors_comp(ids[1:3])
    plot_transect_neighbors_comp(ids[4:6])
    plot_transect_neighbors_comp(ids[7:9])
    plot_transect_neighbors_comp(ids[10:12])
    plot_transect_neighbors_comp(ids[13:15])
    plot_transect_neighbors_comp(ids[c(2, 5, 8)])
    plot_transect_strip_comp(ids, 'comp')
    future_walk(id_comb, ~plot_cross_correlation(.x[1], .x[2]),
                .options=furrr_options(seed=seed))
  } else {
    cat('\nNo summary data found at:', fpath)
  }
  cat('\n', rep('=', 45), sep='')
}

main()