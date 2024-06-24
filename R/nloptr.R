#!/usr/bin/env Rscript

# Load packages and data
cat(rep('~', 45), '\n', sep='')
source('R/functions.R')
load('assets/map_data/map-data.RData')
plan(multicore, workers=availableCores() - 2)

# Nlopt setup
seed <- 42
alg <- 'NLOPT_LN_SBPLX'
ids <- shp_submap$short_name
v_mods <- c('Sph', 'Exp', 'Lin')
fpath <- 'assets/nlopt_data/interpolation-summary.RData'

# Optimize krige models
if (!file.exists(fpath)) {
  nlopt_transects(ids, v_mods, alg)
  nlopt_summary <- summarize_optimal_krige_models(ids)
  interp_diff_summary <- summarize_interpolation_differences(ids)
  interp_accuracy_summary <- summarize_interpolation_accuracy(ids)
  save(list=c('nlopt_summary', 'interp_diff_summary', 'interp_accuracy_summary'), file=fpath)
  plot_nlopt_summary()
  plot_control_point_summary()
  future_walk(ids, plot_transect_buff_comp, .options=furrr_options(seed=seed), .progress=T)
  future_walk(ids, plot_optimal_krige_model, .options=furrr_options(seed=seed), .progress=T)
  future_walk(c('NPA', 'SAM', 'SEA', 'SWP'), plot_interp_accuracy_summary,
              .options=furrr_options(seed=seed), .progress=T)
} else {
  cat('\nOptimal krige model summary found at:\n  ', fpath, '\n', sep='')
  load('assets/nlopt_data/interpolation-summary.RData')
  plot_nlopt_summary()
  plot_control_point_summary()
  future_walk(ids, plot_transect_buff_comp, .options=furrr_options(seed=seed), .progress=T)
  future_walk(ids, plot_optimal_krige_model, .options=furrr_options(seed=seed), .progress=T)
  future_walk(c('NPA', 'SAM', 'SEA', 'SWP'), plot_interp_accuracy_summary,
              .options=furrr_options(seed=seed), .progress=T)
}

cat('\nnloptr.R complete!\n')