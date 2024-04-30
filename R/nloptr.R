#!/usr/bin/env Rscript

# Load packages and data
cat(rep('~', 45), '\n', sep='')
source('R/functions.R')
load('assets/map_data/map-data.RData')

# Nlopt setup
ids <- shp_submap$short_name
alg <- 'NLOPT_LN_NELDERMEAD'
v_mods <- c('Sph', 'Exp', 'Lin')
fpath <- 'assets/nlopt_data/interpolation-summary.RData'

if (!file.exists(fpath)) {
  # Optimize krige models and visualize results
  nlopt_transects(ids, v_mods, alg)

  # Summarize nlopt results
  nlopt_summary <- summarize_optimal_krige_models(ids)
  interp_diff_summary <- summarize_interpolation_differences(ids)
  interp_accuracy_summary <- summarize_interpolation_accuracy(ids)

  # Save nlopt results
  save(list=c('nlopt_summary', 'interp_diff_summary', 'interp_accuracy_summary'), file=fpath)

  # Visualize nlopt results
  walk(ids, plot_transect)
  walk(ids, plot_optimal_krige_model)

} else {
  cat('\nOptimal krige model summary found at:\n', fpath, '\n', sep='')
}

cat('\nnloptr.R complete!\n')