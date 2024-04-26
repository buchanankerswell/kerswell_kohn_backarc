#!/usr/bin/env Rscript

# Load packages and data
cat(rep('~', 45), '\n', sep='')
source('R/functions.R')
load('assets/map_data/map-data.RData')

# Plot global base map
plot_ghfdb_base()

# Set variogram models and submap zones
v_mods <- c('Sph', 'Exp', 'Lin', 'Bes')
ids <- shp_submap$short_name

fpath <- 'assets/nlopt_data/interpolation-summary.RData'
if (!file.exists(fpath)) {
  cat('\nOptimizing krige models after Li et al. (2018) ...\n')
  nlopt_transects(ids, v_mods)
  walk(ids, plot_transect)
  walk(ids, plot_optimal_krige_model)
  cat('\nSummarizing optimal krige models ...\n')
  nlopt_summary <- summarize_optimal_krige_models(ids)
  cat('\nSummarizing interpolation differences ...\n')
  interp_diff_summary <- summarize_interpolation_differences(ids)
  cat('\nSummarizing interpolation accuracies ...\n')
  interp_accuracy_summary <- summarize_interpolation_accuracy(ids)
  save(list=c('nlopt_summary', 'interp_diff_summary', 'interp_accuracy_summary'), file=fpath)
} else {
  cat('\nOptimal krige model summary found at:\n', fpath, '\n')
}

cat('\nnloptr.R complete!\n\n')