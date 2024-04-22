#!/usr/bin/env Rscript

# Load packages and data
cat(rep('~', 45), '\n', sep='')
source('R/functions.R')
load('assets/map_data/map-data.RData')

# Set variogram models and submap zones
v_mods <- c('Sph', 'Exp', 'Lin', 'Bes')
zones <- c('NPA', 'SAM', 'SEA', 'SWP')

walk(zones, ~{
  # Optimize kriging models after Li et al. (2018)
  ids <- shp_submap$short_name[shp_submap$zone == .x]
  nlopt_transects(ids, v_mods)
  # Visualize best variogram solutions
  walk(ids, plot_optimal_variogram)
})

cat('\nnloptr.R complete!\n\n')