#!/usr/bin/env Rscript

# Load packages and data
cat(rep('~', 45), '\n', sep='')
source('R/functions.R')
load('assets/map_data/map-data.RData')

# Set variogram models and submap zones
v_mods <- c('Sph', 'Exp', 'Lin', 'Bes')
zones <- c('NPA', 'SAM', 'SEA', 'SWP')

# Optimize kriging models after Li et al. (2018)
nlopt_transects(shp_submap$short_name, v_mods)

# Visualize best variogram solutions
walk(shp_submap$short_name, plot_optimal_variogram)

cat('\nnloptr.R complete!\n\n')