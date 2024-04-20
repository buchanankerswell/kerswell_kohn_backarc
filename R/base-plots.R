#!/usr/bin/env Rscript

# Load packages and functions
cat(rep('~', 45), '\n', sep='')
source('R/functions.R')
load('assets/map_data/map-data.RData')

# Plot global base map
plot_tglobe_base()

# Plot submap transects
short_names <- shp_submap$short_name
walk(short_names, plot_transect_tglobe)
walk(short_names, plot_transect_sim)

cat('\nbase-plots.R complete!\n\n')