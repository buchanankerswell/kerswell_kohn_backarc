#!/usr/bin/env Rscript

# Load packages and functions
cat(rep('~', 45), '\n', sep='')
source('R/functions.R')
load('assets/map_data/map-data.RData')

# Plot base maps
plot_tglobe_base()
walk(shp_submap$short_name, plot_transect)

cat('\nbase-plots.R complete!\n\n')