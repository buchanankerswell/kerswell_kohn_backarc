#!/usr/bin/env Rscript

source('R/functions.R')
load_map_data('assets/map_data/map-data.RData')
load_nlopt_data('assets/nlopt_data/interpolation-summary.RData')

ids <- shp_submap$short_name