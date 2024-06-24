#!/usr/bin/env Rscript

ids <- shp_submap$short_name

source('R/functions.R')
load('assets/map_data/map-data.RData')
load('assets/nlopt_data/interpolation-summary.RData')

source('R/functions.R')
file.remove('figs/transect_buff/NPA02-transect-buff-comp.png')
plot_transect_buff_comp(ids[2])

plan(multicore, workers=availableCores() - 2)
future_walk(seq(length(ids) - 2), ~{
  ids_triplet <- ids[.x:(.x + 2)]
  plot_transect_neighbors_comp(ids_triplet)
}, .options=furrr_options(seed=seed), .progress=T)