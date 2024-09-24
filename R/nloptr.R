#!/usr/bin/env Rscript

# Load functions and libraries
source('R/functions.R')
load('assets/hf_data/preprocessed-hf-data.RData')

# Testing arguments
args <- commandArgs(trailingOnly=T)
parse_krige_args(args)

# Set parallel computing plan
plan(multicore, workers=n_cores)

# Print settings
cat('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',
    '\n~~~~~                         ~~~~~~',
    '\n~~~~~     nloptr settings     ~~~~~~',
    '\n~~~~~                         ~~~~~~',
    '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',
    '\nAlgorithm       : ', alg,
    '\nMax evaluations : ', max_eval,
    '\nk-fold CV       : ', n_fold,
    '\nInterp weight   : ', iwt,
    '\nVgrm weight     : ', vwt,
    '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')

# Optimize over v_mods for each segment
v_mods <- c('Sph', 'Exp', 'Lin', 'Bes')
seg_names <- extract_RData_object('assets/map_data/preprocessed-map-data.RData', 'seg_names')
optimized_interpolations <-
  tibble(expand.grid(segment=seg_names, v_mod=v_mods, stringsAsFactors=F)) %>%
  arrange(segment, v_mod) %>%
  mutate(opt=future_map2(segment, v_mod,
                         ~optimize_krige(..., shp_hf_crop, alg, max_eval, n_fold, iwt, vwt),
                         .options=furrr_options(seed=42), .progress=T))

# Save nloptr results
save(optimized_interpolations, file="assets/opt_data/nloptr.RData")

cat('\nnloptr.R complete!\n\n')