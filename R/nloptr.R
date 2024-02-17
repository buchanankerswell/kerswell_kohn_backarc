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

# Minimization function
f <- function(segment, v_mod) {
  opt_fun <- function(x) {
    hf_obs <- shp_hf_crop[[segment]]
    cost_function(hf_obs, x[1], x[2], x[3], x[4], v_mod, n_fold, iwt, vwt, segment)
  }
  x0 <- c(3, 50, 3, 10) # Initial search values (cutoff, n_lags, lag_start, n_max)
  lb <- c(1, 30, 1, 2) # Lower bounds
  ub <- c(12, 100, 10, 50) # Upper bounds
  opts <- list(print_level=0, maxeval=max_eval, algorithm=alg, ftol_rel=1e-5)
  nloptr(x0, opt_fun, lb=lb, ub=ub, opts=opts)
}

# Optimize over v_mods for each segment
v_mods <- c('Sph', 'Exp', 'Lin', 'Bes')
seg_names <- extract_RData_object('assets/map_data/preprocessed-map-data.RData', 'seg_names')
optimized_interpolations <-
  tibble(expand.grid(segment=seg_names, v_mod=v_mods, stringsAsFactors=F)) %>%
  arrange(segment, v_mod) %>%
  mutate(opt=future_map2(segment, v_mod, f, .options=furrr_options(seed=42), .progress=T))

# Save nloptr results
save(optimized_interpolations, file="assets/opt_data/nloptr.RData")

cat('\nnloptr.R complete!\n\n')