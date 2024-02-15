#!/usr/bin/env Rscript

# Load functions and libraries
source('R/functions.R')
load('assets/hf_data/preprocessed-hf-data.RData')

# Testing arguments
args <- commandArgs(trailingOnly = TRUE)
parse_krige_args(args)

# Set parallel computing plan
plan(multicore, workers = n_cores)

# Print settings
cat(
  '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',
  '\n~~~~~                         ~~~~~~',
  '\n~~~~~     nloptr settings     ~~~~~~',
  '\n~~~~~                         ~~~~~~',
  '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',
  '\nAlgorithm       : ', alg,
  '\nMax evaluations : ', max_eval,
  '\nk-fold CV       : ', n_fold,
  '\nInterp weight   : ', iwt,
  '\nVgrm weight     : ', vwt,
  '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n'
)

# Minimization function
f <- function(segment, v_mod) {
  # Define cost function
  opt_fun <- function(x) {
    hf_obs <- shp_hf_crop[[segment]]
    cost_function(hf_obs, x[1], x[2], x[3], v_mod, x[4], n_fold, iwt, vwt, segment, T)
  }
  # Initial search values, bounds, and nlopt options (cutoff_prop, n_lags, lag_start, n_max)
  x0 <- c(3, 20, 1, 8)
  lb <- c(1, 15, 1, 5)
  ub <- c(10, 50, 10, 50)
  opts <- list(print_level = 0, maxeval = max_eval, algorithm = alg)
  # Run nloptr to minimize cost function
  nloptr(x0, opt_fun, lb = lb, ub = ub, opts = opts)
}

# Create array of variogram models to optimize for each segment in parallel
v_mods <- c('Sph', 'Exp', 'Lin', 'Bes')
seg_names <- extract_RData_object('assets/map_data/preprocessed-map-data.RData', 'seg_names')
optimized_interpolations <-
  tibble(expand.grid(segment = seg_names, v_mod = v_mods, stringsAsFactors = F)) %>%
  arrange(segment, v_mod) %>%
  mutate(opt = future_map2(segment, v_mod, possibly(f, otherwise = NULL),
                           .progress = T, .options = furrr_options(seed = 42)))

# Save optimized kriging interpolations
save(optimized_interpolations, file = "assets/hf_data/optimize-krige.RData")

cat('\noptimize-krige.R complete!\n\n')