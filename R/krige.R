#!/usr/bin/env Rscript

# Load functions and libraries
source('R/functions.R')
load('assets/map_data/preprocessed-map-data.RData')
load('assets/hf_data/preprocessed-hf-data.RData')
load('assets/opt_data/nloptr.RData')

# Read nloptr trace
cat('\nCompiling nloptr solutions ...')
nloptr_trace <- read_nloptr_trace(paste0('log/nlopt-out-', format(Sys.time(), '%d-%m-%Y')))

# Decode nloptr results and construct variograms
opt_vgrms <- optimized_interpolations %>% drop_na() %>%
  pmap(~decode_opt(shp_hf_crop[[..1]], ..2, ..3))

# Compile solutions
opt_solutions <- optimized_interpolations %>% drop_na() %>%
  mutate(opt_exp_vgrm=map(opt_vgrms, ~.x[[1]]), opt_fit_vgrm=map(opt_vgrms, ~.x[[2]]),
         cost=map_dbl(optimized_interpolations$opt, ~.x$objective))

# Get solutions with minimum cost
best_krige_models <- opt_solutions %>% group_by(segment) %>% slice(which.min(cost))

# Summarise optimized variogram models
nloptr_trace_min <- nloptr_trace %>% group_by(segment, v_mod) %>% slice_min(cost) %>% ungroup()
vgrm_summary <-
  tibble(segment=opt_solutions$segment, v_mod=opt_solutions$v_mod,
         cutoff=map_dbl(opt_solutions$opt, ~.x$solution[1]),
         n_lags=map_dbl(opt_solutions$opt, ~.x$solution[2]),
         lag_start=map_dbl(opt_solutions$opt, ~.x$solution[3]),
         n_max=map_dbl(opt_solutions$opt, ~.x$solution[4]),
         sill=map_dbl(opt_solutions$opt_fit_vgrm, ~.x$psill),
         range=map_dbl(opt_solutions$opt_fit_vgrm, ~.x$range)) %>%
  left_join(nloptr_trace_min, by=c('segment', 'v_mod')) %>%
  left_join(rmse_similarity_interpolation, by=c('segment'))

# Krige segments with optimal variogram models
cat('\nKriging segments with optimal variogram models ...')
krige_params <- list(opt_solutions$segment, opt_solutions$opt_fit_vgrm,
                     map(opt_solutions$opt, ~.x$solution[4]))
shp_interp_krige <- krige_params %>%
  pmap(~Krige(shp_hf_crop[[..1]], ..2, shp_grid_crop[[..1]], ..3, ..1)) %>%
  map(~filter(.x, est_krige >= 0)) %>% set_names(opt_solutions$segment)

# Calculating similarity–krige difference
shp_interp_diff <- shp_interp_krige %>% map(~interp_diff(.x, shp_interp_sim)) %>%
  set_names(opt_solutions$segment)

# Compile solutions
opt_solutions <- opt_solutions %>% mutate(shp_interp_diff, .before=cost)

# Summarise similarity–krige differences
cat('\nComputing interpolation residuals and summarizing ...')
interp_diff_summary <-
  opt_solutions %>%
  pmap_df(~st_set_geometry(..6, NULL) %>%
          mutate(segment=..1, v_mod=..2, cost=..7, .before=est_sim)) %>%
  group_by(segment, v_mod, cost) %>%
  summarise(n=n(), min=round(min(est_diff, na.rm=T)), max=round(max(est_diff, na.rm=T)),
            median=round(median(est_diff, na.rm=T)), IQR=round(IQR(est_diff, na.rm=T)),
            mean=round(mean(est_diff, na.rm=T)), sigma=round(sd(est_diff, na.rm=T)),
            .groups='drop')

# Calculate direct RMSE
vgrm_summary$dir_rmse <-
  map2_dbl(opt_solutions$segment, opt_solutions$shp_interp_diff,
           ~interpolation_rmse(.x, .y, shp_buffer, shp_grid_crop, shp_hf_crop, 'krg'))

# Save results
save(list=c('opt_solutions', 'best_krige_models', 'nloptr_trace', 'vgrm_summary',
            'interp_diff_summary'), file='assets/opt_data/krige.RData')

cat('\nkrige.R complete!\n\n')