#!/usr/bin/env Rscript

# Load functions and libraries
source('R/functions.R')
load('assets/hf_data/optimized-krige-interpolations.R')

# Check for NULL results from errors
if(any(map_lgl(optimized_interpolations$opt, is.null))) {
  cat('\n', rep('~!', 30), sep = '')
  cat(
    '\n\nOptimization failed for',
    sum(map_lgl(optimized_interpolations$opt, is.null)),
    'runs:\n'
  )
  print(optimized_interpolations[map_lgl(optimized_interpolations$opt, is.null),])
  cat('\nContinuing computation for good runs ...')
  cat('\n', rep('~!', 30), sep = '')
}

# Drop missing values
optimized_interpolations <- drop_na(optimized_interpolations)

# Parse nloptr history
read_file <- function(fpath, ...) {
  con <- file(fpath)
  on.exit(close(con))
  readLines(con, ...)
}
raw_trace <- read_file(paste0('log-', Sys.Date()))
segs_trace <-
  raw_trace[grepl('^Segment: ', raw_trace)] %>%
  map_chr(~gsub('Segment: ', '', .x))
vgrm_model_trace <-
  raw_trace[grepl('^Variogram model:', raw_trace)] %>%
  map_chr(~gsub('Variogram model: ', '', .x))
vgrm_weight_trace <-
  raw_trace[grepl('^Variogram weight', raw_trace)] %>%
  map_dbl(~as.numeric(gsub('Variogram weight: ', '', .x)))
vgrm_rmse_trace <-
  raw_trace[grepl('^Variogram rmse', raw_trace)] %>%
  map_dbl(~as.numeric(gsub('Variogram rmse: ', '', .x)))
vgrm_cost_trace <-
  raw_trace[grepl('^Variogram cost', raw_trace)] %>%
  map_dbl(~as.numeric(gsub('Variogram cost: ', '', .x)))
cv_weight_trace <-
  raw_trace[grepl('^Interpolation weight', raw_trace)] %>%
  map_dbl(~as.numeric(gsub('Interpolation weight: ', '', .x)))
cv_rmse_trace <-
  raw_trace[grepl('^Interpolation rmse', raw_trace)] %>%
  map_dbl(~as.numeric(gsub('Interpolation rmse: ', '', .x)))
cv_cost_trace <-
  raw_trace[grepl('^Interpolation cost', raw_trace)] %>%
  map_dbl(~as.numeric(gsub('Interpolation cost: ', '', .x)))
cost_trace <-
  raw_trace[grepl('^Cost', raw_trace)] %>%
  map_dbl(~as.numeric(gsub('Cost: ', '', .x)))
opt_trace <-
  tibble(
    segment = segs_trace,
    v_mod = vgrm_model_trace,
    vgrm_wt = vgrm_weight_trace,
    vgrm_rmse = vgrm_rmse_trace,
    vgrm_cost = vgrm_cost_trace,
    cv_wt = cv_weight_trace,
    cv_rmse = cv_rmse_trace,
    cv_cost = cv_cost_trace,
    cost = cost_trace
  ) %>%
  group_by(segment, v_mod) %>%
  mutate(itr = row_number(), .before = segment) %>%
  ungroup()

cat('\nnloptr trace ...')
cat('\n', rep('-', 40), '\n', sep = '')
print(opt_trace)
cat(rep('-', 40), '\n', sep = '')

# Decode optimization output and construct variograms
cat('\n\n', rep('~', 80), sep = '')
cat('\nDecoding optimized solutions ...')
opt_vgrms <-
  optimized_interpolations %>%
  pmap(
    possibly(
      ~decode_opt(
        model_vgrm = ..2,
        shp_hf = shp_hf_crop[[..1]],
        opt = ..3
      ),
      otherwise = NULL
    )
  )
# Add variograms and cost to tibble for easy manipulation and query
solns <-
  optimized_interpolations %>%
  mutate(
    opt_exp_vgrm = map(opt_vgrms, ~.x[[1]]),
    opt_fit_vgrm = map(opt_vgrms, ~.x[[2]]),
    cost = map_dbl(optimized_interpolations$opt, ~.x$objective)
  )

# Check for NULL results from errors
if(any(map_lgl(solns$opt_exp_vgrm, is.null))) {
  cat('\n', rep('~!', 30), sep = '')
  cat(
    '\n\nVariogram fitting failed during opt decoding for',
    sum(map_lgl(solns$opt_exp_vgrm, is.null)),
    'runs:\n'
  )
  print(solns[map_lgl(solns$opt_exp_vgrm, is.null),])
  cat('\nContinuing computation for good runs ...')
  cat('\n', rep('~!', 30), sep = '')
}

# Get solutions with minimum cost
opt_solns <-
  solns %>%
  drop_na() %>%
  group_by(segment) %>%
  slice(which.min(cost))

# Print lowest-cost solutions
cat('\nBest solutions for each segment:\n')
cat('\n', rep('-', 40), sep = '')
pwalk(opt_solns, ~{
  cat(
    '\nSegment:         ', ..1,
    '\nVariogram cutoff:', ..3$solution[1], 
    '\nNumber of lags:  ', ..3$solution[2],
    '\nLag shift:       ', ..3$solution[3],
    '\nMax points:      ', ..3$solution[4],
    '\nVariogram model: ', ..2,
    '\nCost:            ', ..3$objective,
    '\n'
  )
})
cat('\n', rep('-', 40), '\n', sep = '')

# Drop missing values
solns <- drop_na(solns)

# Summarise variograms
cat('\nVariogram summary:\n')
opt_trace_end <-
  opt_trace %>%
  group_by(segment, v_mod) %>%
  slice_tail() %>%
  ungroup()
vgrm_summary <-
  tibble(
    segment = solns$segment,
    v_mod = solns$v_mod,
    cutoff_prop = map_dbl(solns$opt, ~.x$solution[1]),
    n_lags = map_dbl(solns$opt, ~.x$solution[2]),
    lag_start = map_dbl(solns$opt, ~.x$solution[3]),
    n_max = map_dbl(solns$opt, ~.x$solution[4]),
    sill = map_dbl(solns$opt_fit_vgrm, ~.x$psill),
    range = map_dbl(solns$opt_fit_vgrm, ~.x$range),
  ) %>%
  left_join(opt_trace_end, by = c('segment', 'v_mod')) %>%
  left_join(rmse_luca, by = c('segment'))
print(vgrm_summary)

# Kriging segments
cat('\n', rep('~', 80), sep = '')
cat('\nKriging segments ...\n')
shp_interp_krige <-
  list(
    solns$segment,
    solns$opt_fit_vgrm,
    map(solns$opt, ~.x$solution[4])
  ) %>%
  pmap(
    possibly(
      ~Krige(
       shp_hf = shp_hf_crop[[..1]],
       fitted_vgrm = ..2,
       shp_interp_grid = shp_grid_crop[[..1]],
       n_max = ..3,
       seg_name = ..1
      ),
      otherwise = NULL
    )
  ) %>%
  set_names(solns$segment)

# Filter negative heat flow estimates
shp_interp_krige <- map(shp_interp_krige, ~filter(.x, est_krige >= 0))

# Check for NULL results from errors
if(any(map_lgl(shp_interp_krige, is.null))) {
  cat('\n', rep('~!', 30), sep = '')
  cat(
    '\n\nVariogram fitting failed during opt decoding for',
    sum(map_lgl(shp_interp_krige, is.null)),
    'runs:\n'
  )
  indx <- map_lgl(shp_interp_krige, is.null)
  cat('\nSegments:\n')
  writeLines(solns$segment[indx])
  cat('\nVariogram models:\n')
  print(solns$opt_fit_vgrm[indx])
  cat('\nContinuing computation for good runs ...')
  cat('\n', rep('~!', 30), sep = '')
}

# Remove missing values
solns <-
  solns %>%
  mutate(shp_interp_krige, .before = cost) %>%
  drop_na()

# Calculating interpolation difference
cat('\n', rep('~', 80), sep = '')
cat('\nComputing differences ...\n')
shp_interp_diff <-
  solns %>%
  pmap(
    possibly(
      ~{
        interp_diff(
          shp_interp_krige = ..6,
          shp_interp_sim = shp_interp_luca
        )
      },
      otherwise = NULL
    )
  ) %>% set_names(solns$segment)

# Check for NULL results from errors
if(any(map_lgl(shp_interp_diff, is.null))) {
  cat('\n', rep('~!', 30), sep = '')
  cat(
    '\n\nVariogram fitting failed during opt decoding for',
    sum(map_lgl(shp_interp_diff, is.null)),
    'runs:\n'
  )
  indx <- map_lgl(shp_interp_diff, is.null)
  cat('\nSegments:\n')
  writeLines(solns$segment[indx])
  cat('\nVariogram models:\n')
  print(solns$opt_fit_vgrm[indx])
  cat('\nContinuing computation for good runs ...')
  cat('\n', rep('~!', 30), sep = '')
}

# Remove missing values
solns <-
  solns %>%
  mutate(shp_interp_diff, .before = cost) %>%
  drop_na()

# Summarise interpolation differences
interp_diff_summary <-
  solns %>%
  pmap_df(
    ~st_set_geometry(..7, NULL) %>%
    mutate(segment = ..1, v_mod = ..2, cost = ..8, .before = est_sim)
  ) %>%
  group_by(segment, v_mod, cost) %>%
  summarise(
    n = n(),
    min = round(min(est_diff, na.rm = T)),
    max = round(max(est_diff, na.rm = T)),
    median = round(median(est_diff, na.rm = T)),
    IQR = round(IQR(est_diff, na.rm = T)),
    mean = round(mean(est_diff, na.rm = T)),
    sigma = round(sd(est_diff, na.rm = T)),
    .groups = 'drop'
  )
cat('\nHeat flow difference summary:\n')
print(interp_diff_summary)

# Drop krige column since it is a duplicate of diff column
solns <- select(solns, -shp_interp_krige)

# Calculate direct RMSE for Kriging predictions (to compare with cross-validation RMSE)
vgrm_summary$dir_rmse <-
  map2_dbl(solns$segment, solns$shp_interp_diff, ~itp_rmse(.x, .y, 'krg'))

cat('\n', rep('~', 80), sep = '')
cat('\nSaving to: data/opt.RData', sep = '')
save(
  list = c(
    'solns',
    'opt_solns',
    'opt_trace',
    'vgrm_summary',
    'interp_diff_summary'
  ),
  file = 'data/opt.RData'
)

cat('\nkrige.R complete!\n\n')