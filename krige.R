#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  max.eval <- 100
  alg <- 'NLOPT_GN_DIRECT_L'            # Global search
} else if (length(args) != 0) {
  max.eval <- suppressWarnings(as.numeric(args[1]))
  if(is.na(max.eval)){
    cat('\nPassed non-numeric argument for max iterations!')
    cat('\nDefaulting to 100')
    max.eval <- 100
  }
  if(args[2] == '1') {
    alg <- 'NLOPT_GN_DIRECT_L'          # Global search
  } else if(args[2] == '2') {
    alg <- 'NLOPT_LN_SBPLX'             # Local without gradients
  } else if(args[2] == '3') {
    algorithm <- 'NLOPT_LN_NELDERMEAD'  # Local without gradients
  } else if(args[2] == '4') {
    algorithm <- 'NLOPT_LN_BOBYQA'      # Local without gradients
  } else if(args[2] == '5') {
    algorithm <- 'NLOPT_LN_COBYLA'      # Local without gradients
  } else {
    cat('\nAlgorithm passed is not recognized!')
    cat('\nDefaulting to NLOPT_GN_DIRECT_L')
    alg <- 'NLOPT_GN_DIRECT_L'
  }
}

# Counter
opt.files <- list.files('data/', pattern = 'opt')
cat('\n', rep('~', 60), '\n', sep='')
cat('Found', length(opt.files), 'opt files already in data/\n')

cntr <- length(opt.files) + 1
cat('Saving results to: data/opt', cntr, '.RData', sep = '')

# Load functions and libraries
cat('\nLoading packages and functions')

source('functions.R')
load('data/hf.RData')

# Set parallel computing plan
n.cores <- availableCores()-2
plan(multicore, workers = n.cores)
cat('\nParallel computing across', n.cores, 'cores')

# Define initial values for cost function
# and nloptr settings
cat('\n\n', rep('~', 60), sep='')
cat('\nSetting up constrained nonlinear minimization problem')
cat('\n            with "nloptr" to fit variogram models ...')
cat('\n\nsee https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/')
cat('\nfor information on nloptr and possible algorithms')
# Upper and lower bounds for constrained search
bounds <-
  tibble(
    cutoff.prop = c(1, 20),
    n.lags = c(15, 50),
    lag.start = c(1, 10),
    n.max = c(10, 50)
  )
init.vals <-
  tibble(
    cutoff.prop = runif(1, bounds$cutoff.prop[1], bounds$cutoff.prop[2]),
    n.lags = runif(1, bounds$n.lags[1], bounds$n.lags[2]),
    lag.start = runif(1, bounds$lag.start[1], bounds$lag.start[2]),
    n.max = runif(1, bounds$n.max[1], bounds$n.max[2]),
    algorithm = alg,
    maxeval = max.eval
  )
# Check init values are within bounds
pwalk(init.vals, ~{
  if(
    ..1 < bounds[1,1] |
    ..1 > bounds[2,1] |
    ..2 < bounds[1,2] |
    ..2 > bounds[2,2] |
    ..3 < bounds[1,3] |
    ..3 > bounds[2,3] |
    ..4 < bounds[1,4] |
    ..4 > bounds[2,4]
  ) {
    stop('One or more initial values out of bounds')
  }
})
# Print settings
cat(
  '\n\nnloptr settings:',
  '\nAlgorithm:       ', init.vals$algorithm,
  '\nMax evaluations: ', init.vals$maxeval,
  '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',
  '\n~~~~~~~~~~Initial values~~~~~~~~~~~~',
  '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',
  '\nLag cutoff:     ', init.vals$cutoff.prop,
  '\nNumber of lags: ', init.vals$n.lags,
  '\nLag shift:      ', init.vals$lag.start,
  '\nMax local pairs:', init.vals$n.max, '\n'
)

# Setup cost function and nloptr
f <- function(segment, v.mod) {
  # Define cost function to minimize with
  # input parameters (x)
  opt.fun <- function(x) {
    cost_function(
      shp.hf = shp.hf.crop[[segment]],
      cutoff.prop = x[1],
      n.lags = x[2],
      lag.start = x[3],
      n.max = x[4],
      model.vgrm = v.mod,
      interp.weight = 0.1,
      vgrm.weight = 0.9,
      segment = segment,
      verbose = T
    )
  }
  # Using nloptr optimization algorithms to minimize cost function
  # by tuning parameters (x) within box constrains
  # (up = upper bound, lb = lower bound)
  # with starting values x0
  # see https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/
  # for information on nloptr and possible algorithms
  nloptr(
    x0 = c(
      init.vals$cutoff.prop,
      init.vals$n.lags,
      init.vals$lag.start,
      init.vals$n.max
    ),
    eval_f = opt.fun,
    lb = c(
      bounds$cutoff.prop[1],
      bounds$n.lags[1],
      bounds$lag.start[1],
      bounds$n.max[1]
    ),
    ub = c(
      bounds$cutoff.prop[2],
      bounds$n.lags[2],
      bounds$lag.start[2],
      bounds$n.max[2]
    ),
    opts = list(
      print_level = 1,
      maxeval = init.vals$maxeval,
      algorithm = init.vals$algorithm
    )
  )
}

# Create array of variogram models to optimize
# for each segment and optimize in parallel
tic()
opt.grd <-
  tibble(expand.grid(
    segment = seg.names,
    v.mod = c('Sph', 'Exp', 'Gau', 'Bes', 'Lin', 'Cir'),
    stringsAsFactors = F
  )) %>%
  mutate(
    opt =
      future_map2(
        segment,
        v.mod,
        f,
        .progress = T,
        .options = furrr_options(seed = T)
      )
  )
toc()

# Decode optimization output and construct variograms
cat('\n\n', rep('~', 60), sep='')
cat('\nDecoding optimized solutions ...')
opt.vgrms <-
  opt.grd %>%
  pmap(
    ~decode_opt(
      model.vgrm = ..2,
      shp.hf = shp.hf.crop[[..1]],
      opt = ..3
    )
  )
# Add variograms and cost to tibble for easy manipulation and query
solns <-
  opt.grd %>%
  mutate(
    opt.exp.vgrm = map(opt.vgrms, ~.x[[1]]),
    opt.fit.vgrm = map(opt.vgrms, ~.x[[2]]),
    cost = map_dbl(opt.grd$opt, ~.x$objective)
  )
# Get solutions with minimum cost
opt.solns <-
  solns %>%
  group_by(segment) %>%
  slice(which.min(cost))

# Print lowest-cost solutions
cat('\nBest solutions for each segment:\n')
cat('\n', rep('+', 40), sep='')
pwalk(opt.solns, ~{
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
cat('\n', rep('+', 40), '\n', sep='')

# Summarise variograms
cat('\nVariogram summary:\n')
vgrm.summary <-
  tibble(
    segment = solns$segment,
    v.mod = solns$v.mod,
    cutoff.prop = map_dbl(solns$opt, ~.x$solution[1]),
    n.lags = map_dbl(solns$opt, ~.x$solution[2]),
    lag.start = map_dbl(solns$opt, ~.x$solution[3]),
    n.max = map_dbl(solns$opt, ~.x$solution[4]),
    sill = map_dbl(solns$opt.fit.vgrm, ~.x$psill),
    range = map_dbl(solns$opt.fit.vgrm, ~.x$range),
    itr = map_dbl(solns$opt, ~.x$iterations),
    cost = solns$cost,
  )
print(vgrm.summary)

# Kriging segments
cat('\n', rep('~', 60), sep='')
cat('\nKriging segments ...\n')
shp.interp.krige <-
  list(
    solns$segment,
    solns$opt.fit.vgrm,
    map(solns$opt, ~.x$solution[4])
  ) %>%
  pmap(~
    suppressWarnings(
      Krige(
       shp.hf = shp.hf.crop[[..1]],
       fitted.vgrm = ..2,
       shp.interp.grid = shp.grid.crop[[..1]],
       n.max = ..3,
       seg.name = ..1
      )
    )
  ) %>%
  set_names(solns$segment)

# Calculating interpolation difference
cat('\n', rep('~', 60), sep='')
cat('\nComputing differences ...\n')
shp.interp.diff <-
  pmap(list(shp.interp.krige, solns$segment, solns$v.mod), ~{
    suppressWarnings({
      cat('\nSegment:', ..2)
      cat('\nModel  :', ..3)
      interp_diff(..1, shp.interp.luca)
    })
  })

# Summarise interpolation differences
interp.diff.summary <-
  pmap_df(
    list(
      shp.interp.diff,
      solns$v.mod,
      solns$cost
    ),
    ~st_set_geometry(..1, NULL) %>%
    mutate(v.mod = ..2, cost = ..3, .before = est.sim),
    .id = 'segment'
  ) %>%
  group_by(segment, v.mod, cost) %>%
  drop_na() %>%
  summarise(
    n = n(),
    min = round(min(est.diff, na.rm = T)),
    max = round(max(est.diff, na.rm = T)),
    median = round(median(est.diff, na.rm = T)),
    IQR = round(IQR(est.diff, na.rm = T)),
    mean = round(mean(est.diff, na.rm = T)),
    sigma = round(sd(est.diff, na.rm = T)),
    .groups = 'drop'
  )
cat('\nHeat flow difference summary:\n')
print(interp.diff.summary)

cat('\n', rep('~', 60), sep='')
cat('\nSaving to: data/opt', cntr, '.RData', sep = '')
save(
  list = c(
    'solns',
    'opt.solns',
    'vgrm.summary',
    'shp.interp.diff',
    'interp.diff.summary'
  ),
  file = paste0('data/opt', cntr, '.RData')
)

cat('\n\nDone!')
cat('\n\n', rep('~', 60), '\n', sep='')