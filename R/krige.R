#!/usr/bin/env Rscript

# Capture output
sink(file = paste0('data/log-', Sys.Date()), append = T, type = 'output', split = T)

# Testing arguments
args <- commandArgs(trailingOnly = TRUE)
cat('\n', rep('~', 80), sep = '')
if(length(args) == 0) {
  cat('\nNo arguments passed to krige.R')
  max.eval <- 50
  alg <- 'NLOPT_LN_SBPLX'
  iwt <- 0.5
  vwt <- 0.5
  n.fold <- NULL
  n.cores <- future::availableCores() - 2
  cat('\nUsing defaults')
  cat('\nMax iterations           :', max.eval)
  cat('\nAlgorithm                :', alg)
  cat('\nInterpolation cost weight:', iwt)
  cat('\nVariogram cost weight    :', vwt)
  cat('\nk-folds for CV           :', vwt)
  cat('\nVariogram cost weight    :', vwt)
} else if(length(args) != 0) {
  max.eval <- suppressWarnings(as.numeric(args[1]))
  if(is.na(max.eval)) {
    max.eval <- 50
    cat('\nPassed non-numeric argument for max iterations!')
    cat('\nDefaulting to', max.eval)
  } else if(max.eval <= 0) {
    max.eval <- 50
    cat('\nMax iterations cannot be negative or zero!')
    cat('\nDefaulting to', max.eval)
  }
  alg <- suppressWarnings(as.numeric(args[2]))
  if(is.na(alg)) {
    cat('\nAlgorithm passed is not recognized!')
    cat('\nDefaulting to NLOPT_LN_SBPLX')
    alg <- 'NLOPT_LN_SBPLX'
  } else {
    if(args[2] == 1) {
      alg <- 'NLOPT_GN_DIRECT_L'   # Global search
    } else if(args[2] == 2) {
      alg <- 'NLOPT_LN_SBPLX'      # Local without gradients
    } else if(args[2] == 3) {
      alg <- 'NLOPT_LN_NELDERMEAD' # Local without gradients
    } else if(args[2] == 4) {
      alg <- 'NLOPT_LN_BOBYQA'     # Local without gradients
    } else if(args[2] == 5) {
      alg <- 'NLOPT_LN_COBYLA'     # Local without gradients
    } else {
      cat('\nAlgorithm passed is not recognized!')
      cat('\nDefaulting to NLOPT_LN_SBPLX')
      alg <- 'NLOPT_LN_SBPLX'
    }
  }
  iwt <- suppressWarnings(as.numeric(args[3]))
  if(is.na(iwt)) {
    iwt <- 0.5
    cat('\nPassed non-numeric argument for interpolation weight!')
    cat('\nDefaulting to', iwt)
  }
  vwt <- suppressWarnings(as.numeric(args[4]))
  if(is.na(vwt)) {
    vwt <- 0.5
    cat('\nPassed non-numeric argument for variogram weight!')
    cat('\nDefaulting to', vwt)
  }
  if((iwt+vwt) != 1) {
    iwt <- 0.5
    vwt <- 0.5
    cat('\nInterpolation and variogram weights must add to one!')
    cat('\nDefaulting to', iwt, 'and', vwt)
  }
  n.cores <- suppressWarnings(as.numeric(args[5]))
  if(is.na(n.cores)) {
    n.cores <- future::availableCores() - 2
    cat('\nPassed non-numeric argument for number of cores!')
    cat('\nDefaulting to', n.cores)
  } else if(n.cores > future::availableCores()) {
    n.cores <- future::availableCores() - 2
    cat('\nToo many cores!')
    cat('\nDefaulting to', n.cores)
  } else if(n.cores < 0) {
    n.cores <- future::availableCores() - 2
    cat('\nCannot run negative cores!')
    cat('\nDefaulting to', n.cores)
  }
  n.fold <- suppressWarnings(as.numeric(args[6]))
  if(is.na(n.fold)) {
    n.fold <- NULL
    cat('\nPassed non-numeric argument for k-fold!')
    cat('\nDefaulting to leave-one-out cross-validation')
  } else {
    if(n.fold < 0) {
      n.fold <- NULL
      cat('\nk-fold cannot be negative or zero!')
      cat('\nDefaulting to leave-one-out cross-validation')
    } else if(n.fold == 0) {
      n.fold <- NULL
      cat('\nDefaulting to leave-one-out cross-validation')
    }
  }
}

cat('\nSaving results to: data/opt.RData', sep = '')

# Load functions and libraries
cat('\nLoading packages and functions ...')
source('R/functions.R')
load('data/hf.RData')

# Set parallel computing plan
plan(multicore, workers = n.cores)
cat('\nParallel computing across', n.cores, 'cores')

# Define initial values for cost function
# and nloptr settings
cat('\n', rep('~', 80), sep = '')
cat('\nSetting up constrained nonlinear minimization problem')
cat('\nwith "nloptr" to fit variogram models ...')
cat('\nSee https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/')
cat('\nfor information on nloptr and possible algorithms')

# Upper and lower bounds for constrained search
bounds <-
  tibble(
    cutoff.prop = c(1, 10),
    n.lags = c(15, 50),
    lag.start = c(1, 10),
    n.max = c(5, 50)
  )
init.vals <-
  tibble(
    cutoff.prop = 3,
    n.lags = 20,
    lag.start = 1,
    n.max = 8,
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
  '\n',
  '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',
  '\n~~~~~                         ~~~~~~',
  '\n~~~~~     nloptr settings     ~~~~~~',
  '\n~~~~~                         ~~~~~~',
  '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',
  '\nAlgorithm       : ', init.vals$algorithm,
  '\nMax evaluations : ', init.vals$maxeval,
  '\nk-fold CV       : ', ifelse(is.null(n.fold), 'Leave-one-out', n.fold),
  '\nInterp weight   : ', iwt,
  '\nVgrm weight     : ', vwt,
  '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',
  '\n~~~~~                         ~~~~~~',
  '\n~~~~~      Initial values     ~~~~~~',
  '\n~~~~~                         ~~~~~~',
  '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',
  '\nLag cutoff      :', init.vals$cutoff.prop,
  '\nNumber of lags  :', init.vals$n.lags,
  '\nLag shift       :', init.vals$lag.start,
  '\nMax local pairs :', init.vals$n.max,
  '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',
  '\n~~~~~                         ~~~~~~',
  '\n~~~~~    Heavy Computation!   ~~~~~~',
  '\n~~~~~                         ~~~~~~',
  '\n~~~~~      This may take      ~~~~~~',
  '\n~~~~~     minutes to hours    ~~~~~~',
  '\n~~~~~                         ~~~~~~',
  '\n~~~~~    Have a sip of tea    ~~~~~~',
  '\n~~~~~                         ~~~~~~',
  '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n'
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
      n.fold = n.fold,
      model.vgrm = v.mod,
      interp.weight = iwt,
      vgrm.weight = vwt,
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
      print_level = 0,
      maxeval = init.vals$maxeval,
      algorithm = init.vals$algorithm
    )
  )
}

# Create array of variogram models to optimize
# for each segment in parallel
optimized.interpolations <-
  tibble(expand.grid(
    segment = seg.names,
    v.mod = c('Sph', 'Exp', 'Lin', 'Bes'),
    stringsAsFactors = F
  )) %>%
  arrange(segment, v.mod) %>%
  mutate(
    opt =
      suppressMessages(future_map2(
        segment,
        v.mod,
        possibly(f, otherwise = NULL),
        .progress = T,
        .options = furrr_options(seed = T)
      ))
  )

# Save output
sink()

# Capture output
sink(file = paste0('data/log-', Sys.Date()), append = T, type = 'output', split = T)

# Check for NULL results from errors
if(any(map_lgl(optimized.interpolations$opt, is.null))) {
  cat('\n', rep('~!', 30), sep = '')
  cat(
    '\n\nOptimization failed for',
    sum(map_lgl(optimized.interpolations$opt, is.null)),
    'runs:\n'
  )
  print(optimized.interpolations[map_lgl(optimized.interpolations$opt, is.null),])
  cat('\nContinuing computation for good runs ...')
  cat('\n', rep('~!', 30), sep = '')
}

# Drop missing values
optimized.interpolations <- drop_na(optimized.interpolations)

# Parse nloptr history
read_file <- function(fpath, ...) {
  con <- file(fpath)
  on.exit(close(con))
  readLines(con, ...)
}
raw.trace <- read_file(paste0('data/log-', Sys.Date()))
segs.trace <-
  raw.trace[grepl('^Segment: ', raw.trace)] %>%
  map_chr(~gsub('Segment: ', '', .x))
vgrm.model.trace <-
  raw.trace[grepl('^Variogram model:', raw.trace)] %>%
  map_chr(~gsub('Variogram model: ', '', .x))
vgrm.weight.trace <-
  raw.trace[grepl('^Variogram weight', raw.trace)] %>%
  map_dbl(~as.numeric(gsub('Variogram weight: ', '', .x)))
vgrm.rmse.trace <-
  raw.trace[grepl('^Variogram rmse', raw.trace)] %>%
  map_dbl(~as.numeric(gsub('Variogram rmse: ', '', .x)))
vgrm.cost.trace <-
  raw.trace[grepl('^Variogram cost', raw.trace)] %>%
  map_dbl(~as.numeric(gsub('Variogram cost: ', '', .x)))
cv.weight.trace <-
  raw.trace[grepl('^Interpolation weight', raw.trace)] %>%
  map_dbl(~as.numeric(gsub('Interpolation weight: ', '', .x)))
cv.rmse.trace <-
  raw.trace[grepl('^Interpolation rmse', raw.trace)] %>%
  map_dbl(~as.numeric(gsub('Interpolation rmse: ', '', .x)))
cv.cost.trace <-
  raw.trace[grepl('^Interpolation cost', raw.trace)] %>%
  map_dbl(~as.numeric(gsub('Interpolation cost: ', '', .x)))
cost.trace <-
  raw.trace[grepl('^Cost', raw.trace)] %>%
  map_dbl(~as.numeric(gsub('Cost: ', '', .x)))
opt.trace <-
  tibble(
    segment = segs.trace,
    v.mod = vgrm.model.trace,
    vgrm.wt = vgrm.weight.trace,
    vgrm.rmse = vgrm.rmse.trace,
    vgrm.cost = vgrm.cost.trace,
    cv.wt = cv.weight.trace,
    cv.rmse = cv.rmse.trace,
    cv.cost = cv.cost.trace,
    cost = cost.trace
  ) %>%
  group_by(segment, v.mod) %>%
  mutate(itr = row_number(), .before = segment) %>%
  ungroup()

cat('\nnloptr trace ...')
cat('\n', rep('-', 40), '\n', sep = '')
print(opt.trace)
cat(rep('-', 40), '\n', sep = '')

# Decode optimization output and construct variograms
cat('\n\n', rep('~', 80), sep = '')
cat('\nDecoding optimized solutions ...')
opt.vgrms <-
  optimized.interpolations %>%
  pmap(
    possibly(
      ~decode_opt(
        model.vgrm = ..2,
        shp.hf = shp.hf.crop[[..1]],
        opt = ..3
      ),
      otherwise = NULL
    )
  )
# Add variograms and cost to tibble for easy manipulation and query
solns <-
  optimized.interpolations %>%
  mutate(
    opt.exp.vgrm = map(opt.vgrms, ~.x[[1]]),
    opt.fit.vgrm = map(opt.vgrms, ~.x[[2]]),
    cost = map_dbl(optimized.interpolations$opt, ~.x$objective)
  )

# Check for NULL results from errors
if(any(map_lgl(solns$opt.exp.vgrm, is.null))) {
  cat('\n', rep('~!', 30), sep = '')
  cat(
    '\n\nVariogram fitting failed during opt decoding for',
    sum(map_lgl(solns$opt.exp.vgrm, is.null)),
    'runs:\n'
  )
  print(solns[map_lgl(solns$opt.exp.vgrm, is.null),])
  cat('\nContinuing computation for good runs ...')
  cat('\n', rep('~!', 30), sep = '')
}

# Get solutions with minimum cost
opt.solns <-
  solns %>%
  drop_na() %>%
  group_by(segment) %>%
  slice(which.min(cost))

# Print lowest-cost solutions
cat('\nBest solutions for each segment:\n')
cat('\n', rep('-', 40), sep = '')
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
cat('\n', rep('-', 40), '\n', sep = '')

# Drop missing values
solns <- drop_na(solns)

# Summarise variograms
cat('\nVariogram summary:\n')
opt.trace.end <-
  opt.trace %>%
  group_by(segment, v.mod) %>%
  slice_tail() %>%
  ungroup()
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
  ) %>%
  left_join(opt.trace.end, by = c('segment', 'v.mod')) %>%
  left_join(rmse.luca, by = c('segment'))
print(vgrm.summary)

# Kriging segments
cat('\n', rep('~', 80), sep = '')
cat('\nKriging segments ...\n')
shp.interp.krige <-
  list(
    solns$segment,
    solns$opt.fit.vgrm,
    map(solns$opt, ~.x$solution[4])
  ) %>%
  pmap(
    possibly(
      ~Krige(
       shp.hf = shp.hf.crop[[..1]],
       fitted.vgrm = ..2,
       shp.interp.grid = shp.grid.crop[[..1]],
       n.max = ..3,
       seg.name = ..1
      ),
      otherwise = NULL
    )
  ) %>%
  set_names(solns$segment)

# Filter negative heat flow estimates
shp.interp.krige <- map(shp.interp.krige, ~filter(.x, est.krige >= 0))

# Check for NULL results from errors
if(any(map_lgl(shp.interp.krige, is.null))) {
  cat('\n', rep('~!', 30), sep = '')
  cat(
    '\n\nVariogram fitting failed during opt decoding for',
    sum(map_lgl(shp.interp.krige, is.null)),
    'runs:\n'
  )
  indx <- map_lgl(shp.interp.krige, is.null)
  cat('\nSegments:\n')
  writeLines(solns$segment[indx])
  cat('\nVariogram models:\n')
  print(solns$opt.fit.vgrm[indx])
  cat('\nContinuing computation for good runs ...')
  cat('\n', rep('~!', 30), sep = '')
}

# Remove missing values
solns <-
  solns %>%
  mutate(shp.interp.krige, .before = cost) %>%
  drop_na()

# Calculating interpolation difference
cat('\n', rep('~', 80), sep = '')
cat('\nComputing differences ...\n')
shp.interp.diff <-
  solns %>%
  pmap(
    possibly(
      ~{
        interp_diff(
          shp.interp.krige = ..6,
          shp.interp.sim = shp.interp.luca
        )
      },
      otherwise = NULL
    )
  ) %>% set_names(solns$segment)

# Check for NULL results from errors
if(any(map_lgl(shp.interp.diff, is.null))) {
  cat('\n', rep('~!', 30), sep = '')
  cat(
    '\n\nVariogram fitting failed during opt decoding for',
    sum(map_lgl(shp.interp.diff, is.null)),
    'runs:\n'
  )
  indx <- map_lgl(shp.interp.diff, is.null)
  cat('\nSegments:\n')
  writeLines(solns$segment[indx])
  cat('\nVariogram models:\n')
  print(solns$opt.fit.vgrm[indx])
  cat('\nContinuing computation for good runs ...')
  cat('\n', rep('~!', 30), sep = '')
}

# Remove missing values
solns <-
  solns %>%
  mutate(shp.interp.diff, .before = cost) %>%
  drop_na()

# Summarise interpolation differences
interp.diff.summary <-
  solns %>%
  pmap_df(
    ~st_set_geometry(..7, NULL) %>%
    mutate(segment = ..1, v.mod = ..2, cost = ..8, .before = est.sim)
  ) %>%
  group_by(segment, v.mod, cost) %>%
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

# Drop krige column since it is a duplicate of diff column
solns <- select(solns, -shp.interp.krige)

# Calculate direct RMSE for Kriging predictions (to compare with cross-validation RMSE)
vgrm.summary$dir.rmse <-
  map2_dbl(solns$segment, solns$shp.interp.diff, ~itp_rmse(.x, .y, 'krg'))

cat('\n', rep('~', 80), sep = '')
cat('\nSaving to: data/opt.RData', sep = '')
save(
  list = c(
    'solns',
    'opt.solns',
    'opt.trace',
    'vgrm.summary',
    'interp.diff.summary'
  ),
  file = 'data/opt.RData'
)

cat('\nkrige.R complete!\n\n')
sink()