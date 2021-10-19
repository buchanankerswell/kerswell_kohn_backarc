# Load functions and libraries
cat(rep('~', 60), '\n', sep='')
cat('Loading packages and functions ...\n\n')

source('functions.R')
load('data/hf.RData')

# Set parallel computing plan
plan(multisession, workers = 7)

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
  # Using BOBYQA or COBYLA optimization algorithm to minimize cost function
  # by tuning parameters (x) within box constrains
  # (up = upper bound, lb = lower bound)
  # with starting values x0
  # see https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/
  # for information on nloptr and the BOBYQA algorithm
  nloptr(
    x0 = c(8, 30, 5, 20),
    eval_f = opt.fun,
    lb = c(1, 15, 1, 10),
    ub = c(20, 50, 10, 50),
    opts = list(
      print_level = 1,
#      ftol_rel = 1e-4,
#      xtol_rel = 1e-4,
      maxeval = 500,
#      algorithm = 'NLOPT_LN_BOBYQA'
#      algorithm = 'NLOPT_LN_COBYLA'
#      algorithm = 'NLOPT_GN_DIRECT_L'
      algorithm = 'NLOPT_LN_SBPLX'
    )
  )
}

# Find good variogram models by constrained
# non-linear optimization
cat('\n', rep('~', 60), sep='')
cat('\nOptimizing variogram parameters ...')
cat(
  '\nnloptr settings:',
  '\nAlgorithm:         NLOPT_LN_SBPLX',
  '\nMax evaluations:   500',
  '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',
  '\n~~~~~~~~~Initial values~~~~~~~~~~~',
  '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',
  '\n  Lag cutoff:      8',
  '\n  Number of lags:  30',
  '\n  Lag shift:       5',
  '\n  Max local pairs: 20\n'
)

# Create array of variogram models to optimize for each segment
# and optimize in parallel
tic()
opt.grd <-
  tibble(expand.grid(c('Sph', 'Exp', 'Gau', 'Bes', 'Lin', 'Cir'), seg.names)) %>%
  rename(v.mod = Var1, segment = Var2) %>%
  mutate(across(everything(), as.character)) %>%
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
cat('\n', rep('~', 60), sep='')
cat('\nDecoding optimized solutions ...')
opt.vgrms <-
  opt.grd %>%
  pmap(
    ~decode_opt(
      model.vgrm = ..1,
      shp.hf = shp.hf.crop[[..2]],
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

cat('\nBest solutions for each segment:\n')
cat('\n', rep('+', 40), sep='')
pwalk(opt.solns, ~{
  cat(
    '\nSegment:         ', ..2,
    '\nVariogram cutoff:', ..3$solution[1], 
    '\nNumber of lags:  ', ..3$solution[2],
    '\nLag shift:       ', ..3$solution[3],
    '\nMax points:      ', ..3$solution[4],
    '\nVariogram model: ', ..1,
    '\nCost:            ', ..3$objective,
    '\n'
  )
})
cat('\n', rep('+', 40), '\n', sep='')

# Summarise variograms
cat('\nVariogram summary:\n')
vgrm.summary <-
  tibble(
    segment = opt.solns$segment,
    model.vgrm = opt.solns$v.mod,
    cutoff.prop = map_dbl(opt.solns$opt, ~.x$solution[1]),
    n.lags = map_dbl(opt.solns$opt, ~.x$solution[2]),
    lag.start = map_dbl(opt.solns$opt, ~.x$solution[3]),
    n.max = map_dbl(opt.solns$opt, ~.x$solution[4]),
    sill = map_dbl(opt.solns$opt.fit.vgrm, ~.x$psill),
    range = map_dbl(opt.solns$opt.fit.vgrm, ~.x$range),
    itr = map_dbl(opt.solns$opt, ~.x$iterations),
    cost = opt.solns$cost,
  )
print(vgrm.summary)

# Plot optimum variogram models
cat('\nSaving variogram plots to: figs/vgrms/\n')
dir.create('figs/vgrms', showWarnings = F)
plts <- pmap(solns, ~plot_vgrm(..4, ..5, ..6, ..1, lineCol = 'firebrick'))
list(seq(1,73,6), seq(6,78,6), seg.names) %>%
pwalk(~{
  p <-
    wrap_plots(plts[..1:..2], nrow=3, ncol=2) +
    plot_annotation(title = paste(..3, 'variogram models'))
  ggsave(
    paste0('figs/vgrms/', str_replace_all(..3, ' ', ''), 'Vgrms.png'),
    plot = p,
    device = 'png',
    type = 'cairo',
    width = 6,
    height = 6
  )
  }
)

# Kriging segments
cat('\n', rep('~', 60), sep='')
cat('\nKriging segments ...\n')
shp.interp.krige <-
  list(
    shp.hf.crop,
    opt.solns$opt.fit.vgrm,
    shp.grid.crop,
    map(opt.solns$opt, ~.x$solution[4])
  ) %>%
  pmap(~suppressWarnings(Krige(..1, ..2, ..3, ..4)))

# Calculating interpolation difference
cat('\n', rep('~', 60), sep='')
cat('\nComputing differences ...\n')
shp.interp.diff <-
  shp.interp.krige %>%
  map(~suppressWarnings(interp_diff(.x, shp.interp.luca)))

# Summarise interpolation differences
interp.diff.summary <-
  shp.interp.diff %>%
  map_df(
    ~st_set_geometry(.x, NULL),
    .id = 'segment'
  ) %>%
  group_by(segment) %>%
  drop_na() %>%
  summarise(
    n = n(),
    min = round(min(est.diff, na.rm = T)),
    max = round(max(est.diff, na.rm = T)),
    median = round(median(est.diff, na.rm = T)),
    IQR = round(IQR(est.diff, na.rm = T)),
    mean = round(mean(est.diff, na.rm = T)),
    sigma = round(sd(est.diff, na.rm = T))
  )
cat('\nHeat flow difference summary:\n')
print(interp.diff.summary)

cat('\n', rep('~', 60), sep='')
cat('\nSaving to: data/opt.RData')
save(
  list = c(
    'solns',
    'opt.solns',
    'vgrm.summary',
    'shp.interp.diff',
    'interp.diff.summary'
  ),
  file = 'data/opt2.RData'
)

cat('\n\nDone!')
cat('\n\n', rep('~', 60), sep='')