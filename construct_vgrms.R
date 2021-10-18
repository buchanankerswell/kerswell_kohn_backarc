# Load functions and libraries
cat(rep('~', 60), '\n', sep='')
cat('Loading packages and functions ...\n\n')

source('functions.R')
load('data/hf.RData')

cat(rep('~', 60), sep='')

# Find good variogram models by constrained
# non-linear optimization
cat('\n', rep('~', 60), sep='')
cat('\nOptimizing variogram parameters ...')

# Set parallel computing plan
plan(multisession, workers = availableCores()-1)

library(tictoc)
tick()

solns <-
  tibble(
    expand.grid(
      c(
        'Sph',
        'Exp',
        'Gau',
        'Bes',
        'Lin',
        'Cir'
      ),
      seg.names
    )
  ) %>%
  rename(v.mod = Var1, segment = Var2) %>%
  mutate(across(everything(), as.character)) %>%
  mutate(
    opt =
      future_map2(segment, v.mod, ~{
        # Define cost function to minimize with
        # input parameters (x)
        opt.fun <- function(x) {
          cost_function(
            shp.hf = shp.hf.crop[[..1]],
            cutoff.prop = x[1],
            n.lags = x[2],
            lag.start = x[3],
            n.max = x[4],
            model.vgrm = ..2,
            interp.weight = 0.3,
            vgrm.weight = 0.7,
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
          x0 = c(10, 25, 5, 20),
          eval_f = opt.fun,
          lb = c(1, 15, 1, 10),
          ub = c(20, 50, 20, 50),
          opts = list(
            print_level = 1,
            ftol_rel = 1e-6,
            xtol_abs = 0,
            maxeval = 1000,
            maxtime = 3600,
#             algorithm = 'NLOPT_LN_BOBYQA'
            algorithm = 'NLOPT_LN_COBYLA'
          )
        )
     }, .progress = T)
  )

toc()

cat('\n', rep('~', 60), sep='')

# Decode optimization output and construct variograms
cat('\nDecoding optimized solutions ...')
opt.vgrms <-
  solns %>%
  pmap(
    ~decode_opt(
      model.vgrm = ..1,
      shp.hf = shp.hf.crop[[..2]],
      opt = ..3
    )
  )

# Add variograms and cost to tibble for easy manipulation and query
solns <-
  solns %>%
  mutate(
    opt.exp.vgrm = map(opt.vgrms, ~.x[[1]]),
    opt.fit.vgrm = map(opt.vgrms, ~.x[[2]]),
    cost = map_dbl(solns$opt, ~.x$objective)
  )

# Find best solutions for each SZ segment
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


cat('\nVariogram summary:')

shp.area <- (shp.buffer %>% map(st_area))/e+12 # Mkm^2
shp.length <- (shp.segs %>% map(st_length) %>% as.vector())/1000 # km

vgrm.summary <-
  mutate(
    seg.length = shp.length, # km
    domain.area = shp.area, # Mkm^2
    n.cv = ,
    cv.density = n.cv/domain.area, # n/Mkm^2
    cutoff.prop = ,
    n.lags = ,
    lag.start = ,
    n.max = ,
    model.vgrm = ,
    sill = , # mW/m^-2
    range = , # km
    cv.rmse = ,
    cost = ,
  )

print(variogram.summary)

cat('\nSaving variogram plots to: figs/vgrms/\n')

dir.create('figs/vgrms', showWarnings = F)

plts <- pmap(solns, ~plot_vgrm(..4, ..5, ..6, ..1, lineCol = 'firebrick'))
list(seq(1,73,6), seq(6,78,6), seg.names) %>%
pwalk(~{
  p <-
    wrap_plots(plts[..1:..2], nrow=3, ncol=2) +
    plot_annotation(title = paste(..3, 'variogram models'))
  ggsave(
    paste0('figs/vgrms/', str_replace(..3, ' ', ''), 'Vgrms.png'),
    plot = p,
    device = 'png',
    type = 'cairo',
    width = 6,
    height = 6
  )
  }
)

cat('\n', rep('~', 60), sep='')

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
  pmap(~Krige(..1, ..2, ..3, ..4))

# Calculating interpolation difference
cat('\nCalculating interpolation differences')

shp.interp.diff <-
  shp.interp.krige %>%
  map(~interp_diff(.x, shp.interp.luca))

hf.diff.summary <-
  hf.diff %>%
  map_df(
    ~st_set_geometry(.x, NULL),
    .id = 'segment'
  ) %>%
  group_by(segment) %>%
  summarise(
    n = n(),
    min = round(min(hf.diff)),
    max = round(max(hf.diff)),
    median = round(median(hf.diff)),
    IQR = round(IQR(hf.diff)),
    mean = round(mean(hf.diff)),
    sigma = round(sd(hf.diff))
  )

cat('\nHeat flow difference summary:\n')
print(hf.diff.summary)

cat('\n', rep('~', 60), sep='')
cat('\nDone!')
cat('\n', rep('~', 60), sep='')