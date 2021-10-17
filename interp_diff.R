# Load functions and libraries
cat(rep('~', 60), '\n', sep='')
cat('Loading packages and functions ...\n\n')
source('functions.R')
load('data/hf.RData')
cat(rep('~', 60), sep='')

# Optimizing Kriging
opt.fun <- function(x) {
  cost.function(
    shp.hf = shp.hf.crop[[3]],
    cutoff.prop = x[1],
    n.lags = x[2],
    lag.start = x[3],
    nmax = x[4],
    model.vgrm = x[5],
    verbose = T
  )
}

opt.solver <-
  lbfgs(
    x0 = c(3, 15, 1, 25, 0),
    fn = opt.fun,
    lower = c(1, 15, 1, 20, 0),
    upper = c(10, 115, 5, 500, 3),
    control = list(xtol_rel = 1e-6, maxeval = 2000)
  )

# Calculate experimental variogram
cat('\n', rep('~', 60), sep='')
cat('\nCalculating experimental variograms ...\n')

experimental.vgrms <- map(shp.hf.crop, experimental_vgrm)

# Arbitrary variogram models
cat('\nDefining variogram models')

model.vgrms <- map(seq_along(experimental.vgrms), ~vgm(NA, 'Sph', NA, NA))

# Fit experimental variograms
cat('\nFitting variograms')

fitted.vgrms <-
  map2(experimental.vgrms, model.vgrms, ~{
    suppressWarnings({
      fit.variogram(.x, model = .y, fit.method = 7)})
  })

# Save plot
cat('\nPlotting variogram models and saving to figs/vgrms.png')
cat('\n', rep('~', 60), sep='')

p <-
  list(experimental.vgrms, fitted.vgrms, seg.names) %>%
  pmap(plot_vgrm) %>%
  wrap_plots()

ggsave(
  'figs/vgrms.png',
  plot = p,
  device = 'png',
  type = 'cairo',
  width = 6,
  height = 6
)

cat('\n', rep('~', 60), sep='')
# Calculate refined experimental variogram
cat('\n', rep('~', 60), sep='')
cat('\nCalculating refined experimental variograms ...\n')

# Define refined experimental variogram parameters
cutoff.prop <- c(5, 6, 4, 5, 3, 8, 6, 3, 5, 15, 7, 3, 5)
n.lags <- c(15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15)
lag.start <- c(3, 1, 1, 2, 5, 1, 3, 1, 2, 1, 1, 1, 3)

pwalk(list(seg.names, cutoff.prop, n.lags, lag.start), ~{
  cat(
    '\nSegment:         ', ..1,
    '\nVariogram cutoff:', ..2, 
    '\nNumber of lags:  ', ..3,
    '\nLag shift:       ', ..4,
    '\n'
  )
})

experimental.vgrms <-
  pmap(
    list(
      shp.hf.crop,
      cutoff.prop,
      n.lags,
      lag.start
    ),
    ~experimental_vgrm(..1, ..2, ..3, ..4)
  ) %>%
  set_names(nm = seg.names)

# Arbitrary variogram models
cat('\nDefining variogram models')

model.vgrms <-
  map(
    seq_len(length(experimental.vgrms)),
    ~vgm(NA, 'Sph', NA, NA)
  ) %>%
  set_names(nm = seg.names)

# Fit experimental variograms
cat('\nFitting variograms')

fitted.vgrms <-
  map2(experimental.vgrms, model.vgrms, ~{
    suppressWarnings(fit.variogram(.x, model = .y, fit.method = 7))
  }) %>%
  set_names(seg.names)

# Save plot
cat('\nPlotting variogram models and saving to figs/vgrms.png')
cat('\n', rep('~', 60), sep='')

p <- pmap(list(experimental.vgrms, fitted.vgrms, seg.names), plot_vgrm) %>% wrap_plots()

ggsave(
  'figs/refined_vgrms.png',
  plot = p,
  device = 'png',
  type = 'cairo',
  width = 6,
  height = 6
)

cat('\n', rep('~', 60), sep='')
# Kriging segments
cat('\n', rep('~', 60), sep='')
cat('\nKriging segments ...\n')

shp.interp.krige <-
  list(shp.hf.crop, fitted.vgrms, shp.grid.crop) %>%
  pmap(~Krige(..1, ..2, ..3, nmax=200)) %>%
  set_names(seg.names)

# Calculating interpolation difference
cat('\nCalculating interpolation differences')

shp.interp.diff <-
  shp.interp.krige %>%
  map(~interp_diff(.x, shp.interp.luca)) %>%
  set_names(seg.names)

cat('\n', rep('~', 60), sep='')
cat('\nDone!\n')
cat('\n', rep('~', 60), sep='')

p1 <- 
  ggplot(data=shp.interp.diff[[5]]) +
  geom_sf(aes(color=est.krige)) +
  geom_sf(data=shp.buffer[[5]], color='white', fill=NA) +
  geom_sf(data=shp.hf.crop[[5]], color='deeppink', size=0.5, shape=15)
ggplot(data=shp.interp.diff[[5]])+geom_sf(aes(color=est.similarity)) -> p2
ggplot(data=shp.interp.diff[[5]])+geom_sf(aes(color=est.diff)) -> p3
p1 / p2 / p3