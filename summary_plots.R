# Load functions and libraries
cat(rep('~', 60), '\n', sep='')
cat('Loading packages and functions ...\n\n')

source('functions.R')
load('data/hf.Rdata')
load('data/opt.RData')
dir.create('figs/summary', showWarnings = F)

cat(rep('~', 60), sep='')

# Visualize
cat('\n', rep('~', 60), sep='')
cat('\nVisualizing ...')

# Heat flow summary ridges plot
p1 <-
  shp.hf.crop %>%
  map_df(
    ~st_set_geometry(.x, NULL),
    .id = 'segment'
  ) %>%
  group_by(segment) %>%
  ggplot() +
    stat_density_ridges(
      aes(x=hf, y=segment, group=segment, fill=factor(stat(quantile))),
      size = 0.3,
      scale = 3.5,
      geom='density_ridges_gradient',
      calc_ecdf=T,
      quantiles=4,
      quantile_lines=T
    ) +
    labs(
      x = bquote('Heat flow'~(mWm^-2)),
      y = NULL
    ) +
    scale_x_continuous(limits = c(-25, 250)) +
    scale_fill_viridis_d() +
    scale_y_discrete(
      limits = rev(levels(as.factor(seg.names)))
    ) +
    theme_classic() +
    theme(
      legend.position = 'none',
      plot.margin = margin()
    )
# Save
cat('\nSaving plot to: figs/summary/hf_summary.png')
suppressWarnings(suppressMessages(
  ggsave(
    file = 'figs/summary/hf_summary.png',
    plot = p1,
    device = 'png',
    type = 'cairo',
    width = 6,
    height = 3.5
  )
))
system('open figs/summary/hf_summary.png', wait = F)

# Summarise variogram models
p2 <-
  vgrm.summary %>%
  mutate('range' = range/1000) %>%
  rename(
    'lag~cutoff' = cutoff.prop,
    'number~of~lags' = n.lags,
    'lag~shift' = lag.start,
    'max~local~pairs' = n.max,
    'range~(km)' = range,
    'sill~(mWm^-2)' = sill,
  ) %>%
  pivot_longer(-c(model.vgrm, cost)) %>%
  ggplot(aes(x = cost, y = value)) +
    geom_point() +
    facet_wrap(~name, scales = 'free', labeller = label_parsed) +
    labs(y = NULL, x = bquote('Cost'~(mWm^-2))) +
    theme_classic() +
    theme(
      strip.background = element_rect(fill = 'grey90', color=NA)
    )
# Save
cat('\nSaving plot to: figs/summary/vgrm_summary.png')
suppressWarnings(suppressMessages(
  ggsave(
    file = 'figs/summary/vgrm_summary.png',
    plot = p2,
    device = 'png',
    type = 'cairo',
    width = 6,
    height = 3.5
  )
))
system('open figs/summary/vgrm_summary.png', wait = F)

# Interpolation differences
p3 <-
  shp.interp.diff %>%
  map_df(
    ~st_set_geometry(.x, NULL),
    .id = 'segment'
  ) %>%
  ggplot() +
    stat_density_ridges(
      aes(x=est.diff, y=segment, group=segment, fill=factor(stat(quantile))),
      size = 0.3,
      scale = 1.5,
      geom = 'density_ridges_gradient',
      calc_ecdf = T,
      quantiles = 4,
      quantile_lines = T
    ) +
    labs(
      x = bquote('Interpolation difference'~(mWm^-2)),
      y = NULL
    ) +
    scale_x_continuous(limits = c(-150, 150), breaks = seq(-150, 150, 50)) +
    scale_fill_viridis_d() +
    scale_y_discrete(limits = rev(levels(as.factor(seg.names)))) +
    theme_classic() +
    theme(
      legend.position = 'none',
      plot.margin = margin()
    )
# Save plot
cat('\nSaving plot to: figs/summary/interp_diff_summary.png')
suppressWarnings(suppressMessages(
  ggsave(
    file = 'figs/summary/interp_diff_summary.png',
    plot = p3,
    device = 'png',
    type = 'cairo',
    width = 6,
    height = 3.5
  )
))
system('open figs/summary/interp_diff_summary.png', wait = F)

cat('\nDone!')
cat('\n', rep('~', 60), sep='')