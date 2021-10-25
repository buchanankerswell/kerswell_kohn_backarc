#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop('Provide the optimization data file as the first arugment <opt*>.R', call.=FALSE)
} else if (length(args) == 1) {
  if(!file.exists(args[1])) {
    stop('That <opt*.R> file does not exist in data/', call.=FALSE)
  } else {
    cntr <- as.numeric(regmatches(args[1], regexec('[0-9]', args[1])))
    if(is.na(cntr)){
      cntr <- NULL
    }
  }
}

# Load functions and libraries
cat(rep('~', 60), '\n', sep='')
cat('Loading packages and functions ...\n\n')

source('functions.R')
load('data/hf.Rdata')
load(paste0('data/opt', cntr, '.RData'))
dir.create('figs/summary', showWarnings = F)
dir.create('figs/vgrms', showWarnings = F)
cat('\nSaving plots to: figs/summary/*', cntr, '.RData', sep = '')
cat('\nSaving plots to: figs/vgrms/*', cntr, '.png', sep = '')

# Visualize
cat('\n', rep('~', 60), sep='')
cat('\nVisualizing ...\n')

# Heat flow summary ridges plot
if(!file.exists('figs/summary/hfSummary.png')) {
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
    scale_fill_viridis_d() +
    scale_x_continuous(breaks = seq(0, 250, 50)) +
    scale_y_discrete(
      limits = rev(levels(as.factor(seg.names)))
    ) +
    theme_classic() +
    theme(
      legend.position = 'none',
      plot.margin = margin()
    )
# Save
cat('\nSaving plot to: figs/summary/hfSummary.png')
suppressWarnings(suppressMessages(
  ggsave(
    file = 'figs/summary/hfSummary.png',
    plot = p1,
    device = 'png',
    type = 'cairo',
    width = 6,
    height = 3.5
  )
))
# system('open figs/summary/hfSummary.png', wait = F)
}

# Plot all variogram models
cat('\nSaving plots to: figs/vgrms/')

# Draw composite plots
unique(solns$segment) %>%
walk(~{
  p <-
    filter(solns, segment == .x) %>%
    pmap(~{
      plot_vgrm(
        experimental.vgrm = ..4,
        fitted.vgrm = ..5,
        cost = ..7,
        v.mod = ..2,
        lineCol = 'firebrick'
      )
    }) %>%
  wrap_plots(nrow=3, ncol=2) +
  plot_annotation(title = paste(.x, 'variogram models'))
  ggsave(
    paste0('figs/vgrms/', str_replace_all(.x, ' ', ''), 'Vgrms', cntr, '.png'),
    plot = p,
    device = 'png',
    type = 'cairo',
    width = 6,
    height = 6
  )
  # system(
  #   paste0(
  #     'open figs/vgrms/',
  #     str_replace_all(.x, ' ', ''),
  #     'Vgrms',
  #     cntr,
  #     '.png'
  #   ),
  #   wait = F
  # )
})

# Summarise variogram models
p2 <-
  vgrm.summary %>%
  drop_na() %>%
  mutate('range' = range/1000) %>%
  filter(cost < (median(cost)+(5*IQR(cost))) & range < (median(range)+(5*IQR(range)))) %>%
  rename(
    'lag cutoff' = cutoff.prop,
    'number of lags' = n.lags,
    'lag shift' = lag.start,
    'max local pairs' = n.max,
    'range' = range,
    'sill' = sill,
  ) %>%
  pivot_longer(-c(
    v.mod,
    cost,
    segment,
    itr,
    vgrm.wt,
    vgrm.rmse,
    vgrm.cost,
    cv.wt,
    cv.rmse,
    cv.cost
  )) %>%
  ggplot(aes(x = cost, y = value, shape = v.mod, color = segment)) +
    geom_point(size = 2) +
    facet_wrap(
      ~name,
      ncol = 2,
      scales = 'free'
    ) +
    labs(
      y = NULL,
      x = bquote('Cost'~(mWm^-2)),
      color = 'Segment',
      shape = 'Variogram model'
    ) +
    scale_color_discrete_qualitative('Dark 3') +
    theme_classic() +
    theme(
      legend.box.margin = margin(),
      legend.text = element_text(size = 8),
      plot.margin = margin(),
      strip.background = element_rect(fill = 'grey90', color=NA)
    )
# Save
cat('\nSaving plot to: figs/summary/vgrmSummary', cntr, '.png', sep = '')
suppressWarnings(suppressMessages(
  ggsave(
    file = paste0('figs/summary/vgrmSummary', cntr, '.png'),
    plot = p2,
    device = 'png',
    type = 'cairo',
    width = 6,
    height = 6
  )
))
# system(paste0('open figs/summary/vgrmSummary', cntr, '.png'), wait = F)

# Interpolation differences
dif <-
  solns %>%
  filter(v.mod != 'Gau') %>%
  group_by(segment) %>%
  slice_min(cost)
p3 <-
  dif$shp.interp.diff %>%
  set_names(dif$segment) %>%
  map_df(
    ~st_set_geometry(.x, NULL),
    .id = 'segment'
  ) %>%
  ggplot() +
    stat_density_ridges(
      aes(x=est.diff, y=segment, group=segment, fill=factor(stat(quantile))),
      size = 0.3,
      scale = 2.5,
      geom = 'density_ridges_gradient',
      calc_ecdf = T,
      quantiles = 4,
      quantile_lines = T
    ) +
    labs(
      x = bquote('Interpolation estimate difference'~(mWm^-2)),
      y = NULL
    ) +
    scale_x_continuous(limits = c(-100, 100), breaks = seq(-100, 100, 50)) +
    scale_fill_viridis_d() +
    scale_y_discrete(limits = rev(levels(as.factor(seg.names)))) +
    theme_classic() +
    theme(
      legend.position = 'none',
      plot.margin = margin()
    )
# Save plot
cat('\nSaving plot to: figs/summary/interpDiffSummary', cntr, '.png', sep = '')
suppressWarnings(suppressMessages(
  ggsave(
    file = paste0('figs/summary/interpDiffSummary', cntr, '.png'),
    plot = p3,
    device = 'png',
    type = 'cairo',
    width = 6,
    height = 3.5
  )
))
# system(paste0('open figs/summary/interpDiffSummary', cntr, '.png'), wait = F)

p3b <-
  dif$shp.interp.diff %>%
  set_names(dif$segment) %>%
  map_df(
    ~st_set_geometry(.x, NULL),
    .id = 'segment'
  ) %>%
  ggplot() +
    stat_density_ridges(
      aes(x=sigma.diff, y=segment, group=segment, fill=factor(stat(quantile))),
      size = 0.3,
      scale = 2.5,
      geom = 'density_ridges_gradient',
      calc_ecdf = T,
      quantiles = 4,
      quantile_lines = T
    ) +
    labs(
      x = bquote('Interpolation'~sigma~'difference'~(mWm^-2)),
      y = NULL
    ) +
    scale_x_continuous(limits = c(-100, 100), breaks = seq(-100, 100, 50)) +
    scale_fill_viridis_d() +
    scale_y_discrete(limits = rev(levels(as.factor(seg.names)))) +
    theme_classic() +
    theme(
      legend.position = 'none',
      plot.margin = margin()
    )
# Save plot
cat('\nSaving plot to: figs/summary/interpSigmaDiffSummary', cntr, '.png', sep = '')
suppressWarnings(suppressMessages(
  ggsave(
    file = paste0('figs/summary/interpSigmaDiffSummary', cntr, '.png'),
    plot = p3b,
    device = 'png',
    type = 'cairo',
    width = 6,
    height = 3.5
  )
))
# system(paste0('open figs/summary/interpSigmaDiffSummary', cntr, '.png'), wait = F)
# nloptr trace
p4 <-
  opt.trace %>%
  group_by(segment, v.mod) %>%
  ggplot() +
    geom_path(aes(itr, cost, color = segment, group = segment)) +
    facet_wrap(~v.mod, scales = 'free') +
    labs(x = 'Iteration', y = 'Cost', color = NULL) +
    guides(color = guide_legend(nrow = 3)) +
    theme_classic() +
    theme(
      plot.margin = margin(),
      legend.position = 'bottom',
      legend.justification = 'left',
      legend.box.margin = margin(-10),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 9),
      legend.key.size = unit(0.8, 'lines'),
      strip.background = element_rect(fill = 'grey90')
    )
ggsave(
  file =
  paste0('figs/summary/optTrace', cntr, '.png'),
  plot = p4,
  device = 'png',
  type = 'cairo',
  width = 6,
  height = 3.5
)
# system(paste0('open figs/summary/optTrace', cntr, '.png'), wait = F)

cat('\nDone!\n')