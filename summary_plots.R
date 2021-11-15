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
dir.create('figs/upperPlate', showWarnings = F)
cat('\nSaving plots to: figs/summary/')
cat('\nSaving plots to: figs/vgrms/')

# Visualize
cat('\n', rep('~', 60), sep='')
cat('\nVisualizing ...\n')

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
      quantile_lines=T,
      bandwidth = 10
    ) +
    labs(
      x = bquote('Heat flow'~(mWm^-2)),
      y = NULL
    ) +
    scale_fill_viridis_d() +
    scale_y_discrete(
      limits = rev(levels(as.factor(seg.names)))
    ) +
    theme_classic(base_size = 12) +
    theme(
      panel.background = element_rect(fill = 'grey50', color = NA),
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
})

# Summarise variogram models
p2 <-
  vgrm.summary %>%
  drop_na() %>%
  select(-c(n.obs.sim, rmse.sim)) %>%
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
  ggplot(aes(x = cost, y = value, shape = v.mod)) +
    geom_point(size = 2) +
    facet_wrap(
      ~name,
      ncol = 3,
      scales = 'free'
    ) +
    guides(shape = guide_legend(nrow = 1)) +
    labs(
      y = NULL,
      x = bquote('Cost'~(mWm^-2)),
      shape = 'Variogram model'
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = 'bottom',
      legend.box.margin = margin(-10, 0, 0, 0),
      panel.background = element_rect(fill = 'grey50', color = NA),
      plot.margin = margin()
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
    height = 4.5
  )
))

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
      quantile_lines = T,
      bandwidth = 5,
      rel_min_height = 0.01
    ) +
    labs(
      x = bquote('Interpolation estimate difference'~(mWm^-2)),
      y = NULL
    ) +
    coord_cartesian(xlim = c(-75, 75)) +
    scale_fill_viridis_d() +
    scale_y_discrete(limits = rev(levels(as.factor(seg.names)))) +
    theme_classic(base_size = 12) +
    theme(
      panel.background = element_rect(fill = 'grey50', color = NA),
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
      quantile_lines = T,
      bandwidth = 5,
      rel_min_height = 0.01
    ) +
    labs(
      x = bquote('Interpolation'~sigma~'difference'~(mWm^-2)),
      y = NULL
    ) +
    coord_cartesian(xlim = c(-75, 75)) +
    scale_fill_viridis_d() +
    scale_y_discrete(limits = rev(levels(as.factor(seg.names)))) +
    theme_classic(base_size = 12) +
    theme(
      panel.background = element_rect(fill = 'grey50', color = NA),
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
# nloptr trace
p4 <-
  opt.trace %>%
  group_by(segment, v.mod) %>%
  ggplot() +
    geom_path(aes(itr, cost, group = segment)) +
    facet_wrap(~v.mod, nrow = 2, scales = 'free') +
    labs(x = 'Iteration', y = 'Cost', color = NULL) +
    theme_classic(base_size = 12) +
    theme(
      panel.background = element_rect(fill = 'grey50', color = NA),
      plot.margin = margin()
    )
ggsave(
  file =
  paste0('figs/summary/optTrace', cntr, '.png'),
  plot = p4,
  device = 'png',
  type = 'cairo',
  width = 6,
  height = 4
)

# Cut up segments and visualize heat flow
# within sectors
tibble(
  seg.names = seg.names,
  buf.dir = c('r', 'l', 'r', 'r', 'l', 'l', 'r', 'r', 'r', 'r', 'r', 'r', 'r'),
  seg.num = c(8, 8, 8, 8, 8, 8, 6, 8, 8, 4, 8, 8, 8),
  scale.bar.width = c(5, 3, 4, 4, 3.5, 3, 4, 4, 2.5, 4, 5.5, 3, 2.5),
  scale.bar.height = c(0.05, 0.02, 0.05, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02),
  scale.bar.text = c(3.5, 4, 4, 4, 3.5, 4, 4, 4, 3.5, 3.5, 3.5, 4, 3.5),
  scale.bar.just = c(0.06, 0.02, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.06, 0.05, 0.05),
  scale.bar.position = c('topleft', 'bottomleft', 'bottomleft', 'topleft', 'topleft', 'bottomleft', 'bottomleft', 'bottomleft', 'bottomleft', 'topleft', 'bottomleft', 'bottomleft', 'bottomleft'),
  north.scale = c(0.4, 0.25, 0.25, 0.2, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.3, 0.2, 0.25),
  north.position = c('bottomright', 'topleft', 'topright', 'bottomright', 'bottomright', 'topleft', 'topright', 'topright', 'topright', 'bottomleft', 'topright', 'topleft', 'topright'),
  label.size = c(6, 6, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5)
) %>%
pwalk(~{
  suppressMessages({suppressWarnings({
    split_segment(
      seg.name = ..1,
      buf.dir = ..2,
      seg.num = ..3
    ) %>%
    plot_split_segment(
      running.avg = 5,
      scale.bar.width = ..4,
      scale.bar.height = ..5,
      scale.bar.text = ..6,
      scale.bar.just = ..7,
      scale.bar.position = ..8,
      north.scale = ..9,
      north.position = ..10,
      label.size = ..11
    ) -> p
    ggsave(
      file = paste0('figs/upperPlate/', str_replace_all(..1, ' ', ''), 'UpperPlate.png'),
      plot = p,
      device = 'png',
      type = 'cairo',
      width = 6,
      height = 6
    )
  })})
})

cat('\nDone!\n')