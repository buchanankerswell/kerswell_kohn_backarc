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
if(is.null(cntr)) {
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
  system(
    paste0(
      'open figs/vgrms/',
      str_replace_all(.x, ' ', ''),
      'Vgrms',
      cntr,
      '.png'
    ),
    wait = F
  )
})

# Summarise variogram models
p2 <-
  vgrm.summary %>%
  filter(sill < 10000 & range < 2000000) %>%
  mutate('range' = range/1000) %>%
  rename(
    'lag cutoff' = cutoff.prop,
    'number of lags' = n.lags,
    'lag shift' = lag.start,
    'max local pairs' = n.max,
    'range' = range,
    'sill' = sill,
  ) %>%
  pivot_longer(-c(v.mod, cost, segment, itr)) %>%
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
cat('\nSaving plot to: figs/summary/vgrm_summary', cntr, '.png', sep = '')
suppressWarnings(suppressMessages(
  ggsave(
    file = paste0('figs/summary/vgrm_summary', cntr, '.png'),
    plot = p2,
    device = 'png',
    type = 'cairo',
    width = 6,
    height = 6
  )
))
system(paste0('open figs/summary/vgrm_summary', cntr, '.png'), wait = F)

# Interpolation differences
dif <-
  solns %>%
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
      x = bquote('Interpolation difference'~(mWm^-2)),
      y = NULL
    ) +
    scale_x_continuous(limits = c(-100, 100), breaks = seq(-100, 100, 25)) +
    scale_fill_viridis_d() +
    scale_y_discrete(limits = rev(levels(as.factor(seg.names)))) +
    theme_classic() +
    theme(
      legend.position = 'none',
      plot.margin = margin()
    )
# Save plot
cat('\nSaving plot to: figs/summary/interp_diff_summary', cntr, '.png', sep = '')
suppressWarnings(suppressMessages(
  ggsave(
    file = paste0('figs/summary/interp_diff_summary', cntr, '.png'),
    plot = p3,
    device = 'png',
    type = 'cairo',
    width = 6,
    height = 3.5
  )
))
system(paste0('open figs/summary/interp_diff_summary', cntr, '.png'), wait = F)

cat('\nDone!\n')