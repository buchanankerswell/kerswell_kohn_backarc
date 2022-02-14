#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

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

source('R/functions.R')
load('data/hf.Rdata')
load(paste0('data/opt', cntr, '.RData'))
dir.create('figs/summary', showWarnings = F)
dir.create('figs/vgrms', showWarnings = F)
dir.create('figs/upperPlate', showWarnings = F)
cat('\nSaving plots to: figs/summary/')
cat('\nSaving plots to: figs/vgrms/\n')

# Visualize
cat('\n', rep('~', 60), sep='')
cat('\nVisualizing ...\n')

# Heat flow summary plot
p1 <-
  shp.hf.crop %>%
  map_df(~st_set_geometry(.x, NULL), .id = 'segment') %>%
  group_by(segment) %>%
  ggplot() +
  geom_boxplot(
    aes(x=hf, y=segment, group=segment),
    fill = 'ivory',
    outlier.shape = NA,
    size = 0.5,
    color = 'black'
  ) +
  labs(x = 'milliwatts per square meter', y = NULL) +
  scale_y_discrete(limits = rev(levels(as.factor(seg.names)))) +
  theme_dark(base_size = 12) +
  theme(legend.position = 'none', plot.margin = margin(1, 1, 1, 1))
# Save
cat('\nSaving plot to: figs/summary/hfSummary.png')
suppressWarnings(suppressMessages(
  ggsave(
    file = 'figs/summary/hfSummary.png',
    plot = p1,
    device = 'png',
    type = 'cairo',
    width = 6.5,
    height = 4.34
  )
))

# Plot all variogram models
cat('\nSaving plots to: figs/vgrms/')

# Draw composite plots
unique(solns$segment) %>%
walk(~{
  ylim <-
    c(0,
      round(
        max(map_dbl(map(solns$opt.exp.vgrm[solns$segment == .x], ~range(.x$gamma)), ~.x[2])),
        digits = -1
      )
    )
  xlim <-
    c(0,
      round(
        max(map_dbl(map(solns$opt.exp.vgrm[solns$segment == .x], ~range(.x$dist)), ~.x[2])),
        digits = -1
      )
    )
  p <-
    filter(solns, segment == .x) %>%
    pmap(~{
      plot_vgrm(
        experimental.vgrm = ..4,
        fitted.vgrm = ..5,
        cost = ..7,
        v.mod = ..2,
        lineCol = 'white',
        ylim = ylim,
        xlim = xlim/1e5
      )
    })
  pp <-
    (
      p[[1]] +
      theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
      )
    ) +
    (
      p[[2]] +
      theme(
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()
      )
    ) +
    (p[[3]]) +
    (
      p[[4]] +
      theme(
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()
      )
    ) &
    theme(plot.margin = margin(1, 1, 1, 1))
  ggsave(
    paste0('figs/vgrms/', str_replace_all(.x, ' ', ''), 'Vgrms.png'),
    plot = pp,
    device = 'png',
    type = 'cairo',
    width = 6,
    height = 3
  )
})

# Summarise variogram models
init.vals <-
  tibble(
    'lag cutoff' = 3,
    'number of lags' = 20,
    'lag shift' = 1,
    'max local pairs' = 8,
    'log10 range' = NA,
    'log10 sill' = NA
  ) %>%
  pivot_longer(everything())
p2 <-
  vgrm.summary %>%
  mutate(
    n.lags = n.lags,
    n.max = n.max,
    sill = log10(sill),
    range = log10(range/1000)
  ) %>%
  drop_na() %>%
  select(-c(rmse, dir.rmse)) %>%
  filter(cost < (median(cost)+(5*IQR(cost))) & range < (median(range)+(5*IQR(range)))) %>%
  rename(
    'lag cutoff' = cutoff.prop,
    'number of lags' = n.lags,
    'lag shift' = lag.start,
    'max local pairs' = n.max,
    'log10 range' = range,
    'log10 sill' = sill,
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
  ggplot() +
  geom_point(aes(x = cost*100, y = value, shape = v.mod, group = name), size = 2) +
  facet_wrap(~name, ncol = 2, scales = 'free_y') +
  scale_color_discrete_qualitative('Dark 3') +
  scale_shape_manual(values = 15:18) +
  guides(
    shape =
      guide_legend(
        override.aes = list(size = 3),
        title.position = 'top',
        title.vjust = 1,
        label.position = 'bottom'
      )
  ) +
  labs(y = NULL, x = 'cost x 100', shape = NULL) +
  theme_dark(base_size = 12) +
  theme(
    strip.background = element_rect(color = 'grey50', fill = 'grey95'),
    strip.text = element_text(color = 'black', size = 12, margin = margin(1, 0, 1.5, 0)),
    legend.box.margin = margin(-10, 0, 0, 0),
    plot.margin = margin(1, 1, 1, 1)
  )
p2 <-
  p2 +
  geom_text(
    data = init.vals,
    aes(x = -Inf, y = value, group = name, label = 'inital value'),
    color = 'white',
    size = 3,
    hjust = -0.1,
    vjust = -0.1
  ) +
  geom_hline(
    data = init.vals,
    aes(yintercept = value, group = name),
    color = 'grey90',
    size = 0.8
  ) +
  geom_smooth(
    aes(x = cost*100, y = value, group = name),
    formula = y~x,
    method = 'lm',
    size = 0.5,
    color = 'black',
    fill = 'ivory',
    fullrange = T
  )
# Save
cat('\nSaving plot to: figs/summary/vgrmSummary.png', sep = '')
suppressWarnings(suppressMessages(
  ggsave(
    file = paste0('figs/summary/vgrmSummary.png'),
    plot = p2,
    device = 'png',
    type = 'cairo',
    width = 6.5,
    height = 6.5
  )
))

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
  geom_boxplot(
    aes(x=est.diff, y=segment, group=segment),
    fill = 'ivory',
    outlier.shape = NA,
    size = 0.5,
    color = 'black'
   ) +
  labs(x = 'milliwatts per square meter', y = NULL) +
  coord_cartesian(xlim = c(-50, 50)) +
  scale_y_discrete(limits = rev(levels(as.factor(seg.names)))) +
  theme_dark(base_size = 12) +
  theme(legend.position = 'none', plot.margin = margin(1, 1, 1, 1))
# Save plot
cat('\nSaving plot to: figs/summary/interpDiffSummary.png', sep = '')
suppressWarnings(suppressMessages(
  ggsave(
    file = paste0('figs/summary/interpDiffSummary.png'),
    plot = p3,
    device = 'png',
    type = 'cairo',
    width = 6.5,
    height = 4.34
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
  geom_boxplot(
    aes(x=sigma.diff, y=segment, group=segment),
    fill = 'ivory',
    outlier.shape = NA,
    size = 0.5,
    color = 'black'
  ) +
  labs(x = 'milliwatts per square meter', y = NULL) +
  coord_cartesian(xlim = c(-50, 50)) +
  scale_y_discrete(limits = rev(levels(as.factor(seg.names)))) +
  theme_dark(base_size = 12) +
  theme(legend.position = 'none', plot.margin = margin(1, 1, 1, 1))
# Save plot
cat('\nSaving plot to: figs/summary/interpSigmaDiffSummary.png', sep = '')
suppressWarnings(suppressMessages(
  ggsave(
    file = paste0('figs/summary/interpSigmaDiffSummary.png'),
    plot = p3b,
    device = 'png',
    type = 'cairo',
    width = 6,
    height = 4
  )
))

# nloptr trace
p4 <-
  opt.trace %>%
  group_by(segment, v.mod) %>%
  ggplot() +
  geom_path(aes(itr, cost*100, group = segment)) +
  facet_wrap(~v.mod, nrow = 2) +
  labs(x = 'iteration', y = 'cost x 100', color = NULL) +
  theme_dark(base_size = 12) +
  theme(
    strip.background = element_rect(color = 'grey50', fill = 'grey95'),
    strip.text = element_text(color = 'black', size = 12, margin = margin(1, 0, 1.5, 0)),
    legend.box.margin = margin(-10, 0, 0, 0),
    plot.margin = margin(1, 1, 1, 1)
  )
ggsave(
  file =
  paste0('figs/summary/optTrace.png'),
  plot = p4,
  device = 'png',
  type = 'cairo',
  width = 6,
  height = 4.67
)

if(!file.exists('data/sectors.RData')){
  # Split segments into sectors
  cat('\n\nSplitting segments into sectors ...')
  shp.sectors <-
    tibble(
      seg.names = seg.names,
      buf.dir = c('r', 'l', 'r', 'r', 'l', 'l', 'r', 'r', 'r', 'r', 'r', 'r', 'r'),
      seg.num = c(14, 14, 8, 14, 8, 8, 6, 8, 8, 8, 14, 12, 8),
      sector.exclude =
        list(
          c(3,4,10), c(1,2,3,14), 3, c(1,2,14), NULL, c(2,4), 5,
          c(1,2,7), c(1,3), c(1,4,5,6,7,8), c(2,13,14), 12, 8
        )
    ) %>%
    pmap(~{
      suppressMessages({suppressWarnings({
        split_segment(
          seg.name = ..1,
          buf.dir = ..2,
          seg.num = ..3,
          sector.exclude = ..4
        )
      })})
    }) %>%
    set_names(seg.names)
  # Save
  save(shp.sectors, file = paste0('data/sectors.RData'))
} else {
  load('data/sectors.RData')
}

borders <-
  list(
    # left right top bottom
    c(-0.12, -0.18, 0.5, 0.1), # Alaska Aleutians
    c(0.2, 1.2, -0.13, -0.225), # Andes
    c(-0.2, -0.2, -0.025, -0.25), # Central America
    c(0.45, -0.2, -0.1, -0.17), # Kamchatka Marianas
    c(-0.1, -0.2, -0.05, -0.25), # Kyushu Ryukyu
    c(0.15, -0.28, -0.25, -0.2), # Lesser Antilles
    c(-0.3, 0.1, -0.3, -0.3), # N Philippines
    c(-0.3, -0.15, 0, -0.4), # New Britain Solomon
    c(0.8, -0.2, -0.175, -0.125), # S Philippines
    c(-0.25, -0.25, -0.25, -0.25), # Scotia
    c(-0.1, -0.3, -0.13, -0.15), # Sumatra Banda Sea
    c(0.5, -0.1, -0.125, -0.15), # Tonga New Zealand
    c(-0.2, 0.3, -0.05, -0.25) # Vanuatu
  )

walk2(1:13, borders, ~{
  cat(
    '\nSaving plot to: figs/upperPlate/',
    str_replace_all(names(shp.sectors)[.x], ' ', ''),
    'UpperPlate.png',
    sep = ''
  )
  suppressMessages({suppressWarnings({
    plot_split_segment(
      split.seg = shp.sectors[[.x]],
      running.avg = 3,
      borders = .y
    ) -> p
    ggsave(
      file =
        paste0(
          'figs/upperPlate/',
          str_replace_all(names(shp.sectors)[.x], ' ', ''),
          'UpperPlate.png'
        ),
      plot = p,
      device = 'png',
      type = 'cairo',
      width = 6,
      height = 5.5
    )
  })})
})

cat('\n\nDone!\n\n')