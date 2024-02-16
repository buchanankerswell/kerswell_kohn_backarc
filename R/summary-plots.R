#!/usr/bin/env Rscript

# Load functions and libraries
cat(rep('~', 60), '\n', sep='')
cat('Loading packages and functions ...\n\n')

source('R/functions.R')
load('data/hf.Rdata')
load('data/opt.RData')

# Create directory
dir.create('figs/summary', recursive=T, showWarnings=F)
dir.create('figs/variograms', recursive=T, showWarnings=F)
dir.create('figs/upper_plate', recursive=T, showWarnings=F)

# Visualize
cat('\n', rep('~', 60), sep='')
cat('\nVisualizing ...\n')

# Heat flow summary plot
p1 <-
  shp.hf.crop %>%
  map_df(~st_set_geometry(.x, NULL), .id='segment') %>%
  group_by(segment) %>%
  ggplot() +
  ggtitle('Observed heat flow: thermoglobe') +
  geom_boxplot(
    aes(x=hf, y=reorder(segment, hf, FUN=median), fill=segment, group=segment),
    outlier.shape=NA,
    size=0.5,
    color='black',
    notch=T
  ) +
  labs(x=bquote('heat flow'~(mWm^-2)), y='segment', fill='segment') +
  scale_fill_discrete_qualitative(palette='set 2') +
  theme_dark(base_size=14) +
  theme(
    legend.position='none',
    plot.margin=margin(1, 1, 1, 1),
    panel.grid=element_blank()
  )
# Save
cat('\nSaving plot to: figs/summary/hfSummary.png')
suppressWarnings(suppressMessages(
  ggsave(
    file='figs/summary/hfSummary.png',
    plot=p1,
    device='png',
    type='cairo',
    width=6.5,
    height=4.34,
    dpi=330
  )
))

# Plot all variogram models
cat('\nSaving plots to: figs/variograms/')

# Draw composite plots
unique(solns$segment) %>%
walk(~{
  ylim <-
    c(0,
      round(
        max(map_dbl(map(solns$opt.exp.vgrm[solns$segment == .x], ~range(.x$gamma)), ~.x[2])),
        digits=-1
      )
    )
  xlim <-
    c(0,
      round(
        max(map_dbl(map(solns$opt.exp.vgrm[solns$segment == .x], ~range(.x$dist)), ~.x[2])),
        digits=-1
      )
    )
  p <-
    filter(solns, segment == .x) %>%
    pmap(~{
      plot_vgrm(
        experimental.vgrm=..4,
        fitted.vgrm=..5,
        cost=..7,
        v.mod=..2,
        lineCol='white',
        ylim=sqrt(ylim),
        xlim=xlim/1e3
      )
    })
  pp <-
    (
      p[[1]] +
      theme(
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank()
      )
    ) +
    (
      p[[2]] +
      theme(
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank()
      )
    ) +
    (
      p[[3]]
    ) +
    (
      p[[4]] +
      theme(
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank()
      )
    ) +
    plot_layout(guides='collect') &
    theme(plot.margin=margin(1, 12, 1, 1), legend.position='top')
  ggsave(
    file=paste0('figs/variograms/', str_replace_all(.x, ' ', ''), '-variograms.png'),
    plot=pp,
    device='png',
    type='cairo',
    width=6.5,
    height=5.2,
    dpi=330
  )
})

# Summarise variogram models
param.labels <-
  list(
    'lag cutoff'='lag cutoff',
    'number of lags'='number of lags',
    'lag shift'='lag shift',
    'max local pairs'='max local pairs',
    'range'='range (km)',
    'sill'=bquote(sqrt('sill')~(mWm^-2))
  )
p2 <-
  vgrm.summary %>%
  mutate(
    n.lags=n.lags,
    n.max=n.max,
    sill=sqrt(sill),
    range=range/1e3
  ) %>%
  drop_na() %>%
  select(-c(rmse, dir.rmse)) %>%
  filter(cost < (median(cost)+(5*IQR(cost))) & range < (median(range)+(5*IQR(range)))) %>%
  rename(
    'lag cutoff'=cutoff.prop,
    'number of lags'=n.lags,
    'lag shift'=lag.start,
    'max local pairs'=n.max,
    'range'=range,
    'sill'=sill
  ) %>%
  select(-c(sill, range)) %>%
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
  ggtitle('Kriging parameters selected by iterative optimization') +
  geom_point(
    aes(y=cost, x=value, fill=segment, group=name),
    shape=22,
    stroke=0.5,
    size=2,
    key_glyph='rect'
  ) +
  facet_wrap(
    ~name,
    ncol=2,
    scales='free_x',
    labeller=function(variable, value) {
      return(param.labels[value])
    }
  ) +
  labs(x=NULL, y='cost', color='segment') +
  scale_color_discrete_qualitative('set 2') +
  scale_fill_discrete_qualitative('set 2') +
  theme_dark(base_size=14) +
  theme(
    panel.grid=element_blank(),
    strip.background=element_rect(color=NA, fill=NA),
    strip.text=element_text(color='black', size=14, margin=margin(0, 0, 2, 0)),
    strip.placement='outside',
    legend.box.margin=margin(0, 0, 0, 0),
    legend.key=element_rect(color='black'),
    plot.margin=margin(1, 1, 1, 1)
  )
# Save
cat('\nSaving plot to: figs/summary/variogram-summary.png', sep='')
suppressWarnings(suppressMessages(
  ggsave(
    file=paste0('figs/summary/variogram-summary.png'),
    plot=p2,
    device='png',
    type='cairo',
    width=6.5,
    height=4.34,
    dpi=330
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
    .id='segment'
  ) %>%
  ggplot() +
  ggtitle('Kriging vs. similarity interpolation differences') +
  geom_boxplot(
    aes(
      x=est.diff,
      y=reorder(segment, est.diff, FUN=median),
      fill=segment,
      group=segment
    ),
    outlier.shape=NA,
    size=0.5,
    color='black',
    notch=T
   ) +
  labs(x=bquote('heat flow'~(mWm^-2)), y='segment', fill='segment') +
  scale_fill_discrete_qualitative(palette='set 2') +
  coord_cartesian(xlim=c(-50, 50)) +
  theme_dark(base_size=14) +
  theme(
    legend.position='none',
    plot.margin=margin(1, 1, 1, 1),
    panel.grid=element_blank()
  )
# Save plot
cat('\nSaving plot to: figs/summary/interpDiffSummary.png', sep='')
suppressWarnings(suppressMessages(
  ggsave(
    file=paste0('figs/summary/interpDiffSummary.png'),
    plot=p3,
    device='png',
    type='cairo',
    width=6.5,
    height=4.34,
    dpi=330
  )
))

p3b <-
  dif$shp.interp.diff %>%
  set_names(dif$segment) %>%
  map_df(
    ~st_set_geometry(.x, NULL),
    .id='segment'
  ) %>%
  ggplot() +
  ggtitle('Kriging vs. similarity interpolation differences') +
  geom_boxplot(
    aes(
      x=sigma.diff,
      y=reorder(segment, sigma.diff, FUN=median),
      fill=segment,
      group=segment
    ),
    outlier.shape=NA,
    size=0.5,
    color='black',
    notch=T
  ) +
  labs(x=bquote(sigma~(mWm^-2)), y='segment', fill='segment') +
  scale_fill_discrete_qualitative(palette='set 2') +
  coord_cartesian(xlim=c(-50, 50)) +
  theme_dark(base_size=14) +
  theme(
    legend.position='none',
    plot.margin=margin(1, 1, 1, 1),
    panel.grid=element_blank()
  )
# Save plot
cat('\nSaving plot to: figs/summary/interpSigmaDiffSummary.png', sep='')
suppressWarnings(suppressMessages(
  ggsave(
    file=paste0('figs/summary/interpSigmaDiffSummary.png'),
    plot=p3b,
    device='png',
    type='cairo',
    width=6.5,
    height=4.34,
    dpi=330
  )
))

# nloptr trace
p4 <-
  opt.trace %>%
  group_by(segment, v.mod) %>%
  ggplot() +
  ggtitle('Cost evaluation during iterative optimization') +
  geom_path(aes(itr, cost, group=segment), linewidth=0.2) +
  geom_point(
    aes(itr, cost, group=segment, fill=segment),
    shape=22,
    stroke=0.5,
    key_glyph='rect'
  ) +
  facet_wrap(~v.mod, nrow=2) +
  labs(x='iteration', y='cost', color='segment') +
  scale_fill_discrete_qualitative(palette='set 2') +
  scale_x_continuous(breaks=seq(0, max(opt.trace$itr), 3)) +
  theme_dark(base_size=14) +
  theme(
    strip.background=element_rect(color=NA, fill=NA),
    strip.text=element_text(color='black', size=14, margin=margin()),
    strip.placement='outside',
    legend.box.margin=margin(0, 0, 0, 0),
    legend.key=element_rect(color='black'),
    plot.margin=margin(1, 1, 1, 1),
    panel.grid=element_blank()
  )
ggsave(
  file='figs/summary/optTrace.png',
  plot=p4,
  device='png',
  type='cairo',
  width=6.5,
  height=4.34,
  dpi=330
)

if(!file.exists('data/sectors.RData')){
  # Split segments into sectors
  cat('\n\nSplitting segments into sectors ...')
  shp.sectors <-
    tibble(
      seg.names=seg.names,
      buf.dir=c('r', 'l', 'r', 'r', 'l', 'l', 'r', 'r', 'r', 'r', 'r', 'r', 'r'),
      seg.num=c(14, 14, 8, 14, 8, 8, 6, 8, 8, 8, 14, 12, 8),
      sector.exclude =
        list(
          c(3,4,10), c(1,2,3,14), 3, c(1,2,14), NULL, c(2,4), 5,
          c(1,2,7), c(1,3), c(1,4,5,6,7,8), c(2,13,14), 12, 8
        )
    ) %>%
    pmap(~suppressMessages(suppressWarnings(split_segment(..1, shp_hf_crop, shp_segs,
                                                          shp_volc, ..2, ..3, ..4)))) %>%
    set_names(seg.names)
  # Save
  save(shp.sectors, file=paste0('data/sectors.RData'))
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
  cat('\nSaving plot to: figs/upper_plate/', str_replace_all(names(shp.sectors)[.x], ' ', ''),
      'UpperPlate.png', sep='')
  suppressMessages({suppressWarnings({
    plot_split_segment(shp.sectors[[.x]], 3, .y) -> p
    ggsave(file = paste0('figs/upper_plate/',
                         str_replace_all(names(shp.sectors)[.x], ' ', ''), 'UpperPlate.png'),
           plot=p, device='png', type='cairo', width=6.5, height=6.5, dpi=330)
  })})
})

cat('\nsummary-plots.R complete!\n\n')