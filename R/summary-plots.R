#!/usr/bin/env Rscript

# Load functions and libraries
source('R/functions.R')
load('assets/hf_data/preprocessed-hf-data.Rdata')
load('assets/map_data/preprocessed-map-data.Rdata')
load('assets/opt_data/nloptr.RData')
load('assets/opt_data/krige.RData')

# Create directory
dir.create('figs/summary', recursive=T, showWarnings=F)
dir.create('figs/variograms', recursive=T, showWarnings=F)
dir.create('figs/upper_plate', recursive=T, showWarnings=F)

# Heat flow summary plot
p1 <-
  shp_hf_crop %>%
  map_df(~st_set_geometry(.x, NULL), .id='segment') %>%
  group_by(segment) %>%
  ggplot() +
  ggtitle('Measured Heat Flow') +
  geom_boxplot(aes(x=hf, y=reorder(segment, hf, FUN=median), fill=segment, group=segment),
               size=0.5, color='black', notch=T) +
  labs(x=bquote('Heat Flow'~(mWm^-2)), y='Segment') +
  scale_fill_discrete_qualitative(palette='set 2') +
  theme_bw(base_size=14) +
  theme(legend.position='none', plot.margin=margin(5, 5, 5, 5),
        panel.border=element_rect(linewidth=1.2),
        panel.background=element_rect(fill='grey90'))

cat('\nSaving plot to: figs/summary/hf-summary.png')
ggsave(file='figs/summary/hf-summary.png', plot=p1, width=6.5, height=4.34)

# Draw composite plots
unique(opt_solutions$segment) %>%
  walk(~{
    ylim <-
      c(0, round(max(map_dbl(map(opt_solutions$opt_exp_vgrm[opt_solutions$segment == .x],
                                 ~range(.x$gamma)), ~.x[2])), digits=-1))
    xlim <-
      c(0,round(max(map_dbl(map(opt_solutions$opt_exp_vgrm[opt_solutions$segment == .x],
                                ~range(.x$dist)), ~.x[2])), digits=-1))
    p <-
      filter(opt_solutions, segment == .x) %>%
      pmap(~plot_vgrm(experimental_vgrm=..4, fitted_vgrm=..5, cost=..7, v_mod=..2,
                      ylim=sqrt(ylim), xlim=xlim / 1e3))
    pp <-
      (p[[1]] + theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(),
                      axis.text.x=element_blank())) +
      (p[[2]] + theme(axis.title=element_blank(), axis.ticks=element_blank(),
                      axis.text=element_blank())) +
      (p[[3]]) +
      (p[[4]] + theme(axis.title.y=element_blank(), axis.ticks.y=element_blank(),
                      axis.text.y=element_blank())) +
      plot_layout(guides='collect') & theme(legend.position='top')
    ggsave(file=paste0('figs/variograms/', str_replace_all(.x, ' ', ''), '-variograms.png'),
           plot=pp, width=6.5, height=5.2)
})

# Summarise variogram models
param_labels <- list('lag cutoff'='Lag Cutoff', 'number of lags'='Number of Lags',
                     'lag shift'='Lag Shift', 'max local pairs'='Max Local Pairs',
                     'range'='Range (km)', 'sill'=bquote(sqrt('Sill')~(mWm^-2)))

p2 <-
  vgrm_summary %>%
  mutate(n_lags=n_lags, n_max=n_max, sill=sqrt(sill), range=range / 1e3) %>%
  drop_na() %>%
  select(-c(rmse, dir_rmse)) %>%
  rename('lag cutoff'=cutoff, 'number of lags'=n_lags, 'lag shift'=lag_start,
         'max local pairs'=n_max, 'range'=range, 'sill'=sill) %>%
  select(-c(sill, range)) %>%
  pivot_longer(-c(v_mod, cost, segment, itr, vgrm_wt, vgrm_rmse, vgrm_cost, cv_wt, cv_rmse,
                  cv_cost)) %>%
  ggplot() +
  ggtitle('Optimal Kriging Parameters') +
  geom_point(aes(y=cost, x=value, fill=segment, group=name), shape=22, stroke=0.5, size=2,
             key_glyph='rect') +
  facet_wrap(~name, ncol=2, scales='free_x',
             labeller=labeller(name=param_labels)) +
  labs(x=NULL, y='Cost', color='Segment', fill='Segment') +
  scale_color_discrete_qualitative('set 2') +
  scale_fill_discrete_qualitative('set 2') +
  theme_bw(base_size=14) +
  theme(strip.background=element_rect(color=NA, fill=NA),
        strip.text=element_text(color='black', size=14, margin=margin(0, 0, 2, 0)),
        strip.placement='outside', legend.box.margin=margin(0, 0, 0, 0),
        legend.key=element_rect(color='black'), plot.margin=margin(5, 5, 5, 5),
        panel.border=element_rect(linewidth=1.2),
        panel.background=element_rect(fill='grey90'))

cat('\nSaving plot to: figs/summary/variogram-summary.png', sep='')
ggsave(file=paste0('figs/summary/variogram-summary.png'), plot=p2, width=6.5, height=4.34)

# Interpolation differences
dif <- opt_solutions %>% group_by(segment) %>% slice_min(cost)
p3 <-
  dif$shp_interp_diff %>% set_names(dif$segment) %>%
  map_df(~st_set_geometry(.x, NULL), .id='segment') %>%
  ggplot() +
  ggtitle('Kriging–Similarity Residuals') +
  geom_boxplot(aes(x=est_diff, y=reorder(segment, est_diff, FUN=median), fill=segment,
                   group=segment), size=0.5, color='black', notch=T) +
  labs(x=bquote('Heat Flow'~(mWm^-2)), y='Segment', fill='Segment') +
  scale_fill_discrete_qualitative(palette='set 2') +
  theme_bw(base_size=14) +
  theme(legend.position='none', plot.margin=margin(5, 5, 5, 5),
        panel.border=element_rect(linewidth=1.2),
        panel.background=element_rect(fill='grey90'))

cat('\nSaving plot to: figs/summary/interp-diff-summary.png', sep='')
ggsave(file=paste0('figs/summary/interp-diff-summary.png'), plot=p3, width=6.5, height=4.34)

p3b <-
  dif$shp_interp_diff %>% set_names(dif$segment) %>%
  map_df(~st_set_geometry(.x, NULL), .id='segment') %>%
  ggplot() +
  ggtitle('Kriging–Similarity Uncertainties') +
  geom_boxplot(aes(x=sigma_diff, y=reorder(segment, sigma_diff, FUN=median), fill=segment,
                   group=segment), size=0.5, color='black', notch=T) +
  labs(x=bquote('Uncertainty'~2~sigma(mWm^-2)), y='Segment', fill='Segment') +
  scale_fill_discrete_qualitative(palette='set 2') +
  theme_bw(base_size=14) +
  theme(legend.position='none', plot.margin=margin(5, 5, 5, 5),
        panel.border=element_rect(linewidth=1.2),
        panel.background=element_rect(fill='grey90'))

cat('\nSaving plot to: figs/summary/interp-sigma-diff-summary.png', sep='')
ggsave(file=paste0('figs/summary/interp-sigma-diff-summary.png'), plot=p3b, width=6.5,
       height=4.34,)

# nloptr trace
p4 <-
  nloptr_trace %>% group_by(segment, v_mod) %>%
  ggplot() +
  ggtitle('NLopt Minimization') +
  geom_path(aes(itr, cost, group=segment), linewidth=0.2) +
  geom_point(aes(itr, cost, group=segment, fill=segment), shape=22, stroke=0.5,
             key_glyph='rect') +
  facet_wrap(~v_mod, nrow=2) +
  labs(x='Iteration', y='Cost', color='Segment', fill='Segment') +
  scale_fill_discrete_qualitative(palette='set 2') +
  scale_x_continuous(breaks=c(0, 25, 50, 100)) +
  theme_bw(base_size=14) +
  theme(strip.background=element_rect(color=NA, fill=NA),
        strip.text=element_text(color='black', size=14, margin=margin()),
        strip.placement='outside', legend.box.margin=margin(0, 0, 0, 0),
        legend.key=element_rect(color='black'), plot.margin=margin(5, 5, 5, 5),
        panel.border=element_rect(linewidth=1.2),
        panel.background=element_rect(fill='grey90'))
ggsave(file='figs/summary/opt-trace.png', plot=p4, width=6.5, height=4.34)

if(!file.exists('assets/hf_data/sectors.RData')){
  suppressWarnings({
    # Split segments into sectors
    cat('\n\nSplitting segments into sectors ...')
    shp_sectors <-
      tibble(seg_names=seg_names,
             buf_dir=c('r', 'l', 'r', 'r', 'l', 'l', 'r', 'r', 'r', 'r', 'r', 'r', 'r'),
             seg_num=c(14, 14, 8, 14, 8, 8, 6, 8, 8, 8, 14, 12, 8),
             sector_exclude = list(c(3, 4, 10), c(1, 2, 3, 14), 3, c(1, 2, 14), NULL, c(2, 4),
                                   5, c(1, 2, 7), c(1, 3), c(1, 4, 5, 6, 7, 8), c(2, 13, 14),
                                   12, 8)
             ) %>%
      pmap(~split_segment(..1, shp_hf_crop, shp_segs, shp_volc, ..2, ..3, ..4)) %>%
      set_names(seg_names)
    save(shp_sectors, file=paste0('assets/hf_data/sectors.RData'))
  })
} else {
  load('assets/hf_data/sectors.RData')
}

borders <- list(c(-0.12, -0.18, 0.5, 0.1), # Alaska Aleutians
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
                c(-0.2, 0.3, -0.05, -0.25)) # Vanuatu

walk2(1:13, borders, ~{
  suppressWarnings({
    short_lab <- str_replace_all(names(shp_sectors)[.x], ' ', '-')
    cat('\nSaving plot to: figs/upper_plate/', short_lab, '-upper-plate.png', sep='')
    plot_split_segment(shp_sectors[[.x]], .y) -> p
    ggsave(file = paste0('figs/upper_plate/', short_lab, '-upper-plate.png'), plot=p,
           width=6.5, height=6.5)
  })
})

cat('\nsummary-plots.R complete!\n\n')