# Load functions and libraries
source('functions.R')
load('data/hf.Rdata')

# Define paths and names
list.files('data/diff', pattern = '.RData', full.names = T) -> fpath
purrr::map_chr(
  list.files('data/diff', pattern = '.RData'),
  ~.x %>%
  stringr::str_replace('.RData', '')) -> fname

# Load data
for (i in fpath) load(i)

# Visualize
cat('\n', rep('~', 60), sep='')
cat('\nVisualizing ...')
cat('\n', rep('~', 60), sep='')

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
      size = 0.1,
      geom='density_ridges_gradient',
      calc_ecdf=T,
      quantiles=4,
      quantile_lines=T
    ) +
    labs(
      x = bquote('Heat flow'~(mWm^-2)),
      y = NULL
    ) +
    scale_x_continuous(limits = c(0, 250)) +
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
ggsave(
  file = 'figs/summary/hf_summary.png',
  plot = p1,
  device = 'png',
  type = 'cairo',
  width = 6,
  height = 3.5
)

# Summarise variogram models
p2 <-
  variogram.summary %>%
  select(-c(seg.length, domain.area, model.vgrm)) %>%
  rename(
    'n[CV]' = n.cv,
    'Density[CV]~(n/Mkm^2)' = cv.density,
    'c' = cutoff.prop,
    'n[h]' = n.lags,
    'local~n[max]' = n.max,
    'Range~(km)' = range,
    'Sill~(mWm^-2)' = sill
    'RMSE[CV]~(mWm^-2)' = cv.rmse,
    'Cost~(mWm^-2)' = cost
  ) %>%
  pivot_longer(-c(segment, cv.rmse)) %>%
  ggplot(aes(x = cv.rmse, y = value)) +
    geom_point() +
    facet_wrap(~name, scales = 'free_y', labeller = label_parsed) +
    labs(y = NULL) +
    theme_classic() +
    theme(
      strip.background = element_rect(fill = 'grey90', color=NA)
    )

# Save
cat('\nSaving plot to: figs/summary/variogram_summary.png')
ggsave(
  file = 'figs/summary/variogram_summary.png',
  plot = p2,
  device = 'png',
  type = 'cairo',
  width = 6,
  height = 6
)

# Interpolation difference
# Visualize
d <- hf.diff %>%
bind_rows(.id = 'segment') %>%
mutate('segment' = segment %>% stringr::str_replace_all('_', ' ')) %>%
group_by(segment)

d %>%
ggplot() +
geom_boxplot(
  aes(x = hf.diff, y = reorder(segment, -hf.diff, FUN = median), group = segment),
  width = 0.8,
  notch = T,
  outlier.shape = NA) +
geom_vline(xintercept = 0, size = 0.2) +
labs(
  x = bquote('Heat flow difference'~(mWm^-2)),
  y = NULL) +
scale_x_continuous(breaks = seq(-20, 40, 10), limits = quantile(d$hf.diff, c(0.1, 0.9))) +
theme_classic(base_size = 9) +
theme(strip.background = element_blank()) -> p

d %>%
ggplot() +
stat_density_ridges(
  aes(x=hf.diff, y=segment, group=segment, fill=factor(stat(quantile))),
  size = 0.1,
  geom='density_ridges_gradient',
  calc_ecdf=T,
  quantiles=4,
  quantile_lines=T
) +
geom_vline(xintercept = 0, size = 0.2) +
labs(
  x = bquote('Heat flow difference'~(mWm^-2)),
  y = NULL
) +
scale_x_continuous(limits = c(-80, 80), breaks = seq(-80, 80, 20)) +
scale_fill_viridis_d() +
scale_y_discrete(
  limits = rev(levels(as.factor(seg.names %>% stringr::str_replace_all('_', ' '))))) +
theme_classic() +
theme(
  legend.position='none',
  plot.margin=margin()
) -> p


# Save plot
cat('Saving heat flow difference summary plot to:\nfigs/summary/hf_diff_summary.png')

ggsave(
  file = 'figs/summary/hf_diff_summary.png',
  plot = p,
  device = 'png',
  type = 'cairo',
  width = 6,
  height = 3.5)

cat('\nDone\n')