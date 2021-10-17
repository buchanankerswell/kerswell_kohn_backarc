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

# Bounding Boxes
purrr::map(seg.names,
  ~shp.sa.segs.robin.pacific.buffer %>%
  filter(segment == .x) %>%
  st_bbox() %>%
  bbox_widen(
    crs = proj4.robin.pacific,
    borders = c('top' = 0.1,
    'bottom' = 0.1,
    'left' = 0.1,
    'right' = 0.1))) %>%
purrr::set_names(nm = seg.names) -> shp.box

# Crop data
purrr::map2_df(shp.box, seg.names,
  ~shp.hf %>%
  rename(hf = `heat-flow (mW/m2)`) %>%
  st_crop(.x) %>%
  mutate(segment = .y, .before = country)) -> shp.hf.crop

# Summarize heat flow data
shp.hf.crop %>%
st_set_geometry(NULL) %>%
group_by(segment) %>%
mutate('segment' = segment %>% stringr::str_replace_all('_', ' ')) %>%
rename(Segment = segment) %>%
summarise(
  n = n(),
  Min = round(min(hf)),
  Max = round(max(hf)),
  Median =round(median(hf)),
  IQR = round(IQR(hf)),
  Mean = round(mean(hf)),
  Sigma = round(sd(hf))) -> hf.summary

cat('Heat flow summary:\n')
print(hf.summary)

# Visualize

shp.hf.crop %>%
mutate('segment' = segment %>% stringr::str_replace_all('_', ' ')) %>%
group_by(segment) %>%
ggplot() +
geom_boxplot(
  aes(x = hf, y = segment, group = segment),
  width = 0.5,
  outlier.size = 0.5,
  outlier.color = rgb(0.5, 0.5, 0.5, 0.1)) +
labs(
  x = bquote('Heat flow'~(mWm^-2)),
  y = NULL,
  title = 'Heat flow observations') +
scale_x_continuous(limits = c(0, 250)) +
scale_y_discrete(
  limits = rev(levels(as.factor(seg.names %>% stringr::str_replace_all('_', ' '))))) +
theme_classic(base_size = 9) +
theme(strip.background = element_blank()) -> p

# Save plot
cat('Saving heat flow summary plot to:\nfigs/summary/hf_summary.png\n')

library(ggridges)
shp.hf.crop %>%
mutate('segment' = segment %>% stringr::str_replace_all('_', ' ')) %>%
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
  limits = rev(levels(as.factor(seg.names %>% stringr::str_replace_all('_', ' '))))) +
theme_classic() +
theme(
  legend.position='none',
  plot.margin=margin()
) -> p

ggsave(
  file = 'figs/summary/hf_summary.png',
  plot = p,
  device = 'png',
  type = 'cairo',
  width = 6,
  height = 3.5)

# Summarise variogram models (add rmse's from f.obj)
purrr::map_dbl(shp.box,
  ~{box <- .x %>% st_bbox
    as.numeric((box$xmax-box$xmin)*(box$ymax-box$ymin))
}) -> shp.area

purrr::map_df(fname,
  ~get(.x)$v.mod, .id = 'segment') %>%
  as_tibble() %>%
  mutate('segment' = seg.names %>% stringr::str_replace_all('_', ' ')) %>%
  dplyr::select(segment, model, psill, range) %>%
  rename(
    Segment = segment, Model = model,
    Sill = psill, Range = range) %>%
  mutate(
    'Sill' = round(sqrt(Sill)), # mW/m^-2
    'Range' = round(Range/1000, 1), # km
    'RMSE' = purrr::map_dbl(fname, ~round(attr(get(.x)$k, 'cv.rmse'), 1)),
    'Length' = (shp.sa.segs.robin.pacific %>% st_length() %>% as.vector())/1000, # km
    'Area' = shp.area/10^12, # Mkm^2
    n = hf.summary$n,
    'Density' = n/Area
  ) -> variogram.summary

cat('Variogram summary:\n')
print(variogram.summary)

# Visualize
variogram.summary %>%
ggplot() +
geom_bar(aes(x = Range, y = Segment), stat = 'identity') +
annotate(
  'segment',
  y = 'Kyushu Ryukyu',
  yend = 'Kyushu Ryukyu',
  x = 400,
  xend = 450,
  arrow = arrow(length = unit(2, 'mm')),
  color = 'white') +
annotate(
  'segment',
  y = 'Vanuatu',
  yend = 'Vanuatu',
  x = 400,
  xend = 450,
  arrow = arrow(length = unit(2, 'mm')),
  color = 'white') +
coord_cartesian(xlim = c(0, 450)) +
labs(x = 'Range (km)', y = NULL) +
scale_y_discrete(
  limits = rev(levels(as.factor(seg.names %>% stringr::str_replace_all('_', ' '))))) +
theme_classic(base_size = 9) -> p1

variogram.summary %>%
ggplot(aes(x = Range, y = Length, label = Segment)) +
geom_point() +
geom_text_repel(size = 3, force = 2) +
labs(x = 'Range (km)', y = 'Segment Length (km)') +
theme_classic() -> p2

variogram.summary %>%
ggplot(aes(x = Range, y = n, label = Segment)) +
geom_point() +
geom_text_repel(size = 3, force = 2) +
labs(x = 'Range (km)', y = 'Number of observations') +
theme_classic() -> p3

variogram.summary %>%
mutate(area = shp.area) %>%
ggplot(aes(x = Range, y = Area, label = Segment)) +
  geom_point() +
  geom_text_repel(size = 3, force = 2) +
  labs(x = 'Range (km)', y = bquote('Domain area'~(Mkm^2))) +
  theme_classic() -> p4

variogram.summary %>%
ggplot(aes(x = Range, y = n/Area, label = Segment)) +
  geom_point() +
  geom_text_repel(size = 3, force = 2) +
  labs(x = 'Range (km)', y = bquote('Observation density'~(n/Mkm^2))) +
  theme_classic() -> p5

variogram.summary %>%
rename(
  'Area~(Mkm^2)' = Area,
  'Density~(n/Mkm^2)' = Density,
  'Length~(km)' = Length,
  'Range~(km)' = Range,
  'Sill~(mWm^-2)' = Sill
) %>%
pivot_longer(-c(Segment, RMSE, Model)) %>%
ggplot(aes(x = RMSE, y = value)) +
geom_point() +
facet_wrap(~name, scales = 'free_y', labeller = label_parsed) +
labs(y = NULL) +
theme_classic() +
guides(color=guide_legend(nrow=4)) +
theme(
  strip.background = element_rect(fill = 'grey90', color=NA)
) -> p

# Composition
(p2 + theme(axis.text.x=element_blank(), axis.title.x=element_blank())) +
(p3 + theme(axis.text.x=element_blank(), axis.title.x=element_blank())) +
p4 + p5 +
plot_annotation(tag_levels = 'a') -> p

# Save plot
cat('Saving variogram model summary plot to:\nfigs/summary/variogram_summary.png\n')

ggsave(
  file = 'figs/summary/variogram_summary.png',
  plot = p,
  device = 'png',
  type = 'cairo',
  width = 6,
  height = 3.5)

# Interpolation difference
purrr::map(fname, ~get(.x)$diff %>% st_set_geometry(NULL)) %>%
purrr::set_names(nm = seg.names) -> hf.diff

hf.diff %>%
bind_rows(.id = 'Segment') %>%
group_by(Segment) %>%
mutate('Segment' = Segment %>% stringr::str_replace_all('_', ' ')) %>%
summarise(
  n = n(),
  Min = round(min(hf.diff)),
  Max = round(max(hf.diff)),
  Median = round(median(hf.diff)),
  IQR = round(IQR(hf.diff)),
  Mean = round(mean(hf.diff)),
  Sigma = round(sd(hf.diff))
) -> hf.diff.summary

cat('Heat flow difference summary:\n')
print(hf.diff.summary)

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