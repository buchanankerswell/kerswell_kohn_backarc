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

ggsave(
  file = 'figs/summary/hf_summary.png',
  plot = p,
  device = 'png',
  type = 'cairo',
  width = 4,
  height = 4)

# Summarise variogram models (add rmse's from f.obj)
purrr::map_df(fname,
  ~get(.x)$v.mod, .id = 'segment') %>%
  as_tibble() %>%
  mutate('segment' = seg.names %>% stringr::str_replace_all('_', ' ')) %>%
  dplyr::select(segment, model, psill, range) %>%
  rename(
    Segment = segment, Model = model,
    Sill = psill, Range = range) %>%
  mutate(
    'Sill' = round(sqrt(Sill)),
    'Range' = round(Range/1000, 1),
    'RMSE' = purrr::map_dbl(fname, ~round(attr(get(.x)$k, 'cv.rmse'), 1))) -> variogram.summary

cat('Variogram summary:\n')
print(variogram.summary)

# Visualize
variogram.summary %>%
mutate('Segment' = Segment %>% stringr::str_replace_all('_', ' ')) %>%
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
mutate(
  seg.length = shp.sa.segs.robin.pacific %>%
  st_length() %>%
  as.vector()) %>%
ggplot() +
geom_point(aes(x = Range, y = seg.length/1000)) +
geom_text_repel(
  aes(x = Range, y = seg.length/1000, label = Segment),
  size = 2,
  color = rgb(0, 0, 0, 0.3)) +
labs(x = 'Range (km)', y = 'Segment Length (km)') +
theme_classic(base_size = 9) -> p2

variogram.summary %>%
mutate(n = hf.summary$n) %>%
ggplot() +
geom_point(aes(x = Range, y = n)) +
geom_text_repel(
  aes(x = Range, y = n, label = Segment),
  size = 2,
  color = rgb(0, 0, 0, 0.3)) +
labs(x = 'Range (km)', y = 'Number of observations') +
theme_classic(base_size = 9) -> p3

purrr::map_dbl(
  shp.box,
  ~{box <- .x %>% st_bbox
    as.numeric((box$xmax-box$xmin)*(box$ymax-box$ymin))}) -> shp.area

variogram.summary %>%
mutate(area = shp.area) %>%
ggplot() +
  geom_point(aes(x = Range, y = area/10^12)) +
  geom_text_repel(
    aes(x = Range, y = area/10^12, label = Segment),
    size = 2,
    color = rgb(0, 0, 0, 0.3)) +
  labs(x = 'Range (km)', y = bquote('Domain area'~(km^2%*%10^9))) +
  theme_classic(base_size = 9) -> p4

# Composition
p1 + p2 +
(p4 + theme(axis.title.y = element_text(margin = margin(0, -100, 0, 0)))) +
p3 +
  plot_annotation(
  tag_levels = 'a',
  title = 'Variogram range correlations') -> p

# Save plot
cat('Saving variogram model summary plot to:\nfigs/summary/variogram_summary.png\n')

ggsave(
  file = 'figs/summary/variogram_summary.png',
  plot = p,
  device = 'png',
  type = 'cairo',
  width = 8,
  height = 8)

# Interpolation difference
purrr::map(fname, ~get(.x)$diff %>% st_set_geometry(NULL)) %>%
purrr::set_names(nm = seg.names) -> hf.diff

hf.diff %>%
bind_rows(.id = 'Segment') %>%
group_by(Segment) %>%
mutate('Segment' = Segment %>% stringr::str_replace_all('_', ' ')) %>%
summarise(
  Min = round(min(hf.diff)),
  Max = round(max(hf.diff)),
  Median = round(median(hf.diff)),
  IQR = round(IQR(hf.diff)),
  Mean = round(mean(hf.diff)),
  Sigma = round(sd(hf.diff))) -> hf.diff.summary

cat('Heat flow difference summary:\n')
print(hf.diff.summary)

# Visualize
hf.diff %>%
bind_rows(.id = 'segment') %>%
mutate('segment' = segment %>% stringr::str_replace_all('_', ' ')) %>%
group_by(segment) %>%
ggplot() +
geom_boxplot(
  aes(x = hf.diff, y = segment, group = segment),
  width = 0.5,
  outlier.size = 0.5,
  outlier.color = rgb(0.5, 0.5, 0.5, 0.1)) +
labs(
  x = bquote('Heat flow difference'~(mWm^-2)),
  y = NULL,
  title = 'Prediction difference') +
scale_x_continuous(limits = c(-2*max(hf.diff.summary$IQR), 2*max(hf.diff.summary$IQR))) +
scale_y_discrete(
  limits = rev(levels(as.factor(seg.names %>% stringr::str_replace_all('_', ' '))))) +
theme_classic(base_size = 9) +
theme(strip.background = element_blank()) -> p

# Save plot
cat('Saving heat flow difference summary plot to:\nfigs/summary/hf_diff_summary.png')

ggsave(
  file = 'figs/summary/hf_diff_summary.png',
  plot = p,
  device = 'png',
  type = 'cairo',
  width = 4,
  height = 4)

cat('\nDone\n')