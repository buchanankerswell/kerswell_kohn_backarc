# Load functions and libraries
source('functions.R')
load('data/hf.Rdata')
library(ggridges)

# Define paths and names
list.files('data/diff', pattern = '.RData', full.names = T) -> fpaths
purrr::map_chr(
  list.files('data/diff', pattern = '.RData'),
  ~.x %>%
  stringr::str_replace('.RData', '')) -> fnames

# Load data
for (i in fpaths) load(i)

# Bounding Boxes
purrr::map(
  seg.names,
  ~shp.sa.segs.robin.pacific.buffer %>%
  filter(segment == .x) %>%
  st_bbox() %>%
  bbox_widen(
    crs = proj4.robin.pacific,
    borders = c(
      'top' = 0.1,
      'bottom' = 0.1,
      'left' = 0.1,
      'right' = 0.1))) %>%
purrr::set_names(nm = seg.names) -> shp.box

# Crop data
purrr::map(
  shp.box,
  ~shp.hf %>%
  rename(hf = `heat-flow (mW/m2)`) %>%
  st_crop(.x)) %>%
purrr::set_names(nm = seg.names) -> shp.hf.crop

# Color scale
scale_color_viridis_c(
  option = 'magma',
  limits = c(0, 250),
  na.value = 'white') -> viridis.scale.white
scale_color_viridis_c(
  option = 'magma',
  limits = c(0, 250),
  na.value = 'grey50') -> viridis.scale.grey

# Draw plots
purrr::pmap(list(seg.names, fnames), ~{
  cat('Drawing', ..1, '\n')
  shp.hf <- shp.hf.crop[[..1]]
  shp.cont <- shp.sa.contours.robin.pacific %>% st_crop(shp.hf)
  shp.seg <- shp.sa.segs.robin.pacific %>%
    st_crop(shp.hf) %>%
    mutate('segment' = segment %>% stringr::str_replace_all('_', ' '))
  shp.diff <- get(..2)$diff
  shp.krige <- get(..2)$k
  v.grm <- get(..2)$v.grm
  v.mod <- get(..2)$v.mod
  v.line <- variogramLine(v.mod, maxdist = max(v.grm$dist))
  ggplot() +
  geom_sf(data = shp.diff, aes(color = hf.pred.luca),
          size = 2, shape = 15) +
  geom_sf(data = shp.seg, size = 0.85, color = 'white') +
  geom_sf(data = shp.cont, size = 0.1, color = 'white') +
  labs(color = bquote(mWm^-2),
       title = 'Similarity prediction',
       caption = 'data from Lucazeau (2019)') +
  viridis.scale.white +
  coord_sf(expand = F,
           xlim = c(st_bbox(shp.hf)$xmin, st_bbox(shp.hf)$xmax),
           ylim = c(st_bbox(shp.hf)$ymin, st_bbox(shp.hf)$ymax)) +
  theme_map(font_size = 11) +
  theme(axis.text = element_text(),
        panel.grid = element_line(size = 0.1, color = rgb(1,1,1,0.5)),
        panel.ontop = T) -> p.luca.pred
   ggplot() +
   geom_sf(data = shp.diff, aes(color = hf.pred.krige),
           size = 2, shape = 15) +
   geom_sf(data = shp.seg, size = 0.85, color = 'white') +
   geom_sf(data = shp.cont, size = 0.1, color = 'white') +
   labs(color = bquote(mWm^-2),
        title = 'Krige prediction') +
   viridis.scale.white +
   coord_sf(expand = F,
            xlim = c(st_bbox(shp.hf)$xmin, st_bbox(shp.hf)$xmax),
            ylim = c(st_bbox(shp.hf)$ymin, st_bbox(shp.hf)$ymax)) +
   theme_map(font_size = 11) +
   theme(axis.text = element_text(),
         panel.grid = element_line(size = 0.1, color = rgb(1,1,1,0.5)),
         panel.ontop = T) -> p.krige.pred
   ggplot() +
   geom_sf(data = shp.diff, aes(color = abs(hf.diff)),
           size = 2, shape = 15) +
   geom_sf(data = shp.seg, size = 0.85, color = 'white') +
   geom_sf(data = shp.cont, size = 0.1, color = 'white') +
   labs(color = bquote(mWm^-2),
        title = 'Absolute difference') +
   viridis.scale.white +
   coord_sf(expand = F,
            xlim = c(st_bbox(shp.hf)$xmin, st_bbox(shp.hf)$xmax),
            ylim = c(st_bbox(shp.hf)$ymin, st_bbox(shp.hf)$ymax)) +
   theme_map(font_size = 11) +
   theme(axis.text = element_text(),
         panel.grid = element_line(size = 0.1, color = rgb(1,1,1,0.5)),
         panel.ontop = T) -> p.diff.pred
   ggplot() +
   geom_sf(data = shp.world.robin.pacific, color = NA) +
   geom_sf(data = shp.seg, size = 0.85) +
   geom_sf(data = shp.cont, size = 0.1) +
   geom_sf(data = shp.hf, aes(color = hf), shape = 20, size = 0.8) +
   guides(color = guide_colorbar(barheight = 10)) +
   labs(color = bquote(mWm^-2),
        title = bquote(bold('Observations')~n==.(nrow(shp.hf))),
        caption = 'data from Lucazeau (2019)') +
   viridis.scale.grey +
   coord_sf(expand = F,
            xlim = c(st_bbox(shp.hf)$xmin, st_bbox(shp.hf)$xmax),
            ylim = c(st_bbox(shp.hf)$ymin, st_bbox(shp.hf)$ymax)) +
   theme_map(font_size = 11) +
   theme(axis.text = element_text(),
         panel.grid = element_line(size = 0.1, color = rgb(0,0,0,0.2)),
         panel.ontop = T) -> p.pts
   ggplot() +
   geom_histogram(data = shp.diff, aes(x = hf.diff), binwidth = 2) +
   scale_x_continuous(limits = c(median(shp.diff$hf.diff)-(2*IQR(shp.diff$hf.diff)),
                                 median(shp.diff$hf.diff)+(2*IQR(shp.diff$hf.diff)))) +
   labs(x = bquote(mWm^-2),
        y = 'Frequency',
        title = 'Prediction difference') +
   theme_classic() +
   theme(plot.title = element_text(face = 'bold')) -> p.hist
   d.hf.diff <- shp.diff %>% st_set_geometry(NULL)
   d.hf.diff %>%
   select(-hf.obs.luca) %>%
   rename(
     'Similarity Estimate' = hf.pred.luca,
     'Similarity Error' = sigma.luca,
     'Krige Estimate' = hf.pred.krige,
     'Krige Variance' = sigma.krige,
     'Difference' = hf.diff
   ) %>%
   pivot_longer(everything()) %>%
   ggplot(aes(x=value, y = name, fill=factor(stat(quantile)))) +
   ggridges::stat_density_ridges(
     geom='density_ridges_gradient',
     calc_ecdf=T,
     quantiles=4,
     quantile_lines=T
   ) +
   labs(x = bquote('Heat flow '(mWm^-2))) +
   scale_x_continuous(limits = c(-100,250)) +
   scale_fill_viridis_d() +
   theme_classic() +
   theme(
    axis.title.y = element_blank(),
    legend.position = 'none'
   ) -> p.ridge
   ggplot() +
   geom_point(data = v.grm, aes(x = dist/1000, y = gamma), size = 0.8) +
   geom_line(data = v.line, aes(x = dist/1000, y = gamma)) +
   labs(x = 'Lag (km)', y = 'Semivariance',title = 'Variogram') +
   theme_classic() +
   theme(plot.title = element_text(face = 'bold')) -> p.vgrm
  list('p.luca.pred' = p.luca.pred,
       'p.krige.pred' = p.krige.pred,
       'p.diff.pred' = p.diff.pred,
       'p.pts' = p.pts,
       'p.hist' = p.hist,
       'p.vgrm' = p.vgrm)
}) %>%
purrr::set_names(nm = fnames) -> p

tibble(
  lat = c(2.5, 12.5, 2, 6, 27.5, 15, 5, -60, -57.5, -65, -60, -56.5, 5, 20, -10, -20, -12, -17.5),
  lon = c(258, 255, 270, 263, 310, 287.5, 300, 329, 340, 345, 320, 332.5, 115, 130, 90, 171, 175, 174),
  segment = c('Central America', 'Central America', 'Central America', 'Central America', 'Lesser Antilles', 'Lesser Antilles', 'Lesser Antilles', 'Scotia', 'Scotia', 'Scotia', 'Scotia', 'Scotia', 'Sumatra Banda Sea', 'Sumatra Banda Sea', 'Sumatra Banda Sea', 'Vanuatu', 'Vanuatu', 'Vanuatu'),
  label = c('GTJ', 'EPR', 'CR', 'CP', 'MAR', 'CBP', 'SA', 'ESR', 'TF', 'AP', 'SP', 'SAN', 'SNP', 'PSP', 'AUP', 'NHP', 'BR', 'CWR')) %>%
st_as_sf(coords = c(2, 1), crs = proj4.wgs) %>%
st_transform(proj4.robin.pacific)  -> features

cat('\nDrawing custom compositions')

((p$Central_America$p.luca.pred +
 geom_sf_text(
  data = features %>% filter(segment == 'Central America'),
  aes(label = label), color = 'black') +
 theme(
   plot.caption = element_blank(),
   axis.text.x = element_blank()
 )
) +
(p$Central_America$p.krige.pred + theme(axis.text.x = element_blank()))) /
((p$Central_America$p.diff.pred) + p$Central_America$p.dist) +
plot_layout(guides = 'collect') +
plot_annotation(tag_levels = 'a') &
guides(color = guide_colorbar(barwidth = unit(2, 'in'))) &
theme(
  plot.title = element_text(hjust = 1, size = 10),
  axis.text = element_text(size = 7),
  legend.title = element_text(vjust = 1, size = 10),
  legend.box.margin = margin(),
  legend.text = element_text(size = 8),
  legend.key.height = unit(0.7, 'lines'),
  legend.position = 'bottom',
  legend.justification = 'right'
) -> p1

# Save
ggsave(file = 'figs/diff/custom/Central_America.png',
       plot = p1,
       device = 'png',
       type = 'cairo',
       width = 6,
       height = 5.5
)

(p$Lesser_Antilles$p.luca.pred +
 geom_sf_text(data = features %>% filter(segment == 'Lesser Antilles'),
              aes(label = label), color = 'white') +
 theme(plot.caption = element_blank())) +
(p$Lesser_Antilles$p.krige.pred + theme(axis.text.y = element_blank())) +
(p$Lesser_Antilles$p.diff.pred + theme(axis.text.y = element_blank())) +
plot_layout(guides = 'collect') +
plot_annotation(tag_levels = 'a') &
guides(color = guide_colorbar(barwidth = unit(50, 'mm'))) &
theme(
  plot.title = element_text(hjust = 1, size = 10),
  axis.text = element_text(size = 7),
  legend.title = element_text(vjust = 1, size = 10),
  legend.box.margin = margin(),
  legend.text = element_text(size = 8),
  legend.key.height = unit(0.7, 'lines'),
  legend.position = 'bottom',
  legend.justification = 'right'
) -> p2

# Save
ggsave(file = 'figs/diff/custom/Lesser_Antilles.png',
       plot = p2,
       device = 'png',
       type = 'cairo',
       width = 6,
       height = 3
)

(p$Scotia$p.luca.pred +
 geom_sf_text(data = features %>% filter(segment == 'Scotia'),
              aes(label = label), color = 'white') +
 theme(plot.caption = element_blank())) +
(p$Scotia$p.krige.pred + theme(axis.text.y = element_blank())) +
(p$Scotia$p.diff.pred + theme(axis.text.y = element_blank())) +
plot_layout(guides = 'collect') +
plot_annotation(tag_levels = 'a') &
guides(color = guide_colorbar(barwidth = unit(50, 'mm'))) &
theme(
  plot.title = element_text(hjust = 1, size = 10),
  axis.text = element_text(size = 7),
  legend.title = element_text(vjust = 1, size = 10),
  legend.box.margin = margin(),
  legend.text = element_text(size = 8),
  legend.key.height = unit(0.7, 'lines'),
  legend.position = 'bottom',
  legend.justification = 'right'
) -> p3

# Save
ggsave(file = 'figs/diff/custom/Scotia.png',
       plot = p3,
       device = 'png',
       type = 'cairo',
       width = 6,
       height = 2.5
)

(p$Sumatra_Banda_Sea$p.luca.pred +
 geom_sf_text(data = features %>% filter(segment == 'Sumatra Banda Sea'),
              aes(label = label), color = 'white', size = 3) +
 theme(plot.caption = element_blank())) +
(p$Sumatra_Banda_Sea$p.krige.pred + theme(axis.text.y = element_blank())) +
(p$Sumatra_Banda_Sea$p.diff.pred + theme(axis.text.y = element_blank())) +
plot_layout(guides = 'collect') +
plot_annotation(tag_levels = 'a') &
guides(color = guide_colorbar(barwidth = unit(50, 'mm'))) &
theme(
  plot.title = element_text(hjust = 1, size = 10),
  axis.text = element_text(size = 5),
  legend.title = element_text(vjust = 1, size = 10),
  legend.text = element_text(size = 8),
  legend.box.margin = margin(),
  legend.key.height = unit(0.7, 'lines'),
  legend.position = 'bottom',
  legend.justification = 'right'
) -> p4

# Save
ggsave(file = 'figs/diff/custom/Sumatra_Banda_Sea.png',
       plot = p4,
       device = 'png',
       type = 'cairo',
       width = 6,
       height = 2.5
)

(p$Vanuatu$p.luca.pred + theme(plot.caption = element_blank())) +
(p$Vanuatu$p.krige.pred +
  geom_sf_text(data = features %>% filter(segment == 'Vanuatu'),
              aes(label = label), color = 'white') +
   theme(axis.text.y = element_blank())) +
(p$Vanuatu$p.diff.pred + theme(axis.text.y = element_blank())) +
plot_layout(guides = 'collect') +
plot_annotation(tag_levels = 'a') &
guides(color = guide_colorbar(barwidth = unit(50, 'mm'))) &
theme(
  plot.title = element_text(hjust = 1, size = 10),
  axis.text = element_text(size = 7),
  legend.title = element_text(vjust = 1, size = 10),
  legend.text = element_text(size = 8),
  legend.box.margin = margin(),
  legend.key.height = unit(0.7, 'lines'),
  legend.position = 'bottom',
  legend.justification = 'right'
) -> p5

# Save
ggsave(file = 'figs/diff/custom/Vanuatu.png',
       plot = p5,
       device = 'png',
       type = 'cairo',
       width = 6,
       height = 3
)

cat('\nDone\n')
