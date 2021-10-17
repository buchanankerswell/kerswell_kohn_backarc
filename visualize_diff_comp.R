# Load functions and libraries
source('functions.R')
load('data/hf.Rdata')

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

# Plot sizes
p.widths <- c(11, 11, 8, 8, 8, 8, 8, 8, 8, 8, 11, 8, 8)
p.heights <- c(8, 9, 8, 12, 9, 11, 12, 8.5, 13, 9, 11, 14, 11)

plt <- function(fname, p.height, p.width){
  cat('Plotting', fname, '\n')
  shp.hf <- shp.hf.crop[[fname]]
  shp.cont <- shp.sa.contours.robin.pacific %>% st_crop(shp.hf)
  shp.seg <- shp.sa.segs.robin.pacific %>%
    st_crop(shp.hf) %>%
    mutate('segment' = segment %>% stringr::str_replace_all('_', ' '))
  shp.diff <- get(fname)$diff
  shp.krige <- get(fname)$k
  v.grm <- get(fname)$v.grm
  v.mod <- get(fname)$v.mod
  v.line <- variogramLine(v.mod, maxdist = max(v.grm$dist))
  ggplot() +
  geom_sf(data = shp.diff, aes(color = hf.pred.luca),
          size = 2, shape = 15) +
  geom_sf(data = shp.seg, size = 1.1, color = 'white') +
  geom_sf(data = shp.cont, size = 0.15, color = 'white') +
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
   geom_sf(data = shp.seg, size = 1.1, color = 'white') +
   geom_sf(data = shp.cont, size = 0.15, color = 'white') +
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
   geom_sf(data = shp.seg, size = 1.1, color = 'white') +
   geom_sf(data = shp.cont, size = 0.15, color = 'white') +
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
   geom_sf(data = shp.seg, size = 1.1) +
   geom_sf(data = shp.cont, size = 0.15) +
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
   ggplot() +
   geom_point(data = v.grm, aes(x = dist/1000, y = gamma), size = 0.8) +
   geom_line(data = v.line, aes(x = dist/1000, y = gamma)) +
   labs(x = 'Lag (km)', y = 'Semivariance', title = 'Variogram') +
   theme_classic() +
   theme(plot.title = element_text(face = 'bold')) -> p.vgrm
  # Composition
  (p.luca.pred + theme(legend.position = 'none')) +
  (p.krige.pred + theme(legend.position = 'none')) +
  (p.diff.pred + theme(legend.position = 'none')) +
  p.pts +
  p.hist +
  p.vgrm +
  plot_annotation(
    title = fname %>% stringr::str_replace_all('_', ' '),
    theme = theme(legend.position = 'right',
                  legend.direction = 'vertical')) -> p
  if(fname != 'Andes') {
    p +
    plot_layout(
      nrow = 3,
      ncol = 2,
      widths = 1,
      heights = c(1, 1, 0.5),
      guides = 'collect') -> p.comp
  } else {
    design <-
      '11223344
      55556666'
     p + plot_layout(design = design, heights = c(1, 0.5)) -> p.comp
  }
  # Save
  cat('Saving plot to', paste0('figs/diff/comp/', fname, '.png'), '\n')
  ggsave(file = paste0('figs/diff/comp/', fname, '.png'),
         device = 'png',
         type = 'cairo',
         plot = p.comp,
         height = p.height,
         width = p.width)
}

# Draw plots
purrr::pwalk(list(fnames, p.heights, p.widths), ~{
  cat('Plotting', ..1, '\n')
  shp.hf <- shp.hf.crop[[..1]]
  shp.cont <- shp.sa.contours.robin.pacific %>% st_crop(shp.hf)
  shp.seg <- shp.sa.segs.robin.pacific %>%
    st_crop(shp.hf) %>%
    mutate('segment' = segment %>% stringr::str_replace_all('_', ' '))
  shp.diff <- get(..1)$diff
  shp.krige <- get(..1)$k
  v.grm <- get(..1)$v.grm
  v.mod <- get(..1)$v.mod
  v.line <- variogramLine(v.mod, maxdist = max(v.grm$dist))
  ggplot() +
  geom_sf(data = shp.diff, aes(color = hf.pred.luca),
          size = 2, shape = 15) +
  geom_sf(data = shp.seg, size = 1.1, color = 'white') +
  geom_sf(data = shp.cont, size = 0.15, color = 'white') +
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
   geom_sf(data = shp.seg, size = 1.1, color = 'white') +
   geom_sf(data = shp.cont, size = 0.15, color = 'white') +
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
   geom_sf(data = shp.seg, size = 1.1, color = 'white') +
   geom_sf(data = shp.cont, size = 0.15, color = 'white') +
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
   geom_sf(data = shp.seg, size = 1.1) +
   geom_sf(data = shp.cont, size = 0.15) +
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
   ggplot() +
   geom_point(data = v.grm, aes(x = dist/1000, y = gamma), size = 0.8) +
   geom_line(data = v.line, aes(x = dist/1000, y = gamma)) +
   labs(x = 'Lag (km)', y = 'Semivariance', title = 'Variogram') +
   theme_classic() +
   theme(plot.title = element_text(face = 'bold')) -> p.vgrm
  # Composition
  (p.luca.pred + theme(legend.position = 'none')) +
  (p.krige.pred + theme(legend.position = 'none')) +
  (p.diff.pred + theme(legend.position = 'none')) +
  p.pts +
  p.hist +
  p.vgrm +
  plot_annotation(
    title = ..1 %>% stringr::str_replace_all('_', ' '),
    theme = theme(legend.position = 'right',
                  legend.direction = 'vertical')) -> p
  if(..1 != 'Andes') {
    p +
    plot_layout(
      nrow = 3,
      ncol = 2,
      widths = 1,
      heights = c(1, 1, 0.5),
      guides = 'collect') -> p.comp
  } else {
    design <-
      '11223344
      55556666'
     p + plot_layout(design = design, heights = c(1, 0.5)) -> p.comp
  }
  # Save
  cat('Saving plot to', paste0('figs/diff/comp/', ..1, '.png'), '\n')
  ggsave(file = paste0('figs/diff/comp/', ..1, '.png'),
         device = 'png',
         type = 'cairo',
         plot = p.comp,
         height = ..2,
         width = ..3)
})

cat('\nDone\n')
