#!/usr/bin/env Rscript

# Load packages and functions
cat(rep('~', 45), '\n', sep='')
cat('Loading packages and functions ...\n\n')
source('R/functions.R')
load('data/hf.RData')

# Create directory
dir.create('figs/base', recursive = T, showWarnings = F)

# Visualize
cat('\n', rep('~', 45), sep='')
cat('\nVisualizing ...')

# Global ThermoGlobe buffer
p1 <-
  ggplot() +
  ggtitle('a) Thermoglobe observations') +
  geom_sf(data = shp.relief.world, aes(color = elevation), shape = 15, size = 0.01) +
  scale_color_etopo(guide = 'none') +
  new_scale_color() +
  geom_sf(data = st_union((bind_rows(shp.buffer))), linewidth = 0.3, fill = NA) +
  geom_sf(data = shp.ridge, linewidth = 0.3) +
  geom_sf(data = shp.trench, linewidth = 0.3) +
  geom_sf(data = shp.transform, linewidth = 0.3) +
  geom_sf(data = bind_rows(shp.hf.crop), aes(color = hf), size = 0.3, shape = 20) +
  geom_sf(data = bind_rows(shp.segs), linewidth = 1, color = 'white') +
  scale_color_viridis_c(
    option = 'magma',
    name = bquote(mWm^-2),
    limits = c(0, 250),
    breaks = c(0, 125, 250),
    na.value = 'transparent',
    guide = guide_colorbar(title.vjust = 1, show.limits = T)
  ) +
  scale_x_continuous(breaks = seq(-180, 180, 60)) +
  labs(color = bquote(mWm^-2)) +
  coord_sf() +
  theme_map(font_size = 14)
# Crop countries and similarity to buffer
shp.world.buf <- suppressWarnings(shp.world %>% st_intersection(bind_rows(shp.buffer)))
shp.sim <- suppressWarnings(shp.interp.luca %>% st_intersection(bind_rows(shp.buffer)))
p2 <-
  ggplot() +
  ggtitle('b) Similarity predictions') +
  geom_sf(data = shp.relief.world, aes(color = elevation), shape = 15, size = 0.01) +
  scale_color_etopo(guide = 'none') +
  new_scale_color() +
  geom_sf(data = shp.sim, aes(color = est.sim), size = 0.1, shape = 15) +
  geom_sf(data = shp.world.buf, linewidth = 0.1, fill = 'grey70', alpha = 0.1) +
  geom_sf(data = shp.ridge, linewidth = 0.3) +
  geom_sf(data = shp.trench, linewidth = 0.3) +
  geom_sf(data = shp.transform, linewidth = 0.3) +
  geom_sf(data = bind_rows(shp.segs), linewidth = 1, color = 'white') +
  geom_sf_label_repel(
    data = bind_rows(shp.segs),
    aes(label =
      c('AA', 'AN', 'CA', 'KM', 'KR', 'LA', 'NP', 'NBS', 'SP', 'SC', 'SBS', 'TNZ', 'VN')
    ),
    size = 4,
    fill = rgb(1, 1, 1, 0.9),
    label.padding = unit(0.15, 'lines'),
    label.r = unit(0, 'lines'),
    seed = 19
  ) +
  labs(color = bquote(mWm^-2)) +
  scale_color_viridis_c(
    option = 'magma',
    name = bquote(mWm^-2),
    limits = c(0, 250),
    breaks = c(0, 125, 250),
    na.value = 'transparent',
    guide = guide_colorbar(title.vjust = 1, show.limits = T)
  ) +
  coord_sf() +
  theme_map(font_size = 14)
# Global buffer composition
p3 <-
  (p1 + theme(axis.text = element_blank())) / p2 &
  theme(
    plot.margin = margin(1, 1, 1, 1),
    legend.position = 'top',
    legend.justification = 'right',
    legend.direction = 'horizontal',
    axis.text = element_text(hjust = 1),
    legend.margin = margin(-4, 0, -12, 0),
    legend.box.margin = margin(0, 10, 0, 0),
    legend.key.height = unit(0.125, 'in'),
    legend.key.width = unit(0.2, 'in'),
    legend.title = element_text(vjust = 0, color = 'black', size = 14),
    panel.grid = element_line(size = 0.01, color = 'grey60'),
    plot.title = element_text(vjust = 0, margin = margin(0, 0, -10, 0))
  )
# Save
cat('\nSaving plot to: figs/base/ThermoGlobeBufferComp.png')
suppressWarnings(suppressMessages(
  ggsave(
    file = 'figs/base/ThermoGlobeBufferComp.png',
    plot = p3,
    device = 'png',
    type = 'cairo',
    width = 6.5,
    height = 6.5,
    dpi = 330
  )
))

# Individual Segments
seg.names %>% walk(~{
  # Define map parts
  cnt <- shp.contours[[.x]] # Contour
  seg <- shp.segs[[.x]] # Segment
  buf <- shp.buffer[[.x]] # Buffer
  grd <- shp.grid.crop[[.x]] # Interpolation grid
  bbx <- st_bbox(buf) %>% bbox_widen(borders = c(0.05, 0.05, 0.05, 0.05)) # Bounding box
  world <- suppressWarnings(shp.world %>% st_crop(bbx)) # Countries
  world.buf <- suppressWarnings(world %>% st_intersection(buf)) # Contries within buffer
  volc <- suppressWarnings(shp.volc %>% st_intersection(buf)) # Contries within buffer
  hf <- shp.hf.crop[[.x]] # Heat flow data
  sim <- suppressWarnings(shp.interp.luca %>% st_intersection(buf)) # Similarity interp
  ridge <- shp.ridge.crop[[.x]]
  trench <- shp.trench.crop[[.x]]
  transform <- shp.transform.crop[[.x]]
  relief <- shp.relief.crop[[.x]]
  fts <- shp.fts[shp.fts$segment == .x,]
  pnt.size <- 1
  annt.txt.size <- 7
  base.txt.size <- 14
  # Define map scale 1:50,000 [meters]
  wdth <- (st_bbox(buf)$xmax - st_bbox(buf)$xmin)/5e4
  hght <- (st_bbox(buf)$ymax - st_bbox(buf)$ymin)/5e4
  aspect <- wdth/hght
  # Interpolation domain
  pp1 <- 
    ggplot() +
    ggtitle('a) Observations') +
    geom_sf(data = relief, aes(color = elevation), shape = 15, size = 0.01) +
    scale_color_etopo(guide = 'none') +
    new_scale_color() +
    geom_sf(data = ridge, linewidth = 0.5) +
    geom_sf(data = trench, linewidth = 0.5) +
    geom_sf(data = transform, linewidth = 0.5) +
    geom_sf(data = seg, linewidth = 1, color = 'white') +
    geom_sf(data = hf, aes(color = hf), shape = 15, size = pnt.size) +
    geom_sf(data = volc, size = pnt.size, color = 'white', shape = 18) +
    scale_color_viridis_c(
      option = 'magma',
      name = bquote(mWm^-2),
      limits = c(0, 250),
      breaks = c(0, 125, 250),
      na.value = 'transparent',
      guide = guide_colorbar(title.vjust = 1, show.limits = T)
    ) +
    coord_sf() +
    theme_map(font_size = base.txt.size)
  if(
     .x %in% c('Alaska Aleutians', 'Kamchatka Marianas', 'Tonga New Zealand', 'Vanuatu')
   ) {
    pp1 <-
      pp1 +
      scale_x_continuous(breaks = c(130, 140, 150, 160, 170, 180, -170, -160, -150, -140, -130))
  }
  # Similarity interpolation
  pp2 <- 
    ggplot() +
    ggtitle('b) Predictions') +
    geom_sf(data = relief, aes(color = elevation), shape = 15, size = 0.01) +
    scale_color_etopo(guide = 'none') +
    new_scale_color() +
    geom_sf(data = sim, aes(color = est.sim), size = pnt.size*1.5, shape = 15) +
    geom_sf(data = world.buf, linewidth = 0.1, fill = NA, color = 'grey60') +
    geom_sf(data = ridge, linewidth = 0.5) +
    geom_sf(data = trench, linewidth = 0.5) +
    geom_sf(data = transform, linewidth = 0.5) +
    geom_sf(data = buf, linewidth = 0.1, fill = NA, color = 'grey60') +
    geom_sf(data = seg, linewidth = 1, color = 'white') +
    geom_sf_label(
      data = fts,
      aes(label = label),
      size = annt.txt.size*0.5,
      fill = rgb(1, 1, 1, 0.9),
      label.padding = unit(0.15, 'lines'),
      label.r = unit(0.05, 'lines')
    ) +
    scale_color_viridis_c(
      option = 'magma',
      name = bquote(mWm^-2),
      limits = c(0, 250),
      breaks = c(0, 125, 250),
      na.value = 'transparent',
      guide = guide_colorbar(title.vjust = 1, show.limits = T)
    ) +
    coord_sf() +
    theme_map(font_size = base.txt.size)
  if(
     .x %in% c('Alaska Aleutians', 'Kamchatka Marianas', 'Tonga New Zealand', 'Vanuatu')
   ) {
    pp2 <-
      pp2 +
      scale_x_continuous(
        breaks = c(130, 140, 150, 160, 170, 180, -170, -160, -150, -140, -130)
      )
  }
  if(aspect <= 0.95) {
    # Composition
    p <-
      (pp1 +
        theme(
          axis.text.x =
            element_text(angle = 30, hjust = 1, vjust = 1, margin = margin(5, 0, 0, 0)),
          axis.text.y = element_text(angle = 30, hjust = 1, margin = margin(0, 5, 0, 0))
        )
      ) +
      (pp2 +
       theme(
          axis.text.x =
            element_text(angle = 30, hjust = 1, vjust = 1, margin = margin(5, 0, 0, 0)),
          axis.text.y = element_blank()
       )
      ) &
      theme(
        plot.margin = margin(1, 1, 1, 1),
        panel.grid = element_line(linewidth = 0.01, color = 'grey60'),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.justification = 'center',
        legend.box.margin = margin(),
        legend.text = element_text(margin = margin(t = -4)),
        legend.title = element_text(vjust = 0, color = 'black', size = base.txt.size)
      )
    # Save
    cat('\nSaving plot to: figs/base/', str_replace_all(.x, ' ', ''), 'Comp.png', sep = '')
    suppressWarnings(suppressMessages(
      ggsave(
        file = paste0('figs/base/', str_replace_all(.x, ' ', ''), 'Comp.png'),
        plot = p,
        device = 'png',
        type = 'cairo',
        width = 6.5,
        height = 6.5,
        dpi = 330
      )
    ))
  } else {
    # Composition
    p <-
      (pp1 +
        theme(
          axis.text.x = element_blank(),
          axis.text.y = element_text(angle = 30, hjust = 1)
        )
      ) /
      (pp2 +
        theme(
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
          axis.text.y = element_text(angle = 30, hjust = 1)
        )
      ) &
      theme(
        plot.margin = margin(1, 1, 1, 1),
        panel.grid = element_line(linewidth = 0.01, color = 'grey60'),
        legend.justification = 'center',
        legend.box.margin = margin(),
        legend.text = element_text(margin = margin(t = -4)),
        legend.title = element_text(vjust = 0, color = 'black', size = base.txt.size)
      )
    # Save
    cat('\nSaving plot to: figs/base/', str_replace_all(.x, ' ', ''), 'Comp.png', sep = '')
    suppressWarnings(suppressMessages(
      ggsave(
        file = paste0('figs/base/', str_replace_all(.x, ' ', ''), 'Comp.png'),
        plot = p,
        device = 'png',
        type = 'cairo',
        width = 6.5,
        height = 6.5,
        dpi = 330
      )
    ))
  }
})

cat('\nbase-plots.R complete!\n\n')