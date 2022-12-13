#!/usr/bin/env Rscript

# Capture output
sink(file = paste0('data/log-', Sys.Date()), append = T, type = 'output', split = T)

cat(rep('~', 60), '\n', sep='')
cat('Loading packages and functions ...\n\n')

source('R/functions.R')
load('data/hf.RData')

# Visualize
cat('\n', rep('~', 60), sep='')
cat('\nVisualizing ...')

# Global ThermoGlobe buffer
p1 <-
  ggplot() +
  geom_sf(data = shp.world, size = 0.1, fill = 'grey70', color = 'black') +
  geom_sf(data = st_union((bind_rows(shp.buffer))), size = 0.3, fill = 'black', alpha = 0.3) +
  geom_sf(data = shp.ridge, size = 0.4, color = 'black', alpha = 0.8) +
  geom_sf(data = shp.trench, size = 0.4, color = 'black', alpha = 0.8) +
  geom_sf(data = shp.transform, size = 0.4, color = 'black', alpha = 0.8) +
  geom_sf(data = bind_rows(shp.hf.crop), aes(color = hf), size = 0.3, shape = 20) +
  geom_sf(data = bind_rows(shp.segs), size = 0.8, color = 'white') +
  annotate(
    'label',
    x = -Inf,
    y = -Inf,
    label = 'ThermoGlobe dataset',
    size = 4,
    hjust = 0,
    vjust = 0,
    fill = 'grey90',
    label.padding = unit(0.15, 'lines'),
    label.r = unit(0, 'lines')
  ) +
  scale_color_viridis_c(
    option = 'magma',
    limits = c(0, 250),
    breaks = c(0, 125, 250),
    na.value = 'transparent'
  ) +
  scale_x_continuous(breaks = seq(-180, 180, 60)) +
  labs(color = bquote(mWm^-2)) +
  coord_sf(expand = F) +
  theme_map(font_size = 12) +
  theme(
    axis.text = element_text(),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.direction = 'horizontal',
    legend.box.background = element_rect(fill = rgb(1, 1, 1, 0.8), color = NA),
    legend.box.margin = margin(0, 8, 0, 0),
    legend.margin = margin(),
    legend.key.height = unit(0.125, 'in'),
    legend.key.width = unit(0.15, 'in'),
    legend.title = element_text(vjust = 1, color = 'black', size = 10),
    panel.grid = element_line(size = 0.1, color = 'white'),
    panel.background = element_rect(fill = 'grey50', color = NA),
    plot.margin = margin()
  )

# Crop countries and similarity to buffer
shp.world.buf <- suppressWarnings(shp.world %>% st_intersection(bind_rows(shp.buffer)))
shp.sim <- suppressWarnings(shp.interp.luca %>% st_intersection(bind_rows(shp.buffer)))
p2 <-
  ggplot() +
  geom_sf(data = shp.world, size = 0.1, fill = 'grey70', color = 'black') +
  geom_sf(data = shp.sim, aes(color = est.sim), size = 0.1, shape = 15) +
  geom_sf(data = shp.world.buf, size = 0.1, fill = 'grey70', alpha = 0.1) +
  geom_sf(data = shp.ridge, size = 0.4, color = 'black', alpha = 0.8) +
  geom_sf(data = shp.trench, size = 0.4, color = 'black', alpha = 0.8) +
  geom_sf(data = shp.transform, size = 0.4, color = 'black', alpha = 0.8) +
  geom_sf(data = bind_rows(shp.segs), size = 0.8, color = 'white') +
  geom_sf_label_repel(
    data = bind_rows(shp.segs),
    aes(label =
      c('AA', 'AN', 'CA', 'KM', 'KR', 'LA', 'NP', 'NBS', 'SP', 'SC', 'SBS', 'TNZ', 'VN')
    ),
    size = 3,
    segment.size = 0.3,
    segment.color = 'white',
    segment.curvature = -1e-20,
    arrow = arrow(length = unit(0.015, 'npc')),
    fill = rgb(1, 1, 1, 0.8),
    label.padding = unit(0.15, 'lines'),
    label.r = unit(0, 'lines'),
    force = 6,
    seed = 19
  ) +
  annotate(
    'label',
    x = -Inf,
    y = -Inf,
    label = 'Similarity interpolation (Lucazeau, 2019)',
    hjust = 0,
    vjust = 0,
    size = 4,
    fill = 'grey90',
    label.padding = unit(0.15, 'lines'),
    label.r = unit(0, 'lines')
  ) +
  labs(color = bquote(mWm^-2)) +
  scale_color_viridis_c(
    option = 'magma',
    limits = c(0, 250),
    breaks = c(0, 125, 250),
    na.value = 'transparent'
  ) +
  coord_sf(expand = F) +
  theme_map(font_size = 12) +
  theme(
    axis.text = element_text(),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.direction = 'horizontal',
    legend.box.background = element_rect(fill = rgb(1, 1, 1, 0.8), color = NA),
    legend.box.margin = margin(0, 8, 0, 0),
    legend.margin = margin(),
    legend.key.height = unit(0.125, 'in'),
    legend.key.width = unit(0.15, 'in'),
    legend.title = element_text(vjust = 1, color = 'black', size = 10),
    panel.grid = element_line(size = 0.1, color = 'white'),
    panel.background = element_rect(fill = 'grey50', color = NA),
    plot.margin = margin(1, 1, 1, 1)
  )

# Global buffer composition
p3 <-
  (p1 +
    annotate(
      'label',
      label = 'a',
      x = -Inf,
      y = Inf,
      size = 5,
      hjust = 0,
      vjust = 1,
      fill = 'grey90',
      label.padding = unit(0.15, 'lines'),
      label.r = unit(0, 'in')
    ) +
    theme(axis.text = element_blank())
  ) /
  (p2 +
    annotate(
      'label',
      label = 'b',
      x = -Inf,
      y = Inf,
      size = 5,
      hjust = 0,
      vjust = 1,
      fill = 'grey90',
      label.padding = unit(0.15, 'lines'),
      label.r = unit(0, 'in')
    ) +
    theme(legend.position = 'none')
  ) &
  theme(plot.margin = margin(1, 1, 1, 1))
# Save
cat('\nSaving plot to: figs/base/ThermoGlobeBufferComp.png')
suppressWarnings(suppressMessages(
  ggsave(
    file = 'figs/base/ThermoGlobeBufferComp.png',
    plot = p3,
    device = 'png',
    type = 'cairo',
    width = 6.5,
    height = 6.5
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
  hf <- shp.hf.crop[[.x]]
  sim <- suppressWarnings(shp.interp.luca %>% st_intersection(buf)) # Similarity interp
  ridge <- shp.ridge.crop[[.x]]
  trench <- shp.trench.crop[[.x]]
  transform <- shp.transform.crop[[.x]]
  fts <- shp.fts[shp.fts$segment == .x,]
  grats.x <-
    as_tibble(st_coordinates(shp.grid.crop[[.x]])) %>%
    group_by(Y) %>%
    group_split() %>%
    bind_rows(.id = 'line') %>%
    mutate(line = as.integer(line))
  grats.spacing <-
    mean(
      map_dbl(
        seq_len(max(grats.x$line)-1),
        ~{abs(
          grats.x$Y[grats.x$line == .x][1] -
          grats.x$Y[grats.x$line == .x+1][1]
        )}
      )
    )
  grats <-
    st_make_grid(shp.grid.crop[[.x]], cellsize = c(grats.spacing, grats.spacing)) %>%
    st_intersection(buf)
  pnt.size <- 3
  annt.txt.size <- 4
  base.txt.size <- 12
  # Define map scale 1:50,000 [meters]
  wdth <- (st_bbox(buf)$xmax - st_bbox(buf)$xmin)/5e4
  hght <- (st_bbox(buf)$ymax - st_bbox(buf)$ymin)/5e4
  aspect <- wdth/hght
  if(aspect <= 0.4) {
    const <- 165.1/hght
    p.wdth <- wdth * const * 2
    p.hght <- hght * const
  } else if(aspect > 0.4 & aspect <= 0.95) {
    const <- 82.55/wdth
    p.wdth <- wdth * const * 2
    p.hght <- hght * const
  } else if(aspect >= 2) {
    const <- 165.1/wdth
    p.wdth <- wdth * const
    p.hght <- hght * const * 2
  } else {
    const <- 82.55/hght
    p.wdth <- wdth * const
    p.hght <- hght * const * 2
  }
  # Interpolation domain
  pp1 <- 
    ggplot() +
    geom_sf(data = world, size = 0.1, fill = 'grey70', color = 'black') +
    geom_sf(data = grats, size = 0.1, color = 'black', fill = NA, alpha = 0.5) +
    geom_sf(data = ridge, size = 1, color = 'black', alpha = 0.8) +
    geom_sf(data = trench, size = 1, color = 'black', alpha = 0.8) +
    geom_sf(data = transform, size = 1, color = 'black', alpha = 0.8) +
    geom_sf(data = buf, size = 0.3, fill = NA) +
    geom_sf(data = seg, size = 1.5, color = 'white') +
    geom_sf(data = hf, aes(color = hf), shape = 15, size = pnt.size*0.3) +
    geom_sf(data = volc, size = pnt.size*0.3, color = 'gold', shape = 18) +
    geom_sf_label(
      data = fts,
      aes(label = label),
      size = annt.txt.size,
      fill = rgb(1, 1, 1, 0.8),
      label.padding = unit(0.15, 'lines'),
      label.r = unit(0.05, 'lines')
    ) +
    annotate(
      'label',
      label = 'a',
      x = -Inf,
      y = Inf,
      size = 5,
      hjust = 0,
      vjust = 1,
      fill = 'grey90',
      label.padding = unit(0.15, 'lines'),
      label.r = unit(0, 'in')
    ) +
    scale_color_viridis_c(
      option = 'magma',
      limits = c(0, 250),
      breaks = c(0, 125, 250),
      na.value = 'transparent'
    ) +
    labs(color = bquote(mWm^-2)) +
    coord_sf(expand = F) +
    theme_map(font_size = base.txt.size) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
      axis.text.y = element_text(angle = 30, hjust = 1),
      legend.position = c(1, 1),
      legend.justification = c(1, 1),
      legend.direction = 'horizontal',
      legend.box.background = element_rect(fill = rgb(1, 1, 1, 0.8), color = NA),
      legend.box.margin = margin(1, 8, 1, 2),
      legend.key.height = unit(0.125, 'in'),
      legend.key.width = unit(0.15, 'in'),
      legend.title = element_text(vjust = 1),
      panel.grid = element_line(size = 0.1, color = 'white'),
      panel.background = element_rect(fill = 'grey50', color = NA),
      plot.margin = margin(1, 1, 1, 1)
    )
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
    geom_sf(data = world, size = 0.1, fill = 'grey70', color = 'black') +
    geom_sf(data = sim, aes(color = est.sim), size = pnt.size, shape = 15, show.legend = F) +
    geom_sf(data = world.buf, size = 0.1, color = 'grey80', fill = 'grey70', alpha = 0.1) +
    geom_sf(data = ridge, size = 1, color = 'black', alpha = 0.8) +
    geom_sf(data = trench, size = 1, color = 'black', alpha = 0.8) +
    geom_sf(data = transform, size = 1, color = 'black', alpha = 0.8) +
    geom_sf(data = buf, size = 0.3, fill = NA) +
    geom_sf(data = seg, size = 1.5, color = 'white') +
    annotate(
      'label',
      label = 'b',
      x = -Inf,
      y = Inf,
      size = 5,
      hjust = 0,
      vjust = 1,
      fill = 'grey90',
      label.padding = unit(0.15, 'lines'),
      label.r = unit(0, 'in')
    ) +
    scale_color_viridis_c(
      option = 'magma',
      limits = c(0, 250),
      breaks = c(0, 125, 250),
      na.value = 'transparent'
    ) +
    labs(color = bquote(mWm^-2)) +
    coord_sf(expand = F) +
    theme_map(font_size = base.txt.size) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
      axis.text.y = element_text(angle = 30, hjust = 1),
      legend.position = c(1, 1),
      legend.justification = c(1, 1),
      legend.direction = 'horizontal',
      legend.box.background = element_rect(fill = rgb(1, 1, 1, 0.8), color = NA),
      legend.box.margin = margin(0, 8, 0, 0),
      legend.key.height = unit(0.125, 'in'),
      legend.key.width = unit(0.15, 'in'),
      legend.title = element_text(vjust = 1),
      panel.grid = element_line(size = 0.1, color = 'white'),
      panel.background = element_rect(fill = 'grey50', color = NA),
      plot.margin = margin(1, 1, 1, 1)
    )
  if(
     .x %in% c('Alaska Aleutians', 'Kamchatka Marianas', 'Tonga New Zealand', 'Vanuatu')
   ) {
    pp2 <-
      pp2 +
      scale_x_continuous(breaks = c(130, 140, 150, 160, 170, 180, -170, -160, -150, -140, -130))
  }
  if(aspect <= 0.95) {
    # Composition
    p <-
      (pp1 +
        theme(
          axis.text.y = element_text(angle = 30, hjust = 1),
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
        )
      ) +
      (pp2 +
       theme(
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
       )
      ) &
      theme(
        plot.margin = margin(1, 1, 1, 1),
        panel.grid = element_line(size = 0.1, color = 'white'),
        panel.background = element_rect(fill = 'grey50', color = NA),
        legend.box.margin = margin(),
        legend.text = element_text(margin = margin(t = -4)),
        legend.margin = margin(0, 8, 0, 0)
      )
    # Save
    cat('\nSaving plot to: figs/base/', str_replace_all(.x, ' ', ''), 'Comp.png', sep = '')
    suppressWarnings(suppressMessages(
      ggsave(
        file = paste0('figs/base/', str_replace_all(.x, ' ', ''), 'Comp.png'),
        plot = p,
        device = 'png',
        type = 'cairo',
        width = p.wdth,
        height = p.hght,
        units = 'mm'
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
        panel.grid = element_line(size = 0.1, color = 'white'),
        panel.background = element_rect(fill = 'grey50', color = NA),
        legend.box.margin = margin(),
        legend.text = element_text(margin = margin(t = -4)),
        legend.margin = margin(0, 8, 0, 0)
      )
    # Save
    cat('\nSaving plot to: figs/base/', str_replace_all(.x, ' ', ''), 'Comp.png', sep = '')
    suppressWarnings(suppressMessages(
      ggsave(
        file = paste0('figs/base/', str_replace_all(.x, ' ', ''), 'Comp.png'),
        plot = p,
        device = 'png',
        type = 'cairo',
        width = p.wdth,
        height = p.hght,
        units = 'mm'
      )
    ))
  }
})

cat('\nbase-plots.R complete!\n\n')
sink()