#!/usr/bin/env Rscript

cat(rep('~', 60), '\n', sep='')
cat('Loading packages and functions ...\n\n')

source('functions.R')
load('data/hf.RData')
dir.create('figs/base', showWarnings = F)

# Visualize
cat('\n', rep('~', 60), sep='')
cat('\nVisualizing ...')

# World map
p1 <-
  ggplot(bind_rows(shp.segs)) +
    geom_sf(
      data = shp.world,
      size = 0.1,
      fill = 'grey95'
    ) +
    geom_sf(size = 1.1) +
    geom_sf(data = bind_rows(shp.contours), size = 0.1) +
    geom_sf_label_repel(
      aes(label = segment),
      fill = 'ivory',
      size = 2,
      alpha = 0.8,
      force = 3
    ) +
    coord_sf(expand = F) +
    theme_map(font_size = 8) +
    theme(
      panel.grid = element_line(size = 0.1, color = 'white'),
      panel.background = element_rect(fill = 'grey80', color = NA),
      plot.margin = margin()
  )
# Save 
cat('\nSaving plot to: figs/base/segments.png')
suppressWarnings(suppressMessages(
  ggsave(
    file = 'figs/base/segments.png',
    plot = p1,
    device = 'png',
    type = 'cairo',
    width = 6,
    height = 3
  )
))
system('open figs/base/segments.png', wait = F)

# Global filtered ThermoGlobe
p2 <-
  ggplot() +
    geom_sf(
      data = shp.world,
      size = 0.1,
      fill = 'grey95'
    ) +
    geom_sf(
      data = shp.hf.filtered,
      aes(color = hf),
      size = 0.01,
      shape = 20
    ) +
    v.scale.white +
    labs(color = bquote(mWm^-2)) +
    coord_sf(expand = F) +
    theme_map(font_size = 8) +
    theme(
      axis.text = element_blank(),
      legend.position = c(1, 0),
      legend.justification = c(1, 0),
      legend.direction = 'horizontal',
      legend.box.background = element_rect(fill = rgb(1, 1, 1, 0.7), color = NA),
      legend.box.margin = margin(1, 8, 1, 2),
      legend.key.height = unit(0.1, 'in'),
      legend.key.width = unit(0.3, 'in'),
      legend.title = element_text(vjust = 1),
      panel.grid = element_line(size = 0.1, color = 'white'),
      panel.background = element_rect(fill = 'grey80', color = NA),
      plot.margin = margin()
    )
# Save
cat('\nSaving plot to: figs/base/ThermoGlobeFiltered.png')
suppressWarnings(suppressMessages(
  ggsave(
    file = 'figs/base/ThermoGlobeFiltered.png',
    plot = p2,
    device = 'png',
    type = 'cairo',
    width = 6,
    height = 3
  )
))
system('open figs/base/ThermoGlobeFiltered.png', wait = F)

# Global similarity estimates
p3 <-
  shp.interp.luca %>%
  ggplot() +
    geom_sf(aes(color = est.sim)) +
    geom_sf(
      data = shp.world,
      size = 0.1,
      fill = NA,
      color = 'white'
    ) +
    labs(color = bquote(mWm^-2)) +
    v.scale.white +
    coord_sf(expand = F) +
    theme_map(font_size = 8) +
    theme(
      axis.text = element_blank(),
      legend.position = c(1, 0),
      legend.justification = c(1, 0),
      legend.direction = 'horizontal',
      legend.box.background = element_rect(fill = rgb(1, 1, 1, 0.7), color = NA),
      legend.box.margin = margin(1, 8, 1, 2),
      legend.key.height = unit(0.1, 'in'),
      legend.key.width = unit(0.3, 'in'),
      legend.title = element_text(vjust = 1),
      panel.grid = element_line(size = 0.1, color = 'white'),
      panel.background = element_rect(fill = 'grey80', color = NA),
      plot.margin = margin()
    )
# Save
cat('\nSaving plot to: figs/base/similarityInterp.png')
suppressWarnings(suppressMessages(
  ggsave(
    file = 'figs/base/similarityInterp.png',
    plot = p3,
    device = 'png',
    type = 'cairo',
    width = 6,
    height = 3
  )
))
system('open figs/base/similarityInterp.png', wait = F)

# Global composition
p4 <-
  (p2 + theme(legend.position = 'none', plot.margin = margin(0, 0, 2, 0))) /
  (p3 + theme(plot.margin = margin(2, 0, 0, 0)))
# Save
cat('\nSaving plot to: figs/base/ThermoGlobeComp.png')
suppressWarnings(suppressMessages(
  ggsave(
    file = 'figs/base/ThermoGlobeComp.png',
    plot = p4,
    device = 'png',
    type = 'cairo',
    width = 6,
    height = 6
  )
))
system('open figs/base/ThermoGlobeComp.png', wait = F)

# Global ThermoGlobe buffer
p5 <-
  ggplot() +
    geom_sf(
      data = shp.world,
      size = 0.1,
      fill = 'grey95'
    ) +
    geom_sf(data = bind_rows(shp.buffer), size = 0.1, fill = NA) +
    geom_sf(data = bind_rows(shp.contours), size = 0.05) +
    geom_sf(data = bind_rows(shp.segs), size = 0.5) +
    geom_sf(
      data = bind_rows(shp.hf.crop),
      aes(color = hf),
      size = 0.01,
      shape = 20
    ) +
    v.scale.white +
    labs(color = bquote(mWm^-2)) +
    coord_sf(expand = F) +
    theme_map(font_size = 8) +
    theme(
      axis.text = element_blank(),
      legend.position = c(1, 0),
      legend.justification = c(1, 0),
      legend.direction = 'horizontal',
      legend.box.background = element_rect(fill = rgb(1, 1, 1, 0.7), color = NA),
      legend.box.margin = margin(1, 8, 1, 2),
      legend.key.height = unit(0.1, 'in'),
      legend.key.width = unit(0.3, 'in'),
      legend.title = element_text(vjust = 1),
      panel.grid = element_line(size = 0.1, color = 'white'),
      panel.background = element_rect(fill = 'grey80', color = NA),
      plot.margin = margin()
    )
# Save
cat('\nSaving plot to: figs/base/ThermoGlobeBuffer.png')
suppressWarnings(suppressMessages(
  ggsave(
    file = 'figs/base/ThermoGlobeBuffer.png',
    plot = p5,
    device = 'png',
    type = 'cairo',
    width = 6,
    height = 3
  )
))
system('open figs/base/ThermoGlobeBuffer.png', wait = F)

# Crop countries and similarity to buffer
shp.world.buf <- suppressWarnings(shp.world %>% st_intersection(bind_rows(shp.buffer)))
shp.sim <- suppressWarnings(shp.interp.luca %>% st_intersection(bind_rows(shp.buffer)))
p6 <-
  ggplot() +
    geom_sf(
      data = shp.world,
      size = 0.1,
      fill = 'grey95'
    ) +
    geom_sf(data = shp.sim, aes(color = est.sim), size = 0.1, shape = 15) +
    geom_sf(data = shp.world.buf, size = 0.1, fill = 'grey95', alpha = 0.1) +
    geom_sf(data = bind_rows(shp.contours), size = 0.05) +
    geom_sf(data = bind_rows(shp.segs), size = 0.5, color = 'white') +
    labs(color = bquote(mWm^-2)) +
    v.scale.white +
    coord_sf(expand = F) +
    theme_map(font_size = 8) +
    theme(
      axis.text = element_blank(),
      legend.position = c(1, 0),
      legend.justification = c(1, 0),
      legend.direction = 'horizontal',
      legend.box.background = element_rect(fill = rgb(1, 1, 1, 0.7), color = NA),
      legend.box.margin = margin(1, 8, 1, 2),
      legend.key.height = unit(0.1, 'in'),
      legend.key.width = unit(0.3, 'in'),
      legend.title = element_text(vjust = 1),
      panel.grid = element_line(size = 0.1, color = 'white'),
      panel.background = element_rect(fill = 'grey80', color = NA),
      plot.margin = margin()
    )
# Save
cat('\nSaving plot to: figs/base/similarityBuffer.png')
suppressWarnings(suppressMessages(
  ggsave(
    file = 'figs/base/similarityBuffer.png',
    plot = p6,
    device = 'png',
    type = 'cairo',
    width = 6,
    height = 3
  )
))
system('open figs/base/similarityBuffer.png', wait = F)

# Global buffer composition
p7 <-
  (p5 + theme(legend.position = 'none', plot.margin = margin(0, 0, 2, 0))) /
  (p6 + theme(plot.margin = margin(2, 0, 0, 0)))
# Save
cat('\nSaving plot to: figs/base/ThermoGlobeBufferComp.png')
suppressWarnings(suppressMessages(
  ggsave(
    file = 'figs/base/ThermoGlobeBufferComp.png',
    plot = p7,
    device = 'png',
    type = 'cairo',
    width = 6,
    height = 6
  )
))
system('open figs/base/ThermoGlobeBufferComp.png', wait = F)

# Individual Segments
seg.names %>%
walk(~{
  cat('\nPlotting', .x)
  # Define map parts
  cnt <- shp.contours[[.x]] # Contour
  seg <- shp.segs[[.x]] # Segment
  buf <- shp.buffer[[.x]] # Buffer
  bbx <- st_bbox(buf) %>% bbox_widen(borders = c(0.05, 0.05, 0.05, 0.05)) # Bounding box
  world <- suppressWarnings(shp.world %>% st_crop(bbx)) # Countries
  world.buf <- suppressWarnings(world %>% st_intersection(buf)) # Contries within buffer
  hf <- suppressWarnings(shp.hf.filtered %>% st_crop(bbx)) # Heat flow
  sim <- suppressWarnings(shp.interp.luca %>% st_intersection(buf)) # Similarity interp
  # Define map scale 1:50,000
  wdth <- (st_bbox(buf)$xmax - st_bbox(buf)$xmin)/5e4
  hght <- (st_bbox(buf)$ymax - st_bbox(buf)$ymin)/5e4
  pnt.size <- wdth*hght/2.3e3
  annt.txt.size <- wdth/2.2e1
  base.txt.size <- wdth/1.5e1
  # ThermoGlobe data
  pp1 <- 
    ggplot() +
      geom_sf(data = world, size = 0.1, fill = 'grey95') +
      geom_sf(data = buf, size = 0.3, fill = NA) +
      geom_sf(data = cnt, size = 0.3) +
      geom_sf(data = seg, size = 1.5) +
      geom_sf(data = hf, aes(color = hf), shape = 15, size = pnt.size*0.3) +
      annotate(
        'label',
        x = -Inf,
        y = -Inf,
        label = .x,
        hjust = 0,
        vjust = 0,
        size = annt.txt.size,
        fill = rgb(1, 1, 1, 0.7),
        label.padding = unit(0.15, 'lines'),
        label.r = unit(0.05, 'lines')
      ) +
      v.scale.white +
      labs(color = bquote(mWm^-2)) +
      coord_sf(expand = F) +
      theme_map(font_size = base.txt.size) +
      theme(
        axis.text = element_blank(),
        legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.direction = 'horizontal',
        legend.box.background = element_rect(fill = rgb(1, 1, 1, 0.7), color = NA),
        legend.box.margin = margin(1, 8, 1, 2),
        legend.key.height = unit(hght/30, 'mm'),
        legend.key.width = unit(wdth/10, 'mm'),
        legend.title = element_text(vjust = 1),
        panel.grid = element_line(size = 0.1, color = 'white'),
        panel.background = element_rect(fill = 'grey80', color = NA),
        plot.margin = margin()
      )
  # Save
  cat('\nSaving plot to: figs/base/', str_replace_all(.x, ' ', ''), 'Base.png', sep = '')
  suppressWarnings(suppressMessages(
    ggsave(
      file = paste0('figs/base/', str_replace_all(.x, ' ', ''), 'Base.png'),
      plot = pp1,
      device = 'png',
      type = 'cairo',
      width = wdth,
      height = hght,
      units = 'mm'
    )
  ))
  system(paste0('open figs/base/', str_replace_all(.x, ' ', ''), 'Base.png'), wait = F)
  # Similarity interpolation
  pp2 <- 
    ggplot() +
      geom_sf(data = world, size = 0.1, fill = 'grey95') +
      geom_sf(data = sim, aes(color = est.sim), size = pnt.size, shape = 15) +
      geom_sf(data = world.buf, size = 0.1, fill = 'grey95', alpha = 0.3) +
      geom_sf(data = buf, size = 0.3, fill = NA) +
      geom_sf(data = cnt, size = 0.3, color = 'white') +
      geom_sf(data = seg, size = 1.5, color = 'white') +
      annotate(
        'label',
        x = -Inf,
        y = -Inf,
        label = .x,
        hjust = 0,
        vjust = 0,
        size = annt.txt.size,
        fill = rgb(1, 1, 1, 0.7),
        label.padding = unit(0.15, 'lines'),
        label.r = unit(0.05, 'lines')
      ) +
      v.scale.white +
      labs(color = bquote(mWm^-2)) +
      coord_sf(expand = F) +
      theme_map(font_size = base.txt.size) +
      theme(
        axis.text = element_blank(),
        legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.direction = 'horizontal',
        legend.box.background = element_rect(fill = rgb(1, 1, 1, 0.7), color = NA),
        legend.box.margin = margin(1, 8, 1, 2),
        legend.key.height = unit(hght/30, 'mm'),
        legend.key.width = unit(wdth/10, 'mm'),
        legend.title = element_text(vjust = 1),
        panel.grid = element_line(size = 0.1, color = 'white'),
        panel.background = element_rect(fill = 'grey80', color = NA),
        plot.margin = margin()
      )
  # Save
  cat('\nSaving plot to: figs/base/', str_replace_all(.x, ' ', ''), 'Similarity.png', sep = '')
  suppressWarnings(suppressMessages(
    ggsave(
      file = paste0('figs/base/', str_replace_all(.x, ' ', ''), 'Similarity.png'),
      plot = pp2,
      device = 'png',
      type = 'cairo',
      width = wdth,
      height = hght,
      units = 'mm'
    )
  ))
  system(paste0('open figs/base/', str_replace_all(.x, ' ', ''), 'Similarity.png'), wait = F)
  # ThermoGlobe data without annotation
  pp3 <- 
    ggplot() +
      geom_sf(data = world, size = 0.1, fill = 'grey95') +
      geom_sf(data = buf, size = 0.3, fill = NA) +
      geom_sf(data = cnt, size = 0.3) +
      geom_sf(data = seg, size = 1.5) +
      geom_sf(data = hf, aes(color = hf), shape = 15, size = pnt.size*0.15, show.legend = F) +
      v.scale.white +
      labs(color = bquote(mWm^-2)) +
      coord_sf(expand = F) +
      theme_map(font_size = base.txt.size) +
      theme(
        axis.text = element_blank(),
        legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.direction = 'horizontal',
        legend.box.background = element_rect(fill = rgb(1, 1, 1, 0.7), color = NA),
        legend.box.margin = margin(1, 8, 1, 2),
        legend.key.height = unit(hght/30, 'mm'),
        legend.key.width = unit(wdth/15, 'mm'),
        legend.title = element_text(vjust = 1),
        panel.grid = element_line(size = 0.1, color = 'white'),
        panel.background = element_rect(fill = 'grey80', color = NA),
        plot.margin = margin()
      )
  if(.x != 'Alaska Aleutians') {
    # Composition
    pp4 <-
      (pp3 + theme(plot.margin = margin(0, 2, 0, 0))) +
      (pp2 + theme(plot.margin = margin(0, 0, 0, 2)))
    # Save
    cat('\nSaving plot to: figs/base/', str_replace_all(.x, ' ', ''), 'Comp.png', sep = '')
    suppressWarnings(suppressMessages(
      ggsave(
        file = paste0('figs/base/', str_replace_all(.x, ' ', ''), 'Comp.png'),
        plot = pp4,
        device = 'png',
        type = 'cairo',
        width = wdth*2,
        height = hght,
        units = 'mm'
      )
    ))
    system(paste0('open figs/base/', str_replace_all(.x, ' ', ''), 'Comp.png'), wait = F)
  } else {
    # Composition
    pp4 <-
      (pp3 + theme(plot.margin = margin(0, 0, 2, 0))) /
      (pp2 + theme(plot.margin = margin(2, 0, 0, 0)))
    # Save
    cat('\nSaving plot to: figs/base/', str_replace_all(.x, ' ', ''), 'Comp.png', sep = '')
    suppressWarnings(suppressMessages(
      ggsave(
        file = paste0('figs/base/', str_replace_all(.x, ' ', ''), 'Comp.png'),
        plot = pp4,
        device = 'png',
        type = 'cairo',
        width = wdth,
        height = hght*2,
        units = 'mm'
      )
    ))
    system(paste0('open figs/base/', str_replace_all(.x, ' ', ''), 'Comp.png'), wait = F)
  }
})

cat('\n\n', rep('~', 60), sep='')
cat('\nDone!')
cat('\n', rep('~', 60), '\n', sep='')