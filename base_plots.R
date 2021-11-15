#!/usr/bin/env Rscript

# Check files
if(dir.exists('figs/base') & length(list.files('figs/base')) == 16) {
  cat('\nBase plots already exist')
  cat('\nPassing ...\n')
  quit()
}

cat(rep('~', 60), '\n', sep='')
cat('Loading packages and functions ...\n\n')

source('functions.R')
load('data/hf.RData')
dir.create('figs/base', showWarnings = F)

# Visualize
cat('\n', rep('~', 60), sep='')
cat('\nVisualizing ...')

# Global ThermoGlobe buffer
p1 <-
  ggplot() +
    geom_sf(
      data = shp.world,
      size = 0.1,
      fill = 'grey60'
    ) +
    geom_sf(data = bind_rows(shp.buffer), size = 0.2, fill = 'ivory', alpha = 0.3) +
    geom_sf(data = shp.ridge, size = 0.4, color = 'black', alpha = 0.8) +
    geom_sf(data = shp.trench, size = 0.4, color = 'black', alpha = 0.8) +
    geom_sf(data = shp.transform, size = 0.4, color = 'black', alpha = 0.8) +
    geom_sf(
      data = bind_rows(shp.hf.crop),
      aes(color = hf),
      size = 0.3,
      shape = 20
    ) +
    geom_sf(data = bind_rows(shp.segs), size = 0.8, color = 'white') +
    annotate(
      'label',
      x = -Inf,
      y = -Inf,
      label = 'ThermoGlobe dataset',
      hjust = 0,
      vjust = 0,
      size = 5,
      fill = 'ivory',
      label.padding = unit(0.15, 'lines'),
      label.r = unit(0.05, 'lines')
    ) +
    v.scale.grey +
    labs(color = bquote(mWm^-2)) +
    coord_sf(expand = F) +
    theme_map(font_size = 12) +
    theme(
      axis.text = element_text(),
      legend.position = c(1, 1),
      legend.justification = c(1, 1),
      legend.direction = 'horizontal',
      legend.box.background = element_rect(fill = rgb(1, 1, 0.941), color = NA),
      legend.box.margin = margin(1, 8, 1, 2),
      legend.key.height = unit(0.1, 'in'),
      legend.key.width = unit(0.3, 'in'),
      legend.title = element_text(vjust = 1),
      panel.grid = element_line(size = 0.1, color = 'white'),
      panel.background = element_rect(fill = 'grey50', color = NA),
      plot.margin = margin()
    )
# Save
cat('\nSaving plot to: figs/base/ThermoGlobeBuffer.png')
suppressWarnings(suppressMessages(
  ggsave(
    file = 'figs/base/ThermoGlobeBuffer.png',
    plot = p1,
    device = 'png',
    type = 'cairo',
    width = 6,
    height = 3
  )
))

# Crop countries and similarity to buffer
shp.world.buf <- suppressWarnings(shp.world %>% st_intersection(bind_rows(shp.buffer)))
shp.sim <- suppressWarnings(shp.interp.luca %>% st_intersection(bind_rows(shp.buffer)))
p2 <-
  ggplot() +
    geom_sf(
      data = shp.world,
      size = 0.1,
      fill = 'grey60'
    ) +
    geom_sf(data = shp.sim, aes(color = est.sim), size = 0.1, shape = 15) +
    geom_sf(data = shp.world.buf, size = 0.1, fill = 'grey60', alpha = 0.1) +
    geom_sf(data = shp.ridge, size = 0.4, color = 'black', alpha = 0.8) +
    geom_sf(data = shp.trench, size = 0.4, color = 'black', alpha = 0.8) +
    geom_sf(data = shp.transform, size = 0.4, color = 'black', alpha = 0.8) +
    geom_sf(data = bind_rows(shp.segs), size = 0.8, color = 'white') +
    geom_sf_label_repel(
      data = bind_rows(shp.segs),
      aes(label = c('AA', 'AN', 'CA', 'KM', 'KR', 'LA', 'NP', 'NBS', 'SP', 'SC', 'SBS', 'TNZ', 'VN')),
      fill = 'ivory',
      size = 4,
      alpha = 0.9,
      segment.size = 0.3,
      segment.color = 'white',
      segment.curvature = -1e-20,
      arrow = arrow(length = unit(0.015, 'npc')),
      label.padding = unit(0.15, 'lines'),
      label.r = unit(0.05, 'lines')
    ) +
    annotate(
      'label',
      x = -Inf,
      y = -Inf,
      label = 'Similarity interpolation (Lucazeau, 2019)',
      hjust = 0,
      vjust = 0,
      size = 5,
      fill = 'ivory',
      label.padding = unit(0.15, 'lines'),
      label.r = unit(0.05, 'lines')
    ) +
    labs(color = bquote(mWm^-2)) +
    v.scale.grey +
    coord_sf(expand = F) +
    theme_map(font_size = 12) +
    theme(
      axis.text = element_text(),
      legend.position = c(1, 1),
      legend.justification = c(1, 1),
      legend.direction = 'horizontal',
      legend.box.background = element_rect(fill = rgb(1, 1, 0.941), color = NA),
      legend.box.margin = margin(1, 8, 1, 2),
      legend.key.height = unit(0.1, 'in'),
      legend.key.width = unit(0.3, 'in'),
      legend.title = element_text(vjust = 1),
      panel.grid = element_line(size = 0.1, color = 'white'),
      panel.background = element_rect(fill = 'grey50', color = NA),
      plot.margin = margin()
    )
# Save
cat('\nSaving plot to: figs/base/similarityBuffer.png')
suppressWarnings(suppressMessages(
  ggsave(
    file = 'figs/base/similarityBuffer.png',
    plot = p2,
    device = 'png',
    type = 'cairo',
    width = 6,
    height = 3
  )
))

# Global buffer composition
p3 <-
  (p1 +
    theme(
      axis.text = element_blank(),
      plot.margin = margin(0, 0, 2, 0)
    )
  ) /
  (p2 +
    theme(
      legend.position = 'none',
      plot.margin = margin(2, 0, 0, 0)
    )
  ) +
  plot_annotation(tag_level = 'a') &
  theme(plot.tag = element_text(margin = margin(10, 4, 4, 4), color = 'ivory', size = 14))
# Save
cat('\nSaving plot to: figs/base/ThermoGlobeBufferComp.png')
suppressWarnings(suppressMessages(
  ggsave(
    file = 'figs/base/ThermoGlobeBufferComp.png',
    plot = p3,
    device = 'png',
    type = 'cairo',
    width = 6,
    height = 6
  )
))

# Individual Segments
seg.names %>%
walk(~{
  # Define map parts
  cnt <- shp.contours[[.x]] # Contour
  seg <- shp.segs[[.x]] # Segment
  buf <- shp.buffer[[.x]] # Buffer
  bbx <- st_bbox(buf) %>% bbox_widen(borders = c(0.05, 0.05, 0.05, 0.05)) # Bounding box
  world <- suppressWarnings(shp.world %>% st_crop(bbx)) # Countries
  world.buf <- suppressWarnings(world %>% st_intersection(buf)) # Contries within buffer
  volc <- suppressWarnings(shp.volc %>% st_intersection(buf)) # Contries within buffer
  hf <- shp.hf.crop[[.x]]
  sim <- suppressWarnings(shp.interp.luca %>% st_intersection(buf)) # Similarity interp
  ridge <- shp.ridge.crop[[.x]]
  trench <- shp.trench.crop[[.x]]
  transform <- shp.transform.crop[[.x]]
  pnt.size <- 3
  annt.txt.size <- 5
  base.txt.size <- 12
  # Define map scale 1:50,000
  wdth <- (st_bbox(buf)$xmax - st_bbox(buf)$xmin)/5e4
  hght <- (st_bbox(buf)$ymax - st_bbox(buf)$ymin)/5e4
  if(!(.x %in%
      c(
        'Alaska Aleutians',
        'Central America',
        'Kyushu Ryukyu',
        'Scotia',
        'Vanuatu',
        'Lesser Antilles',
        'N Philippines',
        'Sumatra Banda Sea',
        'New Britain Solomon'
      ))
  ) {
    const <- 76.2/wdth
    wdth <- wdth * const
    hght <- hght * const
  } else {
    const <- 101.6/hght
    wdth <- wdth * const
    hght <- hght * const
  }
  # ThermoGlobe data
  pp1 <- 
    ggplot() +
      geom_sf(data = world, size = 0.1, fill = 'grey60') +
      geom_sf(data = ridge, size = 1, color = 'mediumspringgreen', alpha = 0.8) +
      geom_sf(data = trench, size = 1, color = 'mediumspringgreen', alpha = 0.8) +
      geom_sf(data = transform, size = 1, color = 'mediumspringgreen', alpha = 0.8) +
      geom_sf(data = buf, size = 0.3, fill = 'ivory', alpha = 0.1) +
      geom_sf(data = seg, size = 1.5, color = 'white') +
      geom_sf(data = volc, size = pnt.size*0.5, color = 'gold', shape = 18) +
      geom_sf(data = hf, aes(color = hf), shape = 15, size = pnt.size*0.3) +
      annotate(
        'label',
        x = -Inf,
        y = -Inf,
        label = .x,
        hjust = 0,
        vjust = 0,
        size = annt.txt.size,
        fill = rgb(1, 1, 0.941),
        label.padding = unit(0.15, 'lines'),
        label.r = unit(0.05, 'lines')
      ) +
      v.scale.grey +
      labs(color = bquote(mWm^-2)) +
      coord_sf(expand = F) +
      theme_map(font_size = base.txt.size) +
      theme(
        axis.text = element_text(),
        legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.direction = 'horizontal',
        legend.box.background = element_rect(fill = rgb(1, 1, 0.941), color = NA),
        legend.box.margin = margin(1, 8, 1, 2),
        legend.key.height = unit(0.1, 'in'),
        legend.key.width = unit(0.3, 'in'),
        legend.title = element_text(vjust = 1),
        panel.grid = element_line(size = 0.1, color = 'white'),
        panel.background = element_rect(fill = 'grey50', color = NA),
        plot.margin = margin()
      )
  # Similarity interpolation
  pp2 <- 
    ggplot() +
      geom_sf(data = world, size = 0.1, fill = 'grey60') +
      geom_sf(data = sim, aes(color = est.sim), size = pnt.size, shape = 15) +
      geom_sf(data = world.buf, size = 0.1, fill = 'grey60', alpha = 0.1) +
      geom_sf(data = ridge, size = 1, color = 'mediumspringgreen', alpha = 0.8) +
      geom_sf(data = trench, size = 1, color = 'mediumspringgreen', alpha = 0.8) +
      geom_sf(data = transform, size = 1, color = 'mediumspringgreen', alpha = 0.8) +
      geom_sf(data = buf, size = 0.3, fill = NA) +
      geom_sf(data = seg, size = 1.5, color = 'white') +
      annotate(
        'label',
        x = -Inf,
        y = -Inf,
        label = .x,
        hjust = 0,
        vjust = 0,
        size = annt.txt.size,
        fill = rgb(1, 1, 0.941),
        label.padding = unit(0.15, 'lines'),
        label.r = unit(0.05, 'lines')
      ) +
      v.scale.grey +
      labs(color = bquote(mWm^-2)) +
      coord_sf(expand = F) +
      theme_map(font_size = base.txt.size) +
      theme(
        axis.text = element_text(),
        legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.direction = 'horizontal',
        legend.box.background = element_rect(fill = rgb(1, 1, 0.941), color = NA),
        legend.box.margin = margin(1, 8, 1, 2),
        legend.key.height = unit(0.1, 'in'),
        legend.key.width = unit(0.3, 'in'),
        legend.title = element_text(vjust = 1),
        panel.grid = element_line(size = 0.1, color = 'white'),
        panel.background = element_rect(fill = 'grey50', color = NA),
        plot.margin = margin()
      )
  # ThermoGlobe data without annotation
  pp3 <- 
    ggplot() +
      geom_sf(data = world, size = 0.1, fill = 'grey60') +
      geom_sf(data = ridge, size = 1, color = 'mediumspringgreen', alpha = 0.8) +
      geom_sf(data = trench, size = 1, color = 'mediumspringgreen', alpha = 0.8) +
      geom_sf(data = transform, size = 1, color = 'mediumspringgreen', alpha = 0.8) +
      geom_sf(data = buf, size = 0.3, fill = 'ivory', alpha = 0.1) +
      geom_sf(data = seg, size = 1.5, color = 'white') +
      geom_sf(data = volc, size = pnt.size*0.5, color = 'gold', shape = 18) +
      geom_sf(data = hf, aes(color = hf), shape = 15, size = pnt.size*0.3, show.legend = F) +
      v.scale.grey +
      labs(color = bquote(mWm^-2)) +
      coord_sf(expand = F) +
      theme_map(font_size = base.txt.size) +
      theme(
        axis.text = element_text(),
        legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.direction = 'horizontal',
        legend.box.background = element_rect(fill = rgb(1, 1, 0.941), color = NA),
        legend.box.margin = margin(1, 8, 1, 2),
        legend.key.height = unit(0.1, 'in'),
        legend.key.width = unit(0.3, 'in'),
        legend.title = element_text(vjust = 1),
        panel.grid = element_line(size = 0.1, color = 'white'),
        panel.background = element_rect(fill = 'grey50', color = NA),
        plot.margin = margin()
      )
  if(!(.x %in%
      c(
        'Alaska Aleutians',
        'Central America',
        'Kyushu Ryukyu',
        'Scotia',
        'Vanuatu',
        'Lesser Antilles',
        'N Philippines',
        'Sumatra Banda Sea',
        'New Britain Solomon'
      ))
  ) {
    # Composition
    pp4 <-
      (pp3 +
        theme(
          plot.margin = margin(0, 2, 0, 0),
          plot.tag = element_text(color = 'ivory', margin = margin(10, 0, 0, 35), size = 20)
        )
      ) +
      (pp2 +
        theme(
          axis.text.y = element_blank(),
          plot.margin = margin(0, 0, 0, 2),
          plot.tag = element_text(color = 'ivory', margin = margin(10, 0, 0, 4), size = 20)
        )
      ) +
      plot_annotation(tag_level = 'a')
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
  } else {
    # Composition
    pp4 <-
      (pp3 +
        theme(
          axis.text.x = element_blank(),
          plot.margin = margin(0, 0, 2, 0),
          plot.tag = element_text(color = 'ivory', margin = margin(10, 0, 0, 35), size = 20)
        )
      ) /
      (pp2 +
        theme(
          plot.margin = margin(2, 0, 0, 0),
          plot.tag = element_text(color = 'ivory', margin = margin(10, 0, 0, 35), size = 20)
        )
      ) +
      plot_annotation(tag_level = 'a')
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
  }
})

cat('\n\nDone!\n')