#!/usr/bin/env Rscript

# Capture output
sink(file = paste0('data/log-', Sys.Date()), append = T, type = 'output', split = T)

# Load functions and libraries
cat(rep('~', 60), '\n', sep='')
cat('Loading packages and functions ...\n\n')

source('R/functions.R')
load('data/hf.Rdata')
load('data/opt.RData')
cat('\nSaving plots to: figs/diff/')

# Visualize
cat('\n', rep('~', 60), sep='')
cat('\nVisualizing ...')

# Individual Segments
solns %>%
group_by(segment) %>%
slice_min(cost) %>%
pwalk(~{
  # Define map parts
  cnt <- shp.contours[[..1]] # Contour
  seg <- shp.segs[[..1]] # Segment
  buf <- shp.buffer[[..1]] # Buffer
  bbx <- st_bbox(buf) %>% bbox_widen(borders = c(0.05, 0.05, 0.05, 0.05)) # Bounding box
  world <- suppressWarnings(shp.world %>% st_crop(bbx)) # Countries
  world.buf <- suppressWarnings(world %>% st_intersection(buf)) # Contries within buffer
  volc <- suppressWarnings(shp.volc %>% st_intersection(buf)) # Contries within buffer
  sim <- suppressWarnings(shp.interp.luca %>% st_intersection(buf)) # Similarity interp
  ridge <- shp.ridge.crop[[..1]]
  trench <- shp.trench.crop[[..1]]
  transform <- shp.transform.crop[[..1]]
  dif <- ..6
  fts <- shp.fts[shp.fts$segment == ..1,]
  v.mod <- ..2
  cost <- ..7
  pnt.size <- 3
  annt.txt.size <- 7
  base.txt.size <- 14
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
  # Similarity interpolation
  pp1 <- 
    ggplot() +
    geom_sf(data = world, size = 0.1, fill = 'seashell2', color = 'grey60') +
    geom_sf(data = dif, aes(color = est.sim), size = pnt.size, shape = 15) +
    geom_sf(data = world.buf, size = 0.1, color = 'grey60', fill = NA) +
    geom_sf(data = ridge, size = 1.5, color = 'black', alpha = 0.8) +
    geom_sf(data = trench, size = 1.5, color = 'black', alpha = 0.8) +
    geom_sf(data = transform, size = 1.5, color = 'black', alpha = 0.8) +
    geom_sf(data = buf, size = 0.1, fill = NA, color = 'grey60', alpha = 0.1) +
    geom_sf(data = seg, size = 2, color = 'white') +
    geom_sf(data = volc, size = pnt.size*0.5, color = 'gold', shape = 18) +
    geom_sf_label(
      data = fts,
      aes(label = label),
      size = annt.txt.size*0.5,
      fill = rgb(1, 1, 1, 0.9),
      label.padding = unit(0.15, 'lines'),
      label.r = unit(0.05, 'lines')
    ) +
    annotate(
      'text',
      label = 'a',
      x = -Inf,
      y = Inf,
      size = annt.txt.size,
      hjust = 0,
      vjust = 1
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
    theme_map(font_size = base.txt.size) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
      axis.text.y = element_text(angle = 30, hjust = 1),
      legend.position = c(1, 0),
      legend.justification = c(1, 0),
      legend.direction = 'horizontal',
      legend.box.background = element_rect(fill = rgb(1, 1, 1, 0.9)),
      legend.box.margin = margin(0, 8, 0, 0),
      legend.key.height = unit(0.125, 'in'),
      legend.key.width = unit(0.2, 'in'),
      legend.title = element_text(vjust = 0, color = 'black', size = base.txt.size),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = 'grey60', color = NA),
      plot.margin = margin(1, 1, 1, 1)
    )
  if(
     .x %in% c('Alaska Aleutians', 'Kamchatka Marianas', 'Tonga New Zealand', 'Vanuatu')
   ) {
    pp1 <-
      pp1 +
      scale_x_continuous(breaks = c(130, 140, 150, 160, 170, 180, -170, -160, -150, -140, -130))
  }
  # Krige interpolation
  pp2 <- 
    ggplot() +
    geom_sf(data = world, size = 0.1, fill = 'seashell2', color = 'grey60') +
    geom_sf(data = dif, aes(color = est.krige), size = pnt.size, shape = 15) +
    geom_sf(data = world.buf, size = 0.1, color = 'grey60', fill = NA) +
    geom_sf(data = ridge, size = 1.5, color = 'black', alpha = 0.8) +
    geom_sf(data = trench, size = 1.5, color = 'black', alpha = 0.8) +
    geom_sf(data = transform, size = 1, color = 'black', alpha = 0.8) +
    geom_sf(data = buf, size = 0.1, fill = NA, color = 'grey60', alpha = 0.1) +
    geom_sf(data = seg, size = 2, color = 'white') +
    geom_sf(data = volc, size = pnt.size*0.5, color = 'gold', shape = 18) +
    annotate(
      'label',
      x = Inf,
      y = -Inf,
      label = paste('variogram model:', v.mod),
      hjust = 1,
      vjust = 0,
      fill = rgb(1, 1, 1, 0.9),
      size = annt.txt.size*0.8,
      label.padding = unit(0.15, 'lines'),
      label.r = unit(0.05, 'lines')
    ) +
    annotate(
      'text',
      label = 'b',
      x = -Inf,
      y = Inf,
      size = annt.txt.size,
      hjust = 0,
      vjust = 1
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
    theme_map(font_size = base.txt.size) +
    theme(
      legend.position = 'none',
      panel.background = element_rect(fill = 'grey60', color = NA),
    )
    if(
       .x %in% c('Alaska Aleutians', 'Kamchatka Marianas', 'Tonga New Zealand', 'Vanuatu')
     ) {
      pp2 <-
        pp2 +
        scale_x_continuous(
          breaks = c(130, 140, 150, 160, 170, 180, -170, -160, -150, -140, -130)
        )
    }
  # Save with map scale of 1:5,000,000 [meters]
  cat(
    '\nSaving plot to: figs/diff/',
    str_replace_all(..1, ' ', ''),
    'DiffComp.png',
    sep = ''
  )
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
        legend.box.margin = margin(),
        legend.text = element_text(margin = margin(t = -4)),
        legend.margin = margin(0, 8, 0, 0)
      )
    # Save
    suppressWarnings(suppressMessages(
      ggsave(
        file = paste0('figs/diff/', str_replace_all(..1, ' ', ''), 'DiffComp.png'),
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
        legend.box.margin = margin(),
        legend.text = element_text(margin = margin(t = -4)),
        legend.margin = margin(0, 8, 0, 0)
      )
    # Save
    suppressWarnings(suppressMessages(
      ggsave(
        file = paste0('figs/diff/', str_replace_all(..1, ' ', ''), 'DiffComp.png'),
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

walk(unique(solns$segment), ~{
  # Define mapping scale
  buf <- shp.buffer[[.x]] # Buffer
  bbx <- st_bbox(buf) %>% bbox_widen(borders = c(0.05, 0.05, 0.05, 0.05)) # Bounding box
  # Define map scale 1:50,000 [meters]
  wdth <- (st_bbox(buf)$xmax - st_bbox(buf)$xmin)/5e4
  hght <- (st_bbox(buf)$ymax - st_bbox(buf)$ymin)/5e4
  aspect <- wdth/hght
  if(aspect <= 0.4) {
    const <- 165.1/hght
    p.wdth <- wdth * const * 4
    p.hght <- hght * const
  } else if(aspect > 0.4 & aspect <= 0.95) {
    const <- 82.55/wdth
    p.wdth <- wdth * const * 2
    p.hght <- hght * const * 2
  } else if(aspect >= 2) {
    const <- 165.1/wdth
    p.wdth <- wdth * const
    p.hght <- hght * const * 4
  } else {
    const <- 82.55/hght
    p.wdth <- wdth * const * 2
    p.hght <- hght * const * 2
  }
  plts <-
    solns %>%
    filter(segment == .x) %>%
    pmap(~{
      # Define map parts
      cnt <- shp.contours[[..1]] # Contour
      seg <- shp.segs[[..1]] # Segment
      buf <- shp.buffer[[..1]] # Buffer
      bbx <- st_bbox(buf) %>% bbox_widen(borders = c(0.05, 0.05, 0.05, 0.05)) # Bounding box
      world <- suppressWarnings(shp.world %>% st_crop(bbx)) # Countries
      world.buf <- suppressWarnings(world %>% st_intersection(buf)) # Contries within buffer
      volc <- suppressWarnings(shp.volc %>% st_intersection(buf)) # Contries within buffer
      ridge <- shp.ridge.crop[[..1]]
      trench <- shp.trench.crop[[..1]]
      transform <- shp.transform.crop[[..1]]
      dif <- ..6
      v.mod <- ..2
      cost <- ..7
      # Define points and text sizes
      pnt.size <- 3
      annt.txt.size <- 7
      base.txt.size <- 14
      # Krige interpolation
      p <-
        ggplot() +
        geom_sf(data = world, size = 0.1, fill = 'seashell2', color = 'grey60') +
        geom_sf(data = dif, aes(color = est.krige), size = pnt.size, shape = 15) +
        geom_sf(data = world.buf, size = 0.1, color = 'grey60', fill = NA) +
        geom_sf(data = ridge, size = 1.5, color = 'black', alpha = 0.8) +
        geom_sf(data = trench, size = 1.5, color = 'black', alpha = 0.8) +
        geom_sf(data = transform, size = 1.5, color = 'black', alpha = 0.8) +
        geom_sf(data = buf, size = 0.1, fill = NA, color = 'grey60') +
        geom_sf(data = seg, size = 2, color = 'white') +
        geom_sf(data = volc, size = pnt.size*0.5, color = 'gold', shape = 18) +
        annotate(
          'label',
          x = Inf,
          y = -Inf,
          label = paste('variogram model:', v.mod),
          hjust = 1,
          vjust = 0,
          fill = rgb(1, 1, 1, 0.9),
          size = annt.txt.size*0.8,
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
        theme_map(font_size = base.txt.size) +
        theme(
          legend.position = c(0, 1),
          legend.justification = c(0, 1),
          legend.direction = 'horizontal',
          legend.box.background = element_rect(fill = rgb(1, 1, 1, 0.9)),
          legend.box.margin = margin(1, 10, 1, 1),
          legend.key.height = unit(0.125, 'in'),
          legend.key.width = unit(0.2, 'in'),
          legend.title = element_text(vjust = 0, color = 'black', size = base.txt.size),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = 'grey60', color = NA),
          plot.margin = margin(1, 1, 1, 1)
        )
        if(
           ..1 %in% c('Alaska Aleutians', 'Kamchatka Marianas', 'Tonga New Zealand', 'Vanuatu')
         ) {
          p <-
            p +
            scale_x_continuous(
              breaks = c(130, 140, 150, 160, 170, 180, -170, -160, -150, -140, -130)
            )
        }
        p
    })
  # Save
  cat(
    '\nSaving plot to: figs/diff/',
    str_replace_all(.x, ' ', ''),
    'ModelComparison.png',
    sep = ''
  )
  if(.x == 'Andes') {
    p <-
      (
        (plts[[1]] +
          theme(
            axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
            axis.text.y = element_text(angle = 30, hjust = 1)
          )
        ) |
        (plts[[2]] +
          theme(
            legend.position = 'none',
            axis.text.y = element_blank(),
            axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
          )
        ) |
        (plts[[3]] +
          theme(
            legend.position = 'none',
            axis.text.y = element_blank(),
            axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
          )
        ) |
        (plts[[4]] +
          theme(
            legend.position = 'none',
            axis.text.y = element_blank(),
            axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
          )
        )
      )
  } else if(.x == 'Alaska Aleutians') {
    p <-
      (
        (plts[[1]] +
          theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(angle = 30, hjust = 1)
          )
        ) /
        (plts[[2]] +
          theme(
            legend.position = 'none',
            axis.text.x = element_blank(),
            axis.text.y = element_text(angle = 30, hjust = 1)
          )
        ) /
        (plts[[3]] +
          theme(
            legend.position = 'none',
            axis.text.x = element_blank(),
            axis.text.y = element_text(angle = 30, hjust = 1)
          )
        ) /
        (plts[[4]] +
          theme(
            legend.position = 'none',
            axis.text.y = element_text(angle = 30, hjust = 1),
            axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
          )
        )
      )
  } else {
    p <-
      (
        (plts[[1]] +
          theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(angle = 30, hjust = 1)
          )
        ) +
        (plts[[2]] + theme(legend.position = 'none', axis.text = element_blank())) +
        (plts[[3]] +
          theme(
            legend.position = 'none',
            axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
            axis.text.y = element_text(angle = 30, hjust = 1)
          )
        ) +
        (plts[[4]] +
          theme(
            legend.position = 'none',
            axis.text.y = element_blank(),
            axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
          )
        )
      )
  }
  suppressWarnings(suppressMessages(
    ggsave(
      file = paste0('figs/diff/', str_replace_all(.x, ' ', ''), 'ModelComparison.png'),
      plot = p,
      device = 'png',
      type = 'cairo',
      width = p.wdth,
      height = p.hght,
      units = 'mm'
    )
  ))
})

cat('\ninterpolation-plots.R complete!\n\n')
sink()
