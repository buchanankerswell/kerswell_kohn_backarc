#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

# Test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop('Provide the optimization data file as the first arugment <opt*>.R', call.=FALSE)
} else if (length(args) == 1) {
  if(!file.exists(args[1])) {
    stop('That <opt*.R> file does not exist in data/', call.=FALSE)
  } else {
    cntr <- as.numeric(regmatches(args[1], regexec('[0-9]', args[1])))
    if(is.na(cntr)){
      cntr <- NULL
    }
  }
}

# Load functions and libraries
cat(rep('~', 60), '\n', sep='')
cat('Loading packages and functions ...\n\n')

source('R/functions.R')
load('data/hf.Rdata')
load(paste0('data/opt', cntr, '.RData'))
dir.create('figs/diff', showWarnings = F)
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
  # Similarity interpolation
  pp1 <- 
    ggplot() +
    geom_sf(data = world, size = 0.1, fill = 'grey70', color = 'black') +
    geom_sf(data = dif, aes(color = est.sim), size = pnt.size, shape = 15) +
    geom_sf(data = world.buf, size = 0.1, color = 'grey80', fill = 'grey70', alpha = 0.1) +
    geom_sf(data = ridge, size = 1, color = 'black', alpha = 0.8) +
    geom_sf(data = trench, size = 1, color = 'black', alpha = 0.8) +
    geom_sf(data = transform, size = 1, color = 'black', alpha = 0.8) +
    geom_sf(data = buf, size = 0.3, fill = NA) +
    geom_sf(data = seg, size = 1.5, color = 'white') +
    geom_sf(data = volc, size = pnt.size*0.5, color = 'gold', shape = 18) +
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
      legend.position = c(1, 1),
      legend.justification = c(1, 1),
      legend.direction = 'horizontal',
      legend.box.background = element_rect(color = NA, fill = rgb(1, 1, 1, 0.8)),
      legend.text = element_text(color = 'black'),
      legend.key.height = unit(0.125, 'in'),
      legend.key.width = unit(0.15, 'in'),
      legend.title = element_text(vjust = 1, color = 'black', size = 10)
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
      geom_sf(data = world, size = 0.1, fill = 'grey70', color = 'black') +
      geom_sf(data = dif, aes(color = est.krige), size = pnt.size, shape = 15) +
      geom_sf(data = world.buf, size = 0.1, color = 'grey80', fill = 'grey70', alpha = 0.1) +
      geom_sf(data = ridge, size = 1, color = 'black', alpha = 0.8) +
      geom_sf(data = trench, size = 1, color = 'black', alpha = 0.8) +
      geom_sf(data = transform, size = 1, color = 'black', alpha = 0.8) +
      geom_sf(data = buf, size = 0.3, fill = NA) +
      geom_sf(data = seg, size = 1.5, color = 'white') +
      geom_sf(data = volc, size = pnt.size*0.5, color = 'gold', shape = 18) +
      annotate(
        'label',
        x = -Inf,
        y = -Inf,
        label = paste('model:', v.mod),
        hjust = 0,
        vjust = 0,
        fill = 'grey90',
        size = annt.txt.size,
        label.padding = unit(0.15, 'lines'),
        label.r = unit(0.05, 'lines')
      ) +
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
      theme(legend.position = 'none')
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
        panel.grid = element_line(size = 0.1, color = 'white'),
        panel.background = element_rect(fill = 'grey50', color = NA),
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
        panel.grid = element_line(size = 0.1, color = 'white'),
        panel.background = element_rect(fill = 'grey50', color = NA),
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
      pnt.size <- 2
      annt.txt.size <- 4
      base.txt.size <- 12
      # Krige interpolation
      p <-
        ggplot() +
        geom_sf(data = world, size = 0.1, fill = 'grey70', color = 'black') +
        geom_sf(data = dif, aes(color = est.krige), size = pnt.size, shape = 15) +
        geom_sf(data = world.buf, size = 0.1, color = 'grey80', fill = 'grey70', alpha = 0.1) +
        geom_sf(data = ridge, size = 0.66, color = 'black', alpha = 0.8) +
        geom_sf(data = trench, size = 0.66, color = 'black', alpha = 0.8) +
        geom_sf(data = transform, size = 0.66, color = 'black', alpha = 0.8) +
        geom_sf(data = buf, size = 0.3, fill = NA) +
        geom_sf(data = seg, size = 1, color = 'white') +
        geom_sf(data = volc, size = pnt.size*0.7, color = 'gold', shape = 18) +
        annotate(
          'label',
          x = -Inf,
          y = -Inf,
          label = paste('model:', v.mod),
          hjust = 0,
          vjust = 0,
          size = annt.txt.size,
          fill = 'grey90',
          label.padding = unit(0.15, 'lines'),
          label.r = unit(0.05, 'lines')
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
          axis.text = element_text(color = 'grey20'),
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

cat('\n\nDone!\n')
