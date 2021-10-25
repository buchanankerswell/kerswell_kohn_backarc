#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

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

source('functions.R')
load('data/hf.Rdata')
load(paste0('data/opt', cntr, '.RData'))
dir.create('figs/diff', showWarnings = F)
cat('\nSaving plots to: figs/diff/*', cntr, '.png', sep = '')

# Visualize
cat('\n', rep('~', 60), sep='')
cat('\nVisualizing ...')
# Individual Segments
solns %>%
group_by(segment) %>%
filter(v.mod != 'Gau') %>%
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
  dif <- ..6
  fts <- shp.fts[shp.fts$segment == ..1,]
  v.mod <- ..2
  cost <- ..7
  # Define map scale 1:50,000
  wdth <- (st_bbox(buf)$xmax - st_bbox(buf)$xmin)/5e4
  hght <- (st_bbox(buf)$ymax - st_bbox(buf)$ymin)/5e4
  # Define points and text sizes
  pnt.size <- log(wdth*hght, base = 20)
  annt.txt.size <- log(wdth*hght, base = 12)
  base.txt.size <- log(wdth*hght, base = 3)
  # Similarity interpolation
  pp1 <- 
    ggplot() +
      geom_sf(data = world, size = 0.1, fill = 'grey95') +
      geom_sf(data = dif, aes(color = est.sim), size = pnt.size, shape = 15) +
      geom_sf(data = world.buf, size = 0.1, fill = 'grey95', alpha = 0.1) +
      geom_sf(data = buf, size = 0.3, fill = NA) +
      geom_sf(data = cnt, size = 0.3, color = 'white') +
      geom_sf(data = seg, size = 1.5, color = 'white') +
      geom_sf_label(
        data = fts,
        aes(label = label),
        size = annt.txt.size,
        fill = rgb(1, 1, 0.941, 0.7),
        label.padding = unit(0.15, 'lines'),
        label.r = unit(0.05, 'lines')
      ) +
      annotate(
        'label',
        x = -Inf,
        y = -Inf,
        label = ..1,
        hjust = 0,
        vjust = 0,
        size = annt.txt.size,
        fill = rgb(1, 1, 0.941, 0.7),
        label.padding = unit(0.15, 'lines'),
        label.r = unit(0.05, 'lines')
      ) +
      annotate(
        'label',
        x = Inf,
        y = -Inf,
        label = 'Similarity',
        hjust = 1,
        vjust = 0,
        size = annt.txt.size,
        fill = rgb(1, 1, 0.941, 0.7),
        label.padding = unit(0.15, 'lines'),
        label.r = unit(0.05, 'lines')
      ) +
      v.scale.white +
      labs(color = bquote(mWm^-2)) +
      coord_sf(expand = F) +
      theme_map(font_size = base.txt.size) +
      theme(
        axis.text = element_text(color = 'grey20'),
        legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.direction = 'horizontal',
        legend.box.background = element_rect(fill = rgb(1, 1, 0.941, 0.7), color = NA),
        legend.box.margin = margin(1, 8, 1, 2),
        legend.key.height = unit(hght/30, 'mm'),
        legend.key.width = unit(wdth/12, 'mm'),
        legend.title = element_text(vjust = 1),
        panel.grid = element_line(size = 0.1, color = 'white'),
        panel.background = element_rect(fill = 'grey80', color = NA),
        plot.margin = margin()
      )
  # Krige interpolation
  pp2 <- 
    ggplot() +
      geom_sf(data = world, size = 0.1, fill = 'grey95') +
      geom_sf(data = dif, aes(color = est.krige), size = pnt.size, shape = 15) +
      geom_sf(data = world.buf, size = 0.1, fill = 'grey95', alpha = 0.1) +
      geom_sf(data = buf, size = 0.3, fill = NA) +
      geom_sf(data = cnt, size = 0.3, color = 'white') +
      geom_sf(data = seg, size = 1.5, color = 'white') +
      annotate(
        'label',
        x = Inf,
        y = -Inf,
        label = 'Krige',
        hjust = 1,
        vjust = 0,
        size = annt.txt.size,
        fill = rgb(1, 1, 0.941, 0.7),
        label.padding = unit(0.15, 'lines'),
        label.r = unit(0.05, 'lines')
      ) +
      annotate(
        'label',
        x = -Inf,
        y = -Inf,
        label = paste('Model:', v.mod),
        hjust = 0,
        vjust = 0,
        size = annt.txt.size,
        fill = rgb(1, 1, 0.941, 0.7),
        label.padding = unit(0.15, 'lines'),
        label.r = unit(0.05, 'lines')
      ) +
      v.scale.white +
      labs(color = bquote(mWm^-2)) +
      coord_sf(expand = F) +
      theme_map(font_size = base.txt.size) +
      theme(
        axis.text = element_text(color = 'grey20'),
        legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.direction = 'horizontal',
        legend.box.background = element_rect(fill = rgb(1, 1, 0.941, 0.7), color = NA),
        legend.box.margin = margin(1, 8, 1, 2),
        legend.key.height = unit(hght/30, 'mm'),
        legend.key.width = unit(wdth/12, 'mm'),
        legend.title = element_text(vjust = 1),
        panel.grid = element_line(size = 0.1, color = 'white'),
        panel.background = element_rect(fill = 'grey80', color = NA),
        plot.margin = margin()
      )
  # Interpolation difference
  pp3 <- 
    ggplot() +
      geom_sf(data = world, size = 0.1, fill = 'grey95') +
      geom_sf(data = dif, aes(color = est.diff), size = pnt.size, shape = 15) +
      geom_sf(data = world.buf, size = 0.1, fill = 'grey95', alpha = 0.1) +
      geom_sf(data = buf, size = 0.3, fill = NA) +
      geom_sf(data = cnt, size = 0.3, color = 'white') +
      geom_sf(data = seg, size = 1.5, color = 'white') +
      annotate(
        'label',
        x = Inf,
        y = -Inf,
        label = 'Estimate Difference',
        hjust = 1,
        vjust = 0,
        size = annt.txt.size,
        fill = rgb(1, 1, 0.941, 0.7),
        label.padding = unit(0.15, 'lines'),
        label.r = unit(0.05, 'lines')
      ) +
      scale_color_continuous_diverging(
        palette = 'Berlin',
        limits = c(-100, 125),
        na.value = 'grey80'
      ) +
      labs(color = bquote(mWm^-2)) +
      coord_sf(expand = F) +
      theme_map(font_size = base.txt.size) +
      theme(
        axis.text = element_text(color = 'grey20'),
        legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.direction = 'horizontal',
        legend.box.background = element_rect(fill = rgb(1, 1, 0.941, 0.7), color = NA),
        legend.box.margin = margin(1, 8, 1, 2),
        legend.key.height = unit(hght/30, 'mm'),
        legend.key.width = unit(wdth/12, 'mm'),
        legend.title = element_text(vjust = 1),
        panel.grid = element_line(size = 0.1, color = 'white'),
        panel.background = element_rect(fill = 'grey80', color = NA),
        plot.margin = margin()
      )
  # Uncertainty difference
  pp4 <- 
    ggplot() +
      geom_sf(data = world, size = 0.1, fill = 'grey95') +
      geom_sf(data = dif, aes(color = sigma.diff), size = pnt.size, shape = 15) +
      geom_sf(data = world.buf, size = 0.1, fill = 'grey95', alpha = 0.1) +
      geom_sf(data = buf, size = 0.3, fill = NA) +
      geom_sf(data = cnt, size = 0.3, color = 'white') +
      geom_sf(data = seg, size = 1.5, color = 'white') +
      geom_sf(data = volc, shape = 2, size = pnt.size, color = 'deeppink') +
      annotate(
        'label',
        x = Inf,
        y = -Inf,
        label = 'Uncertainty Difference',
        hjust = 1,
        vjust = 0,
        size = annt.txt.size,
        fill = rgb(1, 1, 0.941, 0.7),
        label.padding = unit(0.15, 'lines'),
        label.r = unit(0.05, 'lines')
      ) +
      scale_color_continuous_diverging(
        palette = 'Berlin',
        limits = c(-125, 100),
        na.value = 'grey80'
      ) +
      labs(color = bquote(mWm^-2)) +
      coord_sf(expand = F) +
      theme_map(font_size = base.txt.size) +
      theme(
        axis.text = element_text(color = 'grey20'),
        legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.direction = 'horizontal',
        legend.box.background = element_rect(fill = rgb(1, 1, 0.941, 0.7), color = NA),
        legend.box.margin = margin(1, 8, 1, 2),
        legend.key.height = unit(hght/30, 'mm'),
        legend.key.width = unit(wdth/12, 'mm'),
        legend.title = element_text(vjust = 1),
        panel.grid = element_line(size = 0.1, color = 'white'),
        panel.background = element_rect(fill = 'grey80', color = NA),
        plot.margin = margin()
      )
  # Save
  cat(
    '\nSaving plot to: figs/diff/',
    str_replace_all(..1, ' ', ''),
    'DiffComp',
    cntr,
    '.png',
    sep = ''
  )
  if(..1 == 'Andes') {
    p <-
      (pp1 +
        (pp2 +
          theme(
            legend.position = 'none',
            axis.text.y = element_blank()
          )
        )
      ) +
      (
        (pp3 +
          theme(
            legend.position = 'none',
            axis.text.y = element_blank()
          )
        ) +
        (pp4 +
          theme(
            axis.text.y = element_blank()
          )
        )
      ) &
      theme(plot.margin = margin(2, 2, 2, 2))
    suppressWarnings(suppressMessages(
      ggsave(
        file = paste0('figs/diff/', str_replace_all(..1, ' ', ''), 'DiffComp', cntr, '.png'),
        plot = p,
        device = 'png',
        type = 'cairo',
        width = wdth*4+(wdth*0.06),
        height = hght,
        units = 'mm'
      )
    ))
    # system(
    #   paste0(
    #     'open figs/diff/',
    #     str_replace_all(..1, ' ', ''),
    #     'DiffComp',
    #     cntr,
    #     '.png'
    #   ), wait = F
    # )
  } else {
    p <-
      (
        (pp1 +
          theme(
            axis.text.x = element_blank()
          )
        ) +
        (pp2 +
          theme(
            legend.position = 'none',
            axis.text = element_blank()
          )
        )
      ) /
      (
        (pp3 +
          theme(
            legend.position = 'none'
          )
        ) +
        (pp4 +
          theme(
            axis.text.y = element_blank()
          )
        )
      ) &
      theme(plot.margin = margin(2, 2, 2, 2))
    suppressWarnings(suppressMessages(
      ggsave(
        file = paste0('figs/diff/', str_replace_all(..1, ' ', ''), 'DiffComp', cntr, '.png'),
        plot = p,
        device = 'png',
        type = 'cairo',
        width = wdth*2,
        height = hght*2,
        units = 'mm'
      )
    ))
    # system(
    #   paste0(
    #     'open figs/diff/',
    #     str_replace_all(..1, ' ', ''),
    #     'DiffComp',
    #     cntr,
    #     '.png'
    #     ),
    #   wait = F
    # )
  }
})

walk(unique(solns$segment), ~{
  # Define mapping scale
  buf <- shp.buffer[[..1]] # Buffer
  bbx <- st_bbox(buf) %>% bbox_widen(borders = c(0.05, 0.05, 0.05, 0.05)) # Bounding box
  wdth <- (st_bbox(buf)$xmax - st_bbox(buf)$xmin)/5e4
  hght <- (st_bbox(buf)$ymax - st_bbox(buf)$ymin)/5e4
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
      dif <- ..6
#      fts <- shp.fts[shp.fts$segment == ..1,]
      v.mod <- ..2
      cost <- ..7
      # Define map scale 1:50,000
      wdth <- (st_bbox(buf)$xmax - st_bbox(buf)$xmin)/5e4
      hght <- (st_bbox(buf)$ymax - st_bbox(buf)$ymin)/5e4
      # Define points and text sizes
      pnt.size <- log(wdth*hght, base = 20)
      annt.txt.size <- log(wdth*hght, base = 12)
      base.txt.size <- log(wdth*hght, base = 3)
      # Krige interpolation
      ggplot() +
        geom_sf(data = world, size = 0.1, fill = 'grey95') +
        geom_sf(data = dif, aes(color = est.krige), size = pnt.size, shape = 15) +
        geom_sf(data = world.buf, size = 0.1, fill = 'grey95', alpha = 0.1) +
        geom_sf(data = buf, size = 0.3, fill = NA) +
        geom_sf(data = cnt, size = 0.3, color = 'white') +
        geom_sf(data = seg, size = 1.5, color = 'white') +
        annotate(
          'label',
          x = Inf,
          y = -Inf,
          label = ..1,
          hjust = 1,
          vjust = 0,
          size = annt.txt.size,
          fill = rgb(1, 1, 0.941, 0.7),
          label.padding = unit(0.15, 'lines'),
          label.r = unit(0.05, 'lines')
        ) +
        annotate(
          'label',
          x = -Inf,
          y = -Inf,
          label = paste('Model:', v.mod),
          hjust = 0,
          vjust = 0,
          size = annt.txt.size,
          fill = rgb(1, 1, 0.941, 0.7),
          label.padding = unit(0.15, 'lines'),
          label.r = unit(0.05, 'lines')
        ) +
        v.scale.white +
        labs(color = bquote(mWm^-2)) +
        coord_sf(expand = F) +
        theme_map(font_size = base.txt.size) +
        theme(
          axis.text = element_text(color = 'grey20'),
          legend.position = c(1, 1),
          legend.justification = c(1, 1),
          legend.direction = 'horizontal',
          legend.box.background = element_rect(fill = rgb(1, 1, 0.941, 0.7), color = NA),
          legend.box.margin = margin(1, 8, 1, 2),
          legend.key.height = unit(hght/30, 'mm'),
          legend.key.width = unit(wdth/12, 'mm'),
          legend.title = element_text(vjust = 1),
          panel.grid = element_line(size = 0.1, color = 'white'),
          panel.background = element_rect(fill = 'grey80', color = NA),
          plot.margin = margin(2, 2, 2, 2)
        )
    })
  # Save
  cat(
    '\nSaving plot to: figs/diff/',
    str_replace_all(.x, ' ', ''),
    'ModelComparison',
    cntr,
    '.png',
    sep = ''
  )
  if(.x == 'Andes') {
    p <-
      plts[[1]] +
      (plts[[2]] + theme(legend.position = 'none', axis.text.y = element_blank())) +
      (plts[[3]] + theme(legend.position = 'none', axis.text.y = element_blank())) +
      (plts[[4]] + theme(legend.position = 'none', axis.text.y = element_blank())) +
      (plts[[5]] + theme(legend.position = 'none', axis.text.y = element_blank())) +
      (plts[[6]] + theme(legend.position = 'none', axis.text.y = element_blank())) +
      plot_layout(ncol = 6)
    suppressWarnings(suppressMessages(
      ggsave(
        file =
          paste0('figs/diff/', str_replace_all(.x, ' ', ''), 'ModelComparison', cntr, '.png'),
        plot = p,
        device = 'png',
        type = 'cairo',
        width = wdth*6,
        height = hght,
        units = 'mm'
      )
    ))
    # system(
    #   paste0(
    #     'open figs/diff/',
    #     str_replace_all(.x, ' ', ''),
    #     'ModelComparison',
    #     cntr,
    #     '.png'
    #   ),
    #   wait = F
    # )
  } else {
    p <-
      (
        (plts[[1]] +
          theme(
            axis.text.x = element_blank()
          )
        ) +
        (plts[[2]] +
          theme(
            legend.position = 'none',
            axis.text = element_blank()
          )
        )
      ) /
      (
        (plts[[3]] +
          theme(
            legend.position = 'none',
            axis.text.x = element_blank()
          )
        ) +
        (plts[[4]] +
          theme(
            legend.position = 'none',
            axis.text = element_blank()
          )
        )
      ) /
      (
        (plts[[5]] +
          theme(
            legend.position = 'none'
          )
        ) +
        (plts[[6]] +
          theme(
            legend.position = 'none',
            axis.text.y = element_blank()
          )
        )
      )
    suppressWarnings(suppressMessages(
      ggsave(
        file =
          paste0('figs/diff/', str_replace_all(.x, ' ', ''), 'ModelComparison', cntr, '.png'),
        plot = p,
        device = 'png',
        type = 'cairo',
        width = wdth*3,
        height = hght*3,
        units = 'mm'
      )
    ))
    # system(
    #   paste0(
    #     'open figs/diff/',
    #     str_replace_all(.x, ' ', ''),
    #     'ModelComparison',
    #     cntr,
    #     '.png'
    #   ),
    #   wait = F
    # )
  }
})

cat('\n\nDone!\n')
