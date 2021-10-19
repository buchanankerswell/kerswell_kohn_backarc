# Load functions and libraries
cat(rep('~', 60), '\n', sep='')
cat('Loading packages and functions ...\n\n')

source('functions.R')
load('data/hf.Rdata')
load('data/opt.RData')
dir.create('figs/diff', showWarnings = F)

# Visualize
cat('\n', rep('~', 60), sep='')
cat('\nVisualizing ...')
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
  volc <- suppressWarnings(shp.volc %>% st_intersection(buf)) # Contries within buffer
  hf <- suppressWarnings(shp.hf.filtered %>% st_crop(bbx)) # Heat flow
  sim <- suppressWarnings(shp.interp.luca %>% st_intersection(buf)) # Similarity interp
  dif <- shp.interp.diff[[.x]] %>% mutate(sigma.diff = sigma.sim - sigma.krige)
  fts <- shp.fts[shp.fts$segment == .x,]
  # Define map scale 1:50,000
  wdth <- (st_bbox(buf)$xmax - st_bbox(buf)$xmin)/5e4
  hght <- (st_bbox(buf)$ymax - st_bbox(buf)$ymin)/5e4
  # Define points and text sizes
  pnt.size <- wdth*hght/2.3e3
  annt.txt.size <- wdth/2.2e1
  base.txt.size <- wdth/1.5e1
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
        fill = rgb(1, 1, 1, 0.7)
      ) +
      annotate(
        'label',
        x = -Inf,
        y = Inf,
        label = 'a',
        hjust = 0,
        vjust = 1,
        size = annt.txt.size*1.5,
        fill = rgb(1, 1, 1, 0.7),
        label.padding = unit(0.15, 'lines'),
        label.r = unit(0.05, 'lines')
      ) +
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
      annotate(
        'label',
        x = Inf,
        y = -Inf,
        label = 'Similarity',
        hjust = 1,
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
        x = -Inf,
        y = Inf,
        label = 'b',
        hjust = 0,
        vjust = 1,
        size = annt.txt.size*1.5,
        fill = rgb(1, 1, 1, 0.7),
        label.padding = unit(0.15, 'lines'),
        label.r = unit(0.05, 'lines')
      ) +
      annotate(
        'label',
        x = Inf,
        y = -Inf,
        label = 'Krige',
        hjust = 1,
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
        x = -Inf,
        y = Inf,
        label = 'c',
        hjust = 0,
        vjust = 1,
        size = annt.txt.size*1.5,
        fill = rgb(1, 1, 1, 0.7),
        label.padding = unit(0.15, 'lines'),
        label.r = unit(0.05, 'lines')
      ) +
      annotate(
        'label',
        x = Inf,
        y = -Inf,
        label = 'Estimate Difference',
        hjust = 1,
        vjust = 0,
        size = annt.txt.size,
        fill = rgb(1, 1, 1, 0.7),
        label.padding = unit(0.15, 'lines'),
        label.r = unit(0.05, 'lines')
      ) +
      br.palette +
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
        x = -Inf,
        y = Inf,
        label = 'd',
        hjust = 0,
        vjust = 1,
        size = annt.txt.size*1.5,
        fill = rgb(1, 1, 1, 0.7),
        label.padding = unit(0.15, 'lines'),
        label.r = unit(0.05, 'lines')
      ) +
      annotate(
        'label',
        x = Inf,
        y = -Inf,
        label = 'Uncertainty Difference',
        hjust = 1,
        vjust = 0,
        size = annt.txt.size,
        fill = rgb(1, 1, 1, 0.7),
        label.padding = unit(0.15, 'lines'),
        label.r = unit(0.05, 'lines')
      ) +
      berlin.palette +
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
  p <-
    (pp1 + pp2) / (pp3 + pp4) &
    theme(plot.margin = margin(2, 2, 2, 2))
  # Save
  cat('\nSaving plot to: figs/diff/', str_replace_all(.x, ' ', ''), 'DiffComp.png', sep = '')
  suppressWarnings(suppressMessages(
    ggsave(
      file = paste0('figs/diff/', str_replace_all(.x, ' ', ''), 'DiffComp.png'),
      plot = p,
      device = 'png',
      type = 'cairo',
      width = wdth*2,
      height = hght*2,
      units = 'mm'
    )
  ))
  system(paste0('open figs/diff/', str_replace_all(.x, ' ', ''), 'DiffComp.png'), wait = F)
})

cat('\n\n', rep('~', 60), sep='')
cat('\nDone!')
cat('\n\n', rep('~', 60), sep='')
