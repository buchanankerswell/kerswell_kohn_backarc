#!/usr/bin/env Rscript

# Load functions and libraries
cat(rep('~', 45), '\n', sep='')
cat('Loading packages and functions ...\n\n')
source('R/functions.R')
load('data/hf.Rdata')
load('data/opt.RData')

# Create directory
dir.create('figs/diff', recursive=T, showWarnings=F)

# Visualize
cat('\n', rep('~', 45), sep='')
cat('\nVisualizing ...')

# Individual Segments
solns %>%
group_by(segment) %>%
slice_min(cost) %>%
pwalk(~{
  # Define map parts
  cnt <- shp.contours[[..1]]
  seg <- shp.segs[[..1]]
  buf <- shp.buffer[[..1]]
  bbx <- st_bbox(buf) %>% bbox_widen(c(0.05, 0.05, 0.05, 0.05))
  world <- suppressWarnings(shp.world %>% st_crop(bbx))
  world.buf <- suppressWarnings(world %>% st_intersection(buf))
  volc <- suppressWarnings(shp.volc %>% st_intersection(buf))
  sim <- suppressWarnings(shp.interp.luca %>% st_intersection(buf))
  ridge <- shp.ridge.crop[[..1]]
  trench <- shp.trench.crop[[..1]]
  transform <- shp.transform.crop[[..1]]
  relief <- shp.relief.crop[[..1]]
  dif <- ..6
  fts <- shp.fts[shp.fts$segment == ..1,]
  v.mod <- ..2
  cost <- ..7
  pnt.size <- 1
  annt.txt.size <- 7
  base.txt.size <- 14
  # Define map scale 1:50,000 [meters]
  wdth <- (st_bbox(buf)$xmax - st_bbox(buf)$xmin)/5e4
  hght <- (st_bbox(buf)$ymax - st_bbox(buf)$ymin)/5e4
  aspect <- wdth/hght
  # Similarity interpolation
  pp1 <- 
    ggplot() +
    ggtitle('a) Similarity') +
    geom_sf(data=relief, aes(color=elevation), shape=15, size=0.01) +
    scale_color_etopo(guide='none') +
    new_scale_color() +
    geom_sf(data=dif, aes(color=est.sim), size=pnt.size*1.5, shape=15) +
    geom_sf(data=world.buf, linewidth=0.1, color='grey60', fill=NA) +
    geom_sf(data=ridge, linewidth=0.5) +
    geom_sf(data=trench, linewidth=0.5) +
    geom_sf(data=transform, linewidth=0.5) +
    geom_sf(data=buf, linewidth=0.1, fill=NA, color='grey60') +
    geom_sf(data=seg, linewidth=1, color='white') +
    geom_sf(data=volc, size=pnt.size, color='white', shape=18) +
    geom_sf_label(
      data=fts,
      aes(label=label),
      size=annt.txt.size*0.5,
      fill=rgb(1, 1, 1, 0.9),
      label.padding=unit(0.15, 'lines'),
      label.r=unit(0.05, 'lines')
    ) +
    scale_color_viridis_c(
      option='magma',
      name=bquote(mWm^-2),
      limits=c(0, 250),
      breaks=c(0, 125, 250),
      na.value='transparent',
      guide=guide_colorbar(title.vjust=1, show.limits=T)
    ) +
    coord_sf() +
    theme_map(font_size=base.txt.size) +
    theme(
      axis.text.x =
        element_text(angle=30, hjust=1, vjust=1, margin=margin(5, 0, 0, 0)),
      axis.text.y=element_text(angle=30, hjust=1, margin=margin(0, 5, 0, 0))
    )
  if(
     .x %in% c('Alaska Aleutians', 'Kamchatka Marianas', 'Tonga New Zealand', 'Vanuatu')
   ) {
    pp1 <-
      pp1 +
      scale_x_continuous(breaks=c(130, 140, 150, 160, 170, 180, -170, -160, -150, -140, -130))
  }
  # Krige interpolation
  pp2 <- 
    ggplot() +
    ggtitle(paste0('b) Kriging: ', v.mod)) +
    geom_sf(data=relief, aes(color=elevation), shape=15, size=0.01) +
    scale_color_etopo(guide='none') +
    new_scale_color() +
    geom_sf(data=dif, aes(color=est.krige), size=pnt.size*1.5, shape=15) +
    geom_sf(data=world.buf, linewidth=0.1, color='grey60', fill=NA) +
    geom_sf(data=ridge, linewidth=0.5) +
    geom_sf(data=trench, linewidth=0.5) +
    geom_sf(data=transform, linewidth=0.5) +
    geom_sf(data=buf, linewidth=0.1, fill=NA, color='grey60') +
    geom_sf(data=seg, linewidth=1, color='white') +
    geom_sf(data=volc, size=pnt.size, color='white', shape=18) +
    scale_color_viridis_c(
      option='magma',
      name=bquote(mWm^-2),
      limits=c(0, 250),
      breaks=c(0, 125, 250),
      na.value='transparent',
      guide=guide_colorbar(title.vjust=1, show.limits=T)
    ) +
    coord_sf() +
    theme_map(font_size=base.txt.size) +
    theme(
      axis.text.x =
        element_text(angle=30, hjust=1, vjust=1, margin=margin(5, 0, 0, 0)),
      axis.text.y=element_text(angle=30, hjust=1, margin=margin(0, 5, 0, 0))
    )
    if(
       .x %in% c('Alaska Aleutians', 'Kamchatka Marianas', 'Tonga New Zealand', 'Vanuatu')
     ) {
      pp2 <-
        pp2 +
        scale_x_continuous(
          breaks=c(130, 140, 150, 160, 170, 180, -170, -160, -150, -140, -130)
        )
    }
  # Save with map scale of 1:5,000,000 [meters]
  cat(
    '\nSaving plot to: figs/diff/',
    str_replace_all(..1, ' ', ''),
    'DiffComp.png',
    sep=''
  )
  if(aspect <= 0.95) {
    # Composition
    p <-
      (pp1 +
        theme(
          axis.text.y=element_text(angle=30, hjust=1),
          axis.text.x=element_text(angle=30, hjust=1, vjust=1)
        )
      ) +
      (pp2 +
       theme(
          axis.text.y=element_blank(),
          axis.text.x=element_text(angle=30, hjust=1, vjust=1)
       )
      ) &
      theme(
        plot.margin=margin(1, 1, 1, 1),
        panel.grid=element_line(linewidth=0.01, color='grey60'),
        legend.position='bottom',
        legend.direction='horizontal',
        legend.justification='center',
        legend.box.margin=margin(),
        legend.text=element_text(margin=margin(t=-4)),
        legend.title=element_text(vjust=0, color='black', size=base.txt.size)
      )
    # Save
    suppressWarnings(suppressMessages(
      ggsave(
        file=paste0('figs/diff/', str_replace_all(..1, ' ', ''), 'DiffComp.png'),
        plot=p,
        device='png',
        type='cairo',
        width=6.5,
        height=6.5,
        units='in'
      )
    ))
  } else {
    # Composition
    p <-
      (pp1 +
        theme(
          axis.text.x=element_blank(),
          axis.text.y=element_text(angle=30, hjust=1)
        )
      ) /
      (pp2 +
        theme(
          axis.text.x=element_text(angle=30, hjust=1, vjust=1),
          axis.text.y=element_text(angle=30, hjust=1)
        )
      ) &
      theme(
        plot.margin=margin(1, 1, 1, 1),
        panel.grid=element_line(linewidth=0.01, color='grey60'),
        legend.justification='center',
        legend.box.margin=margin(),
        legend.text=element_text(margin=margin(t=-4)),
        legend.title=element_text(vjust=0, color='black', size=base.txt.size)
      )
    # Save
    suppressWarnings(suppressMessages(
      ggsave(
        file=paste0('figs/diff/', str_replace_all(..1, ' ', ''), 'DiffComp.png'),
        plot=p,
        device='png',
        type='cairo',
        width=6.5,
        height=6.5,
        units='in'
      )
    ))
  }
})

walk(unique(solns$segment), ~{
  # Define mapping scale
  buf <- shp.buffer[[.x]]
  bbx <- st_bbox(buf) %>% bbox_widen(c(0.05, 0.05, 0.05, 0.05))
  # Define map scale 1:50,000 [meters]
  wdth <- (st_bbox(buf)$xmax - st_bbox(buf)$xmin)/5e4
  hght <- (st_bbox(buf)$ymax - st_bbox(buf)$ymin)/5e4
  aspect <- wdth/hght
  plts <-
    solns %>%
    filter(segment == .x) %>%
    mutate(lab=c('a) ', 'b) ', 'c) ', 'd) ')) %>%
    pmap(~{
      # Define map parts
      cnt <- shp.contours[[..1]]
      seg <- shp.segs[[..1]]
      buf <- shp.buffer[[..1]]
      bbx <- st_bbox(buf) %>% bbox_widen(c(0.05, 0.05, 0.05, 0.05))
      world <- suppressWarnings(shp.world %>% st_crop(bbx))
      world.buf <- suppressWarnings(world %>% st_intersection(buf))
      volc <- suppressWarnings(shp.volc %>% st_intersection(buf))
      ridge <- shp.ridge.crop[[..1]]
      trench <- shp.trench.crop[[..1]]
      transform <- shp.transform.crop[[..1]]
      relief <- shp.relief.crop[[..1]]
      dif <- ..6
      v.mod <- ..2
      cost <- ..7
      # Define points and text sizes
      pnt.size <- 1
      annt.txt.size <- 7
      base.txt.size <- 14
      # Krige interpolation
      p <-
        ggplot() +
        ggtitle(paste0(..8, 'Model: ', v.mod)) +
        geom_sf(data=relief, aes(color=elevation), shape=15, size=0.01) +
        scale_color_etopo(guide='none') +
        new_scale_color() +
        geom_sf(data=dif, aes(color=est.krige), size=pnt.size*1.7, shape=15) +
        geom_sf(data=world.buf, linewidth=0.1, color='grey60', fill=NA) +
        geom_sf(data=ridge, linewidth=0.5) +
        geom_sf(data=trench, linewidth=0.5) +
        geom_sf(data=transform, linewidth=0.5) +
        geom_sf(data=buf, linewidth=0.1, fill=NA, color='grey60') +
        geom_sf(data=seg, linewidth=1, color='white') +
        geom_sf(data=volc, size=pnt.size, color='white', shape=18) +
        scale_color_viridis_c(
          option='magma',
          name=bquote(mWm^-2),
          limits=c(0, 250),
          breaks=c(0, 125, 250),
          na.value='transparent',
          guide=guide_colorbar(title.vjust=1, show.limits=T)
        ) +
        coord_sf() +
        theme_map(font_size=base.txt.size) +
        theme(
          plot.margin=margin(1, 1, 1, 1),
          panel.grid=element_line(linewidth=0.01, color='grey60'),
          legend.justification='center',
          legend.box.margin=margin(),
          legend.text=element_text(margin=margin(t=-4)),
          legend.title=element_text(vjust=0, color='black', size=base.txt.size)
        )
        if(
           ..1 %in% c('Alaska Aleutians', 'Kamchatka Marianas', 'Tonga New Zealand', 'Vanuatu')
         ) {
          p <-
            p +
            scale_x_continuous(
              breaks=c(130, 140, 150, 160, 170, 180, -170, -160, -150, -140, -130)
            )
        }
        p
    })
  # Save
  cat(
    '\nSaving plot to: figs/diff/',
    str_replace_all(.x, ' ', ''),
    'ModelComparison.png',
    sep=''
  )
  p <-
    (
      (plts[[1]] +
        theme(
          axis.text.x=element_blank(),
          axis.text.y=element_text(angle=30, hjust=1, margin=margin(0, 5, 0, 0))
        )
      ) +
      (plts[[2]] + theme(axis.text=element_blank())) +
      (plts[[3]] +
        theme(
          axis.text.x =
            element_text(angle=30, hjust=1, vjust=1, margin=margin(5, 0, 0, 0)),
          axis.text.y=element_text(angle=30, hjust=1, margin=margin(0, 5, 0, 0))
        )
      ) +
      (plts[[4]] +
        theme(
          axis.text.x =
            element_text(angle=30, hjust=1, vjust=1, margin=margin(5, 0, 0, 0)),
          axis.text.y=element_blank()
        )
      )
    )
  suppressWarnings(suppressMessages(
    ggsave(
      file=paste0('figs/diff/', str_replace_all(.x, ' ', ''), 'ModelComparison.png'),
      plot=p,
      device='png',
      type='cairo',
      width=9.75,
      height=9.75,
      units='in'
    )
  ))
})

cat('\ninterpolation-plots.R complete!\n\n')
