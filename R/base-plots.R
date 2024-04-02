#!/usr/bin/env Rscript

# Load packages and functions
source('R/functions.R')
load('assets/map_data/preprocessed-map-data.RData')
load('assets/hf_data/preprocessed-hf-data.RData')

# Create directory
dir.create('figs/base', recursive=T, showWarnings=F)

# Crop countries and similarity to buffer
shp_world_buf <- suppressWarnings(shp_world %>% st_intersection(bind_rows(shp_buffer)))
shp_sim <- suppressWarnings(shp_interp_sim %>% st_intersection(bind_rows(shp_buffer)))

# Define map parts
plate_boundaries <- list(geom_sf(data=shp_ridge, linewidth=0.3),
                         geom_sf(data=shp_trench, linewidth=0.3),
                         geom_sf(data=shp_transform, linewidth=0.3))
relief <- list(geom_sf(data=shp_relief_world, aes(color=elevation), shape=15, size=0.01),
               scale_color_etopo(guide='none'), new_scale_color())
buffer <- geom_sf(data=st_union((bind_rows(shp_buffer))), linewidth=0.3, fill=NA)
buffer_world <- geom_sf(data=shp_world_buf, linewidth=0.1, fill='grey70', alpha=0.1)
buffer_sim <- geom_sf(data=shp_sim, aes(color=est_sim), size=0.1, shape=15)
hf <- geom_sf(data=bind_rows(shp_hf_crop), aes(color=hf), size=0.3, shape=20)
segs <- geom_sf(data=bind_rows(shp_segs), linewidth=1, color='white')
guide <- guide_colorbar(title.vjust=1, show.limits=T)
viridis_colors <- scale_color_viridis_c(option='magma', name=bquote(mWm^-2),
                                            limits=c(0, 250), breaks=c(0, 125, 250),
                                            na.value='transparent', guide=guide)
viridis_scale <- list(viridis_colors, labs(color=bquote(mWm^-2)))
seg_labs <- c('AA', 'AN', 'CA', 'KM', 'KR', 'LA', 'NP', 'NBS', 'SP', 'SC', 'SBS', 'TNZ', 'VN')
seg_labels <- geom_sf_label_repel(data=bind_rows(shp_segs), aes(label=seg_labs), size=4,
                                  fill=rgb(1, 1, 1, 0.9), label.padding=unit(0.15, 'lines'),
                                  label.r=unit(0, 'lines'), seed=42)
map_theme <- list(coord_sf(), theme_map(font_size=14))

# Global ThermoGlobe buffer
p1 <- ggplot() + ggtitle('a) Thermoglobe observations') + relief + buffer + plate_boundaries +
  hf + segs + viridis_scale + scale_x_continuous(breaks=seq(-180, 180, 60)) + map_theme
p2 <- ggplot() + ggtitle('b) Similarity predictions') + relief + buffer_sim + buffer_world +
  plate_boundaries + segs + seg_labels + viridis_scale + map_theme

# Global buffer composition
p3 <-
  (p1 + theme(axis.text=element_blank())) / p2 &
  theme(plot.margin=margin(1, 1, 1, 1), legend.position='top', legend.justification='right',
        legend.direction='horizontal', axis.text=element_text(hjust=1),
        legend.margin=margin(-4, 0, -12, 0), legend.box.margin=margin(0, 10, 0, 0),
        legend.key.height=unit(0.125, 'in'), legend.key.width=unit(0.2, 'in'),
        legend.title=element_text(vjust=0, color='black', size=14),
        panel.grid=element_line(linewidth=0.05, color='grey20'),
        plot.title=element_text(vjust=0, margin=margin(0, 0, -10, 0)))

cat('\nSaving plot to: figs/base/ThermoGlobeBufferComp.png')
ggsave(file='figs/base/ThermoGlobeBufferComp.png', plot=p3, width=6.5, height=6.5)

# Individual Segments
seg_names %>% walk(~{
  # Crop map parts
  seg <- shp_segs[[.x]]
  buf <- shp_buffer[[.x]]
  grd <- shp_grid_crop[[.x]]
  bbx <- st_bbox(buf) %>% bbox_widen(c(0.05, 0.05, 0.05, 0.05))
  world <- suppressWarnings(shp_world %>% st_crop(bbx))
  world_buf <- suppressWarnings(world %>% st_intersection(buf))
  rlf <- shp_relief_crop[[.x]]
  volc <- suppressWarnings(shp_volc %>% st_intersection(buf))
  hf <- shp_hf_crop[[.x]]
  sim <- suppressWarnings(shp_interp_sim %>% st_intersection(buf))
  ridge <- shp_ridge_crop[[.x]]
  trench <- shp_trench_crop[[.x]]
  transform <- shp_transform_crop[[.x]]
  fts <- shp_fts[shp_fts$segment == .x,]

  # Define map scale 1:50,000 [meters]
  wdth <- (st_bbox(buf)$xmax - st_bbox(buf)$xmin) / 5e4
  hght <- (st_bbox(buf)$ymax - st_bbox(buf)$ymin) / 5e4
  aspect <- wdth / hght
  pnt_size <- 1
  annt_txt_size <- 7
  base_txt_size <- 14

  # Define map parts
  plate_boundaries <- list(geom_sf(data=ridge, linewidth=0.5),
                           geom_sf(data=trench, linewidth=0.5),
                           geom_sf(data=transform, linewidth=0.5))
  relief <- list(geom_sf(data=rlf, aes(color=elevation), shape=15, size=0.01),
                 scale_color_etopo(guide='none'), new_scale_color())
  buffer <- geom_sf(data=buf, linewidth=0.1, fill=NA, color='grey60')
  world_buffer <- geom_sf(data=world_buf, linewidth=0.1, fill=NA, color='grey60')
  segs <- geom_sf(data=seg, linewidth=1, color='white')
  volcs <- geom_sf(data=volc, size=pnt_size, color='white', shape=18)
  sim_interp <- geom_sf(data=sim, aes(color=est_sim), size=pnt_size * 1.5, shape=15)
  hf_obs <- geom_sf(data=hf, aes(color=hf), shape=15, size=pnt_size)
  fts_labels <- geom_sf_label(data=fts, aes(label=label), size=annt_txt_size * 0.5,
                              fill=rgb(1, 1, 1, 0.9), label.padding=unit(0.15, 'lines'),
                              label.r=unit(0.05, 'lines'))
  long_scale <- scale_x_continuous(breaks=c(130, 140, 150, 160, 170, 180, -170, -160, -150,
                                            -140, -130))
  wide_scale <- scale_x_continuous(breaks=c(130, 140, 150, 160, 170, 180, -170, -160, -150,
                                            -140, -130))
  map_theme <- list(coord_sf(), theme_map(font_size=base_txt_size))
  short_lab <- str_replace_all(.x, ' ', '')

  # Interpolation domain
  pp1 <- ggplot() + ggtitle('a) Observations') + relief + plate_boundaries + segs + hf_obs +
    volcs + viridis_scale + map_theme

  if(.x %in% c('Alaska Aleutians', 'Kamchatka Marianas', 'Tonga New Zealand', 'Vanuatu')) {
    pp1 <- pp1 + long_scale
  }

  # Similarity interpolation
  pp2 <- ggplot() + ggtitle('b) Predictions') + relief + sim_interp + world_buffer +
    plate_boundaries + buffer + segs + fts_labels + viridis_scale + map_theme

  if(.x %in% c('Alaska Aleutians', 'Kamchatka Marianas', 'Tonga New Zealand', 'Vanuatu')) {
    pp2 <- pp2 + wide_scale
  }

  # Compositions
  if(aspect <= 0.95) {
    p <-
      (pp1 + theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1,
                                              margin=margin(5, 0, 0, 0)),
                   axis.text.y=element_text(angle=30, hjust=1, margin=margin(0, 5, 0, 0)))) +
      (pp2 + theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1,
                                              margin=margin(5, 0, 0, 0)),
                   axis.text.y=element_blank())) &
      theme(plot.margin=margin(1, 1, 1, 1),
            panel.grid=element_line(linewidth=0.05, color='grey20'), legend.position='bottom',
            legend.direction='horizontal', legend.justification='center',
            legend.box.margin=margin(), legend.text=element_text(margin=margin(t=-4)),
            legend.title=element_text(vjust=0, color='black', size=base_txt_size))

    cat('\nSaving plot to: figs/base/', short_lab, 'Comp.png', sep='')
    ggsave(file=paste0('figs/base/', short_lab, 'Comp.png'), plot=p, width=6.5, height=6.5)
  } else {
    p <-
      (pp1 + theme(axis.text.x=element_blank(),
                   axis.text.y=element_text(angle=30, hjust=1))) /
      (pp2 + theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1),
                   axis.text.y=element_text(angle=30, hjust=1))) &
      theme(plot.margin=margin(1, 1, 1, 1),
            panel.grid=element_line(linewidth=0.05, color='grey20'),
            legend.justification='center', legend.box.margin=margin(),
            legend.text=element_text(margin=margin(t=-4)),
            legend.title=element_text(vjust=0, color='black', size=base_txt_size))

    cat('\nSaving plot to: figs/base/', short_lab, 'Comp.png', sep='')
    ggsave(file=paste0('figs/base/', short_lab, 'Comp.png'), plot=p, width=6.5, height=6.5)
  }
})

cat('\nbase-plots.R complete!\n\n')