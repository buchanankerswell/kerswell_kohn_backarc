#!/usr/bin/env Rscript

# Load packages and functions
cat(rep('~', 45), '\n', sep='')
source('R/functions.R')

# Read segment features
shp_fts <-
  read_csv('assets/map_data/segment-feature-labels.csv', show_col_types=F) %>%
  st_as_sf(coords=c(2, 1), crs=wgs) %>% st_transform(eck3)

# Country boundaries
shp_world <- ne_countries(returnclass='sf') %>% slice_dateline()

# Global relief model
cat('\nReading global relief model (ETOPO 2022, NOAA) ...')
cat('\nhttps://www.ncei.noaa.gov/products/etopo-global-relief-model')
if (!dir.exists('assets/map_data/relief')) {dir.create('assets/map_data/relief')}
#shp_relief_world <- get_world_bathy()

# Plate boundaries
cat('\n\nReading UTIG plate boundary data ...')
cat('\nhttp://www-udc.ig.utexas.edu/external/plates/data.htm')
shp_ridge <- st_read('assets/map_data/plates/ridge.gmt', crs=wgs, quiet=T)
shp_trench <- st_read('assets/map_data/plates/trench.gmt', crs=wgs, quiet=T)
shp_transform <- st_read('assets/map_data/plates/transform.gmt', crs=wgs, quiet=T)

# Filter out linestrings with single points
single_pnt_idx <- shp_ridge$geometry %>% map_lgl(~nrow(st_coordinates(.x)) <= 1)

# Fix dateline wrapping and transform
shp_ridge <- shp_ridge[!single_pnt_idx,] %>% slice_dateline()
shp_trench <- shp_trench %>% slice_dateline()
shp_transform <- shp_transform %>% slice_dateline()

# Read Syracuse et al (2006) volcanoes
cat('\n\nReading Syracuse and Abers (2006) volcano data ...')
cat('\nhttps://doi.org/10.1029/2005GC001045\n')
volc <- read_table('assets/map_data/segments/volcanoes.txt',
                   col_type=c('ddddddddddddccf'), na='NULL')
volc2 <- read_table('assets/map_data/segments/volcanoes_nospread.txt',
                    col_type=c('ddddddddddddccf'), na='NULL')

# Change to simple features object and project
shp_volc <-
  bind_rows(volc, volc2) %>% st_as_sf(coords=c('Lon', 'Lat'), crs=wgs) %>% slice_dateline()

# Read submap transects
cat('\nReading Submap data (Lallemand & Heuret, 2017) ...')
cat('\nhttps://submap.gm.umontpellier.fr')
shp_submap <-
  combine_json_to_df(list.files('assets/map_data/submap', full.names=T)) %>%
  rename_all(tolower) %>%
  filter(!is.na(phi)) %>%
  filter(!(grepl('MED', short_name))) %>%
  mutate_all(~ifelse(. == -999, NA, .)) %>%
  mutate(short_name=sub('^(\\D+)(\\d)$', '\\10\\2', short_name)) %>%
  mutate(`next`=sub('^(\\D+)(\\d)$', '\\10\\2', `next`)) %>%
  mutate(previous=sub('^(\\D+)(\\d)$', '\\10\\2', previous)) %>%
  arrange(zone, short_name) %>%
  mutate(id = row_number()) %>%
  mutate(transect=sprintf('LINESTRING(%s %s, %s %s)', lon1, lat1, lon2, lat2)) %>%
  st_as_sf(wkt='transect', crs=wgs) %>%
  slice_dateline()

# Add buffers, bboxes, plate boundaries, volcanoes, and bathy
source('R/functions.R')
x <-
  shp_submap[13:29,] %>%
  mutate(buffer=st_buffer(transect, 5e5, endCapStyle='ROUND')) %>%
  rowwise() %>%
  mutate(bbox=bbox_widen(st_bbox(buffer))) %>%
  mutate(ridge=crop_feature(shp_ridge, bbox)) %>%
  mutate(trench=crop_feature(shp_trench, bbox)) %>%
  mutate(transform=crop_feature(shp_transform, bbox)) %>%
  mutate(volcano=crop_feature(select(shp_volc, geometry), bbox)) %>%
  mutate(bathy=list(get_seg_bathy(bbox))) %>%
  ungroup()

glimpse(x)

ggplot() +
  geom_sf(data=shp_world) +
  geom_sf(data=shp_trench)

x %>%
  ggplot() +
  geom_sf(aes(geometry=bathy)) +
  scale_color_etopo(guide='none') +
  geom_sf(aes(geometry=ridge)) +
  geom_sf(aes(geometry=bbox), fill=NA) +
  geom_sf(aes(geometry=buffer), fill=NA) +
  geom_sf(aes(geometry=transect)) +
#  geom_sf(aes(geometry=trench), color='blue') +
  geom_sf(aes(geometry=volcano), shape=17, color='red')

walk(x$short_name, ~{
  seg <- x %>% filter(short_name == .x)
  if (!is.null(seg$bathy[[1]])) {
    bathy <- list(geom_sf(data=seg$bathy[[1]], aes(color=elev), shape=15, size=0.01),
                  scale_color_etopo(guide='none'), new_scale_color())
    p <- ggplot(seg) + bathy
  } else {
    p <- ggplot(seg)
  }
  p <-
    p +
    geom_sf(aes(geometry=transect), linewidth=1.5) +
    geom_sf(aes(geometry=buffer), fill=NA) +
    geom_sf(aes(geometry=bbox), fill=NA) +
    geom_sf(aes(geometry=volcano), shape=17, color='red') +
    coord_sf() +
    theme_map(font_size=14) +
    theme(plot.margin=margin(1, 1, 1, 1), legend.position='top', legend.justification='right',
          legend.direction='horizontal', axis.text=element_text(hjust=1),
          legend.margin=margin(-4, 0, -12, 0), legend.box.margin=margin(0, 10, 0, 0),
          legend.key.height=unit(0.125, 'in'), legend.key.width=unit(0.2, 'in'),
          legend.title=element_text(vjust=0, color='black', size=14),
          panel.grid=element_line(linewidth=0.05, color='grey20'),
          plot.title=element_text(vjust=0, margin=margin(0, 0, -10, 0)))
  ggsave(file=paste0('~/Downloads/segs/', seg$short_name, '.png'), plot=p, width=6.5,
         height=6.5)
})

# Clean up environment
rm(list=lsf.str())
rm(eck3, wgs, sliver, volc, volc2, shp_volc, shp_ridge, shp_trench, shp_transform,
   single_pnt_idx, files, shp_box, shp_countries, shp_sliver, shp_contours)

# Save
cat('\nSaving data to: assets/map_data/preprocessed-map-data.RData')
save.image('assets/map_data/preprocessed-map-data.RData')

cat('\npreprocess-map-data.R complete!\n')