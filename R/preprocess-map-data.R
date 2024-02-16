#!/usr/bin/env Rscript

# Load packages and functions
cat(rep('~', 45), '\n', sep='')
source('R/functions.R')

# Define map projections
proj4.wgs <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
proj4.rp <- paste0('+proj=robin +lon_0=-155 +lon_wrap=-155 +x_0=0 +y_0=0 +ellps=WGS84 ',
                   '+datum=WGS84 +units=m +no_defs')

# Read segment features
shp_fts <-
  read_csv('assets/map_data/segment-feature-labels.csv', show_col_types=F) %>%
  st_as_sf(coords=c(2, 1), crs=proj4.wgs) %>% st_transform(proj4.rp)

# Create sliver around twenty-five degree long
sliver <- list(rbind(c(25.001, 90), c(25, 90), c(25, -90), c(25.001, -90), c(25.001, 90)))
shp_sliver <- st_polygon(x=sliver) %>% st_sfc() %>% st_set_crs(proj4.wgs)

# Country boundaries
shp_countries <- ne_countries(returnclass='sf')

# Fix dateline wrapping and transform
suppressWarnings({suppressMessages({
  shp_world <-
    shp_countries %>% st_difference(shp_sliver) %>% as_tibble() %>% st_as_sf() %>%
    st_transform(proj4.rp)
})})

# Global relief model
cat('\nReading global relief model (ETOPO 2022, NOAA) ...')
cat('\nhttps://www.ncei.noaa.gov/products/etopo-global-relief-model')
suppressWarnings({suppressMessages({
  path <- 'assets/map_data/relief/'
  shp_relief_world <-
    getNOAA.bathy(180, -180, 90, -90, resolution=15, keep=T, path=path) %>%
    as.SpatialGridDataFrame() %>% st_as_sf() %>% st_transform(proj4.rp) %>%
    st_make_valid() %>% rename(elevation=layer)
})})

# Plate boundaries
cat('\n\nReading UTIG plate boundary data ...')
cat('\nhttp://www-udc.ig.utexas.edu/external/plates/data.htm')
shp_ridge <- st_read('assets/map_data/plates/ridge.gmt', crs=proj4.wgs, quiet=T)
shp_trench <- st_read('assets/map_data/plates/trench.gmt', crs=proj4.wgs, quiet=T)
shp_transform <- st_read('assets/map_data/plates/transform.gmt', crs=proj4.wgs, quiet=T)

# Filter out linestrings with single points
single_pnt_idx <- shp_ridge$geometry %>% map_lgl(~nrow(st_coordinates(.x)) <= 1)

# Fix dateline wrapping and transform
suppressWarnings({suppressMessages({
  shp_ridge <-
    shp_ridge[!single_pnt_idx,] %>% st_wrap_dateline() %>% st_difference(shp_sliver) %>%
    as_tibble() %>% st_as_sf() %>% st_transform(proj4.rp)
  shp_trench <-
    shp_trench %>% st_wrap_dateline() %>% st_difference(shp_sliver) %>% as_tibble() %>%
    st_as_sf() %>% st_transform(proj4.rp)
  shp_transform <-
    shp_transform %>% st_wrap_dateline() %>% st_difference(shp_sliver) %>% as_tibble() %>%
    st_as_sf() %>% st_transform(proj4.rp)
})})

# Read Syracuse et al (2006) volcanoes
cat('\n\nReading Syracuse and Abers (2006) arc segment and volcano data ...')
cat('\nhttps://doi.org/10.1029/2005GC001045\n')
path <- 'assets/map_data/segments/volcanoes.txt'
volc <- read_table(path, col_type=c('ddddddddddddccf'), na='NULL')
path <- 'assets/map_data/segments/volcanoes_nospread.txt'
volc2 <- read_table(path, col_type=c('ddddddddddddccf'), na='NULL')

# Change to simple features object and project
shp_volc <-
  volc %>% bind_rows(volc2) %>% st_as_sf(coords=c('Lon', 'Lat')) %>%
  st_set_crs(proj4.wgs) %>% st_transform(proj4.rp)

# Read syracuse et al 2006 Segments
files <- list.files('assets/map_data/segments/gmts', full.names=TRUE)
seg_names <-
  files %>% str_extract('(?<=gmts\\/)[a-z].*(?=_contours\\.gmt)') %>%
  str_replace_all('_', ' ') %>% str_to_title()

# Read contours and reproject
shp_contours <- files %>% map(~read_latlong(.x, proj4.rp)) %>% set_names(seg_names)

# First contour defines the segment boundaries
shp_segs <- shp_contours %>% map(~.x[1,])

# Draw 1000 km buffer around segment boundaries
buf_dist <- 1e6
shp_buffer <- shp_segs %>% map(~st_buffer(.x, buf_dist, endCapStyle='ROUND'))
shp_box <- shp_buffer %>% map(~st_bbox(.x) %>% bbox_widen(crs=proj4.rp))

# Crop UTIG plate boundaries to bounding boxes
shp_ridge_crop <- suppressWarnings({shp_box %>% map(~shp_ridge %>% st_crop(.x))})
shp_trench_crop <- suppressWarnings({shp_box %>% map(~shp_trench %>% st_crop(.x))})
shp_transform_crop <- suppressWarnings({shp_box %>% map(~shp_transform %>% st_crop(.x))})

# Crop NOAA relief model
cat('\nCropping relief model ...\n')
shp_relief_crop <-
  map(seg_names, ~{
    bbox <-
      if (.x == 'Kamchatka Marianas') {
        c(-179, -2.1, 120, 69.2)
      } else if (.x == 'Alaska Aleutians') {
        c(-130, 39.38, 150, 71.06)
      } else if (.x == 'Scotia') {
        c(-55.16, -73, 0, -44)
      } else if (.x == 'Tonga New Zealand') {
        c(-160.6, -54.86, 158, -2.94)
      } else {
        shp_buffer[[.x]] %>%
        st_bbox() %>%
        bbox_widen(c(0.05, 0.05, 0.05, 0.05)) %>%
        st_transform(proj4.wgs) %>%
        st_bbox()
      }
    suppressWarnings({suppressMessages({
      path <- 'assets/map_data/relief/'
      shp_relief <-
        getNOAA.bathy(bbox[3], bbox[1], bbox[2], bbox[4], resolution=4, keep=T, path=path,
                      antimeridian=ifelse(.x %in% c('Alaska Aleutians', 'Kamchatka Marianas',
                                                    'Tonga New Zealand', 'Vanuatu'), T, F))
      shp_relief <-
        shp_relief %>% as.SpatialGridDataFrame() %>% st_as_sf() %>% st_transform(proj4.rp) %>%
        st_make_valid() %>% rename(elevation=layer)
    })})
  }) %>%
  set_names(seg_names)

# Clean up environment
rm(list=lsf.str())
rm(proj4.rp, proj4.wgs, sliver, volc, volc2, single_pnt_idx, buf_dist, files, shp_box,
   shp_countries, shp_sliver, path)

# Save
cat('\nSaving data to: assets/map_data/preprocessed-map-data.RData')
save.image('assets/map_data/preprocessed-map-data.RData')

cat('\npreprocess-map-data.R complete!\n')