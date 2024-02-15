#!/usr/bin/env Rscript

# Load packages and functions
cat(rep('~', 45), sep = '')
cat('\nLoading packages and functions ...\n')
source('R/functions.R')

# Set seed
set.seed(42)

# Define map projections
proj4.wgs <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
proj4.rp <- paste0('+proj=robin +lon_0=-155 +lon_wrap=-155 +x_0=0 +y_0=0 +ellps=WGS84 ',
                   '+datum=WGS84 +units=m +no_defs')

# Read segment features
shp.fts <- read_csv('assets/map_data/segment-feature-labels.csv', show_col_types = F) %>%
  st_as_sf(coords = c(2, 1), crs = proj4.wgs) %>% st_transform(proj4.rp)

# Create sliver around twenty-five degree long
sliver <- list(rbind(c(25.001, 90), c(25, 90), c(25, -90), c(25.001, -90), c(25.001, 90)))
shp.sliver <- st_polygon(x = sliver) %>% st_sfc() %>% st_set_crs(proj4.wgs)

# Country boundaries
shp.countries <- ne_countries(returnclass = 'sf')

# Fix dateline wrapping and transform
suppressWarnings({suppressMessages({
    shp.world <- shp.countries %>% st_difference(shp.sliver) %>% as_tibble() %>% st_as_sf() %>%
      st_transform(proj4.rp)
})})

# Global relief model from NOAA (ETOPO 2022)
# https://www.ncei.noaa.gov/products/etopo-global-relief-model
cat('\nReading global relief model (ETOPO 2022, NOAA) ...')
cat('\nhttps://www.ncei.noaa.gov/products/etopo-global-relief-model')
suppressWarnings({suppressMessages({
    shp.relief.world <- getNOAA.bathy(180, -180, 90, -90, resolution = 15, keep = T,
                                      path = 'assets/map_data/relief/')
    shp.relief.world <- shp.relief.world %>% as.SpatialGridDataFrame() %>% st_as_sf() %>%
      st_transform(proj4.rp) %>% st_make_valid() %>% rename(elevation = layer)
})})

# Plate boundaries from UTIG
# http://www-udc.ig.utexas.edu/external/plates/data.htm
cat('\n\nReading UTIG plate boundary data ...')
cat('\nhttp://www-udc.ig.utexas.edu/external/plates/data.htm')
shp.ridge <- st_read('assets/map_data/plates/ridge.gmt', crs = proj4.wgs, quiet = T)
shp.trench <- st_read('assets/map_data/plates/trench.gmt', crs = proj4.wgs, quiet = T)
shp.transform <- st_read('assets/map_data/plates/transform.gmt', crs = proj4.wgs, quiet = T)

# Filter out linestrings with single points
single.pnt.idx <- shp.ridge$geometry %>% map_lgl(~nrow(st_coordinates(.x)) <= 1)

# Fix dateline wrapping and transform
suppressWarnings({suppressMessages({
    shp.ridge <- shp.ridge[!single.pnt.idx,] %>% st_wrap_dateline() %>%
      st_difference(shp.sliver) %>% as_tibble() %>% st_as_sf() %>% st_transform(proj4.rp)
    shp.trench <- shp.trench %>% st_wrap_dateline() %>% st_difference(shp.sliver) %>%
      as_tibble() %>% st_as_sf() %>% st_transform(proj4.rp)
    shp.transform <- shp.transform %>% st_wrap_dateline() %>% st_difference(shp.sliver) %>%
      as_tibble() %>% st_as_sf() %>% st_transform(proj4.rp)
})})

# Read Syracuse et al (2006) volcanoes
cat('\n\nReading Syracuse and Abers (2006) arc segment and volcano data ...')
cat('\nhttps://doi.org/10.1029/2005GC001045\n')
volc <- read_table('assets/map_data/segments/volcanoes.txt',
                   col_type = c('ddddddddddddccf'), na = 'NULL')
volc2 <- read_table('assets/map_data/segments/volcanoes_nospread.txt',
                    col_type = c('ddddddddddddccf'), na = 'NULL')

# Change to simple features object and project
shp.volc <- volc %>% bind_rows(volc2) %>% st_as_sf(coords = c('Lon', 'Lat')) %>%
  st_set_crs(proj4.wgs) %>% st_transform(proj4.rp)

# Read syracuse et al 2006 Segments
files <- list.files('assets/map_data/segments/gmts', full.names = TRUE)
seg.names <- files %>% str_extract('(?<=gmts\\/)[a-z].*(?=_contours\\.gmt)') %>%
  str_replace_all('_', ' ') %>% str_to_title()

# Read contours and reproject
shp.contours <- files %>% map(~read_latlong(.x, proj4.rp)) %>% set_names(seg.names)

# First contour defines the segment boundaries
shp.segs <- shp.contours %>% map(~.x[1,])

# Draw 1000 km buffer around segment boundaries
buf.dist <- 1e6
shp.buffer <- shp.segs %>% map(~st_buffer(.x, buf.dist, endCapStyle = 'ROUND'))
shp.box <- shp.buffer %>% map(~st_bbox(.x) %>% bbox_widen(crs = proj4.rp))

# Crop UTIG plate boundaries to bounding boxes
shp.ridge.crop <- suppressWarnings({ shp.box %>% map(~shp.ridge %>% st_crop(.x)) })
shp.trench.crop <- suppressWarnings({ shp.box %>% map(~shp.trench %>% st_crop(.x)) })
shp.transform.crop <- suppressWarnings({ shp.box %>% map(~shp.transform %>% st_crop(.x)) })

# Crop NOAA relief model
shp.relief.crop <-
  map(seg.names, ~{
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
        shp.buffer[[.x]] %>%
        st_bbox() %>%
        bbox_widen(borders = c(0.05, 0.05, 0.05, 0.05)) %>%
        st_transform(proj4.wgs) %>%
        st_bbox()
      }
    suppressWarnings({suppressMessages({
        shp.relief <-
          getNOAA.bathy(bbox[3], bbox[1], bbox[2], bbox[4], resolution = 4, keep = T,
                        path = 'assets/map_data/relief/',
                        antimeridian = ifelse(.x %in% c('Alaska Aleutians',
                                                        'Kamchatka Marianas',
                                                        'Tonga New Zealand',
                                                        'Vanuatu'), T, F))
        cat('\nReprojecting ETOPO relief model for', .x,' ...')
        shp.relief <- shp.relief %>% as.SpatialGridDataFrame() %>% st_as_sf() %>%
          st_transform(proj4.rp) %>% st_make_valid() %>% rename(elevation = layer)
    })})
  }) %>%
  set_names(seg.names)

# Clean up environment
cat('\nCleaning up environment ...')
rm(list = lsf.str())
rm(proj4.rp, proj4.wgs, sliver, volc, volc2, single.pnt.idx, buf.dist, files, shp.box,
   shp.countries, shp.sliver)

# Save
cat('\nSaving data to: assets/map_data/preprocessed-map-data.RData')
save.image('assets/map_data/preprocessed-map-data.RData')

cat('\npreprocess-map-data.R complete!\n')