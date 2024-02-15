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

# Read ThermoGlobe database from Lucazeau (2019)
cat('\nReading ThermoGlobe data from Lucazeau (2019) ...')
hf <- read_csv('assets/hf_data/tglobe.csv', col_type = cols()) %>%
  select(longitude, latitude, `heat-flow (mW/m2)`, code6) %>% rename(hf = `heat-flow (mW/m2)`)

# Filter bad quality data and missing heat flow
cat(
  '\n Removing', nrow(hf[is.na(hf$hf),]), 'NA heat flow values',
  '\n Removing', nrow(hf[hf$code6 == 'D',]), 'poor quality data (code6 = D)',
  '\n Removing', nrow(hf[is.na(hf$latitude) | is.na(hf$longitude),]), 'obs without coords'
)
shp.hf.filtered <- hf %>% filter(!is.na(hf)) %>% filter(code6 != 'D') %>%
  filter(!is.na(longitude)) %>% filter(!is.na(latitude)) %>% filter(hf > 0 & hf <= 250) %>%
  st_as_sf(coords = c(1,2), crs = proj4.wgs) %>% st_transform(proj4.rp)

# Check for duplicate measurements and remove
dup <- sp::zerodist(sf::as_Spatial(shp.hf.filtered))
n.dup <- nrow(dup)
rid <- vector('integer', nrow(dup))

# For duplicate measurements, select the best quality data point,
# if the quality is the same, randomly select one measurement
for(i in 1:nrow(dup)) {
  if(
    shp.hf.filtered$code6[dup[i,1]] != shp.hf.filtered$code6[dup[i,2]] &
    shp.hf.filtered$code6[dup[i,1]] > shp.hf.filtered$code6[dup[i,2]]
  ) {
    rid[i] <- dup[i,1]
  } else if(
    shp.hf.filtered$code6[dup[i,1]] != shp.hf.filtered$code6[dup[i,2]] &
    shp.hf.filtered$code6[dup[i,1]] < shp.hf.filtered$code6[dup[i,2]]
  ) {
    rid[i] <- dup[i,2]
  } else {
    rid[i] <- dup[i,sample(1:2, 1)]
  }
}
shp.hf.filtered <- shp.hf.filtered %>% slice(-rid) %>% select(hf, geometry) %>% rename(hf = hf)
cat('\n Parsed', length(rid), 'duplicates')
cat('\n Filtered ThermoGlobe dataset contains', nrow(shp.hf.filtered), 'observations')

# Load map data
shp.buffer <- extractor_RData('assets/map_data/preprocessed-map-data.RData', 'shp.buffer')
seg.names <- extractor_RData('assets/map_data/preprocessed-map-data.RData', 'seg.names')

# Crop ThermoGlobe data
cat('\nCropping ThermoGlobe data to buffers ...')
shp.hf.crop <- suppressWarnings({
  shp.buffer %>% map(~shp.hf.filtered %>% st_intersection(.x) %>% select(-segment))
})

# Filtered and cropped hf tibble by segment
hf.crop <- shp.hf.crop %>% map_df(~st_set_geometry(.x, NULL), .id = 'segment')

# Read Lucazeau (2019) interpolation
cat('\nReading similarity interpolation from Lucazeau (2019) ...')
interp.luca <- read_delim('assets/hf_data/HFgrid14.csv', delim = ';', col_types = c('ddddd'))
interp.luca <- interp.luca %>%
  rename(lon = longitude, lat = latitude, est.sim = HF_pred, sigma.sim = sHF_pred,
         obs.sim = Hf_obs) %>% filter(est.sim > 0 & est.sim <= 250)
shp.interp.luca <- interp.luca %>% st_as_sf(coords = c(1,2), crs = proj4.wgs) %>%
  st_transform(proj4.rp)

# Extract Lucazeau (2019) 0.5˚ x 0.5˚ grid
shp.grid <- st_geometry(shp.interp.luca)

# Crop Lucazeau (2019) grids
cat('\nCropping interpolation grids to buffers ...')
shp.grid.crop <- shp.buffer %>% map(~shp.grid %>% st_intersection(.x))

# Print info
cat('\n', rep('-', 45), sep = '')
cat('\nNumber of interpolation points by segment:')
walk2(shp.grid.crop, seg.names, ~{cat('\n', ..2, ':', length(..1))})
cat('\n', rep('-', 45), sep = '')

# Calculate regional Similarity RMSE
rmse.luca <- tibble(segment = seg.names,
                    rmse = map_dbl(seg.names,
                                   ~suppressWarnings(itp_rmse(.x, shp.interp.luca, 'sim'))))

# Save global similarity dataset
fname <- 'assets/hf_data/global-similarity-interp.csv'
st_transform(shp.interp.luca, proj4.wgs) %>%
  st_write(fname, layer_options = "GEOMETRY=AS_XY", quiet = T)
read_csv(fname, show_col_types = F) %>% mutate(X = round(X, 3), Y = round(Y, 3)) %>%
  rename(lon = X, lat = Y) %>% write_csv(fname)

# Clean up environment
cat('\nCleaning up environment ...')
rm(list = lsf.str())
rm(proj4.rp, proj4.wgs, i, shp.hf.filtered, shp.grid, dup, rid, interp.luca, seg.names,
   shp.buffer, n.dup, fname)

# Save
cat('\nSaving data to: assets/hf_data/preprocessed-hf-data.RData')
save.image('assets/hf_data/preprocessed-hf-data.RData')

cat('\npreprocess-hf-data.R complete!\n')
