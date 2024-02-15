#!/usr/bin/env Rscript

# Load packages and functions
cat(rep('~', 45), '\n', sep = '')
source('R/functions.R')

# Set seed
set.seed(42)

# Define map projections
proj4_wgs <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
proj4_rp <- paste0('+proj=robin +lon_0=-155 +lon_wrap=-155 +x_0=0 +y_0=0 +ellps=WGS84 ',
                   '+datum=WGS84 +units=m +no_defs')

# Read ThermoGlobe database from Lucazeau (2019)
hf <-
  read_csv('assets/hf_data/tglobe.csv', col_type = cols()) %>%
  select(longitude, latitude, `heat-flow (mW/m2)`, code6) %>%
  rename(hf = `heat-flow (mW/m2)`)

# Filter bad quality data and missing heat flow
shp_hf_filtered <-
  hf %>% filter(!is.na(hf)) %>% filter(code6 != 'D') %>% filter(!is.na(longitude)) %>%
  filter(!is.na(latitude)) %>% filter(hf > 0 & hf <= 250) %>%
  st_as_sf(coords = c(1,2), crs = proj4_wgs) %>%
  st_transform(proj4_rp)

# For duplicate measurements, select the best quality data point,
# if the quality is the same, randomly select one measurement
duplicates <- sp::zerodist(sf::as_Spatial(shp_hf_filtered))
n_duplicates <- nrow(duplicates)
rid <- vector('integer', nrow(duplicates))

for(i in 1:nrow(duplicates)) {
  if(
    shp_hf_filtered$code6[duplicates[i,1]] != shp_hf_filtered$code6[duplicates[i,2]] &
    shp_hf_filtered$code6[duplicates[i,1]] > shp_hf_filtered$code6[duplicates[i,2]]
  ) {
    rid[i] <- duplicates[i,1]
  } else if(
    shp_hf_filtered$code6[duplicates[i,1]] != shp_hf_filtered$code6[duplicates[i,2]] &
    shp_hf_filtered$code6[duplicates[i,1]] < shp_hf_filtered$code6[duplicates[i,2]]
  ) {
    rid[i] <- duplicates[i,2]
  } else {
    rid[i] <- duplicates[i,sample(1:2, 1)]
  }
}

shp_hf_filtered <- shp_hf_filtered %>% slice(-rid) %>% select(hf, geometry) %>% rename(hf = hf)

# Load map data
path <- 'assets/map_data/preprocessed-map-data.RData'
shp_buffer <- extract_RData_object(path, 'shp_buffer')
seg_names <- extract_RData_object(path, 'seg_names')

# Crop ThermoGlobe data
shp_hf_crop <- suppressWarnings({
  shp_buffer %>% map(~shp_hf_filtered %>% st_intersection(.x) %>% select(-segment))
})

# Filtered and cropped hf tibble by segment
hf_crop <- shp_hf_crop %>% map_df(~st_set_geometry(.x, NULL), .id = 'segment')

# Read Lucazeau (2019) interpolation
similarity_interpolation <- read_delim('assets/hf_data/HFgrid14.csv', delim = ';',
                                       col_types = c('ddddd'))
shp_similarity_interpolation <-
  similarity_interpolation %>%
  rename(lon = longitude, lat = latitude, est_sim = HF_pred, sigma_sim = sHF_pred,
         obs_sim = Hf_obs) %>%
  filter(est_sim > 0 & est_sim <= 250) %>%
  st_as_sf(coords = c(1,2), crs = proj4_wgs) %>%
  st_transform(proj4_rp)

# Extract interpolation grid
shp_grid <- st_geometry(shp_similarity_interpolation)

# Crop interpolation grid
shp_grid_crop <- shp_buffer %>% map(~shp_grid %>% st_intersection(.x))

# Calculate regional Similarity RMSE
rmse <-
  seg_names %>% map_dbl(~interpolation_rmse(.x, shp_similarity_interpolation, shp_buffer,
                                            shp_grid_crop, shp_hf_crop, 'sim'))
rmse_similarity_interpolation <- tibble(segment = seg_names, rmse = rmse)

# Write global similarity interpolation to csv
path <- 'assets/hf_data/global-similarity-interp.csv'
st_transform(shp_similarity_interpolation, proj4_wgs) %>%
  st_write(path, layer_options = "GEOMETRY=AS_XY", quiet = T)
read_csv(path, show_col_types = F) %>% mutate(X = round(X, 3), Y = round(Y, 3)) %>%
  rename(lon = X, lat = Y) %>% write_csv(path)

# Clean up environment
rm(list = lsf.str())
rm(proj4_rp, proj4_wgs, i, shp_hf_filtered, shp_grid, duplicates, rid,
   similarity_interpolation, seg_names, shp_buffer, n_duplicates, path)

# Save
cat('\nSaving data to: assets/hf_data/preprocessed-hf-data.RData')
save.image('assets/hf_data/preprocessed-hf-data.RData')

cat('\npreprocess-hf-data.R complete!\n')
