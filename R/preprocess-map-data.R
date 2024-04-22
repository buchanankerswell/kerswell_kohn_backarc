#!/usr/bin/env Rscript

# Load packages
cat(rep('~', 45), '\n', sep='')
source('R/functions.R')

# Read and process segment feature labels
shp_fts <-
  read_csv('assets/map_data/segment-feature-labels.csv', show_col_types=F) %>%
  st_as_sf(coords=c(2, 1), crs=wgs) %>% reproject_center_pacific()

# Read and process country boundaries
shp_world <- ne_countries(returnclass='sf') %>% reproject_center_pacific()

# Read and process global relief model
cat('\nProcessing global relief model (ETOPO 2022, NOAA) ...')
cat('\nhttps://www.ncei.noaa.gov/products/etopo-global-relief-model')
if (!dir.exists('assets/map_data/relief')) {dir.create('assets/map_data/relief')}
shp_relief_world <- get_world_bathy()

# Read and process plate boundaries
cat('\n\nProcessing UTIG plate boundary data ...')
cat('\nhttp://www-udc.ig.utexas.edu/external/plates/data.htm')
shp_ridge <-
  st_read('assets/map_data/plates/ridge.gmt', crs=wgs, quiet=T) %>%
  filter(map_lgl(geometry, ~nrow(st_coordinates(.x)) > 1)) %>%
  reproject_center_pacific()
shp_trench <-
  st_read('assets/map_data/plates/trench.gmt', crs=wgs, quiet=T) %>%
  reproject_center_pacific()
shp_transform <-
  st_read('assets/map_data/plates/transform.gmt', crs=wgs, quiet=T) %>%
  reproject_center_pacific()

# Read and process volcanoes (Syracuse et al., 2006)
cat('\n\nProcessing Syracuse and Abers (2006) volcano data ...')
cat('\nhttps://doi.org/10.1029/2005GC001045\n')
ct <- c('ddddddddddddccf')
shp_volc <- 
  bind_rows(
    read_table('assets/map_data/segments/volcanoes.txt', col_type=ct, na='NULL'),
    read_table('assets/map_data/segments/volcanoes_nospread.txt', col_type=ct, na='NULL')
  ) %>%
  st_as_sf(coords=c('Lon', 'Lat'), crs=wgs) %>%
  reproject_center_pacific()

# Read and process ThermoGlobe database (Lucazeau, 2019)
shp_tglobe <-
  read_csv('assets/hf_data/tglobe.csv', show_col_types=F) %>%
  select(longitude, latitude, `heat-flow (mW/m2)`, code6) %>%
  rename(obs=`heat-flow (mW/m2)`) %>%
  filter(!is.na(obs)) %>% filter(code6 != 'D') %>% filter(!is.na(longitude)) %>%
  filter(!is.na(latitude)) %>% filter(obs > 0 & obs <= 250) %>%
  st_as_sf(coords=c(1,2), crs=wgs) %>% reproject_center_pacific() %>%
  handle_zerodist_obs() %>% rename(tglobe=geometry) %>% select(obs, tglobe)

# Read and process similarity interpolation (Lucazeau, 2019)
shp_sim <-
  read_delim('assets/hf_data/HFgrid14.csv', delim=';', col_types=c('ddddd')) %>%
  rename(lon=longitude, lat=latitude, est_sim=HF_pred, sigma_sim=sHF_pred, obs_sim=Hf_obs) %>%
  filter(est_sim > 0 & est_sim <= 250) %>% st_as_sf(coords=c(1, 2), crs=wgs) %>%
  reproject_center_pacific() %>% rename(similarity=geometry)

# Extract interpolation grid
shp_grid <- st_geometry(shp_sim)

# Read and process submap transects
cat('\nProcessing Submap data (Lallemand & Heuret, 2017) ...')
cat('\nhttps://submap.gm.umontpellier.fr\n')
suppressWarnings({suppressMessages({
  shp_submap <-
    combine_json_to_df(list.files('assets/map_data/submap', full.names=T)) %>%
    rename_all(tolower) %>%
    filter(zone != 'MED') %>%
    filter(!(trench_name %in% c('Sandwich'))) %>%
    mutate_all(~ifelse(. == -999, NA, .)) %>%
    mutate(short_name=sub('^(\\D+)(\\d)$', '\\10\\2', short_name)) %>%
    mutate(`next`=sub('^(\\D+)(\\d)$', '\\10\\2', `next`)) %>%
    mutate(previous=sub('^(\\D+)(\\d)$', '\\10\\2', previous)) %>%
    arrange(zone, short_name) %>%
    mutate(id=row_number()) %>%
    mutate(transect=sprintf('LINESTRING(%s %s, %s %s)', lon1, lat1, lon2, lat2)) %>%
    st_as_sf(wkt='transect', crs=wgs) %>%
    reproject_center_pacific(break_dateline=F) %>%
    select(id, zone, short_name, trench_name, harc1, a, phi, m56_vc, m56_azvc, m56_vsn,
           alphad, dz, dz1, n, n1, tau, tau1)
})})

# Plot base maps
cat('\nDrawing base maps ...')
plot_tglobe_base()
walk(shp_submap$short_name, plot_transect)

# Clean up environment
rm(list=lsf.str())
rm(ct, prj, wgs)

# Save
cat('\n\nSaving data to: assets/map_data/map-data.RData')
save.image('assets/map_data/map-data.RData')

cat('\npreprocess-map-data.R complete!\n')