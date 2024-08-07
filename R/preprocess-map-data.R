#!/usr/bin/env Rscript

main <- function() {
  source('R/functions.R')
  suppressWarnings({suppressMessages({
    if (!file.exists('assets/map_data/map-data.RData')) {
      # Read and process country boundaries
      shp_world <- ne_countries(returnclass='sf') %>% reproject_center_pacific()
    
      # Read and process global relief model
      cat('\nProcessing global relief model (ETOPO 2022, NOAA) ...')
      cat('\n   https://www.ncei.noaa.gov/products/etopo-global-relief-model')
      if (!dir.exists('assets/map_data/relief')) {dir.create('assets/map_data/relief')}
      shp_relief_world <- get_world_bathy()
    
      # Read and process plate boundaries
      cat('\nProcessing UTIG plate boundary data ...')
      cat('\n   http://www-udc.ig.utexas.edu/external/plates/data.htm')
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
      cat('\nProcessing Syracuse and Abers (2006) volcano data ...')
      cat('\n   https://doi.org/10.1029/2005GC001045')
      ct <- c('ddddddddddddccf')
      volc_path1 <- 'assets/map_data/volcanoes/volcanoes.txt'
      volc_path2 <- 'assets/map_data/volcanoes/volcanoes_nospread.txt'
      shp_volc <- bind_rows(read_table(volc_path1, col_type=ct, na='NULL'),
                            read_table(volc_path2, col_type=ct, na='NULL')) %>%
      st_as_sf(coords=c('Lon', 'Lat'), crs=wgs) %>%
      reproject_center_pacific()
    
      # Read and process IHFC database
      cat('\nProcessing IHFC database 2024 release ...')
      cat('\n   https://doi.org/10.5880/fidgeo.2024.014')
      shp_ghf <-
        read_excel('assets/hf_data/IHFC_2024_GHFDB.xlsx', skip=5) %>%
        select(long_EW, lat_NS, q, Quality_Code) %>% rename(obs=q) %>%
        filter(obs > 0 & obs <= 250) %>% st_as_sf(coords=c(1,2), crs=wgs) %>%
        reproject_center_pacific() %>% parse_zerodist_obs() %>%
        rename(ghf=geometry) %>% select(obs, ghf)
    
      # Read and process similarity interpolation (Lucazeau, 2019)
      cat('\nProcessing similarity interpolation (Lucazeau, 2019) ...')
      cat('\n   https://doi.org/10.1029/2019GC008389')
      shp_sim <-
        read_delim('assets/hf_data/HFgrid14.csv', delim=';', col_types=c('ddddd')) %>%
        rename(lon=longitude, lat=latitude, est_sim=HF_pred, sigma_sim=sHF_pred,
               obs_sim=Hf_obs) %>%
        filter(est_sim > 0 & est_sim <= 250) %>% st_as_sf(coords=c(1, 2), crs=wgs) %>%
        reproject_center_pacific() %>% rename(similarity=geometry)
    
      # Extract interpolation grid
      shp_grid <- st_geometry(shp_sim)
    
      # Read and process submap transects
      cat('\nProcessing Submap data (Lallemand & Heuret, 2017) ...')
      cat('\n   https://submap.gm.umontpellier.fr')
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
          select(id, zone, short_name, trench_name, harc1, a, phi, m56_vc, m56_azvc,
                 m56_vsn, alphad, dz, dz1, n, n1, tau, tau1)
      })})
    
      # Clean up environment
      rm(list=lsf.str())
      rm(ct, volc_path1, volc_path2, prj, wgs, seed)
    
      # Save
      cat('\nSaving data to: assets/map_data/map-data.RData')
      save(shp_ghf, shp_grid, shp_relief_world, shp_ridge, shp_sim, shp_submap,
           shp_transform, shp_trench, shp_volc, shp_world,
           file='assets/map_data/map-data.RData')
    } else {
      cat('\nPreprocessed map data found!')
      plot_ghf()
    }
  })})
  cat('\n', rep('=', 45), sep='')
}

main()