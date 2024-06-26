#######################################################
## .0.              Load Libraries               !!! ##
#######################################################
# Quiet loading
sshhh <- function(package_list) {
  suppressWarnings(suppressPackageStartupMessages({
    require(package_list, quietly=T, character.only=T)
  }))
}

# Package list
package_list <- c('tictoc', 'stringr', 'tidyr', 'readr', 'readxl', 'purrr', 'furrr',
                  'tibble', 'dplyr', 'magrittr', 'units', 'ggplot2', 'colorspace', 'metR',
                  'ggrepel', 'ggridges', 'ggnewscale', 'patchwork', 'cowplot', 'ggsflabel',
                  'marmap', 'scales', 'ggspatial', 'gstat', 'rgeos', 'sp', 'sf',
                  'rnaturalearth', 'nloptr', 'zoo', 'jsonlite')

# Load packages quietly
sapply(package_list, sshhh)
rm(package_list, sshhh)

# Turn off S2 geometry
suppressMessages(sf_use_s2(F))

# Set seed
seed <- 42
set.seed(seed)

# Set map projections
wgs <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +no_defs'
prj <- '+proj=eck3 +lon_0=-180 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs'

#######################################################
## .1.         General Helper Functions          !!! ##
#######################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# combine json to df !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
combine_json_to_df <- function(files) {
  map_dfr(files, ~fromJSON(.x))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# bbox extend !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bbox_extend <- function(bbox, ext=rep(0, 4), square=F) {
  b <- bbox
  if (any(ext != 0)) {
    xrange <- b$xmax - b$xmin
    yrange <- b$ymax - b$ymin
    b[1] <- b[1] - (ext[1] * xrange)
    b[3] <- b[3] + (ext[2] * xrange)
    b[2] <- b[2] - (ext[4] * yrange)
    b[4] <- b[4] + (ext[3] * yrange)
  }
  if (square) {
    xrange <- b$xmax - b$xmin
    yrange <- b$ymax - b$ymin
    factor <- xrange / yrange
    if (factor > 1) {
      center_y <- (b[2] + b[4]) / 2
      y_half <- xrange / 2
      b[2] <- center_y - y_half
      b[4] <- center_y + y_half
    } else {
      center_x <- (b[1] + b[3]) / 2
      x_half <- yrange / 2
      b[1] <- center_x - x_half
      b[3] <- center_x + x_half
    }
  }
  box <- c(b$xmin, b$ymax, b$xmin, b$ymin, b$xmax, b$ymin, b$xmax, b$ymax, b$xmin, b$ymax)
  st_polygon(list(matrix(box, ncol=2, byrow=T))) %>% st_sfc(crs=st_crs(bbox))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# reproject center pacific !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
reproject_center_pacific <- function(shp, break_dateline=T, tol=1) {
  if (st_crs(shp) != st_crs(wgs)) {shp <- st_transform(shp, wgs)}
  if (any(st_geometry_type(shp) %in% 'POINT') || !break_dateline) {
    shp %>% st_make_valid() %>% st_transform(prj)
  } else {
    lon0 <- as.numeric(str_extract(prj, '(?<=lon_0=)-?\\d+')) + 360
    suppressWarnings({suppressMessages({
      shp %>% st_make_valid() %>%
        st_wrap_dateline(options=c('WRAPDATELINE=YES', paste0('DATELINEOFFSET=', lon0))) %>%
        st_break_antimeridian(lon_0=lon0, tol=tol) %>% st_transform(prj)
    })})
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# crop feature !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
crop_feature <- function(ft, bbox, crop_within=F, keep_df=F) {
  if (st_crs(ft) != st_crs(bbox)) {ft <- st_transform(ft, st_crs(bbox))}
  if (!crop_within) {
    suppressWarnings({suppressMessages({ft_cropped <- st_crop(ft, bbox)})})
  } else {
    suppressWarnings({suppressMessages({ft_cropped <- st_intersection(ft, bbox)})})
  }
  if (!is.data.frame(ft_cropped)) {l <- length(ft_cropped)} else {l <- nrow(ft_cropped)}
  if (l == 0) {
    if (!keep_df) {NA} else {NULL}
  } else {
    if (!keep_df) {
      suppressWarnings({suppressMessages({st_sfc(st_union(ft_cropped), crs=st_crs(bbox))})})
    } else {ft_cropped}
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get world bathy !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_world_bathy <- function(res=15, path='assets/map_data/relief/') {
  if (!dir.exists(path)) {dir.create(path, recursive=T, showWarnings=F)}
  tryCatch({
    suppressWarnings({suppressMessages({
      getNOAA.bathy(180, -180, 90, -90, res, T, F, path) %>%
        as.SpatialGridDataFrame() %>% st_as_sf(crs=wgs) %>% rename(elev=layer) %>%
        reproject_center_pacific()
    })})
  }, error=function(e) {
    cat('\n!! ERROR occurred in get_world_bathy:\n!!', conditionMessage(e))
    return(NULL)
  })
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get seg bathy !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_seg_bathy <- function(shp, res=2, path='assets/map_data/relief/', tol=1) {
  if (!dir.exists(path)) {dir.create(path, recursive=T, showWarnings=F)}
  tryCatch({
    bbx <- shp %>% st_transform(wgs) %>% st_bbox() %>% round(2)
    if (bbx[3] - bbx[1] > 180) {antim <- T} else {antim <- F}
    getNOAA.bathy(bbx[3], bbx[1], bbx[2], bbx[4], res, T, antim, path) %>%
      as.SpatialGridDataFrame() %>% st_as_sf(crs=wgs) %>% st_make_valid() %>%
      st_transform(prj) %>% rename(elev=layer)
  }, error=function(e) {
    cat('\n!! ERROR occurred in get_seg_bathy:\n!!', conditionMessage(e))
    return(NULL)
  })
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# parse zerodist obs !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
parse_zerodist_obs <- function(df) {
  dup <- zerodist(as_Spatial(df))
  n_dup <- nrow(dup)
  rid <- map_dbl(1:n_dup, ~ {
    i <- .x
    loc1 <- dup[i, 1]
    loc2 <- dup[i, 2]
    locr <- sample(c(loc1, loc2), 1)
    obs1 <- df$obs[loc1]
    obs2 <- df$obs[loc2]
    if (obs1 == obs2) {return(locr)}
    U1 <- str_extract(df$Quality_Code[loc1], '(?<=U).')
    u1 <- if (U1 == 'x') {5} else {as.numeric(U1)}
    U2 <- str_extract(df$Quality_Code[loc2], '(?<=U).')
    u2 <- if (U2 == 'x') {5} else {as.numeric(U2)}
    M1 <- str_extract(df$Quality_Code[loc1], '(?<=M).')
    m1 <- if (M1 == 'x') {5} else {as.numeric(M1)}
    M2 <- str_extract(df$Quality_Code[loc2], '(?<=M).')
    m2 <- if (M2 == 'x') {5} else {as.numeric(M2)}
    sum1 <- sum(u1 + m1, na.rm=T)
    sum2 <- sum(u2 + m2, na.rm=T)
    if (sum1 < sum2) {return(loc2)}
    if (sum2 < sum1) {return(loc1)}
    if (sum2 == sum1) {return(locr)}
  })
  slice(df, -rid)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate interp rmse !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
calculate_interp_rmse <- function(obs, interp, type='similarity') {
  if (!(type %in% c('krige', 'similarity'))) stop('\ninvalid type!')
  if (type == 'similarity') {
    nearest_est <- interp$est_sim[st_nearest_feature(obs, st_geometry(interp))]
    cat('\n', nearest_est)
    sqrt(mean((nearest_est - obs$obs)^2))
  } else if (type == 'krige') {
    interpolation <- interpolation
    nearest_est <- interp$est_krige[st_nearest_feature(obs, st_geometry(interp))]
    sqrt(mean((nearest_est - obs$obs)^2))
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# project obs to transect !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
project_obs_to_transect <- function(transect, shp_obs) {
  if (is.null(shp_obs)) {return(NULL)}
  if (st_crs(transect) != st_crs(shp_obs)) {
    stop('transect and shp_obs crs not the same!')
  }
  if (!any(class(transect) == 'sfc')) {
    stop('\ntransect needs to be an sf object!')
  }
  if (!any(class(shp_obs) == 'sf')) {
    stop('\nUnrecognized shp_obs passed to project_obs_to_transect() !')
  }
  if (!any(names(shp_obs) %in% c('obs', 'est_sim', 'est_krg'))) {
    stop('\nUnrecognized sf object passed to project_obs_to_transect() !')
  }
  if (any(names(shp_obs) == 'ghfdb')) {
    obs <- shp_obs$obs
    sigma <- NA
  } else if (any(names(shp_obs) == 'similarity')) {
    obs <- shp_obs$est_sim
    sigma <- shp_obs$sigma_sim
  } else if (any(names(shp_obs) == 'krige')) {
    obs <- shp_obs$est_krg
    sigma <- shp_obs$sigma_krg
  }
  projected_distances <- st_line_project(transect, st_geometry(shp_obs), normalized=T)
  st_as_sf(st_line_interpolate(transect, projected_distances)) %>% rename(geometry=x) %>%
    mutate(projected_distances=projected_distances, obs=obs, sigma=sigma, .before=geometry)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# fit loess to projected obs !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit_loess_to_projected_obs <- function(df, n=1e3, nmin=10, span=0.65) {
  df <- df %>% st_set_geometry(NULL)
  if (nrow(df) < nmin) {
    cat('\n   Cannot fit loess with less than', nmin, 'projected obs!')
    return(NULL)
  } else {
    mod <- NULL
    while (span <= 0.9 && is.null(mod)) {
      mod <- tryCatch({
        loess(obs~projected_distances, data=df, span=span)
      }, error=function(e) {
        span <- span + 0.05
        cat('\n   Loess failed with', nrow(df), ' obs! Increasing span to ', span)
        NULL
      })
    }
    if (is.null(mod)) {
      cat('\n   Loess failed!')
      return(NULL)
    } else {
      new_dist <- seq(0, 1, length.out = n)
      original_range <- range(df$projected_distances, na.rm = TRUE)
      loess_pred <- predict(mod, newdata=new_dist)
      loess_pred[new_dist < original_range[1] | new_dist > original_range[2]] <- NA
      smooth <- tibble(projected_distances=new_dist, obs=loess_pred)
      return(smooth)
    }
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate cross correlation !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
calculate_cross_correlation <- function(smooth1, smooth2) {
  len1 <- nrow(smooth1)
  len2 <- nrow(smooth2)
  if (len1 < len2) {
    smooth1 <- bind_rows(smooth1, tibble(projected_distances=rep(NA, len2 - len1),
                                         obs=rep(NA, len2 - len1)))
  } else if (len2 < len1) {
    smooth2 <- bind_rows(smooth2, tibble(projected_distances=rep(NA, len1 - len2),
                                         obs=rep(NA, len1 - len2)))
  }
  tryCatch({
    ccf_result <- ccf(smooth1$obs, smooth2$obs, plot=F)
    max_ccf <- max(ccf_result$acf)
    max_lag <- ccf_result$lag[which.max(ccf_result$acf)]
    return(max_ccf)
  }, error=function(e) {
    cat('\n!! ERROR occurred in calculate_cross_correlation:\n!!', conditionMessage(e))
    return(NULL)
  })
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load map data !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load_map_data <- function(fpath) {
  load(fpath, envir=parent.frame())
  if (!exists('shp_submap', envir=parent.frame())) {
    stop('\nMissing map data! Use "load(path/to/map-data.RData)"')
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load nlopt interpolation data !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load_nlopt_interpolation_data <- function(fpath) {
  load(fpath, envir=parent.frame())
  if (!exists('nlopt_summary', envir=parent.frame())) {
    stop('\nMissing nlopt data! Use "load(path/to/interpolation-summary.RData)"')
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# compile transect data !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
compile_transect_data <- function(trans_ids=NULL, lbuff=5e5, sbuff=c(5e4, 1e5, 1.5e5),
                                  fv=NULL, np=NULL) {
  load_map_data('assets/map_data/map-data.RData')
  if (is.null(trans_ids)) {stop('\nMissing submap transect ids!')}
  if (is.numeric(trans_ids) && all(trans_ids >= 1) && all(trans_ids <= nrow(shp_submap))) {
    x <- slice(shp_submap, trans_ids)
  } else if (is.character(trans_ids) && all(trans_ids %in% shp_submap$short_name)) {
    x <- filter(shp_submap, short_name %in% trans_ids)
  } else {
    stop('\nUnrecognized input for trans_ids! Use the transect short_name or row number ...')
  }
  if (nrow(x) == 0) {
    stop('\n', trans_ids, ' not found in shp_submap!')
  }
  suppressWarnings({suppressMessages({
    df <-
      x %>% rowwise() %>%
      mutate(large_buffer=st_buffer(transect, lbuff, endCapStyle='ROUND')) %>%
      mutate(small_buffer1=st_buffer(transect, sbuff[1], endCapStyle='FLAT')) %>%
      mutate(small_buffer2=st_buffer(transect, sbuff[2], endCapStyle='FLAT')) %>%
      mutate(small_buffer3=st_buffer(transect, sbuff[3], endCapStyle='FLAT')) %>%
      mutate(bbox=bbox_extend(st_bbox(large_buffer), square=T)) %>%
      mutate(ridge=crop_feature(shp_ridge, bbox)) %>%
      mutate(trench=crop_feature(shp_trench, bbox)) %>%
      mutate(transform=crop_feature(shp_transform, bbox)) %>%
      mutate(volcano=crop_feature(select(shp_volc, geometry), bbox)) %>%
      mutate(grid=list(crop_feature(shp_grid, large_buffer, T, T))) %>%
      mutate(ghfdb_large_buff=list(crop_feature(shp_ghfdb, large_buffer, T, T))) %>%
      mutate(ghfdb_small_buff1=list(crop_feature(shp_ghfdb, small_buffer1, T, T))) %>%
      mutate(ghfdb_small_buff2=list(crop_feature(shp_ghfdb, small_buffer2, T, T))) %>%
      mutate(ghfdb_small_buff3=list(crop_feature(shp_ghfdb, small_buffer3, T, T))) %>%
      mutate(ghfdb_projected1=list(project_obs_to_transect(transect, ghfdb_small_buff1))) %>%
      mutate(ghfdb_projected2=list(project_obs_to_transect(transect, ghfdb_small_buff2))) %>%
      mutate(ghfdb_projected3=list(project_obs_to_transect(transect, ghfdb_small_buff3))) %>%
      mutate(ghfdb_loess1=list(fit_loess_to_projected_obs(ghfdb_projected1))) %>%
      mutate(ghfdb_loess2=list(fit_loess_to_projected_obs(ghfdb_projected2))) %>%
      mutate(ghfdb_loess3=list(fit_loess_to_projected_obs(ghfdb_projected3))) %>%
      mutate(sim_large_buff=list(crop_feature(shp_sim, large_buffer, T, T))) %>%
      mutate(sim_small_buff1=list(crop_feature(shp_sim, small_buffer1, T, T))) %>%
      mutate(sim_small_buff2=list(crop_feature(shp_sim, small_buffer2, T, T))) %>%
      mutate(sim_small_buff3=list(crop_feature(shp_sim, small_buffer3, T, T))) %>%
      mutate(sim_projected1=list(project_obs_to_transect(transect, sim_small_buff1))) %>%
      mutate(sim_projected2=list(project_obs_to_transect(transect, sim_small_buff2))) %>%
      mutate(sim_projected3=list(project_obs_to_transect(transect, sim_small_buff3))) %>%
      mutate(sim_loess1=list(fit_loess_to_projected_obs(sim_projected1))) %>%
      mutate(sim_loess2=list(fit_loess_to_projected_obs(sim_projected2))) %>%
      mutate(sim_loess3=list(fit_loess_to_projected_obs(sim_projected3))) %>%
      mutate(bathy=list(get_seg_bathy(bbox)))
  })})
  if (!is.null(fv) && !is.null(np)) {
    shp_krg <- Krige(df$ghfdb_large_buff[[1]], fv, df$grid[[1]], np)
    df %>%
      mutate(krg_large_buff=list(crop_feature(shp_krg, large_buffer, T, T))) %>%
      mutate(krg_small_buff1=list(crop_feature(shp_krg, small_buffer1, T, T))) %>%
      mutate(krg_small_buff2=list(crop_feature(shp_krg, small_buffer2, T, T))) %>%
      mutate(krg_small_buff3=list(crop_feature(shp_krg, small_buffer3, T, T))) %>%
      mutate(krg_projected1=list(project_obs_to_transect(transect, krg_small_buff1))) %>%
      mutate(krg_projected2=list(project_obs_to_transect(transect, krg_small_buff2))) %>%
      mutate(krg_projected3=list(project_obs_to_transect(transect, krg_small_buff3))) %>%
      mutate(krg_loess1=list(fit_loess_to_projected_obs(krg_projected1))) %>%
      mutate(krg_loess2=list(fit_loess_to_projected_obs(krg_projected2))) %>%
      mutate(krg_loess3=list(fit_loess_to_projected_obs(krg_projected3))) %>%
      mutate(dff_large_buff=list(interp_diff(krg_large_buff, sim_large_buff))) %>%
      mutate(dff_small_buff1=list(interp_diff(krg_small_buff1, sim_small_buff1))) %>%
      mutate(dff_small_buff2=list(interp_diff(krg_small_buff2, sim_small_buff2))) %>%
      mutate(dff_small_buff3=list(interp_diff(krg_small_buff3, sim_small_buff3))) %>%
      mutate(dff_projected1=list(interp_diff(krg_projected1, sim_projected1))) %>%
      mutate(dff_projected2=list(interp_diff(krg_projected2, sim_projected2))) %>%
      mutate(dff_projected3=list(interp_diff(krg_projected3, sim_projected3))) %>%
      ungroup()
  } else {
    ungroup(df)
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get closest interp obs !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_closest_interp_obs <- function(trans_id=NULL, fv=NULL, np=NULL, thresh=1e4) {
  if (is.null(trans_id)) {stop('\nMissing submap transect ids!')}
  if (!is.null(fv) && !is.null(np)) {
    x <- compile_transect_data(trans_id, fv=fv, np=np)
    itp <- x$krg_large_buff[[1]]
  } else {
    x <- compile_transect_data(trans_id)
    itp <- x$sim_large_buff[[1]]
  }
  grid <- x$grid[[1]]
  obs <- x$ghfdb_large_buff[[1]]
  nearest_obs <- st_nearest_feature(grid, obs)
  nearest_itp <- st_nearest_feature(obs, grid)
  dt_obs <- st_distance(grid, obs[nearest_obs,], by_element=T) < set_units(thresh, 'm')
  dt_itp <- st_distance(obs, grid[nearest_itp,], by_element=T) < set_units(thresh, 'm')
  obs_n <- obs[nearest_obs,][dt_obs,]
  itp_n <- itp[nearest_itp,][dt_itp,]
  if (any(names(itp_n) %in% c('est_sim', 'similarity'))) {
    itp_n %>%
      mutate(obs_ghfdb=obs_n[st_nearest_feature(itp_n, obs_n),]$obs,
             ghfdb=obs_n[st_nearest_feature(itp_n, obs_n),]$ghfdb, .before=similarity) %>%
      select(-c(sigma_sim, obs_sim))
  } else if (any(names(itp_n) %in% c('est_krg', 'krige'))) {
    itp_n %>%
      mutate(obs_ghfdb=obs_n[st_nearest_feature(itp_n, obs_n),]$obs,
             ghfdb=obs_n[st_nearest_feature(itp_n, obs_n),]$ghfdb, .before=krige) %>%
      select(-c(sigma_krg, var_krg))
  } else {
    NULL
  }
}

#######################################################
## .2.             Kriging Functions             !!! ##
#######################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# experimental vgrm !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
experimental_vgrm <- function(shp_hf=NULL, cutoff=3, n_lags=30) {
  if (is.null(shp_hf)) {stop('\nMissing heat flow data!')}
  bbox <- st_bbox(shp_hf)
  bbox_diagonal_distance <- sqrt((bbox$xmax - bbox$xmin)^2 + (bbox$ymax - bbox$ymin)^2)
  lag_cutoff <- as.vector(bbox_diagonal_distance / cutoff)
  bin_width <- lag_cutoff / n_lags
  variogram(obs~1, locations=shp_hf, cutoff=lag_cutoff, width=bin_width)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sample and shuffle equal nfolds !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sample_and_shuffle_equal_nfolds <- function(nrow_data, nfold=10) {
  if (is.null(nrow_data)) {stop('\nMissing nrow data!')}
  if (is.null(nfold)) {nfold <- 10}
  fold_size <- nrow_data %/% nfold
  remainder <- nrow_data %% nfold
  fold <- rep(1:nfold, each=fold_size)
  if (remainder > 0) {fold <- c(fold, sample(1:nfold, remainder))}
  sample(fold)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# cost function !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cost_function <- function(shp_hf=NULL, cutoff=3, n_lags=50, n_max=10, max_dist=1e5,
                          v_mod='Sph', interp_weight=0.5, vgrm_weight=0.5, trans_id=NULL,
                          nfold=NULL) {
  if (is.null(shp_hf)) {stop('\nMissing heat flow data model!')}
  suppressWarnings({
    tryCatch({
      ev <- experimental_vgrm(shp_hf, cutoff, n_lags)
    }, error=function(e) {
      cat('\n!! ERROR occurred in experimental_vgrm:\n!!', conditionMessage(e))
      return(runif(1, 1, 1.5))
    })
    if (nrow(ev) < 2) {return(runif(1, 1, 1.5))}
    tryCatch({
      fv <- fit.variogram(ev, vgm(model=v_mod), fit.method=6)
    }, error=function(e) {
      cat('\n!! ERROR occurred in fit.variogram:\n!!', conditionMessage(e))
      return(runif(1, 1, 1.5))
    })
    if (fv$range < 0) {
      cat('\nVariogram range is negative:', fv$range)
      return(runif(1, 1, 1.5))
    }
    tryCatch({
      k_cv <- krige.cv(obs~1, shp_hf, model=fv, nmax=n_max, maxdist=max_dist, nfold=nfold)
    }, error=function(e) {
      cat('\n!! ERROR occurred in krige.cv:\n!!', conditionMessage(e))
      return(runif(1, 1, 1.5))
    })
    na_sum <- sum(is.na(k_cv$residual))
    na_thresh <- nrow(k_cv) / 2
    if (na_sum != 0) {
      if (na_sum >= na_thresh) {
        cat('\nToo many NAs in krige.cv:', na_sum, '/', na_thresh)
        return(runif(1, 1, 1.5))
      }
    }
    tryCatch({
      k_cv <- k_cv %>% filter(!is.na(residual))
      vgrm_rmse <- sqrt(sqrt(attr(fv, 'SSErr') / nrow(ev)))
      vgrm_sd <- sqrt(sd(ev$gamma, na.rm=T))
      vgrm_cost <- vgrm_weight * vgrm_rmse / vgrm_sd
      interp_rmse <- sqrt(sum(k_cv$residual^2, na.rm=T) / nrow(k_cv))
      interp_sd <- sd(k_cv$var1.pred, na.rm=T)
      interp_cost <- interp_weight * interp_rmse / interp_sd
    }, error=function(e) {
      cat('\n!! ERROR occurred in cost_function:\n!!', conditionMessage(e))
      return(runif(1, 1, 1.5))
    })
  })
  log_dir <- 'assets/nlopt_data/nlopt_itr'
  log_file <- paste0(log_dir, '/nlopt-out-', trans_id, '-', v_mod)
  if (!dir.exists(log_dir)) {dir.create(log_dir, recursive=T, showWarnings=F)}
  if (!file.exists(log_file)) {file.create(log_file, showWarnings=F)}
  sink(log_file, append=T)
  cat('\nTransect:', trans_id)
  cat('\nN obs:', nrow(shp_hf))
  cat('\nCutoff:', cutoff)
  cat('\nN lags:', n_lags)
  cat('\nN pairs:', n_max)
  cat('\nMax dist:', max_dist)
  cat('\nVgrm weight:', vgrm_weight)
  cat('\nVgrm rmse:', vgrm_rmse)
  cat('\nVgrm sd:', vgrm_sd)
  cat('\nVgrm cost:', vgrm_cost)
  cat('\nInt weight:', interp_weight)
  cat('\nInt rmse:', interp_rmse)
  cat('\nInt sd:', interp_sd)
  cat('\nInt cost:', interp_cost)
  cat('\nCost:', vgrm_cost + interp_cost)
  cat('\nVgrm model:', v_mod)
  cat('\n', rep('+', 30), sep='')
  cat('\n', rep('-', 40), sep='')
  sink()
  vgrm_cost + interp_cost
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# decode opt !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
decode_opt <- function(shp_hf=NULL, v_mod=NULL, opt=NULL) {
  suppressWarnings({
    if (is.null(v_mod)) {stop('\nMissing variogram model!')}
    if (is.null(shp_hf)) {stop('\nMissing heat flow data!')}
    if (is.null(opt)) {stop('\nMissing nloptr object!')}
    tryCatch({
      ev <- experimental_vgrm(shp_hf, opt$solution[1], opt$solution[2])
    }, error=function(e) {
      stop('\n!! ERROR occurred in experimental_vgrm:\n!!', conditionMessage(e))
    })
    if (nrow(ev) < 2) {stop('\nExperimental variogram has less than two lags!')}
    if (any(class(ev) == 'try-error')) {stop('\nExperimental variogram error!')}
    tryCatch({
      fv <- fit.variogram(ev, vgm(model=v_mod), fit.method=6)
    }, error=function(e) {
      stop('\n!! ERROR occurred in fit.variogram:\n!!', conditionMessage(e))
    })
  })
  return(list('experimental_vgrm'=as_tibble(ev), 'fitted_vgrm'=fv))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# nlopt krige !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nlopt_krige <- function(trans_id=NULL, v_mod='Sph', alg='NLOPT_LN_COBYLA', max_eval=500,
                        iwt=0.5, vwt=0.5) {
  if (is.null(trans_id)) {stop('\nMissing submap transect id!')}
  nlopt_dir <- 'assets/nlopt_data/nlopt_transects'
  nlopt_id <- paste0('opt-', trans_id, '-', v_mod) 
  nlopt_path <- paste0(nlopt_dir, '/', nlopt_id, '.RData')
  nlopt_itr_path <- paste0('assets/nlopt_data/nlopt_itr/nlopt-out-', trans_id, '-', v_mod)
  if (file.exists(nlopt_path) & file.exists(nlopt_itr_path)) {
    cat('\nOptimized', v_mod, 'kriging model found for submap transect:', trans_id)
    return(invisible())
  }
  if ((!file.exists(nlopt_path) & file.exists(nlopt_itr_path)) |
      (file.exists(nlopt_path) & !file.exists(nlopt_itr_path))) {
    cat('\n!!Optimized', v_mod, 'kriging model failed for submap transect:', trans_id, '!!')
    cat('\n     Remove', nlopt_path, 'and ...')
    cat('\n     Remove', nlopt_itr_path, 'and ...')
    cat('\n     Run "make nlopt" again ...')
    return(invisible())
  }
  if (!dir.exists(nlopt_dir)) {dir.create(nlopt_dir, recursive=T, showWarnings=F)}
  x0 <- c(3, 50, 5, 1e5) # Initial values (cutoff, n_lags, n_max, max_dist)
  lb <- c(1, 30, 2, 5e4) # Lower bound (cutoff, n_lags, n_max, max_dist)
  ub <- c(12, 100, 30, 5e5) # Upper bound (cutoff, n_lags, n_max, max_dist)
  opts <- list(print_level=0, maxeval=max_eval, algorithm=alg, xtol_rel=1e-8, ftol_rel=1e-8)
  x <- compile_transect_data(trans_id)
  obs <- x$ghfdb_large_buff[[1]]
  folds <- sample_and_shuffle_equal_nfolds(nrow(obs))
  nlopt_fun <- function(x) {
    cost_function(obs, x[1], x[2], x[3], x[4], v_mod, iwt, vwt, trans_id, folds)
  }
  tryCatch({
    cat('\n', rep('-', 60), sep='')
    cat('\nOptimizing', v_mod, 'kriging model for submap transect:', trans_id)
    cat('\n', rep('+', 60), sep='')
    cat('\nNLopt algorithm:    ', alg)
    cat('\nMaximum evaluations:', max_eval)
    cat('\nKrige weight:       ', iwt)
    cat('\nVariogram weight:   ', vwt)
    cat('\n                    (cutoff, n_lags, n_max, max_dist)')
    cat('\nInitial parameters: ', x0)
    cat('\nLower bounds:       ', lb)
    cat('\nUpper bounds:       ', ub)
    opt <- nloptr(x0, nlopt_fun, lb=lb, ub=ub, opts=opts)
  }, error=function(e) {
    cat('\n!! ERROR occurred in nlopt_krige:\n!!', conditionMessage(e))
  })
  if (opt$status < 0 | opt$status == 5 | opt$iterations < 10) {
    cat('\n', rep('+', 60), sep='')
    cat('\nNLopt failed to converge:')
    cat('\n   NLopt status:', opt$status)
    cat('\n   NLopt iterations:', opt$iterations)
    cat('\n   ', opt$message)
    cat('\n', rep('-', 60), sep='')
    return(invisible())
  } else {
    cat('\n', rep('+', 60), sep='')
    cat('\nNLopt converged:')
    cat('\n   NLopt status:', opt$status)
    cat('\n   NLopt iterations:', opt$iterations)
    cat('\n   ', opt$message)
    cat('\n', rep('-', 60), sep='')
    opt_decoded <- decode_opt(obs, v_mod, opt)
    assign(str_replace_all(nlopt_id, '-', '_'), opt_decoded)
    save(list=str_replace_all(nlopt_id, '-', '_'), file=nlopt_path)
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# nlopt transects !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nlopt_transects <- function(trans_ids=NULL, v_mods=c('Sph', 'Exp', 'Lin'),
                            alg='NLOPT_LN_COBYLA', max_eval=500, iwt=0.5, vwt=0.5,
                            parallel=T) {
  if (is.null(trans_ids)) {stop('\nMissing submap transect ids!')}
  if (length(list.files('assets/map_data/relief')) < length(trans_ids)) {parallel <- F}
  x <- expand.grid(id=trans_ids, vm=v_mods, stringsAsFactors=F) %>% arrange(id, vm)
  cat('Optimizing krige models after Li et al. (2018) ...\n')
  if (parallel) {
    plan(multisession, workers=availableCores() - 2)
    future_walk2(x$id, x$vm, ~nlopt_krige(..., alg, max_eval, iwt, vwt),
                 .options=furrr_options(seed=seed), .progress=T)
  } else {
    walk2(x$id, x$vm, ~nlopt_krige(..., alg, max_eval, iwt, vwt))
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read nloptr itr !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
read_nloptr_itr <- function(fpath) {
  read_file <- function(path, ...) {
    con <- file(path)
    on.exit(close(con))
    suppressWarnings(readLines(con, ...))
  }
  t <- read_file(fpath)
  trans_id_t <- t[grepl('^Transect: ', t)] %>% map_chr(~gsub('Transect: ', '', .x))
  n_obs_t <- t[grepl('^N obs: ', t)] %>% map_dbl(~as.numeric(gsub('N obs: ', '', .x)))
  cutoff_t <- t[grepl('^Cutoff: ', t)] %>% map_dbl(~as.numeric(gsub('Cutoff: ', '', .x)))
  n_lags_t <- t[grepl('^N lags: ', t)] %>% map_dbl(~as.numeric(gsub('N lags: ', '', .x)))
  n_pairs_t <- t[grepl('^N pairs: ', t)] %>% map_dbl(~as.numeric(gsub('N pairs: ', '', .x)))
  max_d_t <- t[grepl('^Max dist: ', t)] %>% map_dbl(~as.numeric(gsub('Max dist: ', '', .x)))
  vmod_t <- t[grepl('^Vgrm model:', t)] %>% map_chr(~gsub('Vgrm model: ', '', .x))
  vwt_t <- t[grepl('^Vgrm weight', t)] %>% map_dbl(~as.numeric(gsub('Vgrm weight: ', '', .x)))
  vrmse_t <- t[grepl('^Vgrm rmse', t)] %>% map_dbl(~as.numeric(gsub('Vgrm rmse: ', '', .x)))
  vcost_t <- t[grepl('^Vgrm cost', t)] %>% map_dbl(~as.numeric(gsub('Vgrm cost: ', '', .x)))
  iwt_t <- t[grepl('^Int weight', t)] %>% map_dbl(~as.numeric(gsub('Int weight: ', '', .x)))
  irmse_t <- t[grepl('^Int rmse', t)] %>% map_dbl(~as.numeric(gsub('Int rmse: ', '', .x)))
  icost_t <- t[grepl('^Int cost', t)] %>% map_dbl(~as.numeric(gsub('Int cost: ', '', .x)))
  cost_t <- t[grepl('^Cost', t)] %>% map_dbl(~as.numeric(gsub('Cost: ', '', .x)))
  tibble(short_name=trans_id_t, n_obs=n_obs_t, cutoff=cutoff_t, n_lags=n_lags_t,
         n_pairs=n_pairs_t, max_dist=max_d_t, v_mod=vmod_t, vgrm_wt=vwt_t, vgrm_rmse=vrmse_t,
         vgrm_cost=vcost_t, cv_wt=iwt_t, cv_rmse=irmse_t, cv_cost=icost_t, cost=cost_t) %>%
  group_by(short_name, v_mod) %>% mutate(itr=row_number(), .after=short_name) %>% ungroup()
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get optimal krige model !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_optimal_krige_model <- function(trans_id=NULL) {
  if (is.null(trans_id)) {stop('\nMissing submap transect id!')}
  nlopt_dir <- 'assets/nlopt_data/'
  nlopt_itr_dir <- paste0(nlopt_dir, 'nlopt_itr')
  nlopt_trans_dir <- paste0(nlopt_dir, 'nlopt_transects')
  nlopt_itr_paths <- list.files(nlopt_itr_dir, pattern=trans_id, full.names=T)
  nlopt_trans_paths <- list.files(nlopt_trans_dir, pattern=trans_id, full.names=T)
  if (length(nlopt_itr_paths) < 1) {
    cat('\n   No nlopt itr files found for:', trans_id)
    return(NULL)
  }
  if (length(nlopt_trans_paths) < 1) {
    cat('\n   No nlopt itr files found for:', trans_id)
    return(NULL)
  }
  if (length(nlopt_itr_paths) != length(nlopt_trans_paths)) {
    itr_mods <- str_sub(nlopt_itr_paths, start=-3)
    trans_mods <- str_extract(nlopt_trans_paths, '.{3}(?=\\.RData)')
    if (length(trans_mods) < 1) {
      cat('\n   No nlopt transect RData found for:', trans_id)
      return(NULL)
    }
    nlopt_itr_paths <-
      nlopt_itr_paths[str_detect(nlopt_itr_paths, paste(trans_mods, collapse='|'))]
  }
  nlopt_itr <- map_df(nlopt_itr_paths, read_nloptr_itr)
  k_mod <- slice_min(nlopt_itr, cost)
  if (nrow(k_mod) > 1) {k_mod <- slice(k_mod, nrow(k_mod))}
  nlopt_path <- paste0(nlopt_trans_dir, '/opt-', trans_id, '-', k_mod$v_mod, '.RData')
  if (file.exists(nlopt_path)) {
    load(nlopt_path)
    opt_decoded <- get(paste0('opt_', trans_id, '_', k_mod$v_mod))
    list('nlopt_itr'=nlopt_itr, 'opt_krige_mod_summary'=k_mod,
         'experimental_vgrm'=opt_decoded$experimental_vgrm,
         'fitted_vgrm'=opt_decoded$fitted_vgrm)
  } else {
    cat('\n   No nlopt RData found for:', trans_id)
    return(NULL)
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Krige !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Krige <- function(shp_hf=NULL, fv=NULL, shp_grid=NULL, n_max=10, max_dist=1e5) {
  if (is.null(shp_hf)) {stop('\nMissing heat flow data!')}
  if (is.null(fv)) {stop('\nMissing variogram model!')}
  if (is.null(shp_grid)) {stop('\nMissing kriging locations (grid)!')}
  if (is.null(n_max)) {stop('\nMissing max pairs!')}
  if (is.null(max_dist)) {stop('\nMissing max distance!')}
  suppressWarnings({
    tryCatch({
      krige(obs~1, shp_hf, newdata=shp_grid, model=fv, nmax=n_max, maxdist=max_dist,
            debug.level=0) %>% as_tibble() %>% st_as_sf() %>%
        rename(est_krg=var1.pred, var_krg=var1.var, krige=geometry) %>%
        mutate(sigma_krg=sqrt(var_krg), .before=krige) %>%
        mutate(est_krg=ifelse(est_krg > 0 & est_krg <= 250, est_krg, NA),
               var_krg=ifelse(est_krg > 0 & est_krg <= 250, est_krg, NA))
    }, error=function(e) {
      cat('\n!! ERROR occurred in Krige:\n!!', conditionMessage(e))
    })
  })
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# interp diff !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
interp_diff <- function(shp_krg=NULL, shp_sim=NULL) {
  if (is.null(shp_krg)) {stop('\nMissing krige interpolation!')}
  if (is.null(shp_sim)) {stop('\nMissing similarity interpolation!')}
  if (any(names(shp_sim) == 'projected_distances')) {
    shp_sim <-
      shp_sim %>% rename(est_sim=obs, sigma_sim=sigma) %>%
      select(est_sim, sigma_sim, geometry)
    shp_krg <-
      shp_krg %>% rename(est_krg=obs, sigma_krg=sigma) %>%
      select(est_krg, sigma_krg, geometry)
    shp_sim %>% 
      mutate(est_krg=shp_krg$est_krg, sigma_krg=shp_krg$sigma_krg) %>%
      mutate(est_dff=est_sim-est_krg, sigma_dff=sigma_sim-sigma_krg, .before=geometry) %>%
      rename(projected_diff=geometry) %>% select(est_dff, sigma_dff, projected_diff)
  } else {
    shp_sim %>% 
      mutate(est_krg=shp_krg$est_krg, sigma_krg=shp_krg$sigma_krg) %>%
      mutate(est_dff=est_sim-est_krg, sigma_dff=sigma_sim-sigma_krg, .before=similarity) %>%
      rename(diff=similarity) %>% select(est_dff, sigma_dff, diff)
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# summarize optimal krige models !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summarize_optimal_krige_models <- function(trans_ids, parallel=T) {
  f <- function(x) {get_optimal_krige_model(x)$opt_krige_mod_summary}
  cat('\nSummarizing optimal krige models ...\n')
  if (parallel) {
    plan(multisession, workers=availableCores() - 2)
    as_tibble(future_map_dfr(trans_ids, f, .options=furrr_options(seed=seed), .progress=T))
  } else {
    map_df(trans_ids, f)
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# summarize interpolation differences !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summarize_interp_differences <- function(trans_ids, parallel=T) {
  f <- function(x) {
    tryCatch({
      km <- get_optimal_krige_model(x)
      fv <- km$fitted_vgrm
      np <- km$opt_krige_mod_summary$n_pairs
      compile_transect_data(x, fv=fv, np=np)$dff_large_buff[[1]] %>%
        st_set_geometry(NULL) %>%
        summarise(short_name=x, n_grid=n(), n=sum(!is.na(est_dff)),
                  rmse=sqrt(mean(est_dff^2, na.rm=T)), min=min(est_dff, na.rm=T),
                  max=max(est_dff, na.rm=T), med=median(est_dff, na.rm=T),
                  iqr=IQR(est_dff, na.rm=T), mean=mean(est_dff, na.rm=T),
                  sigma=sd(est_dff, na.rm=T))
    }, error=function(e) {
      cat('\n!! ERROR occurred in summarize_interp_differences:\n!!', conditionMessage(e))
      return(NULL)
    })
  }
  cat('\nSummarizing interpolation differences ...\n')
  if (parallel) {
    as_tibble(future_map_dfr(ids, f, .options=furrr_options(seed=seed), .progress=T))
  } else {
    map_df(ids, f)
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# summarize interpolation accuracy !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summarize_interp_accuracy <- function(trans_ids, parallel=T) {
  f <- function(x) {
    suppressWarnings({suppressMessages({
      tryCatch({
        km <- get_optimal_krige_model(x)
        fv <- km$fitted_vgrm
        np <- km$opt_krige_mod_summary$n_pairs
        comps_sim <- get_closest_interp_obs(x)
        comps_krg <- get_closest_interp_obs(x, fv=fv, np=np)
        tibble(short_name=x,
               n_sim=nrow(comps_sim[!is.na(comps_sim$est_sim),]),
               min_obs_sim=min(comps_sim$est_sim - comps_sim$obs_ghfdb, na.rm=T),
               max_obs_sim=max(comps_sim$est_sim - comps_sim$obs_ghfdb, na.rm=T),
               mean_obs_sim=mean(comps_sim$est_sim - comps_sim$obs_ghfdb, na.rm=T),
               sd_obs_sim=sd(comps_sim$est_sim - comps_sim$obs_ghfdb, na.rm=T),
               med_obs_sim=median(comps_sim$est_sim - comps_sim$obs_ghfdb, na.rm=T),
               iqr_obs_sim=IQR(comps_sim$est_sim - comps_sim$obs_ghfdb, na.rm=T),
               rmse_obs_sim=sqrt(mean((comps_sim$est_sim - comps_sim$obs_ghfdb)^2, na.rm=T)),
               n_krg=nrow(comps_krg[!is.na(comps_krg$est_krg),]),
               min_obs_krg=min(comps_krg$est_krg - comps_krg$obs_ghfdb, na.rm=T),
               max_obs_krg=max(comps_krg$est_krg - comps_krg$obs_ghfdb, na.rm=T),
               mean_obs_krg=mean(comps_krg$est_krg - comps_krg$obs_ghfdb, na.rm=T),
               sd_obs_krg=sd(comps_krg$est_krg - comps_krg$obs_ghfdb, na.rm=T),
               med_obs_krg=median(comps_krg$est_krg - comps_krg$obs_ghfdb, na.rm=T),
               iqr_obs_krg=IQR(comps_krg$est_krg - comps_krg$obs_ghfdb, na.rm=T),
               rmse_obs_krg=sqrt(mean((comps_krg$est_krg - comps_krg$obs_ghfdb)^2, na.rm=T)))
      }, error=function(e) {
        cat('\n!! ERROR occurred in summarize_interp_accuracy:\n!!', conditionMessage(e))
        return(NULL)
      })
    })})
  }
  cat('\nSummarizing interpolation accuracies ...\n')
  if (parallel) {
    as_tibble(future_map_dfr(ids, f, .options=furrr_options(seed=seed), .progress=T))
  } else {
    map_df(ids, f)
  }
}

#######################################################
## .3.            Plotting Functions             !!! ##
#######################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot ghfdb !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_ghfdb <- function() {
  load_map_data('assets/map_data/map-data.RData')
  fig_path <- 'figs/ghfdb.png'
  if (!dir.exists('figs')) {dir.create('figs', recursive=T, showWarnings=F)}
  if (file.exists(fig_path)) {
    cat('\n   ', fig_path, ' already exists ...', sep='')
    return(invisible())
  }
  cat('\n   Plotting:', fig_path)
  suppressWarnings({suppressMessages({
    p1 <-
      ggplot() +
      geom_sf(data=shp_relief_world, aes(color=elev), shape=15, size=0.01) +
      scale_color_etopo(guide='none') +
      geom_sf(data=shp_ridge, linewidth=0.3, color='white') +
      geom_sf(data=shp_transform, linewidth=0.3, color='white') +
      geom_sf(data=shp_trench, linewidth=0.3, color='white') +
      geom_sf(data=shp_submap, color='black') +
      ggtitle('Submap Transects') +
      coord_sf(expand=F, lims_method='geometry_bbox') +
      theme_map(font_size=14)
    p2 <-
      ggplot() +
      geom_sf(data=shp_relief_world, aes(color=elev), shape=15, size=0.01) +
      scale_color_etopo(guide='none') +
      new_scale_color() +
      geom_sf(data=shp_ridge, linewidth=0.3, color='white') +
      geom_sf(data=shp_transform, linewidth=0.3, color='white') +
      geom_sf(data=shp_trench, linewidth=0.3, color='white') +
      geom_sf(data=shp_ghfdb, aes(color=obs), size=0.1, shape=20) +
      scale_color_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)),
                            limits=c(0, 250), breaks=c(0, 125, 250), na.value='transparent',
                            guide=guide_colorbar(title.vjust=1, show.limits=T,
                                                 frame.colour='black',
                                                 ticks.colour='black')) +
      ggtitle('Global HF Observations') +
      coord_sf(expand=F, lims_method='geometry_bbox') +
      theme_map(font_size=14)
    p3 <- p1 / p2 &
      theme(plot.margin=margin(), legend.position='bottom',
            legend.justification='center', legend.direction='horizontal',
            axis.text=element_blank(), legend.margin=margin(),
            legend.box.margin=margin(5, 5, 5, 5), legend.key.height=unit(0.5, 'cm'),
            legend.key.width=unit(1, 'cm'),
            legend.title=element_text(vjust=0, color='black', size=14),
            panel.grid=element_line(linewidth=0.05, color='grey20'),
            plot.title=element_text(vjust=0, hjust=0.5, margin=margin(10, 10, 10, 10)))
    ggsave(file=fig_path, plot=p3, width=6.5, height=6.5, bg='white')
  })})
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot transect buff comp !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_transect_buff_comp <- function(trans_id, base_size=20) {
  fig_dir <- 'figs/transect_buff/'
  fig_path1 <- paste0(fig_dir, trans_id, '-transect-buff-ghfdb.png')
  fig_path2 <- paste0(fig_dir, trans_id, '-transect-buff-sim.png')
  fig_path3 <- paste0(fig_dir, trans_id, '-transect-buff-krg.png')
  fig_path_comp <- paste0(fig_dir, trans_id, '-transect-buff-comp.png')
  if (!dir.exists(fig_dir)) {dir.create(fig_dir, recursive=T, showWarnings=F)}
  if (file.exists(fig_path_comp) & file.exists(fig_path1) & file.exists(fig_path2) &
      file.exists(fig_path3)) {
    cat('\n   ', fig_path1, ' already exists ...', sep='')
    cat('\n   ', fig_path2, ' already exists ...', sep='')
    cat('\n   ', fig_path3, ' already exists ...', sep='')
    cat('\n   ', fig_path_comp, ' already exists ...', sep='')
    return(invisible())
  }
  tryCatch({
    cat('\n   Plotting: ', fig_path1, ' ...', sep='')
    cat('\n   Plotting: ', fig_path2, ' ...', sep='')
    cat('\n   Plotting: ', fig_path3, ' ...', sep='')
    cat('\n   Plotting: ', fig_path_comp, ' ...', sep='')
    km <- get_optimal_krige_model(trans_id)
    fv <- km$fitted_vgrm
    np <- km$opt_krige_mod_summary$n_pairs
    x <- compile_transect_data(trans_id, fv=fv, np=np)
    comps_sim <- get_closest_interp_obs(trans_id)
    comps_krg <- get_closest_interp_obs(trans_id, fv=fv, np=np)
    hf_pt_size <- 3
    pt_stroke <- 0.8
    p_title <- paste0('Submap Transect: ', trans_id, ' ', x$trench_name)
    map_theme <- list(
      theme_bw(base_size=base_size),
      theme(plot.margin=margin(5, 5, 5, 5),
            plot.title=element_text(vjust=0, hjust=0.5, margin=margin(0, 5, 10, 5))),
      annotation_scale(location='bl', width_hint=0.33, text_cex=1.6, style='ticks',
                       line_width=4, text_face='bold'),
      annotation_north_arrow(location='bl', which_north='true', pad_x=unit(0.0, 'cm'),
                             height=unit(2, 'cm'), width=unit(2, 'cm'),
                             pad_y=unit(0.5, 'cm'), style=north_arrow_fancy_orienteering)
    )
    profile_theme <-
      list(
        theme_bw(base_size=base_size),
        theme(panel.grid=element_blank(), panel.background=element_rect(fill='grey90'),
              plot.margin=margin(5, 5, 5, 5),
              legend.justification='right', legend.position='inside',
              legend.position.inside=c(0.92, 0.85), legend.direction='horizontal',
              legend.key.height=unit(0.5, 'cm'), legend.key.width=unit(1, 'cm'),
              legend.box.margin=margin(2, 2, 2, 2), legend.margin=margin(),
              legend.title=element_text(vjust=0, size=base_size),
              legend.background=element_blank())
      )
    suppressWarnings({suppressMessages({
      p1 <-
        ggplot(x) +
        geom_sf(data=x$bathy[[1]], aes(color=elev), size=0.5, shape=15) +
        scale_color_etopo(guide='none') +
        new_scale_color() +
        geom_sf(aes(geometry=large_buffer), fill=NA, linewidth=0.5) +
        geom_sf(aes(geometry=ridge), color='white') +
        geom_sf(aes(geometry=transform), color='white') +
        geom_sf(aes(geometry=trench), color='white', linewidth=1.5) +
        geom_sf(aes(geometry=transect), color='black', linewidth=1.5) +
        geom_sf(aes(geometry=volcano), color='black', fill='white', shape=24) +
        geom_sf(data=x$ghfdb_large_buff[[1]], aes(color=obs), shape=20, size=hf_pt_size) +
        geom_sf(aes(geometry=small_buffer1), fill=NA, linewidth=0.5, color='black') +
        geom_sf(aes(geometry=small_buffer2), fill=NA, linewidth=0.5, color='black') +
        geom_sf(aes(geometry=small_buffer3), fill=NA, linewidth=0.5, color='black') +
        ggtitle('Global HF Observations') +
        scale_color_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                              breaks=c(0, 125, 250), na.value='transparent', guide='none') +
        coord_sf(expand=F, lims_method='geometry_bbox') + map_theme
      if (is.null(x$ghfdb_projected1[[1]])) {
        p2 <- ggplot(data=data.frame(), aes())
      } else {
        p2 <-
          ggplot() +
          geom_point(data=x$ghfdb_projected1[[1]], aes(projected_distances, obs), shape=20,
                     color='grey20', size=2) +
          geom_point(data=x$ghfdb_projected2[[1]], aes(projected_distances, obs), shape=20,
                     color='grey20', size=2) +
          geom_point(data=x$ghfdb_projected3[[1]], aes(projected_distances, obs), shape=20,
                     color='grey20', size=2) +
          geom_path(data=x$ghfdb_loess1[[1]], aes(projected_distances, obs), color='black') +
          geom_path(data=x$ghfdb_loess2[[1]], aes(projected_distances, obs), color='black') +
          geom_path(data=x$ghfdb_loess3[[1]], aes(projected_distances, obs), color='black') +
          geom_rug(data=x$ghfdb_loess3[[1]], aes(projected_distances, obs, color=obs),
                   sides='b', length=unit(0.06, 'npc'))
      }
      p2 <- p2 +
        labs(x='Normalized Distance', y=NULL) +
        scale_color_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                              breaks=c(0, 125, 250), na.value='transparent',
                              guide=guide_colorbar(title.vjust=1, show.limits=T,
                                                   frame.colour='black',
                                                   ticks.colour='black')) +
        scale_y_continuous(limits=c(0, 250), breaks=seq(0, 250, 50)) +
        scale_x_continuous(limits=c(0, 1), breaks=c(0, 0.5, 1)) + profile_theme
      p <- (p1 / p2) + plot_layout(widths=1, heights=c(1.5, 1)) &
        theme(plot.title=element_text(size=base_size * 1.5, margin=margin()))
      ggsave(file=fig_path1, plot=p, width=6.5, height=10, dpi=300, bg='white')
      p3 <-
        ggplot(x) +
        geom_sf(data=x$bathy[[1]], aes(color=elev), size=0.5, shape=15) +
        scale_color_etopo(guide='none') +
        new_scale_color() +
        geom_sf(aes(geometry=large_buffer), fill=NA, linewidth=0.5) +
        geom_sf(aes(geometry=ridge), color='white') +
        geom_sf(aes(geometry=transform), color='white') +
        geom_sf(aes(geometry=trench), color='white', linewidth=1.5) +
        geom_sf(aes(geometry=transect), color='black', linewidth=1.5) +
        geom_sf(aes(geometry=volcano), color='black', fill='white', shape=24) +
        geom_sf(data=x$sim_large_buff[[1]], aes(color=est_sim), shape=20, size=hf_pt_size) +
        geom_sf(data=comps_sim[!is.na(comps_sim$est_sim),], aes(fill=est_sim), shape=21,
                size=hf_pt_size, stroke=pt_stroke, color='white') +
        geom_sf(aes(geometry=small_buffer1), fill=NA, linewidth=0.5, color='black') +
        geom_sf(aes(geometry=small_buffer2), fill=NA, linewidth=0.5, color='black') +
        geom_sf(aes(geometry=small_buffer3), fill=NA, linewidth=0.5, color='black') +
        ggtitle('Similarity Interpolation') +
        scale_color_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                              breaks=c(0, 125, 250), na.value='transparent', guide='none') +
        scale_fill_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                              breaks=c(0, 125, 250), na.value='transparent', guide='none') +
        coord_sf(expand=F, lims_method='geometry_bbox') + map_theme
      if (is.null(x$sim_projected1[[1]])) {
        p4 <- ggplot(data=data.frame(), aes())
      } else {
        p4 <-
          ggplot() +
          geom_point(data=x$sim_projected1[[1]], aes(projected_distances, obs), shape=20,
                     color='grey20', size=2) +
          geom_point(data=x$sim_projected2[[1]], aes(projected_distances, obs), shape=20,
                     color='grey20', size=2) +
          geom_point(data=x$sim_projected3[[1]], aes(projected_distances, obs), shape=20,
                     color='grey20', size=2) +
          geom_path(data=x$sim_loess1[[1]], aes(projected_distances, obs), color='black') +
          geom_path(data=x$sim_loess2[[1]], aes(projected_distances, obs), color='black') +
          geom_path(data=x$sim_loess3[[1]], aes(projected_distances, obs), color='black') +
          geom_rug(data=x$sim_loess3[[1]], aes(projected_distances, obs, color=obs),
                   sides='b', length=unit(0.06, 'npc'))
      }
      p4 <- p4 +
        labs(x='Normalized Distance', y=NULL) +
        scale_color_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                              breaks=c(0, 125, 250), na.value='transparent',
                              guide=guide_colorbar(title.vjust=1, show.limits=T,
                                                   frame.colour='black',
                                                   ticks.colour='black')) +
        scale_y_continuous(limits=c(0, 250), breaks=seq(0, 250, 50)) +
        scale_x_continuous(limits=c(0, 1), breaks=c(0, 0.5, 1)) + profile_theme
      p <- (p3 / p4) + plot_layout(widths=1, heights=c(1.5, 1)) &
        theme(plot.title=element_text(size=base_size * 1.5, margin=margin()))
      ggsave(file=fig_path2, plot=p, width=6.5, height=10, dpi=300, bg='white')
      p5 <-
        ggplot(x) +
        geom_sf(data=x$bathy[[1]], aes(color=elev), size=0.5, shape=15) +
        scale_color_etopo(guide='none') +
        new_scale_color() +
        geom_sf(aes(geometry=large_buffer), fill=NA, linewidth=0.5) +
        geom_sf(aes(geometry=ridge), color='white') +
        geom_sf(aes(geometry=transform), color='white') +
        geom_sf(aes(geometry=trench), color='white', linewidth=1.5) +
        geom_sf(aes(geometry=transect), color='black', linewidth=1.5) +
        geom_sf(aes(geometry=volcano), color='black', fill='white', shape=24) +
        geom_sf(data=x$krg_large_buff[[1]], aes(color=est_krg), shape=20, size=hf_pt_size) +
        geom_sf(data=comps_krg[!is.na(comps_krg$est_krg),], aes(fill=est_krg), shape=21,
                size=hf_pt_size, stroke=pt_stroke, color='white') +
        geom_sf(aes(geometry=small_buffer1), fill=NA, linewidth=0.5, color='black') +
        geom_sf(aes(geometry=small_buffer2), fill=NA, linewidth=0.5, color='black') +
        geom_sf(aes(geometry=small_buffer3), fill=NA, linewidth=0.5, color='black') +
        ggtitle('Krige Interpolation') +
        scale_color_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                              breaks=c(0, 125, 250), na.value='transparent', guide='none') +
        scale_fill_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                              breaks=c(0, 125, 250), na.value='transparent', guide='none') +
        coord_sf(expand=F, lims_method='geometry_bbox') + map_theme
      if (is.null(x$krg_projected1[[1]])) {
        p6 <- ggplot(data=data.frame(), aes())
      } else {
        p6 <-
          ggplot() +
          geom_point(data=x$krg_projected1[[1]], aes(projected_distances, obs), shape=20,
                     color='grey20', size=2) +
          geom_point(data=x$krg_projected2[[1]], aes(projected_distances, obs), shape=20,
                     color='grey20', size=2) +
          geom_point(data=x$krg_projected3[[1]], aes(projected_distances, obs), shape=20,
                     color='grey20', size=2) +
          geom_path(data=x$krg_loess1[[1]], aes(projected_distances, obs), color='black') +
          geom_path(data=x$krg_loess2[[1]], aes(projected_distances, obs), color='black') +
          geom_path(data=x$krg_loess3[[1]], aes(projected_distances, obs), color='black') +
          geom_rug(data=x$krg_loess3[[1]], aes(projected_distances, obs, color=obs),
                   sides='b', length=unit(0.06, 'npc'))
      }
      p6 <- p6 +
        labs(x='Normalized Distance', y=NULL) +
        scale_color_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                              breaks=c(0, 125, 250), na.value='transparent',
                              guide=guide_colorbar(title.vjust=1, show.limits=T,
                                                   frame.colour='black',
                                                   ticks.colour='black')) +
        scale_y_continuous(limits=c(0, 250), breaks=seq(0, 250, 50)) +
        scale_x_continuous(limits=c(0, 1), breaks=c(0, 0.5, 1)) + profile_theme
      p <- (p5 / p6) + plot_layout(widths=1, heights=c(1.5, 1)) &
        theme(plot.title=element_text(size=base_size * 1.5, margin=margin()))
      ggsave(file=fig_path3, plot=p, width=6.5, height=10, dpi=300, bg='white')
      p7 <-
        (p1 + p5 + p3) / (p2 + p6 + p4) +
        plot_layout(widths=1, heights=c(1.5, 1)) +
        plot_annotation(title=p_title,
                        tag_levels=list(c('a)', 'b)', 'c)', '  ', '  ', '  ')),
                        theme=theme(plot.title=element_text(size=base_size * 1.5,
                                                            margin=margin()))) &
        theme(plot.tag=element_text(size=base_size * 1.5, margin=margin(0, 0, -10, 0)))
      ggsave(file=fig_path_comp, plot=p7, width=19.5, height=11, dpi=300, bg='white')
    })})
  }, error=function(e) {
    cat('\n!! ERROR occurred in plot_transect:\n!!', conditionMessage(e))
  })
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot transect neighbors comp !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_transect_neighbors_comp <- function(trans_ids, base_size=20) {
  if (length(trans_ids) < 3) {stop('\nNeed at least 3 transect ids!')}
  if (length(trans_ids) > 3) {trans_ids <- trans_ids[1:3]}
  fig_dir <- 'figs/transect_neighbors/'
  trans_id_lab <- paste0(trans_ids[1], '-', trans_ids[2], '-', trans_ids[3])
  fig_paths <- list(paste0(fig_dir, trans_id_lab, '-transect-neighbors-ghfdb.png'),
                    paste0(fig_dir, trans_id_lab, '-transect-neighbors-sim.png'),
                    paste0(fig_dir, trans_id_lab, '-transect-neighbors-krg.png'))
  if (!dir.exists(fig_dir)) {dir.create(fig_dir, recursive=T, showWarnings=F)}
  if (file.exists(fig_paths[[1]]) & file.exists(fig_paths[[2]]) &
      file.exists(fig_paths[[3]])) {
    cat('\n   ', fig_paths[[1]], ' already exists ...', sep='')
    cat('\n   ', fig_paths[[2]], ' already exists ...', sep='')
    cat('\n   ', fig_paths[[3]], ' already exists ...', sep='')
    return(invisible())
  }
  f <- function(trans_id) {
    tryCatch({
      km <- get_optimal_krige_model(trans_id)
      fv <- km$fitted_vgrm
      np <- km$opt_krige_mod_summary$n_pairs
      x <- compile_transect_data(trans_id, fv=fv, np=np)
      comps_sim <- get_closest_interp_obs(trans_id)
      comps_krg <- get_closest_interp_obs(trans_id, fv=fv, np=np)
      hf_pt_size <- 3
      pt_stroke <- 0.8
      map_theme <- list(
        theme_bw(base_size=base_size),
        theme(plot.margin=margin(5, 5, 5, 5),
              plot.title=element_text(vjust=0, hjust=0.5, margin=margin(0, 5, 10, 5))),
        annotation_scale(location='bl', width_hint=0.33, text_cex=1.6, style='ticks',
                         line_width=4, text_face='bold'),
        annotation_north_arrow(location='bl', which_north='true', pad_x=unit(0.0, 'cm'),
                               height=unit(2, 'cm'), width=unit(2, 'cm'),
                               pad_y=unit(0.5, 'cm'), style=north_arrow_fancy_orienteering)
      )
      profile_theme <-
        list(
          theme_bw(base_size=base_size),
          theme(panel.grid=element_blank(), panel.background=element_rect(fill='grey90'),
                plot.margin=margin(5, 5, 5, 5),
                legend.justification='right', legend.position='inside',
                legend.position.inside=c(0.92, 0.85), legend.direction='horizontal',
                legend.key.height=unit(0.5, 'cm'), legend.key.width=unit(1, 'cm'),
                legend.box.margin=margin(2, 2, 2, 2), legend.margin=margin(),
                legend.title=element_text(vjust=0, size=base_size),
                legend.background=element_blank())
        )
      suppressWarnings({suppressMessages({
        p1 <-
          ggplot(x) +
          geom_sf(data=x$bathy[[1]], aes(color=elev), size=0.5, shape=15) +
          scale_color_etopo(guide='none') +
          new_scale_color() +
          geom_sf(aes(geometry=large_buffer), fill=NA, linewidth=0.5) +
          geom_sf(aes(geometry=ridge), color='white') +
          geom_sf(aes(geometry=transform), color='white') +
          geom_sf(aes(geometry=trench), color='white', linewidth=1.5) +
          geom_sf(aes(geometry=transect), color='black', linewidth=1.5) +
          geom_sf(aes(geometry=volcano), color='black', fill='white', shape=24) +
          geom_sf(data=x$ghfdb_large_buff[[1]], aes(color=obs), shape=20, size=hf_pt_size) +
          geom_sf(aes(geometry=small_buffer1), fill=NA, linewidth=0.5, color='black') +
          geom_sf(aes(geometry=small_buffer2), fill=NA, linewidth=0.5, color='black') +
          geom_sf(aes(geometry=small_buffer3), fill=NA, linewidth=0.5, color='black') +
          ggtitle('Global HF Observations') +
          scale_color_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                                breaks=c(0, 125, 250), na.value='transparent',
                                guide='none') +
          coord_sf(expand=F, lims_method='geometry_bbox') + map_theme
        if (is.null(x$ghfdb_projected1[[1]])) {
          p2 <- ggplot(data=data.frame(), aes())
        } else {
          p2 <-
            ggplot() +
            geom_point(data=x$ghfdb_projected1[[1]], aes(projected_distances, obs), shape=20,
                       color='grey20', size=2) +
            geom_point(data=x$ghfdb_projected2[[1]], aes(projected_distances, obs), shape=20,
                       color='grey20', size=2) +
            geom_point(data=x$ghfdb_projected3[[1]], aes(projected_distances, obs), shape=20,
                       color='grey20', size=2) +
            geom_path(data=x$ghfdb_loess1[[1]], aes(projected_distances, obs),
                      color='black') +
            geom_path(data=x$ghfdb_loess2[[1]], aes(projected_distances, obs),
                      color='black') +
            geom_path(data=x$ghfdb_loess3[[1]], aes(projected_distances, obs),
                      color='black') +
            geom_rug(data=x$ghfdb_loess3[[1]], aes(projected_distances, obs, color=obs),
                     sides='b', length=unit(0.06, 'npc'))
        }
        p2 <- p2 +
          labs(x='Normalized Distance', y=NULL) +
          scale_color_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                                breaks=c(0, 125, 250), na.value='transparent',
                                guide=guide_colorbar(title.vjust=1, show.limits=T,
                                                     frame.colour='black',
                                                     ticks.colour='black')) +
          scale_y_continuous(limits=c(0, 250), breaks=seq(0, 250, 50)) +
          scale_x_continuous(limits=c(0, 1), breaks=c(0, 0.5, 1)) + profile_theme
        p3 <-
          ggplot(x) +
          geom_sf(data=x$bathy[[1]], aes(color=elev), size=0.5, shape=15) +
          scale_color_etopo(guide='none') +
          new_scale_color() +
          geom_sf(aes(geometry=large_buffer), fill=NA, linewidth=0.5) +
          geom_sf(aes(geometry=ridge), color='white') +
          geom_sf(aes(geometry=transform), color='white') +
          geom_sf(aes(geometry=trench), color='white', linewidth=1.5) +
          geom_sf(aes(geometry=transect), color='black', linewidth=1.5) +
          geom_sf(aes(geometry=volcano), color='black', fill='white', shape=24) +
          geom_sf(data=x$sim_large_buff[[1]], aes(color=est_sim), shape=20,
                  size=hf_pt_size) +
          geom_sf(data=comps_sim[!is.na(comps_sim$est_sim),], aes(fill=est_sim),
                  shape=21, size=hf_pt_size, stroke=pt_stroke, color='white') +
          geom_sf(aes(geometry=small_buffer1), fill=NA, linewidth=0.5, color='black') +
          geom_sf(aes(geometry=small_buffer2), fill=NA, linewidth=0.5, color='black') +
          geom_sf(aes(geometry=small_buffer3), fill=NA, linewidth=0.5, color='black') +
          ggtitle('Similarity Interpolation') +
          scale_color_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                                breaks=c(0, 125, 250), na.value='transparent',
                                guide='none') +
          scale_fill_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                                breaks=c(0, 125, 250), na.value='transparent',
                                guide='none') +
          coord_sf(expand=F, lims_method='geometry_bbox') + map_theme
        if (is.null(x$sim_projected1[[1]])) {
          p4 <- ggplot(data=data.frame(), aes())
        } else {
          p4 <-
            ggplot() +
            geom_point(data=x$sim_projected1[[1]], aes(projected_distances, obs), shape=20,
                       color='grey20', size=2) +
            geom_point(data=x$sim_projected2[[1]], aes(projected_distances, obs), shape=20,
                       color='grey20', size=2) +
            geom_point(data=x$sim_projected3[[1]], aes(projected_distances, obs), shape=20,
                       color='grey20', size=2) +
            geom_path(data=x$sim_loess1[[1]], aes(projected_distances, obs), color='black') +
            geom_path(data=x$sim_loess2[[1]], aes(projected_distances, obs), color='black') +
            geom_path(data=x$sim_loess3[[1]], aes(projected_distances, obs), color='black') +
            geom_rug(data=x$sim_loess3[[1]], aes(projected_distances, obs, color=obs),
                     sides='b', length=unit(0.06, 'npc'))
        }
        p4 <- p4 +
          labs(x='Normalized Distance', y=NULL) +
          scale_color_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                                breaks=c(0, 125, 250), na.value='transparent',
                                guide=guide_colorbar(title.vjust=1, show.limits=T,
                                                     frame.colour='black',
                                                     ticks.colour='black')) +
          scale_y_continuous(limits=c(0, 250), breaks=seq(0, 250, 50)) +
          scale_x_continuous(limits=c(0, 1), breaks=c(0, 0.5, 1)) + profile_theme
        p5 <-
          ggplot(x) +
          geom_sf(data=x$bathy[[1]], aes(color=elev), size=0.5, shape=15) +
          scale_color_etopo(guide='none') +
          new_scale_color() +
          geom_sf(aes(geometry=large_buffer), fill=NA, linewidth=0.5) +
          geom_sf(aes(geometry=ridge), color='white') +
          geom_sf(aes(geometry=transform), color='white') +
          geom_sf(aes(geometry=trench), color='white', linewidth=1.5) +
          geom_sf(aes(geometry=transect), color='black', linewidth=1.5) +
          geom_sf(aes(geometry=volcano), color='black', fill='white', shape=24) +
          geom_sf(data=x$krg_large_buff[[1]], aes(color=est_krg), shape=20,
                  size=hf_pt_size) +
          geom_sf(data=comps_krg[!is.na(comps_krg$est_krg),], aes(fill=est_krg),
                  shape=21, size=hf_pt_size, stroke=pt_stroke, color='white') +
          geom_sf(aes(geometry=small_buffer1), fill=NA, linewidth=0.5, color='black') +
          geom_sf(aes(geometry=small_buffer2), fill=NA, linewidth=0.5, color='black') +
          geom_sf(aes(geometry=small_buffer3), fill=NA, linewidth=0.5, color='black') +
          ggtitle('Krige Interpolation') +
          scale_color_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                                breaks=c(0, 125, 250), na.value='transparent',
                                guide='none') +
          scale_fill_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                                breaks=c(0, 125, 250), na.value='transparent',
                                guide='none') +
          coord_sf(expand=F, lims_method='geometry_bbox') + map_theme
        if (is.null(x$krg_projected1[[1]])) {
          p6 <- ggplot(data=data.frame(), aes())
        } else {
          p6 <-
            ggplot() +
            geom_point(data=x$krg_projected1[[1]], aes(projected_distances, obs), shape=20,
                       color='grey20', size=2) +
            geom_point(data=x$krg_projected2[[1]], aes(projected_distances, obs), shape=20,
                       color='grey20', size=2) +
            geom_point(data=x$krg_projected3[[1]], aes(projected_distances, obs), shape=20,
                       color='grey20', size=2) +
            geom_path(data=x$krg_loess1[[1]], aes(projected_distances, obs), color='black') +
            geom_path(data=x$krg_loess2[[1]], aes(projected_distances, obs), color='black') +
            geom_path(data=x$krg_loess3[[1]], aes(projected_distances, obs), color='black') +
            geom_rug(data=x$krg_loess3[[1]], aes(projected_distances, obs, color=obs),
                     sides='b', length=unit(0.06, 'npc'))
        }
        p6 <- p6 +
          labs(x='Normalized Distance', y=NULL) +
          scale_color_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                                breaks=c(0, 125, 250), na.value='transparent',
                                guide=guide_colorbar(title.vjust=1, show.limits=T,
                                                     frame.colour='black',
                                                     ticks.colour='black')) +
          scale_y_continuous(limits=c(0, 250), breaks=seq(0, 250, 50)) +
          scale_x_continuous(limits=c(0, 1), breaks=c(0, 0.5, 1)) + profile_theme
        return(list(p1, p2, p3, p4, p5, p6))
      })})
    }, error=function(e) {
      cat('\n!! ERROR occurred in plot_transect:\n!!', conditionMessage(e))
    })
  }
  plan(multisession, workers=availableCores() - 2)
  plots <- future_map(trans_ids, f, .options=furrr_options(seed=seed)) %>% reduce(c)
  suppressWarnings({suppressMessages({
    p_title <- paste0('Submap Transect: ', trans_id_lab, ' ')
    p1 <-
      (plots[[1]] + plots[[7]] + plots[[13]]) / (plots[[2]] + plots[[8]] + plots[[14]]) +
      plot_layout(widths=1, heights=c(1.5, 1)) +
      plot_annotation(title=p_title, tag_levels=list(c('a)', 'b)', 'c)', '  ', '  ', '  ')),
                      theme=theme(plot.title=element_text(size=base_size * 1.5,
                                                          margin=margin()))) &
      theme(plot.tag=element_text(size=base_size * 1.5, margin=margin(0, 0, -10, 0)))
    cat('\n   Plotting: ', fig_paths[[1]], ' ...', sep='')
    ggsave(file=fig_paths[[1]], plot=p1, width=19.5, height=11, dpi=300, bg='white')
    p2 <-
      (plots[[3]] + plots[[9]] + plots[[15]]) / (plots[[4]] + plots[[10]] + plots[[16]]) +
      plot_layout(widths=1, heights=c(1.5, 1)) +
      plot_annotation(title=p_title, tag_levels=list(c('a)', 'b)', 'c)', '  ', '  ', '  ')),
                      theme=theme(plot.title=element_text(size=base_size * 1.5,
                                                          margin=margin()))) &
      theme(plot.tag=element_text(size=base_size * 1.5, margin=margin(0, 0, -10, 0)))
    cat('\n   Plotting: ', fig_paths[[2]], ' ...', sep='')
    ggsave(file=fig_paths[[2]], plot=p2, width=19.5, height=11, dpi=300, bg='white')
    p3 <-
      (plots[[5]] + plots[[11]] + plots[[17]]) / (plots[[6]] + plots[[12]] + plots[[18]]) +
      plot_layout(widths=1, heights=c(1.5, 1)) +
      plot_annotation(title=p_title, tag_levels=list(c('a)', 'b)', 'c)', '  ', '  ', '  ')),
                      theme=theme(plot.title=element_text(size=base_size * 1.5,
                                                          margin=margin()))) &
      theme(plot.tag=element_text(size=base_size * 1.5, margin=margin(0, 0, -10, 0)))
    cat('\n   Plotting: ', fig_paths[[3]], ' ...', sep='')
    ggsave(file=fig_paths[[3]], plot=p3, width=19.5, height=11, dpi=300, bg='white')
  })})
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot transect strip comp !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_transect_strip_comp <- function(trans_ids=NULL, fname=NULL, buff=1, base_size=14) {
  if (is.null(trans_ids)) {stop('\nMissing submap transect ids!')}
  if (is.null(fname)) {fname <- 'test'}
  if (buff == 1) {
    select_cols <- c('short_name', 'ghfdb_loess1', 'krg_loess1', 'sim_loess1')
  } else if (buff == 2) {
    select_cols <- c('short_name', 'ghfdb_loess2', 'krg_loess2', 'sim_loess2')
  } else if (buff == 3) {
    select_cols <- c('short_name', 'ghfdb_loess3', 'krg_loess3', 'sim_loess3')
  } else {
    select_cols <- c('short_name', 'ghfdb_loess1', 'krg_loess1', 'sim_loess1')
  }
  fig_dir <- 'figs/summary/'
  if (!dir.exists(fig_dir)) {dir.create(fig_dir, recursive=T, showWarnings=F)}
  fig_path <- paste0(fig_dir, 'strip-', fname, '.png')
  if (file.exists(fig_path)) {
    cat('\n   ', fig_path, ' already exists ...', sep='')
    return(invisible())
  }
  tryCatch({
    cat('\n   Plotting: ', fig_path, ' ...', sep='')
    f <- function(x) {
      km <- get_optimal_krige_model(x)
      fv <- km$fitted_vgrm
      np <- km$opt_krige_mod_summary$n_pairs
      compile_transect_data(x, fv=fv, np=np)
    }
    plan(multisession, workers=availableCores() - 2)
    df <-
      future_map(trans_ids, f, .options=furrr_options(seed=seed)) %>%
      reduce(rbind) %>%
      st_set_geometry(NULL) %>%
      select(all_of(select_cols)) %>%
      pivot_longer(-c(short_name)) %>%
      unnest(value)
    p <-
      ggplot(df) +
      geom_vline(aes(xintercept=projected_distances, color=obs)) +
      scale_x_continuous(breaks=c(0, 1)) +
      labs(x='Normalized Distance', y=NULL) +
      facet_grid(vars(short_name), vars(name)) +
      scale_color_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                            breaks=c(0, 125, 250), na.value='transparent', guide='none') +
      theme_bw(base_size=base_size) +
      theme(panel.grid=element_blank(), panel.background=element_rect(fill='grey90'),
            plot.margin=margin(5, 5, 5, 5),
            legend.justification='right', legend.position='inside',
            legend.position.inside=c(0.92, 0.85), legend.direction='horizontal',
            legend.key.height=unit(0.5, 'cm'), legend.key.width=unit(1, 'cm'),
            legend.box.margin=margin(2, 2, 2, 2), legend.margin=margin(),
            legend.title=element_text(vjust=0, size=base_size),
            legend.background=element_blank())
    ggsave(file=fig_path, plot=p, width=6.5, height=length(trans_ids) * 0.75, dpi=300,
           bg='white')
  }, error=function(e) {
    cat('\n!! ERROR occurred in plot_transect_strip_comp:\n!!', conditionMessage(e))
  })
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot cross correlation !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_cross_correlation <- function(trans_id1, trans_id2, min_ccf=0.8, base_size=14) {
  fig_dir <- 'figs/transect_xcorr/'
  fig_path <- paste0(fig_dir, trans_id1, '-', trans_id2, '-transect-xcorr.png')
  if (!dir.exists(fig_dir)) {dir.create(fig_dir, recursive=T, showWarnings=F)}
  if (file.exists(fig_path)) {
    cat('\n   ', fig_path, ' already exists ...', sep='')
    return(invisible())
  }
  tryCatch({
    cat('\n   Plotting: ', fig_path, ' ...', sep='')
    kmx <- get_optimal_krige_model(trans_id1)
    kmy <- get_optimal_krige_model(trans_id2)
    fvx <- kmx$fitted_vgrm
    fvy <- kmy$fitted_vgrm
    npx <- kmx$opt_krige_mod_summary$n_pairs
    npy <- kmy$opt_krige_mod_summary$n_pairs
    x <- compile_transect_data(trans_id1, fv=fvx, np=npx)
    y <- compile_transect_data(trans_id2, fv=fvy, np=npy)
    color_mapping <- setNames(c('darkorange', 'navy'), c(trans_id1, trans_id2))
    p_theme <- list(
      theme_bw(base_size=base_size),
      theme(panel.grid=element_blank(), panel.background=element_rect(fill='grey90'),
            plot.margin=margin(5, 5, 5, 5), plot.title=element_text(hjust=0.5), 
            legend.justification='right', legend.position='inside',
            legend.position.inside=c(0.92, 0.85), legend.direction='horizontal',
            legend.key.height=unit(0.5, 'cm'), legend.key.width=unit(1, 'cm'),
            legend.box.margin=margin(2, 2, 2, 2), legend.margin=margin(),
            legend.title=element_text(vjust=0, size=base_size),
            legend.background=element_blank())
    )
    smooth1_ghfdb <- x$ghfdb_loess3[[1]]
    points1_ghfdb <- x$ghfdb_projected3[[1]]
    smooth2_ghfdb <- y$ghfdb_loess3[[1]]
    points2_ghfdb <- y$ghfdb_projected3[[1]]
    ccf_result_ghfdb <- ccf(smooth1_ghfdb$obs, smooth2_ghfdb$obs, na.action=na.pass)
    max_ccf_ghfdb <- max(ccf_result_ghfdb$acf)
    max_lag_ghfdb <- ccf_result_ghfdb$lag[which(ccf_result_ghfdb$acf == max_ccf_ghfdb)]
    smooth1_sim <- x$sim_loess3[[1]]
    points1_sim <- x$sim_projected3[[1]]
    smooth2_sim <- y$sim_loess3[[1]]
    points2_sim <- y$sim_projected3[[1]]
    ccf_result_sim <- ccf(smooth1_sim$obs, smooth2_sim$obs, na.action=na.pass)
    max_ccf_sim <- max(ccf_result_sim$acf)
    max_lag_sim <- ccf_result_sim$lag[which(ccf_result_sim$acf == max_ccf_sim)]
    smooth1_krg <- x$krg_loess3[[1]]
    points1_krg <- x$krg_projected3[[1]]
    smooth2_krg <- y$krg_loess3[[1]]
    points2_krg <- y$krg_projected3[[1]]
    ccf_result_krg <- ccf(smooth1_krg$obs, smooth2_krg$obs, na.action=na.pass)
    max_ccf_krg <- max(ccf_result_krg$acf)
    max_lag_krg <- ccf_result_krg$lag[which(ccf_result_krg$acf == max_ccf_krg)]
    if (max_ccf_ghfdb > min_ccf | max_ccf_sim > min_ccf | max_ccf_krg > min_ccf) {
      p1 <-
        ggplot() +
        geom_point(data=points1_ghfdb, aes(projected_distances, obs, color=trans_id1),
                   na.rm=T, shape=20, size=0.7) +
        geom_point(data=points2_ghfdb, aes(projected_distances, obs, color=trans_id2),
                   na.rm=T, shape=20, size=0.7) +
        geom_path(data=smooth1_ghfdb, aes(projected_distances, obs, color=trans_id1),
                  na.rm=T) +
        geom_path(data=smooth2_ghfdb, aes(projected_distances, obs, color=trans_id2),
                  na.rm=T) +
        labs(x='Normalized Distance', y=bquote('Q'~(mWm^-2))) +
        ggtitle('Global Heat Flow') +
        scale_color_manual(name=NULL, values=color_mapping) + p_theme
      p2 <-
        ggplot(tibble(acf=ccf_result_ghfdb$acf, lag=ccf_result_ghfdb$lag)) +
        geom_hline(yintercept=min_ccf, linetype=2) +
        geom_hline(yintercept=-min_ccf, linetype=2) +
        geom_path(aes(lag, acf)) +
        geom_point(data=tibble(acf=max_ccf_ghfdb, lag=max_lag_ghfdb), aes(lag, acf)) +
        labs(x='Lag', y='Correlation') +
        scale_y_continuous(limits=c(-1, 1), breaks=seq(-1, 1, 0.5)) + p_theme
      p3 <-
        ggplot() +
        geom_point(data=points1_sim, aes(projected_distances, obs, color=trans_id1),
                   na.rm=T, shape=20, size=0.7) +
        geom_point(data=points2_sim, aes(projected_distances, obs, color=trans_id2),
                   na.rm=T, shape=20, size=0.7) +
        geom_path(data=smooth1_sim, aes(projected_distances, obs, color=trans_id1),
                  na.rm=T) +
        geom_path(data=smooth2_sim, aes(projected_distances, obs, color=trans_id2),
                  na.rm=T) +
        labs(x='Normalized Distance', y=bquote('Q'~(mWm^-2))) +
        ggtitle('Similarity') +
        scale_color_manual(name=NULL, values=color_mapping) + p_theme
      p4 <-
        ggplot(tibble(acf=ccf_result_sim$acf, lag=ccf_result_sim$lag)) +
        geom_hline(yintercept=min_ccf, linetype=2) +
        geom_hline(yintercept=-min_ccf, linetype=2) +
        geom_path(aes(lag, acf)) +
        geom_point(data=tibble(acf=max_ccf_sim, lag=max_lag_sim), aes(lag, acf)) +
        labs(x='Lag', y='Correlation') +
        scale_y_continuous(limits=c(-1, 1), breaks=seq(-1, 1, 0.5)) + p_theme
      p5 <-
        ggplot() +
        geom_point(data=points1_krg, aes(projected_distances, obs, color=trans_id1),
                   na.rm=T, shape=20, size=0.7) +
        geom_point(data=points2_krg, aes(projected_distances, obs, color=trans_id2),
                   na.rm=T, shape=20, size=0.7) +
        geom_path(data=smooth1_krg, aes(projected_distances, obs, color=trans_id1),
                  na.rm=T) +
        geom_path(data=smooth2_krg, aes(projected_distances, obs, color=trans_id2),
                  na.rm=T) +
        labs(x='Normalized Distance', y=bquote('Q'~(mWm^-2))) +
        ggtitle('Krige') +
        scale_color_manual(name=NULL, values=color_mapping) + p_theme
      p6 <-
        ggplot(tibble(acf=ccf_result_krg$acf, lag=ccf_result_krg$lag)) +
        geom_hline(yintercept=min_ccf, linetype=2) +
        geom_hline(yintercept=-min_ccf, linetype=2) +
        geom_path(aes(lag, acf)) +
        geom_point(data=tibble(acf=max_ccf_krg, lag=max_lag_krg), aes(lag, acf)) +
        labs(x='Lag', y='Correlation') +
        scale_y_continuous(limits=c(-1, 1), breaks=seq(-1, 1, 0.5)) + p_theme
      p_title <- paste0('Cross-correlation: ', trans_id1, ' ', trans_id2)
      p_rax <- theme(axis.title.y=element_blank(), axis.text.y=element_blank())
      p <-
        (p1 + (p3 + p_rax) + (p5 + p_rax)) / (p2 + (p4 + p_rax) + (p6 + p_rax)) +
        plot_annotation(title=p_title,
                        tag_levels=list(c('a)', 'b)', 'c)', '  ', '  ', '  ')),
                        theme=theme(plot.title=element_text(size=base_size * 1.5,
                                                            margin=margin()))) &
        theme(plot.tag=element_text(size=base_size * 1.5, margin=margin(0, 0, -10, 0)))
      ggsave(file=fig_path, plot=p, width=13, height=5, dpi=300, bg='white')
    } else {
      cat('\n   No significant cross-correlations at 95% CI ...')
    }
  }, error=function(e) {
    cat('\n!! ERROR occurred in plot_cross_correlation:\n!!', conditionMessage(e))
  })
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot optimal krige model !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_optimal_krige_model <- function(trans_id=NULL, base_size=22) {
  if (is.null(trans_id)) {stop('\nMissing submap transect ids!')}
  fig_dir <- 'figs/nlopt/'
  fig_path <- paste0(fig_dir, trans_id, '-opt-krige.png')
  if (!dir.exists(fig_dir)) {dir.create(fig_dir, recursive=T, showWarnings=F)}
  if (file.exists(fig_path)) {
    cat('\n   ', fig_path, ' already exists ...', sep='')
    return(invisible())
  }
  tryCatch({
    cat('\n   Plotting:', fig_path)
    opt_krige_mod <- get_optimal_krige_model(trans_id)
    opt_kmod <- opt_krige_mod$opt_krige_mod_summary
    nlopt_itr <- opt_krige_mod$nlopt_itr
    ev <- opt_krige_mod$experimental_vgrm
    fv <- opt_krige_mod$fitted_vgrm
    fv_line <- variogramLine(fv, maxdist=max(ev$dist))
    p_theme <-
      list(theme_bw(base_size=base_size),
           theme(panel.grid=element_blank(), panel.background=element_rect(fill='grey90'),
                 plot.margin=margin(5, 5, 5, 5),
                 plot.title=element_text(vjust=0, hjust=0.5, margin=margin(10, 10, 10, 10)),
                 legend.position='none'))
    suppressWarnings({suppressMessages({
      p0 <-
        ggplot(filter(nlopt_itr, v_mod == opt_kmod$v_mod)) +
        geom_path(aes(itr, max_dist / 1e3), color='forestgreen', linewidth=1) +
        geom_point(data=opt_kmod, aes(itr, max_dist / 1e3), color='forestgreen', size=3) +
        geom_label_repel(data=opt_kmod, size=base_size * 0.28,
                         aes(itr, max_dist / 1e3, label=round(max_dist / 1e3))) +
        annotate('text', x=Inf, y=Inf, label='Variogram Max Distance', hjust=1.05, vjust=1.5,
                 size=base_size * 0.35) +
        labs(x='NLopt Iteration', y='Distance (km)', color=NULL) + p_theme
      p1 <-
        ggplot(filter(nlopt_itr, v_mod == opt_kmod$v_mod)) +
        geom_path(aes(itr, n_pairs), color='firebrick', linewidth=1) +
        geom_point(data=opt_kmod, aes(itr, n_pairs), color='firebrick', size=3) +
        geom_label_repel(data=opt_kmod, aes(itr, n_pairs, label=round(n_pairs)),
                         size=base_size * 0.28) +
        annotate('text', x=Inf, y=Inf, label='Variogram Pairs', hjust=1.05, vjust=1.5,
                 size=base_size * 0.35) +
        labs(x='NLopt Iteration', y='Pairs', color=NULL) + p_theme
      p2 <-
        ggplot(filter(nlopt_itr, v_mod == opt_kmod$v_mod)) +
        geom_path(aes(itr, cutoff), color='orchid4', linewidth=1) +
        geom_point(data=opt_kmod, aes(itr, cutoff), color='orchid4', size=3) +
        geom_label_repel(data=opt_kmod, aes(itr, cutoff, label=round(cutoff)),
                         size=base_size * 0.28) +
        annotate('text', x=Inf, y=Inf, label='Variogram Cutoff', hjust=1.05, vjust=1.5,
                 size=base_size * 0.35) +
        labs(x='NLopt Iteration', y='Cutoff', color=NULL) + p_theme
      p3 <-
        ggplot(filter(nlopt_itr, v_mod == opt_kmod$v_mod)) +
        geom_path(aes(itr, n_lags), color='saddlebrown', linewidth=1) +
        geom_point(data=opt_kmod, aes(itr, n_lags), color='saddlebrown', size=3) +
        geom_label_repel(data=opt_kmod, aes(itr, n_lags, label=round(n_lags)),
                         size=base_size * 0.28) +
        annotate('text', x=Inf, y=Inf, label='Variogram Lags', hjust=1.05, vjust=1.5,
                 size=base_size * 0.35) +
        labs(x='NLopt Iteration', y='Lags', color=NULL) + p_theme
      p4 <-
        ggplot(filter(nlopt_itr, v_mod == opt_kmod$v_mod)) +
        geom_path(aes(itr, cost), color='black', linewidth=1) +
        geom_point(data=opt_kmod, aes(itr, cost), color='black', shape=20, size=5) +
        geom_label_repel(data=opt_kmod, aes(itr, cost, label=paste0(round(cost, 3))),
                         size=base_size * 0.28) +
        annotate('text', x=Inf, y=Inf, label='Cost Function', hjust=1.05, vjust=1.5,
                 size=base_size * 0.35) +
        labs(x='NLopt Iteration', y='Cost', color=NULL) + p_theme
      p5 <- 
        ggplot(ev) +
        geom_point(aes(x=dist / 1e3, y=sqrt(gamma)), shape=20) +
        geom_line(data=fv_line, aes(x=dist / 1e3, y=sqrt(gamma)), linewidth=1) +
        annotate('text', x=Inf, y=-Inf, label=paste0('Optimal model (', opt_kmod$v_mod, ')'),
                 hjust=1.05, vjust=-0.5, size=base_size * 0.35) +
        labs(x='Lag Distance (km)', y=bquote('Variance'~(mWm^-2))) + p_theme
      p6 <- (p0 + p1) / (p2 + p3) / (p4 + p5) +
        plot_annotation(
          title=paste0('Kriging Optimization: ', unique(nlopt_itr$short_name)),
          theme=theme(plot.title=element_text(size=base_size * 1.2, hjust=0.5))) +
        plot_layout(widths=1, heights=1) &
        theme(plot.tag = element_text(size=base_size * 1.5))
      ggsave(file=fig_path, plot=p6, width=13, height=10, dpi=300, bg='white')
    })})
  }, error=function(e) {
    cat('\n!! ERROR occurred in plot_optimal_variogram:\n!!', conditionMessage(e))
  })
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot interp accuracy summary !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_interp_accuracy_summary <- function(submap_zone=NULL, base_size=22) {
  load_nlopt_interpolation_data('assets/nlopt_data/interpolation-summary.RData')
  if (is.null(submap_zone)) {stop('\nMissing submap zone!')}
  fig_dir <- 'figs/summary/'
  if (!dir.exists(fig_dir)) {dir.create(fig_dir, recursive=T, showWarnings=F)}
  f <- function(x, y) {
    if (file.exists(y)) {
      cat('\n   ', y, ' already exists ...', sep='')
      return(invisible())
    }
    tryCatch({
      cat('\n   Plotting:', y)
      p_theme <- list(
        theme_bw(base_size=14),
        theme(panel.grid=element_blank(), panel.background=element_rect(fill='grey90'),
              axis.text.x=element_text(angle=40, hjust=1), legend.position='bottom',
              plot.title=element_text(hjust=0.5)))
      outliers_diff <-
        interp_diff_summary$short_name[which(interp_diff_summary$mean + 2 *
                                             interp_diff_summary$sigma > 300)]
      outliers_acc <-
        interp_accuracy_summary$short_name[which(interp_accuracy_summary$rmse_obs_sim > 300 |
                                                 interp_accuracy_summary$rmse_obs_krg > 300)]
      outliers <- unique(c(outliers_diff, outliers_acc))
      df_acc_rmse <-
        interp_accuracy_summary %>%
        filter(!is.na(rmse_obs_krg)) %>%
        pivot_longer(-c(short_name)) %>%
        mutate(method=ifelse(str_detect(name, 'sim'), 'sim', 'krg'),
               name=str_split(name, '_', simplify=T)[,1],
               zone=str_sub(short_name, 1, 3)) %>%
        rename(metric=name) %>% filter(zone == submap_zone) %>%
        filter(short_name %in% x & !(short_name %in% outliers)) %>%
        mutate(method=ifelse(method == 'sim', 'Similarity', 'Krige')) %>%
        filter(metric == 'rmse') %>% rename(rmse=value) %>% select(-metric) %>%
        left_join(select(nlopt_summary, short_name, n_obs), by='short_name') %>%
        left_join(select(interp_accuracy_summary, short_name, n_krg), by='short_name') %>%
        rename(n_control=n_krg) %>%
        left_join(select(interp_diff_summary, short_name, n_grid), by='short_name') %>%
        left_join(select(interp_diff_summary, short_name, n), by='short_name') %>%
        mutate(n=ifelse(method == 'Similarity', n_grid, n)) %>%
        rename(n_itp=n) %>% mutate(coverage=n_itp / n_grid * 100)
      suppressWarnings({suppressMessages({
        p0 <-
          nlopt_summary %>%
          filter(short_name %in% x & !(short_name %in% outliers)) %>%
          ggplot() +
          geom_col(aes(short_name, n_obs, fill='Observation'), color='black') +
          geom_col(data=filter(df_acc_rmse, method == 'Krige'),
                   aes(short_name, n_control, fill='Control point'), color='black') +
          scale_fill_manual(values=c('Observation'='forestgreen',
                                     'Control point'='darkred')) +
          labs(x=NULL, y='N', color=NULL, fill='Observation type') +
          ggtitle('Heat flow observations') + p_theme +
          theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
        p1 <-
          interp_diff_summary %>%
          filter(short_name %in% x & !(short_name %in% outliers)) %>%
          ggplot() +
          geom_hline(yintercept=0) +
          geom_crossbar(aes(short_name, mean, ymin=mean - 2 * sigma,
                            ymax=mean + 2 * sigma)) +
          labs(x=NULL, y=bquote('Difference'~(mWm^-2)), color=NULL) +
          ggtitle('Point-by-Point Interpolation Differences') + p_theme +
          theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
        p2 <-
          interp_diff_summary %>%
          filter(short_name %in% x & !(short_name %in% outliers)) %>%
          ggplot() +
          geom_col(aes(short_name, 100, fill='Similarity'), color='black') +
          geom_col(aes(short_name, n / n_grid * 100, fill='Krige'), color='black') +
          scale_fill_manual(values=c('Krige'='darkorange', 'Similarity'='navy'),
                            guide='none') +
          labs(x=NULL, y='Coverage (%)', fill='Method') +
          ggtitle('Spatial Coverage') + p_theme +
          theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                legend.position='none')
        p3 <-
          df_acc_rmse %>%
          ggplot() +
          geom_hline(yintercept=0) +
          geom_col(aes(short_name, rmse, fill=method), color='black', position='identity') +
          scale_fill_manual(values=c('Krige'='darkorange', 'Similarity'='navy')) +
          labs(x=NULL, y=bquote('RMSE'~(mWm^-2)), fill='Method') +
          ggtitle('Interpolation Accuracies') + p_theme
        p4 <- p0 / p1 / p2 / p3 + plot_layout(guides='collect') &
          theme(legend.position='bottom')
        ggsave(file=y, plot=p4, width=13, height=10, dpi=300, bg='white')
      })})
    }, error=function(e) {
      cat('\n!! ERROR occurred in plot_interp_accuracy_summary:\n!!', conditionMessage(e))
    })
  }
  x <- nlopt_summary$short_name[str_detect(nlopt_summary$short_name, submap_zone)]
  if (submap_zone %in% c('SAM', 'SEA')) {
    midpoint <- length(x) %/% 2
    walk2(list(x[1:midpoint], x[(midpoint + 1):length(x)]), c(1, 2), ~{
      f(.x, paste0(fig_dir, submap_zone, '-interpolation-accuracy-summary', .y, '.png'))
    })
  } else {
    f(x, paste0(fig_dir, submap_zone, '-interpolation-accuracy-summary.png'))
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot nlopt summary !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_nlopt_summary <- function(base_size=14) {
  load_nlopt_interpolation_data('assets/nlopt_data/interpolation-summary.RData')
  fig_dir <- 'figs/summary/'
  fig_path <- paste0(fig_dir, 'nlopt-summary.png')
  if (!dir.exists(fig_dir)) {dir.create(fig_dir, recursive=T, showWarnings=F)}
  if (file.exists(fig_path)) {
    cat('\n   ', fig_path, ' already exists ...', sep='')
    return(invisible())
  }
  tryCatch({
    cat('\n   Plotting:', fig_path)
    outliers_diff <-
      interp_diff_summary$short_name[which(interp_diff_summary$mean + 2 *
                                           interp_diff_summary$sigma > 300)]
    outliers_acc <-
      interp_accuracy_summary$short_name[which(interp_accuracy_summary$rmse_obs_sim > 300 |
                                               interp_accuracy_summary$rmse_obs_krg > 300)]
    outliers <- unique(c(outliers_diff, outliers_acc))
    suppressWarnings({suppressMessages({
      p1 <-
        nlopt_summary %>%
        filter(!(short_name %in% outliers)) %>%
        mutate(zone=str_sub(short_name, 1, 3)) %>%
        select(-c(vgrm_wt, vgrm_cost, vgrm_rmse, cv_wt, cv_cost, cv_rmse)) %>%
        pivot_longer(-c(short_name, zone, v_mod, cost)) %>%
        ggplot() +
        geom_point(aes(value, cost, fill=zone), shape=21, color='black') +
        labs(x=NULL, y='Cost', color='Zone') +
        scale_fill_manual(values=c('NPA'='darkorange', 'SAM'='navy', 'SEA'='darkred',
                                   'SWP'='forestgreen')) +
        facet_wrap(~name, scales='free_x') +
        theme_bw(base_size=14) +
        theme(panel.grid=element_blank(), panel.background=element_rect(fill='grey90'),
              strip.background=element_blank())
      ggsave(file=fig_path, plot=p1, width=6.5, height=4, dpi=300, bg='white')
    })})
  }, error=function(e) {
    cat('\n!! ERROR occurred in plot_nlopt_summary:\n!!', conditionMessage(e))
  })
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot control point summary !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_control_point_summary <- function(base_size=14) {
  load_nlopt_interpolation_data('assets/nlopt_data/interpolation-summary.RData')
  fig_dir <- 'figs/summary/'
  fig_path <- paste0(fig_dir, 'control-point-summary.png')
  if (!dir.exists(fig_dir)) {dir.create(fig_dir, recursive=T, showWarnings=F)}
  if (file.exists(fig_path)) {
    cat('\n   ', fig_path, ' already exists ...', sep='')
    return(invisible())
  }
  p_theme <- list(
    theme_bw(base_size=14),
    theme(panel.grid=element_blank(), panel.background=element_rect(fill='grey90'),
          plot.title=element_text(hjust=0.5)))
  df <-
    interp_accuracy_summary %>%
    filter(!is.na(rmse_obs_krg)) %>%
    pivot_longer(-c(short_name)) %>%
    mutate(method=ifelse(str_detect(name, 'sim'), 'sim', 'krg'),
           name=str_split(name, '_', simplify=T)[,1],
           zone=str_sub(short_name, 1, 3)) %>%
    rename(metric=name) %>%
    mutate(method=ifelse(method == 'sim', 'Similarity', 'Krige')) %>%
    filter(metric == 'rmse') %>% rename(rmse=value) %>% select(-metric) %>%
    left_join(select(nlopt_summary, short_name, n_obs), by='short_name') %>%
    left_join(select(interp_accuracy_summary, short_name, n_krg), by='short_name') %>%
    rename(n_control=n_krg) %>%
    left_join(select(interp_diff_summary, short_name, n_grid), by='short_name') %>%
    left_join(select(interp_diff_summary, short_name, n), by='short_name') %>%
    mutate(n=ifelse(method == 'Similarity', n_grid, n)) %>%
    rename(n_itp=n) %>% mutate(coverage=n_control / n_grid * 100)
  tryCatch({
    cat('\n   Plotting:', fig_path)
    suppressWarnings({suppressMessages({
      p1 <-
        ggplot(df) +
        geom_point(aes(coverage, rmse, fill=zone), shape=21, color='black') +
        scale_fill_manual(values=c('NPA'='darkorange', 'SAM'='navy', 'SEA'='darkred',
                                   'SWP'='forestgreen')) +
        facet_wrap(~method, scales='free_x') +
        labs(x='Control coverage (%)', y='RMSE', fill='Zone') + p_theme +
        theme(axis.title.x=element_blank())
      p2 <-
        ggplot(filter(df, coverage <= 20)) +
        geom_point(aes(coverage, rmse, fill=zone), shape=21, color='black') +
        scale_fill_manual(values=c('NPA'='darkorange', 'SAM'='navy', 'SEA'='darkred',
                                   'SWP'='forestgreen'), guide='none') +
        facet_wrap(~method, scales='free_x') +
        labs(x='Control coverage (%)', y='RMSE', fill='Zone') + p_theme +
        theme(axis.title.x=element_blank())
      p3 <-
        ggplot(filter(df, coverage <= 5)) +
        geom_point(aes(coverage, rmse, fill=zone), shape=21, color='black') +
        scale_fill_manual(values=c('NPA'='darkorange', 'SAM'='navy', 'SEA'='darkred',
                                   'SWP'='forestgreen'), guide='none') +
        facet_wrap(~method, scales='free_x') +
        labs(x='Control coverage (%)', y='RMSE', fill='Zone') + p_theme
      p4 <- p1 / p2 / p3 + plot_layout(guides='collect')
      ggsave(file=fig_path, plot=p4, width=6.5, height=8, dpi=300, bg='white')
    })})
  }, error=function(e) {
    cat('\n!! ERROR occurred in plot_control_point_summary:\n!!', conditionMessage(e))
  })
}
