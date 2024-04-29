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
sf_use_s2(F)

# Set seed
set.seed(42)

# Set map projections
wgs <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +no_defs'
prj <- '+proj=eck3 +lon_0=-180 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs'

#######################################################
## .1.         General Helper Functions          !!! ##
#######################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get RData object !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
extract_RData_object <- function(file, object) {
  E <- new.env()
  load(file=file, envir=E)
  return(get(object, envir=E, inherits=F))
}

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
    cat('\nAn error occurred in get_world_bathy:\n', conditionMessage(e))
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
    cat('\nAn error occurred in get_seg_bathy:\n', conditionMessage(e))
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
# compile transect data !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
compile_transect_data <- function(trans_ids=NULL, lbuff=5e5, sbuff=5e4, fv=NULL, np=NULL) {
  if (is.null(trans_ids)) {stop('\nMissing submap transect ids!')}
  if (!exists('shp_submap', envir=.GlobalEnv)) {
    stop('\nMissing map data! Use "load(path/to/map-data.RData)"')
  }
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
      mutate(small_buffer=st_buffer(transect, sbuff, endCapStyle='FLAT')) %>%
      mutate(bbox=bbox_extend(st_bbox(large_buffer), square=T)) %>%
      mutate(ridge=crop_feature(shp_ridge, bbox)) %>%
      mutate(trench=crop_feature(shp_trench, bbox)) %>%
      mutate(transform=crop_feature(shp_transform, bbox)) %>%
      mutate(volcano=crop_feature(select(shp_volc, geometry), bbox)) %>%
      mutate(grid=list(crop_feature(shp_grid, large_buffer, T, T))) %>%
      mutate(ghfdb_large_buff=list(crop_feature(shp_ghfdb, large_buffer, T, T))) %>%
      mutate(ghfdb_small_buff=list(crop_feature(shp_ghfdb, small_buffer, T, T))) %>%
      mutate(ghfdb_projected=list(project_obs_to_transect(transect, ghfdb_small_buff))) %>%
      mutate(sim_large_buff=list(crop_feature(shp_sim, large_buffer, T, T))) %>%
      mutate(sim_small_buff=list(crop_feature(shp_sim, small_buffer, T, T))) %>%
      mutate(sim_projected=list(project_obs_to_transect(transect, sim_small_buff))) %>%
      mutate(bathy=list(get_seg_bathy(bbox)))
  })})
  if (!is.null(fv) && !is.null(np)) {
    shp_krg <- Krige(df$ghfdb_large_buff[[1]], fv, df$grid[[1]], np)
    df %>%
      mutate(krg_large_buff=list(crop_feature(shp_krg, large_buffer, T, T))) %>%
      mutate(krg_small_buff=list(crop_feature(shp_krg, small_buffer, T, T))) %>%
      mutate(krg_projected=list(project_obs_to_transect(transect, krg_small_buff))) %>%
      mutate(dff_large_buff=list(interp_diff(krg_large_buff, sim_large_buff))) %>%
      mutate(dff_small_buff=list(interp_diff(krg_small_buff, sim_small_buff))) %>%
      mutate(dff_projected=list(interp_diff(krg_projected, sim_projected))) %>%
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
# cost function !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cost_function <- function(shp_hf=NULL, cutoff=3, n_lags=50, n_max=10, max_dist=1e5,
                          v_mod='Sph', n_fold=NULL, interp_weight=0.5, vgrm_weight=0.5,
                          trans_id=NULL) {
  if (is.null(shp_hf)) {stop('\nMissing heat flow data model!')}
  if (is.null(n_fold)) {n_fold <- nrow(shp_hf)}
  if (!is.null(n_fold)) {if (n_fold > 0 & n_fold <= 1) {n_fold <- nrow(shp_hf) * n_fold}}
  suppressWarnings({
    tryCatch({
      ev <- experimental_vgrm(shp_hf, cutoff, n_lags)
    }, error=function(e) {
      cat('\nAn error occurred in experimental_vgrm:\n', conditionMessage(e))
      return(3)
    })
    if (nrow(ev) < 2) {return(3)}
    tryCatch({
      fv <- fit.variogram(ev, vgm(model=v_mod), fit.method=6)
    }, error=function(e) {
      cat('\nAn error occurred in fit.variogram:\n', conditionMessage(e))
      return(3)
    })
    if (fv$range < 0) {
      cat('\nVariogram range is negative:', fv$range)
      return(3)
    }
    tryCatch({
      k_cv <- krige.cv(obs~1, shp_hf, model=fv, nmax=n_max, maxdist=max_dist, nfold=n_fold)
    }, error=function(e) {
      cat('\nAn error occurred in krige.cv:\n', conditionMessage(e))
      return(3)
    })
    na_sum <- sum(is.na(k_cv$residual))
    na_thresh <- nrow(k_cv) / 8
    if (na_sum != 0) {
      if (na_sum >= na_thresh) {
        cat('\nToo many NAs in krige.cv:', na_sum, '/', na_thresh)
        return(3)
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
      cat('\nAn error occurred in cost_function:\n', conditionMessage(e))
      return(3)
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
      stop('\nAn error occurred in experimental_vgrm:\n', conditionMessage(e))
    })
    if (nrow(ev) < 2) {stop('\nExperimental variogram has less than two lags!')}
    if (any(class(ev) == 'try-error')) {stop('\nExperimental variogram error!')}
    tryCatch({
      fv <- fit.variogram(ev, vgm(model=v_mod), fit.method=6)
    }, error=function(e) {
      stop('\nAn error occurred in fit.variogram:\n', conditionMessage(e))
    })
  })
  return(list('experimental_vgrm'=as_tibble(ev), 'fitted_vgrm'=fv))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# nlopt krige !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nlopt_krige <- function(trans_id=NULL, v_mod='Sph', alg='NLOPT_LN_COBYLA', max_eval=500,
                        n_fold=10, iwt=0.5, vwt=0.5) {
  if (is.null(trans_id)) {stop('\nMissing submap transect id!')}
  nlopt_dir <- 'assets/nlopt_data/nlopt_transects'
  nlopt_id <- paste0('opt-', trans_id, '-', v_mod) 
  nlopt_path <- paste0(nlopt_dir, '/', nlopt_id, '.RData')
  nlopt_itr_path <- paste0('assets/nlopt_data/nlopt_itr/nlopt-out-', trans_id, '-', v_mod)
  if (file.exists(nlopt_path) & file.exists(nlopt_itr_path)) {
    cat('\nOptimized', v_mod, 'kriging model found for submap transect:', trans_id)
    return(invisible())
  }
  if (!file.exists(nlopt_path) & file.exists(nlopt_itr_path)) {
    cat('\nOptimized', v_mod, 'kriging model failed for submap transect:', trans_id)
    cat('\nRemove', nlopt_itr_path, 'and run "make nlopt" again ...')
    return(invisible())
  }
  if (!dir.exists(nlopt_dir)) {dir.create(nlopt_dir, recursive=T, showWarnings=F)}
  x0 <- c(3, 50, 10, 1e5) # Initial values (cutoff, n_lags, n_max, max_dist)
  lb <- c(1, 30, 2, 5e4) # Lower bound (cutoff, n_lags, n_max, max_dist)
  ub <- c(12, 100, 50, 5e5) # Upper bound (cutoff, n_lags, n_max, max_dist)
  opts <-
    list(print_level=0, maxeval=max_eval, algorithm=alg, xtol_rel=1e-15, ftol_rel=1e-10)
  x <- compile_transect_data(trans_id)
  obs <- x$ghfdb_large_buff[[1]]
  nlopt_fun <- function(x) {
    cost_function(obs, x[1], x[2], x[3], x[4], v_mod, n_fold, iwt, vwt, trans_id)
  }
  tryCatch({
    cat('\nOptimizing', v_mod, 'kriging model for submap transect:', trans_id)
    opt <- nloptr(x0, nlopt_fun, lb=lb, ub=ub, opts=opts)
  }, error=function(e) {
    cat('\nAn error occurred in nlopt_krige:\n', conditionMessage(e))
  })
  if (opt$status != 0) {
    cat('\n', rep('+', 60), sep='')
    cat('\nNlopt failed to converge:')
    cat('\n', rep('-', 60), sep='')
    cat('\nNloptr status:', opt$status)
    cat('\nNloptr iterations:', opt$iterations)
    cat('\nNloptr message:\n', opt$message)
    cat('\n', rep('+', 60), sep='')
    return(invisible())
  }
  opt_decoded <- decode_opt(obs, v_mod, opt)
  assign(str_replace_all(nlopt_id, '-', '_'), opt_decoded)
  save(list=str_replace_all(nlopt_id, '-', '_'), file=nlopt_path)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# nlopt transects !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nlopt_transects <- function(trans_ids=NULL, v_mods=c('Sph', 'Exp', 'Lin', 'Bes'),
                            alg='NLOPT_LN_COBYLA', max_eval=500, n_fold=10,
                            iwt=0.5, vwt=0.5, parallel=T) {
  if (is.null(trans_ids)) {stop('\nMissing submap transect ids!')}
  if (length(list.files('assets/map_data/relief')) < length(trans_ids)) {parallel <- F}
  x <- expand.grid(id=trans_ids, vm=v_mods, stringsAsFactors=F) %>% arrange(id, vm)
  if (parallel) {
    plan(multicore, workers=availableCores() - 2)
    future_walk2(x$id, x$vm, ~nlopt_krige(..., alg, max_eval, n_fold, iwt, vwt),
                 .options=furrr_options(seed=42), .progress=T)
  } else {
    walk2(x$id, x$vm, ~nlopt_krige(..., alg, max_eval, n_fold, iwt, vwt))
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
    cat('\nNo nlopt itr files found for:', trans_id)
    return(NULL)
  }
  if (length(nlopt_itr_paths) != length(nlopt_trans_paths)) {
    itr_mods <- str_sub(nlopt_itr_paths, start=-3)
    trans_mods <- str_extract(nlopt_trans_paths, ".{3}(?=\\.RData)")
    nlopt_itr_paths <-
      nlopt_itr_paths[str_detect(nlopt_itr_paths, paste(trans_mods, collapse="|"))]
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
    cat('\nNo nlopt RData found for:', trans_id)
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
      cat('\nAn error occurred in Krige:\n', conditionMessage(e))
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
  if (parallel) {
    plan(multicore, workers=availableCores() - 2)
    as_tibble(future_map_dfr(trans_ids, f, .options=furrr_options(seed=42), .progress=T))
  } else {
    map_df(trans_ids, f)
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# summarize interpolation differences !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summarize_interpolation_differences <- function(trans_ids, parallel=T) {
  f <- function(x) {
    tryCatch({
      km <- get_optimal_krige_model(x)
      fv <- km$fitted_vgrm
      np <- km$opt_krige_mod_summary$n_pairs
      compile_transect_data(x, fv=fv, np=np)$dff_large_buff[[1]] %>%
        st_set_geometry(NULL) %>%
        summarise(short_name=x, n=sum(!is.na(est_dff)), rmse=sqrt(mean(est_dff^2, na.rm=T)),
                  min=min(est_dff, na.rm=T), max=max(est_dff, na.rm=T),
                  med=median(est_dff, na.rm=T), iqr=IQR(est_dff, na.rm=T),
                  mean=mean(est_dff, na.rm=T), sigma=sd(est_dff, na.rm=T))
    }, error=function(e) {
      cat('\nAn error occurred in interp_diff_summary:\n', conditionMessage(e))
      return(NULL)
    })
  }
  if (parallel) {
    as_tibble(future_map_dfr(ids, f, .options=furrr_options(seed=42), .progress=T))
  } else {
    map_df(ids, f)
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# summarize interpolation accuracy !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summarize_interpolation_accuracy <- function(trans_ids, parallel=T) {
  f <- function(x) {
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
      cat('\nAn error occurred in interp_accuracy_summary:\n', conditionMessage(e))
      return(NULL)
    })
  }
  if (parallel) {
    as_tibble(future_map_dfr(ids, f, .options=furrr_options(seed=42), .progress=T))
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
plot_ghfdb_base <- function() {
  fig_path <- 'figs/ghfdb.png'
  if (!dir.exists('figs')) {dir.create('figs', recursive=T, showWarnings=F)}
  if (file.exists(fig_path)) {
    cat('\n', fig_path, ' already exists ...', sep='')
    return(invisible())
  }
  cat('\nPlotting:', fig_path)
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
# plot transect !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_transect <- function(trans_id, base_size=20) {
  fig_dir <- 'figs/transect/'
  fig_path <- paste0(fig_dir, trans_id, '-transect.png')
  if (!dir.exists(fig_dir)) {dir.create(fig_dir, recursive=T, showWarnings=F)}
  if (file.exists(fig_path)) {
    cat('\n', fig_path, ' already exists ...', sep='')
    return(invisible())
  }
  tryCatch({
    cat('\nPlotting:', fig_path)
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
              legend.justification='right', legend.position=c(0.92, 0.85),
              legend.direction='horizontal', legend.key.height=unit(0.5, 'cm'),
              legend.key.width=unit(1, 'cm'), legend.box.margin=margin(2, 2, 2, 2),
              legend.margin=margin(), legend.title=element_text(vjust=0, size=base_size),
              legend.background=element_blank())
      )
    suppressWarnings({suppressMessages({
      p1 <-
        ggplot(x) +
        geom_sf(data=x$bathy[[1]], aes(color=elev), size=0.5, shape=15) +
        scale_color_etopo(guide='none') +
        new_scale_color() +
        geom_sf(aes(geometry=large_buffer), fill=NA, linewidth=0.5) +
        geom_sf(aes(geometry=small_buffer), fill=NA, linewidth=0.5) +
        geom_sf(aes(geometry=ridge), color='white') +
        geom_sf(aes(geometry=transform), color='white') +
        geom_sf(aes(geometry=trench), color='white', linewidth=1.5) +
        geom_sf(aes(geometry=transect), color='black', linewidth=1.5) +
        geom_sf(aes(geometry=volcano), color='black', fill='white', shape=24) +
        geom_sf(data=x$ghfdb_large_buff[[1]], aes(color=obs), shape=20, size=hf_pt_size) +
        ggtitle('Global HF Observations') +
        scale_color_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                              breaks=c(0, 125, 250), na.value='transparent', guide='none') +
        coord_sf(expand=F, lims_method='geometry_bbox') + map_theme
      if (is.null(x$ghfdb_projected[[1]])) {
        p2 <- ggplot(data=data.frame(), aes())
      } else {
        p2 <-
          ggplot(x$ghfdb_projected[[1]]) +
          geom_smooth(aes(projected_distances, obs), method='gam', color='black') +
          geom_point(aes(projected_distances, obs), shape=20, color='grey20', size=3) +
          geom_rug(aes(projected_distances, obs, color=obs), length=unit(0.06, "npc"))
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
        geom_sf(aes(geometry=small_buffer), fill=NA, linewidth=0.5) +
        geom_sf(aes(geometry=ridge), color='white') +
        geom_sf(aes(geometry=transform), color='white') +
        geom_sf(aes(geometry=trench), color='white', linewidth=1.5) +
        geom_sf(aes(geometry=transect), color='black', linewidth=1.5) +
        geom_sf(aes(geometry=volcano), color='black', fill='white', shape=24) +
        geom_sf(data=x$sim_large_buff[[1]], aes(color=est_sim), shape=20, size=hf_pt_size) +
        geom_sf(data=comps_sim[!is.na(comps_sim$est_sim),], aes(fill=est_sim), shape=21,
                size=hf_pt_size, stroke=pt_stroke, color='white') +
        ggtitle('Similarity Interpolation') +
        scale_color_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                              breaks=c(0, 125, 250), na.value='transparent', guide='none') +
        scale_fill_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                              breaks=c(0, 125, 250), na.value='transparent', guide='none') +
        coord_sf(expand=F, lims_method='geometry_bbox') + map_theme
      if (is.null(x$sim_projected[[1]])) {
        p4 <- ggplot(data=data.frame(), aes())
      } else {
        p4 <-
          ggplot(x$sim_projected[[1]]) +
          geom_smooth(aes(projected_distances, obs), method='gam', color='black') +
          geom_point(aes(projected_distances, obs), shape=20, color='grey20', size=3) +
          geom_rug(aes(projected_distances, obs, color=obs), length=unit(0.06, "npc"))
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
        geom_sf(aes(geometry=small_buffer), fill=NA, linewidth=0.5) +
        geom_sf(aes(geometry=ridge), color='white') +
        geom_sf(aes(geometry=transform), color='white') +
        geom_sf(aes(geometry=trench), color='white', linewidth=1.5) +
        geom_sf(aes(geometry=transect), color='black', linewidth=1.5) +
        geom_sf(aes(geometry=volcano), color='black', fill='white', shape=24) +
        geom_sf(data=x$krg_large_buff[[1]], aes(color=est_krg), shape=20, size=hf_pt_size) +
        geom_sf(data=comps_krg[!is.na(comps_krg$est_krg),], aes(fill=est_krg), shape=21,
                size=hf_pt_size, stroke=pt_stroke, color='white') +
        ggtitle('Krige Interpolation') +
        scale_color_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                              breaks=c(0, 125, 250), na.value='transparent', guide='none') +
        scale_fill_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                              breaks=c(0, 125, 250), na.value='transparent', guide='none') +
        coord_sf(expand=F, lims_method='geometry_bbox') + map_theme
      if (is.null(x$krg_projected[[1]])) {
        p6 <- ggplot(data=data.frame(), aes())
      } else {
        p6 <-
          ggplot(x$krg_projected[[1]]) +
          geom_smooth(aes(projected_distances, obs), method='gam', color='black') +
          geom_point(aes(projected_distances, obs), shape=20, color='grey20', size=3) +
          geom_rug(aes(projected_distances, obs, color=obs), length=unit(0.06, "npc"))
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
      p_title <- paste0('Submap Transect: ', trans_id, ' ', x$trench_name)
      p7 <-
        (p1 + p5 + p3) / (p2 + p6 + p4) +
        plot_layout(widths=1, heights=c(1.5, 1)) +
        plot_annotation(title=p_title, tag_levels = list(c('a)', 'b)', 'c)', '', '', '')),
                        theme=theme(plot.title=element_text(size=base_size * 1.5,
                                                            margin=margin()))) &
        theme(plot.tag=element_text(size=base_size * 1.5, margin=margin(0, 0, -45, 0)))
      ggsave(file=fig_path, plot=p7, width=19.5, height=11, dpi=300, bg='white')
    })})
  }, error=function(e) {
    cat('\nAn error occurred in plot_transect:\n', conditionMessage(e))
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
    cat('\n', fig_path, ' already exists ...', sep='')
    return(invisible())
  }
  tryCatch({
    cat('\nPlotting:', fig_path)
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
    p0 <-
      ggplot(filter(nlopt_itr, v_mod == opt_kmod$v_mod)) +
      geom_path(aes(itr, max_dist / 1e3), color='forestgreen', linewidth=1) +
      geom_point(data=opt_kmod, aes(itr, max_dist / 1e3), color='forestgreen', size=3) +
      geom_label_repel(data=opt_kmod, aes(itr, max_dist / 1e3, label=round(max_dist / 1e3)),
                       size=base_size * 0.28) +
      annotate('text', x=Inf, y=Inf, label='Variogram Max Distance', hjust=1.05, vjust=1.5,
               size=base_size * 0.35) +
      labs(x='Nlopt Iteration', y='Distance (km)', color=NULL) + p_theme
    p1 <-
      ggplot(filter(nlopt_itr, v_mod == opt_kmod$v_mod)) +
      geom_path(aes(itr, n_pairs), color='firebrick', linewidth=1) +
      geom_point(data=opt_kmod, aes(itr, n_pairs), color='firebrick', size=3) +
      geom_label_repel(data=opt_kmod, aes(itr, n_pairs, label=round(n_pairs)),
                       size=base_size * 0.28) +
      annotate('text', x=Inf, y=Inf, label='Variogram Pairs', hjust=1.05, vjust=1.5,
               size=base_size * 0.35) +
      labs(x='Nlopt Iteration', y='Pairs', color=NULL) + p_theme
    p2 <-
      ggplot(filter(nlopt_itr, v_mod == opt_kmod$v_mod)) +
      geom_path(aes(itr, cutoff), color='orchid4', linewidth=1) +
      geom_point(data=opt_kmod, aes(itr, cutoff), color='orchid4', size=3) +
      geom_label_repel(data=opt_kmod, aes(itr, cutoff, label=round(cutoff)),
                       size=base_size * 0.28) +
      annotate('text', x=Inf, y=Inf, label='Variogram Cutoff', hjust=1.05, vjust=1.5,
               size=base_size * 0.35) +
      labs(x='Nlopt Iteration', y='Cutoff', color=NULL) + p_theme
    p3 <-
      ggplot(filter(nlopt_itr, v_mod == opt_kmod$v_mod)) +
      geom_path(aes(itr, n_lags), color='saddlebrown', linewidth=1) +
      geom_point(data=opt_kmod, aes(itr, n_lags), color='saddlebrown', size=3) +
      geom_label_repel(data=opt_kmod, aes(itr, n_lags, label=round(n_lags)),
                       size=base_size * 0.28) +
      annotate('text', x=Inf, y=Inf, label='Variogram Lags', hjust=1.05, vjust=1.5,
               size=base_size * 0.35) +
      labs(x='Nlopt Iteration', y='Lags', color=NULL) + p_theme
    p4 <-
      ggplot(filter(nlopt_itr, v_mod == opt_kmod$v_mod)) +
      geom_path(aes(itr, vgrm_cost), color='darkorange', linewidth=1) +
      geom_point(data=opt_kmod, aes(itr, vgrm_cost), color='darkorange', size=3) +
      geom_label_repel(data=opt_kmod, aes(itr, vgrm_cost, label=round(vgrm_cost, 3)),
                       size=base_size * 0.28) +
      geom_path(aes(itr, cv_cost), color='navy', linewidth=1) +
      geom_point(data=opt_kmod, aes(itr, cv_cost), color='navy', shape=20, size=5) +
      geom_label_repel(data=opt_kmod, aes(itr, cv_cost, label=paste0(round(cv_cost, 3))),
                       size=base_size * 0.28) +
      annotate('text', x=Inf, y=Inf, label='Cost Function', hjust=1.05, vjust=1.5,
               size=base_size * 0.35) +
      labs(x='Nlopt Iteration', y='Cost', color=NULL) + p_theme
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
  }, error=function(e) {
    cat('\nAn error occurred in plot_optimal_variogram:\n', conditionMessage(e))
  })
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot interp accuracy summary !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_interp_accuracy_summary <- function(submap_zone=NULL, base_size=22) {
  if (is.null(submap_zone)) {stop('\nMissing submap zone!')}
  fig_dir <- 'figs/summary/'
  data_path <- 'assets/nlopt_data/interpolation-summary.RData'
  if (!dir.exists(fig_dir)) {dir.create(fig_dir, recursive=T, showWarnings=F)}
  if (!exists('interp_accuracy_summary', envir=.GlobalEnv)) {
    stop('\nMissing nlopt summary data! Use "load(path/to/interpolation-summary.RData)"')
  }
  f <- function(x, y) {
    if (file.exists(y)) {
      cat('\n', y, ' already exists ...', sep='')
      return(invisible())
    }
    tryCatch({
      cat('\nPlotting:', y)
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
      p1 <-
        interp_diff_summary %>%
        left_join(select(nlopt_summary, short_name, n_obs), by='short_name') %>%
        filter(short_name %in% x & !(short_name %in% outliers)) %>%
        ggplot() +
        geom_hline(yintercept=0) +
        geom_crossbar(aes(short_name, mean, ymin=mean - 2 * sigma, ymax=mean + 2 * sigma)) +
        labs(x=NULL, y=bquote('Difference'~(mWm^-2))) +
        ggtitle('Point-by-Point Interpolation Differences') + p_theme
      p2 <-
        interp_accuracy_summary %>%
        filter(!is.na(rmse_obs_krg)) %>%
        pivot_longer(-c(short_name)) %>%
        mutate(method=ifelse(str_detect(name, 'sim'), 'sim', 'krg'),
               name=str_split(name, '_', simplify=T)[,1],
               zone=str_sub(short_name, 1, 3)) %>%
        rename(metric=name) %>%
        filter(zone == submap_zone) %>%
        filter(short_name %in% x & !(short_name %in% outliers)) %>%
        mutate(method=ifelse(method == 'sim', 'Similarity', 'Krige')) %>%
        group_by(method) %>%
        filter(metric == 'rmse') %>%
        ggplot() +
        geom_col(aes(short_name, value, fill=method), color='black', position='dodge') +
        scale_fill_manual(values=c('Krige'='darkorange', 'Similarity'='navy')) +
        labs(x=NULL, y=bquote('RMSE'~(mWm^-2)), fill='Method') +
        ggtitle('Interpolation Accuracies') + p_theme
      p3 <- p1 / p2
      ggsave(file=y, plot=p3, width=13, height=6.5, dpi=300, bg='white')
    }, error=function(e) {
      cat('\nAn error occurred in plot_interp_accuracy_summary:\n', conditionMessage(e))
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
  fig_dir <- 'figs/summary/'
  fig_path <- paste0(fig_dir, 'nlopt-summary.png')
  data_path <- 'assets/nlopt_data/interpolation-summary.RData'
  if (!dir.exists(fig_dir)) {dir.create(fig_dir, recursive=T, showWarnings=F)}
  if (!exists('nlopt_summary', envir=.GlobalEnv)) {
    stop('\nMissing nlopt summary data! Use "load(path/to/interpolation-summary.RData)"')
  }
  tryCatch({
    cat('\nPlotting:', fig_path)
    outliers_diff <-
      interp_diff_summary$short_name[which(interp_diff_summary$mean + 2 *
                                           interp_diff_summary$sigma > 300)]
    outliers_acc <-
      interp_accuracy_summary$short_name[which(interp_accuracy_summary$rmse_obs_sim > 300 |
                                               interp_accuracy_summary$rmse_obs_krg > 300)]
    outliers <- unique(c(outliers_diff, outliers_acc))
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
  }, error=function(e) {
    cat('\nAn error occurred in plot_nlopt_summary:\n', conditionMessage(e))
  })
}
