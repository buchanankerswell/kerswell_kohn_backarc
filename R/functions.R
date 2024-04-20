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
package_list <- c('tictoc', 'stringr', 'tidyr', 'readr', 'purrr', 'furrr', 'tibble', 'dplyr',
                  'magrittr', 'ggplot2', 'colorspace', 'metR', 'ggrepel', 'ggridges',
                  'ggnewscale', 'patchwork', 'cowplot', 'ggsflabel', 'marmap', 'scales',
                  'ggspatial', 'gstat', 'rgeos', 'sp', 'sf', 'rnaturalearth', 'nloptr', 'zoo',
                  'jsonlite')

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
    ft_cropped <- st_crop(ft, bbox)
  } else {
    ft_cropped <- st_intersection(ft, bbox)
  }
  if (!is.data.frame(ft_cropped)) {l <- length(ft_cropped)} else {l <- nrow(ft_cropped)}
  if (l == 0) {NA} else {
    if (!keep_df) {st_sfc(st_union(ft_cropped), crs=st_crs(bbox))} else {ft_cropped}
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get world bathy !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_world_bathy <- function(res=15, path='assets/map_data/relief/') {
  if (!dir.exists(path)) {dir.create(path)}
  tryCatch({
    suppressWarnings({suppressMessages({
      getNOAA.bathy(180, -180, 90, -90, res, T, F, path) %>%
        as.SpatialGridDataFrame() %>% st_as_sf(crs=wgs) %>% rename(elev=layer) %>%
        reproject_center_pacific()
    })})
  }, error = function(e) {
    cat('An error occurred in get_world_bathy:', conditionMessage(e), '\n')
    return(NULL)
  })
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get seg bathy !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_seg_bathy <- function(shp, res=2, path='assets/map_data/relief/', tol=1) {
  if (!dir.exists(path)) {dir.create(path)}
  tryCatch({
    bbx <- shp %>% st_transform(wgs) %>% st_bbox() %>% round(2)
    if (bbx[3] - bbx[1] > 180) {antim <- T} else {antim <- F}
    getNOAA.bathy(bbx[3], bbx[1], bbx[2], bbx[4], res, T, antim, path) %>%
      as.SpatialGridDataFrame() %>% st_as_sf(crs=wgs) %>% st_make_valid() %>%
      st_transform(prj) %>% rename(elev=layer)
  }, error = function(e) {
    cat('An error occurred in get_seg_bathy:', conditionMessage(e), '\n')
    return(NULL)
  })
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# handle zerodist obs !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
handle_zerodist_obs <- function(df) {
  dup <- zerodist(as_Spatial(df))
  n_dup <- nrow(dup)
  rid <- map_dbl(1:n_dup, ~ {
    i <- .x
    if (
      df$code6[dup[i, 1]] != df$code6[dup[i, 2]] &
      df$code6[dup[i, 1]] > df$code6[dup[i, 2]]
    ) {
      return(dup[i, 1])
    } else if (
      df$code6[dup[i, 1]] != df$code6[dup[i, 2]] &
      df$code6[dup[i, 1]] < df$code6[dup[i, 2]]
    ) {
      return(dup[i, 2])
    } else {
      return(dup[i, sample(1:2, 1)])
    }
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
# compile transect data !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
compile_transect_data <- function(trans_id, large_buff=5e5, small_buff=5e4) {
  if (!exists('shp_submap', envir=.GlobalEnv)) {
    stop('\nMissing map data! Use "load(path/to/map-data.RData"')
  }
  if (is.numeric(trans_id) && all(trans_id >= 1) && all(trans_id <= nrow(shp_submap))) {
    x <- slice(shp_submap, trans_id)
  } else if (is.character(trans_id) && all(trans_id %in% shp_submap$short_name)) {
    x <- filter(shp_submap, short_name %in% trans_id)
  } else {
    stop('\nUnrecognized input for trans_id! Use the transect short_name or row number ...')
  }
  if (nrow(x) == 0) {
    stop('\n', trans_id, ' not found in shp_submap!')
  }
  cat('\nCompiling transect data for:', trans_id)
  suppressWarnings({suppressMessages({
    x %>% rowwise() %>%
      mutate(large_buffer=st_buffer(transect, large_buff, endCapStyle='ROUND')) %>%
      mutate(small_buffer=st_buffer(transect, small_buff, endCapStyle='FLAT')) %>%
      mutate(bbox=bbox_extend(st_bbox(large_buffer), square=T)) %>%
      mutate(grid=crop_feature(shp_grid, bbox)) %>%
      mutate(ridge=crop_feature(shp_ridge, bbox)) %>%
      mutate(trench=crop_feature(shp_trench, bbox)) %>%
      mutate(transform=crop_feature(shp_transform, bbox)) %>%
      mutate(volcano=crop_feature(select(shp_volc, geometry), bbox)) %>%
      mutate(tglobe_large_buff=list(crop_feature(shp_tglobe, large_buffer, T, T))) %>%
      mutate(tglobe_small_buff=list(crop_feature(shp_tglobe, small_buffer, T, T))) %>%
      mutate(tglobe_projected=list(project_obs_to_transect(transect, tglobe_small_buff))) %>%
      mutate(sim=list(crop_feature(shp_sim, large_buffer, T, T))) %>%
      mutate(bathy=list(get_seg_bathy(bbox))) %>%
      ungroup()
  })})
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# project obs to transect !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
project_obs_to_transect <- function(transect, tglobe_small_buff) {
  if (st_crs(transect) != st_crs(tglobe_small_buff)) {
    stop('transect and tglobe_small_buff crs not the same!')
  }
  if (!any(class(transect) == 'sfc')) {
    stop('\ntransect needs to be an sf object!')
  }
  if (!any(class(tglobe_small_buff) == 'sf')) {
    stop('\nUnrecognized tglobe_small_buff passed to project_obs_to_transect() !')
  }
  if (!any(names(tglobe_small_buff) == 'obs')) {
    stop('\nUnrecognized sf object passed to project_obs_to_transect() !')
  }
  projected_distances <- st_line_project(transect, tglobe_small_buff$tglobe, normalized=T)
  st_as_sf(st_line_interpolate(transect, projected_distances)) %>%
    mutate(projected_distances=projected_distances, obs=tglobe_small_buff$obs)
}


#######################################################
## .2.             Kriging Functions             !!! ##
#######################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# parse krige args !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
parse_krige_args <- function(args) {
  if (length(args) == 0) {
    max_eval <- 30
    alg <- 'NLOPT_LN_SBPLX'
    iwt <- 0.5
    vwt <- 0.5
    n_fold <- NULL
    n_cores <- availableCores() - 2
  } else if (length(args) != 0) {
    max_eval <- suppressWarnings(as.numeric(args[1]))
    if (is.na(max_eval)) {
      max_eval <- 50
    } else if (max_eval <= 0) {
      max_eval <- 50
    }
    alg <- suppressWarnings(as.numeric(args[2]))
    if (is.na(alg)) {
      alg <- 'NLOPT_LN_SBPLX'
    } else {
      if (args[2] == 1) {
        alg <- 'NLOPT_GN_DIRECT'
      } else if (args[2] == 2) {
        alg <- 'NLOPT_GN_DIRECT_L'
      } else if (args[2] == 3) {
        alg <- 'NLOPT_LN_SBPLX'
      } else if (args[2] == 4) {
        alg <- 'NLOPT_LN_NELDERMEAD'
      } else if (args[2] == 5) {
        alg <- 'NLOPT_LN_BOBYQA'
      } else if (args[2] == 6) {
        alg <- 'NLOPT_LN_COBYLA'
      } else {
        alg <- 'NLOPT_LN_COBYLA'
      }
    }
    iwt <- suppressWarnings(as.numeric(args[3]))
    if (is.na(iwt)) {
      iwt <- 0.5
    }
    vwt <- suppressWarnings(as.numeric(args[4]))
    if (is.na(vwt)) {
      vwt <- 0.5
    }
    if ((iwt + vwt) != 1) {
      iwt <- 0.5
      vwt <- 0.5
    }
    n_fold <- suppressWarnings(as.numeric(args[5]))
    if (is.na(n_fold)) {
      n_fold <- NULL
    } else {
      if (n_fold < 0) {
        n_fold <- NULL
      } else if (n_fold == 0) {
        n_fold <- NULL
      }
    }
    n_cores <- suppressWarnings(as.numeric(args[6]))
    if (is.na(n_cores)) {
      n_cores <- availableCores() - 2
    } else if (n_cores > availableCores()) {
      n_cores <- availableCores() - 2
    } else if (n_cores < 0) {
      n_cores <- availableCores() - 2
    }
  }
  # Load variables into current env
  assign('max_eval', max_eval, envir=.GlobalEnv)
  assign('alg', alg, envir=.GlobalEnv)
  assign('iwt', iwt, envir=.GlobalEnv)
  assign('vwt', vwt, envir=.GlobalEnv)
  assign('n_fold', n_fold, envir=.GlobalEnv)
  assign('n_cores', n_cores, envir=.GlobalEnv)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# optimize krige !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Minimization function
optimize_krige <- function(segment, v_mod, hf, alg, max_eval, n_fold, iwt, vwt) {
  # Initial values (cutoff, n_lags, lag_start, n_max)
  x0 <- c(3, 50, 3, 10)
  lb <- c(1, 30, 1, 2)
  ub <- c(12, 100, 10, 50)
  opts <- list(print_level=0, maxeval=max_eval, algorithm=alg, ftol_rel=1e-5)
  hf_obs <- shp_hf_crop[[segment]]
  opt_fun <- function(x) {
    cost_function(hf_obs, x[1], x[2], x[3], x[4], v_mod, n_fold, iwt, vwt, segment)
  }
  nloptr(x0, opt_fun, lb=lb, ub=ub, opts=opts)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# experimental vgrm !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
experimental_vgrm <- function(shp_hf=NULL, cutoff=3, n_lags=30, lag_start=1) {
  if (is.null(shp_hf)) {stop('\nMissing heat flow data!')}
  bbox <- st_bbox(shp_hf)
  bbox_diagonal_distance <- sqrt((bbox$xmax - bbox$xmin)^2 + (bbox$ymax - bbox$ymin)^2)
  lag_cutoff <- as.vector(bbox_diagonal_distance / cutoff)
  bin_width <- lag_cutoff / n_lags
  shifted_cutoff <- bin_width * (n_lags + lag_start - 1)
  variogram(hf~1, locations=shp_hf, cutoff=shifted_cutoff, width=bin_width) %>%
    slice(floor(lag_start):n())
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# cost function !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cost_function <- function(shp_hf, cutoff=3, n_lags=50, lag_start=3, n_max=10,
                          model_vgrm='Sph', n_fold=NULL, interp_weight=0.5, vgrm_weight=0.5,
                          segment=NULL, verbose=T) {
  if (is.null(shp_hf)) {stop('\nMissing heat flow data model!')}
  if (is.null(n_fold)) {n_fold <- nrow(shp_hf)}
  if (n_fold > 0 & n_fold <= 1) {n_fold <- nrow(shp_hf) * n_fold}
  suppressWarnings({
    experimental_vgrm <- try(experimental_vgrm(shp_hf, cutoff, n_lags, lag_start))
  })
  if (nrow(experimental_vgrm) < 2) {
    if (verbose) {
      cat('\nExperimental variogram has less than two lags!')
      cat('\nReturning arbitrarily high cost')
      return(runif(1, 1, 1.5))
    }
  }
  if (any(class(experimental_vgrm) == 'try-error')) {
    if (verbose) {
      cat('\nExperimental variogram error!')
      cat('\nReturning arbitrarily high cost')
      return(runif(1, 1, 1.5))
    }
  }
  suppressWarnings({
    fitted_vgrm <- try(fit.variogram(experimental_vgrm, vgm(model=model_vgrm), fit.method=6))
  })
  if (any(class(fitted_vgrm) == 'try-error')) {
    if (verbose) {
      cat('\nVariogram fitting error!\n')
      cat('\nReturning arbitrarily high cost')
      return(runif(1, 1, 1.5))
    }
  }
  if (verbose) {
    cat('\n', rep('+', 30), sep='')
    cat('\nComputing cross-validation')
  }
  suppressWarnings({
    k_cv <- try(krige.cv(hf~1, shp_hf, model=fitted_vgrm, nmax=n_max, nfold=n_fold))
  })
  if (any(class(k_cv) == 'try-error')) {
    if (verbose) {
      cat('\nCross-validation error!')
      cat('\nReturning arbitrarily high cost\n')
      return(runif(1, 1, 1.5))
    }
  }
  if (sum(is.na(k_cv)) != 0) {
    if (verbose) {
      cat('\nCross-validation produced NAs!')
      if (sum(is.na(k_cv$residual)) >= nrow(k_cv) / 2) {
        cat('\nCross-validation produced too many NAs!')
        cat('\nReturning arbitrarily high cost\n')
        return(runif(1, 1, 1.5))
      } else {
        cat('\nComputing cost despite', sum(is.na(k_cv$residual)), '/', nrow(k_cv), 'NAs')
      }
    }
  }
  k_cv <- k_cv %>% filter(!is.na(residual))
  # Calculating cost after Li et al. (2018)
  suppressWarnings({
    vgrm_rmse <- sqrt(sqrt(attr(fitted_vgrm, 'SSErr') / nrow(experimental_vgrm)))
    vgrm_sd <- sqrt(sd(experimental_vgrm$gamma, na.rm=T))
    vgrm_cost <- vgrm_weight * vgrm_rmse / vgrm_sd
    interp_rmse <- sqrt(sum(k_cv$residual^2, na.rm=T) / nrow(k_cv))
    interp_sd <- sd(k_cv$var1.pred, na.rm=T)
    interp_cost <- interp_weight * interp_rmse / interp_sd
  })
  if (verbose) {
    sink(paste0('log/nlopt-out-', format(Sys.Date(), '%d-%m-%Y')), append=T, split=T)
    cat('\nSegment:', segment)
    cat('\nCutoff proportion:', cutoff)
    cat('\nNumber of lags:', n_lags)
    cat('\nLag start:', lag_start)
    cat('\nMax pairs:', n_max)
    cat('\nVariogram weight:', vgrm_weight)
    cat('\nVariogram rmse:', vgrm_rmse)
    cat('\nVariogram sd:', vgrm_sd)
    cat('\nVariogram cost:', vgrm_cost)
    cat('\nInterpolation weight:', interp_weight)
    cat('\nInterpolation rmse:', interp_rmse)
    cat('\nInterpolation sd:', interp_sd)
    cat('\nInterpolation cost:', interp_cost)
    cat('\nCost:', vgrm_cost + interp_cost)
    cat('\nVariogram model:', model_vgrm)
    cat('\n')
    print(fitted_vgrm)
    cat(rep('+', 30), sep='')
    cat('\n', rep('-', 40), sep='')
    sink()
  }
  return(vgrm_cost + interp_cost)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read nloptr trace !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
read_nloptr_trace <- function(fpath) {
  read_file <- function(path, ...) {
    con <- file(path)
    on.exit(close(con))
    suppressWarnings(readLines(con, ...))
  }
  t <- read_file(fpath)
  segs_t <- t[grepl('^Segment: ', t)] %>% map_chr(~gsub('Segment: ', '', .x))
  vgrm_model_t <- t[grepl('^Variogram model:', t)] %>%
    map_chr(~gsub('Variogram model: ', '', .x))
  vgrm_weight_t <- t[grepl('^Variogram weight', t)] %>%
    map_dbl(~as.numeric(gsub('Variogram weight: ', '', .x)))
  vgrm_rmse_t <- t[grepl('^Variogram rmse', t)] %>%
    map_dbl(~as.numeric(gsub('Variogram rmse: ', '', .x)))
  vgrm_cost_t <- t[grepl('^Variogram cost', t)] %>%
    map_dbl(~as.numeric(gsub('Variogram cost: ', '', .x)))
  cv_weight_t <- t[grepl('^Interpolation weight', t)] %>%
    map_dbl(~as.numeric(gsub('Interpolation weight: ', '', .x)))
  cv_rmse_t <- t[grepl('^Interpolation rmse', t)] %>%
    map_dbl(~as.numeric(gsub('Interpolation rmse: ', '', .x)))
  cv_cost_t <- t[grepl('^Interpolation cost', t)] %>%
    map_dbl(~as.numeric(gsub('Interpolation cost: ', '', .x)))
  cost_t <- t[grepl('^Cost', t)] %>% map_dbl(~as.numeric(gsub('Cost: ', '', .x)))
  opt_t <- tibble(segment=segs_t, v_mod=vgrm_model_t, vgrm_wt=vgrm_weight_t,
                  vgrm_rmse=vgrm_rmse_t, vgrm_cost=vgrm_cost_t, cv_wt=cv_weight_t,
                  cv_rmse=cv_rmse_t, cv_cost=cv_cost_t, cost=cost_t)
  opt_t %>% group_by(segment, v_mod) %>% mutate(itr=row_number(), .before=segment) %>%
    ungroup()
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Krige !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Krige <- function(shp_hf=NULL, fitted_vgrm=NULL, shp_interp_grid=NULL, n_max=10,
                  seg_name=NULL) {
  if (is.null(shp_hf)) {stop('\nMissing heat flow data!')}
  if (is.null(fitted_vgrm)) {stop('\nMissing variogram model!')}
  if (is.null(shp_interp_grid)) {stop('\nMissing kriging locations (grid)!')}
  if (is.null(n_max)) {stop('\nNumer of max local point-pairs!')}
  suppressWarnings({
    k <- try(krige(hf~1, shp_hf, newdata=shp_interp_grid, model=fitted_vgrm,
                   nmax=n_max, debug.level=0))
    if (any(class(k) == 'try-error')) {stop('\nVariogram fitting error!')}
    k %>% as_tibble() %>% st_as_sf() %>% rename(est_krige=var1.pred, var_krige=var1.var) %>%
      mutate(sigma_krige=sqrt(var_krige), .before=geometry)
  })
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# decode opt !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
decode_opt <- function(shp_hf=NULL, model_vgrm=NULL, opt=NULL) {
  suppressWarnings({
    if (is.null(model_vgrm)) {stop('\nMissing variogram model!')}
    if (is.null(shp_hf)) {stop('\nMissing heat flow data!')}
    if (is.null(opt)) {stop('\nMissing nloptr object!')}
    experimental_vgrm <-
      try(experimental_vgrm(shp_hf, opt$solution[1], opt$solution[2], opt$solution[3]))
    if (nrow(experimental_vgrm) < 2) {
      stop('\nExperimental variogram has less than two lags!')
    }
    if (any(class(experimental_vgrm) == 'try-error')) {
      stop('\nExperimental variogram error!')
    }
    fitted_vgrm <- try(fit.variogram(experimental_vgrm, vgm(model=model_vgrm), fit.method=6))
    if (any(class(fitted_vgrm) == 'try-error')) {
      print(fitted_vgrm)
      stop('\nVariogram fitting error!')
    }
  })
  return(list('experimental_vgrm'=experimental_vgrm, 'fitted_vgrm'=fitted_vgrm))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# interp diff !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
interp_diff <- function(shp_interp_krige, shp_interp_sim) {
  if (is.null(shp_interp_krige)) {stop('\nMissing krige data!')}
  if (is.null(shp_interp_sim)) {stop('\nMissing similarity data!')}
  suppressWarnings(shp_interp_sim %>% st_intersection(shp_interp_krige) %>%
                   mutate(est_diff=est_sim - est_krige, sigma_diff=sigma_sim - sigma_krige,
                          .before=geometry))
}

#######################################################
## .3.            Plotting Functions             !!! ##
#######################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot tglobe !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_tglobe_base <- function() {
  fig_path <- 'figs/tglobe.png'
  if (!dir.exists('figs')) {
    dir.create('figs', recursive=T, showWarnings=F)
  }
  if (file.exists(fig_path)) {
    cat('\n', fig_path, 'already exists ...')
    return(invisible())
  }
  cat('\nplotting: tglobe base')
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
      geom_sf(data=shp_tglobe, aes(color=obs), size=0.1, shape=20) +
      scale_color_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)),
                            limits=c(0, 250), breaks=c(0, 125, 250), na.value='transparent',
                            guide=guide_colorbar(title.vjust=1, show.limits=T,
                                                 frame.colour='black',
                                                 ticks.colour='black')) +
      ggtitle('Thermoglobe Observations') +
      coord_sf(expand=F, lims_method='geometry_bbox') +
      theme_map(font_size=14)
    p3 <- p1 / p2 &
      theme(plot.margin=margin(), legend.position='bottom',
            legend.justification='center', legend.direction='horizontal',
            axis.text=element_blank(), legend.margin=margin(),
            legend.box.margin=margin(5, 5, 5, 5), legend.key.height=unit(0.5, 'cm'),
            legend.key.width=unit(0.6, 'cm'),
            legend.title=element_text(vjust=0, color='black', size=14),
            panel.grid=element_line(linewidth=0.05, color='grey20'),
            plot.title=element_text(vjust=0, hjust=0.5, margin=margin(10, 10, 10, 10)))
    ggsave(file=fig_path, plot=p3, width=6.5, height=6.5, bg='white')
  })})
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot transect tglobe !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_transect_tglobe <- function(short_name, base_size=20) {
  fig_path <- paste0('figs/transect_tglobe/', short_name,'-tglobe.png')
  if (!dir.exists('figs/transect_tglobe')) {
    dir.create('figs/transect_tglobe', recursive=T, showWarnings=F)
  }
  if (file.exists(fig_path)) {
    cat('\n', fig_path, 'already exists ...')
    return(invisible())
  }
  cat('\nCompiling map data for:', short_name)
  x <- compile_transect_data(short_name)
  cat('\nplotting: ', short_name)
  suppressWarnings({suppressMessages({
    p1 <-
      ggplot(x) +
      geom_sf(data=x$bathy[[1]], aes(color=elev), size=0.5, shape=15) +
      scale_color_etopo(name='Elevation (m)',
                        labels=label_number(scale_cut=cut_short_scale()),
                        guide=guide_colorbar(title.vjust=1, show.limits=T,
                                             frame.colour='black', ticks.colour='black')) +
      new_scale_color() +
      geom_sf(aes(geometry=large_buffer), fill=NA, linewidth=0.5) +
      geom_sf(aes(geometry=small_buffer), fill=NA, linewidth=0.5) +
      geom_sf(aes(geometry=ridge), color='white') +
      geom_sf(aes(geometry=transform), color='white') +
      geom_sf(aes(geometry=trench), color='white', linewidth=1.5) +
      geom_sf(aes(geometry=transect), color='black', linewidth=1.5) +
      geom_sf(aes(geometry=volcano), color='black', fill='white', shape=24) +
      geom_sf(data=x$tglobe_large_buff[[1]], aes(color=obs), shape=20) +
      scale_color_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                            breaks=c(0, 125, 250), na.value='transparent', guide='none') +
      xlab('Longitude') + ylab('Latitude') +
      ggtitle(paste0('Submap ', short_name, ': ', x$trench_name)) +
      coord_sf(expand=F, lims_method='geometry_bbox') +
      theme_bw(base_size=base_size) +
      theme(plot.margin=margin(5, 5, 5, 5), axis.text=element_text(hjust=1),
            panel.grid=element_line(linewidth=0.05, color='grey20'),
            plot.title=element_text(vjust=0, hjust=0.5, margin=margin(0, 0, 10, 10)),
            legend.justification='center', legend.position='bottom',
            legend.direction='horizontal', legend.key.height=unit(0.5, 'cm'),
            legend.key.width=unit(0.9, 'cm'), legend.box.margin=margin(5, 5, 5, 5),
            legend.margin=margin(), legend.title=element_text(vjust=0, size=base_size)) +
      annotation_scale(location='bl', width_hint=0.33, text_cex=1, style='ticks',
                       line_width=2.5, text_face='bold') +
      annotation_north_arrow(location='bl', which_north='true', pad_x=unit(0.0, 'cm'),
                             pad_y=unit(0.5, 'cm'), style=north_arrow_fancy_orienteering)
    p2 <-
      ggplot(x$tglobe_projected[[1]]) +
      geom_smooth(aes(projected_distances, obs), method='loess', formula=y~x, color='black') +
      geom_point(aes(projected_distances, obs), shape=20, color='grey20') +
      geom_rug(aes(projected_distances, obs, color=obs), length=unit(0.06, "npc")) +
      labs(x='Normalized Distance', y=bquote('Q'~(mWm^-2))) +
      scale_color_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                            breaks=c(0, 125, 250), na.value='transparent',
                            guide=guide_colorbar(title.vjust=1, show.limits=T,
                                                 frame.colour='black',
                                                 ticks.colour='black')) +
      scale_y_continuous(limits=c(0, 250), breaks=seq(0, 250, 50)) +
      scale_x_continuous(limits=c(0, 1), breaks=c(0, 0.5, 1)) +
      ggtitle('Projected Observations') +
      theme_bw(base_size=base_size) +
      theme(plot.margin=margin(5, 5, 5, 5),
            panel.grid=element_line(linewidth=0.05, color='grey20'),
            plot.title=element_text(vjust=0, hjust=0.5, margin=margin(0, 0, 10, 10)),
            legend.justification='center', legend.position='bottom',
            legend.direction='horizontal', legend.key.height=unit(0.5, 'cm'),
            legend.key.width=unit(0.9, 'cm'), legend.box.margin=margin(5, 5, 5, 5),
            legend.margin=margin(), legend.title=element_text(vjust=0, size=base_size))
    p3 <- p1 + p2 + plot_annotation(tag_levels = 'a', tag_suffix = ')') &
      theme(plot.tag = element_text(size=base_size*1.5))
    ggsave(file=fig_path, plot=p3, width=13, height=6.5, dpi=300, bg='white')
  })})
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot transect sim !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_transect_sim <- function(short_name) {
  fig_path <- paste0('figs/transect_sim/', short_name,'-sim.png')
  if (!dir.exists('figs/transect_sim')) {
    dir.create('figs/transect_sim', recursive=T, showWarnings=F)
  }
  if (file.exists(fig_path)) {
    cat('\n', fig_path, 'already exists ...')
    return(invisible())
  }
  cat('\nCompiling map data for:', short_name)
  x <- compile_transect_data(short_name)
  cat('\nplotting: ', short_name)
  suppressWarnings({suppressMessages({
    p <-
      ggplot(x) +
      geom_sf(data=x$bathy[[1]], aes(color=elev), size=0.5, shape=15) +
      scale_color_etopo(name='Elevation (m)',
                        guide=guide_colorbar(title.vjust=1, show.limits=T),
                        labels=label_number(scale_cut=cut_short_scale())) +
      new_scale_color() +
      geom_sf(aes(geometry=large_buffer), fill=NA, linewidth=0.5) +
      geom_sf(aes(geometry=ridge), color='white') +
      geom_sf(aes(geometry=transform), color='white') +
      geom_sf(aes(geometry=trench), color='white', linewidth=1.5) +
      geom_sf(aes(geometry=transect), color='black', linewidth=1.5) +
      geom_sf(aes(geometry=volcano), color='black', fill='white', shape=24) +
      geom_sf(data=x$sim[[1]], aes(color=est_sim), shape=20) +
      scale_color_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                            breaks=c(0, 125, 250), na.value='transparent', guide='none') +
      xlab('Longitude') + ylab('Latitude') +
      ggtitle(paste0('Submap ', short_name, ': ', x$trench_name)) +
      coord_sf(expand=F, lims_method='geometry_bbox') +
      theme_bw(base_size=14) +
      theme(plot.margin=margin(5, 5, 5, 5), legend.position='bottom',
            legend.justification='center', legend.direction='horizontal',
            axis.text=element_text(hjust=1), legend.margin=margin(),
            legend.box.margin=margin(5, 5, 5, 5), legend.key.height=unit(0.5, 'cm'),
            legend.key.width=unit(0.75, 'cm'), legend.text=element_text(size=base_size*0.694),
            legend.title=element_text(vjust=0, color='black', size=14),
            panel.grid=element_line(linewidth=0.05, color='grey20'),
            plot.title=element_text(vjust=0, hjust=0.5, margin=margin(0, 0, 10, 10))) +
      annotation_scale(location='bl', width_hint=0.33, text_cex=1, style='ticks',
                       line_width=2.5, text_face='bold') +
      annotation_north_arrow(location='bl', which_north='true', pad_x=unit(0.0, 'cm'),
                             pad_y=unit(0.5, 'cm'), style=north_arrow_fancy_orienteering)
    ggsave(file=fig_path, plot=p, width=6.5, height=6.5, dpi=300, bg='white')
  })})
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot vgrm !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_vgrm <- function(experimental_vgrm, fitted_vgrm=NULL, cost=NULL, v_mod=NULL,
                      line_col='black', ylim=NULL, xlim=NULL) {
  if (is.null(experimental_vgrm)) {stop('\nMissing experimental variogram!')}
  sill <- round(sqrt(fitted_vgrm[[2]]))
  range <- round(fitted_vgrm[[3]] / 1e3)
  p <- 
    ggplot() +
    geom_segment(aes(x=range, xend=range, y=Inf, yend=sill, color='range'),
                 arrow=arrow(angle=10, length=unit(0.1, 'inches'), type='closed'),
                 linewidth=0.8) +
    geom_segment(aes(x=-Inf, xend=range, y=sill, yend=sill, color='sill'),
                 arrow=arrow(angle=10, length=unit(0.1, 'inches'), type='closed'),
                 linewidth=0.8) +
    labs(x='lag distance (km)', y=bquote(sqrt(hat(gamma)[(h)])~(mWm^-2))) +
    scale_color_discrete_qualitative(name=NULL, palette='set 2') +
    new_scale_color() +
    coord_cartesian(ylim=ylim, xlim=xlim) +
    theme_bw(base_size=14) +
    theme(plot.margin=margin(5, 5, 5, 5), panel.background=element_rect(fill='grey90'),
          panel.border=element_rect(linewidth=1.2))
  plt <- tryCatch(
    {
      p +
      geom_point(data=experimental_vgrm, aes(x=dist / 1e3, y=sqrt(gamma)), shape=20) +
      geom_line(data=variogramLine(fitted_vgrm, maxdist=max(experimental_vgrm$dist)),
                aes(x=dist / 1e3, y=sqrt(gamma), color='variogram model')) +
      annotate('text', label=paste0('(', v_mod, ') cost: ', round(cost, 4)),
               x=xlim[2] / 2, y=0, vjust=-0.2, hjust=0.5) +
      scale_color_manual(name=NULL, values=c('experimental variogram'='black',
                                             'variogram model'=line_col))
    },
    error=function(cond) {
      exp_plt <- p +
      geom_point(data=experimental_vgrm, aes(x=dist / 1e3, y=sqrt(gamma)), shape=20)
      return(exp_plt)
    }
  )
  return(plt)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot split segment !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_split_segment <- function(split_seg, ext=c(0.1, 0.1, 0.1, 0.1), running_avg=3) {
  seg_name <- split_seg$interp[[1]]$segment[1]
  seg_num <- as.numeric(unique(bind_rows(split_seg$pnts)$split_fID))
  bx <- bbox_extend(st_bbox(st_buffer(st_combine(split_seg$seg), dist=5e5)), ext=ext)
  seg <- bind_rows(split_seg$seg)
  ridge <- shp_ridge_crop[[seg_name]]
  trench <- shp_trench_crop[[seg_name]]
  transform <- shp_transform_crop[[seg_name]]
  relief <- shp_relief_crop[[seg_name]]
  volc <- split_seg$volc
  wdth <- range(st_bbox(bind_rows(split_seg$buf))[c('xmin', 'xmax')])/1e3
  rng <- bind_rows(split_seg$interp) %>% st_set_geometry(NULL) %>%
    select(est_sim, est_krige, split_fID) %>% summarise(min_krige=min(est_krige),
                                                        max_krige=max(est_krige),
                                                        min_sim=min(est_sim),
                                                        max_sim=max(est_sim))
  pnt_size <- 1
  p0 <-
    ggplot() +
    geom_sf(data=bx, color=NA, fill=NA) +
    geom_sf(data=relief, aes(color=elevation), shape=15, size=0.25) +
    scale_color_etopo(guide='none') +
    new_scale_color() +
    geom_sf(data=ridge, linewidth=0.5) +
    geom_sf(data=trench, linewidth=0.5) +
    geom_sf(data=transform, linewidth=0.5) +
    geom_sf(data=seg, linewidth=1, color='white') +
    geom_sf(data=bind_rows(split_seg$buf),
            aes(fill=factor(split_fID, levels=seg_num[order(seg_num)])),
            linewidth=0.5, color='black', show.legend=F) +
    geom_sf(data=bind_rows(split_seg$pnts),
            aes(fill=factor(split_fID, levels=seg_num[order(seg_num)]),
                group=factor(split_fID, levels=seg_num[order(seg_num)])),
            size=pnt_size, shape=22, stroke=0.5, show.legend=F) +
    geom_sf(data=bind_rows(volc), color='white', shape=18, size=pnt_size) +
    annotate('label', label='a', x=-Inf, y=Inf, size=7, hjust=0, vjust=1, fill='white',
             label.padding=unit(0.02, 'in'), label.r=unit(0, 'in')) +
    coord_sf(xlim=c(st_bbox(bx)$xmin, st_bbox(bx)$xmax),
             ylim=c(st_bbox(bx)$ymin, st_bbox(bx)$ymax), label_axes='--EN') +
    scale_color_discrete_qualitative('set 2') +
    scale_fill_discrete_qualitative('set 2') +
    guides(fill=guide_legend(nrow=1, override.aes=list(alpha=1, size=8),
                             title.position='top', title.vjust=1, label.position='bottom')) +
    theme_map(font_size=9) +
    theme(plot.margin=margin(5, 5, 5, 5),
          axis.text.x=element_text(color='grey40', angle=30, hjust=0, vjust=0),
          axis.text.y=element_text(color='grey40', angle=30, hjust=0),
          panel.grid=element_line(size=0.01, color='grey60'))
    if (seg_name %in% c('Alaska Aleutians', 'Kamchatka Marianas', 'Tonga New Zealand',
                        'Vanuatu')) {
      p0 <- p0 + scale_x_continuous(breaks=c(80, 90, 100, 110, 120, 130, 140, 150, 160, 170,
                                             180, -170, -160, -150, -140, -130))
    }
  p1 <-
    bind_rows(split_seg$interp) %>%
    st_set_geometry(NULL) %>%
    select(est_sim, est_krige, split_fID) %>%
    rename(Similarity=est_sim, Kriging=est_krige) %>%
    pivot_longer(-split_fID) %>%
    group_by(name) %>%
    filter(split_fID %in% seg_num) %>%
    ggplot() +
    geom_boxplot(aes(y=value, x=factor(split_fID, levels=seg_num[order(seg_num)]),
                     fill=factor(split_fID, levels=seg_num[order(seg_num)]), color=name,
                     outlier.color=name), outlier.size=pnt_size * 0.2, size=0.5) +
    annotate('label', label='b', x=-Inf, y=Inf, size=7, hjust=0, vjust=1, fill='white',
             label.padding=unit(0.02, 'in'), label.r=unit(0, 'in')) +
    guides(color=guide_legend(nrow=1, override.aes=list(alpha=1, fill=NA),
                              title.position='top', title.vjust=1, label.position='left'),
           fill='none') +
    coord_cartesian(ylim=c(min(rng$min_krige, rng$min_sim),
                           max(rng$max_krige, rng$max_sim))) +
    scale_y_continuous(position='right') +
    scale_fill_discrete_qualitative(palette='set 2') +
    scale_color_manual(values=c('black', 'black')) +
    labs(y=bquote('heat flow'~(mWm^-2)), x=NULL, color='interpolation method', fill=NULL) +
    theme_bw(base_size=14) +
    theme(legend.title=element_text(margin=margin(0, 0, -5, 0)),
          plot.margin=margin(5, 5, 5, 5), panel.background=element_rect(fill='grey90'),
          panel.border=element_rect(linewidth=1.2))
  labl <- tibble(name=c('Kriging', 'Similarity'), lbl=c('c', 'd'))
  p2 <-
    bind_rows(split_seg$interp) %>%
    st_set_geometry(NULL) %>%
    filter(distance_from_seg <= 500000) %>%
    select(est_sim, est_krige, split_fID, distance_from_seg) %>%
    arrange(distance_from_seg) %>%
    mutate(Similarity=rollmean(est_sim, running_avg, fill=NA),
           Kriging=rollmean(est_krige, running_avg, fill=NA)) %>%
    select(-c(est_sim, est_krige)) %>%
    pivot_longer(-c(split_fID, distance_from_seg)) %>%
    group_by(name) %>%
    ggplot() +
    geom_point(data=bind_rows(volc),
               aes(distance_from_seg / 1e3, min(rng$min_krige, rng$min_sim)), color='black',
               size=pnt_size, shape=18, key_glyph='rect') +
    geom_point(data=bind_rows(split_seg$pnts),
               aes(distance_from_seg / 1e3, hf,
                   fill=factor(split_fID, levels=seg_num[order(seg_num)]),
                   group=factor(split_fID, levels=seg_num[order(seg_num)])),
               size=pnt_size, shape=22, key_glyph='rect') +
    geom_smooth(aes(distance_from_seg / 1e3, value,
                    color=factor(split_fID, levels=seg_num[order(seg_num)]),
                    group=factor(split_fID, levels=seg_num[order(seg_num)])),
                fill=NA, method='loess', formula = y ~ x, alpha=0.1, size=1, se=T,
                key_glyph='rect', show.legend=F) +
    geom_label(data=labl, aes(label=lbl), x=-Inf, y=Inf, size=7, hjust=0, vjust=1,
               fill='white', label.padding=unit(0.02, 'in'), label.r=unit(0, 'in')) +
    labs(x='distance from trench (km)', y=bquote('heat flow'~(mWm^-2)), fill='sector') +
    guides(fill=guide_legend(nrow=1, override.aes=list(alpha=1, size=5), title.position='top',
                             title.vjust=1, label.position='bottom')) +
    coord_cartesian(ylim=c(min(rng$min_krige, rng$min_sim),
                           max(rng$max_krige, rng$max_sim))) +
    scale_color_discrete_qualitative('set 2') +
    scale_fill_discrete_qualitative('set 2') +
    facet_wrap(~name, ncol=2) +
    theme_bw(base_size=14) +
    theme(legend.title=element_text(margin=margin(0, 0, -5, 0)),
          strip.background=element_rect(color=NA, fill=NA),
          strip.text=element_text(color='black', size=14, margin=margin(0, 0, 2, 0)),
          strip.placement='outside', legend.box.margin=margin(0, 0, 0, 0),
          legend.key=element_rect(color='black'), legend.spacing.x=unit(0, 'mm'),
          plot.margin=margin(5, 5, 5, 5), panel.background=element_rect(fill='grey90'),
          panel.border=element_rect(linewidth=1.2))
  p <-
    ((p0 | p1) / p2) + plot_annotation(title='Comparing heat flow interpolations by sector') +
    plot_layout(widths=1, guides='collect') &
    theme(plot.title=element_text(size=18), plot.margin=margin(5, 5, 5, 5),
          legend.box.margin=margin(0, 0, 0, 0), legend.margin=margin(),
          legend.position='bottom', legend.direction='horizontal',
          legend.justification='left')
  p
}
