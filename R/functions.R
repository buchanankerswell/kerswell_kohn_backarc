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
  if (l == 0) {
    if(!keep_df) {NA} else {NULL}
  } else {
    if (!keep_df) {st_sfc(st_union(ft_cropped), crs=st_crs(bbox))} else {ft_cropped}
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
    cat('\nAn error occurred in get_world_bathy:', conditionMessage(e), '\n')
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
    cat('\nAn error occurred in get_seg_bathy:', conditionMessage(e), '\n')
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
    stop('\nMissing map data! Use "load(path/to/map-data.RData)"')
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
  cat('\nCompiling transect data for submap transect:', trans_id)
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
      mutate(sim_large_buff=list(crop_feature(shp_sim, large_buffer, T, T))) %>%
      mutate(sim_small_buff=list(crop_feature(shp_sim, small_buffer, T, T))) %>%
      mutate(sim_projected=list(project_obs_to_transect(transect, sim_small_buff))) %>%
      mutate(bathy=list(get_seg_bathy(bbox))) %>%
      ungroup()
  })})
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
  if (!any(names(shp_obs) %in% c('obs', 'est_sim'))) {
    stop('\nUnrecognized sf object passed to project_obs_to_transect() !')
  }
  if (any(names(shp_obs) == 'tglobe')) {
    obs <- shp_obs$obs
  } else if (any(names(shp_obs) == 'similarity')) {
    obs <- shp_obs$est_sim
  } else if (any(names(shp_obs) == 'krige')) {
    obs <- shp_obs$est_krg
  }
  projected_distances <- st_line_project(transect, st_geometry(shp_obs), normalized=T)
  st_as_sf(st_line_interpolate(transect, projected_distances)) %>%
    mutate(projected_distances=projected_distances, obs=obs)
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
cost_function <- function(shp_hf=NULL, cutoff=3, n_lags=50, n_max=10, v_mod='Sph',
                          n_fold=NULL, interp_weight=0.5, vgrm_weight=0.5, trans_id=NULL) {
  if (is.null(shp_hf)) {stop('\nMissing heat flow data model!')}
  if (is.null(n_fold)) {n_fold <- nrow(shp_hf)}
  if (!is.null(n_fold)) {if (n_fold > 0 & n_fold <= 1) {n_fold <- nrow(shp_hf) * n_fold}}
  suppressWarnings({
    tryCatch({
      ev <- experimental_vgrm(shp_hf, cutoff, n_lags)
    }, error=function(e) {
      cat('\nAn error occurred in experimental_vgrm:', conditionMessage(e), '\n')
      return(runif(1, 1, 1.5))
    })
    if (nrow(ev) < 2) {return(runif(1, 1, 1.5))}
    tryCatch({
      fv <- fit.variogram(ev, vgm(model=v_mod), fit.method=6)
    }, error=function(e) {
      cat('\nAn error occurred in fit.variogram:', conditionMessage(e), '\n')
      return(runif(1, 1, 1.5))
    })
    tryCatch({
      k_cv <- krige.cv(obs~1, shp_hf, model=fv, nmax=n_max, nfold=n_fold)
    }, error=function(e) {
      cat('\nAn error occurred in krige.cv:', conditionMessage(e), '\n')
      return(runif(1, 1, 1.5))
    })
    if (sum(is.na(k_cv)) != 0) {
      if (sum(is.na(k_cv$residual)) >= nrow(k_cv) / 2) {return(runif(1, 1, 1.5))}
    }
  })
  k_cv <- k_cv %>% filter(!is.na(residual))
  suppressWarnings({
    vgrm_rmse <- sqrt(sqrt(attr(fv, 'SSErr') / nrow(ev)))
    vgrm_sd <- sqrt(sd(ev$gamma, na.rm=T))
    vgrm_cost <- vgrm_weight * vgrm_rmse / vgrm_sd
    interp_rmse <- sqrt(sum(k_cv$residual^2, na.rm=T) / nrow(k_cv))
    interp_sd <- sd(k_cv$var1.pred, na.rm=T)
    interp_cost <- interp_weight * interp_rmse / interp_sd
  })
  log_dir <- 'assets/opt_data/nlopt_traces'
  log_file <- paste0(log_dir, '/nlopt-out-', trans_id, '-', v_mod)
  if (!dir.exists(log_dir)) {dir.create(log_dir, recursive=T, showWarnings=F)}
  if (!file.exists(log_file)) {file.create(log_file, showWarnings=F)}
  sink(log_file, append=T)
  cat('\nTransect:', trans_id)
  cat('\nCutoff proportion:', cutoff)
  cat('\nNumber of lags:', n_lags)
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
  cat('\nVariogram model:', v_mod)
  cat('\n')
  print(fv)
  cat(rep('+', 30), sep='')
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
      stop('\nAn error occurred in experimental_vgrm:', conditionMessage(e), '\n')
    })
    if (nrow(ev) < 2) {stop('\nExperimental variogram has less than two lags!')}
    if (any(class(ev) == 'try-error')) {stop('\nExperimental variogram error!')}
    tryCatch({
      fv <- fit.variogram(ev, vgm(model=v_mod), fit.method=6)
    }, error=function(e) {
      stop('\nAn error occurred in fit.variogram:', conditionMessage(e), '\n')
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
  opt_dir <- 'assets/opt_data'
  opt_id <- paste0('opt-', trans_id, '-', v_mod) 
  opt_path <- paste0(opt_dir, '/', opt_id, '.RData')
  if (file.exists(opt_path)) {
    cat('\nOptimized', v_mod, 'kriging model found for submap transect:', trans_id)
    return(invisible())
  }
  if (!dir.exists(opt_dir)) {dir.create(opt_dir, recursive=T, showWarnings=F)}
  x0 <- c(3, 50, 10) # Initial values (cutoff, n_lags, n_max)
  lb <- c(1, 30, 2) # Lower bound (cutoff, n_lags, n_max)
  ub <- c(12, 100, 50) # Upper bound (cutoff, n_lags, n_max)
  opts <- list(print_level=0, maxeval=max_eval, algorithm=alg, xtol_rel=1e-5, ftol_rel=1e-5)
  x <- compile_transect_data(trans_id)
  obs <- x$tglobe_large_buff[[1]]
  nlopt_fun <- function(x) {
    cost_function(obs, x[1], x[2], x[3], v_mod, n_fold, iwt, vwt, trans_id)
  }
  tryCatch({
    cat('\nOptimizing', v_mod, 'kriging model for submap transect:', trans_id)
    opt <- nloptr(x0, nlopt_fun, lb=lb, ub=ub, opts=opts)
    opt_decoded <- decode_opt(obs, v_mod, opt)
    assign(str_replace_all(opt_id, '-', '_'), opt_decoded)
    save(list=str_replace_all(opt_id, '-', '_'), file=opt_path)
  }, error=function(e) {
    cat('\nAn error occurred in nlopt_krige:', conditionMessage(e), '\n')
  })
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# nlopt transects !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nlopt_transects <- function(trans_ids=NULL, v_mods=c('Sph', 'Exp', 'Lin', 'Bes'),
                            alg='NLOPT_LN_COBYLA', max_eval=500, n_fold=10,
                            iwt=0.5, vwt=0.5, parallel=T) {
  if (is.null(trans_ids)) {stop('\nMissing submap transect ids!')}
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
# read nloptr trace !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
read_nloptr_trace <- function(fpath) {
  read_file <- function(path, ...) {
    con <- file(path)
    on.exit(close(con))
    suppressWarnings(readLines(con, ...))
  }
  t <- read_file(fpath)
  trans_id_t <- t[grepl('^Transect: ', t)] %>% map_chr(~gsub('Transect: ', '', .x))
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
  opt_t <- tibble(short_name=trans_id_t, v_mod=vgrm_model_t, vgrm_wt=vgrm_weight_t,
                  vgrm_rmse=vgrm_rmse_t, vgrm_cost=vgrm_cost_t, cv_wt=cv_weight_t,
                  cv_rmse=cv_rmse_t, cv_cost=cv_cost_t, cost=cost_t)
  opt_t %>% group_by(short_name, v_mod) %>% mutate(itr=row_number(), .before=short_name) %>%
    ungroup()
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Krige !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Krige <- function(shp_hf=NULL, fv=NULL, shp_interp_grid=NULL, n_max=10,
                  seg_name=NULL) {
  if (is.null(shp_hf)) {stop('\nMissing heat flow data!')}
  if (is.null(fv)) {stop('\nMissing variogram model!')}
  if (is.null(shp_interp_grid)) {stop('\nMissing kriging locations (grid)!')}
  if (is.null(n_max)) {stop('\nNumer of max local point-pairs!')}
  suppressWarnings({
    k <- try(krige(hf~1, shp_hf, newdata=shp_interp_grid, model=fv,
                   nmax=n_max, debug.level=0))
    if (any(class(k) == 'try-error')) {stop('\nVariogram fitting error!')}
    k %>% as_tibble() %>% st_as_sf() %>% rename(est_krige=var1.pred, var_krige=var1.var) %>%
      mutate(sigma_krige=sqrt(var_krige), .before=geometry)
  })
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
  if (!dir.exists('figs')) {dir.create('figs', recursive=T, showWarnings=F)}
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
# plot transect !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_transect <- function(short_name, base_size=20) {
  fig_dir <- 'figs/transect/'
  fig_path <- paste0(fig_dir, short_name, '-tglobe-sim.png')
  if (!dir.exists(fig_dir)) {dir.create(fig_dir, recursive=T, showWarnings=F)}
  if (file.exists(fig_path)) {cat('\n', fig_path, 'already exists ...'); return(invisible())}
  x <- compile_transect_data(short_name)
  cat('\nPlotting:', fig_path)
  suppressWarnings({suppressMessages({
    map_theme <- list(
      theme(plot.margin=margin(2, 2, 2, 2),
            panel.grid=element_line(linewidth=0.05, color='grey20'),
            plot.title=element_text(vjust=0, hjust=0.5, margin=margin(0, 0, 10, 0)),
            legend.justification='center', legend.position='bottom',
            legend.direction='horizontal', legend.key.height=unit(0.5, 'cm'),
            legend.key.width=unit(1, 'cm'), legend.box.margin=margin(2, 2, 2, 2),
            legend.margin=margin(), legend.title=element_text(vjust=0, size=base_size)),
      annotation_scale(location='bl', width_hint=0.33, text_cex=1, style='ticks',
                       line_width=2.5, text_face='bold'),
      annotation_north_arrow(location='bl', which_north='true', pad_x=unit(0.0, 'cm'),
                             pad_y=unit(0.5, 'cm'), style=north_arrow_fancy_orienteering)
    )
    profile_theme <-
      theme(plot.margin=margin(2, 2, 2, 2),
            panel.grid=element_line(linewidth=0.05, color='grey20'),
            plot.title=element_text(vjust=0, hjust=0.5, margin=margin(0, 0, 10, 0)),
            legend.justification='center', legend.position='bottom',
            legend.direction='horizontal', legend.key.height=unit(0.5, 'cm'),
            legend.key.width=unit(1, 'cm'), legend.box.margin=margin(2, 2, 2, 2),
            legend.margin=margin(), legend.title=element_text(vjust=0, size=base_size))

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
      theme_bw(base_size=base_size) + map_theme
    if (is.null(x$tglobe_projected[[1]])) {
      p2 <- ggplot(data=data.frame(), aes())
    } else {
      p2 <-
        ggplot(x$tglobe_projected[[1]]) +
        geom_smooth(aes(projected_distances, obs), method='loess', formula=y~x,
                    color='black') +
        geom_point(aes(projected_distances, obs), shape=20, color='grey20') +
        geom_rug(aes(projected_distances, obs, color=obs), length=unit(0.06, "npc"))
    }
    p2 <- p2 +
      labs(x='Normalized Distance', y=bquote('Q'~(mWm^-2))) +
      scale_color_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                            breaks=c(0, 125, 250), na.value='transparent',
                            guide=guide_colorbar(title.vjust=1, show.limits=T,
                                                 frame.colour='black',
                                                 ticks.colour='black')) +
      scale_y_continuous(limits=c(0, 250), breaks=seq(0, 250, 50)) +
      scale_x_continuous(limits=c(0, 1), breaks=c(0, 0.5, 1)) +
      ggtitle('Thermoglobe Observations') +
      theme_bw(base_size=base_size) + profile_theme
    p3 <-
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
      geom_sf(data=x$sim_large_buff[[1]], aes(color=est_sim), shape=20) +
      scale_color_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                            breaks=c(0, 125, 250), na.value='transparent', guide='none') +
      xlab('Longitude') + ylab('Latitude') +
      ggtitle(paste0('Submap ', short_name, ': ', x$trench_name)) +
      coord_sf(expand=F, lims_method='geometry_bbox') +
      theme_bw(base_size=base_size) + map_theme
    if (is.null(x$sim_projected[[1]])) {
      p4 <- ggplot(data=data.frame(), aes())
    } else {
      p4 <-
        ggplot(x$sim_projected[[1]]) +
        geom_smooth(aes(projected_distances, obs), method='loess', formula=y~x,
                    color='black') +
        geom_point(aes(projected_distances, obs), shape=20, color='grey20') +
        geom_rug(aes(projected_distances, obs, color=obs), length=unit(0.06, "npc"))
    }
    p4 <- p4 +
      labs(x='Normalized Distance', y=bquote('Q'~(mWm^-2))) +
      scale_color_viridis_c(option='magma', name=bquote('Q'~(mWm^-2)), limits=c(0, 250),
                            breaks=c(0, 125, 250), na.value='transparent',
                            guide=guide_colorbar(title.vjust=1, show.limits=T,
                                                 frame.colour='black',
                                                 ticks.colour='black')) +
      scale_y_continuous(limits=c(0, 250), breaks=seq(0, 250, 50)) +
      scale_x_continuous(limits=c(0, 1), breaks=c(0, 0.5, 1)) +
      ggtitle('Similarity Interpolation') +
      theme_bw(base_size=base_size) + profile_theme
    p5 <-
      (p1 + p2 & theme(legend.position='none', axis.title.x=element_blank())) /
      (p3 + p4) +
      plot_annotation(tag_levels = 'a', tag_suffix = ')') &
      plot_layout(widths=1, heights=1) &
      theme(plot.tag = element_text(size=base_size * 1.5))
    ggsave(file=fig_path, plot=p5, width=13, height=11, dpi=300, bg='white')
  })})
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot optimal variogram !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_optimal_variogram <- function(trans_id=NULL, base_size=14) {
  if (is.null(trans_id)) {stop('\nMissing submap transect ids!')}
  fig_dir <- 'figs/opt_variogram/'
  fig_path <- paste0(fig_dir, trans_id, '-opt-variogram.png')
  if (!dir.exists(fig_dir)) {dir.create(fig_dir, recursive=T, showWarnings=F)}
  if (file.exists(fig_path)) {cat('\n', fig_path, 'already exists ...'); return(invisible())}
  x <- compile_transect_data(trans_id)
  shp_hf <- x$tglobe_large_buff[[1]]
  cat('\nPlotting:', fig_path)
  opt_dir <- 'assets/opt_data/'
  trace_dir <- paste0(opt_dir, '/nlopt_traces')
  trace_paths <- list.files(trace_dir, pattern=trans_id, full.names=T)
  opt_trace <- map_df(trace_paths, read_nloptr_trace)
  abs_min_cost <- opt_trace %>% slice_min(cost)
  best_vmod <- abs_min_cost$v_mod
  load(paste0(opt_dir, 'opt-', trans_id, '-', best_vmod, '.RData'))
  opt_decoded <- get(paste0('opt_', trans_id, '_', best_vmod))
  ev <- opt_decoded$experimental_vgrm
  fv <- opt_decoded$fitted_vgrm
  fv_line <- variogramLine(fv, maxdist=max(ev$dist))
  p_theme <-
    theme(plot.margin=margin(2, 2, 2, 2),
          panel.grid=element_line(linewidth=0.05, color='grey20'),
          plot.title=element_text(vjust=0, hjust=0.5, margin=margin(0, 0, 10, 0)),
          legend.position='none')
  p1 <-
    ggplot(filter(opt_trace, v_mod == best_vmod)) +
    geom_path(aes(itr, vgrm_cost), color='firebrick', linewidth=1) +
    geom_point(data=abs_min_cost, aes(itr, vgrm_cost), color='firebrick', shape=20, size=5) +
    geom_label_repel(data=abs_min_cost, aes(itr, vgrm_cost, label=round(vgrm_cost, 3)),
                     label.padding=0.12) +
    geom_path(aes(itr, cv_cost), color='navy', linewidth=1) +
    geom_point(data=abs_min_cost, aes(itr, cv_cost), color='navy', shape=20, size=5) +
    geom_label_repel(data=abs_min_cost, aes(itr, cv_cost, label=paste0(round(cv_cost, 3))),
                     label.padding=0.12) +
    ggtitle(paste0('Submap: ', unique(opt_trace$short_name), ' Kriging Optimization')) +
    labs(x='Iteration', y='Cost', color=NULL) +
    theme_bw(base_size=base_size) + p_theme
  p2 <- 
    ggplot(ev) +
    geom_point(aes(x=dist / 1e3, y=sqrt(gamma)), shape=19) +
    geom_line(data=fv_line, aes(x=dist / 1e3, y=sqrt(gamma))) +
    ggtitle(paste0('Best model: ', v_mod, ' (', round(abs_min_cost$cost, 3),
                   ' total cost)')) +
    labs(x='Lag Distance (km)', y=bquote('Variance'~(mWm^-2))) +
    theme_bw(base_size=base_size) + p_theme
  p3 <- p1 / p2 +
    plot_annotation(tag_level='a', tag_suffix=')') +
    plot_layout(widths=1, heights=1) &
    theme(plot.tag = element_text(size=base_size * 1.5))
  ggsave(file=fig_path, plot=p3, width=6.5, height=6.5, dpi=300, bg='white')
}
