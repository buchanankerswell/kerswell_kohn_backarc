#######################################################
## .0.              Load Libraries               !!! ##
#######################################################
# Quiet loading
sshhh <- function(package_list) {
  suppressWarnings(suppressPackageStartupMessages({
    require(package_list, quietly = T, character.only = T)
  }))
}

# Package list
package_list <- c('tictoc', 'stringr', 'tidyr', 'readr', 'purrr', 'furrr', 'tibble', 'dplyr',
                  'magrittr', 'ggplot2', 'colorspace', 'metR', 'ggrepel', 'ggridges',
                  'ggnewscale', 'patchwork', 'cowplot', 'ggsflabel', 'marmap', 'gstat',
                  'rgeos', 'sf', 'stars', 'rnaturalearth', 'nloptr')

# Load packages quietly
sapply(package_list, sshhh)
rm(package_list, sshhh)

# Don't allow sf to use google's s2 library for spherical geometry
# This reverts to using the GEOS library instead which is what sf used before 1.0 release
sf_use_s2(F)

#######################################################
## .1.         General Helper Functions          !!! ##
#######################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get RData object !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
extract_RData_object <- function(file, object) {
  E <- new.env()
  load(file = file, envir = E)
  return(get(object, envir = E, inherits = F))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# interpolation rmse !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
interpolation_rmse <- function(seg_name, interpolation = NULL, shp_buffer = NULL,
                               shp_grid_crop = NULL, shp_hf_crop = NULL, type = 'sim') {
  if (is.null(interpolation)) {stop('need interpolation sf object!')}
  if (is.null(shp_buffer)) {stop('need buffer sf object!')}
  if (is.null(shp_grid_crop)) {stop('need grid_crop sf object!')}
  if (is.null(shp_hf_crop)) {stop('need hf_crop sf object!')}
  buf <- shp_buffer[[seg_name]]
  grd <- shp_grid_crop[[seg_name]]
  obs <- shp_hf_crop[[seg_name]]
  if (type == 'sim') {
    interpolation <- suppressWarnings(st_intersection(interpolation, buf))
    nearest_est_sim <- interpolation$est_sim[st_nearest_feature(obs, grd)]
    return(sqrt(mean((nearest_est_sim - obs$hf)^2)))
  } else if (type == 'krg') {
    interpolation <- interpolation
    nearest_est_krige <- interpolation$est_krige[st_nearest_feature(obs, grd)]
    return(sqrt(mean((nearest_est_krige - obs$hf)^2)))
  } else {
    stop('invalid type!')
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# bbox widen !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bbox_widen <- function(bbox, crs = NULL, borders = c(0, 0, 0, 0)) {
  if (is.null(bbox)) stop('\nMissing bounding box!')
  if (is.null(crs)) {crs <- st_crs(bbox)}
  b <- bbox
  xrange <- b$xmax - b$xmin
  yrange <- b$ymax - b$ymin
  b[1] <- b[1] - (borders[1] * xrange)
  b[3] <- b[3] + (borders[2] * xrange)
  b[4] <- b[4] + (borders[3] * yrange)
  b[2] <- b[2] - (borders[4] * yrange)
  box <- c(b$xmin, b$ymax, b$xmin, b$ymin, b$xmax, b$ymin, b$xmax, b$ymax, b$xmin, b$ymax)
  shp_box <- st_polygon(list(matrix(box, ncol = 2, byrow = T))) %>% st_sfc(crs = crs)
  return(shp_box)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read latlong !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
read_latlong <- function(file = NULL, crs = NULL) {
  if (is.null(file)) {stop('\nMissing filename!')}
  if (is.null(crs)) {stop('\nMissing coordinate reference system!')}
  pattern <- '(?<=gmts\\/)[a-z].*(?=_contours\\.gmt)'
  seg_name <- file %>% str_extract(pattern) %>% str_replace_all('_', ' ') %>% str_to_title()
  proj <- '+proj=longlat +lon_wrap=180 +ellps=WGS84 +datum=WGS84 +no_defs'
  st_read(file, crs = proj, quiet = T) %>% st_transform(crs) %>% as_tibble() %>%
    st_as_sf() %>% mutate(segment = seg_name, .before = geometry)
}

#######################################################
## .2.             Kriging Functions             !!! ##
#######################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# parse krige args !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
parse_krige_args <- function(args) {
  if (length(args) == 0) {
    cat('\nNo arguments passed to krige.R')
    cat('\nUsing defaults')
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
        alg <- 'NLOPT_GN_DIRECT_L'   # Global search
      } else if (args[2] == 2) {
        alg <- 'NLOPT_LN_SBPLX'      # Local without gradients
      } else if (args[2] == 3) {
        alg <- 'NLOPT_LN_NELDERMEAD' # Local without gradients
      } else if (args[2] == 4) {
        alg <- 'NLOPT_LN_BOBYQA'     # Local without gradients
      } else if (args[2] == 5) {
        alg <- 'NLOPT_LN_COBYLA'     # Local without gradients
      } else {
        alg <- 'NLOPT_LN_SBPLX'
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
  assign('max_eval', max_eval, envir = .GlobalEnv)
  assign('alg', alg, envir = .GlobalEnv)
  assign('iwt', iwt, envir = .GlobalEnv)
  assign('vwt', vwt, envir = .GlobalEnv)
  assign('n_fold', n_fold, envir = .GlobalEnv)
  assign('n_cores', n_cores, envir = .GlobalEnv)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Krige !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Krige <-
  function(
    shp_hf = NULL,
    fitted_vgrm = NULL,
    shp_interp_grid = NULL,
    n_max = 20,
    seg_name = NULL
  ) {
  # Check for missing arguments
  if (is.null(shp_hf)) stop('\nMissing heat flow data!')
  if (is.null(fitted_vgrm)) stop('\nMissing variogram model!')
  if (is.null(shp_interp_grid)) stop('\nMissing kriging locations (grid)!')
  if (is.null(n_max)) stop('\nNumer of max local point-pairs!')
  # Print grid and parameters info
  cat(
    '\nKriging:', seg_name,
    '\nObservations:', nrow(shp_hf),
    '\nGrid size:', length(shp_interp_grid),
    '\nMax nearest points:', n_max,
    '\nVariogram params:\n'
  )
  print(fitted_vgrm)
  # Kriging
  k <-
    try(
      suppressWarnings(
        krige(
          formula = hf~1,
          locations = shp_hf,
          newdata = shp_interp_grid,
          model = fitted_vgrm,
          nmax = n_max,
          debug.level = 0
        ) %>%
        as_tibble() %>%
        st_as_sf() %>%
        rename(
          est_krige = var1.pred,
          var_krige = var1.var
        ) %>%
        mutate(
          sigma_krige = sqrt(var_krige),
          .before = geometry
        )
      )
    )
  # Check for error during fitting
  if (any(class(k) == 'try-error')) {
    # Print grid and parameters info
    cat(
      '\nKriging:', seg_name,
      '\nObservations:', nrow(shp_hf),
      '\nGrid size:', length(shp_interp_grid),
      '\nMax nearest points:', n_max,
      '\nVariogram params:\n'
    )
    print(fitted_vgrm)
    stop('\nVariogram fitting error!')
  }
  return(k)
}

# Take difference
interp_diff <- function(shp_interp_krige, shp_interp_sim) {
  # Check for missing arguments
  if (is.null(shp_interp_krige)) stop('\nMissing krige data!')
  if (is.null(shp_interp_sim)) stop('\nMissing similarity data!')
  # Crop Similarity estimates to Kriging bounding box
  dif <-
    suppressWarnings(
      shp_interp_sim %>%
      st_intersection(shp_interp_krige) %>%
      mutate(
        est_diff = est_sim - est_krige,
        sigma_diff = sigma_sim - sigma_krige,
        .before = geometry
      )
    )
  return(dif)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# experimental vgrm !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
experimental_vgrm <- function(shp_hf = NULL, cutoff_prop = 3, n_lags = 20, lag_start = 1) {
  # Check for missing arguments
  if (is.null(shp_hf)) stop('\nMissing heat flow data!')
  # Get bounding box
  bbox <- st_bbox(shp_hf)
  # Calculate corner-to-corner distance of bbox
  bbox_diagonal_distance <-
    sqrt((bbox$xmax-bbox$xmin)^2 + (bbox$ymax-bbox$ymin)^2)
  # Lag cutoff is a proportion of the corner-to-corner bounding box distance
  # note: defaults to the bounding box diagonal divided by three
  lag_cutoff <- as.vector(bbox_diagonal_distance/cutoff_prop)
  bin_width <- lag_cutoff/n_lags
  # Calculate new cutoff for shifting variogram lag_start to the left
  shifted_cutoff <- bin_width*(n_lags+lag_start-1)
  # Return shifted experimental variogram
  variogram(hf~1, locations = shp_hf, cutoff = shifted_cutoff, width = bin_width) %>%
    slice(floor(lag_start):n())
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# cost function !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cost_function <- function(shp_hf, cutoff_prop = 3, n_lags = 30, lag_start = 3,
                          model_vgrm = 'Sph', n_max = 20, n_fold = NULL, interp_weight = 0.5,
                          vgrm_weight = 0.5, segment = NULL, verbose = T) {
  if (is.null(shp_hf)) {stop('\nMissing heat flow data model!')}
  if (is.null(n_fold)) {n_fold <- nrow(shp_hf)}
  if (n_fold > 0 & n_fold <= 1) {n_fold <- nrow(shp_hf) * n_fold}
  experimental_vgrm <-
    try(experimental_vgrm(shp_hf = shp_hf, cutoff_prop = cutoff_prop, n_lags = n_lags,
                          lag_start = lag_start), silent = T)
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
  fitted_vgrm <-
    try(fit.variogram(object = experimental_vgrm, fit.method = 7, debug.level = 0,
                      model = vgm(model = model_vgrm)), silent = T)
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
  k_cv <- 
    try(krige.cv(formula = hf~1, locations = shp_hf, model = fitted_vgrm, nmax = n_max,
                 nfold = n_fold, verbose = F))
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
      if (sum(is.na(k_cv$residual)) >= nrow(k_cv)/2) {
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
  iwt <- interp_weight
  vwt <- vgrm_weight
  vgrm_rmse <- sqrt(attr(fitted_vgrm, "SSErr") / nrow(experimental_vgrm))
  vgrm_cost <- vwt * sqrt(vgrm_rmse / sd(experimental_vgrm$gamma, na.rm = T))
  interp_rmse <- sqrt(sum(k_cv$residual^2, na.rm = T) / nrow(k_cv))
  interp_cost <- iwt * interp_rmse / sd(k_cv$var1.pred, na.rm = T)
  # Return cost = vgrm_error + interp_error
  if (verbose) {
    cat('\nSegment:', segment)
    cat('\nCutoff proportion:', cutoff_prop)
    cat('\nNumber of lags:', n_lags)
    cat('\nLag start:', lag_start)
    cat('\nMax pairs:', n_max)
    cat('\nVariogram weight:', vwt)
    cat('\nVariogram rmse:', vgrm_rmse)
    cat('\nVariogram cost:', vgrm_cost)
    cat('\nInterpolation weight:', iwt)
    cat('\nInterpolation rmse:', interp_rmse)
    cat('\nInterpolation cost:', interp_cost)
    cat('\nCost:', vgrm_cost + interp_cost)
    cat('\nVariogram model:', model_vgrm)
    cat('\n')
    print(fitted_vgrm)
    cat('\n', rep('-', 40), '\n', sep='')
  }
  return(vgrm_cost + interp_cost)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# decode opt !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
decode_opt <- function(model_vgrm = NULL, shp_hf = NULL, opt = NULL) {
  # Check for missing arguments
  if (is.null(model_vgrm)) stop('\nMissing variogram model!')
  if (is.null(shp_hf)) stop('\nMissing heat flow data!')
  if (is.null(opt)) stop('\nMissing nloptr object!')
  # Compute experimental vgrm
  experimental_vgrm <-
    try(experimental_vgrm(shp_hf, cutoff_prop = opt$solution[1], n_lags = opt$solution[2],
                          lag_start = opt$solution[3]), silent = T)
  # Need more than one lag to fit variogram
  if (nrow(experimental_vgrm) < 2) {
    stop('\nExperimental variogram has less than two lags!')
  }
  if (any(class(experimental_vgrm) == 'try-error')) {
    stop('\nExperimental variogram error!')
  }
  # Fit experimental variogram
  fitted_vgrm <-
    try(fit.variogram(object = experimental_vgrm, fit.method = 7, debug.level = 0,
                      model = vgm(model = model_vgrm)), silent = T)
  # Check for error during fitting
  if (any(class(fitted_vgrm) == 'try-error')) {
    print(fitted_vgrm)
    stop('\nVariogram fitting error!')
  }
  return(list('experimental_vgrm' = experimental_vgrm, 'fitted_vgrm' = fitted_vgrm))
}

#######################################################
## .3.            Plotting Functions             !!! ##
#######################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# split lines !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
split_lines <- function(input_lines, max_length, id = NULL) {
  geom_column <- attr(input_lines, "sf_column")
  input_crs <- st_crs(input_lines)
  input_lines[["geom_len"]] <- st_length(input_lines[[geom_column]])
  attr(input_lines[["geom_len"]], "units") <- NULL
  input_lines[["geom_len"]] <- as.numeric(input_lines[["geom_len"]])
  too_long <- 
    input_lines %>%
    select(all_of(id), all_of(geom_column), geom_len) %>%
    filter(geom_len >= max_length)
  rm(input_lines) # just to control memory usage in case this is big.
  too_long <-
    too_long %>%
    mutate(
      pieces = ceiling(geom_len / max_length),
      piece_len = (geom_len / pieces),
      fID = 1:nrow(too_long)
    )
  split_points <-
    st_set_geometry(too_long, NULL)[rep(seq_len(nrow(too_long)), too_long[["pieces"]]),]
  split_points <-
    split_points %>%
    mutate(split_fID = row.names(split_points)) %>%
    select(-geom_len, -pieces) %>%
    group_by(fID) %>%
    mutate(ideal_len = cumsum(piece_len)) %>%
    ungroup()
  coords <- data.frame(st_coordinates(too_long[[geom_column]]))
  rm(too_long)
  coords <- rename(coords, fID = L1) %>% mutate(nID = 1:nrow(coords))
  split_nodes <-
    coords %>%
    group_by(fID) %>%
    # First calculate cumulative length by feature.
    mutate(len  = sqrt(((X - (lag(X)))^2) + (((Y - (lag(Y)))^2)))) %>%
    mutate(len = ifelse(is.na(len), 0, len)) %>%
    mutate(len = cumsum(len)) %>%
    # Now join nodes to split points -- this generates all combinations.
    left_join(select(split_points, fID, ideal_len, split_fID), by = "fID") %>%
    # Calculate the difference between node-wise distance and split-point distance.
    mutate(diff_len = abs(len - ideal_len)) %>%
    # regroup by the new split features.
    group_by(split_fID) %>%
    # filter out na then grab the min distance
    filter(!is.na(diff_len) & diff_len == min(diff_len)) %>%
    ungroup() %>%
    # Grab the start node for each geometry -- the end node of the geometry before it.
    mutate(
      start_nID = lag(nID),
      # need to move the start node one for new features.
      new_feature = fID - lag(fID, default = -1),
      start_nID = ifelse(new_feature == 1, start_nID + 1, start_nID)
    ) %>%
    # Clean up the mess
    select(fID, split_fID, start_nID, stop_nID = nID, -diff_len, -ideal_len, -len, -X, -Y)
  split_nodes$start_nID[1] <- 1
  split_points <- 
    split_points %>%
    left_join(select(split_nodes, split_fID, start_nID, stop_nID), by = "split_fID")
  new_line <- function(start_stop, coords) {
    st_linestring(as.matrix(coords[start_stop[1]:start_stop[2], c("X", "Y")]))
  }
  split_lines <-
    apply(
      as.matrix(split_points[c("start_nID", "stop_nID")]),
      MARGIN = 1,
      FUN = new_line,
      coords = coords
    )
  split_lines <-
    st_sf(split_points[c(id, "split_fID")], geometry = st_sfc(split_lines, crs = input_crs))
  return(split_lines)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# split segment !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
split_segment <-
  function(
    seg_name,
    buf_dir = 'l',
    seg_num = 6,
    buf_len = 5e5,
    sector_exclude = NULL
  ) {
  pnts <- shp_hf_crop[[seg_name]]
  seg <- shp_segs[[seg_name]]
  buf <- st_buffer(seg, dist = buf_len, endCapStyle = 'ROUND')
  volc <- shp_volc
  split_seg <- split_lines(seg, as.numeric(st_length(seg)/seg_num))
  split_buf <-
    map(1:seg_num, ~
      st_buffer(
        split_seg[.x,],
        dist = ifelse(buf_dir == 'l', -1*buf_len, 1*buf_len),
        singleSide = T
      ) %>%
      st_intersection(buf)
    )
  pnts_buf <-
    map(1:seg_num, ~
      st_intersection(pnts, split_buf[[.x]]) %>%
      mutate(
        distance_from_seg = as.vector(st_distance(split_seg[.x,]$geometry, geometry)),
        .before = geometry
      )
    )
  volc_buf <-
    map(1:seg_num, ~
      st_intersection(volc, split_buf[[.x]]) %>%
      mutate(
        distance_from_seg = as.vector(st_distance(split_seg[.x,]$geometry, geometry)),
        .before = geometry
      )
    )
  best_mod <-
    solns %>%
    filter(segment == seg_name) %>%
    slice_min(cost)
  interp <- best_mod[['shp_interp_diff']][[1]]
  interp_buf <-
    map(1:seg_num, ~
      st_intersection(interp, split_buf[[.x]]) %>%
      mutate(
        distance_from_seg = as.vector(st_distance(split_seg[.x,]$geometry, geometry)),
        .before = geometry
      )
    )
  if (!is.null(sector_exclude)) {
    split_buf <- split_buf[-sector_exclude]
    pnts_buf <- pnts_buf[-sector_exclude]
    volc_buf <- volc_buf[-sector_exclude]
    interp_buf <- interp_buf[-sector_exclude]
  }
  return(
    list(
      'seg' = split_seg,
      'buf' = split_buf,
      'pnts' = pnts_buf,
      'volc' = volc_buf,
      'interp' = interp_buf
    )
  )
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot split segment !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_split_segment <-
  function(
    split_seg,
    running_avg = 3,
    borders =
      c(
        'left' = 0.1,
        'right' = 0.1,
        'top' = 0.1,
        'bottom' = 0.1
      )
  ) {
  seg_name <- split_seg$interp[[1]]$segment[1]
  seg_num <- as.numeric(unique(bind_rows(split_seg$pnts)$split_fID))
  bx <- bbox_widen(st_bbox(st_buffer(st_combine(split_seg$seg), dist = 5e5)), borders = borders)
  seg <- bind_rows(split_seg$seg)
  ridge <- shp_ridge_crop[[seg_name]]
  trench <- shp_trench_crop[[seg_name]]
  transform <- shp_transform_crop[[seg_name]]
  relief <- shp_relief_crop[[seg_name]]
  volc <- split_seg$volc
  wdth <- range(st_bbox(bind_rows(split_seg$buf))[c('xmin', 'xmax')])/1e3
  rng <-
    bind_rows(split_seg$interp) %>%
    st_set_geometry(NULL) %>%
    select(est_sim, est_krige, split_fID) %>%
    summarise(
      min_krige = min(est_krige),
      max_krige = max(est_krige),
      min_sim = min(est_sim),
      max_sim = max(est_sim)
    )
  pnt_size <- 1
  p0 <-
    ggplot() +
    geom_sf(data = bx, color = NA, fill = NA) +
    geom_sf(data = relief, aes(color = elevation), shape = 15, size = 0.25) +
    scale_color_etopo(guide = 'none') +
    new_scale_color() +
    geom_sf(data = ridge, linewidth = 0.5) +
    geom_sf(data = trench, linewidth = 0.5) +
    geom_sf(data = transform, linewidth = 0.5) +
    geom_sf(data = seg, linewidth = 1, color = 'white') +
    geom_sf(
      data = bind_rows(split_seg$buf),
      aes(fill = factor(split_fID, levels = seg_num[order(seg_num)])),
      linewidth = 0.5,
      color = 'black',
      show_legend = F
    ) +
    geom_sf(
      data = bind_rows(split_seg$pnts),
      aes(
        fill = factor(split_fID, levels = seg_num[order(seg_num)]),
        group = factor(split_fID, levels = seg_num[order(seg_num)])
      ),
      size = pnt_size,
      shape = 22,
      stroke = 0.5,
      show_legend = F
    ) +
    geom_sf(data = bind_rows(volc), color = 'white', shape = 18, size = pnt_size) +
    annotate(
      'label',
      label = 'a',
      x = -Inf,
      y = Inf,
      size = 7,
      hjust = 0,
      vjust = 1,
      fill = 'grey90',
      label.padding = unit(0.02, 'in'),
      label.r = unit(0, 'in')
    ) +
    coord_sf(
      xlim = c(st_bbox(bx)$xmin, st_bbox(bx)$xmax),
      ylim = c(st_bbox(bx)$ymin, st_bbox(bx)$ymax),
      label_axes = '--EN'
    ) +
    scale_color_discrete_qualitative('set 2') +
    scale_fill_discrete_qualitative('set 2') +
    guides(
      fill = 
        guide_legend(
          nrow = 1,
          override.aes = list(alpha = 1, size = 8),
          title.position = 'top',
          title.vjust = 1,
          label.position = 'bottom'
        )
    ) +
    theme_map(font_size = 9) +
    theme(
      plot.margin = margin(1, 1, 1, 5),
      axis.text.x = element_text(color = 'grey40', angle = 30, hjust = 0, vjust = 0),
      axis.text.y = element_text(color = 'grey40', angle = 30, hjust = 0),
      panel.grid = element_line(size = 0.01, color = 'grey60')
    )
    if (
       seg_name %in% c('Alaska Aleutians', 'Kamchatka Marianas', 'Tonga New Zealand', 'Vanuatu')
     ) {
      p0 <-
        p0 +
        scale_x_continuous(
          breaks = c(
            80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, -170, -160, -150, -140, -130
          )
        )
    }
  p1 <-
    bind_rows(split_seg$interp) %>%
    st_set_geometry(NULL) %>%
    select(est_sim, est_krige, split_fID) %>%
    rename(Similarity = est_sim, Kriging = est_krige) %>%
    pivot_longer(-split_fID) %>%
    group_by(name) %>%
    filter(split_fID %in% seg_num) %>%
    ggplot() +
    geom_boxplot(
      aes(
        y = value,
        x = factor(split_fID, levels = seg_num[order(seg_num)]),
        fill = factor(split_fID, levels = seg_num[order(seg_num)]),
        color = name,
        outlier.color = name
      ),
      outlier.size = pnt_size * 0.2,
      size = 0.5
    ) +
    annotate(
      'label',
      label = 'b',
      x = -Inf,
      y = Inf,
      size = 7,
      hjust = 0,
      vjust = 1,
      fill = 'grey90',
      label.padding = unit(0.02, 'in'),
      label.r = unit(0, 'in')
    ) +
    guides(
      color =
        guide_legend(
          nrow = 1,
          override.aes = list(alpha = 1, fill = 'grey50'),
          title.position = 'top',
          title.vjust = 1,
          label.position = 'left'
        ),
      fill = 'none'
    ) +
    coord_cartesian(
      ylim = c(min(rng$min_krige, rng$min_sim), max(rng$max_krige, rng$max_sim))
    ) +
    scale_y_continuous(position = 'right') +
    scale_fill_discrete_qualitative(palette = 'set 2') +
    scale_color_manual(values = c('black', 'white')) +
    labs(y = bquote('heat flow'~(mWm^-2)), x = NULL, color = 'interpolation method', fill = NULL) +
    theme_dark(base_size = 14) +
    theme(
      legend.title = element_text(margin = margin(0, 0, -5, 0)),
      plot.margin = margin(1, 1, 1, 1),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = 'grey50')
    )
  labl <-
    tibble(
      name = c('Kriging', 'Similarity'),
      lbl = c('c', 'd')
    )
  p2 <-
    bind_rows(split_seg$interp) %>%
    st_set_geometry(NULL) %>%
    filter(distance_from_seg <= 500000) %>%
    select(est_sim, est_krige, split_fID, distance_from_seg) %>%
    arrange(distance_from_seg) %>%
    mutate(
      Similarity = rollmean(est_sim, running.avg, fill = NA),
      Kriging = rollmean(est_krige, running.avg, fill = NA)
    ) %>%
    select(-c(est_sim, est_krige)) %>%
    pivot_longer(-c(split_fID, distance_from_seg)) %>%
    group_by(name) %>%
    ggplot() +
    geom_point(
      data = bind_rows(volc),
      aes(distance_from_seg/1e3, min(rng$min_krige, rng$min_sim)),
      color = 'white',
      size = pnt_size,
      shape = 18,
      key_glyph = 'rect'
    ) +
    geom_point(
      data = bind_rows(split_seg$pnts),
      aes(
        distance_from_seg/1e3,
        hf,
        fill = factor(split_fID, levels = seg_num[order(seg_num)]),
        group = factor(split_fID, levels = seg_num[order(seg_num)])
      ),
      size = pnt_size,
      shape = 22,
      key_glyph = 'rect'
    ) +
    geom_smooth(
      aes(
        distance_from_seg/1e3,
        value,
        color = factor(split_fID, levels = seg_num[order(seg_num)]),
        group = factor(split_fID, levels = seg_num[order(seg_num)])
      ),
      fill = NA,
      method = 'loess',
      #span = 0.95,
      alpha = 0.1,
      size = 1,
      se = T,
      key_glyph = 'rect',
      show.legend = F
    ) +
    geom_label(
      data = labl,
      aes(label = lbl),
      x = -Inf,
      y = Inf,
      size = 7,
      hjust = 0,
      vjust = 1,
      fill = 'grey90',
      label.padding = unit(0.02, 'in'),
      label.r = unit(0, 'in')
    ) +
    labs(
      x = 'distance from trench (km)',
      y = bquote('heat flow'~(mWm^-2)),
      fill = 'sector'
    ) +
    guides(
      fill = 
        guide_legend(
          nrow = 1,
          override.aes = list(alpha = 1, size = 5),
          title.position = 'top',
          title.vjust = 1,
          label.position = 'bottom'
        )
    ) +
    coord_cartesian(
      ylim = c(min(rng$min_krige, rng$min_sim), max(rng$max_krige, rng$max_sim))
    ) +
    scale_color_discrete_qualitative('set 2') +
    scale_fill_discrete_qualitative('set 2') +
    facet_wrap(~name, ncol = 2) +
    theme_dark(base_size = 14) +
    theme(
      legend.title = element_text(margin = margin(0, 0, -5, 0)),
      strip.background = element_rect(color = NA, fill = NA),
      strip.text = element_text(color = 'black', size = 14, margin = margin(0, 0, 2, 0)),
      strip.placement = 'outside',
      legend.box.margin = margin(0, 0, 0, 0),
      legend.key = element_rect(color = 'black'),
      legend.spacing.x = unit(0, 'mm'),
      plot.margin = margin(1, 1, 1, 1),
      panel.grid = element_blank()
    )
  p <-
    ((p0 | p1) / p2) +
    plot_annotation(title = 'Comparing heat flow interpolations by sector') +
    plot_layout(widths = 1, guides = 'collect') &
    theme(
      plot.title = element_text(size = 18),
      plot.margin = margin(5, 5, 5, 5),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.margin = margin(),
      legend.position = 'bottom',
      legend.direction = 'horizontal',
      legend.justification = 'left'
    )
  p
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot vgrm !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_vgrm <-
  function(
    experimental_vgrm,
    fitted_vgrm = NULL,
    cost = NULL,
    v_mod = NULL,
    lineCol = 'white',
    ylim = NULL,
    xlim = NULL
  ) {
  # Check for missing arguments
  if (is.null(experimental_vgrm)) stop('\nMissing experimental variogram!')
  sill <- round(sqrt(fitted_vgrm[[2]]))
  range <- round(fitted_vgrm[[3]]/1e3)
  p <- 
    ggplot() +
    geom_segment(
      aes(x = range, xend = range, y = Inf, yend = sill, color = 'range'),
      arrow = arrow(angle = 10, length = unit(0.05, 'inches'), type = 'closed')
    ) +
    geom_segment(
      aes(x = -Inf, xend = range, y = sill, yend = sill, color = 'sill'),
      arrow = arrow(angle = 10, length = unit(0.05, 'inches'), type = 'closed')
    ) +
    labs(x = 'lag distance (km)', y = bquote(sqrt(hat(gamma)[(h)])~(mWm^-2))) +
    scale_color_discrete_qualitative(name = NULL, palette = 'set 2') +
    new_scale_color() +
    coord_cartesian(ylim = ylim, xlim = xlim) +
    theme_dark(base_size = 14) +
    theme(panel.grid = element_blank())
  plt <- tryCatch(
    {
      p +
      geom_point(
        data = experimental_vgrm,
        aes(x = dist/1e3, y = sqrt(gamma), color = 'experimental variogram'),
        shape = 20,
        key_glyph = 'path'
      ) +
      geom_path(
        data = experimental_vgrm,
        aes(x = dist/1e3, y = sqrt(gamma), color = 'experimental variogram'),
        key_glyph = 'path'
      ) +
      geom_line(
        data = variogramLine(fitted_vgrm, maxdist = max(experimental_vgrm$dist)),
        aes(x = dist/1e3, y = sqrt(gamma), color = 'model variogram'),
        key_glyph = 'path'
      ) +
      annotate(
        'text',
        label = paste0(
          'fitted variogram: ', v_mod,
          '\nsill: ', sill, ', ',
          'range: ', range, ', ',
          'cost: ', round(cost, 4)
        ),
        x = xlim[2]/2,
        y = -Inf,
        vjust = -0.2,
        hjust = 0.5
      ) +
      scale_color_manual(
        name = NULL,
        values = c('experimental variogram' = 'black', 'model variogram' = lineCol)
      )
    },
    error=function(cond) {
      # If there was an error with the variogram fit ...
      # plot only the experimental variogram
      cat('Somthing is up with variogram model!')
      cat('\nPlotting experimental variogram only\n')
      exp_plt <- p +
      geom_point(
        data = experimental_vgrm,
        aes(x = dist/1e3, y = sqrt(gamma)),
        shape = 20
      )
      return(exp_plt)
    }
  )
  return(plt)
}
