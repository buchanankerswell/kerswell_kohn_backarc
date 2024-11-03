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
package_list <- c(
  "tictoc", "stringr", "tidyr", "readr", "readxl", "purrr", "furrr",
  "tibble", "dplyr", "magrittr", "units", "ggplot2", "colorspace", "metR",
  "ggrepel", "ggridges", "ggnewscale", "patchwork", "cowplot", "ggsflabel",
  "marmap", "scales", "ggspatial", "gstat", "sp", "sf", "rnaturalearth",
  "nloptr", "zoo", "jsonlite"
)

# Load packages quietly
sapply(package_list, sshhh)
rm(package_list, sshhh)

# Turn off S2 geometry
suppressMessages(sf_use_s2(F))

# Set seed
seed <- 42
set.seed(seed)

# Set map projections
wgs <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +no_defs"
prj <- "+proj=eck4 +lon_0=-180 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"

#######################################################
## .1.         General Helper Functions          !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read json files to df !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read multiple JSON files into a single data frame
read_json_files_to_df <- function(files) {
  safely_read_json <- safely(fromJSON)
  results <- map(files, safely_read_json)
  data_list <- map(results, "result") %>% compact()
  if (length(data_list) == 0) {
    return(tibble())
  }
  bind_rows(data_list)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# bbox extend !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bbox_extend <- function(bbox, ext = rep(0, 4), square = F) {
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
  st_polygon(list(matrix(box, ncol = 2, byrow = T))) %>% st_sfc(crs = st_crs(bbox))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# reproject center pacific !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
reproject_center_pacific <- function(shp, break_dateline = T, tol = 1) {
  if (st_crs(shp) != st_crs(wgs)) {
    shp <- st_transform(shp, wgs)
  }
  if (any(st_geometry_type(shp) %in% "POINT") || !break_dateline) {
    shp %>%
      st_make_valid() %>%
      st_transform(prj)
  } else {
    lon0 <- as.numeric(str_extract(prj, "(?<=lon_0=)-?\\d+")) + 360
    suppressWarnings({
      suppressMessages({
        shp %>%
          st_make_valid() %>%
          st_wrap_dateline(options = c("WRAPDATELINE=YES", paste0("DATELINEOFFSET=", lon0))) %>%
          st_break_antimeridian(lon_0 = lon0, tol = tol) %>%
          st_transform(prj)
      })
    })
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# reproject wgs !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
reproject_wgs <- function(shp) {
  if (st_crs(shp) != st_crs(wgs)) {
    shp <- st_transform(shp, wgs)
    shp %>% st_make_valid()
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# crop feature !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
crop_feature <- function(ft, bbox, crop_within = F, keep_df = F) {
  if (st_crs(ft) != st_crs(bbox)) {
    ft <- st_transform(ft, st_crs(bbox))
  }
  if (!crop_within) {
    suppressWarnings({
      suppressMessages({
        ft_cropped <- st_crop(ft, bbox)
      })
    })
  } else {
    suppressWarnings({
      suppressMessages({
        ft_cropped <- st_intersection(ft, bbox)
      })
    })
  }
  if (!is.data.frame(ft_cropped)) {
    l <- length(ft_cropped)
  } else {
    l <- nrow(ft_cropped)
  }
  if (l == 0) {
    if (!keep_df) {
      NA
    } else {
      NULL
    }
  } else {
    if (!keep_df) {
      suppressWarnings({
        suppressMessages({
          st_sfc(st_union(ft_cropped), crs = st_crs(bbox))
        })
      })
    } else {
      ft_cropped
    }
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get world bathy !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_world_bathy <- function(res = 15, path = "assets/map_data/relief/") {
  if (!dir.exists(path)) {
    dir.create(path, recursive = T, showWarnings = F)
  }
  tryCatch(
    {
      suppressWarnings({
        suppressMessages({
          getNOAA.bathy(180, -180, 90, -90, res, T, F, path) %>%
            as.SpatialGridDataFrame() %>%
            st_as_sf(crs = wgs) %>%
            rename(elev = layer) %>%
            reproject_center_pacific()
        })
      })
    },
    error = function(e) {
      cat("\n!! ERROR occurred in get_world_bathy:\n!!", conditionMessage(e))
      return(NULL)
    }
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get seg bathy !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_seg_bathy <- function(shp, res = 2, path = "assets/map_data/relief/", tol = 1) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = T, showWarnings = F)
  }
  tryCatch(
    {
      bbx <- shp %>%
        st_transform(wgs) %>%
        st_bbox() %>%
        round(2)
      if (bbx[3] - bbx[1] > 180) {
        antim <- T
      } else {
        antim <- F
      }
      getNOAA.bathy(bbx[3], bbx[1], bbx[2], bbx[4], res, T, antim, path) %>%
        as.SpatialGridDataFrame() %>%
        st_as_sf(crs = wgs) %>%
        st_make_valid() %>%
        st_transform(prj) %>%
        rename(elev = layer)
    },
    error = function(e) {
      cat("\n!! ERROR occurred in get_seg_bathy:\n!!", conditionMessage(e))
      return(NULL)
    }
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# parse zerodist obs !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    if (obs1 == obs2) {
      return(locr)
    }
    U1 <- str_extract(df$Quality_Code[loc1], "(?<=U).")
    u1 <- if (U1 == "x") {
      5
    } else {
      as.numeric(U1)
    }
    U2 <- str_extract(df$Quality_Code[loc2], "(?<=U).")
    u2 <- if (U2 == "x") {
      5
    } else {
      as.numeric(U2)
    }
    M1 <- str_extract(df$Quality_Code[loc1], "(?<=M).")
    m1 <- if (M1 == "x") {
      5
    } else {
      as.numeric(M1)
    }
    M2 <- str_extract(df$Quality_Code[loc2], "(?<=M).")
    m2 <- if (M2 == "x") {
      5
    } else {
      as.numeric(M2)
    }
    sum1 <- sum(u1 + m1, na.rm = T)
    sum2 <- sum(u2 + m2, na.rm = T)
    if (sum1 < sum2) {
      return(loc2)
    }
    if (sum2 < sum1) {
      return(loc1)
    }
    if (sum2 == sum1) {
      return(locr)
    }
  })
  slice(df, -rid)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate interp rmse !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
calculate_interp_rmse <- function(obs, interp, type = "similarity") {
  if (!(type %in% c("krige", "similarity"))) stop("\ninvalid type!")
  if (type == "similarity") {
    nearest_est <- interp$est_sim[st_nearest_feature(obs, st_geometry(interp))]
    cat("\n", nearest_est)
    sqrt(mean((nearest_est - obs$obs)^2))
  } else if (type == "krige") {
    interpolation <- interpolation
    nearest_est <- interp$est_krige[st_nearest_feature(obs, st_geometry(interp))]
    sqrt(mean((nearest_est - obs$obs)^2))
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# project obs to transect !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
project_obs_to_transect <- function(transect, shp_obs) {
  if (is.null(shp_obs)) {
    return(NULL)
  }
  if (st_crs(transect) != st_crs(shp_obs)) {
    stop("transect and shp_obs crs not the same!")
  }
  if (!any(class(transect) == "sfc")) {
    stop("\ntransect needs to be an sf object!")
  }
  if (!any(class(shp_obs) == "sf")) {
    stop("\nUnrecognized shp_obs passed to project_obs_to_transect() !")
  }
  if (!any(names(shp_obs) %in% c("obs", "est_sim", "est_krg"))) {
    stop("\nUnrecognized sf object passed to project_obs_to_transect() !")
  }
  if (any(names(shp_obs) == "ghf")) {
    obs <- shp_obs$obs
    sigma <- NA
  } else if (any(names(shp_obs) == "similarity")) {
    obs <- shp_obs$est_sim
    sigma <- shp_obs$sigma_sim
  } else if (any(names(shp_obs) == "krige")) {
    obs <- shp_obs$est_krg
    sigma <- shp_obs$sigma_krg
  }
  projected_distances <- st_line_project(transect, st_geometry(shp_obs), normalized = T)
  st_as_sf(st_line_interpolate(transect, projected_distances)) %>%
    rename(geometry = x) %>%
    mutate(projected_distances = projected_distances, obs = obs, sigma = sigma, .before = geometry)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# fit loess to projected obs !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit_loess_to_projected_obs <- function(df, n = 1e3, nmin = 10, span = 0.65) {
  if (is.null(df)) {
    return(NULL)
  }
  df <- df %>% st_set_geometry(NULL)
  if (nrow(df) < nmin) {
    cat("\n   Cannot fit loess with less than", nmin, "projected obs!")
    return(NULL)
  } else {
    mod <- NULL
    while (span <= 0.9 && is.null(mod)) {
      mod <- tryCatch(
        {
          loess(obs ~ projected_distances, data = df, span = span)
        },
        error = function(e) {
          NULL
        }
      )
      if (is.null(mod)) {
        span <- span + 0.05
        cat("\n   Loess failed with", nrow(df), "obs! Increasing span to", span)
      }
    }
    if (is.null(mod) || is.null(mod$fitted)) {
      cat("\n   Loess fitting failed!")
      return(NULL)
    }
    new_dist <- seq(0, 1, length.out = n)
    original_range <- range(df$projected_distances, na.rm = TRUE)
    loess_pred <- tryCatch(
      {
        predict(mod, newdata = new_dist)
      },
      error = function(e) {
        cat("\n   Prediction failed!")
        return(NULL)
      }
    )
    if (is.null(loess_pred)) {
      return(NULL)
    }
    loess_pred[new_dist < original_range[1] | new_dist > original_range[2]] <- NA
    smooth <- tibble(projected_distances = new_dist, obs = loess_pred)
    return(smooth)
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate cross correlation !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
calculate_cross_correlation <- function(smooth1, smooth2) {
  len1 <- nrow(smooth1)
  len2 <- nrow(smooth2)
  if (len1 < len2) {
    smooth1 <- bind_rows(smooth1, tibble(
      projected_distances = rep(NA, len2 - len1),
      obs = rep(NA, len2 - len1)
    ))
  } else if (len2 < len1) {
    smooth2 <- bind_rows(smooth2, tibble(
      projected_distances = rep(NA, len1 - len2),
      obs = rep(NA, len1 - len2)
    ))
  }
  tryCatch(
    {
      ccf_result <- ccf(smooth1$obs, smooth2$obs, plot = F)
      max_ccf <- max(ccf_result$acf)
      max_lag <- ccf_result$lag[which.max(ccf_result$acf)]
      return(max_ccf)
    },
    error = function(e) {
      cat("\n!! ERROR occurred in calculate_cross_correlation:\n!!", conditionMessage(e))
      return(NULL)
    }
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load map data !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load_map_data <- function(fpath) {
  load(fpath, envir = parent.frame())
  if (!exists("shp_submap", envir = parent.frame())) {
    stop('\nMissing map data! Use "load(path/to/map-data.RData)"')
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load nlopt interpolation data !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load_nlopt_interpolation_data <- function(fpath) {
  load(fpath, envir = parent.frame())
  if (!exists("nlopt_summary", envir = parent.frame())) {
    stop('\nMissing nlopt data! Use "load(path/to/interpolation-summary.RData)"')
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generate global rects !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
generate_global_rects <- function(spacing = 10) {
  lon_breaks <- seq(-180, 180, by = spacing)
  lat_breaks <- seq(-90, 90, by = spacing)
  grid <- expand.grid(
    lon_min = lon_breaks[-length(lon_breaks)],
    lat_min = lat_breaks[-length(lat_breaks)]
  )
  grid$lon_max <- grid$lon_min + spacing
  grid$lat_max <- grid$lat_min + spacing
  f <- function(lon_min, lat_min, lon_max, lat_max) {
    mat <- matrix(c(
      lon_min, lat_min, lon_min, lat_max, lon_max,
      lat_max, lon_max, lat_min, lon_min, lat_min
    ), ncol = 2, byrow = TRUE)
    st_polygon(list(mat))
  }
  polys <- pmap(grid, f)
  st_as_sf(tibble(rect = st_sfc(polys, crs = wgs))) %>%
    mutate(id = row_number(), short_name = paste0("rect", row_number()), .before = rect)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# compile rect data !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
compile_rect_data <- function(rect_ids = NULL, fv = NULL, np = NULL) {
  load_map_data("assets/map_data/map-data.RData")
  if (is.null(rect_ids)) {
    stop("\nMissing rect ids!")
  }
  if (all(rect_ids >= 1) && all(rect_ids <= nrow(shp_global_rects))) {
    x <- slice(shp_global_rects, rect_ids)
  } else if (is.character(rect_ids) && all(rect_ids %in% shp_global_rects$short_name)) {
    x <- filter(shp_global_rects, short_name %in% rect_ids)
  } else {
    stop("\nUnrecognized input for rect_ids! Use the rect short_name or id number ...")
  }
  if (nrow(x) == 0) {
    stop("\n", rect_ids, " not found in shp_global_rects!")
  }
  suppressWarnings({
    suppressMessages({
      df <-
        x %>%
        rowwise() %>%
        mutate(ridge = crop_feature(shp_ridge, rect), .before = rect) %>%
        mutate(trench = crop_feature(shp_trench, rect), .before = rect) %>%
        mutate(transform = crop_feature(shp_transform, rect), .before = rect) %>%
        mutate(volcano = crop_feature(select(shp_volc, geometry), rect), .before = rect) %>%
        mutate(grid = list(crop_feature(shp_grid, rect, T, T)), .before = rect) %>%
        mutate(ghf_rect = list(crop_feature(shp_ghf_raw_wgs, rect, T, T)), .before = rect) %>%
        mutate(sim_rect = list(crop_feature(shp_sim_raw_wgs, rect, T, T)), .before = rect) %>%
        mutate(bathy = list(get_seg_bathy(rect)), .before = rect)
    })
  })
  if (!is.null(fv) && !is.null(np)) {
    shp_krg <- Krige(df$ghf_rect[[1]], fv, df$grid[[1]], np)
    df %>%
      mutate(krg_rect = list(crop_feature(shp_krg, rect, T, T))) %>%
      mutate(dff_rect = list(interp_diff(krg_rect, sim_rect))) %>%
      ungroup()
  } else {
    ungroup(df)
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# compile transect data !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
compile_trans_data <- function(trans_ids = NULL, lbuff = 5e5, sbuff = c(5e4, 1e5, 1.5e5),
                               fv = NULL, np = NULL) {
  load_map_data("assets/map_data/map-data.RData")
  if (is.null(trans_ids)) {
    stop("\nMissing submap transect ids!")
  }
  if (is.numeric(trans_ids) && all(trans_ids >= 1) && all(trans_ids <= nrow(shp_submap))) {
    x <- slice(shp_submap, trans_ids)
  } else if (is.character(trans_ids) && all(trans_ids %in% shp_submap$short_name)) {
    x <- filter(shp_submap, short_name %in% trans_ids)
  } else {
    stop("\nUnrecognized input for trans_ids! Use the transect short_name or id number ...")
  }
  if (nrow(x) == 0) {
    stop("\n", trans_ids, " not found in shp_submap!")
  }
  suppressWarnings({
    suppressMessages({
      df <-
        x %>%
        rowwise() %>%
        mutate(large_buffer = st_buffer(transect, lbuff, endCapStyle = "ROUND")) %>%
        mutate(small_buffer1 = st_buffer(transect, sbuff[1], endCapStyle = "FLAT")) %>%
        mutate(small_buffer2 = st_buffer(transect, sbuff[2], endCapStyle = "FLAT")) %>%
        mutate(small_buffer3 = st_buffer(transect, sbuff[3], endCapStyle = "FLAT")) %>%
        mutate(bbox = bbox_extend(st_bbox(large_buffer), square = T)) %>%
        mutate(ridge = crop_feature(shp_ridge, bbox)) %>%
        mutate(trench = crop_feature(shp_trench, bbox)) %>%
        mutate(transform = crop_feature(shp_transform, bbox)) %>%
        mutate(volcano = crop_feature(select(shp_volc, geometry), bbox)) %>%
        mutate(grid = list(crop_feature(shp_grid, large_buffer, T, T))) %>%
        mutate(ghf_large_buff = list(crop_feature(shp_ghf, large_buffer, T, T))) %>%
        mutate(ghf_small_buff1 = list(crop_feature(shp_ghf, small_buffer1, T, T))) %>%
        mutate(ghf_small_buff2 = list(crop_feature(shp_ghf, small_buffer2, T, T))) %>%
        mutate(ghf_small_buff3 = list(crop_feature(shp_ghf, small_buffer3, T, T))) %>%
        mutate(ghf_projected1 = list(project_obs_to_transect(transect, ghf_small_buff1))) %>%
        mutate(ghf_projected2 = list(project_obs_to_transect(transect, ghf_small_buff2))) %>%
        mutate(ghf_projected3 = list(project_obs_to_transect(transect, ghf_small_buff3))) %>%
        mutate(ghf_loess1 = list(fit_loess_to_projected_obs(ghf_projected1))) %>%
        mutate(ghf_loess2 = list(fit_loess_to_projected_obs(ghf_projected2))) %>%
        mutate(ghf_loess3 = list(fit_loess_to_projected_obs(ghf_projected3))) %>%
        mutate(sim_large_buff = list(crop_feature(shp_sim, large_buffer, T, T))) %>%
        mutate(sim_small_buff1 = list(crop_feature(shp_sim, small_buffer1, T, T))) %>%
        mutate(sim_small_buff2 = list(crop_feature(shp_sim, small_buffer2, T, T))) %>%
        mutate(sim_small_buff3 = list(crop_feature(shp_sim, small_buffer3, T, T))) %>%
        mutate(sim_projected1 = list(project_obs_to_transect(transect, sim_small_buff1))) %>%
        mutate(sim_projected2 = list(project_obs_to_transect(transect, sim_small_buff2))) %>%
        mutate(sim_projected3 = list(project_obs_to_transect(transect, sim_small_buff3))) %>%
        mutate(sim_loess1 = list(fit_loess_to_projected_obs(sim_projected1))) %>%
        mutate(sim_loess2 = list(fit_loess_to_projected_obs(sim_projected2))) %>%
        mutate(sim_loess3 = list(fit_loess_to_projected_obs(sim_projected3))) %>%
        mutate(bathy = list(get_seg_bathy(bbox)))
    })
  })
  if (!is.null(fv) && !is.null(np)) {
    shp_krg <- Krige(df$ghf_large_buff[[1]], fv, df$grid[[1]], np)
    df %>%
      mutate(krg_large_buff = list(crop_feature(shp_krg, large_buffer, T, T))) %>%
      mutate(krg_small_buff1 = list(crop_feature(shp_krg, small_buffer1, T, T))) %>%
      mutate(krg_small_buff2 = list(crop_feature(shp_krg, small_buffer2, T, T))) %>%
      mutate(krg_small_buff3 = list(crop_feature(shp_krg, small_buffer3, T, T))) %>%
      mutate(krg_projected1 = list(project_obs_to_transect(transect, krg_small_buff1))) %>%
      mutate(krg_projected2 = list(project_obs_to_transect(transect, krg_small_buff2))) %>%
      mutate(krg_projected3 = list(project_obs_to_transect(transect, krg_small_buff3))) %>%
      mutate(krg_loess1 = list(fit_loess_to_projected_obs(krg_projected1))) %>%
      mutate(krg_loess2 = list(fit_loess_to_projected_obs(krg_projected2))) %>%
      mutate(krg_loess3 = list(fit_loess_to_projected_obs(krg_projected3))) %>%
      mutate(dff_large_buff = list(interp_diff(krg_large_buff, sim_large_buff))) %>%
      mutate(dff_small_buff1 = list(interp_diff(krg_small_buff1, sim_small_buff1))) %>%
      mutate(dff_small_buff2 = list(interp_diff(krg_small_buff2, sim_small_buff2))) %>%
      mutate(dff_small_buff3 = list(interp_diff(krg_small_buff3, sim_small_buff3))) %>%
      mutate(dff_projected1 = list(interp_diff(krg_projected1, sim_projected1))) %>%
      mutate(dff_projected2 = list(interp_diff(krg_projected2, sim_projected2))) %>%
      mutate(dff_projected3 = list(interp_diff(krg_projected3, sim_projected3))) %>%
      ungroup()
  } else {
    ungroup(df)
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get closest interp obs !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_closest_interp_obs <- function(trans_id = NULL, fv = NULL, np = NULL, thresh = 1e4) {
  if (is.null(trans_id)) {
    stop("\nMissing submap transect ids!")
  }
  if (!is.null(fv) && !is.null(np)) {
    x <- compile_trans_data(trans_id, fv = fv, np = np)
    itp <- x$krg_large_buff[[1]]
  } else {
    x <- compile_trans_data(trans_id)
    itp <- x$sim_large_buff[[1]]
  }
  grid <- x$grid[[1]]
  obs <- x$ghf_large_buff[[1]]
  nearest_obs <- st_nearest_feature(grid, obs)
  nearest_itp <- st_nearest_feature(obs, grid)
  dt_obs <- st_distance(grid, obs[nearest_obs, ], by_element = T) < set_units(thresh, "m")
  dt_itp <- st_distance(obs, grid[nearest_itp, ], by_element = T) < set_units(thresh, "m")
  obs_n <- obs[nearest_obs, ][dt_obs, ]
  itp_n <- itp[nearest_itp, ][dt_itp, ]
  if (any(names(itp_n) %in% c("est_sim", "similarity"))) {
    itp_n %>%
      mutate(
        obs_ghf = obs_n[st_nearest_feature(itp_n, obs_n), ]$obs,
        ghf = obs_n[st_nearest_feature(itp_n, obs_n), ]$ghf, .before = similarity
      ) %>%
      select(-c(sigma_sim, obs_sim))
  } else if (any(names(itp_n) %in% c("est_krg", "krige"))) {
    itp_n %>%
      mutate(
        obs_ghf = obs_n[st_nearest_feature(itp_n, obs_n), ]$obs,
        ghf = obs_n[st_nearest_feature(itp_n, obs_n), ]$ghf, .before = krige
      ) %>%
      select(-c(sigma_krg, var_krg))
  } else {
    NULL
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# safe execute !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
safe_execute <- function(expr) {
  error_message <- paste0("\nError in expression: ", deparse(substitute(expr)))
  tryCatch(expr, error = function(e) {
    warning(error_message, "\n", conditionMessage(e), "\n")
    return(NULL)
  })
}

#######################################################
## .2.             Kriging Functions             !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# fit variogram model !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit_variogram_model <- function(ev, v_mod) {
  if (nrow(ev) < 2) {
    return(NULL)
  }
  fv <- safe_execute(fit.variogram(ev, vgm(model = v_mod), fit.method = 6))
  if (!is.null(fv) && fv$range < 0) {
    warning("Variogram range is negative: ", fv$range)
    return(NULL)
  }
  return(fv)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# krige cv !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
krige_cv <- function(shp_hf, fv, n_max, max_dist, nfold) {
  k_cv <- safe_execute(
    krige.cv(obs ~ 1, shp_hf, model = fv, nmax = n_max, maxdist = max_dist, nfold = nfold)
  )
  return(k_cv)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make experimental vgrm !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
make_experimental_vgrm <- function(shp_hf = NULL, cutoff = 3, n_lags = 30) {
  if (is.null(shp_hf)) {
    stop("\nMissing heat flow data!")
  }
  bbox <- st_bbox(shp_hf)
  bbox_diagonal_distance <- sqrt((bbox$xmax - bbox$xmin)^2 + (bbox$ymax - bbox$ymin)^2)
  lag_cutoff <- as.vector(bbox_diagonal_distance / cutoff)
  bin_width <- lag_cutoff / n_lags
  vgrm <- safe_execute(
    variogram(obs ~ 1, locations = shp_hf, cutoff = lag_cutoff, width = bin_width)
  )
  return(vgrm)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sample and shuffle equal nfolds !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sample_and_shuffle_equal_nfolds <- function(nrow_data, nfold = 10) {
  if (is.null(nrow_data)) {
    stop("\nMissing nrow data!")
  }
  if (is.null(nfold)) {
    nfold <- 10
  }
  fold_size <- nrow_data %/% nfold
  remainder <- nrow_data %% nfold
  fold <- rep(1:nfold, each = fold_size)
  if (remainder > 0) {
    fold <- c(fold, sample(1:nfold, remainder))
  }
  sample(fold)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# log optimization !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log_optimization <- function(log_file, details) {
  if (!dir.exists(dirname(log_file))) {
    dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  }
  sink(log_file, append = TRUE)
  cat(rep("-", 60), "\n", sep = "")
  details_str <- paste(names(details), details, sep = ": ", collapse = "\n")
  cat(details_str, "\n", sep = "")
  sink()
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# cost function !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cost_function <- function(shp_hf, cutoff = 3, n_lags = 50, n_max = 10, max_dist = 1e5,
                          v_mod = "Sph", interp_weight = 0.5, vgrm_weight = 0.5, id = NULL,
                          nfold = NULL, na_thresh_prop = 0.5) {
  if (is.null(shp_hf)) stop("Missing heat flow data model!")
  ## TODO: How to handle cases where nlopt converges even though
  ##       none of the iterations get to log_optimization??
  rand_high_cost <- 1e6 + runif(1, 0, 1e5)
  ev <- make_experimental_vgrm(shp_hf, cutoff, n_lags)
  if (is.null(ev)) {
    return(rand_high_cost)
  }
  fv <- fit_variogram_model(ev, v_mod)
  if (is.null(fv)) {
    return(rand_high_cost)
  }
  k_cv <- krige_cv(shp_hf, fv, n_max, max_dist, nfold)
  if (is.null(k_cv)) {
    return(rand_high_cost)
  }
  na_thresh <- nrow(k_cv) * na_thresh_prop
  if (sum(is.na(k_cv$residual)) >= na_thresh) {
    warning("Too many NAs in krige.cv: ", sum(is.na(k_cv$residual)), "/", na_thresh)
    return(rand_high_cost)
  }
  vgrm_rmse <- sqrt(sqrt(attr(fv, "SSErr") / nrow(ev)))
  vgrm_sd <- sqrt(sd(ev$gamma, na.rm = TRUE))
  vgrm_cost <- vgrm_weight * vgrm_rmse / vgrm_sd
  interp_rmse <- sqrt(sum(k_cv$residual^2, na.rm = TRUE) / nrow(k_cv))
  interp_sd <- sd(k_cv$var1.pred, na.rm = TRUE)
  interp_cost <- interp_weight * interp_rmse / interp_sd
  cost_value <- vgrm_cost + interp_cost
  if (vgrm_sd == 0 || interp_sd == 0) {
    warning("Vgrm or interpolation std is zero")
    return(rand_high_cost)
  }
  log_optimization(
    paste0("assets/nlopt_data/iterations/", id, "-", v_mod),
    c(
      "ID" = id, "N obs" = nrow(shp_hf), "Cutoff" = cutoff,
      "N lags" = n_lags, "N pairs" = n_max, "Max dist" = max_dist,
      "Vgrm weight" = vgrm_weight, "Vgrm rmse" = vgrm_rmse, "Vgrm sd" = vgrm_sd,
      "Vgrm cost" = vgrm_cost, "Int weight" = interp_weight,
      "Int rmse" = interp_rmse, "Int sd" = interp_sd, "Int cost" = interp_cost,
      "Cost" = cost_value, "Vgrm model" = v_mod
    )
  )
  return(cost_value)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# decode opt !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
decode_opt <- function(shp_hf = NULL, v_mod = NULL, opt = NULL) {
  suppressWarnings({
    if (is.null(v_mod)) {
      stop("\nMissing variogram model!")
    }
    if (is.null(shp_hf)) {
      stop("\nMissing heat flow data!")
    }
    if (is.null(opt)) {
      stop("\nMissing nloptr object!")
    }
    ev <- make_experimental_vgrm(shp_hf, opt$solution[1], opt$solution[2])
    if (is.null(ev)) {
      return(NULL)
    }
    fv <- fit_variogram_model(ev, v_mod)
    if (is.null(fv)) {
      return(NULL)
    }
  })
  return(list("experimental_vgrm" = as_tibble(ev), "fitted_vgrm" = fv))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# nlopt krige !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nlopt_krige <- function(id = NULL, type = "trans", v_mod = "Sph", alg = "NLOPT_LN_COBYLA",
                        max_eval = 500, iwt = 0.5, vwt = 0.5) {
  if (is.null(id)) {
    stop("\nMissing id!")
  }
  nlopt_dir <- "assets/nlopt_data"
  model_dir <- paste0(nlopt_dir, "/krige_models")
  iter_dir <- paste0(nlopt_dir, "/iterations")
  nlopt_id <- paste0(id, "-", v_mod)
  nlopt_path <- paste0(model_dir, "/", nlopt_id, ".RData")
  nlopt_itr_path <- paste0(iter_dir, "/", id, "-", v_mod)
  if (file.exists(nlopt_path) & file.exists(nlopt_itr_path)) {
    cat("\nOptimized", v_mod, "kriging model found for", type, ":", id)
    return(invisible())
  }
  if ((!file.exists(nlopt_path) & file.exists(nlopt_itr_path)) |
    (file.exists(nlopt_path) & !file.exists(nlopt_itr_path))) {
    cat("\n!!Optimized", v_mod, "kriging model failed for", type, ":", id, "!!")
    cat("\n     Remove", nlopt_path, "and ...")
    cat("\n     Remove", nlopt_itr_path, "and ...")
    cat('\n     Run "make nlopt" again ...')
    return(invisible())
  }
  if (!dir.exists(model_dir)) {
    dir.create(model_dir, recursive = T, showWarnings = F)
  }
  if (type == "rect") {
    x <- compile_rect_data(id)
    obs <- x$ghf_rect[[1]]
    na_thresh_prop <- 1
    x0 <- c(3, 50, 5, 0.5)
    lb <- c(1, 30, 2, 0.25)
    ub <- c(12, 100, 30, 0.75)
    xtol <- 1e-8
    ftol <- 1e-8
  } else if (type == "trans") {
    x <- compile_trans_data(id)
    obs <- x$ghf_large_buff[[1]]
    na_thresh_prop <- 0.5
    x0 <- c(3, 50, 5, 1e5)
    lb <- c(1, 30, 2, 5e4)
    ub <- c(12, 100, 30, 5e5)
    xtol <- 1e-8
    ftol <- 1e-8
  }
  opts <- list(print_level = 0, maxeval = max_eval, algorithm = alg, xtol_rel = xtol, ftol_rel = ftol)
  folds <- sample_and_shuffle_equal_nfolds(nrow(obs))
  nlopt_fun <- function(x) {
    cost_function(obs, x[1], x[2], x[3], x[4], v_mod, iwt, vwt, id, folds, na_thresh_prop)
  }
  tryCatch(
    {
      cat("\n", rep("-", 60), sep = "")
      cat("\nOptimizing", v_mod, "kriging model for:", id)
      cat("\n", rep("+", 60), sep = "")
      cat("\nNLopt algorithm:    ", alg)
      cat("\nMaximum evaluations:", max_eval)
      cat("\nKrige weight:       ", iwt)
      cat("\nVariogram weight:   ", vwt)
      cat("\n                    (cutoff, n_lags, n_max, max_dist)")
      cat("\nInitial parameters: ", x0)
      cat("\nLower bounds:       ", lb)
      cat("\nUpper bounds:       ", ub)
      opt <- nloptr(x0, nlopt_fun, lb = lb, ub = ub, opts = opts)
    },
    error = function(e) {
      cat("\n!! ERROR occurred in nlopt_krige:\n!!", conditionMessage(e))
    }
  )
  if (opt$status < 0 | opt$status == 5 | opt$iterations < 10) {
    cat("\n", rep("+", 60), sep = "")
    cat("\nNLopt failed to converge:")
    cat("\n   NLopt status:", opt$status)
    cat("\n   NLopt iterations:", opt$iterations)
    cat("\n   ", opt$message)
    cat("\n", rep("-", 60), sep = "")
    return(invisible())
  } else {
    cat("\n", rep("+", 60), sep = "")
    cat("\nNLopt converged:")
    cat("\n   NLopt status:", opt$status)
    cat("\n   NLopt iterations:", opt$iterations)
    cat("\n   ", opt$message)
    cat("\n", rep("-", 60), sep = "")
    opt_decoded <- decode_opt(obs, v_mod, opt)
    assign(str_replace_all(nlopt_id, "-", "_"), opt_decoded)
    save(list = str_replace_all(nlopt_id, "-", "_"), file = nlopt_path)
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# nlopt rects !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nlopt_rects <- function(rect_ids = NULL, v_mods = c("Sph", "Exp", "Lin"), alg = "NLOPT_LN_COBYLA",
                        max_eval = 500, iwt = 0.5, vwt = 0.5, parallel = T) {
  if (is.null(rect_ids)) {
    stop("\nMissing global rect ids!")
  }
  if (length(list.files("assets/map_data/relief")) < length(rect_ids)) {
    parallel <- F
  }
  x <- expand.grid(id = rect_ids, vm = v_mods, stringsAsFactors = F) %>% arrange(id, vm)
  cat("\nOptimizing krige models after Li et al. (2018) ...\n")
  if (parallel) {
    plan(multisession, workers = availableCores() - 2)
    future_walk2(x$id, x$vm, ~ nlopt_krige_rects(..., alg, max_eval, iwt, vwt),
      .options = furrr_options(seed = seed), .progress = T
    )
  } else {
    walk2(x$id, x$vm, ~ nlopt_krige_rects(..., alg, max_eval, iwt, vwt))
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# nlopt transects !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nlopt_transects <- function(trans_ids = NULL, v_mods = c("Sph", "Exp", "Lin"),
                            alg = "NLOPT_LN_COBYLA", max_eval = 500, iwt = 0.5, vwt = 0.5,
                            parallel = T) {
  if (is.null(trans_ids)) {
    stop("\nMissing submap transect ids!")
  }
  if (length(list.files("assets/map_data/relief")) < length(trans_ids)) {
    parallel <- F
  }
  x <- expand.grid(id = trans_ids, vm = v_mods, stringsAsFactors = F) %>% arrange(id, vm)
  cat("\nOptimizing krige models after Li et al. (2018) ...\n")
  if (parallel) {
    plan(multisession, workers = availableCores() - 2)
    future_walk2(x$id, x$vm, ~ nlopt_krige_trans(..., alg, max_eval, iwt, vwt),
      .options = furrr_options(seed = seed), .progress = T
    )
  } else {
    walk2(x$id, x$vm, ~ nlopt_krige_trans(..., alg, max_eval, iwt, vwt))
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read nloptr itr !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
read_nloptr_itr <- function(fpath) {
  read_file <- function(path, ...) {
    con <- file(path)
    on.exit(close(con))
    suppressWarnings(readLines(con, ...))
  }
  extract_value <- function(lines, pattern, as_numeric = TRUE) {
    matched_lines <- lines[grepl(pattern, lines)]
    clean_values <- map_chr(matched_lines, ~ gsub(pattern, "", .x))
    if (as_numeric) {
      return(map_dbl(clean_values, ~ as.numeric(.x)))
    } else {
      return(clean_values)
    }
  }
  t <- read_file(fpath)
  id_t <- extract_value(t, "^ID: ", FALSE)
  vmod_t <- extract_value(t, "^Vgrm model: ", FALSE)
  n_obs_t <- extract_value(t, "^N obs: ")
  cutoff_t <- extract_value(t, "^Cutoff: ")
  n_lags_t <- extract_value(t, "^N lags: ")
  n_pairs_t <- extract_value(t, "^N pairs: ")
  max_d_t <- extract_value(t, "^Max dist: ")
  vwt_t <- extract_value(t, "^Vgrm weight: ")
  vrmse_t <- extract_value(t, "^Vgrm rmse: ")
  vcost_t <- extract_value(t, "^Vgrm cost: ")
  iwt_t <- extract_value(t, "^Int weight: ")
  irmse_t <- extract_value(t, "^Int rmse: ")
  icost_t <- extract_value(t, "^Int cost: ")
  cost_t <- extract_value(t, "^Cost: ")
  tibble(
    short_name = id_t, n_obs = n_obs_t, cutoff = cutoff_t, n_lags = n_lags_t,
    n_pairs = n_pairs_t, max_dist = max_d_t, v_mod = vmod_t, vgrm_wt = vwt_t, vgrm_rmse = vrmse_t,
    vgrm_cost = vcost_t, cv_wt = iwt_t, cv_rmse = irmse_t, cv_cost = icost_t, cost = cost_t
  ) %>%
    group_by(short_name, v_mod) %>%
    mutate(itr = row_number(), .after = short_name) %>%
    ungroup()
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get optimal krige model !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_optimal_krige_model <- function(id = NULL) {
  if (is.null(id)) {
    stop("\nMissing id!")
  }
  nlopt_dir <- "assets/nlopt_data/"
  itr_dir <- paste0(nlopt_dir, "iterations")
  model_dir <- paste0(nlopt_dir, "krige_models")
  itr_paths <- list.files(itr_dir, pattern = id, full.names = T)
  model_paths <- list.files(model_dir, pattern = id, full.names = T)
  if (length(itr_paths) < 1) {
    cat("\n   No nlopt itr files found for:", id)
    return(NULL)
  }
  if (length(model_paths) < 1) {
    cat("\n   No nlopt itr files found for:", id)
    return(NULL)
  }
  if (length(itr_paths) != length(model_paths)) {
    itr_mods <- str_sub(itr_paths, start = -3)
    krige_mods <- str_extract(model_paths, ".{3}(?=\\.RData)")
    if (length(krige_mods) < 1) {
      cat("\n   No nlopt RData found for:", id)
      return(NULL)
    }
    itr_paths <- itr_paths[str_detect(itr_paths, paste(krige_mods, collapse = "|"))]
  }
  nlopt_itr <- map_df(itr_paths, read_nloptr_itr)
  k_mod <- slice_min(nlopt_itr, cost)
  if (nrow(k_mod) > 1) {
    k_mod <- slice(k_mod, nrow(k_mod))
  }
  nlopt_path <- paste0(model_dir, "/", id, "-", k_mod$v_mod, ".RData")
  if (file.exists(nlopt_path)) {
    load(nlopt_path)
    opt_decoded <- get(paste0(id, "_", k_mod$v_mod))
    list(
      "nlopt_itr" = nlopt_itr, "opt_krige_mod_summary" = k_mod,
      "experimental_vgrm" = opt_decoded$experimental_vgrm,
      "fitted_vgrm" = opt_decoded$fitted_vgrm
    )
  } else {
    cat("\n   No nlopt RData found for:", id)
    return(NULL)
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Krige !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Krige <- function(shp_hf = NULL, fv = NULL, shp_grid = NULL, n_max = 10, max_dist = 1e5) {
  if (is.null(shp_hf)) {
    stop("\nMissing heat flow data!")
  }
  if (is.null(fv)) {
    stop("\nMissing variogram model!")
  }
  if (is.null(shp_grid)) {
    stop("\nMissing kriging locations (grid)!")
  }
  if (is.null(n_max)) {
    stop("\nMissing max pairs!")
  }
  if (is.null(max_dist)) {
    stop("\nMissing max distance!")
  }
  suppressWarnings({
    tryCatch(
      {
        krige(obs ~ 1, shp_hf,
          newdata = shp_grid, model = fv, nmax = n_max, maxdist = max_dist,
          debug.level = 0
        ) %>%
          as_tibble() %>%
          st_as_sf() %>%
          rename(est_krg = var1.pred, var_krg = var1.var, krige = geometry) %>%
          mutate(sigma_krg = sqrt(var_krg), .before = krige) %>%
          mutate(
            est_krg = ifelse(est_krg > 0 & est_krg <= 250, est_krg, NA),
            var_krg = ifelse(est_krg > 0 & est_krg <= 250, est_krg, NA)
          )
      },
      error = function(e) {
        cat("\n!! ERROR occurred in Krige:\n!!", conditionMessage(e))
      }
    )
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# interp diff !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
interp_diff <- function(shp_krg = NULL, shp_sim = NULL) {
  if (is.null(shp_krg)) {
    stop("\nMissing krige interpolation!")
  }
  if (is.null(shp_sim)) {
    stop("\nMissing similarity interpolation!")
  }
  if (any(names(shp_sim) == "projected_distances")) {
    shp_sim <-
      shp_sim %>%
      rename(est_sim = obs, sigma_sim = sigma) %>%
      select(est_sim, sigma_sim, geometry)
    shp_krg <-
      shp_krg %>%
      rename(est_krg = obs, sigma_krg = sigma) %>%
      select(est_krg, sigma_krg, geometry)
    shp_sim %>%
      mutate(est_krg = shp_krg$est_krg, sigma_krg = shp_krg$sigma_krg) %>%
      mutate(est_dff = est_sim - est_krg, sigma_dff = sigma_sim - sigma_krg, .before = geometry) %>%
      rename(projected_diff = geometry) %>%
      select(est_dff, sigma_dff, projected_diff)
  } else {
    shp_sim %>%
      mutate(est_krg = shp_krg$est_krg, sigma_krg = shp_krg$sigma_krg) %>%
      mutate(est_dff = est_sim - est_krg, sigma_dff = sigma_sim - sigma_krg, .before = similarity) %>%
      rename(diff = similarity) %>%
      select(est_dff, sigma_dff, diff)
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# summarize optimal krige models !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summarize_optimal_krige_models <- function(trans_ids, parallel = T) {
  f <- function(x) {
    get_optimal_krige_model(x)$opt_krige_mod_summary
  }
  cat("\nSummarizing optimal krige models ...\n")
  if (parallel) {
    plan(multisession, workers = availableCores() - 2)
    as_tibble(future_map_dfr(trans_ids, f, .options = furrr_options(seed = seed), .progress = T))
  } else {
    map_df(trans_ids, f)
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# summarize interpolation differences !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summarize_interp_differences <- function(trans_ids, parallel = T) {
  f <- function(x) {
    tryCatch(
      {
        km <- get_optimal_krige_model(x)
        fv <- km$fitted_vgrm
        np <- km$opt_krige_mod_summary$n_pairs
        compile_trans_data(x, fv = fv, np = np)$dff_large_buff[[1]] %>%
          st_set_geometry(NULL) %>%
          summarise(
            short_name = x, n_grid = n(), n = sum(!is.na(est_dff)),
            rmse = sqrt(mean(est_dff^2, na.rm = T)), min = min(est_dff, na.rm = T),
            max = max(est_dff, na.rm = T), med = median(est_dff, na.rm = T),
            iqr = IQR(est_dff, na.rm = T), mean = mean(est_dff, na.rm = T),
            sigma = sd(est_dff, na.rm = T)
          )
      },
      error = function(e) {
        cat("\n!! ERROR occurred in summarize_interp_differences:\n!!", conditionMessage(e))
        return(NULL)
      }
    )
  }
  cat("\nSummarizing interpolation differences ...\n")
  if (parallel) {
    plan(multisession, workers = availableCores() - 2)
    as_tibble(future_map_dfr(trans_ids, f, .options = furrr_options(seed = seed), .progress = T))
  } else {
    map_df(trans_ids, f)
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# summarize interpolation accuracy !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summarize_interp_accuracy <- function(trans_ids, parallel = T) {
  f <- function(x) {
    suppressWarnings({
      suppressMessages({
        tryCatch(
          {
            km <- get_optimal_krige_model(x)
            fv <- km$fitted_vgrm
            np <- km$opt_krige_mod_summary$n_pairs
            comps_sim <- get_closest_interp_obs(x)
            comps_krg <- get_closest_interp_obs(x, fv = fv, np = np)
            tibble(
              short_name = x,
              n_sim = nrow(comps_sim[!is.na(comps_sim$est_sim), ]),
              min_obs_sim = min(comps_sim$est_sim - comps_sim$obs_ghf, na.rm = T),
              max_obs_sim = max(comps_sim$est_sim - comps_sim$obs_ghf, na.rm = T),
              mean_obs_sim = mean(comps_sim$est_sim - comps_sim$obs_ghf, na.rm = T),
              sd_obs_sim = sd(comps_sim$est_sim - comps_sim$obs_ghf, na.rm = T),
              med_obs_sim = median(comps_sim$est_sim - comps_sim$obs_ghf, na.rm = T),
              iqr_obs_sim = IQR(comps_sim$est_sim - comps_sim$obs_ghf, na.rm = T),
              rmse_obs_sim = sqrt(mean((comps_sim$est_sim - comps_sim$obs_ghf)^2, na.rm = T)),
              n_krg = nrow(comps_krg[!is.na(comps_krg$est_krg), ]),
              min_obs_krg = min(comps_krg$est_krg - comps_krg$obs_ghf, na.rm = T),
              max_obs_krg = max(comps_krg$est_krg - comps_krg$obs_ghf, na.rm = T),
              mean_obs_krg = mean(comps_krg$est_krg - comps_krg$obs_ghf, na.rm = T),
              sd_obs_krg = sd(comps_krg$est_krg - comps_krg$obs_ghf, na.rm = T),
              med_obs_krg = median(comps_krg$est_krg - comps_krg$obs_ghf, na.rm = T),
              iqr_obs_krg = IQR(comps_krg$est_krg - comps_krg$obs_ghf, na.rm = T),
              rmse_obs_krg = sqrt(mean((comps_krg$est_krg - comps_krg$obs_ghf)^2, na.rm = T))
            )
          },
          error = function(e) {
            cat("\n!! ERROR occurred in summarize_interp_accuracy:\n!!", conditionMessage(e))
            return(NULL)
          }
        )
      })
    })
  }
  cat("\nSummarizing interpolation accuracies ...\n")
  if (parallel) {
    as_tibble(future_map_dfr(trans_ids, f, .options = furrr_options(seed = seed), .progress = T))
  } else {
    map_df(trans_ids, f)
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write ghf raw !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_ghf_raw <- function() {
  save_dir <- "assets/hf_data"
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = T, showWarnings = F)
  }
  tryCatch(
    {
      fpath <- paste0(save_dir, "/ghf-raw.csv")
      if (file.exists(fpath)) {
        cat("\nRaw HF observations csv found!")
        return(invisible())
      }
      cat("\n   Processing and writing raw IHFC_2024GHFDB.xlsx ...\n")
      read_excel("assets/hf_data/IHFC_2024_GHFDB.xlsx", skip = 5) %>%
        select(long_EW, lat_NS, q, elevation, Quality_Code, ID) %>%
        rename(lon = long_EW, lat = lat_NS, obs = q, elev = elevation, quality = Quality_Code, id = ID) %>%
        mutate(across(where(is.numeric), ~ formatC(.x, format = "f", digits = 3))) %>%
        write_csv(fpath)
    },
    error = function(e) {
      cat("\n!! ERROR occurred in write_ghf_raw:\n!!", conditionMessage(e))
      return(NULL)
    }
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write ghf transects !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_ghf_transects <- function(trans_ids) {
  save_dir <- "assets/hf_data/ghf_observations_submap_transects"
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = T, showWarnings = F)
  }
  f <- function(x) {
    tryCatch(
      {
        fpath <- paste0(save_dir, "/", x, "-ghf.csv")
        if (file.exists(fpath)) {
          cat("\nKriging predictions csv found for submap transect:", x)
          return(invisible())
        }
        cat("\n   Processing:", x)
        trans_df <- compile_trans_data(x)
        trans_df$ghf_large_buff[[1]] %>%
          reproject_wgs() %>%
          mutate(lon = st_coordinates(.)[, 1], lat = st_coordinates(.)[, 2], .before = obs) %>%
          mutate(across(where(is.numeric), ~ formatC(.x, format = "f", digits = 3))) %>%
          st_set_geometry(NULL) %>%
          mutate(
            transect_id = trans_df$id, transect_zone = trans_df$zone,
            transect_short_name = trans_df$short_name, .before = obs
          ) %>%
          write_csv(fpath)
      },
      error = function(e) {
        cat("\n!! ERROR occurred in write_ghf_transects:\n!!", conditionMessage(e))
        return(NULL)
      }
    )
  }
  cat("\nWriting heat flow observations ...")
  walk(trans_ids, f)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write sim raw !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_sim_raw <- function() {
  save_dir <- "assets/hf_data"
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = T, showWarnings = F)
  }
  tryCatch(
    {
      fpath <- paste0(save_dir, "/sim-raw.csv")
      if (file.exists(fpath)) {
        cat("\nRaw Similarity csv found!")
        return(invisible())
      }
      cat("\n   Processing and writing raw HFgrid14.csv ...\n")
      read_delim("assets/hf_data/HFgrid14.csv", delim = ";", col_types = c("ddddd")) %>%
        rename(
          lon = longitude, lat = latitude, est_sim = HF_pred, sigma_sim = sHF_pred,
          obs_sim = Hf_obs
        ) %>%
        mutate(across(where(is.numeric), ~ formatC(.x, format = "f", digits = 3))) %>%
        write_csv(fpath)
    },
    error = function(e) {
      cat("\n!! ERROR occurred in write_sim_raw:\n!!", conditionMessage(e))
      return(NULL)
    }
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write sim predictions transects !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_sim_predictions_transects <- function(trans_ids) {
  save_dir <- "assets/hf_data/similarity_predictions_submap_transects"
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = T, showWarnings = F)
  }
  f <- function(x) {
    tryCatch(
      {
        fpath <- paste0(save_dir, "/", x, "-sim.csv")
        if (file.exists(fpath)) {
          cat("\nKriging predictions csv found for submap transect:", x)
          return(invisible())
        }
        cat("\n   Processing:", x)
        trans_df <- compile_trans_data(x)
        trans_df$sim_large_buff[[1]] %>%
          reproject_wgs() %>%
          mutate(lon = st_coordinates(.)[, 1], lat = st_coordinates(.)[, 2], .before = est_sim) %>%
          mutate(across(where(is.numeric), ~ formatC(.x, format = "f", digits = 3))) %>%
          st_set_geometry(NULL) %>%
          mutate(
            transect_id = trans_df$id, transect_zone = trans_df$zone,
            transect_short_name = trans_df$short_name, .before = est_sim
          ) %>%
          write_csv(fpath)
      },
      error = function(e) {
        cat("\n!! ERROR occurred in write_sim_predictions_transects:\n!!", conditionMessage(e))
        return(NULL)
      }
    )
  }
  cat("\nWriting Similarity predictions ...")
  walk(trans_ids, f)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write krg predictions transects !!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_krg_predictions_transects <- function(trans_ids) {
  save_dir <- "assets/hf_data/kriging_predictions_submap_transects"
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = T, showWarnings = F)
  }
  f <- function(x) {
    tryCatch(
      {
        fpath <- paste0(save_dir, "/", x, "-krg.csv")
        if (file.exists(fpath)) {
          cat("\nKriging predictions csv found for submap transect:", x)
          return(invisible())
        }
        cat("\n   Processing:", x)
        km <- get_optimal_krige_model(x)
        fv <- km$fitted_vgrm
        np <- km$opt_krige_mod_summary$n_pairs
        trans_df <- compile_trans_data(x, fv = fv, np = np)
        trans_df$krg_large_buff[[1]] %>%
          reproject_wgs() %>%
          mutate(lon = st_coordinates(.)[, 1], lat = st_coordinates(.)[, 2], .before = est_krg) %>%
          mutate(across(where(is.numeric), ~ formatC(.x, format = "f", digits = 3))) %>%
          st_set_geometry(NULL) %>%
          mutate(
            transect_id = trans_df$id, transect_zone = trans_df$zone,
            transect_short_name = trans_df$short_name, .before = est_krg
          ) %>%
          write_csv(fpath)
      },
      error = function(e) {
        cat("\n!! ERROR occurred in write_krg_predictions_transects:\n!!", conditionMessage(e))
        return(NULL)
      }
    )
  }
  cat("\nWriting Krige predictions ...")
  walk(trans_ids, f)
}

