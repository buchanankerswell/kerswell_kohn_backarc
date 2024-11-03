#######################################################
## Process spatial datasets                          ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
normalize_lon <- function(lon, to360 = FALSE) {
  if (!to360) {
    ((lon + 180) %% 360) - 180
  } else {
    ((lon %% 360) + 360) %% 360
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
modular_lon_span <- function(xmin, xmax) {
  lon1_360 <- normalize_lon(xmin, to360 = TRUE)
  lon2_360 <- normalize_lon(xmax, to360 = TRUE)
  delta <- (lon2_360 - lon1_360) %% 360
  if (delta > 180) 360 - delta else delta
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
crosses_antimeridian <- function(xmin, xmax) {
  lon1 <- normalize_lon(xmin, to360 = TRUE)
  lon2 <- normalize_lon(xmax, to360 = TRUE)
  delta <- modular_lon_span(lon1, lon2)
  (lon1 > lon2) && (delta < 180 && delta > 0)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
modular_mid_lon <- function(xmin, xmax) {
  lon1_360 <- normalize_lon(xmin, to360 = TRUE)
  lon2_360 <- normalize_lon(xmax, to360 = TRUE)
  delta <- (lon2_360 - lon1_360) %% 360
  if (delta <= 180) {
    mid_360 <- (lon1_360 + delta / 2) %% 360
  } else {
    delta2 <- 360 - delta
    mid_360 <- (lon2_360 + delta2 / 2) %% 360
  }
  normalize_lon(mid_360, to360 = FALSE)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
determine_projection <- function(sf_obj, proj_type = c("ortho", "aeqd", "albers", "longlat")) {
  with_error_handling({
    bbox <- st_bbox(sf_obj)

    lon0 <- modular_mid_lon(bbox["xmin"], bbox["xmax"])
    lat0 <- (bbox["ymin"] + bbox["ymax"]) / 2
    proj_type <- match.arg(proj_type)

    if (proj_type == "ortho") {
      proj_string <- paste0("+proj=ortho +lon_0=", lon0, " +lat_0=", lat0, " +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs")
    } else if (proj_type == "aeqd") {
      proj_string <- paste0("+proj=aeqd +lon_0=", lon0, " +lat_0=", lat0, " +datum=WGS84 +units=m +no_defs")
    } else if (proj_type == "albers") {
      lat1 <- bbox["ymin"]
      lat2 <- bbox["ymax"]
      if (lat1 == lat2) {
        lat1 <- lat1 - 1
        lat2 <- lat2 + 1
      }
      proj_string <- paste0("+proj=aea +lat_1=", lat1, " +lat_2=", lat2, " +lat_0=", lat0, " +lon_0=", lon0, " +datum=WGS84 +units=m +no_defs")
    } else if (proj_type == "longlat") {
      proj_string <- "+proj=longlat +datum=WGS84 +no_defs"
    } else {
      proj_string <- "+proj=longlat +datum=WGS84 +no_defs"
    }

    list(wkt = st_crs(proj_string), lon0 = lon0, lat0 = lat0)
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
crop_to_bbox <- function(sf_obj, sf_bbox, crop_within_region = TRUE, keep_df = TRUE) {
  with_error_handling({
    if (st_crs(sf_obj) != st_crs(sf_bbox)) sf_obj <- st_transform(sf_obj, st_crs(sf_bbox))

    sf_obj <- st_make_valid(sf_obj)
    sf_bbox <- st_make_valid(sf_bbox)

    sf_obj_cropped <- if (!crop_within_region) st_crop(sf_obj, sf_bbox) else st_intersection(sf_obj, sf_bbox)
    sf_obj_length <- if (!is.data.frame(sf_obj_cropped)) length(sf_obj_cropped) else nrow(sf_obj_cropped)

    if (sf_obj_length == 0) {
      invisible()
    } else {
      if (!keep_df) st_sfc(st_union(sf_obj_cropped), crs = st_crs(sf_bbox)) else sf_obj_cropped
    }
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
parse_zerodist_observations <- function(df, data_col, qc_col) {
  with_error_handling({
    duplicate_observations <- zerodist(as_Spatial(df))
    n_duplicate_observaitons <- nrow(duplicate_observations)

    locations_to_remove <- map_dbl(1:n_duplicate_observaitons, ~ {
      location1 <- duplicate_observations[.x, 1]
      location2 <- duplicate_observations[.x, 2]

      location_random <- sample(c(location1, location2), 1)

      obs1 <- df[[data_col]][location1]
      obs2 <- df[[data_col]][location2]

      if (obs1 == obs2) {
        return(location_random)
      }

      u1_str <- str_extract(df[[qc_col]][location1], "(?<=U).")
      u1 <- if (u1_str == "x") {
        5
      } else {
        as.numeric(u1_str)
      }

      u2_str <- str_extract(df[[qc_col]][location2], "(?<=U).")
      u2 <- if (u2_str == "x") {
        5
      } else {
        as.numeric(u2_str)
      }

      m1_str <- str_extract(df[[qc_col]][location1], "(?<=M).")
      m1 <- if (m1_str == "x") {
        5
      } else {
        as.numeric(m1_str)
      }

      m2_str <- str_extract(df[[qc_col]][location2], "(?<=M).")
      m2 <- if (m2_str == "x") {
        5
      } else {
        as.numeric(m2_str)
      }

      total_quality_location1 <- sum(u1 + m1, na.rm = TRUE)
      total_quality_location2 <- sum(u2 + m2, na.rm = TRUE)

      if (total_quality_location1 < total_quality_location2) {
        return(location2)
      }
      if (total_quality_location2 < total_quality_location1) {
        return(location1)
      }
      if (total_quality_location2 == total_quality_location1) {
        return(location_random)
      }
    })

    slice(df, -locations_to_remove)
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
build_spatraster_from_matrix <- function(mat, lons, lats) {
  with_error_handling({
    lons <- as.numeric(lons)
    lats <- as.numeric(lats)

    ord_lon <- order(lons)
    ord_lat <- order(lats, decreasing = TRUE)

    mat <- mat[ord_lon, ord_lat, drop = FALSE]
    lons <- lons[ord_lon]
    lats <- lats[ord_lat]

    r <- rast(t(mat), extent = c(min(lons), max(lons), min(lats), max(lats)), crs = "EPSG:4326")

    crosses <- crosses_antimeridian(min(lons), max(lons))
    if (crosses) {
      lons_360 <- normalize_lon(lons, to360 = crosses)
      ord_lon <- order(lons_360)
      mat <- mat[ord_lon, , drop = FALSE]
      lons_360 <- lons_360[ord_lon]
      r <- rast(t(mat), extent = c(min(lons_360), max(lons_360), min(lats), max(lats)), crs = "EPSG:4326")
    }

    return(r)
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
to_spatraster <- function(obj, data_col = NULL) {
  with_error_handling({
    if (inherits(obj, "Raster") || inherits(obj, "SpatRaster")) {
      return(rast(obj))
    }

    if (inherits(obj, "stars")) {
      return(rast(obj))
    }

    if (inherits(obj, "bathy")) {
      bathy_mat <- as.matrix(obj)
      lons <- as.numeric(rownames(bathy_mat))
      lats <- as.numeric(colnames(bathy_mat))

      r <- build_spatraster_from_matrix(bathy_mat, lons, lats)
      return(r)
    }

    if (inherits(obj, "sf")) {
      if (is.null(data_col)) stop("For sf input, provide 'data_col' to rasterize")

      coords <- st_coordinates(obj)

      coord_lon <- round(coords[, 1], digits = 3)
      coord_lat <- round(coords[, 2], digits = 3)
      lons <- sort(unique(coord_lon))
      lats <- sort(unique(coord_lat))

      sf_mat <- matrix(NA_real_, nrow = length(lons), ncol = length(lats), dimnames = list(lons, lats))

      lon_idx <- match(coord_lon, lons)
      lat_idx <- match(coord_lat, lats)

      sf_mat[cbind(lon_idx, lat_idx)] <- obj[[data_col]]

      r <- build_spatraster_from_matrix(sf_mat, lons, lats)
      return(r)
    }

    stop("Unsupported input type: must be Raster, SpatRaster, stars, bathy, or sf")
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_noaa_relief <- function(lon1, lon2, lat1, lat2, data_dir, out_dir, arc_seconds = 30, res = NULL, add_geoid = FALSE, sea_level_offset = NULL) {
  with_error_handling({
    if (is.null(res)) res_tag <- "native-res" else paste0("res-", res)
    if (is.null(sea_level_offset)) sea_level_tag <- "0m-offset" else paste0(sea_level_offset, "m-offset")
    fname <- sprintf("relief-wgs84-lon(%.2f,%.2f)-lat(%.2f,%.2f)-%.0f-arc-seconds-%s-%s.tif", lon1, lon2, lat1, lat2, arc_seconds, res_tag, sea_level_tag)
    cache_file <- file.path(out_dir, fname)
    if (file.exists(cache_file)) {
      return(rast(cache_file))
    }

    if (arc_seconds == 30) {
      base_relief_url <- "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/data/30s/30s_surface_elev_gtif/ETOPO_2022_v1_30s_N90W180_surface.tif"
      base_geoid_url <- "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/data/30s/30s_geoid_gtif/ETOPO_2022_v1_30s_N90W180_geoid.tif"
    } else if (arc_seconds == 60) {
      base_relief_url <- "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/data/60s/60s_surface_elev_gtif/ETOPO_2022_v1_60s_N90W180_surface.tif"
      base_geoid_url <- "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/data/60s/60s_geoid_gtif/ETOPO_2022_v1_60s_N90W180_geoid.tif"
    } else {
      base_relief_url <- "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/data/60s/60s_surface_elev_gtif/ETOPO_2022_v1_60s_N90W180_surface.tif"
      base_geoid_url <- "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/data/60s/60s_geoid_gtif/ETOPO_2022_v1_60s_N90W180_geoid.tif"
    }

    base_relief_file <- file.path(data_dir, basename(base_relief_url))
    if (!file.exists(base_relief_file)) {
      message(paste0("Downloading ETOPO2022 global relief model (", arc_seconds, " arc-second)..."))
      old_timeout <- getOption("timeout")
      options(timeout = max(2400, old_timeout))
      download.file(base_relief_url, base_relief_file, mode = "wb", quiet = FALSE)
      options(timeout = old_timeout)
    }

    base_geoid_file <- file.path(data_dir, basename(base_geoid_url))
    if (!file.exists(base_geoid_file)) {
      message(paste0("Downloading ETOPO2022 global geoid model (", arc_seconds, " arc-second)..."))
      old_timeout <- getOption("timeout")
      options(timeout = max(2400, old_timeout))
      download.file(base_geoid_url, base_geoid_file, mode = "wb", quiet = FALSE)
      options(timeout = old_timeout)
    }

    r_geoid <- rast(base_geoid_file)
    r_relief <- rast(base_relief_file)
    pixel_width <- xres(r_relief)
    ymin <- min(lat1, lat2)
    ymax <- max(lat1, lat2)
    crosses <- crosses_antimeridian(lon1, lon2)

    if (!crosses) {
      r_relief <- crop(r_relief, ext(min(lon1, lon2), max(lon1, lon2), ymin, ymax))
      r_geoid <- crop(r_geoid, ext(min(lon1, lon2), max(lon1, lon2), ymin, ymax))
    } else {
      r_relief <- rotate(r_relief)
      r_geoid <- rotate(r_geoid)

      lon1_360 <- normalize_lon(lon1, to360 = crosses)
      lon2_360 <- normalize_lon(lon2, to360 = crosses)

      rotated_extent <- ext(min(lon1_360, lon2_360), max(lon1_360, lon2_360), ymin, ymax)
      r_relief <- crop(r_relief, rotated_extent)
      r_geoid <- crop(r_geoid, rotated_extent)
    }

    if (add_geoid) r_relief <- r_relief + r_geoid
    if (!is.null(sea_level_offset) && is.numeric(sea_level_offset)) r_relief <- r_relief + sea_level_offset

    if (!is.null(res)) {
      target_deg <- res / 60
      fact <- round(target_deg / pixel_width)
      if (fact > 1) r_relief <- aggregate(r_relief, fact = fact, fun = "mean", na.rm = TRUE)
    }

    crs(r_relief) <- "EPSG:4326"

    create_dir(out_dir)
    writeRaster(r_relief, cache_file, overwrite = TRUE)

    return(r_relief)
  })
}
