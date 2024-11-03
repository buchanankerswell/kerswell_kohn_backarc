#######################################################
## Preprocess spatial datasets                       ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
download_data_from_osf <- function(data_dir) {
  with_error_handling({
    if (!file.exists(file.path(data_dir, "depth-profile-data.csv"))) {
      if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
      osf_retrieve_node("ca6zu") |>
        osf_ls_files(pattern = "data") |>
        osf_download(path = data_dir, conflicts = "overwrite", recurse = TRUE)
    } else {
      cat(" -- Found data: ", data_dir, "\n", sep = "")
    }
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fetch_coastline <- function(scale = 10) {
  with_error_handling({
    cat(" -> Processing world map\n")
    cat("    Fetched from: rnaturalearth coastlines (", scale, "m)\n", sep = "")
    shp <- ne_coastline(scale = scale, returnclass = "sf") |> st_transform(4326)
    shp <- shp[st_length(shp) > set_units(1e5, "m"), ]
    st_simplify(shp, dTolerance = 0.05, preserveTopology = TRUE)
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fetch_plate_boundaries <- function(path) {
  with_error_handling({
    cat(" -> Processing UTIG plate boundary data\n")
    cat("    Fetched from: http://www-udc.ig.utexas.edu/external/plates/data.htm\n")
    cat("    Reading data: ", path, "\n", sep = "")
    st_read(path, crs = 4326, quiet = TRUE) |> filter(map_lgl(geometry, ~ nrow(st_coordinates(.x)) > 1))
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fetch_volcanoes <- function(dir) {
  with_error_handling({
    cat(" -> Processing Syracuse and Abers (2006) volcano data\n")
    cat("    Fetched from: https://doi.org/10.1029/2005GC001045\n")

    paths <- file.path(dir, c("volcanoes.txt", "volcanoes_nospread.txt"))
    walk(paths, ~ cat("    Reading data: ", .x, "\n", sep = ""))

    map_dfr(paths, ~ read_table(.x, col_type = "ddddddddddddccf", na = "NULL")) |>
      st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fetch_ihfc <- function(path, limits = NULL) {
  with_error_handling({
    cat(" -> Processing IHFC database 2024 release\n")
    cat("    Fetched from: https://doi.org/10.5880/fidgeo.2024.014\n")
    cat("    Reading data: ", path, "\n", sep = "")
    shp <-
      read_excel(path, skip = 5) |>
      select(long_EW, lat_NS, q, elevation, Quality_Code, ID) |>
      rename(lon = long_EW, lat = lat_NS, ihfc2024_obs = q, id = ID, quality = Quality_Code)

    if (!is.null(limits) && length(limits) == 2 && is.numeric(limits)) {
      cat("    Filtering data: heatflow [", limits[1], ", ", limits[2], "] mW/m^2\n", sep = "")
      shp <- shp |> filter(ihfc2024_obs > limits[1] & ihfc2024_obs <= limits[2])
    }

    shp |>
      st_as_sf(coords = c(1, 2), crs = 4326, remove = FALSE) |>
      parse_zerodist_observations("ihfc2024_obs", "quality") |>
      rename(ihfc2024_geometry = geometry)
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fetch_lucazeau <- function(path, limits = NULL) {
  with_error_handling({
    cat(" -> Processing similarity interpolation (Lucazeau, 2019)\n")
    cat("    Fetched from: https://doi.org/10.1029/2019GC008389\n")
    cat("    Reading data: ", path, "\n", sep = "")
    shp <-
      read_delim(path, delim = ";", col_types = c("ddddd")) |>
      rename(lon = longitude, lat = latitude, lucazeau2019_sim_est = HF_pred, lucazeau2019_sim_sigma = sHF_pred, lucazeau2019_obs = Hf_obs)

    if (!is.null(limits) && length(limits) == 2 && is.numeric(limits)) {
      cat("    Filtering data: heatflow [", limits[1], ", ", limits[2], "] mW/m^2\n", sep = "")
      shp <- shp |> filter(lucazeau2019_sim_est > 0 & lucazeau2019_sim_est <= 250)
    }

    shp |>
      st_as_sf(coords = c(1, 2), crs = 4326, remove = FALSE) |>
      rename(lucazeau2019_geometry = geometry)
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fetch_submap_transects <- function(dir) {
  with_error_handling({
    cat(" -> Processing Submap data (Lallemand & Heuret, 2017)\n")
    cat("    Fetched from: https://submap.gm.umontpellier.fr\n")
    cat("    Reading data: ", dir, "\n", sep = "")
    compile_jsons_to_df(list.files(dir, full.names = TRUE)) |>
      rename_all(tolower) |>
      filter(zone != "MED") |>
      filter(!(trench_name %in% c("Sandwich"))) |>
      mutate_all(~ ifelse(. == -999, NA, .)) |>
      mutate(short_name = sub("^(\\D+)(\\d)$", "\\10\\2", short_name)) |>
      mutate(`next` = sub("^(\\D+)(\\d)$", "\\10\\2", `next`)) |>
      mutate(previous = sub("^(\\D+)(\\d)$", "\\10\\2", previous)) |>
      arrange(zone, short_name) |>
      mutate(geometry = sprintf("LINESTRING(%s %s, %s %s)", lon1, lat1, lon2, lat2)) |>
      st_as_sf(wkt = "geometry", crs = 4326, remove = FALSE) |>
      rename(submap_transect_zone = zone, submap_transect_id = short_name, submap_transect_geometry = geometry) |>
      select(submap_transect_zone, submap_transect_id, trench_name, harc1, a, phi, m56_vc, m56_azvc, m56_vsn, alphad, dz, dz1, n, n1, tau, tau1)
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_submap_hulls <- function(sf_submap, buffer_size = 4e5, type = "concave", concavity = 6) {
  with_error_handling({
    cat(" -> Drawing ", type, " hulls around submap transect buffers\n", sep = "")
    sf_use_s2(FALSE)

    make_ids <- function(prefix, nums, width = 2) {
      paste0(prefix, str_pad(nums, width, pad = "0"))
    }

    match_list <- list(
      "NPA_SET1" = make_ids("NPA", 1:4),
      "NPA_SET2" = make_ids("NPA", 5:10),
      "NPA_SET3" = make_ids("NPA", 11:14),
      "NPA_SET4" = make_ids("NPA", 15:22),
      "NPA_SET5" = make_ids("NPA", 23:30),
      "NPA_SET6" = make_ids("NPA", 31:36),
      "NPA_SET7" = make_ids("NPA", 37:42),
      "NPA_SET8" = make_ids("NPA", 43:47),
      "SAM_SET1" = make_ids("SAM", 1:8),
      "SAM_SET2" = make_ids("SAM", 9:15),
      "SAM_SET3" = make_ids("SAM", 16:21),
      "SAM_SET4" = make_ids("SAM", 22:33),
      "SAM_SET5" = make_ids("SAM", 34:40),
      "SAM_SET6" = make_ids("SAM", 41:46),
      "SAM_SET7" = make_ids("SAM", 47:51),
      "SAM_SET8" = make_ids("SAM", 52:56),
      "SAM_SET9" = make_ids("SAM", 57:61),
      "SAM_SET10" = make_ids("SAM", 62:66),
      "SEA_SET1" = make_ids("SEA", 1:7),
      "SEA_SET2" = make_ids("SEA", 8:12),
      "SEA_SET3" = make_ids("SEA", 13:17),
      "SEA_SET4" = make_ids("SEA", c(18:22, 29:30)),
      "SEA_SET5" = make_ids("SEA", c(23:28, 31:34)),
      "SEA_SET6" = make_ids("SEA", 35:54),
      "SEA_SET7" = make_ids("SEA", 55:60),
      "SEA_SET8" = make_ids("SEA", 61:66),
      "SEA_SET9" = make_ids("SEA", 67:72),
      "SEA_SET10" = make_ids("SEA", 73:77),
      "SEA_SET11" = make_ids("SEA", 78:86),
      "SEA_SET12" = make_ids("SEA", 87:90),
      "SWP_SET1" = make_ids("SWP", 1:9),
      "SWP_SET2" = make_ids("SWP", 10:16),
      "SWP_SET3" = make_ids("SWP", 17:26),
      "SWP_SET4" = make_ids("SWP", 27:32)
    )

    id_to_group <- imap_dfr(match_list, ~ tibble(submap_transect_id = .x, submap_transect_set = .y)) |> unnest(submap_transect_set)
    sbmp <- sf_submap |>
      left_join(id_to_group, by = "submap_transect_id") |>
      relocate(submap_transect_set, .before = harc1)
    submap_transect_sets <- unique(na.omit(sbmp$submap_transect_set))

    hulls <- map(submap_transect_sets, ~ {
      sbmp_set <- filter(sbmp, submap_transect_set == .x)

      if (nrow(sbmp_set) == 0) {
        return(st_as_sf(tibble(submap_transect_set = .x, submap_transect_ids = NA_character_, hull_geometry = NULL), crs = 4326))
      }

      proj <- determine_projection(sbmp_set, proj_type = "aeqd")
      sbmp_proj <- sbmp_set |>
        st_transform(proj$wkt) |>
        rowwise() |>
        mutate(buff_proj = st_buffer(st_transform(submap_transect_geometry, proj$wkt), buffer_size, endCapStyle = "FLAT")) |>
        ungroup()

      shp <- c(st_make_valid(sbmp_proj$buff_proj)) |>
        st_combine() |>
        st_union()

      if (type == "convex") {
        hull <- st_convex_hull(shp) |> st_transform(crs = 4326)
      } else if (type == "concave") {
        pts <- st_cast(shp, "POINT") |> st_make_valid()
        hull <- concaveman(pts, concavity = concavity, length_threshold = 0) |> st_transform(crs = 4326)
      } else {
        hull <- st_convex_hull(shp) |> st_transform(crs = 4326)
      }

      tibble(submap_transect_set = .x, submap_transect_ids = paste0(sbmp_set$submap_transect_id, collapse = ","), geometry = hull) |>
        st_as_sf(crs = 4326)
    })

    sf_use_s2(TRUE)
    bind_rows(hulls) |>
      filter(!st_is_empty(geometry)) |>
      rename(hull_geometry = geometry)
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_ihfc2024 <- function(shp, out_path) {
  with_error_handling({
    if (check_data_path(out_path)) {
      return(invisible())
    }
    cat(" -> Processing IHFC 2024 observations globally\n")

    data <- shp |>
      st_set_geometry(NULL) |>
      mutate(across(where(is.numeric), ~ formatC(.x, format = "f", digits = 3)))

    cat("    Writing data: ", out_path, "\n", sep = "")
    write_csv(data, out_path)
    invisible()
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_lucazeau2019 <- function(shp, out_path) {
  with_error_handling({
    if (check_data_path(out_path)) {
      return(invisible())
    }
    cat(" -> Processing Lucazeau (2019) similarity interpolation globally\n")

    data <- shp |>
      st_set_geometry(NULL) |>
      mutate(across(where(is.numeric), ~ formatC(.x, format = "f", digits = 3)))

    cat("    Writing data: ", out_path, "\n", sep = "")
    write_csv(data, out_path)
    invisible()
  })
}
