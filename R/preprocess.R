#!/usr/bin/env Rscript

#######################################################
## Preprocess spatial datasets                       ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 3) {
    cat("    --------------------------------------------------\n")
    cat(" !! Usage: Rscript preprocess.R [util_dir] [data_dir] [out_dir]\n")
    return(invisible())
  }

  util_dir <- args[1]
  data_dir <- args[2]
  out_dir <- args[3]

  lapply(list.files(util_dir, pattern = "\\.R$", full.names = TRUE), source)

  dirs <- list(
    hf = file.path(data_dir, "heatflow"),
    plates = file.path(data_dir, "mapping", "plates"),
    volc = file.path(data_dir, "mapping", "volcanoes"),
    submap = file.path(data_dir, "mapping", "submap"),
    etopo = file.path(data_dir, "mapping", "relief"),
    relief = file.path(out_dir, "relief")
  )

  files <- list(
    hf = c("IHFC_2024_GHFDB.xlsx", "HFgrid14.csv", "goutorbe2011-similarity.csv"),
    plates = c("ridge.gmt", "transform.gmt", "trench.gmt"),
    volc = c("volcanoes.txt", "volcanoes_nospread.txt"),
    etopo = c(
      "ETOPO_2022_v1_60s_N90W180_surface.tif",
      "ETOPO_2022_v1_60s_N90W180_geoid.tif",
      "ETOPO_2022_v1_30s_N90W180_surface.tif",
      "ETOPO_2022_v1_30s_N90W180_geoid.tif"
    )
  )

  paths <- map2(files, names(files), ~ file.path(dirs[[.y]], .x))
  data_exist <- all(map_lgl(paths, ~ all(file.exists(.x)))) && dir.exists(dirs$submap)

  cat("    --------------------------------------------------\n")
  cat("    Fetching data from OSF repo: https://osf.io/ca6zu\n")
  cat("    --------------------------------------------------\n")
  if (!data_exist) {
    options(timeout = 1200)
    download_data_from_osf(data_dir)
  } else {
    imap(paths, ~ cat(" -- Found ", .y, " files:\n", paste0("    ", .x, collapse = "\n"), "\n", sep = ""))
    cat(" -- Found submap directory: ", dirs$submap, " (", length(list.files(dirs$submap)), " files)\n", sep = "")
  }

  cat("    --------------------------------------------------\n")
  cat("    Processing map data\n")
  cat("    --------------------------------------------------\n")
  out_path <- file.path(out_dir, "map-data.RData")

  if (!file.exists(out_path)) {
    create_path_dir(out_path)

    sf_coast <- fetch_coastline()
    sf_ridge <- fetch_plate_boundaries(paths$plates[1])
    sf_trench <- fetch_plate_boundaries(paths$plates[2])
    sf_transform <- fetch_plate_boundaries(paths$plates[3])
    sf_volc <- fetch_volcanoes(dirs$volc)
    sf_ihfc_raw <- fetch_ihfc(paths$hf[1])
    sf_ihfc <-  filter(sf_ihfc_raw, ihfc2024_obs > 0 & ihfc2024_obs <= 250)
    sf_sim <- fetch_lucazeau(paths$hf[2], c(0, 250))
    sf_grid <- st_geometry(sf_sim)
    sf_submap <- fetch_submap_transects(dirs$submap)
    sf_hull <- draw_submap_hulls(sf_submap)

    cat(" -> Saving data: ", out_path, "\n", sep = "")
    save(sf_coast, sf_ridge, sf_trench, sf_transform, sf_volc, sf_ihfc, sf_ihfc_raw, sf_sim, sf_grid, sf_submap, sf_hull, file = out_path)


    out_path_ihfc <- file.path(out_dir, "heatflow", "ihfc2024-obs.csv")
    out_path_sim <- file.path(out_dir, "heatflow", "lucazeau2019-sim.csv")

    write_ihfc2024(sf_ihfc_raw, out_path_ihfc)
    write_lucazeau2019(sf_sim, out_path_sim)
  } else {
    cat(" -- Found data: ", out_path, "\n", sep = "")
  }
}

if (!interactive() && (sys.nframe() == 0 || identical(environment(), globalenv()))) {
  main()
}
