#!/usr/bin/env Rscript

#######################################################
## Visualize                                         ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 4) {
    cat("    --------------------------------------------------\n")
    cat(" !! Usage: Rscript visualize.R [util_dir] [data_dir] [out_dir] [fig_dir]\n")
    return(invisible())
  }

  util_dir <- args[1]
  data_dir <- args[2]
  out_dir <- args[3]
  fig_dir <- args[4]

  lapply(list.files(util_dir, pattern = "\\.R$", full.names = TRUE), source)

  cat("    --------------------------------------------------\n")
  cat("    Visualizing global heatflow database\n")
  cat("    --------------------------------------------------\n")
  map_data <- file.path(out_dir, "map-data.RData")
  load_data(map_data)
  draw_global_dataset_composition(data_dir, out_dir, fig_dir)

  cat("    --------------------------------------------------\n")
  cat("    Visualizing nlopt summary\n")
  cat("    --------------------------------------------------\n")
  nlopt_data <- file.path(out_dir, "nlopt-data.RData")
  load_data(nlopt_data)
  draw_nlopt_summary(out_dir, fig_dir)
  draw_interpolation_accuracy_summary(out_dir, fig_dir)

  cat("    --------------------------------------------------\n")
  cat("    Visualizing interpolations\n")
  cat("    --------------------------------------------------\n")
  submap_transect_sets <- unique(sf_hull$submap_transect_set)
  draw_submap_transect_sets(data_dir, out_dir, fig_dir, submap_transect_sets)
}

if (!interactive() && (sys.nframe() == 0 || identical(environment(), globalenv()))) {
  main()
}
