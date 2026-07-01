#!/usr/bin/env Rscript

#######################################################
## Krige submap transects                            ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 4) {
    cat("    --------------------------------------------------\n")
    cat(" !! Usage: Rscript test.R [util_dir] [data_dir] [out_dir] [fig_dir]\n")
    return(invisible())
  }

  util_dir <- args[1]
  data_dir <- args[2]
  out_dir <- args[3]
  fig_dir <- args[4]

  lapply(list.files(util_dir, pattern = "\\.R$", full.names = TRUE), source)

  map_data <- file.path(out_dir, "map-data.RData")
  load_data(map_data)

  submap_transect_sets <- unique(sf_hull$submap_transect_set)
  process_submap_transect_sets(data_dir, out_dir, submap_transect_sets, quiet = TRUE)
}

if (!interactive() && (sys.nframe() == 0 || identical(environment(), globalenv()))) main()
