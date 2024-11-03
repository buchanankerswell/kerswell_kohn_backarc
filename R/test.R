#!/usr/bin/env Rscript

#######################################################
## Testing                                           ##
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

  # optimize_krige_model(data_dir, out_dir, "NPA_SET2", "Sph")
  draw_submap_transect_sets(data_dir, out_dir, fig_dir, c("NPA_SET2"))

}

if (!interactive() && (sys.nframe() == 0 || identical(environment(), globalenv()))) {
  main()
}
