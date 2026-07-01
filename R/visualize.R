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

  map_data <- file.path(out_dir, "map-data.RData")
  load_data(map_data)
  draw_global_dataset_summary(data_dir, out_dir, fig_dir, quiet = TRUE)

  nlopt_data <- file.path(out_dir, "post-processed.RData")
  load_data(nlopt_data)
  draw_nlopt_summary(fig_dir, nlopt_summary, quiet = TRUE)

  submap_transect_sets <- unique(sf_hull$submap_transect_set)
  draw_submap_transect_sets(data_dir, out_dir, fig_dir, submap_transect_sets, quiet = TRUE)
  draw_point_by_point_summary(data_dir, out_dir, fig_dir, submap_transect_sets, quiet = TRUE)

  draw_global_ccf_summary(fig_dir, global_ccf_summary, method_ccf_summary, quiet = TRUE)
  draw_method_ccf_summary(fig_dir, global_ccf_summary, method_ccf_summary, quiet = TRUE)

  walk(submap_transect_sets, ~ draw_case_study_composition(data_dir, out_dir, fig_dir, .x, quiet = TRUE))
}

if (!interactive() && (sys.nframe() == 0 || identical(environment(), globalenv()))) {
  main()
}
