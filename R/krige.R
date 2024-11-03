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

  cat("    --------------------------------------------------\n")
  cat("    Kriging submap transects (after Li et al., 2018)\n")
  cat("    --------------------------------------------------\n")
  map_data <- file.path(out_dir, "map-data.RData")
  nlopt_data <- file.path(out_dir, "nlopt-data.RData")

  if (file.exists(nlopt_data)) {
    cat(" -- Found data: ", nlopt_data, "\n", sep = "")
    return(invisible())
  }

  load_data(map_data)
  submap_transect_sets <- unique(sf_hull$submap_transect_set)
  process_submap_transect_sets(data_dir, out_dir, submap_transect_sets, nprocs = 8)

  cat("    --------------------------------------------------\n")
  cat("    Summarizing results\n")
  cat("    --------------------------------------------------\n")
  nlopt_summary <- summarize_optimal_krige_models(out_dir, submap_transect_sets)
  interp_diff_summary <- summarize_interpolation_differences(data_dir, out_dir, submap_transect_sets)
  interp_accuracy_summary <- summarize_interpolation_accuracies(data_dir, out_dir, submap_transect_sets)

  if (!dir.exists(dirname(nlopt_data))) dir.create(dirname(nlopt_data), recursive = TRUE, showWarnings = FALSE)

  cat("    Writing data: ", nlopt_data, "\n", sep = "")
  save(nlopt_summary, interp_diff_summary, interp_accuracy_summary, file = nlopt_data)

  cat("    --------------------------------------------------\n")
  cat("    Writing heatflow data\n")
  cat("    --------------------------------------------------\n")
  suppressWarnings({
    write_submap_transect_set_ihfc2024_obs(data_dir, out_dir, submap_transect_sets)
    write_submap_transect_set_lucazeau2019_sim(data_dir, out_dir, submap_transect_sets)
    write_submap_transect_set_kerswell2025_krg(data_dir, out_dir, submap_transect_sets)
  })
}

if (!interactive() && (sys.nframe() == 0 || identical(environment(), globalenv()))) main()
