#!/usr/bin/env Rscript

#######################################################
## Postprocess                                       ##
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

  submap_transect_sets <- unique(sf_hull$submap_transect_set)
  write_ihfc2024_summary_tables(data_dir, out_dir, submap_transect_sets, quiet = TRUE)
  write_submap_transect_ihfc2024_summary_tables(data_dir, out_dir, submap_transect_sets, quiet = TRUE)

  post_processed_data <- file.path(out_dir, "post-processed.RData")
  if (file.exists(post_processed_data)) {
    load_data(post_processed_data)
  } else {
    nlopt_summary <- summarize_optimal_krige_models(out_dir, submap_transect_sets, quiet = TRUE)
    interp_est_diff_summary <- summarize_interpolation_est_differences(data_dir, out_dir, submap_transect_sets, quiet = TRUE)
    interp_sigma_diff_summary <- summarize_interpolation_sigma_differences(data_dir, out_dir, submap_transect_sets, quiet = TRUE)
    interp_accuracy_summary <- summarize_interpolation_accuracies(data_dir, out_dir, submap_transect_sets, quiet = TRUE)

    global_ccf_summary <- summarize_global_ccf(data_dir, out_dir, submap_transect_sets, quiet = TRUE)
    method_ccf_summary <- summarize_method_ccf(data_dir, out_dir, submap_transect_sets, quiet = TRUE)

    enriched_data <- add_transect_set_regime_classifications(global_ccf_summary, method_ccf_summary, quiet = TRUE)
    global_ccf_summary <- enriched_data$global
    method_ccf_summary <- enriched_data$method

    case_ccf_summary <- summarize_case_study_ccf(data_dir, out_dir, submap_transect_sets, quiet = TRUE)

    if (!dir.exists(dirname(post_processed_data))) dir.create(dirname(post_processed_data), recursive = TRUE, showWarnings = FALSE)
    save(
      nlopt_summary,
      interp_est_diff_summary,
      interp_sigma_diff_summary,
      interp_accuracy_summary,
      global_ccf_summary,
      method_ccf_summary,
      case_ccf_summary,
      file = post_processed_data
    )
  }

  write_submap_transect_set_ihfc2024_obs_csv(data_dir, out_dir, submap_transect_sets, quiet = TRUE)
  write_submap_transect_set_lucazeau2019_sim_csv(data_dir, out_dir, submap_transect_sets, quiet = TRUE)
  write_submap_transect_set_kerswell2025_krg_csv(data_dir, out_dir, submap_transect_sets, quiet = TRUE)

  write_nlopt_summary_tables(out_dir, nlopt_summary, quiet = TRUE)
  write_est_diff_summary_tables(out_dir, interp_est_diff_summary, quiet = TRUE)
  write_sigma_diff_summary_tables(out_dir, interp_sigma_diff_summary, quiet = TRUE)
  write_similarity_residual_summary_tables(out_dir, interp_accuracy_summary, quiet = TRUE)
  write_krige_residual_summary_tables(out_dir, interp_accuracy_summary, quiet = TRUE)
  write_global_ccf_summary_tables(out_dir, global_ccf_summary, quiet = TRUE)
  write_method_ccf_summary_tables(out_dir, method_ccf_summary, quiet = TRUE)
  write_case_study_ccf_summary_tables(out_dir, case_ccf_summary, quiet = TRUE)
}

if (!interactive() && (sys.nframe() == 0 || identical(environment(), globalenv()))) {
  main()
}
