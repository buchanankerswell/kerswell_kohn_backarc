#!/usr/bin/env Rscript

#######################################################
## Sync data from OSF repo ca6zu                     ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 3) {
    cat("    --------------------------------------------------\n")
    cat(" !! Usage: Rscript visualize.R [util_dir] [data_dir] [sync_dir]\n")
    return(invisible())
  }

  util_dir <- args[1]
  data_dir <- args[2]
  sync_dir <- args[3]

  lapply(list.files(util_dir, pattern = "\\.R$", full.names = TRUE), source)

  if (sync_dir == "download") {
    download_simulation_results_from_osf(data_dir)
  } else if (sync_dir == "upload") {
    upload_simulation_results_to_osf(data_dir)
  } else {
    cat(" !! Error: unknown sync command '", sync_dir, "'. Use 'upload' or 'download'.\n", sep = "")
  }
}

if (!interactive() && (sys.nframe() == 0 || identical(environment(), globalenv()))) {
  main()
}
