#!/usr/bin/env Rscript

#######################################################
## Prepare R environment                         !!! ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
main <- function() {
  if (!requireNamespace("renv", quietly = TRUE)) {
    install.packages("renv", quiet = TRUE)
  }

  lockfile <- "renv.lock"

  if (file.exists(lockfile)) {
    renv::restore(prompt = FALSE)
  } else {
    cran_pkgs <- c(
      "tidyverse", "readxl", "jsonlite", "rlang", "units", "nloptr", "lhs", "zoo",
      "mgcv", "quantreg", "furrr", "concaveman", "sp", "sf", "stars", "terra",
      "gstat", "tidyterra", "marmap", "ggOceanMaps", "scales", "ggsci", "ggspatial",
      "ggrepel", "ggnewscale", "ggridges", "ggpattern", "cowplot", "patchwork",
      "knitr", "osfr"
    )

    github_pkgs <- c("yutannihilation/ggsflabel", "ropensci/rnaturalearthhires")

    suppressWarnings({
      renv::install(c(cran_pkgs, github_pkgs))
      renv::snapshot(prompt = FALSE)
    })
  }
}

if (!interactive() && (sys.nframe() == 0 || identical(environment(), globalenv()))) {
  main()
}
