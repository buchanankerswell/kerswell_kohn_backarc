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
      "tidyverse", "furrr", "nloptr", "lhs", "zoo", "concaveman", "sp", "sf", "gstat", "tidyterra",
      "marmap", "ggOceanMaps", "ggsci", "ggspatial", "ggrepel", "ggnewscale", "cowplot", "patchwork"
    )

    github_pkgs <- c("yutannihilation/ggsflabel", "ropensci/rnaturalearthhires")

    suppressWarnings({
      renv::install(c(cran_pkgs, github_pkgs))
      renv::snapshot(packages = c(cran_pkgs, "ggsflabel"), prompt = FALSE)
    })
  }
}

if (!interactive() && (sys.nframe() == 0 || identical(environment(), globalenv()))) {
  main()
}
