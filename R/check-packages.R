#!/usr/bin/env Rscript

# Quiet loading
sshhh <- function(p) {
  suppressWarnings(suppressPackageStartupMessages({
    require(p, quietly = T, character.only = TRUE)
  }))
}

# Install dependencies
check_dependencies <- function(packages) {
  pkgs <- unlist(list(packages))
  req <- unlist(lapply(pkgs, sshhh))
  need <- pkgs[req == FALSE]
  if (length(need) > 0) {
    install.packages(need)
    lapply(need, sshhh)
  }
  req <- unlist(lapply(pkgs, sshhh))
  need <- pkgs[req == FALSE]
  if (length(need) > 0) {
    cat('\n', rep('~', 45), sep = '')
    cat('\nFailed to install packages:\n')
    writeLines(need)
    cat('\nFor troubleshooting tips see:')
    cat('\nhttps://github.com/buchanankerswell/kerswell_kohn_backarc')
    cat('\n', rep('~', 45), '\n', sep = '')
    stop()
  }
}

# Package list
package_list <- c('tictoc', 'stringr', 'tidyr', 'readr', 'readxl', 'purrr', 'furrr',
                  'tibble', 'dplyr', 'magrittr', 'units', 'ggplot2', 'colorspace', 'metR',
                  'ggrepel', 'ggridges', 'ggnewscale', 'patchwork', 'cowplot', 'ggsflabel',
                  'marmap', 'scales', 'ggspatial', 'gstat', 'rgeos', 'sp', 'sf',
                  'rnaturalearth', 'nloptr', 'zoo', 'jsonlite')

# Check dependencies and install missing packages
cat('Checking required packages ...')
check_dependencies(package_list)

# Print session info
cat('\n', rep('~', 45), '\n', sep='')
sessionInfo()

cat('\nAll required packages installed and available!')
cat('\ncheck-packages.R complete!\n')