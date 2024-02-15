#!/usr/bin/env Rscript

# Quiet loading
sshhh <- function(p) {
  suppressWarnings(suppressPackageStartupMessages({
    require(p, quietly = T, character.only = TRUE)
  }))
}

# Install dependencies
check_dependencies <- function(packages) {
  # Try loading required packages
  pkgs <- unlist(list(packages))
  req <- unlist(lapply(pkgs, sshhh))
  need <- pkgs[req == FALSE]
  cat('Checking for required R packages:')
  cat('\n', pkgs)

  # Try installing required packages
  if(length(need) > 0) {
    install.packages(need)
    lapply(need, sshhh)
  }

  # Check installs
  req <- unlist(lapply(pkgs, sshhh))
  need <- pkgs[req == FALSE]

  if(length(need) > 0) {
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
package.list <- c('tictoc', 'stringr', 'tidyr', 'readr', 'purrr', 'furrr', 'tibble', 'dplyr',
                  'magrittr', 'ggplot2', 'metR', 'colorspace', 'ggrepel', 'ggridges',
                  'ggnewscale', 'patchwork', 'cowplot', 'ggsflabel', 'marmap', 'mapproj',
                  'gstat', 'rgeos', 'sf', 'stars', 'rnaturalearth', 'nloptr', 'future', 'zoo',
                  'sp')

# Check dependencies and install missing packages
check_dependencies(package.list)

# Print session info
cat('\n', rep('~', 45), '\n', sep='')
sessionInfo()

cat('\nAll required packages installed and available!')
cat('\ncheck-packages.R complete!\n')