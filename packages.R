#!/usr/bin/env Rscript

# Quiet loading
sshhh <- function(p) {
  suppressWarnings(
    suppressPackageStartupMessages({
      cat('\nChecking required package:', p)
      require(p, quietly = T, character.only = TRUE)
    })
  )
}

# Install dependencies method
using <- function(...) {
  cat('\n', rep('~', 60), sep = '')
  # Try loading required packages
  pkgs <- unlist(list(...))
  req <- unlist(lapply(pkgs, sshhh))
  need <- pkgs[req == FALSE]
  # Try installing required packages
  if(length(need) > 0) {
    install.packages(need)
    lapply(need, sshhh)
  }
  # Check installs
  req <- unlist(lapply(pkgs, sshhh))
  need <- pkgs[req == FALSE]
  if(length(need) > 0) {
    cat('\n', rep('~', 60), sep = '')
    cat('\nFailed to install packages:\n')
    writeLines(need)
    cat('\nFor troubleshooting tips see:')
    cat('\nhttps://github.com/buchanankerswell/kerswell_kohn_backarc')
    cat('\n', rep('~', 60), '\n', sep = '')
    stop()
  }
}

# Package list
package.list <- c(
  'colorspace',
  'tictoc',
  'stringr',
  'tidyr',
  'readr',
  'purrr',
  'furrr',
  'tibble',
  'dplyr',
  'magrittr',
  'ggplot2',
  'ggrepel',
  'ggridges',
  'ggsflabel',
  'patchwork',
  'cowplot',
  'gstat',
  'rgeos',
  'sf',
  'rnaturalearth',
  'nloptr',
  'future',
  'zoo',
  'sp'
)

using(package.list)

cat('\n\nAll required packages installed and available!')
cat('\n', rep('~', 60), '\n\n', sep = '')