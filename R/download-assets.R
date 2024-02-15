#!/usr/bin/env Rscript

# Set download timeout to 10min
options(timeout = 600)

# Download data from osf
# https://osf.io/ca6zu/files/osfstorage
assets.url <- 'https://files.osf.io/v1/resources/ca6zu/providers/osfstorage/65cddc3e6c2a400199187f6b/?zip=&_gl=1*l3s9hi*_ga*NjQ1MzM5ODAxLjE3MDA2NTIyMTc.*_ga_YE9BMGGWX8*MTcwNzk4OTgzNC4yMi4xLjE3MDc5OTAyOTguNjAuMC4w'

# Download .zip file
cat('Downloading assets from osf: https://osf.io/ca6zu/files/osfstorage')
download.file(assets.url, 'assets.zip', quiet = T)

# Extract
unzip('assets.zip', exdir = 'assets')

# Remove .zip file
file.remove('assets.zip')

cat('\nAll assets downloaded and available!')
cat('download-assets.R complete!\n')