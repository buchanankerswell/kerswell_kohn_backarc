#!/usr/bin/env Rscript

# Capture output
sink(file = paste0('data/log-', Sys.Date()), append = T, type = 'output', split = T)
# Set download timeout to 10min
cat(rep('~', 80), sep='')
cat('\nSetting download timeout to 10 minutes ...')
options(timeout = 600)
# Download data from osf
# https://osf.io/ca6zu/files/osfstorage
data.url <- 
  'https://files.osf.io/v1/resources/ca6zu/providers/osfstorage/639362fd084e0f013dcdaba1/?zip='
# Download .zip file
cat('\nDownloading data from osf ...')
cat('\nurl: https://osf.io/ca6zu/files/osfstorage')
cat('\nThis may take > 1 minute ...')
download.file(data.url, 'data.zip', quiet = T)
# Extract
unzip('data.zip', exdir = 'data')
# Remove .zip file
file.remove('data.zip')
# Write log
cat('download-data.R complete!\n')
sink()