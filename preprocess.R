#!/usr/bin/env Rscript

# Capture output
sink(file = 'data/preprocess_log', type = 'output', split = T)

cat(rep('~', 60), '\n', sep='')
sessionInfo()
cat('\n', rep('~', 60), '\n', sep='')

# Load packages and functions
cat('Loading packages and functions ...\n\n')

source('functions.R')

# Color scales
v.scale.white <-
  scale_color_viridis_c(
    option = 'magma',
    limits = c(0, 250),
    na.value = 'white'
  )
v.scale.grey <-
  scale_color_viridis_c(
    option = 'magma',
    limits = c(0, 250),
    na.value = 'grey50'
  )
v.scale.trans <-
  scale_color_viridis_c(
    option = 'magma',
    limits = c(0, 250),
    na.value = 'transparent'
  )

# Projections
# WGS84
proj4.wgs <-
  '+proj=longlat +lon_wrap=180 +ellps=WGS84 +datum=WGS84 +no_defs'

# Robinson Pacific centered
proj4.rp <-
  '+proj=robin +lon_0=-155 +lon_wrap=-155 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'

# Print Projections
cat('\n', rep('~', 60), '\n', sep='')
cat(
  'Defining coordinate reference system ...\n',
  '\nWGS:\n',
  proj4.wgs,
  '\nRobinson Pacific Centered:\n',
  proj4.rp
)

# SZ features to include in visualizations
shp.fts <-
  tibble(
    lat =
      c(2.5, 12.5, 2, 6, 15, 5, -60, -57.5,
        -65, -57.5, -56.5, 5, 15, -10, -16.5, -15, -18),
    lon =
      c(258, 255, 270, 263, 287.5, 300, 329, 340,
        345, 322, 332.5, 115, 130, 90, 172, 176, 176),
    segment =
      c('Central America', 'Central America', 'Central America',
        'Central America', 'Lesser Antilles',
        'Lesser Antilles', 'Scotia', 'Scotia', 'Scotia', 'Scotia',
        'Scotia', 'Sumatra Banda Sea', 'Sumatra Banda Sea',
        'Sumatra Banda Sea', 'Vanuatu', 'Vanuatu', 'Vanuatu'),
    label =
      c('GTJ', 'EPR', 'CR', 'CP', 'CBP', 'SA', 'ESR', 'TF',
        'AP', 'SP', 'SAN', 'SNP', 'PSP', 'AUP', 'NHP', 'BR', 'CWR')
  ) %>%
  st_as_sf(coords = c(2, 1), crs = proj4.wgs) %>%
  st_transform(proj4.rp)

# Country Boundaries

# Create sliver around twenty-five degree long to
# cut country boundaries so they don't project
# as straight lines across the globe
shp.sliver <-
  st_polygon(
    x = list(rbind(
      c(25.001, 90),
      c(25, 90),
      c(25, -90),
      c(25.001, -90),
      c(25.001, 90))
    )
  ) %>%
  st_sfc() %>%
  st_set_crs(proj4.wgs)

# Trim the countries boundaries overlapping at zero degree long and reproject
shp.countries <-
  ne_countries(returnclass = 'sf')

suppressWarnings({
  suppressMessages({
    shp.world <-
      shp.countries %>%
      st_difference(shp.sliver) %>%
      as_tibble() %>%
      st_as_sf() %>%
      st_transform(proj4.rp)
  })
})

cat('\n\n', rep('~', 60), '\n', sep='')
cat('Reading UTIG plate boundary data ...\n')
# Plate boundaries from UTIG
# http://www-udc.ig.utexas.edu/external/plates/data.htm
shp.ridge <- st_read('data/utig/ridge.gmt', crs = proj4.wgs, quiet = T)
shp.trench <- st_read('data/utig/trench.gmt', crs = proj4.wgs, quiet = T)
shp.transform <- st_read('data/utig/transform.gmt', crs = proj4.wgs, quiet = T)

# Filter out linestrings with single points
single.pnt.idx <-
  shp.ridge$geometry %>%
  map_lgl(~nrow(st_coordinates(.x)) <= 1)

# Fix dateline wrapping and transform
suppressWarnings({
  suppressMessages({
    shp.ridge <-
      shp.ridge[!single.pnt.idx,] %>%
      st_wrap_dateline() %>%
      st_difference(shp.sliver) %>%
      as_tibble() %>%
      st_as_sf() %>%
      st_transform(proj4.rp)
    shp.trench <-
      shp.trench %>%
      st_wrap_dateline() %>%
      st_difference(shp.sliver) %>%
      as_tibble() %>%
      st_as_sf() %>%
      st_transform(proj4.rp)
    shp.transform <-
      shp.transform %>%
      st_wrap_dateline() %>%
      st_difference(shp.sliver) %>%
      as_tibble() %>%
      st_as_sf() %>%
      st_transform(proj4.rp)
  })
})

# Read Syracuse et al (2006) volcanoes
cat('\n\n', rep('~', 60), '\n', sep='')
cat('Reading Syracuse and Abers (2006) data ...\n')

volc <-
  read_table(
    'data/sa2006/volcanoes.txt',
    col_type = c('ddddddddddddccf'),
    na = 'NULL'
  )

volc.no.spread <-
  read_table(
    'data/sa2006/volcanoes_nospread.txt',
    col_type = c('ddddddddddddccf'),
    na = 'NULL'
  )

volc <- volc %>% bind_rows(volc.no.spread)

# Change to simple features object and project
shp.volc <-
  volc %>%
  st_as_sf(coords = c('Lon', 'Lat')) %>%
  st_set_crs(proj4.wgs) %>%
  st_transform(proj4.rp)

cat('\nFound', nrow(volc), 'volcanoes\n')

# Read syracuse et al 2006 Segments
files <-
  list.files('data/sa2006/gmts', full.names = TRUE)

seg.names <-
  files %>%
  str_extract('(?<=gmts\\/)[a-z].*(?=_contours\\.gmt)') %>%
  str_replace_all('_', ' ') %>%
  str_to_title()

cat('\nFound segments:', seg.names, sep = '\n')

# Read contours and project to robinson pacific centered
shp.contours <-
  files %>%
  map(~read_latlong(.x, proj4.rp)) %>%
  set_names(seg.names)

# First contour defines the segment boundaries
shp.segs <-
  shp.contours %>%
  map(~.x[1,])

# Read ThermoGlobe database from Lucazeau (2019)
cat('\n', rep('~', 60), sep='')
cat('\nReading ThermoGlobe data from Lucazeau (2019) ...\n')

# Read ThermoGlobe data from Lucazeau (2019)
hf <- 
  read_csv('data/tglobe/tglobe.csv', col_type = cols()) %>%
  select(longitude, latitude, `heat-flow (mW/m2)`, code6) %>%
  rename(hf = `heat-flow (mW/m2)`)

cat(
  '\nRemoving', nrow(hf[is.na(hf$hf),]), 'NA heat flow values',
  '\nRemoving', nrow(hf[hf$code6 == 'D',]), 'poor quality data (code6 = D)',
  '\nRemoving', nrow(hf[is.na(hf$latitude) | is.na(hf$longitude),]), 'obs without coords\n'
)

# Transform to sf object
shp.hf <-
  hf %>%
  filter(!is.na(longitude)) %>% # remove no lat
  filter(!is.na(latitude)) %>% # remove no long
  st_as_sf(coords = c(1,2), crs = proj4.wgs) %>% # make sf object
  st_transform(proj4.rp)

# Filter bad quality data and missing heat flow
shp.hf.filtered <-
  hf %>% 
  filter(!is.na(hf)) %>% # remove NA
  filter(code6 != 'D') %>% # remove poor quality data
  filter(!is.na(longitude)) %>% # remove no lat
  filter(!is.na(latitude)) %>% # remove no long
  filter(hf > 0 & hf <= 250) %>%
  st_as_sf(coords = c(1,2), crs = proj4.wgs) %>% # make sf object
  st_transform(proj4.rp)

# Check for duplicate measurements and remove
dup <- sp::zerodist(sf::as_Spatial(shp.hf.filtered))
n.dup <- nrow(dup)

cat('\nParsing', nrow(dup), 'duplicate measurements:\n')
cat('\n', rep('+', 40), sep='')
cat(
  '\nIf x is better quality than y, keep x',
  '\notherwise randomly keep either x or y'
)
cat('\n', rep('+', 40), sep='')

# For duplicate measurements, select the best quality data point,
# if the quality is the same, randomly select one measurement
rid <- vector('integer', nrow(dup))

for(i in 1:nrow(dup)) {
  if(
    shp.hf.filtered$code6[dup[i,1]] != shp.hf.filtered$code6[dup[i,2]] &
    shp.hf.filtered$code6[dup[i,1]] > shp.hf.filtered$code6[dup[i,2]]
  ) {
    rid[i] <- dup[i,1]
  } else if(
    shp.hf.filtered$code6[dup[i,1]] != shp.hf.filtered$code6[dup[i,2]] &
    shp.hf.filtered$code6[dup[i,1]] < shp.hf.filtered$code6[dup[i,2]]
  ) {
    rid[i] <- dup[i,2]
  } else {
    rid[i] <- dup[i,sample(1:2, 1)]
  }
}

cat('\n\nParsed', length(rid), 'duplicates')

shp.hf.filtered <-
  shp.hf.filtered %>%
  slice(-rid) %>%
  select(hf, geometry) %>%
  rename(hf = hf)

cat('\n\nFinal ThermoGlobe dataset contains', nrow(shp.hf.filtered), 'observations')

# Read Lucazeau (2019) interpolation
cat('\n\n', rep('~', 60), sep='')
cat('\nReading similarity interpolation from Lucazeau (2019) ...')

interp.luca <-
  read_delim(
    'data/tglobe/HFgrid14.csv',
	  delim = ';',
	  col_types = c('ddddd')
  ) %>%
  rename(
    lon = longitude,
    lat = latitude,
    est.sim = HF_pred,
    sigma.sim = sHF_pred,
    obs.sim = Hf_obs
  ) %>%
  filter(est.sim > 0 & est.sim <= 250)

# Make simple feature object and project to robinson pacific centered
shp.interp.luca <-
  interp.luca %>% 
  st_as_sf(coords = c(1,2), crs = proj4.wgs) %>% 
  st_transform(proj4.rp)

# Extract 0.5˚ x 0.5˚ grid from Lucazeau (2019)
cat('\n\nExtracting similarity interpolation grid locations')

shp.grid <- st_geometry(shp.interp.luca)

cat('\nTotal similarity grid size:', length(shp.grid))

# Defining interpolation domain
cat('\n\n', rep('~', 60), sep='')
cat('\nDefining interpolation domain ...\n')

# Draw 1000 km buffer around segment boundaries
buf.dist <- 1e6

cat('\nDrawing', buf.dist/1000, 'km buffers around segments')

shp.buffer <-
  shp.segs %>%
  map(~st_buffer(.x, buf.dist, endCapStyle = 'ROUND'))

shp.box <-
  shp.buffer %>%
  map(~st_bbox(.x) %>% bbox_widen(crs = proj4.rp))

# Crop UTIG plate boundaries to bounding boxes
cat('\nCropping UTIG plate boundaries to segment bounding boxes')

shp.ridge.crop <-
  suppressWarnings({
    shp.box %>%
    map(~shp.ridge %>% st_crop(.x))
  })
shp.trench.crop <-
  suppressWarnings({
    shp.box %>%
    map(~shp.trench %>% st_crop(.x))
  })
shp.transform.crop <-
  suppressWarnings({
    shp.box %>%
    map(~shp.transform %>% st_crop(.x))
  })

# Crop ThermoGlobe data to bounding boxes
cat('\nCropping ThermoGlobe data to buffers')

shp.hf.crop <-
  suppressWarnings({
    shp.buffer %>%
    map(~shp.hf.filtered %>% st_intersection(.x) %>% select(-segment))
  })

cat('\nNumber of observations by segment:')
cat('\n\n', rep('+', 40), sep='')
walk2(shp.hf.crop, seg.names, ~{
  cat('\n', ..2, ':', nrow(..1))
})
cat('\n', rep('+', 40), sep='')

# Filtered and cropped hf tibble by segment
hf.crop <-
  shp.hf.crop %>%
  map_df(
    ~st_set_geometry(.x, NULL),
    .id = 'segment'
  )

# Cropped grids from Lucazeau (2019) to bounding boxes
cat('\n\nCropping interpolation grids to buffers')

shp.grid.crop <-
  shp.buffer %>%
  map(~shp.grid %>% st_intersection(.x))

cat('\nNumber of Kriging estimates by segment:')
cat('\n\n', rep('+', 40), sep='')
walk2(shp.grid.crop, seg.names, ~{
  cat('\n', ..2, ':', length(..1))
})
cat('\n', rep('+', 40), sep='')

# Calculate regional Similarity RMSE
cat('\n\nCalculating Similarity RMSE')
rmse.luca <-
  tibble(
    segment = seg.names,
    rmse =
      map_dbl(seg.names, ~suppressWarnings(
        itp_rmse(.x, shp.interp.luca, 'sim')
      ))
  )
# Summarize heat flow data
cat('\n\nHeat flow summary:\n')
hf.summary <-
  shp.hf.crop %>%
  map_df(
    ~st_set_geometry(.x, NULL),
    .id = 'segment'
  ) %>%
  group_by(segment) %>%
  summarise(
    n = n(),
    min = round(min(hf)),
    max = round(max(hf)),
    median =round(median(hf)),
    IQR = round(IQR(hf)),
    mean = round(mean(hf)),
    sd = round(sd(hf))
  )

print(hf.summary)

# Clean up environment
cat('\n', rep('~', 60), sep='')
cat('\nCleaning up environment ...')

rm(list = lsf.str())

rm(
   files,
   volc,
   proj4.rp,
   proj4.wgs,
   i,
   shp.box,
   shp.sliver,
   shp.hf,
   shp.hf.filtered,
   shp.grid,
   shp.countries,
   volc.no.spread,
   dup,
   rid,
   buf.dist,
   interp.luca,
   single.pnt.idx
)

# Save
cat('\nSaving data to: data/hf.RData')
save.image('data/hf.RData')

cat('\nDone!\n')
sink()
