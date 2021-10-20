#!/usr/bin/env Rscript

# Load packages and functions
cat(rep('~', 60), '\n', sep='')
cat('Loading packages and functions ...\n\n')

source('functions.R')

# Color scales
v.scale.white <-
  scale_color_viridis_c(
    option = 'magma',
    limits = c(0, 250),
    na.value = rgb(0.99, 0.99, 0.90)
  )
v.scale.grey <-
  scale_color_viridis_c(
    option = 'magma',
    limits = c(0, 250),
    na.value = 'grey50'
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
      c(2.5, 12.5, 2, 6, 27.5, 15, 5, -60, -57.5,
        -65, -60, -56.5, 5, 20, -10, -20, -12, -17.5),
    lon =
      c(258, 255, 270, 263, 310, 287.5, 300, 329, 340,
        345, 320, 332.5, 115, 130, 90, 171, 175, 174),
    segment =
      c('Central America', 'Central America', 'Central America',
        'Central America', 'Lesser Antilles', 'Lesser Antilles',
        'Lesser Antilles', 'Scotia', 'Scotia', 'Scotia', 'Scotia',
        'Scotia', 'Sumatra Banda Sea', 'Sumatra Banda Sea',
        'Sumatra Banda Sea', 'Vanuatu', 'Vanuatu', 'Vanuatu'),
    label =
      c('GTJ', 'EPR', 'CR', 'CP', 'MAR', 'CBP', 'SA', 'ESR', 'TF',
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
      c(25.000001, 90),
      c(25, 90),
      c(25, -90),
      c(25.000001, -90),
      c(25.000001, 90))
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
hf <- read_csv('data/tglobe/tglobe.csv', col_type = cols())

cat(
  '\nRemoving', nrow(hf[is.na(hf$`heat-flow (mW/m2)`),]), 'NA heat flow values',
  '\nRemoving', nrow(hf[hf$code6 == 'D',]), 'poor quality data (code6 = D)',
  '\nRemoving', nrow(hf[is.na(hf$latitude) | is.na(hf$longitude),]), 'obs without coords\n'
)

# Transform to sf object
shp.hf <-
  hf %>%
  filter(!is.na(longitude)) %>% # remove no lat
  filter(!is.na(latitude)) %>% # remove no long
  st_as_sf(coords = c(9,10), crs = proj4.wgs) %>% # make sf object
  st_transform(proj4.rp)

# Filter bad quality data and missing heat flow
shp.hf.filtered <-
  hf %>% 
  filter(!is.na(`heat-flow (mW/m2)`)) %>% # remove NA
  filter(code6 != 'D') %>% # remove poor quality data
  filter(!is.na(longitude)) %>% # remove no lat
  filter(!is.na(latitude)) %>% # remove no long
  st_as_sf(coords = c(9,10), crs = proj4.wgs) %>% # make sf object
  st_transform(proj4.rp)

# Check for duplicate measurements and remove
dup <- sp::zerodist(sf::as_Spatial(shp.hf.filtered))

cat('\nParsing', length(dup), 'duplicate measurements:\n')
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
  select(`heat-flow (mW/m2)`, geometry) %>%
  rename(hf = `heat-flow (mW/m2)`)

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
    lon = longiyude,
    lat = latitude,
    est.sim = HF_pred,
    sigma.sim = sHF_pred,
    obs.sim = Hf_obs
  )

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
buf.dist <- 1000000

cat('\nDrawing', buf.dist/1000, 'km buffers around segments')

shp.buffer <-
  shp.segs %>%
  map(~st_buffer(.x, buf.dist, endCapStyle = 'ROUND'))

# shp.box <-
#   shp.buffer %>%
#   map(~st_bbox(.x) %>% bbox_widen(crs = proj4.rp))

# Crop ThermoGlobe data to bounding boxes
cat('\nCropping ThermoGlobe data to buffers')

shp.hf.crop <-
  suppressWarnings({
#     shp.box %>%
    shp.buffer %>%
#     map(~shp.hf.filtered %>% st_crop(.x))
    map(~shp.hf.filtered %>% st_intersection(.x) %>% select(-segment))
  })

cat('\nNumber of observations by segment:')
cat('\n\n', rep('+', 40), sep='')
walk2(shp.hf.crop, seg.names, ~{
  cat('\n', ..2, ':', nrow(..1))
})
cat('\n', rep('+', 40), sep='')

# Cropped grids from Lucazeau (2019) to bounding boxes
cat('\n\nCropping interpolation grids to buffers')

shp.grid.crop <-
#   shp.box %>%
  shp.buffer %>%
#   map(~shp.grid %>% st_crop(.x))
  map(~shp.grid %>% st_intersection(.x))

cat('\nNumber of Kriging estimates by segment:')
cat('\n\n', rep('+', 40), sep='')
walk2(shp.grid.crop, seg.names, ~{
  cat('\n', ..2, ':', length(..1))
})
cat('\n', rep('+', 40), sep='')

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
   hf,
   proj4.rp,
   proj4.wgs,
   i,
   shp.sliver,
   shp.hf,
   shp.grid,
   shp.countries,
   volc.no.spread,
   dup,
   rid,
   buf.dist,
   interp.luca
)

# Save
cat('\n\nSaving data to: data/hf.RData')
save.image('data/hf.RData')

cat('\n\nDone!')
cat('\n', rep('~', 60), '\n', sep='')
