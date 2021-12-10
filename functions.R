# Load packages
# Quiet loading
sshhh <- function(p) {
  suppressWarnings(
    suppressPackageStartupMessages(
      library(p, quietly = T, character.only=TRUE)
    )
  )
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
  'nloptr'
)

# auto-load quietly
sapply(package.list, sshhh)
rm(package.list, sshhh)

# Don't allow sf to use google's s2 library
# for spherical geometry
# This reverts to using the GEOS library instead
# which is what sf used before 1.0 release
sf_use_s2(FALSE)

# Calculate Similarity rmse by inverse distance weighting
# interpolation of observations to nearby grid points
sim_rmse <- function(seg.name, maxdist = 5.5e4, idp = 2, plot = F){
  buf <- shp.buffer[[seg.name]]
  grd <- shp.grid.crop[[seg.name]]
  obs <- shp.hf.crop[[seg.name]]
  sim <- st_intersection(shp.interp.luca, buf)
  obs.itp <-
    idw(
      formula = hf~1,
      locations = obs,
      newdata = grd,
      maxdist = maxdist,
      idp = idp,
      debug.level = 0
    ) %>%
    mutate(hf = var1.pred) %>%
    select(hf, geometry)
  idx <- !is.na(obs.itp$hf)
  if(plot){
    seg <- shp.segs[[seg.name]]
    world <- shp.world %>% st_crop(buf)
    ridge <- shp.ridge.crop[[seg.name]] %>% st_crop(buf)
    trench <- shp.trench.crop[[seg.name]] %>% st_crop(buf)
    transform <- shp.transform.crop[[seg.name]] %>% st_crop(buf)
    p <-
      ggplot() +
        geom_sf(data = world, size = 0.1, fill = 'grey60') +
        geom_sf(data = buf, fill = NA) +
        geom_sf(data = ridge, size = 0.5, alpha = 0.8) +
        geom_sf(data = trench, size = 0.5, alpha = 0.8) +
        geom_sf(data = transform, size = 0.5, alpha = 0.8) +
        geom_sf(data = seg, size = 2) +
        geom_sf(data = obs.itp, aes(color = hf), size = 1, shape = 19) +
        geom_sf(data = obs, color = 'white', shape = 19, size = 0.1) +
        labs(color = bquote(mWm^-2), fill = bquote(mWm^-2)) +
        scale_fill_viridis_c(
          option = 'magma',
          limits = c(0, 250),
          na.value = 'transparent'
        ) +
        scale_color_viridis_c(
          option = 'magma',
          limits = c(0, 250),
          na.value = 'transparent'
        ) +
        theme_map(font_size = 10) +
        theme(
          plot.tag = element_text(face = 'bold', size = 14),
          axis.text = element_text(),
          axis.text.x = element_text(angle = 30),
          panel.grid = element_line(size = 0.1, color = 'white'),
          panel.background = element_rect(fill = 'grey50', color = NA),
          plot.margin = margin()
        )
    print(p)
  }
  return(sqrt(sum((sim[idx,]$est.sim - obs.itp[idx,]$hf)^2)/length(idx)))
}

# Draw a widened box from a st_bbox object
bbox_widen <-
  function(
    bbox,
    crs = NULL,
    borders =
      c(
        'left' = 0,
        'right' = 0,
        'top' = 0,
        'bottom' = 0
      )
  ) {
  # Check for missing arguments
  if(is.null(bbox)) stop('\nMissing bounding box!')
  if(is.null(crs)){
    crs <- st_crs(bbox)
  }
  b <- bbox # current bounding box
  xrange <- b$xmax - b$xmin # range of x values
  yrange <- b$ymax - b$ymin # range of y values
  b[1] <- b[1] - (borders[1] * xrange) # xmin - left
  b[3] <- b[3] + (borders[2] * xrange) # xmax - right
  b[4] <- b[4] + (borders[3] * yrange) # ymax - top
  b[2] <- b[2] - (borders[4] * yrange) # ymin - bottom
  box <- st_polygon(
    list(
      matrix(
        c(
          b$xmin,
          b$ymax,
          b$xmin,
          b$ymin,
          b$xmax,
          b$ymin,
          b$xmax,
          b$ymax,
          b$xmin,
          b$ymax
        ),
        ncol = 2,
        byrow = TRUE
      )
    )
  ) %>%
  st_sfc(crs = crs)
  return(box)
}

# Read gmt files, wrap the dateline to avoid plotting horizontal lines on map,
# make into tibble, add segment names, transform projection, and bind into one sf object
read_latlong <- function(file = NULL, crs = NULL) {
  # Check for missing arguments
  if(is.null(file)) stop('\nMissing filename!')
  if(is.null(crs)) stop('\nMissing coordinate reference system!')
  # Parse seg.name from filename
  seg.name <-
    file %>%
    str_extract('(?<=gmts\\/)[a-z].*(?=_contours\\.gmt)') %>%
    str_replace_all('_', ' ') %>%
    str_to_title()
  # Read txt file and project to WGS
  st_read(
    file,
    crs = '+proj=longlat +lon_wrap=180 +ellps=WGS84 +datum=WGS84 +no_defs',
    quiet = TRUE
  ) %>%
  # Transform projection
  st_transform(crs) %>%
  # Make into tibble for better printing
  tibble::as_tibble() %>%
  # Make back into sf object
  st_as_sf() %>%
  # Add segment column with segment name from file
  mutate(segment = seg.name, .before = geometry)
}

# Calculating experimental variograms
experimental_vgrm <-
  function(
    shp.hf = NULL,
    cutoff.prop = 3,
    n.lags = 20,
    lag.start = 1
  ) {
  # Check for missing arguments
  if(is.null(shp.hf)) stop('\nMissing heat flow data!')
  # Get bounding box
  bbox <- st_bbox(shp.hf)
  # Calculate corner-to-corner distance of bbox
  bbox.diagonal.distance <-
    sqrt((bbox$xmax-bbox$xmin)^2 + (bbox$ymax-bbox$ymin)^2)
  # Lag cutoff is a proportion of the corner-to-corner bounding box distance
  # note: defaults to the bounding box diagonal divided by three
  lag.cutoff <- as.vector(bbox.diagonal.distance/cutoff.prop)
  bin.width <- lag.cutoff/n.lags
  # Calculate new cutoff for shifting variogram lag.start to the left
  shifted.cutoff <- bin.width*(n.lags+lag.start-1)
  # Return shifted experimental variogram
  variogram(
    hf~1,
    locations = shp.hf,
    cutoff = shifted.cutoff,
    width = bin.width
  ) %>%
  slice(floor(lag.start):n())
}

# Wrapper function around gstat krige method
Krige <-
  function(
    shp.hf = NULL,
    fitted.vgrm = NULL,
    shp.interp.grid = NULL,
    n.max = 20,
    seg.name = NULL
  ) {
  # Check for missing arguments
  if(is.null(shp.hf)) stop('\nMissing heat flow data!')
  if(is.null(fitted.vgrm)) stop('\nMissing variogram model!')
  if(is.null(shp.interp.grid)) stop('\nMissing kriging locations (grid)!')
  if(is.null(n.max)) stop('\nNumer of max local point-pairs!')
  # Print grid and parameters info
  cat(
    '\nKriging:', seg.name,
    '\nObservations:', nrow(shp.hf),
    '\nGrid size:', length(shp.interp.grid),
    '\nMax nearest points:', n.max,
    '\nVariogram params:\n'
  )
  print(fitted.vgrm)
  # Kriging
  k <-
    try(
      suppressWarnings(
        krige(
          formula = hf~1,
          locations = shp.hf,
          newdata = shp.interp.grid,
          model = fitted.vgrm,
          nmax = n.max,
          debug.level = 0
        ) %>%
        as_tibble() %>%
        st_as_sf() %>%
        rename(
          est.krige = var1.pred,
          var.krige = var1.var
        ) %>%
        mutate(
          sigma.krige = sqrt(var.krige),
          .before = geometry
        )
      )
    )
  # Check for error during fitting
  if(any(class(k) == 'try-error')){
    # Print grid and parameters info
    cat(
      '\nKriging:', seg.name,
      '\nObservations:', nrow(shp.hf),
      '\nGrid size:', length(shp.interp.grid),
      '\nMax nearest points:', n.max,
      '\nVariogram params:\n'
    )
    print(fitted.vgrm)
    stop('\nVariogram fitting error!')
  }
  return(k)
}

# Take difference
interp_diff <- function(shp.interp.krige, shp.interp.sim) {
  # Check for missing arguments
  if(is.null(shp.interp.krige)) stop('\nMissing krige data!')
  if(is.null(shp.interp.sim)) stop('\nMissing similarity data!')
  # Crop Similarity estimates to Kriging bounding box
  dif <-
    suppressWarnings(
      shp.interp.sim %>%
      st_intersection(shp.interp.krige) %>%
      mutate(
        est.diff = est.sim - est.krige,
        sigma.diff = sigma.sim - sigma.krige,
        .before = geometry
      )
    )
  return(dif)
}

# Decode output of nloptr method and
# construct optimized variogram models
decode_opt <- function(model.vgrm = NULL, shp.hf = NULL, opt = NULL) {
  # Check for missing arguments
  if(is.null(model.vgrm)) stop('\nMissing variogram model!')
  if(is.null(shp.hf)) stop('\nMissing heat flow data!')
  if(is.null(opt)) stop('\nMissing nloptr object!')
  # Compute experimental vgrm
  experimental.vgrm <-
    try(
      experimental_vgrm(
        shp.hf,
        cutoff.prop = opt$solution[1],
        n.lags = opt$solution[2],
        lag.start = opt$solution[3]
      ),
      silent = T
    )
  # Need more than one lag to fit variogram
  if(nrow(experimental.vgrm) < 2) {
    stop('\nExperimental variogram has less than two lags!')
  }
  if(any(class(experimental.vgrm) == 'try-error')) {
    stop('\nExperimental variogram error!')
  }
  # Fit experimental variogram
  fitted.vgrm <-
    try(
      fit.variogram(
        object = experimental.vgrm,
        fit.method = 7,
        debug.level = 0,
        model = vgm(model = model.vgrm)
      ),
      silent = T
    )
  # Check for error during fitting
  if(any(class(fitted.vgrm) == 'try-error')){
    print(fitted.vgrm)
    stop('\nVariogram fitting error!')
  }
  return(
    list(
      'experimental.vgrm' = experimental.vgrm,
      'fitted.vgrm' = fitted.vgrm
    )
  )
}

# Cost function for Kriging estimates
#   note: after Li et al 2018
#   see https://doi.org/10.1016/j.cageo.2018.07.011
cost_function <- 
  function(
    shp.hf,
    cutoff.prop = 3,
    n.lags = 30,
    lag.start = 3,
    model.vgrm = 'Sph',
    n.max = 20,
    n.fold = NULL,
    interp.weight = 0.5,
    vgrm.weight = 0.5,
    segment = NULL,
    verbose = T
  ) {
  # Check for missing arguments
  if(is.null(shp.hf)) stop('\nMissing heat flow data model!')
  if(is.null(n.fold)) n.fold <- nrow(shp.hf)
  if(n.fold > 0 & n.fold <= 1) n.fold <- nrow(shp.hf)*n.fold
  # Calculate experimental variogram
  experimental.vgrm <-
    try(
      experimental_vgrm(
        shp.hf = shp.hf,
        cutoff.prop = cutoff.prop,
        n.lags = n.lags,
        lag.start = lag.start
      ),
      silent = T
    )
  # Need more than one lag to fit variogram
  if(nrow(experimental.vgrm) < 2) {
    if(verbose) {
      cat('\nExperimental variogram has less than two lags!')
      cat('\nReturning arbitrarily high cost')
      return(runif(1, 1, 1.5))
    }
  }
  if(any(class(experimental.vgrm) == 'try-error')) {
    if(verbose) {
      cat('\nExperimental variogram error!')
      cat('\nReturning arbitrarily high cost')
      return(runif(1, 1, 1.5))
    }
  }
  # Fit experimental variogram
  fitted.vgrm <-
    try(
      fit.variogram(
        object = experimental.vgrm,
        fit.method = 7,
        debug.level = 0,
        model = vgm(model = model.vgrm)
      ),
      silent = T
    )
  # Check for error during fitting
  if(any(class(fitted.vgrm) == 'try-error')){
    if(verbose) {
      cat('\nVariogram fitting error!\n')
      cat('\nReturning arbitrarily high cost')
      return(runif(1, 1, 1.5))
    }
  }
  # Compute n-fold cross validation
  if(verbose){
    cat('\n', rep('+', 30), sep='')
    cat('\nComputing cross-validation')
  }
  k.cv <- 
    try(
      krige.cv(
        formula = hf~1,
        locations = shp.hf,
        model = fitted.vgrm,
        nmax = n.max,
        nfold = n.fold,
        verbose = F
      ),
      silent = T
    )
  # Handle errors
  # Return arbitrarily high cost if
  # k.cv throws an error to keep
  # minimization algorithm searching
  if(any(class(k.cv) == 'try-error')) {
    if(verbose){
      cat('\nCross-validation error!')
      cat('\nReturning arbitrarily high cost\n')
      return(runif(1, 1, 1.5))
    }
  }
  if(sum(is.na(k.cv)) != 0) {
    if(verbose){
      cat('\nCross-validation produced NAs!')
      if(sum(is.na(k.cv$residual)) >= nrow(k.cv)/2){
        cat('\nCross-validation produced too many NAs!')
        cat('\nReturning arbitrarily high cost\n')
        return(runif(1, 1, 1.5))
      } else {
        cat('\nComputing cost despite', sum(is.na(k.cv$residual)), '/', nrow(k.cv), 'NAs')
      }
    }
  }
  k.cv <- k.cv %>% filter(!is.na(residual))
  # Calculating cost function after Li et al 2018
  # https://doi.org/10.1016/j.cageo.2018.07.011
  # Weights
  iwt <- interp.weight # interpolation error
  vwt <- vgrm.weight # variogram fit error
  # Calculate variogram fit cost
  #   equation: sqrt( vgrm.weight * RMSE / sd(experimental.vgrm) ) [mWm^-2]
  #   note: The weighted RMSE is minimized during the fitting process with weights
  #   defined in fit.variogram method, see table 4.2 in http://www.gstat.org/gstat.pdf
  vgrm.rmse <-
    sqrt(attr(fitted.vgrm, "SSErr") / nrow(experimental.vgrm))
  # Important! Take square root of variogram RMSE to get in units of mWm^-2
  vgrm.cost <-
    sqrt(vwt * vgrm.rmse / sd(experimental.vgrm$gamma, na.rm = T))
  # Calculate interpolation cost similarly to vgrm.cost
  #   equation: interp.weight * RMSE / sd(cross-validation estimates) [mWm^-2]
  interp.rmse <-
    sqrt(sum(k.cv$residual^2, na.rm = T) / nrow(k.cv))
  interp.cost <-
    iwt * interp.rmse / sd(k.cv$var1.pred, na.rm = T)
  # Return cost = vgrm.error + interp.error
  if(verbose) {
    cat('\nSegment:', segment)
    cat('\nCutoff proportion:', cutoff.prop)
    cat('\nNumber of lags:', n.lags)
    cat('\nLag start:', lag.start)
    cat('\nMax pairs:', n.max)
    cat('\nVariogram weight:', vwt)
    cat('\nVariogram rmse:', vgrm.rmse)
    cat('\nVariogram cost:', vgrm.cost)
    cat('\nInterpolation weight:', iwt)
    cat('\nInterpolation rmse:', interp.rmse)
    cat('\nInterpolation cost:', interp.cost)
    cat('\nCost:', vgrm.cost + interp.cost)
    cat('\nVariogram model:', model.vgrm)
    cat('\n')
    print(fitted.vgrm)
    cat('\n', rep('+', 30), '\n', sep='')
  }
  return(vgrm.cost + interp.cost)
}

# Splits lines longer than a given threshold into the minimum number
# of pieces to all be under the given threshold
# From: https://gist.github.com/dblodgett-usgs/cf87392c02d73f1b7d16153d2b66a8f3
split_lines <- function(input_lines, max_length, id = NULL) {
  geom_column <- attr(input_lines, "sf_column")

  input_crs <- sf::st_crs(input_lines)

  input_lines[["geom_len"]] <- sf::st_length(input_lines[[geom_column]])

  attr(input_lines[["geom_len"]], "units") <- NULL
  input_lines[["geom_len"]] <- as.numeric(input_lines[["geom_len"]])

  too_long <- 
    input_lines %>%
    select(all_of(id), all_of(geom_column), geom_len) %>%
    filter(geom_len >= max_length)

  rm(input_lines) # just to control memory usage in case this is big.

  too_long <-
    too_long %>%
    mutate(
      pieces = ceiling(geom_len / max_length),
      piece_len = (geom_len / pieces),
      fID = 1:nrow(too_long)
    )

  split_points <-
    sf::st_set_geometry(too_long, NULL)[rep(seq_len(nrow(too_long)), too_long[["pieces"]]),]

  split_points <-
    split_points %>%
    mutate(split_fID = row.names(split_points)) %>%
    select(-geom_len, -pieces) %>%
    group_by(fID) %>%
    mutate(ideal_len = cumsum(piece_len)) %>%
    ungroup()

  coords <- data.frame(sf::st_coordinates(too_long[[geom_column]]))
  rm(too_long)

  coords <- rename(coords, fID = L1) %>% mutate(nID = 1:nrow(coords))

  split_nodes <-
    coords %>%
    group_by(fID) %>%
    # First calculate cumulative length by feature.
    mutate(len  = sqrt(((X - (lag(X)))^2) + (((Y - (lag(Y)))^2)))) %>%
    mutate(len = ifelse(is.na(len), 0, len)) %>%
    mutate(len = cumsum(len)) %>%
    # Now join nodes to split points -- this generates all combinations.
    left_join(select(split_points, fID, ideal_len, split_fID), by = "fID") %>%
    # Calculate the difference between node-wise distance and split-point distance.
    mutate(diff_len = abs(len - ideal_len)) %>%
    # regroup by the new split features.
    group_by(split_fID) %>%
    # filter out na then grab the min distance
    filter(!is.na(diff_len) & diff_len == min(diff_len)) %>%
    ungroup() %>%
    # Grab the start node for each geometry -- the end node of the geometry before it.
    mutate(start_nID = lag(nID),
           # need to move the start node one for new features.
           new_feature = fID - lag(fID, default = -1),
           start_nID = ifelse(new_feature == 1, start_nID + 1, start_nID)) %>%
    # Clean up the mess
    select(fID, split_fID, start_nID, stop_nID = nID, -diff_len, -ideal_len, -len, -X, -Y)

  split_nodes$start_nID[1] <- 1

  split_points <- 
    split_points %>%
    left_join(select(split_nodes, split_fID, start_nID, stop_nID), by = "split_fID")

  new_line <- function(start_stop, coords) {
    sf::st_linestring(as.matrix(coords[start_stop[1]:start_stop[2], c("X", "Y")]))
  }

  split_lines <- apply(as.matrix(split_points[c("start_nID", "stop_nID")]),
                        MARGIN = 1, FUN = new_line, coords = coords)

  split_lines <- st_sf(split_points[c(id, "split_fID")], geometry = st_sfc(split_lines, crs = input_crs))

  return(split_lines)
}

# Splits SZ segment buffers into equal parts and crops data points and
# interpolations by intersection for each buffer segment
split_segment <-
  function(
    seg.name,
    buf.dir = 'l',
    seg.num = 6,
    buf.len = 5e5,
    sector.exclude = NULL
  ) {
  pnts <- shp.hf.crop[[seg.name]]
  seg <- shp.segs[[seg.name]]
  buf <- st_buffer(seg, dist = buf.len, endCapStyle = 'ROUND')
  volc <- shp.volc
  split.seg <- split_lines(seg, as.numeric(st_length(seg)/seg.num))
  split.buf <-
    map(1:seg.num, ~
      st_buffer(
        split.seg[.x,],
        dist = ifelse(buf.dir == 'l', -1*buf.len, 1*buf.len),
        singleSide = T
      ) %>%
      st_intersection(buf)
    )
  pnts.buf <-
    map(1:seg.num, ~
      st_intersection(pnts, split.buf[[.x]]) %>%
      mutate(
        distance.from.seg = as.vector(st_distance(split.seg[.x,]$geometry, geometry)),
        .before = geometry
      )
    )
  volc.buf <-
    map(1:seg.num, ~
      st_intersection(volc, split.buf[[.x]]) %>%
      mutate(
        distance.from.seg = as.vector(st_distance(split.seg[.x,]$geometry, geometry)),
        .before = geometry
      )
    )
  best.mod <-
    solns %>%
    filter(segment == seg.name & v.mod != 'Gau') %>%
    slice_min(cost)
  interp <- best.mod[['shp.interp.diff']][[1]]
  interp.buf <-
    map(1:seg.num, ~
      st_intersection(interp, split.buf[[.x]]) %>%
      mutate(
        distance.from.seg = as.vector(st_distance(split.seg[.x,]$geometry, geometry)),
        .before = geometry
      )
    )
  if(!is.null(sector.exclude)) {
    split.buf <- split.buf[-sector.exclude]
    pnts.buf <- pnts.buf[-sector.exclude]
    volc.buf <- volc.buf[-sector.exclude]
    interp.buf <- interp.buf[-sector.exclude]
  }
  if(0 %in% map_dbl(pnts.buf, nrow)) {
    split.buf <- split.buf[!(map_dbl(pnts.buf, nrow) %in% 0)]
    pnts.buf <- pnts.buf[!(map_dbl(pnts.buf, nrow) %in% 0)]
    volc.buf <- volc.buf[-sector.exclude]
    interp.buf <- interp.buf[!(map_dbl(pnts.buf, nrow) %in% 0)]
  }
  return(
    list(
      'seg' = split.seg,
      'buf' = split.buf,
      'pnts' = pnts.buf,
      'volc' = volc.buf,
      'interp' = interp.buf
    )
  )
}

# Plot split segment
plot_split_segment <-
  function(
    split.seg,
    running.avg = 3,
    borders =
      c(
        'left' = 0.1,
        'right' = 0.1,
        'top' = 0.1,
        'bottom' = 0.1
      )
  ) {
  seg.name <- split.seg$interp[[1]]$segment[1]
  seg.num <- as.numeric(unique(bind_rows(split.seg$pnts)$split_fID))
  bx <- bbox_widen(st_bbox(st_buffer(st_combine(split.seg$seg), dist = 5e5)), borders = borders)
  world <- shp.world %>% st_crop(bx)
  ridge <- shp.ridge.crop[[seg.name]] %>% st_crop(bx)
  trench <- shp.trench.crop[[seg.name]] %>% st_crop(bx)
  transform <- shp.transform.crop[[seg.name]] %>% st_crop(bx)
  volc <- split.seg$volc
  wdth <- range(st_bbox(bind_rows(split.seg$buf))[c('xmin', 'xmax')])/1e3
  y.lim <- 
    c(
      median(bind_rows(split.seg$interp)$est.sim)
      - 2*IQR(bind_rows(split.seg$interp)$est.sim),
      median(bind_rows(split.seg$interp)$est.sim)
      + 2*IQR(bind_rows(split.seg$interp)$est.sim),
      median(bind_rows(split.seg$interp)$est.krige)
      - 2*IQR(bind_rows(split.seg$interp)$est.krige),
      median(bind_rows(split.seg$interp)$est.krige)
      + 2*IQR(bind_rows(split.seg$interp)$est.krige)
    )
  if(range(y.lim)[1] < 0) {
    y.lim[y.lim <= 0] <- 0
  }
  med.sim <- map_dbl(split.seg$interp, ~median(.x$est.sim))
  range.med.sim <- range(med.sim)
  med.krige <- map_dbl(split.seg$interp, ~median(.x$est.krige))
  range.med.krige <- range(med.krige)
  p0 <-
    ggplot() +
      geom_sf(data = world, size = 0.1, fill = 'grey60') +
      geom_sf(data = ridge, size = 0.5, alpha = 0.8) +
      geom_sf(data = trench, size = 0.5, alpha = 0.8) +
      geom_sf(data = transform, size = 0.5, alpha = 0.8) +
      geom_sf(data = split.seg$seg, size = 2) +
      geom_sf(
        data = bind_rows(split.seg$pnts),
        aes(fill = split_fID, group = split_fID),
        shape = 22,
        size = 1,
        show.legend = F
      ) +
      geom_sf(
        data = bind_rows(volc),
        color = 'gold',
        shape = 18
      ) +
      geom_sf(
        data = bind_rows(split.seg$buf),
        aes(color = split_fID),
        size = 1,
        fill = NA,
        show.legend = F
      ) +
      theme_map(font_size = 10) +
      theme(
        plot.tag = element_text(face = 'bold', size = 14),
        axis.text = element_text(),
        axis.text.x = element_text(angle = 30),
        panel.grid = element_line(size = 0.1, color = 'white'),
        panel.background = element_rect(fill = 'grey50', color = NA),
        plot.margin = margin()
      )
  p1 <-
    bind_rows(split.seg$interp) %>%
    st_set_geometry(NULL) %>%
    select(est.sim, est.krige, split_fID) %>%
    rename(Similarity = est.sim, Krige = est.krige) %>%
    pivot_longer(-split_fID) %>%
    group_by(name) %>%
    ggplot() +
    annotate(
      'rect',
      xmin = range.med.krige[1],
      xmax = range.med.krige[2],
      ymin = -Inf,
      ymax = Inf,
      fill = 'grey20',
      color = 'black',
      alpha = 0.8
    ) +
    annotate(
      'rect',
      xmin = range.med.sim[1],
      xmax = range.med.sim[2],
      ymin = -Inf,
      ymax = Inf,
      fill = 'ivory',
      color = 'black',
      alpha = 0.8
    ) +
    stat_density_ridges(
      aes(x = value, y = split_fID, fill = name, linetype = name),
      quantile_lines = T,
      quantiles = 2,
      rel_min_height = 0.03
    ) +
    coord_cartesian(xlim = range(y.lim)) +
    scale_fill_manual(values = c('grey20', 'ivory')) +
    labs(x = bquote('Heat Flow'~(mWm^-2)), y = 'Sector', linetype = NULL, fill = NULL) +
    theme_classic(base_size = 10) +
    theme(
      legend.box.margin = margin(),
      legend.background = element_rect(fill = NA, color = NA),
      legend.key.size = unit(0.8, 'lines'),
      panel.background = element_rect(fill = 'grey50'),
      axis.title.x = element_text(vjust = 5)
    )
  p2 <-
    bind_rows(split.seg$interp) %>%
    st_set_geometry(NULL) %>%
    filter(distance.from.seg <= 500000) %>%
    select(est.sim, est.krige, split_fID, distance.from.seg) %>%
    arrange(distance.from.seg) %>%
    mutate(
      Similarity = zoo::rollmean(est.sim, running.avg, fill = NA),
      Krige = zoo::rollmean(est.krige, running.avg, fill = NA)
    ) %>%
    select(-c(est.sim, est.krige)) %>%
    pivot_longer(-c(split_fID, distance.from.seg)) %>%
    group_by(name) %>%
    ggplot() +
      geom_point(
        data = bind_rows(volc),
        aes(distance.from.seg/1e3, range(y.lim)[1]),
        color = 'gold',
        shape = 18
      ) +
      geom_point(
        data = bind_rows(split.seg$pnts),
        aes(distance.from.seg/1e3, hf, fill = split_fID, group = split_fID),
        size = 1.5,
        shape = 22,
        alpha = 0.8
      ) +
      geom_point(
        aes(distance.from.seg/1e3, value, color = split_fID, group = split_fID),
        size = 0.3,
        shape = 3,
        alpha = 0.8,
        show.legend = F
      ) +
      geom_smooth(
        aes(distance.from.seg/1e3, value, color = split_fID, group = split_fID),
        method = 'loess',
        span = 0.95,
        size = 1,
        se = F,
        show.legend = F
      ) +
      labs(
        x = 'Distance from Trench (km)',
        y = bquote('Heat Flow'~(mWm^-2)),
        fill = 'Sector'
      ) +
      guides(
        fill = guide_legend(
          nrow = 1,
          override.aes = list(color = NA, size = 5, alpha = 1)
        )
      ) +
      scale_fill_discrete_qualitative(palette = 'Dark 3') +
      coord_cartesian(ylim = range(y.lim)) +
      facet_wrap(~name) +
      theme_classic(base_size = 10) +
      theme(
        panel.background = element_rect(fill = 'grey50'),
        legend.key.size = unit(0.8, 'lines'),
        legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.dir = 'horizontal',
        legend.background = element_rect(fill = NA)
      )
  pp1 <- p0 + p1 + plot_layout(widths = 1)
  p <-
    pp1 / p2 +
    plot_layout(widths = c(1, 2), guides = 'collect') +
    plot_annotation(
      tag_level = 'a',
      caption = seg.name
    ) &
    theme(
      plot.margin = margin(2, 2, 2, 2),
      plot.tag = element_text(face = 'bold', vjust = 1),
      plot.title = element_text(hjust = 0.5, vjust = 0),
      legend.position = 'bottom'
    )
  p
}

# Plot variogram
plot_vgrm <-
  function(
    experimental.vgrm,
    fitted.vgrm = NULL,
    cost = NULL,
    v.mod = NULL,
    lineCol = 'deeppink'
  ){
  # Check for missing arguments
  if(is.null(experimental.vgrm)) stop('\nMissing experimental variogram!')
  p <- 
    ggplot() +
    labs(
      x = 'Lag Distance (km)',
      y = bquote('Semivariance'~(mWm^-2))
    ) +
    theme_classic() +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  plt <- tryCatch(
    {
      p +
      geom_line(
        data = variogramLine(fitted.vgrm, maxdist = max(experimental.vgrm$dist)),
        aes(x = dist/1000, y = gamma),
        color = lineCol
      ) +
      geom_point(
        data = experimental.vgrm,
        aes(x = dist/1000, y = gamma),
        size = 0.8,
        shape = 19
      ) +
      annotate(
        'label',
        x = -Inf,
        y = Inf,
        label = paste('Model:', v.mod, '\nCost:', round(cost, 3)),
        label.padding = unit(0.1, 'lines'),
        label.r = unit(0, 'lines'),
        alpha = 0.8,
        fill = 'grey80',
        size = 3,
        vjust = 1.1,
        hjust = -0.1
      )
    },
    error=function(cond) {
      # If there was an error with the variogram fit ...
      # plot only the experimental variogram
      cat('Somthing is up with variogram model!')
      cat('\nPlotting experimental variogram only\n')
      exp.plt <- p +
      geom_point(
        data = experimental.vgrm,
        aes(x = dist/1000, y = gamma),
        size = 0.8,
        shape = 19
      )
      return(exp.plt)
    }
  )
  return(plt)
}
