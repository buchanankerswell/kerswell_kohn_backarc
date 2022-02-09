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
  'colorspace',
  'metR',
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

# Calculate Similarity rmse
itp_rmse <- function(seg.name, itp = NULL, type = 'sim'){
  if(is.null(itp)) {
    stop('need interpolation sf object')
  }
  buf <- shp.buffer[[seg.name]]
  grd <- shp.grid.crop[[seg.name]]
  obs <- shp.hf.crop[[seg.name]]
  if(type == 'sim') {
    itp <- suppressWarnings(st_intersection(itp, buf))
    # Find nearest Similarity estimate for each observation
    nearest.est.sim <- itp$est.sim[st_nearest_feature(obs, grd)]
    return(sqrt(mean((nearest.est.sim - obs$hf)^2)))
  } else if(type == 'krg') {
    itp <- itp
    # Find nearest Kriging estimate for each observation
    nearest.est.krige <- itp$est.krige[st_nearest_feature(obs, grd)]
    return(sqrt(mean((nearest.est.krige - obs$hf)^2)))
  } else {
    stop('invalid type!')
  }
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
    vwt * sqrt(vgrm.rmse / sd(experimental.vgrm$gamma, na.rm = T))
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
    mutate(
      start_nID = lag(nID),
      # need to move the start node one for new features.
      new_feature = fID - lag(fID, default = -1),
      start_nID = ifelse(new_feature == 1, start_nID + 1, start_nID)
    ) %>%
    # Clean up the mess
    select(fID, split_fID, start_nID, stop_nID = nID, -diff_len, -ideal_len, -len, -X, -Y)
  split_nodes$start_nID[1] <- 1
  split_points <- 
    split_points %>%
    left_join(select(split_nodes, split_fID, start_nID, stop_nID), by = "split_fID")
  new_line <- function(start_stop, coords) {
    sf::st_linestring(as.matrix(coords[start_stop[1]:start_stop[2], c("X", "Y")]))
  }
  split_lines <-
    apply(
      as.matrix(split_points[c("start_nID", "stop_nID")]),
      MARGIN = 1,
      FUN = new_line,
      coords = coords
    )
  split_lines <-
    st_sf(split_points[c(id, "split_fID")], geometry = st_sfc(split_lines, crs = input_crs))
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
    filter(segment == seg.name) %>%
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
  seg <- bind_rows(split.seg$seg)
  world <- shp.world
  ridge <- shp.ridge.crop[[seg.name]]
  trench <- shp.trench.crop[[seg.name]]
  transform <- shp.transform.crop[[seg.name]]
  volc <- split.seg$volc
  wdth <- range(st_bbox(bind_rows(split.seg$buf))[c('xmin', 'xmax')])/1e3
  rng <-
    bind_rows(split.seg$interp) %>%
    st_set_geometry(NULL) %>%
    select(est.sim, est.krige, split_fID) %>%
    summarise(
      min.krige = min(est.krige),
      max.krige = max(est.krige),
      min.sim = min(est.sim),
      max.sim = max(est.sim)
    )

  p0 <-
    ggplot() +
    geom_sf(data = bx, color = NA, fill = NA) +
    geom_sf(data = world, size = 0.1, fill = 'grey70', color = 'black') +
    geom_sf(data = ridge, size = 0.5, alpha = 0.8) +
    geom_sf(data = trench, size = 0.5, alpha = 0.8) +
    geom_sf(data = transform, size = 0.5, alpha = 0.8) +
    geom_sf(data = seg, size = 2) +
    geom_sf(
      data = bind_rows(split.seg$buf),
      aes(
        fill = factor(split_fID, levels = seg.num[order(seg.num)]),
        color = factor(split_fID, levels = seg.num[order(seg.num)])
      ),
      size = 0.5,
      alpha = 0.5,
      show.legend = F
    ) +
    geom_sf(
      data = bind_rows(split.seg$pnts),
      aes(
        fill = factor(split_fID, levels = seg.num[order(seg.num)]),
        group = factor(split_fID, levels = seg.num[order(seg.num)])
      ),
      size = 0.7,
      shape = 22,
      show.legend = F
    ) +
    geom_sf(
      data = bind_rows(volc),
      color = 'gold',
      shape = 18,
      size = 0.8
    ) +
    annotate(
      'label',
      label = 'b',
      x = -Inf,
      y = Inf,
      size = 5,
      hjust = 0,
      vjust = 1,
      color = 'grey40',
      fill = 'white',
      label.padding = unit(0.02, 'in'),
      label.r = unit(0, 'in')
    ) +
    coord_sf(
      xlim = c(st_bbox(bx)$xmin, st_bbox(bx)$xmax),
      ylim = c(st_bbox(bx)$ymin, st_bbox(bx)$ymax),
      label_axes = 'EN--'
      ) +
    scale_color_discrete_qualitative('Dark 3') +
    scale_fill_discrete_qualitative('Dark 3') +
    guides(
      fill = 
        guide_legend(
          nrow = 1,
          override.aes = list(alpha = 1, size = 8),
          title.position = 'top',
          title.vjust = 1,
          label.position = 'bottom'
        )
    ) +
    theme_map(font_size = 9) +
    theme(
      axis.text.x = element_text(color = 'grey40', angle = 30, hjust = 0, vjust = 0),
      axis.text.y = element_text(color = 'grey40', angle = 30, hjust = 0),
      panel.grid = element_line(size = 0.1, color = 'white'),
      panel.background = element_rect(fill = 'grey50', color = NA)
    )
    if(
       seg.name %in% c('Alaska Aleutians', 'Kamchatka Marianas', 'Tonga New Zealand', 'Vanuatu')
     ) {
      p0 <-
        p0 +
        scale_x_continuous(
          breaks = c(
            80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, -170, -160, -150, -140, -130
          )
        )
    }
  p1 <-
    bind_rows(split.seg$interp) %>%
    st_set_geometry(NULL) %>%
    select(est.sim, est.krige, split_fID) %>%
    rename(Similarity = est.sim, Kriging = est.krige) %>%
    pivot_longer(-split_fID) %>%
    group_by(name) %>%
    filter(split_fID %in% seg.num) %>%
    ggplot() +
    geom_boxplot(
      aes(
        y = value,
        x = factor(split_fID, levels = seg.num[order(seg.num)]),
        fill = name
      ),
      color = 'black',
      outlier.size = 0.1,
      outlier.color = 'grey30',
      size = 0.5
    ) +
    annotate(
      'label',
      label = 'a',
      x = -Inf,
      y = Inf,
      size = 5,
      hjust = 0,
      vjust = 1,
      color = 'grey40',
      fill = 'white',
      label.padding = unit(0.02, 'in'),
      label.r = unit(0, 'in')
    ) +
    guides(
      fill =
        guide_legend(
          nrow = 1,
          override.aes = list(alpha = 1),
          title.position = 'top',
          title.vjust = 1,
          label.position = 'bottom'
        )
    ) +
    coord_cartesian(
      ylim = c(min(rng$min.krige, rng$min.sim), max(rng$max.krige, rng$max.sim))
    ) +
    scale_x_discrete(position = 'top') +
    scale_fill_manual(values = c('ivory', 'cornflowerblue')) +
    labs(y = bquote(mWm^-2), x = NULL, fill = 'method') +
    theme_dark(base_size = 12) +
    theme(
      axis.title.x = element_text(margin = margin(t = -10)),
      legend.title = element_text(margin = margin(0, 0, -5, 0)),
      legend.key.height = unit(0.31, 'in')
    )
  p2 <-
    bind_rows(split.seg$interp) %>%
    st_set_geometry(NULL) %>%
    filter(distance.from.seg <= 500000) %>%
    select(est.sim, est.krige, split_fID, distance.from.seg) %>%
    arrange(distance.from.seg) %>%
    mutate(
      Similarity = zoo::rollmean(est.sim, running.avg, fill = NA),
      Kriging = zoo::rollmean(est.krige, running.avg, fill = NA)
    ) %>%
    select(-c(est.sim, est.krige)) %>%
    pivot_longer(-c(split_fID, distance.from.seg)) %>%
    group_by(name) %>%
    ggplot() +
    geom_point(
      data = bind_rows(volc),
      aes(distance.from.seg/1e3, min(rng$min.krige, rng$min.sim)),
      color = 'gold',
      shape = 18
    ) +
    geom_point(
      data = bind_rows(split.seg$pnts),
      aes(
        distance.from.seg/1e3,
        hf,
        fill = factor(split_fID, levels = seg.num[order(seg.num)]),
        group = factor(split_fID, levels = seg.num[order(seg.num)])
      ),
      size = 0.8,
      shape = 22,
      alpha = 0.3
    ) +
    geom_smooth(
      aes(
        distance.from.seg/1e3,
        value,
        color = factor(split_fID, levels = seg.num[order(seg.num)]),
        group = factor(split_fID, levels = seg.num[order(seg.num)])
      ),
      fill = 'ivory',
      method = 'loess',
      span = 0.95,
      alpha = 0.1,
      size = 1,
      se = T,
      show.legend = F
    ) +
    annotate(
      'label',
      label = 'c',
      x = -Inf,
      y = Inf,
      size = 5,
      hjust = 0,
      vjust = 1,
      color = 'grey40',
      fill = 'white',
      label.padding = unit(0.02, 'in'),
      label.r = unit(0, 'in')
    ) +
    labs(
      x = 'kilometers from trench',
      y = bquote(mWm^-2),
      fill = 'sector'
    ) +
    guides(
      fill = 
        guide_legend(
          nrow = 1,
          override.aes = list(alpha = 1, size = 8),
          title.position = 'top',
          title.vjust = 1,
          label.position = 'bottom'
        )
    ) +
    coord_cartesian(
      ylim = c(min(rng$min.krige, rng$min.sim), max(rng$max.krige, rng$max.sim))
    ) +
    scale_color_discrete_qualitative('Dark 3') +
    scale_fill_discrete_qualitative('Dark 3') +
    facet_wrap(~name, ncol = 2) +
    theme_dark(base_size = 12) +
    theme(
      legend.title = element_text(margin = margin(0, 0, -5, 0)),
      strip.background = element_rect(color = 'grey50', fill = 'grey95'),
      strip.text = element_text(color = 'black', size = 12, margin = margin(1, 0, 1.5, 0)),
      legend.spacing.x = unit(0, 'mm')
    )
  p <-
    ((p1 | p0) / p2) +
    plot_layout(guides = 'collect') &
    theme(
      plot.margin = margin(1, 1, 1, 1),
      legend.box.margin = margin(t = -5),
      legend.margin = margin(),
      legend.position = 'bottom',
      legend.direction = 'horizontal',
      legend.justification = 'left'
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
    lineCol = 'deeppink',
    ylim = NULL,
    xlim = NULL
  ){
  # Check for missing arguments
  if(is.null(experimental.vgrm)) stop('\nMissing experimental variogram!')
  p <- 
    ggplot() +
    labs(x = 'kilometer x 100', y = 'semivariance') +
    coord_cartesian(ylim = ylim, xlim = xlim) +
    theme_dark(base_size = 12) +
    theme(
      plot.margin = margin(1, 1, 1, 1)
    )
  plt <- tryCatch(
    {
      p +
      geom_point(
        data = experimental.vgrm,
        aes(x = dist/1e5, y = gamma),
        shape = 20
      ) +
      geom_line(
        data = variogramLine(fitted.vgrm, maxdist = max(experimental.vgrm$dist)),
        aes(x = dist/1e5, y = gamma),
        color = lineCol
      ) +
      annotate(
        'label',
        x = -Inf,
        y = Inf,
        label = paste('model:', v.mod, '\ncost:', round(cost, 4)),
        label.padding = unit(0.1, 'lines'),
        label.r = unit(0, 'lines'),
        alpha = 0.8,
        fill = 'grey80',
        size = 3,
        vjust = 1,
        hjust = 0
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
        aes(x = dist/1e5, y = gamma),
        shape = 20
      )
      return(exp.plt)
    }
  )
  return(plt)
}
