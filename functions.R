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
      select(-obs.sim) %>%
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
        fill = 'grey90',
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
