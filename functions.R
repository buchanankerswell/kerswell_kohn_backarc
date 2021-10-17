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
  'stringr',
  'tidyr',
  'readr',
  'purrr',
  'tibble',
  'dplyr',
  'magrittr',
  'ggsflabel',
  'ggplot2',
  'ggrepel',
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

# Don't allow sf to use google's s2 library
# for spherical geometry
# This reverts to using the GEOS library instead
# which is what sf used before 1.0 release
sf_use_s2(FALSE)

# Draw a widened box from a st_bbox object
bbox_widen <-
  function(
    bbox,
    crs,
    borders =
      c(
        'left' = 0,
        'right' = 0,
        'top' = 0,
        'bottom' = 0
      )
  ) {
  b <- bbox # current bounding box
  xrange <- b$xmax - b$xmin # range of x values
  yrange <- b$ymax - b$ymin # range of y values
  b[1] <- b[1] - (borders['left'] * xrange) # xmin - left
  b[3] <- b[3] + (borders['right'] * xrange) # xmax - right
  b[2] <- b[2] - (borders['bottom'] * yrange) # ymin - bottom
  b[4] <- b[4] + (borders['top'] * yrange) # ymax - top
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
read_latlong <- function(file, crs) {
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
    shp.hf,
    cutoff.prop = 3,
    n.lags = 20,
    lag.start = 1
  ) {
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
  v <-
    variogram(
      hf~1,
      locations = shp.hf,
      cutoff = shifted.cutoff,
      width = bin.width
    )
  # Return shifted variogram
  return(v[lag.start:nrow(v),])
}

# Wrapper function around gstat krige method
Krige <- function(shp.hf, fitted.vgrm, shp.interp.grid, nmax=200) {
  # Print grid and parameters info
  cat(
    '\nKriging ...',
    '\nObservations:', nrow(shp.hf),
    '\nGrid size:', length(shp.interp.grid),
    '\nMax nearest points:', nmax,
    '\nVariogram params:\n'
  )
  print(fitted.vgrm)
  # Kriging
  suppressMessages({
    krige(
      formula = hf~1,
      locations = shp.hf,
      newdata = shp.interp.grid,
      model = fitted.vgrm,
      nmax = nmax,
      debug.level = -1
    ) %>%
    as_tibble() %>%
    st_as_sf() %>%
    rename(
      est = var1.pred,
      variance = var1.var
    ) %>%
    mutate(sigma = sqrt(variance), .before = geometry)
  })
}

# Take difference
interp_diff <- function(shp.interp.krige, shp.interp.compare) {
  suppressWarnings({
    # Crop Similarity estimates to Kriging bounding box
    shp.interp.compare %>%
    select(-obs) %>%
    st_intersection(shp.interp.krige) %>%
    # Rename
    rename(
      est.similarity = est,
      sigma.similarity = sigma
    ) %>%
    # Add Kriging estimates and compute difference
    mutate(
      est.krige = shp.interp.krige$est,
      sigma.krige = shp.interp.krige$sigma,
      est.diff = est.similarity - est.krige,
      .before = geometry
    )
  })
}

# Cost function for Kriging estimates
#   note: after Li et al 2018
#   see https://doi.org/10.1016/j.cageo.2018.07.011
cost.function <- 
  function(
    shp.hf,
    cutoff.prop = 3,
    n.lags = 20,
    lag.start = 1,
    model.vgrm = 0,
    nmax = 200,
    nfold = 10,
    verbose = F
  ) {
  # Calculate experimental variogram
  experimental.vgrm <-
    experimental_vgrm(
      shp.hf = shp.hf,
      cutoff.prop = cutoff.prop,
      n.lags = n.lags,
      lag.start = lag.start
    )
  # Fit experimental variogram with a model variogram
  fitted.vgrm <- NULL
  k.cv <- NA
  itr <- 1
  while(is.null(fitted.vgrm) | sum(is.na(k.cv)) != 0) {
    if(itr > 10) {
      cat('\nVariogram model:', v.mod)
      cat('\nFitted variogram model:\n')
      print(fitted.vgrm)
      cat('\nCross-validation NAs:', sum(is.na(k.cv)), '\n')
      stop('Failed to compute cost after ', itr-1, ' iterations')
    }
    itr <- itr + 1
    ft <- try(
      fit.variogram(
        object = experimental.vgrm,
        fit.method = 7,
        debug.level = 0,
        model = vgm(NA, v.mod, NA, NA)
      ),
      silent = T
    )
    if(any(class(ft) == 'try-error')){
      if(verbose){
        cat('\nCould not fit variogram:', v.mod)
      }
      v.mod <- sample(c('Sph', 'Gau', 'Exp'), 1)
      if(verbose){
        cat('\nTrying new variogram model:', v.mod)
      }
    } else {
      if(verbose){
        cat('\nFitted variogram model:\n')
        print(ft)
      }
      fitted.vgrm <- ft
      # Kriging with n-fold cross validation
      if(verbose){
        cat('Computing cross-validation\n')
      }
      k.cv <- try({
        krige.cv(
          formula = hf~1,
          locations = shp.hf,
          model = fitted.vgrm,
          nmax = nmax,
          nfold = nfold,
          verbose = ifelse(verbose, T, F)
        )
      })
      if(any(class(k.cv) == 'try-error')) {
        k.cv <- NA
        if(verbose){
          cat('\nCould not fit variogram:', v.mod)
        }
        v.mod <- sample(c('Sph', 'Gau', 'Exp'), 1)
        if(verbose){
          cat('\nTrying new variogram model:', v.mod)
        }
      }
    }
  }
  # Calculating cost function after Li et al 2018
  # https://doi.org/10.1016/j.cageo.2018.07.011
  # Weights
  interp.weight <- 0.6 # interpolation error
  vgrm.weight <- 0.4 # variogram fit error
  # Calculate variogram fit error
  #   equation: vgrm.weight * RMSE / sd(experimental.vgrm) [mWm^-2]
  #   note: The weighted RMSE is minimized during the fitting process with weights
  #   defined in fit.variogram method, see table 4.2 in http://www.gstat.org/gstat.pdf
  vgrm.error <-
    vgrm.weight *
    sqrt(attr(fitted.vgrm, "SSErr") / nrow(experimental.vgrm)) /
    sd(sqrt(experimental.vgrm$gamma))
  if(verbose){
    cat('\nVariogram error:', vgrm.error)
  }
  # Calculate interpolation error similarly to vgrm.error
  #   equation: interp.weight * RMSE / sd(cross-validation estimates)
  interp.error <-
    interp.weight *
    sqrt(sum(k.cv$residual^2) / nrow(k.cv)) /
    sd(k.cv$var1.pred)
  if(verbose){
    cat('\nInterpolation error:', interp.error)
  }
  # Return cost = vgrm.error + interp.error
  if(verbose){
    cat('\nCost:', vgrm.error + interp.error, '\n')
  }
  return(vgrm.error + interp.error)
}

# Plot variogram
plot_vgrm <- function(experimental.vgrm, fitted.vgrm, seg.name){
  p <- 
    ggplot() +
    theme_classic() +
    labs(x = NULL, title = seg.name) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 9),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    )
  tc <- tryCatch(
    {
      p +
      geom_line(
        data = variogramLine(fitted.vgrm, maxdist = max(experimental.vgrm$dist)),
        aes(x = dist/1000, y = gamma),
        color = 'deeppink'
      ) +
      geom_point(data = experimental.vgrm, aes(x = dist/1000, y = gamma))
    },
    error=function(cond) {
      # If there was an error with the variogram fit ...
      # plot only the experimental variogram
      message(paste('Somthing is up with', seg.name, 'variogram model:'))
      message(cond)
      message('\nPlotting experimental variogram without fitted variogram model ...\n')

      plt <- p +
      geom_point(data = experimental.vgrm, aes(x = dist/1000, y = gamma))

      return(plt)
    }
  )
  return(tc)
}

