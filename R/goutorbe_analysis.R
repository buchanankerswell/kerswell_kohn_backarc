#!/usr/bin/env Rscript

# Load packages and functions
source('R/functions.R')
load('data/hf.RData')
dir.create('figs/goutorbe2011_param', showWarnings = F)

# Encoding function for categorical variables
encode_ordinal <- function(x, order = unique(x)) {
  x <- as.numeric(factor(x, levels = order, exclude = NULL))
  x
}

# Map projections
# WGS84
proj4.wgs <-
  '+proj=longlat +lon_wrap=180 +ellps=WGS84 +datum=WGS84 +no_defs'
# Robinson Pacific centered
proj4.rp <-
  '+proj=robin +lon_0=-155 +lon_wrap=-155 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'

# Read goutorbe2011 data
g <- read_csv('data/goutorbe2011.csv', show_col_types = F)
names(g) <- str_to_lower(names(g))

# Tidy data
g.full <-
  g %>%
  mutate(
    `age ocean (999 on continents)` =
      ifelse(`age ocean (999 on continents)` > 900, NA, `age ocean (999 on continents)`),
    `distance to young rift` =
      ifelse(
        `distance to young rift` < 0,
        abs(`distance to young rift`),
        `distance to young rift`
      )
  ) %>%
  rename(`age ocean` = `age ocean (999 on continents)`)
g.collocated <- g.full %>% filter(!is.na(`observed heat flow`))

# Make into sf object
shp.g.full <-
  st_as_sf(g.full, coords = c(1,2), crs = proj4.wgs) %>%
  st_transform(proj4.rp)
shp.g.collocated <-
  st_as_sf(g.collocated, coords = c(1,2), crs = proj4.wgs) %>%
  st_transform(proj4.rp)

# Compute global point density
dns <-
  MASS::kde2d(
  map_dbl(shp.g.collocated$geometry, ~.x[1]),
  map_dbl(shp.g.collocated$geometry, ~.x[2]),
    n = 100,
    lims = c(
      c(st_bbox(shp.world)$xmin, st_bbox(shp.world)$xmax),
      c(st_bbox(shp.world)$ymin, st_bbox(shp.world)$ymax)
    )
  )
dns.collocated <-
  expand.grid(dns$x, dns$y) %>%
  as_tibble() %>%
  rename(x = Var1, y = Var2) %>%
  mutate(
    dens = as.vector(dns$z),
    cnt = nrow(shp.g.collocated) / sum(dens) * dens
  )

# Units for goutorbe2011 dataset
unts <- c(
  'microwatts per cubic meter', 'microwatts per cubic meter', 'kilometer',
  'kilometer', 'grams per cubic centimeter', 'meter', 'mega annum', 'kilometer',
  'kilometer/kilometer', 'kilometer', 'mega annum', 'mega annum', '', '', '',
  'kilometers from ridge', 'kilometer', 'kilometer', 'kilometers from young rift', 'kilometer',
  'milliwatts per square meter', 'milliwatts per square meter',
  'milliwatts per square meter', 'milliwatts per square meter',
  'milliwatts per square meter', ''
)

# Plot global density
p1 <-
  ggplot() +
  geom_sf(
    data = shp.world,
    size = 0.3,
    fill = 'grey70',
    color = 'black'
  ) +
  geom_sf(
    data = st_union((bind_rows(shp.buffer))),
    size = 0.3, fill = NA, color = rgb(1, 1, 1, 0.5)
  ) +
  geom_sf(
    data = shp.g.collocated,
    shape = 15,
    size = 0.3
  ) +
  geom_sf(data = shp.ridge, size = 0.4, color = 'black', alpha = 0.8) +
  geom_sf(data = shp.trench, size = 0.4, color = 'black', alpha = 0.8) +
  geom_sf(data = shp.transform, size = 0.4, color = 'black', alpha = 0.8) +
  geom_sf(data = st_union(bind_rows(shp.segs)), size = 0.8, color = 'white') +
  annotate(
    'label',
    label = 'a',
    x = -Inf,
    y = Inf,
    size = 5,
    hjust = 0,
    vjust = 1,
    fill = 'grey90',
    label.padding = unit(0.02, 'in'),
    label.r = unit(0, 'in')
  ) +
  coord_sf(expand = F) +
  theme_map(font_size = 12) +
  theme(
    panel.grid = element_line(size = 0.1, color = 'white'),
    panel.background = element_rect(fill = 'grey50', color = NA),
    plot.margin = margin(1, 1, 1, 1),
    legend.position = 'bottom',
    legend.justification = 'left',
    legend.spacing.x = unit(0, 'in'),
    legend.box.margin = margin(),
    legend.margin = margin()
  )
p2 <-
  ggplot() +
  geom_contour_fill(
    data = dns.collocated,
    aes(x, y, z = cnt)
  ) +
  geom_contour2(
    data = dns.collocated,
    aes(x, y, z = cnt),
    size = 0.3,
    color = rgb(1, 1, 1, 0.1)
  ) +
  scale_fill_continuous_sequential(
    'Turku',
    breaks = MakeBreaks(bins = 7),
    name = 'observation count density',
    guide =
      guide_colorstrip(
        inside = T,
        title.position = 'top',
        title.vjust = 1,
        barwidth = unit(2.5, 'in')
      )
  ) +
  geom_sf(
    data = shp.world,
    size = 0.3,
    fill = rgb(0, 0, 0, 0.3),
    color = rgb(1, 1, 1, 0.1)
  ) +
  geom_sf(
    data = shp.g.full,
    shape = 15,
    size = 0.3,
    alpha = 0.1
  ) +
  geom_sf(
    data = st_union((bind_rows(shp.buffer))),
    size = 0.3, fill = NA, color = rgb(1, 1, 1, 0.5)
  ) +
  geom_sf(data = shp.ridge, size = 0.4, color = 'black', alpha = 0.8) +
  geom_sf(data = shp.trench, size = 0.4, color = 'black', alpha = 0.8) +
  geom_sf(data = shp.transform, size = 0.4, color = 'black', alpha = 0.8) +
  geom_sf(data = st_union(bind_rows(shp.segs)), size = 0.8, color = 'white') +
  annotate(
    'label',
    label = 'b',
    x = -Inf,
    y = Inf,
    size = 5,
    hjust = 0,
    vjust = 1,
    fill = 'grey90',
    label.padding = unit(0.02, 'in'),
    label.r = unit(0, 'in')
  ) +
  coord_sf(expand = F) +
  theme_map(font_size = 12) +
  theme(
    panel.grid = element_line(size = 0.1, color = 'white'),
    panel.background = element_rect(fill = 'grey50', color = NA),
    plot.margin = margin(1, 1, 1, 1),
    legend.position = 'bottom',
    legend.justification = 'left',
    legend.spacing.x = unit(0, 'in'),
    legend.box.margin = margin(),
    legend.margin = margin()
  )
p <- p1 / p2
ggsave(
  'figs/goutorbe2011_param/global_dens.png', plot = p, type = 'cairo',
  width = 6, height = 9, dpi = 330
)

# Plot collocated maps
plts1 <-
  map2(names(shp.g.full)[1:20], unts[1:20], ~{
    if(!(.x %in% c('geometry', 'segment'))) {
      cat('\nPlotting collocated grids for ', .x, sep = '')
      shp.g.collocated <- shp.g.collocated[!is.na(shp.g.collocated[[.x]]),]
      if(!is.numeric(shp.g.collocated[[.x]])) {
        shp.g.collocated[[.x]] <- encode_ordinal(shp.g.collocated[[.x]])
      }
      p <-
        ggplot() +
        geom_sf(data = shp.world, size = 0.3, fill = 'grey70', color = 'black') +
        geom_sf(data = shp.ridge, size = 0.4, color = 'black', alpha = 0.8) +
        geom_sf(data = shp.trench, size = 0.4, color = 'black', alpha = 0.8) +
        geom_sf(data = shp.transform, size = 0.4, color = 'black', alpha = 0.8) +
        geom_sf(
          data = shp.g.collocated[.x],
          aes(color = get(.x)),
          shape = 15,
          size = 0.3,
          show.legend = F
        ) +
        geom_sf(data = st_union(bind_rows(shp.segs)), size = 0.8, color = 'white') +
        geom_sf(
          data = st_union((bind_rows(shp.buffer))),
          size = 0.3, fill = NA, color = rgb(1, 1, 1, 0.5)
        ) +
        scale_color_continuous_sequential(
          'Turku',
          breaks = MakeBreaks(bins = 7),
          na.value = 'transparent',
          guide =
            guide_colorstrip(
              inside = T,
              title.position = 'top',
              title.vjust = 1,
              barwidth = unit(3, 'in')
            )
        ) +
        labs(color = .y) +
        annotate(
          'label',
          label = 'a',
          x = -Inf,
          y = Inf,
          size = 5,
          hjust = 0,
          vjust = 1,
          fill = 'grey90',
          label.padding = unit(0.02, 'in'),
          label.r = unit(0, 'in')
        ) +
        coord_sf(expand = F) +
        theme_map(font_size = 12) +
        theme(
          panel.grid = element_line(size = 0.1, color = 'white'),
          panel.background = element_rect(fill = 'grey50', color = NA),
          plot.margin = margin(1, 1, 1, 1)
        )
    } else {
      p <- NULL
    }
    p
  })

# Plot full grid maps
plts2 <-
  map2(names(shp.g.full)[1:20], unts[1:20], ~{
    if(!(.x %in% c('geometry', 'segment'))) {
      cat('\nPlotting full grids for ', .x, sep = '')
      shp.g.full <- shp.g.full[!is.na(shp.g.full[[.x]]),]
      if(!is.numeric(shp.g.full[[.x]])) {
        shp.g.full[[.x]] <- encode_ordinal(shp.g.full[[.x]])
      }
      p <-
        ggplot() +
        geom_sf(
          data = shp.g.full[.x],
          aes(color = get(.x)),
          shape = 15,
          size = 0.3
        ) +
        geom_sf(
          data = shp.world,
          size = 0.3,
          fill = NA,
          color = 'black'
        ) +
        geom_sf(
          data = st_union((bind_rows(shp.buffer))),
          size = 0.3, fill = NA, color = rgb(1, 1, 1, 0.5)
        ) +
        geom_sf(data = shp.ridge, size = 0.4, color = 'black', alpha = 0.8) +
        geom_sf(data = shp.trench, size = 0.4, color = 'black', alpha = 0.8) +
        geom_sf(data = shp.transform, size = 0.4, color = 'black', alpha = 0.8) +
        geom_sf(data = st_union(bind_rows(shp.segs)), size = 0.8, color = 'white') +
        annotate(
          'label',
          label = 'b',
          x = -Inf,
          y = Inf,
          size = 5,
          hjust = 0,
          vjust = 1,
          fill = 'grey90',
          label.padding = unit(0.02, 'in'),
          label.r = unit(0, 'in')
        ) +
        scale_color_continuous_sequential(
          'Turku',
          breaks = MakeBreaks(bins = 7),
          na.value = 'transparent',
          guide =
            guide_colorstrip(
              inside = T,
              title.position = 'top',
              title.vjust = 1,
              barwidth = unit(3, 'in')
            )
        ) +
        labs(color = .y) +
        coord_sf(expand = F) +
        theme_map(font_size = 12) +
        theme(
          panel.grid = element_line(size = 0.1, color = 'white'),
          panel.background = element_rect(fill = 'grey50', color = NA),
          plot.margin = margin(1, 1, 1, 1),
          legend.position = 'bottom',
          legend.justification = 'left',
          legend.box.margin = margin(t = -10)
        )
    } else {
      p <- NULL
    }
    p
  })

# Plot parameter densities
plts3 <-
  map2(names(shp.g.full)[1:20], unts[1:20], ~{
    if(!(.x %in% c('geometry', 'segment'))) {
      cat('\nPlotting densities for ', .x, sep = '')
      if(!is.numeric(shp.g.full[[.x]]) | .x %in% c('Up mantle velocity structure (class)')) {
        X <- encode_ordinal(g.full[[.x]][!is.na(g.full[[.x]])])
        Y <- encode_ordinal(g.collocated[[.x]][!is.na(g.collocated[[.x]])])
        d <-
          as_tibble(table(X)) %>%
          rename(type = X, full.count = n) %>%
          mutate(
            type = as.integer(type),
            full.count = full.count/length(X),
            collocated.count = as_tibble(table(Y))$n/length(X),
            diff.count = full.count - collocated.count,
            'undersampling' =
              ifelse(diff.count > 0, full.count-collocated.count, NA),
            'oversampling' =
              ifelse(diff.count < 0, full.count-collocated.count, NA)
          ) %>%
          rename(
            'global' = full.count,
            'hf obs' = collocated.count
          ) %>%
          select(-diff.count)
        p <-
          d %>%
          arrange(type) %>%
          mutate(type = factor(type, levels = type)) %>%
          pivot_longer(-type) %>%
          ggplot() +
          geom_col(
            aes(type, value, fill = name),
            color = NA,
            position = 'dodge'
          ) +
          annotate(
            'label',
            label = 'c',
            x = -Inf,
            y = Inf,
            size = 5,
            hjust = 0,
            vjust = 1,
            fill = 'grey90',
            label.padding = unit(0.02, 'in'),
            label.r = unit(0, 'in')
          ) +
          scale_y_continuous(breaks = c(0), position = 'right') +
          scale_fill_manual(
            values =
              c(
                'oversampling' = 'firebrick',
                'undersampling' = 'navy',
                'global' = 'grey80',
                'hf obs' = 'grey20'
              ),
            guide =
              guide_legend(
                override.aes = list(size = 2),
                title.position = 'top',
                title.vjust = 1,
                nrow = 1
              )
          ) +
          labs(x = .y, y = 'frequency', color = NULL, fill = NULL) +
          theme_dark(base_size = 12) +
          theme(
            plot.margin = margin(1, 1, 1, 1),
            legend.box.margin = margin(),
            legend.margin = margin()
          )
      } else {
        X <- g.full[[.x]][!is.na(g.full[[.x]])]
        Y <- g.collocated[[.x]][!is.na(g.collocated[[.x]])]
        dns.full <-
          tibble(
            x =
              density(
                X,
                from = range(c(X, Y), na.rm = T)[1],
                to = range(c(X, Y), na.rm = T)[2],
                na.rm = T,
                n = 1000
              )$x,
            dns =
              density(
                X,
                from = range(c(X, Y), na.rm = T)[1],
                to = range(c(X, Y), na.rm = T)[2],
                na.rm = T,
                n = 1000
              )$y / length(X)
          )
        dns.collocated <-
          tibble(
            x =
              density(
                X,
                from = range(c(X, Y), na.rm = T)[1],
                to = range(c(X, Y), na.rm = T)[2],
                na.rm = T,
                n = 1000
              )$x,
            dns =
              density(
                Y,
                from = range(c(X, Y), na.rm = T)[1],
                to = range(c(X, Y), na.rm = T)[2],
                na.rm = T,
                n = 1000
              )$y / length(X)
          )
        dns.diff <-
          tibble(
            x =
              density(
                X,
                from = range(c(X, Y), na.rm = T)[1],
                to = range(c(X, Y), na.rm = T)[2],
                na.rm = T,
                n = 1000
              )$x,
            dns =
              (density(
                X,
                from = range(c(X, Y), na.rm = T)[1],
                to = range(c(X, Y), na.rm = T)[2],
                na.rm = T,
                n = 1000
              )$y / length(X)) -
              (density(
                Y,
                from = range(c(X, Y), na.rm = T)[1],
                to = range(c(X, Y), na.rm = T)[2],
                na.rm = T,
                n = 1000
              )$y / length(X))
          )
        p <-
          ggplot() +
          geom_ribbon(
            data = filter(dns.diff, dns <= 0),
            aes(x, ymin = dns, ymax = 0, fill = 'oversampling'),
            size = 0.1
          ) +
          geom_ribbon(
            data = filter(dns.diff, dns >= 0),
            aes(x, ymin = 0, ymax = dns, fill = 'undersampling'),
            size = 0.1
          ) +
          geom_path(
            data = dns.full,
            aes(x, dns, color = 'global'),
            size = 1
          ) +
          geom_path(
            data = dns.collocated,
            aes(x, dns, color = 'hf obs'),
            size = 1
          ) +
          annotate(
            'label',
            label = 'c',
            x = -Inf,
            y = Inf,
            size = 5,
            hjust = 0,
            vjust = 1,
            fill = 'grey90',
            label.padding = unit(0.02, 'in'),
            label.r = unit(0, 'in')
          ) +
          scale_y_continuous(breaks = c(0), position = 'right') +
          scale_color_grey(
            guide =
              guide_legend(
                override.aes = list(size = 2),
                title.position = 'top',
                title.vjust = 1,
                nrow = 1
              )
          ) +
          scale_fill_manual(
            values = c(
              'oversampling' = 'firebrick',
              'undersampling' = 'navy'
            ),
            guide =
              guide_legend(
                override.aes = list(size = 2),
                title.position = 'top',
                title.vjust = 1,
                nrow = 1
              )
          ) +
          labs(x = .y, y = 'density', color = NULL, fill = NULL) +
          theme_dark(base_size = 12) +
          theme(
            plot.margin = margin(1, 1, 1, 1),
            legend.box.margin = margin(),
            legend.box = 'vertical',
            legend.margin = margin()
          )
      }
    } else {
      p <- NULL
    }
    p
  })

# Save composite plots
pwalk(list(plts1, plts2, plts3, names(shp.g.full)[1:20]), ~{
  if(!(..4 %in% c('geometry', 'segment'))) {
    cat('\nComposing plots for ', ..4, sep = '')
    p <-
      (..1 /
      (..2 + theme(plot.title = element_blank()))) |
      (..3 + theme(plot.title = element_blank())) +
      plot_layout(guides = 'collect', widths = 1) &
      theme(
        legend.position = 'bottom',
        legend.justification = 'left',
        legend.box.margin = margin(t = -10),
        plot.margin = margin(1, 1, 1, 1)
      )
      suppressWarnings(
        ggsave(
        paste0('figs/goutorbe2011_param/', str_replace_all(..4, ' ', '_'), '.png'),
        plot = p, type = 'cairo', width = 8, height = 4.9
      )
    )
  }
})

cat('\n\nDone!\n\n')