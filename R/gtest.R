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
# Tidy data
g.full <-
  g %>%
  mutate(
    `Age ocean (999 on continents)` =
      ifelse(`Age ocean (999 on continents)` > 900, NA, `Age ocean (999 on continents)`)
  ) %>%
  rename(`Age ocean` = `Age ocean (999 on continents)`)
g.colocated <- g.full %>% filter(!is.na(`Observed heat flow`))
# Make into sf object
shp.g.full <-
  st_as_sf(g.full, coords = c(1,2), crs = proj4.wgs) %>%
  st_transform(proj4.rp)
shp.g.colocated <-
  st_as_sf(g.colocated, coords = c(1,2), crs = proj4.wgs) %>%
  st_transform(proj4.rp)
# Compute global point density
dns <-
  MASS::kde2d(
  map_dbl(shp.g.colocated$geometry, ~.x[1]),
  map_dbl(shp.g.colocated$geometry, ~.x[2]),
    n = 100,
    lims = c(
      c(st_bbox(shp.world)$xmin, st_bbox(shp.world)$xmax),
      c(st_bbox(shp.world)$ymin, st_bbox(shp.world)$ymax)
    )
  )
dns.colocated <-
  expand.grid(dns$x, dns$y) %>%
  as_tibble() %>%
  rename(x = Var1, y = Var2) %>%
  mutate(
    dens = as.vector(dns$z),
    cnt = nrow(shp.g.colocated) / sum(dens) * dens
  )
# Units for goutorbe2011 dataset
unts <- c(
  'microwatts per cubic meter', 'microwatts per cubic meter', 'kilometer',
  'kilometer', 'grams per cubic centimeter', 'meter', 'mega annum', 'kilometer',
  'kilometer/kilometer', 'kilometer', 'mega annum', 'mega annum', '', '', '',
  'kilometer', 'kilometer', 'kilometer', 'kilometer', 'kilometer',
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
    data = shp.g.colocated,
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
    color = 'grey40',
    fill = 'white',
    label.padding = unit(0.02, 'in'),
    label.r = unit(0, 'in')
  ) +
  ggtitle('Global coverage of heat flow observations') +
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
    data = dns.colocated,
    aes(x, y, z = cnt)
  ) +
  geom_contour2(
    data = dns.colocated,
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
        barwidth = unit(2, 'in')
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
    color = 'grey40',
    fill = 'white',
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

# Plot colocated maps
plts1 <-
  map2(names(shp.g.full)[1:20], unts[1:20], ~{
    if(!(.x %in% c('geometry', 'segment'))) {
      shp.g.colocated <- shp.g.colocated[!is.na(shp.g.colocated[[.x]]),]
      if(!is.numeric(shp.g.colocated[[.x]])) {
        shp.g.colocated[[.x]] <- encode_ordinal(shp.g.colocated[[.x]])
      }
      p <-
        ggplot() +
        geom_sf(data = shp.world, size = 0.3, fill = 'grey70', color = 'black') +
        geom_sf(data = shp.ridge, size = 0.4, color = 'black', alpha = 0.8) +
        geom_sf(data = shp.trench, size = 0.4, color = 'black', alpha = 0.8) +
        geom_sf(data = shp.transform, size = 0.4, color = 'black', alpha = 0.8) +
        geom_sf(
          data = shp.g.colocated[.x],
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
              barwidth = unit(2, 'in')
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
          color = 'grey40',
          fill = 'white',
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
          color = 'grey40',
          fill = 'white',
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
              barwidth = unit(2, 'in')
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
      if(!is.numeric(shp.g.full[[.x]]) | .x %in% c('Up mantle velocity structure (class)')) {
        X <- encode_ordinal(g.full[[.x]][!is.na(g.full[[.x]])])
        Y <- encode_ordinal(g.colocated[[.x]][!is.na(g.colocated[[.x]])])
        d <-
          as_tibble(table(X)) %>%
          rename(type = X, full.count = n) %>%
          mutate(
            type = as.integer(type),
            colocated.count = as_tibble(table(Y))$n,
            diff.count = full.count - colocated.count,
            'undersampling' =
              ifelse(diff.count > 0, full.count-colocated.count, NA),
            'oversampling' =
              ifelse(diff.count < 0, full.count-colocated.count, NA)
          ) %>%
          rename(
            'full grid (Goutorbe et al., 2011)' = full.count,
            'colocated w/ hf observations' = colocated.count
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
            color = 'grey40',
            fill = 'white',
            label.padding = unit(0.02, 'in'),
            label.r = unit(0, 'in')
          ) +
          scale_y_continuous(breaks = c(0), position = 'right') +
          scale_fill_manual(
            values =
              c(
                'oversampling' = 'firebrick',
                'undersampling' = 'navy',
                'full grid (Goutorbe et al., 2011)' = 'grey80',
                'colocated w/ hf observations' = 'grey20'
              ),
            guide =
              guide_legend(
                override.aes = list(size = 2),
                title.position = 'top',
                title.vjust = 1,
                nrow = 2
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
        Y <- g.colocated[[.x]][!is.na(g.colocated[[.x]])]
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
        dns.colocated <-
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
            aes(x, dns, color = 'full grid (Goutorbe et al., 2011)'),
            size = 1
          ) +
          geom_path(
            data = dns.colocated,
            aes(x, dns, color = 'colocated w/ hf observations'),
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
            color = 'grey40',
            fill = 'white',
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
                nrow = 2
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
                nrow = 2
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
    p <-
      (..1 /
      (..2 + theme(plot.title = element_blank()))) |
      (..3 + theme(plot.title = element_blank())) +
      plot_layout(guides = 'collect') &
      theme(
        legend.position = 'bottom',
        legend.justification = 'left',
        legend.box.margin = margin(t = -10),
        plot.margin = margin(1, 1, 1, 1)
      )
    ggsave(
      paste0('figs/goutorbe2011_param/', str_replace_all(..4, ' ', '_'), '.png'),
      plot = p, type = 'cairo', width = 6.5, height = 4
    )
  }
})