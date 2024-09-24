#!/usr/bin/env Rscript

# Load packages and functions
source('R/functions.R')

# Define map projections
proj4.wgs <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
proj4.rp <- paste0('+proj=robin +lon_0=-155 +lon_wrap=-155 +x_0=0 +y_0=0 +ellps=WGS84 ',
                   '+datum=WGS84 +units=m +no_defs')

load('data/hf.RData')

# Create directory
dir.create('figs/global_density', recursive=T, showWarnings=F)

# Encoding function for categorical variables
encode_ordinal <- function(x, order=unique(x)) {
  x <- as.numeric(factor(x, levels=order, exclude=NULL))
  x
}

# Read goutorbe2011 data
g <- read_csv('data/g11-sim.csv', show_col_types=F)
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
  rename(`age ocean`=`age ocean (999 on continents)`)
g.colocated <- g.full %>% filter(!is.na(`observed heat flow`))

# Make into sf object
shp.g.full <-
  st_as_sf(g.full, coords=c(1,2), crs=proj4.wgs) %>%
  st_transform(proj4.rp)
shp.g.colocated <-
  st_as_sf(g.colocated, coords=c(1,2), crs=proj4.wgs) %>%
  st_transform(proj4.rp)

# Compute global point density
dns <-
  MASS::kde2d(
  map_dbl(shp.g.colocated$geometry, ~.x[1]),
  map_dbl(shp.g.colocated$geometry, ~.x[2]),
    n=100,
    lims=c(
      c(st_bbox(shp.world)$xmin, st_bbox(shp.world)$xmax),
      c(st_bbox(shp.world)$ymin, st_bbox(shp.world)$ymax)
    )
  )
dns.colocated <-
  expand.grid(dns$x, dns$y) %>%
  as_tibble() %>%
  rename(x=Var1, y=Var2) %>%
  mutate(
    dens=as.vector(dns$z),
    cnt=nrow(shp.g.colocated) / sum(dens) * dens
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
  ggtitle('a) Sediment thickness') +
  geom_sf(data=shp.relief.world, aes(color=elevation), shape=15, size=0.01) +
  scale_color_etopo(guide='none') +
  new_scale_color() +
  geom_sf(data=st_union((bind_rows(shp.buffer))), linewidth=0.3, fill=NA) +
  geom_sf(data=shp.ridge, linewidth=0.3) +
  geom_sf(data=shp.trench, linewidth=0.3) +
  geom_sf(data=shp.transform, linewidth=0.3) +
  geom_sf(
    data=shp.g.colocated,
    aes(color=`sediment thickness`),
    shape=15,
    size=0.3
  ) +
  geom_sf(data=st_union(bind_rows(shp.segs)), linewidth=1, color='white') +
  scale_color_viridis_c(
    option='viridis',
    name='km',
    limits=c(
      round(min(shp.g.colocated[['sediment thickness']])),
      round(max(shp.g.colocated[['sediment thickness']]))
    ),
    breaks=c(
      round(min(shp.g.colocated[['sediment thickness']])),
      round(max(shp.g.colocated[['sediment thickness']]))
    ),
    na.value='transparent',
    guide=guide_colorbar(title.vjust=1, show.limits=T)
  ) +
  coord_sf() +
  theme_map(font_size=14)
p2 <-
  ggplot() +
  ggtitle('b) Observational frequency') +
  geom_sf(data=shp.relief.world, aes(color=elevation), shape=15, size=0.01) +
  scale_color_etopo(guide='none') +
  new_scale_color() +
  geom_contour_fill(data=dns.colocated, aes(x, y, z=cnt), alpha=0.6) +
  geom_contour2(data=dns.colocated, aes(x, y, z=cnt), size=0.3) +
  scale_fill_viridis_c(
    option='viridis',
    name='frequency',
    limits=c(min(dns.colocated$cnt), max(dns.colocated$cnt)),
    breaks=c(min(dns.colocated$cnt), max(dns.colocated$cnt)),
    labels=c('low', 'high'),
    na.value='transparent',
    guide=guide_colorbar(title.vjust=1, show.limits=T)
  ) +
  geom_sf(data=st_union((bind_rows(shp.buffer))), linewidth=0.3, fill=NA) +
  geom_sf(data=shp.ridge, linewidth=0.3) +
  geom_sf(data=shp.trench, linewidth=0.3) +
  geom_sf(data=shp.transform, linewidth=0.3) +
  geom_sf(data=shp.g.colocated, shape=15, size=0.3) +
  geom_sf(data=st_union(bind_rows(shp.segs)), linewidth=1, color='white') +
  coord_sf() +
  theme_map(font_size=14)
p <-
  (p1 + theme(axis.text=element_blank())) / p2 &
  theme(
    plot.margin=margin(1, 1, 1, 1),
    legend.position='top',
    legend.justification='right',
    legend.direction='horizontal',
    axis.text=element_text(hjust=1),
    legend.margin=margin(-4, 0, -12, 0),
    legend.box.margin=margin(0, 10, 0, 0),
    legend.key.height=unit(0.125, 'in'),
    legend.key.width=unit(0.2, 'in'),
    legend.title=element_text(vjust=0, color='black', size=14),
    panel.grid=element_line(size=0.01, color='grey60'),
    plot.title=element_text(vjust=0, margin=margin(0, 0, -10, 0))
  )
cat('\nSaving plot to: figs/global_density/global-dens.png')
ggsave(
  file='figs/global_density/global-dens.png',
  plot=p,
  device='png',
  type='cairo',
  width=6.5,
  height=6.5,
  dpi=330
)

# Plot colocated maps
plts1 <-
  map2(names(shp.g.full)[1:20], unts[1:20], ~{
    if(!(.x %in% c('geometry', 'segment'))) {
      cat('\nPlotting colocated grids for ', .x, sep='')
      shp.g.colocated <- shp.g.colocated[!is.na(shp.g.colocated[[.x]]),]
      if(!is.numeric(shp.g.colocated[[.x]])) {
        shp.g.colocated[[.x]] <- encode_ordinal(shp.g.colocated[[.x]])
      }
      p <-
        ggplot() +
        ggtitle(paste0('a) ', .x)) +
        geom_sf(data=shp.relief.world, aes(color=elevation), shape=15, size=0.01) +
        scale_color_etopo(guide='none') +
        new_scale_color() +
        geom_sf(data=st_union((bind_rows(shp.buffer))), linewidth=0.3, fill=NA) +
        geom_sf(data=shp.ridge, linewidth=0.3) +
        geom_sf(data=shp.trench, linewidth=0.3) +
        geom_sf(data=shp.transform, linewidth=0.3) +
        geom_sf(
          data=shp.g.colocated[.x],
          aes(color=get(.x)),
          shape=15,
          size=0.3,
          show.legend=F
        ) +
        geom_sf(data=st_union(bind_rows(shp.segs)), linewidth=1, color='white') +
        scale_color_viridis_c(
          option='viridis',
          name=.y,
          limits=c(
            round(min(shp.g.colocated[[.x]])),
            round(max(shp.g.colocated[[.x]]))
          ),
          breaks=c(
            round(min(shp.g.colocated[[.x]])),
            round(max(shp.g.colocated[[.x]]))
          ),
          na.value='transparent',
          guide=guide_colorbar(title.vjust=1, show.limits=T)
        ) +
        labs(color=.y) +
        coord_sf() +
        theme_map(font_size=14) +
        theme(
          plot.margin=margin(1, 1, 1, 1),
          legend.position='top',
          legend.justification='right',
          legend.direction='horizontal',
          axis.text=element_text(hjust=1),
          legend.margin=margin(-4, 0, -12, 0),
          legend.box.margin=margin(0, 10, 0, 0),
          legend.key.height=unit(0.125, 'in'),
          legend.key.width=unit(0.2, 'in'),
          legend.title=element_text(vjust=0, color='black', size=14),
          panel.grid=element_line(size=0.01, color='grey60'),
          plot.title=element_text(vjust=0, margin=margin(0, 0, -10, 0))
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
      cat('\nPlotting full grids for ', .x, sep='')
      shp.g.full <- shp.g.full[!is.na(shp.g.full[[.x]]),]
      if(!is.numeric(shp.g.full[[.x]])) {
        shp.g.full[[.x]] <- encode_ordinal(shp.g.full[[.x]])
      }
      p <-
        ggplot() +
        ggtitle(paste0('b) ', .x)) +
        geom_sf(
          data=shp.g.full[.x],
          aes(color=get(.x)),
          shape=15,
          size=0.3
        ) +
        geom_sf(data=shp.relief.world, aes(color=elevation), shape=15, size=0.01) +
        scale_color_etopo(guide='none') +
        new_scale_color() +
        geom_sf(data=st_union((bind_rows(shp.buffer))), linewidth=0.3, fill=NA) +
        geom_sf(data=shp.ridge, linewidth=0.3) +
        geom_sf(data=shp.trench, linewidth=0.3) +
        geom_sf(data=shp.transform, linewidth=0.3) +
        geom_sf(data=st_union(bind_rows(shp.segs)), linewidth=1, color='white') +
        scale_color_viridis_c(
          option='viridis',
          name=.y,
          limits=c(
            round(min(shp.g.colocated[[.x]])),
            round(max(shp.g.colocated[[.x]]))
          ),
          breaks=c(
            round(min(shp.g.colocated[[.x]])),
            round(max(shp.g.colocated[[.x]]))
          ),
          na.value='transparent',
          guide=guide_colorbar(title.vjust=1, show.limits=T)
        ) +
        labs(color=.y) +
        coord_sf() +
        theme_map(font_size=14) +
        theme(
          plot.margin=margin(1, 1, 1, 1),
          legend.position='top',
          legend.justification='right',
          legend.direction='horizontal',
          axis.text=element_text(hjust=1),
          legend.margin=margin(-4, 0, -12, 0),
          legend.box.margin=margin(0, 10, 0, 0),
          legend.key.height=unit(0.125, 'in'),
          legend.key.width=unit(0.2, 'in'),
          legend.title=element_text(vjust=0, color='black', size=14),
          panel.grid=element_line(size=0.01, color='grey60'),
          plot.title=element_text(vjust=0, margin=margin(0, 0, -10, 0))
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
      cat('\nPlotting densities for ', .x, sep='')
      if(!is.numeric(shp.g.full[[.x]]) | .x %in% c('Up mantle velocity structure (class)')) {
        X <- encode_ordinal(g.full[[.x]][!is.na(g.full[[.x]])])
        Y <- encode_ordinal(g.colocated[[.x]][!is.na(g.colocated[[.x]])])
        d <-
          as_tibble(table(X)) %>%
          rename(type=X, full.count=n) %>%
          mutate(
            type=as.integer(type),
            full.count=full.count/length(X),
            colocated.count=as_tibble(table(Y))$n/length(X),
            diff.count=full.count - colocated.count,
            'undersampling' =
              ifelse(diff.count > 0, full.count-colocated.count, NA),
            'oversampling' =
              ifelse(diff.count < 0, full.count-colocated.count, NA)
          ) %>%
          rename(
            'global'=full.count,
            'hf obs'=colocated.count
          ) %>%
          select(-diff.count)
        p <-
          d %>%
          arrange(type) %>%
          mutate(type=factor(type, levels=type)) %>%
          pivot_longer(-type) %>%
          ggplot() +
          geom_col(
            aes(type, value, fill=name),
            color=NA,
            position='dodge'
          ) +
          annotate(
            'label',
            label='c',
            x=-Inf,
            y=Inf,
            size=5,
            hjust=0,
            vjust=1,
            fill='grey90',
            label.padding=unit(0.02, 'in'),
            label.r=unit(0, 'in')
          ) +
          scale_y_continuous(breaks=c(0), position='right') +
          scale_fill_manual(
            values =
              c(
                'oversampling'='firebrick',
                'undersampling'='navy',
                'global'='grey80',
                'hf obs'='grey20'
              ),
            guide =
              guide_legend(
                override.aes=list(size=2),
                title.position='top',
                title.vjust=1,
                nrow=1
              )
          ) +
          labs(x=.y, y='frequency', color=NULL, fill=NULL) +
          theme_dark(base_size=12) +
          theme(
            plot.margin=margin(1, 1, 1, 1),
            legend.box.margin=margin(),
            legend.margin=margin()
          )
      } else {
        X <- g.full[[.x]][!is.na(g.full[[.x]])]
        Y <- g.colocated[[.x]][!is.na(g.colocated[[.x]])]
        dns.full <-
          tibble(
            x =
              density(
                X,
                from=range(c(X, Y), na.rm=T)[1],
                to=range(c(X, Y), na.rm=T)[2],
                na.rm=T,
                n=1000
              )$x,
            dns =
              density(
                X,
                from=range(c(X, Y), na.rm=T)[1],
                to=range(c(X, Y), na.rm=T)[2],
                na.rm=T,
                n=1000
              )$y / length(X)
          )
        dns.colocated <-
          tibble(
            x =
              density(
                X,
                from=range(c(X, Y), na.rm=T)[1],
                to=range(c(X, Y), na.rm=T)[2],
                na.rm=T,
                n=1000
              )$x,
            dns =
              density(
                Y,
                from=range(c(X, Y), na.rm=T)[1],
                to=range(c(X, Y), na.rm=T)[2],
                na.rm=T,
                n=1000
              )$y / length(X)
          )
        dns.diff <-
          tibble(
            x =
              density(
                X,
                from=range(c(X, Y), na.rm=T)[1],
                to=range(c(X, Y), na.rm=T)[2],
                na.rm=T,
                n=1000
              )$x,
            dns =
              (density(
                X,
                from=range(c(X, Y), na.rm=T)[1],
                to=range(c(X, Y), na.rm=T)[2],
                na.rm=T,
                n=1000
              )$y / length(X)) -
              (density(
                Y,
                from=range(c(X, Y), na.rm=T)[1],
                to=range(c(X, Y), na.rm=T)[2],
                na.rm=T,
                n=1000
              )$y / length(X))
          )
        p <-
          ggplot() +
          geom_ribbon(
            data=filter(dns.diff, dns <= 0),
            aes(x, ymin=dns, ymax=0, fill='oversampling'),
            size=0.1
          ) +
          geom_ribbon(
            data=filter(dns.diff, dns >= 0),
            aes(x, ymin=0, ymax=dns, fill='undersampling'),
            size=0.1
          ) +
          geom_path(
            data=dns.full,
            aes(x, dns, color='global'),
            size=1
          ) +
          geom_path(
            data=dns.colocated,
            aes(x, dns, color='hf obs'),
            size=1
          ) +
          annotate(
            'label',
            label='c',
            x=-Inf,
            y=Inf,
            size=5,
            hjust=0,
            vjust=1,
            fill='grey90',
            label.padding=unit(0.02, 'in'),
            label.r=unit(0, 'in')
          ) +
          scale_y_continuous(breaks=c(0), position='right') +
          scale_color_grey(
            guide =
              guide_legend(
                override.aes=list(size=2),
                title.position='top',
                title.vjust=1,
                nrow=1
              )
          ) +
          scale_fill_manual(
            values=c(
              'oversampling'='firebrick',
              'undersampling'='navy'
            ),
            guide =
              guide_legend(
                override.aes=list(size=2),
                title.position='top',
                title.vjust=1,
                nrow=1
              )
          ) +
          labs(x=.y, y='density', color=NULL, fill=NULL) +
          theme_dark(base_size=12) +
          theme(
            plot.margin=margin(1, 1, 1, 1),
            legend.box.margin=margin(),
            legend.box='vertical',
            legend.margin=margin()
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
    cat('\nComposing plots for ', ..4, sep='')
    p <-
      (..1 /
      (..2 + theme(plot.title=element_blank()))) |
      (..3 + theme(plot.title=element_blank())) +
      plot_layout(guides='collect', widths=1) &
      theme(
        legend.position='bottom',
        legend.justification='left',
        legend.box.margin=margin(t=-10),
        plot.margin=margin(1, 1, 1, 1)
      )
      suppressWarnings(
        ggsave(
        paste0('figs/global_density/', str_replace_all(..4, ' ', '-'), '.png'),
        plot=p, type='cairo', width=8, height=4.9
      )
    )
  }
})

cat('\ngoutorbe-analysis.R complete!\n\n')