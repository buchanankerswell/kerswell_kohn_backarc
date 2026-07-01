#!/usr/bin/env Rscript

# Load packages and functions
source('R/functions.R')
load('assets/map_data/map-data.RData')

# Create directory
dir.create('figs/global_density', recursive=T, showWarnings=F)

# Encoding function for categorical variables
encode_ordinal <- function(x, order=unique(x)) {
  x <- as.numeric(factor(x, levels=order, exclude=NULL))
  x
}

# Read goutorbe2011 data
g <- read_csv('assets/hf_data/g11-sim.csv', show_col_types=F)
names(g) <- str_to_lower(names(g))

# Tidy data
g_full <-
  g %>%
  mutate(
    `age ocean (999 on continents)` =
      ifelse(`age ocean (999 on continents)` > 900, NA, `age ocean (999 on continents)`),
    `distance to young rift` =
      ifelse(
        `distance to young rift` < 0,
        abs(`distance to young rift`),
        `distance to young rift`
      ),
    `thermo-tectonic age` = encode_ordinal(`thermo-tectonic age`),
    `rift type` = encode_ordinal(`rift type`)
  ) %>%
  rename(`age ocean`=`age ocean (999 on continents)`)
g_colocated <- g_full %>% filter(!is.na(`observed heat flow`))

# Make into sf object
shp_g_full <-
  st_as_sf(g_full, coords=c(1,2), crs=wgs) %>%
  reproject_center_pacific()
shp_g_colocated <-
  st_as_sf(g_colocated, coords=c(1,2), crs=wgs) %>%
  reproject_center_pacific()

# Compute global point density
dns <-
  MASS::kde2d(
  map_dbl(shp_g_colocated$geometry, ~.x[1]),
  map_dbl(shp_g_colocated$geometry, ~.x[2]),
    n=100,
    lims=c(
      c(st_bbox(shp_world)$xmin, st_bbox(shp_world)$xmax),
      c(st_bbox(shp_world)$ymin, st_bbox(shp_world)$ymax)
    )
  )
density_colocated <-
  expand.grid(dns$x, dns$y) %>%
  as_tibble() %>%
  rename(x=Var1, y=Var2) %>%
  mutate(
    dens=as.vector(dns$z),
    cnt=nrow(shp_g_colocated) / sum(dens) * dens
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
  geom_sf(data=shp_relief_world, aes(color=elev), shape=15, size=0.01) +
  scale_color_etopo(guide='none') +
  new_scale_color() +
  geom_sf(data=shp_ridge, linewidth=0.3, color='white') +
  geom_sf(data=shp_transform, linewidth=0.3, color='white') +
  geom_sf(data=shp_trench, linewidth=0.3, color='white') +
  geom_sf(data=shp_g_colocated, aes(color=`sediment thickness`), size=0.1, shape=20) +
  geom_sf(data=st_union(bind_rows(shp_submap)), linewidth=1) +
  scale_color_viridis_c(option='viridis', name='km',
                        limits=c(
                          round(min(shp_g_colocated[['sediment thickness']])),
                          round(max(shp_g_colocated[['sediment thickness']]))
                        ),
                        breaks=c(
                          round(min(shp_g_colocated[['sediment thickness']])),
                          round(max(shp_g_colocated[['sediment thickness']]))
                        ),
                        na.value='transparent',
                        guide=guide_colorbar(title.vjust=1, show.limits=T,
                                             frame.colour='black',
                                             ticks.colour='black')) +
  ggtitle('Sediment Thickness') +
  coord_sf(expand=F, lims_method='geometry_bbox') +
  theme_map(font_size=14)
p2 <-
  ggplot() +
  geom_sf(data=shp_relief_world, aes(color=elev), shape=15, size=0.01) +
  scale_color_etopo(guide='none') +
  new_scale_color() +
  geom_contour_fill(data=density_colocated, aes(x, y, z=cnt), alpha=0.6) +
  geom_contour2(data=density_colocated, aes(x, y, z=cnt), size=0.3) +
  scale_fill_viridis_c(option='viridis', name='density',
                       limits=c(min(density_colocated$cnt), max(density_colocated$cnt)),
                       breaks=c(min(density_colocated$cnt), max(density_colocated$cnt)),
                       labels=c('low', 'high'),
                       na.value='transparent',
                       guide=guide_colorbar(title.vjust=1, show.limits=T,
                                            frame.colour='black',
                                            ticks.colour='black')) +
  geom_sf(data=shp_ridge, linewidth=0.3, color='white') +
  geom_sf(data=shp_transform, linewidth=0.3, color='white') +
  geom_sf(data=shp_trench, linewidth=0.3, color='white') +
  geom_sf(data=shp_g_colocated, size=0.1, shape=20) +
  geom_sf(data=st_union(bind_rows(shp_submap)), linewidth=1) +
  ggtitle('Observational Density') +
  coord_sf() +
  theme_map(font_size=14)
p <- p1 / p2 &
   theme(plot.margin=margin(), legend.position='bottom',
         legend.justification='center', legend.direction='horizontal',
         axis.text=element_blank(), legend.margin=margin(),
         legend.box.margin=margin(5, 5, 5, 5), legend.key.height=unit(0.5, 'cm'),
         legend.key.width=unit(1, 'cm'),
         legend.title=element_text(vjust=0, color='black', size=14),
         panel.grid=element_line(linewidth=0.05, color='grey20'),
         plot.title=element_text(vjust=0, hjust=0.5, margin=margin(10, 10, 10, 10)))
fig_path <- 'figs/global_density/global-dens.png'
cat('\nSaving plot to: ', fig_path)
ggsave(file=fig_path, plot=p, width=6.5, height=6.5, bg='white')

# Plot colocated maps
plts1 <-
  map2(names(shp_g_full)[1:20], unts[1:20], ~{
    if(!(.x %in% c('geometry', 'segment'))) {
      cat('\nPlotting colocated grids for ', .x, sep='')
      shp_g_colocated <- shp_g_colocated[!is.na(shp_g_colocated[[.x]]),]
      p <-
        ggplot() +
        ggtitle(paste0('a) ', .x)) +
        geom_sf(data=shp_relief_world, aes(color=elev), shape=15, size=0.01) +
        scale_color_etopo(guide='none') +
        new_scale_color() +
        geom_sf(data=shp_ridge, linewidth=0.3) +
        geom_sf(data=shp_trench, linewidth=0.3) +
        geom_sf(data=shp_transform, linewidth=0.3) +
        geom_sf(
          data=shp_g_colocated[.x],
          aes(color=get(.x)),
          shape=15,
          size=0.3,
          show.legend=F
        ) +
        geom_sf(data=st_union(bind_rows(shp_submap)), linewidth=1) +
        scale_color_viridis_c(
          option='viridis',
          name=.y,
          limits=c(
            round(min(shp_g_colocated[[.x]])),
            round(max(shp_g_colocated[[.x]]))
          ),
          breaks=c(
            round(min(shp_g_colocated[[.x]])),
            round(max(shp_g_colocated[[.x]]))
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
          panel.grid=element_line(linewidth=0.01, color='grey60'),
          plot.title=element_text(vjust=0, margin=margin(0, 0, -10, 0))
        )
    } else {
      p <- NULL
    }
    p
  })

# Plot full grid maps
plts2 <-
  map2(names(shp_g_full)[1:20], unts[1:20], ~{
    if(!(.x %in% c('geometry', 'segment'))) {
      cat('\nPlotting_full grids for ', .x, sep='')
      shp_g_full <- shp_g_full[!is.na(shp_g_full[[.x]]),]
      p <-
        ggplot() +
        ggtitle(paste0('b) ', .x)) +
        geom_sf(data=shp_g_full[.x], aes(color=get(.x)), shape=15, size=0.3) +
        geom_sf(data=shp_relief_world, aes(color=elev), shape=15, size=0.01) +
        scale_color_etopo(guide='none') +
        new_scale_color() +
        geom_sf(data=shp_ridge, linewidth=0.3) +
        geom_sf(data=shp_trench, linewidth=0.3) +
        geom_sf(data=shp_transform, linewidth=0.3) +
        geom_sf(data=st_union(bind_rows(shp_submap)), linewidth=1) +
        scale_color_viridis_c(
          option='viridis',
          name=.y,
          limits=c(
            round(min(shp_g_colocated[[.x]])),
            round(max(shp_g_colocated[[.x]]))
          ),
          breaks=c(
            round(min(shp_g_colocated[[.x]])),
            round(max(shp_g_colocated[[.x]]))
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
          panel.grid=element_line(linewidth=0.01, color='grey60'),
          plot.title=element_text(vjust=0, margin=margin(0, 0, -10, 0))
        )
    } else {
      p <- NULL
    }
    p
  })

# Plot parameter densities
plts3 <-
  map2(names(shp_g_full)[1:20], unts[1:20], ~{
    if(!(.x %in% c('geometry', 'segment'))) {
      cat('\nPlotting densities for ', .x, sep='')
      if(!is.numeric(shp_g_full[[.x]]) | .x %in% c('Up mantle velocity structure (class)')) {
        X <- g_full[[.x]]
        Y <- g_colocated[[.x]]
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
          annotate('label', label='c', x=-Inf, y=Inf, size=5, hjust=0, vjust=1, fill='grey90',
                   label.padding=unit(0.02, 'in'), label.r=unit(0, 'in')) +
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
        X <- g_full[[.x]][!is.na(g_full[[.x]])]
        Y <- g_colocated[[.x]][!is.na(g_colocated[[.x]])]
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
        density_colocated <-
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
            linewidth=0.1
          ) +
          geom_ribbon(
            data=filter(dns.diff, dns >= 0),
            aes(x, ymin=0, ymax=dns, fill='undersampling'),
            linewidth=0.1
          ) +
          geom_path(
            data=dns.full,
            aes(x, dns, color='global'),
            linewidth=1
          ) +
          geom_path(
            data=density_colocated,
            aes(x, dns, color='hf obs'),
            linewidth=1
          ) +
          annotate('label', label='c', x=-Inf, y=Inf, size=5, hjust=0, vjust=1, fill='grey90',
                   label.padding=unit(0.02, 'in'), label.r=unit(0, 'in')) +
          scale_y_continuous(breaks=c(0), position='right') +
          scale_color_grey(
            guide = guide_legend(override.aes=list(size=2), title.position='top',
                                 title.vjust=1, nrow=1)
          ) +
          scale_fill_manual(
            values=c(
              'oversampling'='firebrick',
              'undersampling'='navy'
            ),
            guide = guide_legend(override.aes=list(size=2), title.position='top',
                                 title.vjust=1, nrow=1)
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
pwalk(list(plts1, plts2, plts3, names(shp_g_full)[1:20]), ~{
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