# Script to compute Figure 5 from the manuscript, and Figures S14, S20 to S24 from the supplement
# - Bubble plot of cross-link fractions (CLFs) and node degrees (NDGs) 
#   of paleoclimate networks of ACER AP and simulated TSF, TS, PR of TraCE for the raw signal.
#   + Principal Component Analysis (PCA) on a map for the first principal component for ACER AP and TraCE TS
# - Analogous bubble plot (S20) and PCA maps (S14) for ACER and TraCE networks on millennial time scales
# - Analogous bubble plots for HadCM3 (S21), LOVECLIM (S22), and the emulated response of TraCE TSF (S24) and HadCM3 TSF (S23)
#   to TS, PR, CO2, and combined forcing.

#DONE & TESTED

## initialize data (skips automatically if done previously in the session, see `main_pcns.R`) ----
source('main_pcns.R')
source('main_pcn_meas.R')
source('main_pca.R')
source('main_emulation_pcn_meas.R')


## prepare maps for matrix axes ----
## read continent data and define labels and positions
cnts <- read_sf(dsn = file.path(DIR_DATASETS_SUPPLEMENTARY, 'esri_data/World_Continents/qgis_modified_shp/continents.shp'),
                layer = "continents")
cnt2plt <- list('Africa', 'Asia', c('Australia', 'Oceania'), 'North America', 'South America', 'Europe') # 'Antarctica'
cntcrs <- list(Africa = '+proj=longlat +datum=WGS84 +no_defs',
               Asia = '+proj=lcc +lat_1=30 +lat_2=62 +lat_0=0 +lon_0=105 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs',
               'Australia Oceania' = '+proj=lcc +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=134 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs ',
               'North America' = '+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
               'South America' = '+proj=lcc +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs',
               Europe = '+proj=lcc +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs')
cntlb <- list(Africa = 'AF',
              Asia = 'AS',
              'Australia Oceania' = 'OC',
              'North America' = 'NA',
              'South America' = 'SA',
              Europe = 'EU')

## plot function for the continent strip ----
plot_continents_row <- function(direction, label.position, list_only = F) {
  if (!label.position %in% c('aligned', 'jittered')) stop('label.position must be `aligned` or `jittered')
  if (!direction %in% c('vertical', 'horizontal')) stop('direction must be `vertical` or `horizontal')
  
  if (label.position == 'jittered') {
    cntx <- list(Africa = 0.3, Asia = 0.8,'Australia Oceania' = 0.75,'North America' = 0.22,'South America' = 0.7,Europe = 0.75)
    cnty <- list(Africa = 0.3, Asia = 0.42,'Australia Oceania' = 0.6,'North America' = 0.55,'South America' = 0.2,Europe = 0.15)
    xlim <- c(0,1); ylim <- c(0,1)
  } else if (label.position == 'aligned') {
    if (direction == 'horizontal') {
      cntx <- list(Africa = 0.5, Asia = 0.5,'Australia Oceania' = 0.5,'North America' = 0.5,'South America' = 0.5,Europe = 0.5)
      cnty <- list(Africa = 1.1, Asia = 1.1,'Australia Oceania' = 1.1,'North America' = 1.1,'South America' = 1.1,Europe = 1.1)
      xlim <- c(0,1); ylim <- c(0,1.3)
    } else if (direction == 'vertical') {
      cntx <- list(Africa = 1.1, Asia = 1.1,'Australia Oceania' = 1.1,'North America' = 1.1,'South America' = 1.1,Europe = 1.1)
      cnty <- list(Africa = 0.5, Asia = 0.5,'Australia Oceania' = 0.5,'North America' = 0.5,'South America' = 0.5,Europe = 0.5)
      xlim <- c(0,1.3); ylim <- c(0,1)
    }
  }
  cntsp <- list(Africa = 1, Asia = 1.25,'Australia Oceania' = 1,'North America' = 1.3,'South America' = 1.1,Europe = 1.15)
  
  cntplt <- lapply(1:length(cnt2plt), function(i) {
    cntnm <- cnt2plt[[i]]
    cntnm2 <- paste(cntnm, collapse = ' ')
    cnt <- cnts %>%
      filter(CONTINENT %in% cntnm) %>% 
      as(., 'Spatial') %>% 
      gSimplify(., cntsp[[cntnm2]]) %>% 
      spTransform(., CRS(cntcrs[[cntnm2]])) %>% 
      as(., 'sf')
    plt <- ggplot(cnt) + 
      geom_sf(fill = 'grey85', color = 'grey30', size = 0.5) + 
      coord_sf() + 
      theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
            plot.background = element_blank(), panel.background = element_blank())
    plt <- ggdraw(plt, xlim = xlim, ylim = ylim) + 
      annotate('label', x = cntx[[cntnm2]], y = cnty[[cntnm2]], label = cntlb[[cntnm2]],
               size = GLOBAL_FONT_SIZE - 8, fontface = GLOBAL_FONT_FACE_TITLE, 
               label.size = 1)
  }) %>% 
    setNames(lapply(cnt2plt, function(cntnm) {paste(cntnm, collapse = ' ')}))
  
  if (list_only) {
    return(cntplt) 
  } else {
    if (direction == 'horizontal') {
      cntlw <- list(Africa = 0.2, Asia = 0.2,'Australia Oceania' = 0.2,Europe = 0.2,'North America' = 0.2,'South America' = 0.2)
      cntlh <- list(Africa = 0.2, Asia = 0.2,'Australia Oceania' = 0.2,Europe = 0.2,'North America' = 0.2,'South America' = 0.2)
      cntlx <- list(Africa = 0.05, Asia = 0.21,'Australia Oceania' = 0.365,Europe = 0.52,'North America' = 0.68,'South America' = 0.84)
      cntly <- list(Africa = 0.1, Asia = 0.1,'Australia Oceania' = 0.1,Europe = 0.1,'North America' = 0.1,'South America' = 0.1)
    } else if (direction == 'vertical') {
      cntlw <- list(Africa = 0.23, Asia = 0.23,'Australia Oceania' = 0.23,Europe = 0.23,'North America' = 0.23,'South America' = 0.23)
      cntlh <- list(Africa = 0.18, Asia = 0.18,'Australia Oceania' = 0.18,Europe = 0.18,'North America' = 0.18,'South America' = 0.18)
      cntly <- list(Africa = 0.05, Asia = 0.21,'Australia Oceania' = 0.3675,Europe = 0.525,'North America' = 0.68,'South America' = 0.835)
      cntlx <- list(Africa = 0.1, Asia = 0.1,'Australia Oceania' = 0.1,Europe = 0.1,'North America' = 0.1,'South America' = 0.1)
    }
  }
  
  plt <- ggdraw()
  for (i in 1:length(cntplt)) {
    nm <- names(cntplt)[i]
    plt <- plt + 
      draw_plot(cntplt[[nm]], x  = cntlx[[nm]], y = cntly[[nm]], width = cntlw[[nm]], height = cntlh[[nm]])
  }
  return(plt)
}

## plot functions for the bubble matrixs----
### 2 by 2 (first row raw signal, second row millennial-scale signal)
plot_nw_clf_and_ndg_matr2x2_nosig <- function(datclf_raw, datclf_mil, datndg_raw, datndg_mil, suffix_pairs,
                                              nw_meas_clf = 'pcross_link',
                                              nw_meas_ndg = 'node_degree',
                                              regions = REGIONS,
                                              legend = TRUE, leg_option = 'share', no_borders = FALSE) {
  
  plts <- lapply(c('raw', 'mil'), function(tscale) {
    datclf <- get(paste0('datclf_',tscale))
    datndg <- get(paste0('datndg_',tscale))
    plts <- lapply(1:length(suffix_pairs[[tscale]]), function(i){
      suffixes <- suffix_pairs[[tscale]][[i]]
      print(suffixes)
      # prepare clf data
      datclf <- datclf %>% 
        filter(cross_reg1 %in% unique(regions$region) & cross_reg2 %in% unique(regions$region)) %>% 
        mutate(cross_reg1 = substr(cross_reg1,1,2),cross_reg2 = substr(cross_reg2,1,2))
      datclf <- datclf %>% select(window, cross_reg1, cross_reg2, 
                                  !!!syms(mapply(function(nm, pr) return(paste0(pr, nm)), pr = c(paste0(nw_meas_clf, c('.p','.n'), '.')), nm = sort(rep(x = suffixes, 2)), USE.NAMES = FALSE)))
      datclf <- lapply(c(paste0(nw_meas_clf, c('.p','.n'))), function(c) {
        return(datclf %>% select(window, cross_reg1, cross_reg2, one_of(paste0(c, '.', suffixes))) %>% gather(key = 'type', value = !!sym(paste0(c, '_val')), -window, -cross_reg1, -cross_reg2) %>%
                 rowwise() %>% mutate(type = str_extract(type, suffixes) %>% .[!is.na(.)]) %>% ungroup())
      }) %>% plyr::join_all(., by = c('window', 'cross_reg1', 'cross_reg2', 'type'), type = 'inner') %>% as_tibble() %>% mutate_at(vars(window), as.character) %>% 
        rowwise %>% 
        rename(p = paste0(nw_meas_clf, '.p_val'), n = paste0(nw_meas_clf, '.n_val')) %>% 
        mutate(type = str_replace(type, '\\.', ' ')) %>% ungroup()
      datclf$type <- factor(datclf$type, levels = str_replace(suffixes, '\\.', ' '))
      
      regs_abv <- tibble(reg = distinct(bind_rows(select(datclf,cross_reg1),select(datclf,cross_reg2) %>% rename(cross_reg1 = cross_reg2)))$cross_reg1) %>% 
        rownames_to_column(var = 'reg_n') %>% 
        mutate_at(vars(reg_n), as.numeric)
      
      scaler_clf <- 0.5
      diag_start <- 0.5
      diag_end <- length(unique(c(as.vector(datclf$cross_reg1), as.vector(datclf$cross_reg2)))) + 0.5
      
      datclf <- datclf %>% 
        mutate_at(vars(cross_reg1,cross_reg2), as.character) %>% 
        inner_join(regs_abv %>% rename(cross_reg1 = reg, reg1_n = reg_n)) %>% 
        inner_join(regs_abv %>% rename(cross_reg2 = reg, reg2_n = reg_n)) %>% 
        mutate(sum = p+n, sum_scaled = sum*scaler_clf) %>% # mapping to srqt(sum) done in stat = 'pie' (consider legend below)
        filter(sum != 0) %>% 
        select(type,reg1_n,reg2_n,p,n,sum,sum_scaled)
      
      datclf <- datclf %>% 
        group_by(type) %>% 
        mutate(reg1_n_sort = case_when(type == suffixes[1] & reg1_n > reg2_n ~ reg2_n,
                                       type == suffixes[1] & reg1_n < reg2_n ~ reg1_n,
                                       type == suffixes[2] & reg1_n > reg2_n ~ reg1_n,
                                       type == suffixes[2] & reg1_n < reg2_n ~ reg2_n,
                                       TRUE ~ reg1_n),
               reg2_n_sort = case_when(type == suffixes[1] & reg1_n > reg2_n ~ reg1_n,
                                       type == suffixes[1] & reg1_n < reg2_n ~ reg2_n,
                                       type == suffixes[2] & reg1_n > reg2_n ~ reg2_n,
                                       type == suffixes[2] & reg1_n < reg2_n ~ reg1_n,
                                       TRUE ~ reg2_n)) %>% 
        rename('+' = p, '-' = n)
      
      # prepare ndg data
      datndg <- datndg %>% select(window, site_id,
                                  !!!syms(mapply(function(nm, pr) return(paste0(pr, nm)), pr = c(paste0(nw_meas_ndg, c('','.p','.n'), '.')), nm = sort(rep(x = suffixes, 3)), USE.NAMES = FALSE)))
      datndg <- lapply(c(paste0(nw_meas_ndg, c('', '.p','.n'))), function(c) {
        return(datndg %>% select(window, site_id, one_of(paste0(c, '.', suffixes))) %>% gather(key = 'type', value = !!sym(paste0(c)), -window, -site_id) %>%
                 rowwise() %>% mutate(type = str_extract(type, suffixes) %>% .[!is.na(.)]) %>% ungroup())
      }) %>% plyr::join_all(., by = c('window', 'site_id', 'type'), type = 'inner') %>% as_tibble() %>% mutate_at(vars(window), as.character) %>% 
        rowwise %>%
        gather(key = 'ndg_type', value = 'ndg_val', !!!syms(paste0(paste0(nw_meas_ndg, c('', '.p','.n'))))) %>% 
        rowwise %>%
        mutate(ndg_type = case_when(str_detect(ndg_type, '.p') ~ 'p', str_detect(ndg_type, '.n') ~ 'n', TRUE ~ 'sum')) %>% 
        #rename(p = paste0(nw_meas_ndg, '.p_val'), n = paste0(nw_meas_ndg, '.n_val')) %>% 
        filter(type %in% suffixes)
      
      scaler_ndg <- 0.05
      
      datndg <- datndg %>% 
        inner_join(regions, by = 'site_id') %>% 
        mutate(region = substr(region,1,2)) %>% 
        spread(ndg_type,ndg_val) %>% 
        group_by(region, type, window) %>% 
        summarise(sum = mean(sum, na.rm = T), p = mean(p, na.rm = T), n = mean(n, na.rm = T)) %>% 
        ungroup() %>% 
        inner_join(regs_abv %>% rename(region = reg), by = 'region') %>% 
        mutate(sum_scaled = sqrt(sum)*scaler_ndg) 
      numndg <- datndg %>% 
        #filter(ndg_type == 'node_degree') %>% 
        group_by(region) %>% 
        count() %>% mutate_at(vars(n), as.character)
      datndg <- datndg %>% 
        mutate(frac_p = case_when(p != 0 & n != 0 ~ p/sum,
                                  p == 0 & n != 0 ~ 0,
                                  p != 0 & n == 0 ~ 1,
                                  T ~ 0.5),
               frac_n = case_when(p != 0 & n != 0 ~ n/sum,
                                  p == 0 & n != 0 ~ 1,
                                  p != 0 & n == 0 ~ 0,
                                  T ~ 0.5)) %>% 
        gather(key = 'sgn', value = 'single',p,n) %>% 
        mutate(start = case_when(type == suffixes[1] & sgn == 'p' ~ 1/4*pi - frac_p*pi,
                                 type == suffixes[2] & sgn == 'p' ~ 1/4*pi,
                                 type == suffixes[1] & sgn == 'n' ~ -3/4*pi,
                                 type == suffixes[2] & sgn == 'n' ~ 5/4*pi - frac_n*pi),
               end = case_when(type == suffixes[1] & sgn == 'p' ~ 1/4*pi,
                               type == suffixes[2] & sgn == 'p' ~ 1/4*pi + frac_p*pi,
                               type == suffixes[1] & sgn == 'n' ~ -3/4*pi + frac_n*pi,
                               type == suffixes[2] & sgn == 'n' ~ 5/4*pi)) %>% 
        mutate(type = str_replace_all(type, '_', ' ')) %>% ungroup() %>% #mutate(shapetype = str_replace(str_extract(type, '\\w+\\s'), ' ', '')) %>% ungroup()
        mutate(sgn = case_when(sgn == 'n' ~ '-',
                               sgn == 'p' ~ '+', 
                               TRUE ~ sgn))
      
      datndg$type <- factor(datndg$type, levels = str_replace_all(suffixes, '_', ' '))
      
      # plot
      # clf pies first
      plot <-  ggplot() + 
        geom_point(data = expand.grid(regs_abv$reg_n,regs_abv$reg_n), mapping = aes(x = Var1, y = Var2), size = 0.5, color = 'black') + 
        scatterpie::geom_scatterpie(aes(x = reg1_n_sort, y = reg2_n_sort, r = sum_scaled), data = datclf, cols = c('+', '-'), color = NA, sorted_by_radius = F)
      
      lgd <- ggplot() + 
        geom_scatterpie_legend_custom(radius = c(0, sqrt(0.5)*scaler_clf, sqrt(0.1)*scaler_clf, scaler_clf), n = 4,
                                      x = diag_end + 2.6, y = diag_end/2 + 0.25, labeller = function(r) return((r/scaler_clf)**2),
                                      size = GLOBAL_FONT_SIZE - 9)
      # ndg pies second
      plot <- plot + 
        ggforce::geom_arc_bar(aes(x0 = reg_n, y0 = reg_n, r = sum_scaled, start = start, end = end, fill = sgn, r0 = 0), data = datndg, color = NA) + 
        annotate('segment', x = diag_start, y = diag_start, xend = diag_end, yend = diag_end, color = 'grey30')
      
      lgd <- lgd + 
        geom_scatterpie_legend_custom(radius = c(0, sqrt(16)*scaler_ndg, sqrt(64)*scaler_ndg), n = 3, x = diag_end + 2.6, y = diag_end/2-1, labeller = function(r) return((r/scaler_ndg)**2), #datndg$sum_scaled[datndg$sum_scaled != 0]
                                      size = GLOBAL_FONT_SIZE - 9) + 
        annotate('text', x = diag_end + 1, y = diag_end/2 + 0.15, label = 'CLF', vjust = 0, hjust = 0, size = GLOBAL_FONT_SIZE - 8, fontface = GLOBAL_FONT_FACE_TITLE) +
        annotate('text', x = diag_end + 1, y = diag_end/2-1.15, label = 'NDG', vjust = 0, hjust = 0, size = GLOBAL_FONT_SIZE - 8, fontface = GLOBAL_FONT_FACE_TITLE) +
        annotate('rect', xmin = diag_end + 0.85, ymin = diag_end/2-1.6, xmax = diag_end + 4, ymax = diag_end / 2 + 1, fill = NA, 
                 color = {if(no_borders) {'white'} else {'black'}}, size = 1)
      
      mrg <- c(0, 0, 0, 0)
      
      plot <- plot +
        scale_alpha_manual(values = list('TRUE' = 1, 'FALSE' = 0.65), guide = FALSE) +
        scale_fill_manual(values = rev(c('#ED3537', '#1188d6')),
                          guide = {if (legend) {guide_legend('SIGN', title.position = 'left')} else {FALSE}}) + #, override.aes = list(shape = 22))) + 
        scale_x_continuous(breaks = regs_abv$reg_n, labels = regs_abv$reg) + #, sec.axis = dup_axis()) + 
        scale_y_continuous(breaks = regs_abv$reg_n, labels = regs_abv$reg) + #, sec.axis = dup_axis()) + 
        global_title_and_axis() + 
        theme(legend.direction = 'horizontal', legend.box = 'vertical', legend.position = c(1.25,0), legend.justification = c(0.5,0),
              #axis.text.x.bottom = element_text(angle = 0, hjust = 0.5, size = GLOBAL_FONT_SIZE-4),
              #axis.text.x.top = element_text(angle = 0, hjust = 0.5, size = GLOBAL_FONT_SIZE-4),
              #axis.text.y.left = element_text(angle = 0, hjust = 1, size = GLOBAL_FONT_SIZE-4),
              #axis.text.y.right = element_text(angle = 0, hjust = 0, size = GLOBAL_FONT_SIZE-4),
              text = element_text(size = GLOBAL_FONT_SIZE - 9, face = GLOBAL_FONT_FACE_TITLE),
              axis.text.x = element_text(angle = 0, hjust = 0.5, size = GLOBAL_FONT_SIZE-2),
              axis.text.y = element_text(angle = 0, hjust = 1, size = GLOBAL_FONT_SIZE-2),
              legend.title = element_text(size = GLOBAL_FONT_SIZE - 5),
              plot.margin = unit(mrg, 'cm'), axis.title = element_blank(),
              panel.grid.minor = element_blank())
      if (no_borders) {
        plot <- plot + 
          theme(panel.border = element_blank(), legend.box.background = element_blank(), legend.background = element_blank())
      }
      plot <- plot + 
        coord_equal(clip = 'off', xlim = c(diag_start,diag_end), ylim = c(diag_start,diag_end))
      return(list(plot = plot, leg_circ = lgd))
    })
  })
  
  leg_glob <- ggdraw(plts[[1]][[2]]$leg_circ + 
                       coord_equal(expand = F) +
                       theme(axis.text = element_blank(), axis.title = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
                             plot.background = element_blank(), panel.background = element_blank()), xlim = c(-0,1), c(-0.25,1)) + 
    draw_grob(get_legend(plts[[1]][[2]]$plot + 
                           theme(legend.justification = 'center',
                                 legend.key.width = unit((GLOBAL_FONT_SIZE/2 + 2) * 2.5, units = 'pt'),
                                 legend.margin = margin(t = GLOBAL_FONT_SIZE/2 + 2, r = GLOBAL_FONT_SIZE/2 + 2, b = GLOBAL_FONT_SIZE/2 + 2, l = GLOBAL_FONT_SIZE/2 + 2),
                                 legend.title = element_text(size = GLOBAL_FONT_SIZE + 10, face = GLOBAL_FONT_FACE_TITLE),
                                 text = element_text(size = GLOBAL_FONT_SIZE + 10))),
              x = -0.87, y = -0.1, width = 1.1, height = 0.24)
  
  
  if (length(suffix_pairs[[2]]) == 2) {
    plts <- list(plts[[1]][[1]]$plot + 
                   theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = 'none',
                         plot.margin = unit(c(0,-0.5,0,0), 'cm')),
                 plts[[1]][[2]]$plot + 
                   theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = 'none',
                         plot.margin = unit(c(0,-3,0,0), 'cm')),
                 plts[[2]][[1]]$plot + 
                   theme(legend.position = 'none', axis.text.x = element_blank(), axis.text.y = element_blank(),
                         plot.margin = unit(c(0,-0.5,0,0), 'cm')),
                 plts[[2]][[2]]$plot + 
                   theme(axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = 'none',
                         plot.margin = unit(c(0,-3,0,0), 'cm')))
  } else if (length(suffix_pairs[[2]]) == 1) {
    plts <- list(plts[[1]][[1]]$plot + 
                   theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = 'none',
                         plot.margin = unit(c(0,-0.5,0,0), 'cm')),
                 plts[[1]][[2]]$plot + 
                   theme(axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = 'none',
                         plot.margin = unit(c(0,-3,0,0), 'cm')),
                 plts[[2]][[1]]$plot + 
                   theme(legend.position = 'none', axis.text.x = element_blank(), axis.text.y = element_blank(),
                         plot.margin = unit(c(0,-0.5,0,0), 'cm')))
  }
  
  plts[sapply(plts, is.null)] <- NULL
  if (leg_option == 'share') {
    plt <- ggdraw(xlim = c(-0.195, 1.24), ylim = c(-0.145,1.05)) + 
      draw_plot(plot_grid(plotlist = plts, nrow = 2, ncol = 2, align = 'hv', axis = 'tblr', 
                          labels = 'AUTO', label_size = GLOBAL_FONT_SIZE + 9, label_fontfamily = GLOBAL_FONT_FAMILY, label_fontface = GLOBAL_FONT_FACE_TITLE,
                          label_x = 0.03, label_y = 0.99)) + 
      draw_plot(leg_glob, x = 0.98, y = 0.6, width = 0.29, height = 0.31, scale = 0.75) + 
      annotate('text', x = 1.015, y = 0.98, size = GLOBAL_FONT_SIZE - 6, fontface = GLOBAL_FONT_FACE_TITLE, label = 'RAW', hjust = 0) + 
      annotate('text', x = 1.015, y = 0.48, size = GLOBAL_FONT_SIZE - 6, fontface = GLOBAL_FONT_FACE_TITLE, label = 'MIL', hjust = 0) + 
      annotate('text', x = 1.12, y = 0.31, label = str_replace_all(str_to_upper(suffix_pairs[['raw']][[1]][1]), '_', '\n'), vjust = -0.07, hjust = 1.07, size = GLOBAL_FONT_SIZE - 9, fontface = GLOBAL_FONT_FACE_TITLE) + 
      annotate('text', x = 1.12, y = 0.31, label = str_replace_all(str_to_upper(suffix_pairs[['raw']][[1]][2]), '_', '\n'), vjust = 1.07, hjust = -0.07, size = GLOBAL_FONT_SIZE - 9, fontface = GLOBAL_FONT_FACE_TITLE) + 
      annotate('segment', x = 1.07, y = 0.26, xend = 1.17, yend = 0.36) + 
      annotate('text', x = 1.03, y = 0.385, size = GLOBAL_FONT_SIZE - 9, fontface = GLOBAL_FONT_FACE_TITLE, label = 'A & C', hjust = 0, vjust = 0) + 
      annotate('rect', xmin = 1.03, ymin = 0.24, xmax = 1.21, ymax = 0.38, fill = NA, color = 'black') +
      annotate('text', x = 1.12, y = 0.13, label = str_replace_all(str_to_upper(suffix_pairs[['mil']][[2]][1]), '_', '\n'), vjust = -0.07, hjust = 1.07, size = GLOBAL_FONT_SIZE - 9, fontface = GLOBAL_FONT_FACE_TITLE) + 
      annotate('text', x = 1.12, y = 0.13, label = str_replace_all(str_to_upper(suffix_pairs[['mil']][[2]][2]), '_', '\n'), vjust = 1.07, hjust = -0.07, size = GLOBAL_FONT_SIZE - 9, fontface = GLOBAL_FONT_FACE_TITLE) + 
      annotate('segment', x = 1.07, y = 0.08, xend = 1.17, yend = 0.18) + 
      annotate('text', x = 1.03, y = 0.205, size = GLOBAL_FONT_SIZE - 9, fontface = GLOBAL_FONT_FACE_TITLE, label = 'B & D', hjust = 0, vjust = 0) + 
      annotate('rect', xmin = 1.03, ymin = 0.06, xmax = 1.21, ymax = 0.2, fill = NA, color = 'black')
  } else if (leg_option == 'not_shared') {
    plt <- ggdraw(xlim = c(-0.195, 1.24), ylim = c(-0.145,1.05)) + 
      draw_plot(plot_grid(plotlist = plts, nrow = 2, ncol = 2, align = 'hv', axis = 'tblr', 
                          labels = 'AUTO', label_size = GLOBAL_FONT_SIZE + 9, label_fontfamily = GLOBAL_FONT_FAMILY, label_fontface = GLOBAL_FONT_FACE_TITLE,
                          label_x = 0.03, label_y = 0.99)) + 
      draw_plot(leg_glob, x = 0.98, y = 0.6, width = 0.29, height = 0.31, scale = 0.75) + 
      annotate('text', x = 1.015, y = 0.98, size = GLOBAL_FONT_SIZE - 6, fontface = GLOBAL_FONT_FACE_TITLE, label = 'RAW', hjust = 0) + 
      annotate('text', x = 1.015, y = 0.48, size = GLOBAL_FONT_SIZE - 6, fontface = GLOBAL_FONT_FACE_TITLE, label = 'MIL', hjust = 0) + 
      annotate('text', x = 0.66, y = 0.28, label = str_replace_all(str_to_upper(suffix_pairs[['raw']][[1]][1]), '_', '\n'), vjust = -0.07, hjust = 1.07, size = GLOBAL_FONT_SIZE - 9, fontface = GLOBAL_FONT_FACE_TITLE) + 
      annotate('text', x = 0.66, y = 0.28, label = str_replace_all(str_to_upper(suffix_pairs[['raw']][[1]][2]), '_', '\n'), vjust = 1.07, hjust = -0.07, size = GLOBAL_FONT_SIZE - 9, fontface = GLOBAL_FONT_FACE_TITLE) + 
      annotate('segment', x = 0.61, y = 0.23, xend = 0.71, yend = 0.33) + 
      annotate('text', x = 0.57, y = 0.355, size = GLOBAL_FONT_SIZE - 9, fontface = GLOBAL_FONT_FACE_TITLE, label = 'A', hjust = 0, vjust = 0) + 
      annotate('rect', xmin = 0.57, ymin = 0.21, xmax = 0.75, ymax = 0.35, fill = NA, color = 'black') +
      annotate('text', x = 0.86, y = 0.28, label = str_replace_all(str_to_upper(suffix_pairs[['raw']][[2]][1]), '_', '\n'), vjust = -0.07, hjust = 1.07, size = GLOBAL_FONT_SIZE - 9, fontface = GLOBAL_FONT_FACE_TITLE) + 
      annotate('text', x = 0.86, y = 0.28, label = str_replace_all(str_to_upper(suffix_pairs[['raw']][[2]][2]), '_', '\n'), vjust = 1.07, hjust = -0.07, size = GLOBAL_FONT_SIZE - 9, fontface = GLOBAL_FONT_FACE_TITLE) + 
      annotate('segment', x = 0.81, y = 0.23, xend = 0.91, yend = 0.33) + 
      annotate('text', x = 0.77, y = 0.355, size = GLOBAL_FONT_SIZE - 9, fontface = GLOBAL_FONT_FACE_TITLE, label = 'B', hjust = 0, vjust = 0) + 
      annotate('rect', xmin = 0.77, ymin = 0.21, xmax = 0.95, ymax = 0.35, fill = NA, color = 'black') +
      annotate('text', x = 0.66, y = 0.1, label = str_replace_all(str_to_upper(suffix_pairs[['mil']][[1]][1]), '_', '\n'), vjust = -0.07, hjust = 1.07, size = GLOBAL_FONT_SIZE - 9, fontface = GLOBAL_FONT_FACE_TITLE) + 
      annotate('text', x = 0.66, y = 0.1, label = str_replace_all(str_to_upper(suffix_pairs[['mil']][[1]][2]), '_', '\n'), vjust = 1.07, hjust = -0.07, size = GLOBAL_FONT_SIZE - 9, fontface = GLOBAL_FONT_FACE_TITLE) + 
      annotate('segment', x = 0.61, y = 0.05, xend = 0.71, yend = 0.15) + 
      annotate('text', x = 0.57, y = 0.175, size = GLOBAL_FONT_SIZE - 9, fontface = GLOBAL_FONT_FACE_TITLE, label = 'C', hjust = 0, vjust = 0) + 
      annotate('rect', xmin = 0.57, ymin = 0.03, xmax = 0.75, ymax = 0.17, fill = NA, color = 'black')
  } else if (leg_option == 'internal') {
    plt <- ggdraw(xlim = c(-0.195, 1.24), ylim = c(-0.145,1.05)) + 
      draw_plot(plot_grid(plotlist = plts, nrow = 2, ncol = 2, align = 'hv', axis = 'tblr', 
                          labels = paste0(LETTERS[1:4], ' | ', 
                                          c(paste0(str_split(suffix_pairs[['raw']][[1]][1], '_')[[1]][1], ' ', str_to_upper(str_split(suffix_pairs[['raw']][[1]][1], '_')[[1]][2])), 
                                            paste0(str_split(suffix_pairs[['raw']][[2]][1], '_')[[1]][1], ' ', str_to_upper(str_split(suffix_pairs[['raw']][[2]][1], '_')[[1]][2])),
                                            paste0(str_split(suffix_pairs[['mil']][[1]][1], '_')[[1]][1], ' ', str_to_upper(str_split(suffix_pairs[['mil']][[1]][1], '_')[[1]][2])),
                                            paste0(str_split(suffix_pairs[['mil']][[2]][1], '_')[[1]][1], ' ', str_to_upper(str_split(suffix_pairs[['mil']][[2]][1], '_')[[1]][2])))),
                          label_size = GLOBAL_FONT_SIZE + 9, label_fontfamily = GLOBAL_FONT_FAMILY, label_fontface = GLOBAL_FONT_FACE_TITLE,
                          label_x = 0.04, label_y = 0.996, hjust = 0.0)) + 
      draw_plot(leg_glob, x = 0.98, y = 0.6, width = 0.29, height = 0.31, scale = 0.75) + 
      annotate('text', x = 1.015, y = 0.98, size = GLOBAL_FONT_SIZE - 6, fontface = GLOBAL_FONT_FACE_TITLE, label = 'RAW', hjust = 0) + 
      annotate('text', x = 1.015, y = 0.48, size = GLOBAL_FONT_SIZE - 6, fontface = GLOBAL_FONT_FACE_TITLE, label = 'MIL', hjust = 0) + 
      annotate('segment', x = 0.008, y = 0.5, xend = 1.008, yend = 0.5, size = 3) + 
      annotate('text', x = 0.47, y = 0.524, label = paste0(if_else(is.na(suffix_pairs[['raw']][[1]][2]), '', str_split(suffix_pairs[['raw']][[1]][2], '_')[[1]][1]), ' ', if_else(is.na(suffix_pairs[['raw']][[1]][2]), '', str_to_upper(str_split(suffix_pairs[['raw']][[1]][2], '_')[[1]][2]))), vjust = 1.0, hjust = 1.0, size = GLOBAL_FONT_SIZE - 8, fontface = GLOBAL_FONT_FACE_TITLE) +
      annotate('text', x = 0.97, y = 0.524, label = paste0(str_split(suffix_pairs[['raw']][[2]][2], '_')[[1]][1], ' ', str_to_upper(str_split(suffix_pairs[['raw']][[2]][2], '_')[[1]][2])), vjust = 1.0, hjust = 1.0, size = GLOBAL_FONT_SIZE - 8, fontface = GLOBAL_FONT_FACE_TITLE) + 
      annotate('text', x = 0.47, y = 0.024, label = paste0(if_else(is.na(suffix_pairs[['mil']][[1]][2]), '', str_split(suffix_pairs[['mil']][[1]][2], '_')[[1]][1]), ' ', if_else(is.na(suffix_pairs[['mil']][[1]][2]), '', str_to_upper(str_split(suffix_pairs[['mil']][[1]][2], '_')[[1]][2]))), vjust = 1.0, hjust = 1.0, size = GLOBAL_FONT_SIZE - 8, fontface = GLOBAL_FONT_FACE_TITLE) + 
      annotate('text', x = 0.97, y = 0.024, label = paste0(str_split(suffix_pairs[['mil']][[2]][2], '_')[[1]][1], ' ', str_to_upper(str_split(suffix_pairs[['mil']][[2]][2], '_')[[1]][2])), vjust = 1.0, hjust = 1.0, size = GLOBAL_FONT_SIZE - 8, fontface = GLOBAL_FONT_FACE_TITLE)
  } else if (leg_option == 'not_shared_internal') {
    plt <- ggdraw(xlim = c(-0.195, 1.24), ylim = c(-0.145,1.05)) + 
      draw_plot(plot_grid(plotlist = plts, nrow = 2, ncol = 2, align = 'hv', axis = 'tblr', 
                          labels = paste0(LETTERS[1:4], ' | ', 
                                          c(paste0(str_split(suffix_pairs[['raw']][[1]][1], '_')[[1]][1], ' ', str_to_upper(str_split(suffix_pairs[['raw']][[1]][1], '_')[[1]][2])), 
                                            paste0(str_split(suffix_pairs[['raw']][[2]][1], '_')[[1]][1], ' ', str_to_upper(str_split(suffix_pairs[['raw']][[2]][1], '_')[[1]][2])),
                                            paste0(str_split(suffix_pairs[['mil']][[1]][1], '_')[[1]][1], ' ', str_to_upper(str_split(suffix_pairs[['mil']][[1]][1], '_')[[1]][2])))),
                          label_size = GLOBAL_FONT_SIZE + 9, label_fontfamily = GLOBAL_FONT_FAMILY, label_fontface = GLOBAL_FONT_FACE_TITLE,
                          label_x = 0.04, label_y = 0.996, hjust = 0.0)) + 
      draw_plot(leg_glob, x = 0.98, y = 0.6, width = 0.29, height = 0.31, scale = 0.75) + 
      annotate('text', x = 1.015, y = 0.98, size = GLOBAL_FONT_SIZE - 6, fontface = GLOBAL_FONT_FACE_TITLE, label = 'RAW', hjust = 0) + 
      annotate('text', x = 1.015, y = 0.48, size = GLOBAL_FONT_SIZE - 6, fontface = GLOBAL_FONT_FACE_TITLE, label = 'MIL', hjust = 0) + 
      annotate('segment', x = 0, y = 0.54, xend = 1, yend = 0.54, size = 3) + 
      annotate('text', x = 0.47, y = 0.524, label = paste0(str_split(suffix_pairs[['raw']][[1]][2], '_')[[1]][1], ' ', str_to_upper(str_split(suffix_pairs[['raw']][[1]][2], '_')[[1]][2])), vjust = 1.0, hjust = 1.0, size = GLOBAL_FONT_SIZE - 8, fontface = GLOBAL_FONT_FACE_TITLE) +
      annotate('text', x = 0.97, y = 0.524, label = paste0(str_split(suffix_pairs[['raw']][[2]][2], '_')[[1]][1], ' ', str_to_upper(str_split(suffix_pairs[['raw']][[2]][2], '_')[[1]][2])), vjust = 1.0, hjust = 1.0, size = GLOBAL_FONT_SIZE - 8, fontface = GLOBAL_FONT_FACE_TITLE) + 
      annotate('text', x = 0.47, y = 0.024, label = paste0(str_split(suffix_pairs[['mil']][[1]][2], '_')[[1]][1], ' ', str_to_upper(str_split(suffix_pairs[['mil']][[1]][2], '_')[[1]][2])), vjust = 1.0, hjust = 1.0, size = GLOBAL_FONT_SIZE - 8, fontface = GLOBAL_FONT_FACE_TITLE)
  }
  
  return(list(plt, leg_glob))
}

plot_nw_clf_and_ndg_matr1x2_nosig <- function(datclf, datndg, suffix_pairs,
                                              nw_meas_clf = 'pcross_link',
                                              nw_meas_ndg = 'node_degree',
                                              tscales = c('raw'),
                                              regions = REGIONS, plts_sep = FALSE,
                                              legend = TRUE, leg_option = 'share', no_borders = FALSE) {
  
  plts <- lapply(tscales, function(tscale) {
    plts <- lapply(1:length(suffix_pairs[[tscale]]), function(i){
      suffixes <- suffix_pairs[[tscale]][[i]]
      print(suffixes)
      # prepare clf data
      datclf <- datclf %>% 
        filter(cross_reg1 %in% unique(regions$region) & cross_reg2 %in% unique(regions$region)) %>% 
        mutate(cross_reg1 = substr(cross_reg1,1,2),cross_reg2 = substr(cross_reg2,1,2))
      datclf <- datclf %>% select(window, cross_reg1, cross_reg2, 
                                  !!!syms(mapply(function(nm, pr) return(paste0(pr, nm)), pr = c(paste0(nw_meas_clf, c('.p','.n'), '.')), nm = sort(rep(x = suffixes, 2)), USE.NAMES = FALSE)))
      datclf <- lapply(c(paste0(nw_meas_clf, c('.p','.n'))), function(c) {
        return(datclf %>% select(window, cross_reg1, cross_reg2, one_of(paste0(c, '.', suffixes))) %>% gather(key = 'type', value = !!sym(paste0(c, '_val')), -window, -cross_reg1, -cross_reg2) %>%
                 rowwise() %>% mutate(type = str_extract(type, suffixes) %>% .[!is.na(.)]) %>% ungroup())
      }) %>% plyr::join_all(., by = c('window', 'cross_reg1', 'cross_reg2', 'type'), type = 'inner') %>% as_tibble() %>% mutate_at(vars(window), as.character) %>% 
        rowwise %>% 
        rename(p = paste0(nw_meas_clf, '.p_val'), n = paste0(nw_meas_clf, '.n_val')) %>% 
        mutate(type = str_replace(type, '\\.', ' ')) %>% ungroup()
      datclf$type <- factor(datclf$type, levels = str_replace(suffixes, '\\.', ' '))
      
      regs_abv <- tibble(reg = distinct(bind_rows(select(datclf,cross_reg1),select(datclf,cross_reg2) %>% rename(cross_reg1 = cross_reg2)))$cross_reg1) %>% 
        rownames_to_column(var = 'reg_n') %>% 
        mutate_at(vars(reg_n), as.numeric)
      
      scaler_clf <- 0.5
      diag_start <- 0.5
      diag_end <- length(unique(c(as.vector(datclf$cross_reg1), as.vector(datclf$cross_reg2)))) + 0.5
      
      datclf <- datclf %>% 
        mutate_at(vars(cross_reg1,cross_reg2), as.character) %>% 
        inner_join(regs_abv %>% rename(cross_reg1 = reg, reg1_n = reg_n)) %>% 
        inner_join(regs_abv %>% rename(cross_reg2 = reg, reg2_n = reg_n)) %>% 
        mutate(sum = p+n, sum_scaled = sum*scaler_clf) %>% # mapping to srqt(sum) done in stat = 'pie' (consider legend below)
        filter(sum != 0) %>% 
        select(type,reg1_n,reg2_n,p,n,sum,sum_scaled)
      
      datclf <- datclf %>% 
        group_by(type) %>% 
        mutate(reg1_n_sort = case_when(type == suffixes[1] & reg1_n > reg2_n ~ reg2_n,
                                       type == suffixes[1] & reg1_n < reg2_n ~ reg1_n,
                                       type == suffixes[2] & reg1_n > reg2_n ~ reg1_n,
                                       type == suffixes[2] & reg1_n < reg2_n ~ reg2_n,
                                       TRUE ~ reg1_n),
               reg2_n_sort = case_when(type == suffixes[1] & reg1_n > reg2_n ~ reg1_n,
                                       type == suffixes[1] & reg1_n < reg2_n ~ reg2_n,
                                       type == suffixes[2] & reg1_n > reg2_n ~ reg2_n,
                                       type == suffixes[2] & reg1_n < reg2_n ~ reg1_n,
                                       TRUE ~ reg2_n)) %>% 
        rename('+' = p, '-' = n)
      
      # prepare ndg data
      datndg <- datndg %>% select(window, site_id,
                                  !!!syms(mapply(function(nm, pr) return(paste0(pr, nm)), pr = c(paste0(nw_meas_ndg, c('','.p','.n'), '.')), nm = sort(rep(x = suffixes, 3)), USE.NAMES = FALSE)))
      datndg <- lapply(c(paste0(nw_meas_ndg, c('', '.p','.n'))), function(c) {
        return(datndg %>% select(window, site_id, one_of(paste0(c, '.', suffixes))) %>% gather(key = 'type', value = !!sym(paste0(c)), -window, -site_id) %>%
                 rowwise() %>% mutate(type = str_extract(type, suffixes) %>% .[!is.na(.)]) %>% ungroup())
      }) %>% plyr::join_all(., by = c('window', 'site_id', 'type'), type = 'inner') %>% as_tibble() %>% mutate_at(vars(window), as.character) %>% 
        rowwise %>%
        gather(key = 'ndg_type', value = 'ndg_val', !!!syms(paste0(paste0(nw_meas_ndg, c('', '.p','.n'))))) %>% 
        rowwise %>%
        mutate(ndg_type = case_when(str_detect(ndg_type, '.p') ~ 'p', str_detect(ndg_type, '.n') ~ 'n', TRUE ~ 'sum')) %>% 
        #rename(p = paste0(nw_meas_ndg, '.p_val'), n = paste0(nw_meas_ndg, '.n_val')) %>% 
        filter(type %in% suffixes)
      
      scaler_ndg <- 0.05
      
      datndg <- datndg %>% 
        inner_join(regions, by = 'site_id') %>% 
        mutate(region = substr(region,1,2)) %>% 
        spread(ndg_type,ndg_val) %>% 
        group_by(region, type, window) %>% 
        summarise(sum = mean(sum, na.rm = T), p = mean(p, na.rm = T), n = mean(n, na.rm = T)) %>% 
        ungroup() %>% 
        inner_join(regs_abv %>% rename(region = reg), by = 'region') %>% 
        mutate(sum_scaled = sqrt(sum)*scaler_ndg) 
      numndg <- datndg %>% 
        #filter(ndg_type == 'node_degree') %>% 
        group_by(region) %>% 
        count() %>% mutate_at(vars(n), as.character)
      datndg <- datndg %>% 
        mutate(frac_p = case_when(p != 0 & n != 0 ~ p/sum,
                                  p == 0 & n != 0 ~ 0,
                                  p != 0 & n == 0 ~ 1,
                                  T ~ 0.5),
               frac_n = case_when(p != 0 & n != 0 ~ n/sum,
                                  p == 0 & n != 0 ~ 1,
                                  p != 0 & n == 0 ~ 0,
                                  T ~ 0.5)) %>% 
        gather(key = 'sgn', value = 'single',p,n) %>% 
        mutate(start = case_when(type == suffixes[1] & sgn == 'p' ~ 1/4*pi - frac_p*pi,
                                 type == suffixes[2] & sgn == 'p' ~ 1/4*pi,
                                 type == suffixes[1] & sgn == 'n' ~ -3/4*pi,
                                 type == suffixes[2] & sgn == 'n' ~ 5/4*pi - frac_n*pi),
               end = case_when(type == suffixes[1] & sgn == 'p' ~ 1/4*pi,
                               type == suffixes[2] & sgn == 'p' ~ 1/4*pi + frac_p*pi,
                               type == suffixes[1] & sgn == 'n' ~ -3/4*pi + frac_n*pi,
                               type == suffixes[2] & sgn == 'n' ~ 5/4*pi)) %>% 
        mutate(type = str_replace_all(type, '_', ' ')) %>% ungroup() %>% #mutate(shapetype = str_replace(str_extract(type, '\\w+\\s'), ' ', '')) %>% ungroup()
        mutate(sgn = case_when(sgn == 'n' ~ '-',
                               sgn == 'p' ~ '+', 
                               TRUE ~ sgn))
      
      datndg$type <- factor(datndg$type, levels = str_replace_all(suffixes, '_', ' '))
      
      # plot
      # clf pies first
      plot <-  ggplot() + 
        geom_point(data = expand.grid(regs_abv$reg_n,regs_abv$reg_n), mapping = aes(x = Var1, y = Var2), size = 0.5, color = 'black') + 
        scatterpie::geom_scatterpie(aes(x = reg1_n_sort, y = reg2_n_sort, r = sum_scaled), data = datclf, cols = c('+', '-'), color = NA, sorted_by_radius = F)
      
      lgd <- ggplot() + 
        geom_scatterpie_legend_custom(radius = c(0, sqrt(0.5)*scaler_clf, sqrt(0.1)*scaler_clf, scaler_clf), n = 4,
                                      x = diag_end + 2.6, y = diag_end/2 + 0.25, labeller = function(r) return((r/scaler_clf)**2),
                                      size = GLOBAL_FONT_SIZE - 9)
      
      # ndg pies second
      plot <- plot + 
        ggforce::geom_arc_bar(aes(x0 = reg_n, y0 = reg_n, r = sum_scaled, start = start, end = end, fill = sgn, r0 = 0), data = datndg, color = NA) + 
        annotate('segment', x = diag_start, y = diag_start, xend = diag_end, yend = diag_end, color = 'grey30')
      
      lgd <- lgd + 
        geom_scatterpie_legend_custom(radius = c(0, sqrt(16)*scaler_ndg, sqrt(64)*scaler_ndg), n = 3, x = diag_end + 2.6, y = diag_end/2-1, labeller = function(r) return((r/scaler_ndg)**2), #datndg$sum_scaled[datndg$sum_scaled != 0]
                                      size = GLOBAL_FONT_SIZE - 9) + 
        annotate('text', x = diag_end + 1, y = diag_end/2 + 0.15, label = 'CLF', vjust = 0, hjust = 0, size = GLOBAL_FONT_SIZE - 8, fontface = GLOBAL_FONT_FACE_TITLE) +
        annotate('text', x = diag_end + 1, y = diag_end/2-1.15, label = 'NDG', vjust = 0, hjust = 0, size = GLOBAL_FONT_SIZE - 8, fontface = GLOBAL_FONT_FACE_TITLE)
      
      if (!plts_sep) {
        lgd <- lgd + 
          annotate('rect', xmin = diag_end + 0.85, ymin = diag_end/2-1.6, xmax = diag_end + 4, ymax = diag_end / 2 + 1, fill = NA, 
                   color = {if(no_borders) {'white'} else {'black'}}, size = 1)
      } else {
        lgd <- lgd + 
          annotate('rect', xmin = diag_end + 0.85, ymin = diag_end/2-1.6, xmax = diag_end + 4, ymax = diag_end / 2 + 1, fill = NA, 
                   color = 'white', size = 1)
      }
      
      mrg <- c(0, 0, 0, 0)
      if (plts_sep) sign_tlt <- 'correlation\nsign' else sign_tlt <- 'SIGN'
      
      #print(regs_abv)
      plot <- plot +
        scale_alpha_manual(values = list('TRUE' = 1, 'FALSE' = 0.65), guide = FALSE) +
        scale_fill_manual(values = rev(c('#ED3537', '#1188d6')),
                          guide = {if (legend & plts_sep) {guide_legend(sign_tlt, title.position = 'top')} else if (legend) {guide_legend(sign_tlt, title.position = 'left')} else {FALSE}}) + #, override.aes = list(shape = 22))) + 
        scale_x_continuous(breaks = regs_abv$reg_n, labels = regs_abv$reg) + #, sec.axis = dup_axis()) + 
        scale_y_continuous(breaks = regs_abv$reg_n, labels = regs_abv$reg) + #, sec.axis = dup_axis()) + 
        global_title_and_axis() + 
        theme(legend.direction = 'horizontal', legend.box = 'vertical', legend.position = c(1.25,0), legend.justification = c(0.5,0),
              text = element_text(size = GLOBAL_FONT_SIZE - 9, face = GLOBAL_FONT_FACE_TITLE),
              axis.text.x = element_text(angle = 0, hjust = 0.5, size = GLOBAL_FONT_SIZE-2),
              axis.text.y = element_text(angle = 0, hjust = 1, size = GLOBAL_FONT_SIZE-2),
              legend.title = element_text(size = GLOBAL_FONT_SIZE - 5),
              plot.margin = unit(mrg, 'cm'), axis.title = element_blank(),
              panel.grid.minor = element_blank())
      if (no_borders) {
        plot <- plot + 
          theme(panel.border = element_blank(), legend.box.background = element_blank(), legend.background = element_blank())
      }
      plot <- plot + 
        coord_equal(clip = 'off', xlim = c(diag_start,diag_end), ylim = c(diag_start,diag_end))
      return(list(plot = plot, leg_circ = lgd))
    })
  })
  
  if (plts_sep) {
    leg1 <- plts[[1]][[1]]$leg_circ + 
      coord_equal(expand = F) +
      theme(axis.text = element_blank(), axis.title = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
            plot.background = element_blank(), panel.background = element_blank(),
            legend.box.background = element_blank(), legend.background = element_blank())
    
    leg2 <- get_legend(plts[[1]][[2]]$plot + 
                             theme(legend.justification = 'center',
                                   legend.key.width = unit((GLOBAL_FONT_SIZE/2 + 2) * 2.5, units = 'pt'),
                                   legend.margin = margin(t = GLOBAL_FONT_SIZE/2 + 2, r = GLOBAL_FONT_SIZE/2 + 2, b = GLOBAL_FONT_SIZE/2 + 2, l = GLOBAL_FONT_SIZE/2 + 2),
                                   legend.title = element_text(size = GLOBAL_FONT_SIZE + 10, face = GLOBAL_FONT_FACE_TITLE),
                                   text = element_text(size = GLOBAL_FONT_SIZE + 10),
                                   legend.box.background = element_blank(), legend.background = element_blank()))
    
    plts <- list(plts[[1]][[1]]$plot + 
                   theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = 'none',
                         plot.margin = unit(c(0,0,0,0), 'cm')),
                 plts[[1]][[2]]$plot + 
                   theme(axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = 'none',
                         plot.margin = unit(c(0,0,0,0), 'cm')))
    
    plts[sapply(plts, is.null)] <- NULL
    
    return(list(plts = plts, leg1 = leg1, leg2 = leg2))
  } else {
    # only leg_option == 'not_shared_internal'
    leg_glob <- ggdraw(plts[[1]][[2]]$leg_circ + 
                         coord_equal(expand = F) +
                         theme(axis.text = element_blank(), axis.title = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
                               plot.background = element_blank(), panel.background = element_blank()), xlim = c(-0,1), c(-0.25,1)) + 
      draw_grob(get_legend(plts[[1]][[2]]$plot + 
                             theme(legend.justification = 'center',
                                   legend.key.width = unit((GLOBAL_FONT_SIZE/2 + 2) * 2.5, units = 'pt'),
                                   legend.margin = margin(t = GLOBAL_FONT_SIZE/2 + 2, r = GLOBAL_FONT_SIZE/2 + 2, b = GLOBAL_FONT_SIZE/2 + 2, l = GLOBAL_FONT_SIZE/2 + 2),
                                   legend.title = element_text(size = GLOBAL_FONT_SIZE + 10, face = GLOBAL_FONT_FACE_TITLE),
                                   text = element_text(size = GLOBAL_FONT_SIZE + 10))),
                x = -0.87, y = -0.1, width = 1.1, height = 0.24)
    
    plts <- list(plts[[1]][[1]]$plot + 
                   theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = 'none',
                         plot.margin = unit(c(0,-0.5,0,0), 'cm')),
                 plts[[1]][[2]]$plot + 
                   theme(axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = 'none',
                         plot.margin = unit(c(0,-3,0,0), 'cm')))
    
    plts[sapply(plts, is.null)] <- NULL
    
    
    plt <- ggdraw(xlim = c(-0.195, 1.24), ylim = c(-0.145,1.05)) + 
      draw_plot(plot_grid(plotlist = plts, nrow = 1, ncol = 2, align = 'hv', axis = 'tblr', 
                          labels = paste0(LETTERS[1:2], ' | ', 
                                          c(paste0(if_else(is.na(suffix_pairs[[tscales[1]]][[1]][1]), '', str_split(suffix_pairs[[tscales[1]]][[1]][1], '_')[[1]][1]), ' ', if_else(is.na(suffix_pairs[[tscales[1]]][[1]][1]), '',str_to_upper(str_split(suffix_pairs[[tscales[1]]][[1]][1], '_')[[1]][2]))), 
                                            paste0(str_split(suffix_pairs[[tscales[1]]][[2]][1], '_')[[1]][1], ' ', str_to_upper(str_split(suffix_pairs[[tscales[1]]][[2]][1], '_')[[1]][2])))),
                          label_size = GLOBAL_FONT_SIZE + 9, label_fontfamily = GLOBAL_FONT_FAMILY, label_fontface = GLOBAL_FONT_FACE_TITLE,
                          label_x = 0.04, label_y = 0.761, hjust = 0.0)) + 
      draw_plot(leg_glob, x = 0.98, y = 0.35, width = 0.29, height = 0.31, scale = 0.75) + 
      annotate('text', x = 0.47, y = 0.267, label = paste0(if_else(is.na(suffix_pairs[[tscales[1]]][[1]][2]), '', str_split(suffix_pairs[[tscales[1]]][[1]][2], '_')[[1]][1]), ' ', if_else(is.na(suffix_pairs[[tscales[1]]][[1]][2]), '', str_to_upper(str_split(suffix_pairs[[tscales[1]]][[1]][2], '_')[[1]][2]))), vjust = 1.0, hjust = 1.0, size = GLOBAL_FONT_SIZE - 8, fontface = GLOBAL_FONT_FACE_TITLE) +
      annotate('text', x = 0.97, y = 0.267, label = paste0(str_split(suffix_pairs[[tscales[1]]][[2]][2], '_')[[1]][1], ' ', str_to_upper(str_split(suffix_pairs[[tscales[1]]][[2]][2], '_')[[1]][2])), vjust = 1.0, hjust = 1.0, size = GLOBAL_FONT_SIZE - 8, fontface = GLOBAL_FONT_FACE_TITLE)
    
    return(list(plt, leg_glob))
  }
}

plot_nw_clf_and_ndg_matr1x1_nosig <- function(datclf, datndg, suffix_pairs,
                                              nw_meas_clf = 'pcross_link',
                                              nw_meas_ndg = 'node_degree',
                                              tscales = 'raw',
                                              regions = REGIONS, plts_sep = FALSE,
                                              legend = TRUE, leg_option = 'share', no_borders = FALSE) {
  
  plts <- lapply(tscales, function(tscale) {
    plts <- lapply(1:length(suffix_pairs[[tscale]]), function(i){
      suffixes <- suffix_pairs[[tscale]][[i]]
      print(suffixes)
      # prepare clf data
      datclf <- datclf %>% 
        filter(cross_reg1 %in% unique(regions$region) & cross_reg2 %in% unique(regions$region)) %>% 
        mutate(cross_reg1 = substr(cross_reg1,1,2),cross_reg2 = substr(cross_reg2,1,2))
      datclf <- datclf %>% select(window, cross_reg1, cross_reg2, 
                                  !!!syms(mapply(function(nm, pr) return(paste0(pr, nm)), pr = c(paste0(nw_meas_clf, c('.p','.n'), '.')), nm = sort(rep(x = suffixes, 2)), USE.NAMES = FALSE)))
      datclf <- lapply(c(paste0(nw_meas_clf, c('.p','.n'))), function(c) {
        return(datclf %>% select(window, cross_reg1, cross_reg2, one_of(paste0(c, '.', suffixes))) %>% gather(key = 'type', value = !!sym(paste0(c, '_val')), -window, -cross_reg1, -cross_reg2) %>%
                 rowwise() %>% mutate(type = str_extract(type, suffixes) %>% .[!is.na(.)]) %>% ungroup())
      }) %>% plyr::join_all(., by = c('window', 'cross_reg1', 'cross_reg2', 'type'), type = 'inner') %>% as_tibble() %>% mutate_at(vars(window), as.character) %>% 
        rowwise %>% 
        rename(p = paste0(nw_meas_clf, '.p_val'), n = paste0(nw_meas_clf, '.n_val')) %>% 
        mutate(type = str_replace(type, '\\.', ' ')) %>% ungroup()
      datclf$type <- factor(datclf$type, levels = str_replace(suffixes, '\\.', ' '))
      
      regs_abv <- tibble(reg = distinct(bind_rows(select(datclf,cross_reg1),select(datclf,cross_reg2) %>% rename(cross_reg1 = cross_reg2)))$cross_reg1) %>% 
        rownames_to_column(var = 'reg_n') %>% 
        mutate_at(vars(reg_n), as.numeric)
      
      scaler_clf <- 0.5
      diag_start <- 0.5
      diag_end <- length(unique(c(as.vector(datclf$cross_reg1), as.vector(datclf$cross_reg2)))) + 0.5
      
      datclf <- datclf %>% 
        mutate_at(vars(cross_reg1,cross_reg2), as.character) %>% 
        inner_join(regs_abv %>% rename(cross_reg1 = reg, reg1_n = reg_n)) %>% 
        inner_join(regs_abv %>% rename(cross_reg2 = reg, reg2_n = reg_n)) %>% 
        mutate(sum = p+n, sum_scaled = sum*scaler_clf) %>% # mapping to srqt(sum) done in stat = 'pie' (consider legend below)
        filter(sum != 0) %>% 
        select(type,reg1_n,reg2_n,p,n,sum,sum_scaled)
      
      datclf <- datclf %>% 
        group_by(type) %>% 
        mutate(reg1_n_sort = case_when(type == suffixes[1] & reg1_n > reg2_n ~ reg2_n,
                                       type == suffixes[1] & reg1_n < reg2_n ~ reg1_n,
                                       type == suffixes[2] & reg1_n > reg2_n ~ reg1_n,
                                       type == suffixes[2] & reg1_n < reg2_n ~ reg2_n,
                                       TRUE ~ reg1_n),
               reg2_n_sort = case_when(type == suffixes[1] & reg1_n > reg2_n ~ reg1_n,
                                       type == suffixes[1] & reg1_n < reg2_n ~ reg2_n,
                                       type == suffixes[2] & reg1_n > reg2_n ~ reg2_n,
                                       type == suffixes[2] & reg1_n < reg2_n ~ reg1_n,
                                       TRUE ~ reg2_n)) %>% 
        rename('+' = p, '-' = n)
      
      # prepare ndg data
      datndg <- datndg %>% select(window, site_id,
                                  !!!syms(mapply(function(nm, pr) return(paste0(pr, nm)), pr = c(paste0(nw_meas_ndg, c('','.p','.n'), '.')), nm = sort(rep(x = suffixes, 3)), USE.NAMES = FALSE)))
      datndg <- lapply(c(paste0(nw_meas_ndg, c('', '.p','.n'))), function(c) {
        return(datndg %>% select(window, site_id, one_of(paste0(c, '.', suffixes))) %>% gather(key = 'type', value = !!sym(paste0(c)), -window, -site_id) %>%
                 rowwise() %>% mutate(type = str_extract(type, suffixes) %>% .[!is.na(.)]) %>% ungroup())
      }) %>% plyr::join_all(., by = c('window', 'site_id', 'type'), type = 'inner') %>% as_tibble() %>% mutate_at(vars(window), as.character) %>% 
        rowwise %>%
        gather(key = 'ndg_type', value = 'ndg_val', !!!syms(paste0(paste0(nw_meas_ndg, c('', '.p','.n'))))) %>% 
        rowwise %>%
        mutate(ndg_type = case_when(str_detect(ndg_type, '.p') ~ 'p', str_detect(ndg_type, '.n') ~ 'n', TRUE ~ 'sum')) %>% 
        filter(type %in% suffixes)
      
      scaler_ndg <- 0.05
      
      datndg <- datndg %>% 
        inner_join(regions, by = 'site_id') %>% 
        mutate(region = substr(region,1,2)) %>% 
        spread(ndg_type,ndg_val) %>% 
        group_by(region, type, window) %>% 
        summarise(sum = mean(sum, na.rm = T), p = mean(p, na.rm = T), n = mean(n, na.rm = T)) %>% 
        ungroup() %>% 
        inner_join(regs_abv %>% rename(region = reg), by = 'region') %>% 
        mutate(sum_scaled = sqrt(sum)*scaler_ndg) 
      numndg <- datndg %>% 
        group_by(region) %>% 
        count() %>% mutate_at(vars(n), as.character)
      datndg <- datndg %>% 
        mutate(frac_p = case_when(p != 0 & n != 0 ~ p/sum,
                                  p == 0 & n != 0 ~ 0,
                                  p != 0 & n == 0 ~ 1,
                                  T ~ 0.5),
               frac_n = case_when(p != 0 & n != 0 ~ n/sum,
                                  p == 0 & n != 0 ~ 1,
                                  p != 0 & n == 0 ~ 0,
                                  T ~ 0.5)) %>% 
        gather(key = 'sgn', value = 'single',p,n) %>% 
        mutate(start = case_when(type == suffixes[1] & sgn == 'p' ~ 1/4*pi - frac_p*pi,
                                 type == suffixes[2] & sgn == 'p' ~ 1/4*pi,
                                 type == suffixes[1] & sgn == 'n' ~ -3/4*pi,
                                 type == suffixes[2] & sgn == 'n' ~ 5/4*pi - frac_n*pi),
               end = case_when(type == suffixes[1] & sgn == 'p' ~ 1/4*pi,
                               type == suffixes[2] & sgn == 'p' ~ 1/4*pi + frac_p*pi,
                               type == suffixes[1] & sgn == 'n' ~ -3/4*pi + frac_n*pi,
                               type == suffixes[2] & sgn == 'n' ~ 5/4*pi)) %>% 
        mutate(type = str_replace_all(type, '_', ' ')) %>% ungroup() %>% #mutate(shapetype = str_replace(str_extract(type, '\\w+\\s'), ' ', '')) %>% ungroup()
        mutate(sgn = case_when(sgn == 'n' ~ '-',
                               sgn == 'p' ~ '+', 
                               TRUE ~ sgn))
      
      datndg$type <- factor(datndg$type, levels = str_replace_all(suffixes, '_', ' '))
      
      # plot
      # clf pies first
      plot <-  ggplot() + 
        geom_point(data = expand.grid(regs_abv$reg_n,regs_abv$reg_n), mapping = aes(x = Var1, y = Var2), size = 0.5, color = 'black') + 
        scatterpie::geom_scatterpie(aes(x = reg1_n_sort, y = reg2_n_sort, r = sum_scaled), data = datclf, cols = c('+', '-'), color = NA, sorted_by_radius = F)
      
      lgd <- ggplot() + 
        geom_scatterpie_legend_custom(radius = c(0, sqrt(0.5)*scaler_clf, sqrt(0.1)*scaler_clf, scaler_clf), n = 4,
                                      x = diag_end + 2.6, y = diag_end/2 + 0.25, labeller = function(r) return((r/scaler_clf)**2),
                                      size = GLOBAL_FONT_SIZE - 9)
      
      # ndg pies second
      plot <- plot + 
        ggforce::geom_arc_bar(aes(x0 = reg_n, y0 = reg_n, r = sum_scaled, start = start, end = end, fill = sgn, r0 = 0), data = datndg, color = NA) + 
        annotate('segment', x = diag_start, y = diag_start, xend = diag_end, yend = diag_end, color = 'grey30')
      
      lgd <- lgd + 
        geom_scatterpie_legend_custom(radius = c(0, sqrt(16)*scaler_ndg, sqrt(64)*scaler_ndg), n = 3, x = diag_end + 2.6, y = diag_end/2-1, labeller = function(r) return((r/scaler_ndg)**2), #datndg$sum_scaled[datndg$sum_scaled != 0]
                                      size = GLOBAL_FONT_SIZE - 9) + 
        annotate('text', x = diag_end + 1, y = diag_end/2 + 0.15, label = 'CLF', vjust = 0, hjust = 0, size = GLOBAL_FONT_SIZE - 8, fontface = GLOBAL_FONT_FACE_TITLE) +
        annotate('text', x = diag_end + 1, y = diag_end/2-1.15, label = 'NDG', vjust = 0, hjust = 0, size = GLOBAL_FONT_SIZE - 8, fontface = GLOBAL_FONT_FACE_TITLE)
      
      if (!plts_sep) {
        lgd <- lgd + 
          annotate('rect', xmin = diag_end + 0.85, ymin = diag_end/2-1.6, xmax = diag_end + 4, ymax = diag_end / 2 + 1, fill = NA, 
                   color = {if(no_borders) {'white'} else {'black'}}, size = 1)
      } else {
        lgd <- lgd + 
          annotate('rect', xmin = diag_end + 0.85, ymin = diag_end/2-1.6, xmax = diag_end + 4, ymax = diag_end / 2 + 1, fill = NA, 
                   color = 'white', size = 1)
      }
      
      mrg <- c(0, 0, 0, 0)
      if (plts_sep) sign_tlt <- 'correlation\nsign' else sign_tlt <- 'SIGN'
      
      #print(regs_abv)
      plot <- plot +
        scale_alpha_manual(values = list('TRUE' = 1, 'FALSE' = 0.65), guide = FALSE) +
        scale_fill_manual(values = rev(c('#ED3537', '#1188d6')),
                          guide = {if (legend & plts_sep) {guide_legend(sign_tlt, title.position = 'top')} else if (legend) {guide_legend(sign_tlt, title.position = 'left')} else {FALSE}}) + #, override.aes = list(shape = 22))) + 
        scale_x_continuous(breaks = regs_abv$reg_n, labels = regs_abv$reg) + 
        scale_y_continuous(breaks = regs_abv$reg_n, labels = regs_abv$reg) + 
        global_title_and_axis() + 
        theme(legend.direction = 'horizontal', legend.box = 'vertical', legend.position = c(1.25,0), legend.justification = c(0.5,0),
              text = element_text(size = GLOBAL_FONT_SIZE - 9, face = GLOBAL_FONT_FACE_TITLE),
              axis.text.x = element_text(angle = 0, hjust = 0.5, size = GLOBAL_FONT_SIZE-2),
              axis.text.y = element_text(angle = 0, hjust = 1, size = GLOBAL_FONT_SIZE-2),
              legend.title = element_text(size = GLOBAL_FONT_SIZE - 5),
              plot.margin = unit(mrg, 'cm'), axis.title = element_blank(),
              panel.grid.minor = element_blank())
      if (no_borders) {
        plot <- plot + 
          theme(panel.border = element_blank(), legend.box.background = element_blank(), legend.background = element_blank())
      }
      plot <- plot + 
        coord_equal(clip = 'off', xlim = c(diag_start,diag_end), ylim = c(diag_start,diag_end))
      return(list(plot = plot, leg_circ = lgd))
    })
  })
  
  if (plts_sep) {
    leg1 <- plts[[1]][[1]]$leg_circ + 
      coord_equal(expand = F) +
      theme(axis.text = element_blank(), axis.title = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
            plot.background = element_blank(), panel.background = element_blank(),
            legend.box.background = element_blank(), legend.background = element_blank())
    
    leg2 <- get_legend(plts[[1]][[1]]$plot + 
                         theme(legend.justification = 'center',
                               legend.key.width = unit((GLOBAL_FONT_SIZE/2 + 2) * 2.5, units = 'pt'),
                               legend.margin = margin(t = GLOBAL_FONT_SIZE/2 + 2, r = GLOBAL_FONT_SIZE/2 + 2, b = GLOBAL_FONT_SIZE/2 + 2, l = GLOBAL_FONT_SIZE/2 + 2),
                               legend.title = element_text(size = GLOBAL_FONT_SIZE + 10, face = GLOBAL_FONT_FACE_TITLE),
                               text = element_text(size = GLOBAL_FONT_SIZE + 10),
                               legend.box.background = element_blank(), legend.background = element_blank()))
    
    plts <- list(plts[[1]][[1]]$plot + 
                   theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = 'none',
                         plot.margin = unit(c(0,0,0,0), 'cm')))
    
    plts[sapply(plts, is.null)] <- NULL
    
    return(list(plts = plts, leg1 = leg1, leg2 = leg2))
  } else {
    # only leg_option == 'not_shared_internal'
    leg_glob <- ggdraw(plts[[1]][[1]]$leg_circ + 
                         coord_equal(expand = F) +
                         theme(axis.text = element_blank(), axis.title = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
                               plot.background = element_blank(), panel.background = element_blank()), xlim = c(-0,1), c(-0.25,1)) + 
      draw_grob(get_legend(plts[[1]][[1]]$plot + 
                             theme(legend.justification = 'center',
                                   legend.key.width = unit((GLOBAL_FONT_SIZE/2 + 2) * 2.5, units = 'pt'),
                                   legend.margin = margin(t = GLOBAL_FONT_SIZE/2 + 2, r = GLOBAL_FONT_SIZE/2 + 2, b = GLOBAL_FONT_SIZE/2 + 2, l = GLOBAL_FONT_SIZE/2 + 2),
                                   legend.title = element_text(size = GLOBAL_FONT_SIZE + 10, face = GLOBAL_FONT_FACE_TITLE),
                                   text = element_text(size = GLOBAL_FONT_SIZE + 10))),
                x = -0.87, y = -0.1, width = 1.1, height = 0.24)
    
    plts <- plts[[1]][[1]]$plot + 
                   theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = 'none',
                         plot.margin = unit(c(0,-0.5,0,0), 'cm'))
    
    plt <- ggdraw(xlim = c(-0.195, 1.24), ylim = c(-0.145,1.05)) + 
      draw_plot(plot_grid(NULL, plts, nrow = 1, ncol = 2, align = 'hv', axis = 'tblr')) + #, use ncol = 2 here for correct scaling 
      draw_plot(leg_glob, x = 0.98, y = 0.35, width = 0.29, height = 0.31, scale = 0.75) 
    
    return(list(plt, leg_glob))
  }
}


## plot function for the PCA loading map with PC time series inset ----
### combined plot function
plot_pca_map <- function(pcares,pcval_lims=NULL,ldg_lims=NULL,pc=1,clb_sep=TRUE,plts_sep=FALSE,clb_ttl_pos = 'top') {
  # data for inset of PC1
  pcts <- tibble(time = pcares$time_axis,
                 pcval = pcares$pca$x[,pc],
                 pccomp = rep('PC1',dim(pcares$pca$x)[1]))
  
  # data for loadings of PC1
  pcld <- tibble(lon = pcares$lon,
                 lat = pcares$lat,
                 loading = pcares$pca$rotation[,pc[1]],
                 pccomp = rep('PC1',dim(pcares$pca$rotation)[1]))
  locs <- as_tibble(project(cbind(pcld$lon, pcld$lat), 
                            proj = GLOBAL_CRS[['robinson']])) %>% 
    rename(lon = V1, lat = V2)
  pcld$lon <- locs$lon
  pcld$lat <- locs$lat
  
  ## turn PC1 for increase over deglaciation (sign of PCs are arbitrary)
  if (pcts$pcval[which.min(pcts$time)] < pcts$pcval[which.max(pcts$time)]) {
    pcts <- pcts %>% 
      mutate(pcval = -1*pcval)
    pcld <- pcld %>% 
      mutate(loading = -1*loading)
  }
  
  # plot for time series inset
  pltts <- ggplot(data = pcts,
                  mapping = aes(x = time/1e3, y = pcval, color = pccomp)) + 
    geom_line(alpha = 0.75, size = 1.5) + 
    scale_color_manual(values = rev(chronosphere::ipccPrec(11)),#viridis(10)[c(3,5,7,9)], 
                       guide = guide_legend(title = 'Component', title.position = 'top', override.aes = list(size = 4), nrow = 2)) +
    scale_y_continuous(limits = pcval_lims) + 
    labs(x = 'time [ka BP]', y = {if(length(pc) > 1) {'PC [a.u.]'} else {paste0('PC',pc,' [a.u.]')}}) + 
    global_title_and_axis() + 
    coord_cartesian(clip = 'off', xlim = c(min(pcts$time/1e3),max(pcts$time/1e3)), ylim = c(min(pcts$pcval),max(pcts$pcval))) + 
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          legend.position = {if(length(pc) > 1) {'bottom'} else {'none'}}, 
          legend.direction = 'horizontal',
          panel.grid.minor = element_blank(),
          #panel.grid.major.y = element_blank(),
          #panel.grid.major.x = element_line(color = 'grey50'),
          panel.spacing.y = unit(0, 'lines'), 
          panel.border = my_border(), #element_blank(), 
          plot.background = element_rect(size = 1, fill = 'grey98', color = 'black'), 
          #plot.background = element_blank(),
          panel.background = element_rect(fill = 'grey98'),#, color = 'black'),
          strip.background = element_blank(), strip.text.y = element_blank(),
          plot.margin = unit(c(1,1,1,1), 'mm'),
          title = element_text(size = GLOBAL_FONT_SIZE + 6), text = element_text(size = GLOBAL_FONT_SIZE + 6)) + 
    annotate('label', label = paste0('expl. var.: ',paste(round(pcares$expl_var[pc],2)*100,collapse = ', '),'%'), 
             x = max(pcts$time/1e3), y = max(pcts$pcval), vjust = 0.65, hjust = 0.91, label.size = 0,
             fill = 'grey98', label.padding = unit(0.15, "lines"),
             fontface = GLOBAL_FONT_FACE_TITLE, size = GLOBAL_FONT_SIZE - 11)
  
  # map plot of loadings
  pltld <- shapefile_map_plot_(ggobj = ggplot(), projection = 'robinson', grat_labs = list(x = c(-180, -90, 0, 90, 180), y = c(-90, -60, -30, 0, 30, 60, 90)),
                               graticules = 30, no_bathy_clr = c('grey55', 'grey98'), #bg_clr = rev(c('#96c1e3', '#e8e8d1', '#fcfcfc')),#'#cae1f4', '#e8e8e8',# c('grey75', 'grey55', 'grey98'),
                               lgm_ice = FALSE, cb_friendly = TRUE, bg = NULL, zoom = c(-180, -65, 180, 65),
                               yaxislclip = -90, xaxisrclip = 180, xaxislclip = -180, latlabs = 'outside') + 
    geom_point(data = pcld, mapping = aes(x = lon, y = lat, fill = loading), size = 5, alpha=I(0.6), stroke = 1, shape = 21) + 
    scale_fill_gradientn(colours = chronosphere::ipccPrec(11), 
                         limits = ldg_lims,
                         guide = guide_colorbar(paste0('PC',pc[1],'\nloading [a.u.]'), direction = 'horizontal', title.position = clb_ttl_pos,
                                                barheight = 1.5, barwidth = 16))
  if (clb_sep) leg1 <- get_legend((pltld + theme(legend.box.background = element_blank(), legend.background = element_blank(),
                                                 title = element_text(size = GLOBAL_FONT_SIZE + 9), text = element_text(size = GLOBAL_FONT_SIZE + 9))))
  
  pltld <- pltld + 
    theme(legend.position = {if(clb_sep) {'none'} else {'right'}}, panel.border = element_blank())
  
  # extract legends
  #leg1 <- get_legend(pltld)
  #if (length(pc)>1) leg2 <- get_legend(pltts)
  
  # combine map with inset
  if (plts_sep) {
    plt <- list(ts = pltts, map = pltld)
  } else {
    plt <- ggdraw(pltld) + 
      draw_plot(pltts, x = 0.07, y = 0.24, width = 0.18, height = 0.3)
  }
  
  if (!clb_sep) return(plt) else return(list(plt = plt, clb = leg1))
}


## Figure 5: TRACE ----
### CLF and NDG bubble matrices and continent strips
trc_suff <- list(
  'raw' = list(c('ACER_ap','TraCE_tsf'), c('TraCE_pr','TraCE_ts')),
  'mil' = list(c('ACER_ap','TraCE_tsf'), c('TraCE_pr','TraCE_ts')),
  'unused' = c('TraCE_vc')
)

matrs_cmb <- plot_nw_clf_and_ndg_matr1x2_nosig(ACERtrc_clp_raw %>% rename_if(., str_detect(names(.), 'TRACE_apsb'), ~ str_replace(., 'TRACE_apsb', 'TRACE_tsf')) %>% 
                                                     rename_if(., str_detect(names(.), 'TRACE'), ~ str_replace(., 'TRACE', 'TraCE')),
                                                   ACERtrc_ndg_raw %>% rename_if(., str_detect(names(.), 'TRACE_apsb'), ~ str_replace(., 'TRACE_apsb', 'TRACE_tsf')) %>% 
                                                     rename_if(., str_detect(names(.), 'TRACE'), ~ str_replace(., 'TRACE', 'TraCE')),
                                                   plts_sep = TRUE,
                                                   suffix_pairs = trc_suff[1], legend = T, leg_option = 'not_shared_internal')

cntpltx <- plot_continents_row('horizontal', 'aligned',list_only = TRUE)
cntplty <- plot_continents_row('vertical', 'aligned',list_only = TRUE)

# for saving the individual bubble plot
#ggsave(plot = matrs_cmb, filename = file.path(DIR_FIGURES, 'CLF_NDG_ACERTRACE_raw_updnorm_correctedNA_woanom_compare.pdf'), 
#       width = 54, height = 45, units = 'cm', device = cairo_pdf) # height = 45

### maps and insets of PCA results
pcval_lims <- lapply(pca_results_filt, function(p) max(abs(p$pca$x[,1]))) %>% 
  unlist(use.names = FALSE) %>% max()
pcval_lims <- c(-1*pcval_lims,pcval_lims) %>% round(.,2) %>% add(.,c(-0.01,0.01))
ldg_lims <- lapply(pca_results_filt, function(p) max(abs(p$pca$rotation[,1]))) %>% 
  unlist(use.names = FALSE) %>% max()
ldg_lims <- c(-1*ldg_lims,ldg_lims) %>% round(.,2) %>% add(.,c(-0.01,0.01))
clb_poss <- c('top','left','left','top')
pca_maps <- lapply(1:length(pca_results_filt), function(i) plot_pca_map(pcares = pca_results_filt[[i]], 
                                                                   pcval_lims = pcval_lims, ldg_lims = ldg_lims, 
                                                                   plts_sep = TRUE, clb_ttl_pos = clb_poss[[i]])) %>%
  setNames(names(pca_results_filt)) # TS, PR already for S14 below


### combine step by step
cntlx_o <- list(Africa = 0.09, Asia = 0.1385,'Australia Oceania' = 0.188,Europe = 0.238,'North America' = 0.2855,'South America' = 0.335)
cntly_o <- list(Africa = 0.485, Asia = 0.575,'Australia Oceania' = 0.655,Europe = 0.74,'North America' = 0.82,'South America' = 0.9)

plt_matr_row <- plot_grid(NULL,matrs_cmb[[1]][[1]],NULL,matrs_cmb[[1]][[2]],NULL,
                          nrow = 1, rel_widths = c(0,0.5,-0.02,0.5,0.01), axis = 'tblr', align = 'hv')

plt_map_row <- plot_grid(NULL,pca_maps$ACER_ap$plt$map,NULL,pca_maps$TRACE_apsb$plt$map,NULL,
                         rel_widths = c(0.0,0.5,-0.02,0.5,0.0),
                         align = 'hv', axis = 'tblr', nrow = 1)

plt_cmb <- plot_grid(plt_matr_row,
                     NULL,
                     plt_map_row,
                     ncol = 1,
                     rel_heights = c(0.52,0.1,0.38))

for (i in 1:length(cntpltx)) {
  nm <- names(cntpltx)[i]
  plt_cmb <- plt_cmb + 
    draw_plot(cntpltx[[nm]], x  = cntlx_o[[nm]]-0.004, y = 0.39, width = 0.09, height = 0.09) + 
    draw_plot(cntpltx[[nm]], x  = cntlx_o[[nm]]+0.4845-0.004, y = 0.39, width = 0.09, height = 0.09)
}

for (i in 1:length(cntplty)) {
  nm <- names(cntpltx)[i]
  plt_cmb <- plt_cmb + 
    draw_plot(cntpltx[[nm]], x  = 0.024, y = cntly_o[[nm]], width = 0.09, height = 0.09) 
}

plt_cmb_a <- plt_cmb + 
  draw_plot_label('A', hjust = 0, size = GLOBAL_FONT_SIZE + 16) + 
  draw_plot_label('B', hjust = 0, y = 0.39, size = GLOBAL_FONT_SIZE + 16) + 
  draw_plot(matrs_cmb$leg1, 0.3949, 0.7538, 0.2, 0.2) + 
  draw_grob(matrs_cmb$leg2, 0.1836, 0.71, 0.23, 0.23) + 
  draw_grob(pca_maps$TRACE_apsb$clb, 0.381, 0.48, 0.23, 0.23) + 
  draw_grob(grid::rectGrob(gp = grid::gpar(fill = NA, lineheight = 10, fontsize = 20)),
            x = 0.42, y = 0.53, width = 0.1514, height = 0.42) + 
  draw_plot(pca_maps$ACER_ap$plt$ts, 0.008,0.028,0.13,0.14) + 
  draw_plot(pca_maps$TRACE_apsb$plt$ts, 0.496,0.028,0.13,0.14) + 
  draw_plot_label('ACER AP', hjust = 0, vjust = 0, x = 0.1, y = 0.97, size = GLOBAL_FONT_SIZE + 16) + 
  draw_plot_label('TraCE TSF', hjust = 0, vjust = 0, x = 0.32, y = 0.486, size = GLOBAL_FONT_SIZE + 16) + 
  draw_plot_label('TraCE PR', hjust = 0, vjust = 0, x = 0.585, y = 0.97, size = GLOBAL_FONT_SIZE + 16) +
  draw_plot_label('TraCE TS', hjust = 0, vjust = 0, x = 0.805, y = 0.486, size = GLOBAL_FONT_SIZE + 16) + 
  draw_plot_label('ACER AP', hjust = 0, x = 0.1, y = 0.39, size = GLOBAL_FONT_SIZE + 16) + 
  draw_plot_label('TraCE TSF', hjust = 0, x = 0.585, y = 0.39, size = GLOBAL_FONT_SIZE + 16)

ggsave(plt_cmb_a,
       filename = file.path(DIR_FIGURES, 'CLF_NDG_PCAfilt_ACERTRACE_raw_updnorm_correctedNA_woanom.pdf'),
       width = 64, height = 40, units = 'cm', device = cairo_pdf)
ggsave(plt_cmb_a,
       filename = file.path(DIR_FIGURES, 'CLF_NDG_PCAfilt_ACERTRACE_raw_updnorm_correctedNA_woanom.png'),
       width = 64, height = 40, units = 'cm', dpi = 100)


## Figure S20: TraCE millennial time scales ----
cntpltx <- plot_continents_row('horizontal', 'aligned', list_only = FALSE)
cntplty <- plot_continents_row('vertical', 'aligned', list_only = FALSE)

matrs_cmb <- plot_nw_clf_and_ndg_matr1x2_nosig(ACERtrc_clp_mil %>% rename_if(., str_detect(names(.), 'TRACE_apsb'), ~ str_replace(., 'TRACE_apsb', 'TRACE_tsf')) %>% 
                                                     rename_if(., str_detect(names(.), 'TRACE'), ~ str_replace(., 'TRACE', 'TraCE')),
                                                   ACERtrc_ndg_mil %>% filter(site_id <= 100) %>% rename_if(., str_detect(names(.), 'TRACE_apsb'), ~ str_replace(., 'TRACE_apsb', 'TRACE_tsf')) %>% 
                                                     rename_if(., str_detect(names(.), 'TRACE'), ~ str_replace(., 'TRACE', 'TraCE')),
                                                   tscales = 'mil',
                                                   suffix_pairs = trc_suff[2], legend = T, leg_option = 'not_shared_internal')[[1]]

matrs_cmb <- matrs_cmb + 
  draw_plot(cntpltx, x = -0.0025, y = 0.1, width = 0.48, height = 0.48, scale = 1) + 
  draw_plot(cntpltx, x = 0.4975, y = 0.1, width = 0.48, height = 0.48, scale = 1) + 
  draw_plot(cntplty, x = -0.197, y = 0.24, width = 0.6, height = 0.49, scale = 1)
ggsave(plot = matrs_cmb, filename = file.path(DIR_FIGURES, 'CLF_NDG_ACERTRACE_mil2-8ka_updnorm_correctedNA_woanom_compare.pdf'), 
       width = 54, height = 45, units = 'cm', device = cairo_pdf)


## Figure S21: HadCM3 ----
bbctrc_comb_suff <- list(
  'raw' = list(c('HadCM3_vc','HadCM3_tsf'), c('HadCM3_pr','HadCM3_ts'))
)

# save data for quicker changes to plot
bbctrc_comb_suff <- list('raw' = list(c('HadCM3_tsf'), c('HadCM3_pr','HadCM3_ts')))

matrs_cmb <- plot_nw_clf_and_ndg_matr1x2_nosig(ACERhcm_clp_raw %>% 
                                                 rename_if(., str_detect(names(.), 'BBC_ap'), ~ str_replace(., 'BBC_ap', 'BBC_tsf')) %>% 
                                                 rename_if(., str_detect(names(.), 'BBC'), ~ str_replace(., 'BBC', 'HadCM3')),
                                               ACERhcm_ndg_raw %>% 
                                                 rename_if(., str_detect(names(.), 'BBC_ap'), ~ str_replace(., 'BBC_ap', 'BBC_tsf')) %>% 
                                                 rename_if(., str_detect(names(.), 'BBC'), ~ str_replace(., 'BBC', 'HadCM3')),
                                               suffix_pairs = bbctrc_comb_suff, legend = T, leg_option = 'not_shared_internal')[[1]] # 'center'

matrs_cmb <- matrs_cmb + 
  draw_plot(cntpltx, x = -0.0025, y = 0.1, width = 0.48, height = 0.48, scale = 1) + 
  draw_plot(cntpltx, x = 0.4975, y = 0.1, width = 0.48, height = 0.48, scale = 1) + 
  draw_plot(cntplty, x = -0.197, y = 0.24, width = 0.6, height = 0.49, scale = 1)
ggsave(plot = matrs_cmb, filename = file.path(DIR_FIGURES, 'CLF_NDG_ACERBRIDGE_raw_updnorm_correctedNA_woanom_compare.pdf'), 
       width = 54, height = 45, units = 'cm', device = cairo_pdf)



## Figure S22: LOVECLIM ----
lc_suff <- list(
  'raw' = list(c('ACER_ap'),c('LOVECLIM_pr','LOVECLIM_ts')),
  'mil' = list(c('ACER_ap'),c('LOVECLIM_pr','LOVECLIM_ts'))
)

lc_matrs_cmb <- plot_nw_clf_and_ndg_matr2x2_nosig(ACERlc_clp_raw %>% rename_if(., str_detect(names(.), 'LC_'), ~ str_replace(., 'LC_', 'LOVECLIM_')),
                                                  ACERlc_clp_mil %>% rename_if(., str_detect(names(.), 'LC_'), ~ str_replace(., 'LC_', 'LOVECLIM_')),
                                                  ACERlc_ndg_raw %>% rename_if(., str_detect(names(.), 'LC_'), ~ str_replace(., 'LC_', 'LOVECLIM_')),
                                                  ACERlc_ndg_mil %>% rename_if(., str_detect(names(.), 'LC_'), ~ str_replace(., 'LC_', 'LOVECLIM_')),
                                                  suffix_pairs = lc_suff[c(1,2)], legend = T, leg_option = 'internal')[[1]] # 'center'

lc_matrs_cmb <- lc_matrs_cmb + 
  draw_plot(cntpltx, x = -0.0025, y = -0.145, width = 0.48, height = 0.48, scale = 1) + 
  draw_plot(cntpltx, x = 0.4975, y = -0.145, width = 0.48, height = 0.48, scale = 1) + 
  draw_plot(cntplty, x = -0.195, y = 0.495, width = 0.6, height = 0.48, scale = 1) + 
  draw_plot(cntplty, x = -0.195, y = -0.005, width = 0.6, height = 0.48, scale = 1)
ggsave(plot = lc_matrs_cmb, filename = file.path(DIR_FIGURES, 'CLF_NDG_ACERLCLIM_noVC_rawmil2-8ka_correctedNA_woanom_compare.pdf'), 
       width = 54, height = 45, units = 'cm', device = cairo_pdf)


## Figure S24: emulated TraCE ----
suff_em <- list(
  'raw' = list(c('EM TSF_all', 'EM TSF_ts'),c('EM TSF_co2','EM TSF_pr'))
)

matrs_cmb <- plot_nw_clf_and_ndg_matr1x2_nosig(ACERtrc_em_clp %>% rename_if(., str_detect(names(.), 'force'), ~ str_replace(., 'force', 'EM TSF_')),
                                                   ACERtrc_em_ndg %>% rename_if(., str_detect(names(.), 'force'), ~ str_replace(., 'force', 'EM TSF_')),
                                                   tscales = 'raw',
                                                   suffix_pairs = suff_em, legend = T, leg_option = 'not_shared_internal')[[1]]

matrs_cmb <- matrs_cmb + 
  draw_plot(cntpltx, x = -0.0025, y = 0.1, width = 0.48, height = 0.48, scale = 1) + 
  draw_plot(cntpltx, x = 0.4975, y = 0.1, width = 0.48, height = 0.48, scale = 1) + 
  draw_plot(cntplty, x = -0.197, y = 0.24, width = 0.6, height = 0.49, scale = 1)
ggsave(plot = matrs_cmb, filename = file.path(DIR_FIGURES, 'CLF_NDG_TRACE_em_updnorm_correctedNA_woanom_compare.pdf'), 
       width = 54, height = 45, units = 'cm', device = cairo_pdf)


## Figure S23: emulated HadCM3 ----
matrs_cmb <- plot_nw_clf_and_ndg_matr1x2_nosig(ACERhcm_em_clp %>% rename_if(., str_detect(names(.), 'force'), ~ str_replace(., 'force', 'EM TSF_')),
                                                   ACERhcm_em_ndg %>% rename_if(., str_detect(names(.), 'force'), ~ str_replace(., 'force', 'EM TSF_')),
                                                   tscales = 'raw',
                                                   suffix_pairs = suff_em, legend = T, leg_option = 'not_shared_internal')[[1]]

matrs_cmb <- matrs_cmb + 
  draw_plot(cntpltx, x = -0.0025, y = 0.1, width = 0.48, height = 0.48, scale = 1) + 
  draw_plot(cntpltx, x = 0.4975, y = 0.1, width = 0.48, height = 0.48, scale = 1) + 
  draw_plot(cntplty, x = -0.197, y = 0.24, width = 0.6, height = 0.49, scale = 1)
ggsave(plot = matrs_cmb, filename = file.path(DIR_FIGURES, 'CLF_NDG_HadCM3_em_updnorm_correctedNA_woanom_compare.pdf'), 
       width = 54, height = 45, units = 'cm', device = cairo_pdf)


## Figure S14: PCA results for TraCE TS and PR for the raw signals ----
### use data and plotting code prepared above similar to Figure 5
plt_map_row <- plot_grid(NULL,pca_maps$TRACE_ts$plt$map,NULL,pca_maps$TRACE_pr$plt$map,NULL,
                         rel_widths = c(0.0,0.5,-0.02,0.5,0.0),
                         align = 'hv', axis = 'tblr', nrow = 1)

plt_cmb <- plot_grid(NULL,
                     plt_map_row,
                     NULL,
                     ncol = 1,
                     rel_heights = c(0.01,0.38,0.05))

plt_cmb_a <- plt_cmb + 
  draw_grob(pca_maps$TRACE_pr$clb, 0.381, -0.04, 0.23, 0.23) + 
  draw_grob(grid::rectGrob(gp = grid::gpar(fill = NA, lineheight = 10, fontsize = 20)),
            x = 0.3765, y = 0.001, width = 0.24, height = 0.15) + 
  draw_plot(pca_maps$TRACE_ts$plt$ts, 0.008,0.165,0.13,0.335) + 
  draw_plot(pca_maps$TRACE_pr$plt$ts, 0.496,0.165,0.13,0.335) + 
  draw_plot_label('TraCE TS', hjust = 0, vjust = 0, x = 0.1, y = 0.94, size = GLOBAL_FONT_SIZE + 16) + 
  draw_plot_label('TraCE PR', hjust = 0, vjust = 0, x = 0.585, y = 0.94, size = GLOBAL_FONT_SIZE + 16)

ggsave(plt_cmb_a,
       filename = file.path(DIR_FIGURES, 'PCAfilt_TRACEtspr_raw_updnorm_correctedNA_woanom.pdf'),
       width = 64, height = 17, units = 'cm', device = cairo_pdf)
ggsave(plt_cmb_a,
       filename = file.path(DIR_FIGURES, 'PCAfilt_TRACEtspr_raw_updnorm_correctedNA_woanom.png'),
       width = 64, height = 17, units = 'cm', dpi = 100)

