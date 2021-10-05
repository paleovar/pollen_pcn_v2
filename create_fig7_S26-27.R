# Script to compute 
#   Figure 7 from the manuscript
#   -  Impact of several site-wise statistics on node degree and the fraction of positive node degrees
#      for the raw signals of ACER AP paleoclimate network.
#   - sensitivity of link density on site type (terrestrial vs marine records; 
#     analogous to Figure 4 from the manuscript)
#   Supplementary Figure S27
#   -  Same as Figure 7 but for ACER AP network based on millennial-scale variability
#   Supplementary Figure S26
#   -  Map of covered time period within 22-6 ka BP window (is based on time coverage data used for
#      Figures 7 and S27 as well)

#DONE, NOT TESTED YET

# load/compute network measures (need node degrees only here)
source('main_pcn_meas.R')
source('main_explained_variance.R')
source('main_pcn_sensitivity.R')

# prepare data to investigate relations
## NDG ~ resolution
## NDG ~ explained variance
## NDG ~ covered time span
## NDG ~ mean distance to other (used) records
## NDG ~ mean age uncertainty (if available)
# for raw signal and signal filtered for millennial time scales ----
trc_suff_dat <- list(
  'raw' = list(c('ACER_ap','TRACE_tsf'), c('TRACE_pr','TRACE_ts')),
  'mil' = list(c('ACER_ap','TRACE_tsf'), c('TRACE_pr','TRACE_ts')),
  'unused' = c()
)

ACERhere <- PCNdata() %>% 
  mix_sample_dating_and_err() %>% 
  harmonized_pollen_to_db() %>% 
  arboreal_pollen_to_db()

get_data_nw_ndg <- function(datndg_raw, datndg_mil, suffix_pairs,
                            nw_meas_ndg = 'node_degree',
                            regions = REGIONS) {
  
  plts <- lapply(c('raw', 'mil'), function(tscale) {
    datndg <- get(paste0('datndg_',tscale))
    plts <- lapply(1:length(suffix_pairs[[tscale]]), function(i){
      suffixes <- suffix_pairs[[tscale]][[i]]
      print(suffixes)
      
      # prepare ndg data
      datndg <- datndg %>% select(window, site_id,
                                  !!!syms(mapply(function(nm, pr) return(paste0(pr, nm)), pr = c(paste0(nw_meas_ndg, c('','.p','.n'), '.')), #, 'diff_sig.'
                                                 nm = sort(rep(x = suffixes, 3)), USE.NAMES = FALSE)))
      datndg <- lapply(c(paste0(nw_meas_ndg, c('', '.p','.n'))), function(c) { # , 'diff_sig'
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
        spread(ndg_type,ndg_val) %>% 
        mutate(frac_p = case_when(p != 0 & n != 0 ~ p/sum,
                                  p == 0 & n != 0 ~ 0,
                                  p != 0 & n == 0 ~ 1,
                                  T ~ 0.5),
               frac_n = case_when(p != 0 & n != 0 ~ n/sum,
                                  p == 0 & n != 0 ~ 1,
                                  p != 0 & n == 0 ~ 0,
                                  T ~ 0.5)) %>% 
        gather(key = 'sgn', value = 'single',p,n) %>% 
        mutate(type = str_replace_all(type, '_', ' ')) %>% ungroup() %>%
        mutate(sgn = case_when(sgn == 'n' ~ '-',
                               sgn == 'p' ~ '+', 
                               TRUE ~ sgn))
      
      datndg$type <- factor(datndg$type, levels = str_replace_all(suffixes, '_', ' '))
      
      return(datndg)
    })
    out <- list()
    out[[tscale]] <- plts
    return(out)
  })
  return(plts)
}

datACERtrc_ndg <- get_data_nw_ndg(ACERtrc_ndg_raw %>% rename_if(., str_detect(names(.), 'TRACE_apsb'), ~ str_replace(., 'TRACE_apsb', 'TRACE_tsf')),
                                  ACERtrc_ndg_mil %>% filter(site_id <= 100) %>% rename_if(., str_detect(names(.), 'TRACE_apsb'), ~ str_replace(., 'TRACE_apsb', 'TRACE_tsf')),
                                  suffix_pairs = trc_suff_dat[c(1,2)]) %>% 
  lapply(., function(x) x[[1]] %>% bind_rows()) %>% 
  setNames(c('raw', 'mil'))

# explained variance
ndg_scatterplt_dat_raw <- inner_join(expl_var_raw, datACERtrc_ndg$raw)
ndg_scatterplt_dat_mil <- inner_join(expl_var_mil, datACERtrc_ndg$mil)

# pos, neg vals
ndg_scatterplt_dat_raw$sgn <- factor(ndg_scatterplt_dat_raw$sgn, levels = c('+', '-'))
ndg_scatterplt_dat_mil$sgn <- factor(ndg_scatterplt_dat_mil$sgn, levels = c('+', '-'))

# resolution
res_data_6_22 <- compute_age_interval_stats(samples = ACERhere$sample_dating, ranges = list(list(start = 6000, stop = 22000))) %>% 
  select(site_id, med_smp_res)

ndgres_scatterplt_dat_raw <- inner_join(res_data_6_22, datACERtrc_ndg$raw)
ndgres_scatterplt_dat_mil <- inner_join(res_data_6_22, datACERtrc_ndg$mil)

# pos, neg vals
ndgres_scatterplt_dat_raw$sgn <- factor(ndgres_scatterplt_dat_raw$sgn, levels = c('+', '-'))
ndgres_scatterplt_dat_mil$sgn <- factor(ndgres_scatterplt_dat_mil$sgn, levels = c('+', '-'))

# length of covered time
tcover_dat <- compute_age_interval_stats(samples = ACERhere$sample_dating, ranges = list(list(start = 6000, stop = 22000))) %>% 
  mutate(tcovered = oldest_dp_in_interv - youngest_dp_in_interv) %>% 
  select(site_id, tcovered)

ndgtcover_scatterplt_dat_raw <- inner_join(tcover_dat, datACERtrc_ndg$raw)
ndgtcover_scatterplt_dat_mil <- inner_join(tcover_dat, datACERtrc_ndg$mil)

# mean distance to other records (only used ones)
p4sCRS_dist <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
ACERdist_raw <- ACERhere$sites %>% 
  filter(site_id %in% datACERtrc_ndg$raw$site_id) %>% 
  arrange(site_id)
ACERdist_raw_site_id <- ACERdist_raw %>% select(site_id)
ACERdist_raw <- ACERdist_raw %>% 
  select(long, lat) %>% 
  SpatialPoints()
proj4string(ACERdist_raw) <- p4sCRS_dist

ACERdist_raw <- raster::pointDistance(ACERdist_raw,ACERdist_raw,lonlat = T,allpairs = T) #%>% apply(.,1,sort)
dimACERdist <- dim(ACERdist_raw)[1]
ACERdist_raw <- ACERdist_raw %>% apply(.,1,mean) #ACERdist_raw[2:dim(ACERdist_raw)[1],] %>% apply(.,1,mean)
ACERdist_raw <- ACERdist_raw/1000. * dimACERdist / (dimACERdist - 1) # correct mean for self-distance
ACERdist_raw <- tibble(site_id = ACERdist_raw_site_id$site_id, dist = ACERdist_raw)

ACERdist_mil <- ACERhere$sites %>% 
  filter(site_id %in% datACERtrc_ndg$mil$site_id) %>% 
  arrange(site_id)
ACERdist_mil_site_id <- ACERdist_mil %>% select(site_id)
ACERdist_mil <- ACERdist_mil %>% 
  select(long, lat) %>% 
  SpatialPoints()
proj4string(ACERdist_mil) <- p4sCRS_dist

ACERdist_mil <- raster::pointDistance(ACERdist_mil,ACERdist_mil,lonlat = T,allpairs = T) #%>% apply(.,1,sort)
dimACERdist <- dim(ACERdist_mil)[1]
ACERdist_mil <- ACERdist_mil %>% apply(.,1,mean) #ACERdist_mil[2:dim(ACERdist_mil)[1],] %>% apply(.,1,mean)
ACERdist_mil <- ACERdist_mil/1000.* dimACERdist / (dimACERdist - 1) # correct mean for self-distance
ACERdist_mil <- tibble(site_id = ACERdist_mil_site_id$site_id, dist = ACERdist_mil)

ndgdist_scatterplt_dat_raw <- inner_join(ACERdist_raw, datACERtrc_ndg$raw)
ndgdist_scatterplt_dat_mil <- inner_join(ACERdist_mil, datACERtrc_ndg$mil)

# mean age uncertainty (if available, otherwise is NA and ignored for the plot)
uncert_data <- ACERhere$sample_dating %>% 
  select(site_id, sample_id, mixed_age_min, mixed_age, mixed_age_max) %>% 
  filter(mixed_age >= 6000. & mixed_age <= 22000.) %>% 
  filter(mixed_age_max > 0 & mixed_age_min > 0) %>% 
  mutate(uncert = (mixed_age_max - mixed_age_min) / 2.)%>% 
  group_by(site_id) %>% 
  summarise(uncert = mean(uncert,na.rm = T)) %>% 
  filter(uncert > 0.)

ndguct_scatterplt_dat_raw <- inner_join(uncert_data, datACERtrc_ndg$raw)
ndguct_scatterplt_dat_mil <- inner_join(uncert_data, datACERtrc_ndg$mil)

# combine data for NDG scatter plots of ACER TSF ----
dat_ndg_scatter <- bind_rows(plyr::join_all(list(datACERtrc_ndg$raw %>% filter(type == 'ACER ap') %>% select(site_id, sum, frac_p) %>% distinct(), 
                                                 expl_var_raw, 
                                                 res_data_6_22 %>% select(site_id, med_smp_res) %>% distinct(), 
                                                 tcover_dat, 
                                                 ACERdist_raw, 
                                                 uncert_data)) %>% 
                               mutate(tscale = 'RAW'),
                             plyr::join_all(list(datACERtrc_ndg$mil %>% filter(type == 'ACER ap') %>% select(site_id, sum, frac_p) %>% distinct(), 
                                                 expl_var_mil, 
                                                 res_data_6_22 %>% select(site_id, med_smp_res) %>% distinct(), 
                                                 tcover_dat, 
                                                 ACERdist_raw, 
                                                 uncert_data)) %>% 
                               mutate(tscale = 'MIL')) %>% 
  as_tibble() %>% 
  mutate(expl_var = 100. * expl_var, tcovered = tcovered / 1000., dist = dist / 1000., med_smp_res = med_smp_res / 100., uncert = uncert / 100.) %>% 
  rename('mean dist. [1000 km]' = dist, 'expl. var. [%]' = expl_var, 'med. ISR [100 yrs]' = med_smp_res, 'time covered [Kyr]' = tcovered, 'mean unct. [100 yrs]' = uncert) %>% 
  gather(key = 'measure_name', value = 'measure_val', -site_id, -sum, -frac_p, -tscale) %>% 
  rename(NDG = sum, 'frac. pos. NDG' = frac_p) %>% 
  gather(key = 'NDG_type', value = 'NDG_val', NDG, !!sym('frac. pos. NDG'))

dat_ndg_scatter$NDG_type <- factor(dat_ndg_scatter$NDG_type, levels = c('NDG', 'frac. pos. NDG'))
dat_ndg_scatter$measure_name <- factor(dat_ndg_scatter$measure_name, levels = c('med. ISR [100 yrs]', 'expl. var. [%]', 'mean unct. [100 yrs]', 'time covered [Kyr]', 'mean dist. [1000 km]'))

# compute linear coefficient of determination, plot and save scatter plots
# Figure 7, panel A: impact of statistics on NDG for raw signal ----
dat_ndg_scatter_raw <- dat_ndg_scatter %>% 
  filter(tscale == 'RAW')

dat_rsq_raw <- dat_ndg_scatter_raw %>% 
  group_by(tscale, NDG_type, measure_name) %>% 
  filter(!is.na(NDG_val) & !is.na(measure_val)) %>% 
  summarise(rsq = cor(NDG_val, measure_val) ^ 2) %>% 
  rowwise() %>% 
  mutate(rsq = paste0("R^2*' = '*",as.character(round(rsq, digits = 3)))) %>% # if used with non-embedded fonts, i.e. w/o cairo_pdf below: paste0("R^2~'='~",as.character(round(rsq, digits = 3)))
  ungroup()

plt <- ggplot(dat_ndg_scatter_raw, aes(x = measure_val, y = NDG_val)) + 
  geom_point(color = GLOBAL_BLUE_DDARK, shape = 21, fill = GLOBAL_BLUE_DARK, size = 2.4) + 
  geom_label(dat = dat_rsq_raw, x = Inf, y = Inf, mapping = aes(label = rsq), parse = T, hjust = 1.04, vjust = 1.12, size = 4.5) + 
  facet_grid(rows = vars(NDG_type), cols = vars(measure_name), switch = 'both', scales = 'free') + 
  global_title_and_axis() + 
  theme(axis.title = element_blank(),
        axis.text = element_text(size = GLOBAL_FONT_SIZE),
        strip.text.x = element_text(hjust = 0.5, size = GLOBAL_FONT_SIZE),
        strip.text.y = element_text(size = GLOBAL_FONT_SIZE),
        strip.background.y = element_blank(),
        strip.placement = 'outside',
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(r = 9, t = 4))

# to save separately, uncomment these lines
#ggsave(plot = plt, width = 33, height = 14, units = 'cm', 
#       filename = file.path(DIR_FIGURES, 'NDGscatterplt_ACERapraw.pdf'), device = cairo_pdf)
#ggsave(plot = plt, width = 33, height = 14, units = 'cm', 
#       filename = file.path(DIR_FIGURES, 'NDGscatterplt_ACERapraw.png'))

# Figure 7, panel B: sensitivity of link density to site type ----
ld_lims_stype <- round(max(abs(ACERsens_stype_err_sgn$link_dens)),3) + 0.001
ACERsens_stat <- ACERsens_stype_err_sgn %>% mutate(signal = paste0(signal, '\n(n = ',nnode,')')) %>% group_by(signal,ld_sgn) %>% 
  summarise(median = median(link_dens), qtl = quantile(link_dens, probs = 0.05), qth = quantile(link_dens, probs = 0.95))
ld_rank_lcl <- ld_rank %>% full_join(ACERsens_stype_err_sgn %>% distinct(signal,nnode)) %>% 
  mutate(nnode = if_else(is.na(nnode), ACERsens_stype_err_sgn %>% distinct(signal,nnode) %>% .$nnode %>% sum(), nnode)) %>% 
  mutate(signal = paste0(signal, '\n(n = ',nnode,')'))

plt_ld_rank <- ggplot(mapping = aes(x = signal)) + #color = ld_sgn, 
  geom_col(data = ld_rank_lcl, mapping = aes(y = link_dens, fill = ld_sgn), width = 0.8) + 
  geom_errorbar(data = ACERsens_stat, mapping = aes(ymin = qtl, ymax = qth), color = 'black', width = 0.07, size = 0.6) +
  geom_point(data = ACERsens_stat, mapping = aes(y = median), color = 'black', size = 2.5) + 
  scale_fill_manual(values = rev(GLOBAL_PN_COLORS), guide = guide_legend('correlation sign', title.position = 'left')) + 
  scale_color_manual(values = rev(GLOBAL_PN_COLORS), guide = FALSE) + 
  scale_y_continuous(limits = c(-1*ld_lims_stype,ld_lims_stype)) + 
  labs(y = 'LD') + 
  global_title_and_axis() + 
  theme(legend.direction = 'horizontal', axis.title.x = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = 'bottom', axis.text = element_text(size = GLOBAL_FONT_SIZE))
# to save separately, uncomment these lines
#ggsave(plt_ld_rank, filename = file.path(DIR_FIGURES, 'ACER_AP_LD_vs_stype.pdf'),
#       width = 16, height = 14, units = 'cm')

# Figure 7: combine panels and save plot ----
# (geom_point warnings arise from NAs in uncert data if uncertainties are not provided)
plt_cmb <- plot_grid(plt,
                     NULL,
                     plot_grid(NULL,plt_ld_rank,NULL,nrow = 3, rel_heights = c(-0.0035,0.8,-0.005)),
                     ncol = 3,
                     rel_widths = c(0.72,0.0,0.28),
                     axis = 'tb', 
                     labels = c('A','','B'), label_size = GLOBAL_FONT_SIZE + 4)

ggsave(plt_cmb, filename = file.path(DIR_FIGURES, 'NDGscatterplt_LDvsstype_ACERapraw.pdf'),
       width = 46, height = 14, units = 'cm')


# Figure S27: impact of statistics on NDG for millennial scale signal ----
dat_ndg_scatter_mil <- dat_ndg_scatter %>% 
  filter(tscale == 'MIL')

dat_rsq_mil <- dat_ndg_scatter_mil %>% 
  group_by(tscale, NDG_type, measure_name) %>% 
  filter(!is.na(NDG_val) & !is.na(measure_val)) %>% 
  summarise(rsq = cor(NDG_val, measure_val) ^ 2) %>% 
  rowwise() %>% 
  mutate(rsq = paste0("R^2*' = '*",as.character(round(rsq, digits = 3)))) %>% # if used with non-embedded fonts, i.e. w/o cairo_pdf below: paste0("R^2~'='~",as.character(round(rsq, digits = 3)))
  ungroup()


plt <- ggplot(dat_ndg_scatter_mil, aes(x = measure_val, y = NDG_val)) + 
  geom_point(color = GLOBAL_BLUE_DDARK, shape = 21, fill = GLOBAL_BLUE_DARK, size = 2.4) + 
  geom_label(dat = dat_rsq_mil, x = Inf, y = Inf, mapping = aes(label = rsq), parse = T, hjust = 1.04, vjust = 1.12, size = 4.5) + 
  facet_grid(rows = vars(NDG_type), cols = vars(measure_name), switch = 'both', scales = 'free') + 
  global_title_and_axis() + 
  theme(axis.title = element_blank(),
        axis.text = element_text(size = GLOBAL_FONT_SIZE),
        strip.text.x = element_text(hjust = 0.5, size = GLOBAL_FONT_SIZE),
        strip.text.y = element_text(size = GLOBAL_FONT_SIZE),
        strip.background.y = element_blank(),
        strip.placement = 'outside',
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(r = 9, t = 4))

ggsave(plot = plt, width = 33, height = 14, units = 'cm', 
       filename = file.path(DIR_FIGURES, 'NDGscatterplt_ACERapmil.pdf'), device = cairo_pdf)
ggsave(plot = plt, width = 33, height = 14, units = 'cm', 
       filename = file.path(DIR_FIGURES, 'NDGscatterplt_ACERapmil.png'))


# Figure S26: time coverage map ----
# define plot function
plot_ACER_sites_tcovered_map <- function(dat) {
  projection <- 'robinson'
  zoom <- c(-180, -60, 180, 75) 
  graticules <- 30
  
  sites <- arrange(dat, site_id)
  sites_transformed <- as_tibble(project(cbind(sites$long, sites$lat), proj = GLOBAL_CRS[[projection]])) %>% 
    rename(long_new = V1, lat_new = V2) %>% 
    bind_cols(select(sites, site_id))
  
  sites <- sites %>% 
    full_join(sites_transformed,
              by = 'site_id') %>% 
    select(-lat, -long) %>% 
    rename(long = long_new, lat = lat_new)
  if (!is.null(zoom)) {
    yaxislclip <- zoom[2]; xaxislclip <- zoom[1]; xaxisrclip <- zoom[3]
  } else {
    yaxislclip <- -90; xaxislclip <- -180; xaxisrclip <- 180
  }
  
  data <- inner_join(sites, dat %>% select(site_id, tcovered))
  
  map_plot <- shapefile_map_plot_(ggobj = ggplot(), projection = projection, grat_labs = list(x = c(-180, -90, 0, 90, 180), y = c(-90, -60, -30, 0, 30, 60, 90)),
                                  graticules = graticules, 
                                  lgm_ice = F, cb_friendly = TRUE, zoom = zoom,
                                  yaxislclip = yaxislclip, xaxisrclip = xaxisrclip, xaxislclip = xaxislclip)
  map_plot <- map_plot + 
    geom_point(data = data,
               mapping = aes(x = long, y = lat,
                             fill = tcovered / 1000.),
               size = 5,
               shape = 21,
               alpha = I(0.6),
               stroke = 1) +
    scale_fill_fermenter(palette = 'Spectral', direction = 1, guide = guide_colorsteps(title = 'Covered time period\n(in 6-22 ka BP) [Kyr]', order = 3, show.limits = T, ticks = T, direction = 'horizontal',
                                                                                       barheight = 0.75, barwidth = 10, title.position = 'top')) + 
    theme(panel.background = element_blank(), 
          axis.text = element_blank(), 
          panel.grid = element_blank(), 
          axis.ticks = element_blank(), 
          panel.border = element_blank())
  
  return(map_plot)
}

# reuse above data on time coverage, plot, and save
plt <- plot_ACER_sites_tcovered_map(inner_join(ACERhere$sites %>% 
                                                 select(site_id, lat, long), tcover_dat) %>%
                                      inner_join(datACERtrc_ndg$raw %>%
                                                   select(site_id) %>% 
                                                   distinct()))
ggsave(file.path(DIR_FIGURES, 'map_tcovered.pdf'), 
       width = 38, height = 22, units = 'cm', dpi = 'print', limitsize = FALSE, device = cairo_pdf)
ggsave(file.path(DIR_FIGURES, 'map_tcovered.png'), 
       width = 38, height = 22, units = 'cm', dpi = 'print', limitsize = FALSE)

