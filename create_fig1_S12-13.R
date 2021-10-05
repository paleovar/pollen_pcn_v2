# Script to create Figure 1 from the manuscript
# - on a world map: individual mean ISR (6-22ka), explained variance, exemplary time series at ACER site locations
# and Figures S12 and S13 from the Supplement
# - explained variances of ACER AP when filtered for orbital scales (>= 8kyrs) and millennial time scales (8-22kyrs)

FILE_FIG1DATA_CACHE <- file.path(DIR_CACHE,'pre_processed_fig1_data.RData')

# DONE & TESTED

## compute variances explained by AP signal (can jump to next code section when using pre-computed .RData object)
source('main_explained_variance.R')

## check if cached data available ----
FIG1DATA_CACHED <- FALSE
if (file.exists(FILE_FIG1DATA_CACHE)) {
  message('Using pre-processed data, therefore skipping aggregating them from scratch')
  load(FILE_FIG1DATA_CACHE)
  FIG1DATA_CACHED <- TRUE
} else {
  message('no matching .RData file found, attempting to compute explained variance from scratch')
}

## ACER data extended with NGRIP and EPICA Dome C ----
ACERe <- PCNdata() %>% 
  harmonized_pollen_to_db() %>% 
  arboreal_pollen_to_db() %>% 
  extend_dataset_in_db(epica_ngrip_cores %>% mutate(age = age * 1e3), 
                           'd18O',
                           'arb_pollen_data', 
                           'arboreal_pollen', 
                           'AP',
                           'ACERarbpoll_NGECice',
                           'ACERarbpoll_NGECice_data',
                           'proxy', 
                           'proxy_id', 
                           'proxy_type', 
                           tibble(site_name = c('ngrip', 'epica'), 
                                  lat = c(75.1, -75.0), 
                                  long = c(-42.32, 125.0)))

## load TraCE LGM land sea mask and ICE-5G ice sheets as background (used within TRACE-21ka) ----
### TraCE lsm from PFT file
if (!FIG1DATA_CACHED) {
  nc <- nc_open(file.path(DIR_MODEL_DATA,'trace21ka/processed/trace.01-36.22000BP.PFTFRAC_decavg_400BCE.nc'))
  trace_pfts <- ncvar_get(nc, 'PFTFRAC')
  nc_close(nc)
  trace_lsm_lgm <- trace_pfts[1,,,1]
  trace_lsm_lgm[!is.nan(trace_lsm_lgm)] <- 1
  trace_lsm_lgm[is.nan(trace_lsm_lgm)] <- 0
  temp <- trace_lsm_lgm[1:48,]
  trace_lsm_lgm[1:48,] <- trace_lsm_lgm[49:96,]
  trace_lsm_lgm[49:96,] <- temp
  rm(temp)
  
  ### ICE-5G ice sheets
  nc <- nc_open('/stacywork/fredweasley/data/Peltier_ICE-5G/ice5g_v1.2_21.0k_3.75deg_remapbil.nc')
  ice5g_lgm <- ncvar_get(nc, 'sftgif')
  nc_close(nc)
  ice5g_lgm[ice5g_lgm == 100] <- 1
  ice5g_lgm[is.na(ice5g_lgm)] <- 0
  temp <- ice5g_lgm[1:48,]
  ice5g_lgm[1:48,] <- ice5g_lgm[49:96,]
  ice5g_lgm[49:96,] <- temp
  ice5g_evo <- ice5g_lgm 
  trace_lsm_lgm[ice5g_evo != 0] <- ice5g_evo[ice5g_evo != 0] + 2
}

### map plot as background in panel B (raw signals) ----
plt <- plot_ACER_sites_on_map2(sites = ACERe$sites, ACERe$sample_dating, expl_var_raw,
                               lgm_ice = F,
                               names = F, site_id = F,
                               ISR = T, ISR_stat = 'med', ISR_window = c(6000,22000),
                               save_plot = list(activate=F), 
                               legend_inside = F, #force_shapefiles = T, 
                               projection = 'robinson', 
                               bg = trace_lsm_lgm,#ice5g_lgm,
                               #zoom = c(-180, -60, 180, 75), 
                               graticules = 30,
                               cats = list(
                                 'smp_per_10k_in_interv' = c(0, 25, 28, Inf),
                                 'mean_smp_res' = c(0, 250, 500, Inf),
                                 'med_smp_res' = c(0, 250, 500, Inf),
                                 'ev_res_in_interv' = c(0, 0.6, 0.8, 0.9, Inf)
                               )
) + 
  theme(panel.background = element_blank(), 
        axis.text = element_blank(), 
        panel.grid = element_blank(), 
        axis.ticks = element_blank(), 
        panel.border = element_blank())

### save map only
ggsave(plot = plt, filename = file.path(DIR_FIGURES, 'fig1_map_dots_only_raw.pdf'), 
       width = 38, height = 22, units = 'cm', dpi = 'print', limitsize = FALSE)

### map plot (orbital signal) ----
plt_orb <- plot_ACER_sites_on_map2(sites = ACERe$sites, ACERe$sample_dating, expl_var_orb,
                                   lgm_ice = F,
                                   names = F, site_id = F,
                                   ISR = T, ISR_stat = 'med', ISR_window = c(6000,22000),
                                   save_plot = list(activate=F), 
                                   legend_inside = F, #force_shapefiles = T, 
                                   projection = 'robinson', 
                                   bg = NULL,#trace_lsm_lgm,
                                   #zoom = c(-180, -60, 180, 75), 
                                   graticules = 30,
                                   cats = list(
                                     'smp_per_10k_in_interv' = c(0, 25, 28, Inf),
                                     'mean_smp_res' = c(0, 250, 500, Inf),
                                     'med_smp_res' = c(0, 250, 500, Inf),
                                     'ev_res_in_interv' = c(0, 0.6, 0.8, 0.9, Inf)
                                   ),
                                   expl_var_tscale = 'orb'
) + 
  theme(panel.background = element_blank(), 
        axis.text = element_blank(), 
        panel.grid = element_blank(), 
        axis.ticks = element_blank(), 
        panel.border = element_blank())

### save map only
ggsave(plot = plt_orb, filename = file.path(DIR_FIGURES, 'fig1_map_dots_only_orb.pdf'), 
       width = 38, height = 22, units = 'cm', dpi = 'print', limitsize = FALSE)
ggsave(plot = plt_orb, filename = file.path(DIR_FIGURES, 'fig1_map_dots_only_orb.png'), 
       width = 38, height = 22, units = 'cm', dpi = 150, limitsize = FALSE)


### map plot (millennial signal) ----
plt_mil <- plot_ACER_sites_on_map2(sites = ACERe$sites, ACERe$sample_dating, expl_var_mil,
                                   lgm_ice = F,
                                   names = F, site_id = F,
                                   ISR = T, ISR_stat = 'med', ISR_window = c(6000,22000),
                                   save_plot = list(activate=F), 
                                   legend_inside = F, #force_shapefiles = T, 
                                   projection = 'robinson', 
                                   bg = NULL,
                                   #zoom = c(-180, -60, 180, 75), 
                                   graticules = 30,
                                   cats = list(
                                     'smp_per_10k_in_interv' = c(0, 25, 28, Inf),
                                     'mean_smp_res' = c(0, 250, 500, Inf),
                                     'med_smp_res' = c(0, 250, 500, Inf),
                                     'ev_res_in_interv' = c(0, 0.6, 0.8, 0.9, Inf)
                                   ),
                                   expl_var_tscale = 'mil'
) + 
  theme(panel.background = element_blank(), 
        axis.text = element_blank(), 
        panel.grid = element_blank(), 
        axis.ticks = element_blank(), 
        panel.border = element_blank())

### save map only
ggsave(plot = plt_mil, filename = file.path(DIR_FIGURES, 'fig1_map_dots_only_mil.pdf'), 
       width = 38, height = 22, units = 'cm', dpi = 'print', limitsize = FALSE)
ggsave(plot = plt_mil, filename = file.path(DIR_FIGURES, 'fig1_map_dots_only_mil.png'), 
       width = 38, height = 22, units = 'cm', dpi = 150, limitsize = FALSE)


## tseries plots as insets in panel B ----
### generate ACERtrc list first in main_simulation_data.R (incl ACER AP)
source('main_simulation_data.R')

### create object for inset time series data
yshifts <- tibble(data = c('ACER AP', 'TRACE PR', 'TRACE TS', 'TRACE TSF'), yshift = rev(seq(0,15,5)))
ACERtrcdat <- lapply(1:length(ACERtrc), function(i){
  nm <- names(ACERtrc)[i]
  dat <- ACERtrc[[nm]]$arb_pollen_data %>% 
    rename(!!sym(paste0(nm, '_pcnt_arb_pollen')) := pcnt_arb_pollen) %>% 
    select(-arb_pollen_data_id)
  dat
}) %>% 
  plyr::join_all(., type = 'inner', by = c('site_id', 'sample_id')) %>% 
  as_tibble() %>% 
  gather(key = 'data', value = 'pcnt_arb_pollen', paste0(names(ACERtrc), '_pcnt_arb_pollen')) %>% 
  mutate(data = str_replace(data,'_pcnt_arb_pollen','')) %>% 
  group_by(site_id, data) %>% 
  mutate(pcnt_arb_pollen = (pcnt_arb_pollen-mean(pcnt_arb_pollen))/sd(pcnt_arb_pollen)) %>% 
  inner_join(ACERtrc$ACER_ap$sample_dating %>% select(mixed_age, site_id, sample_id), by = c('site_id', 'sample_id')) %>% 
  ungroup() %>% 
  mutate(data = str_replace(data, '_', ' ')) %>% 
  filter(mixed_age >= 6000 & mixed_age <= 22000) %>% 
  filter(data != 'TRACE vc') %>% 
  rowwise() %>% 
  mutate(data = if_else(data == 'TRACE apsb', 'TRACE tsf', data)) %>% 
  mutate(data = toupper(data)) %>% 
  ungroup() %>% 
  inner_join(yshifts, by = 'data') %>% 
  mutate(pcnt_arb_pollen = pcnt_arb_pollen + yshift) %>% 
  mutate(data = str_replace(data, 'TRACE', 'TraCE'))

### time series plot function
plot_ts_single <- function(dat, grouper = 'data', leg = T){
  plt <- ggplot(data = dat,
                mapping = aes(x = mixed_age/1e3, y = pcnt_arb_pollen, color = !!sym(grouper))) + 
    geom_line(alpha = 0.75, size = 0.5) + 
    #scale_color_viridis_d(guide = guide_legend(title = 'Dataset', title.position = 'top', 
    #                                           override.aes = list(size = 4), nrow = 2)) +
    scale_color_manual(values = viridis(10)[c(3,5,7,9)], 
                       guide = guide_legend(title = 'Dataset', title.position = 'top', 
                                            override.aes = list(size = 4), nrow = 2)) +
    #scale_linetype_manual(guide = guide_legend(title = 'dataset', title.position = 'top'), values = c('solid', 'twodash', 'longdash', 'dotdash')) + 
    labs(x = 'Age / ka BP') + 
    #facet_grid(rows = vars(!!sym(grouper))) + 
    global_title_and_axis() + 
    coord_cartesian(xlim = c(6,22), ylim = c(-2.5,17)) + 
    theme(legend.position = {if(leg) {'bottom'} else {'none'}}, legend.direction = 'horizontal',
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          #panel.grid.major.x = element_line(color = 'grey50'),
          panel.spacing.y = unit(0, 'lines'), 
          panel.border = element_blank(), 
          plot.background = element_rect(color = 'black', size = 1, fill = 'grey98'), panel.background = element_rect(fill = 'grey98'),
          strip.background = element_blank(), strip.text.y = element_blank(),
          plot.margin = unit(c(1,1,1,1), 'mm'),
          title = element_text(size = GLOBAL_FONT_SIZE - 10), text = element_text(size = GLOBAL_FONT_SIZE - 4))
  return(plt)
}

### styling tweaks
tplt <- lapply(c(37, 69, 4,15,21,39), function(i) plot_ts_single(dat = filter(ACERtrcdat, site_id == i), leg = F))
tlgd <- (plot_ts_single(dat = filter(ACERtrcdat, site_id == 37), leg = T) + 
           theme(title = element_text(size = GLOBAL_FONT_SIZE), text = element_text(size = GLOBAL_FONT_SIZE))) %>% 
  get_legend(.)
tplx <- c(0.355,0.18,0.85,0.65,0.42,0.1)
tply <- c(0.56,0.18,0.55,0.18,0.18,0.53)
tplw <- rep(0.1,6)
tplh <- rep(0.18, 6)
tlbl <- c(37,69,4,15,21,39)
tlbx <- tplx + 0.01
tlby <- tply + 0.185
tseg <- tibble(
  id = c(rep(seq(1,6), 1)),
  x = rep(c(tplx[1:2] + tplw[1:2], tplx[3], tplx[4:6] + tplw[4:6]), 1), 
  y = c(tply + tplh), #tply
  xend = rep(c(0.535,0.3325,0.79,0.9025,0.5825,0.24), 1),
  yend = rep(c(0.725,0.3,0.635,0.295,0.45,0.69), 1)
)

### combine map, time series insets
cmbplt <- ggdraw(plt + theme(legend.position = 'none'), xlim = c(0,1), ylim = c(0,1))
for (i in 1:length(tplt)) {
  cmbplt <- cmbplt + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend), data = tseg %>% filter(id == i), lineend = "round") +
    draw_plot(tplt[i][[1]], x = tplx[i], y = tply[i], width = tplw[i], height = tplh[i]) + 
    annotate('text', label = 'ka BP', x = tplx[i], y = tply[i], hjust = -0.1, vjust = -0.5,
             fontface = GLOBAL_FONT_FACE_TITLE, size = GLOBAL_FONT_SIZE - 13.5) + 
    annotate('text', label = c('10','15','20'), x = tplx[i] + c(0.03,0.0555,0.084), y = tply[i], hjust = -0.1, vjust = -0.5,
             fontface = GLOBAL_FONT_FACE_TEXT, size = GLOBAL_FONT_SIZE - 13.5) + 
    annotate('label', label = tlbl[i], x = tlbx[i], y = tlby[i],
             fill = 'white', fontface = GLOBAL_FONT_FACE_TITLE, size = GLOBAL_FONT_SIZE - 11.5, label.size = 0.5) #draw_plot_label(tlbl[i], x = tlbx[i], y = tlby[i])
}
cmbplt <- cmbplt +
  annotate('label', label = c('NGRIP', 'EDC'), x = c(0.4425, 0.695), y = c(0.905, 0.105), hjust = c(0,1),
           fill = 'white', fontface = GLOBAL_FONT_FACE_TITLE, size = GLOBAL_FONT_SIZE - 11.5, label.size = 0.5)

lgd <- plot_grid(NULL,get_legend(plt),NULL,tlgd,nrow = 1, rel_heights = c(1,1,1,1), rel_widths = c(0.16,1,-0.475,1), align = 'h')

cmbplt <- plot_grid(cmbplt,
                    lgd,
                    rel_heights = c(1, 0.15), ncol = 1)

### for separate saving
#ggsave(file.path(DIR_FIGURES, 'fig1_map_ts_update21.pdf'), 
#       width = 38, height = 22, units = 'cm', dpi = 'print', limitsize = FALSE)


## data for panel A: global mean time series of the simulations & orbital and CO2 forcing ----
### load CO2 forcing from Köhler et al. 2017 ----
if (!FIG1DATA_CACHED) {
  CO2_forc <- paleodata_windowing(load_forcing_dataset()$co2_smooth,-22000,-6000)

### load model 2m temperature ----
  model <- c("trace",
             "loveclim",
             "hadcm3")
  sim_names <- c("TraCE-Full",
                 "LOVECLIM_full",
                 "BBC")
  sim_dir <- c(file.path(DIR_MODEL_DATA,"trace21ka/trace.01-36.22000BP.cam2.TREFHT.22000BP_decavg_400BCE.nc"),
               file.path(DIR_MODEL_DATA,"LOVECLIM/deglaciation_ns/annual_surface_temperature_at_2_meter.nc"))
  
  sim_data <- list()
  sim_gmst <- list()
  
  for (i in 1:length(sim_dir)) {
    if (model[i] == "trace") {
      nc <- nc_open(sim_dir[i])
      tas_sim <- ncvar_get(nc,"TREFHT")[c(1:96,1:96),,]
      lon_sim <- ncvar_get(nc,"lon")
      lon_sim <- c(lon_sim-360,lon_sim)
      lat_sim <- ncvar_get(nc,"lat")
      time_sim <- ncvar_get(nc,"time")*1000
    }
    if (model[i] == "loveclim") {
      nc <- nc_open(sim_dir[i])
      tas_sim <- ncvar_get(nc,"T2M")[c(34:64,1:65,2:32),,]
      lon_sim <- ncvar_get(nc,"LONN31_33")
      lon_sim <- c(lon_sim[34:64]-360,lon_sim,lon_sim[2:32]+360)
      lat_sim <- ncvar_get(nc,"LAT")
      decs <- seq(1,dim(tas_sim)[3],by=10)
      if (nc$dim$TAX$units == "YRS AFTER 18K BP") {
        #time_sim <- ncvar_get(nc,"TAX")-18000
        time_sim <- -18000 + decs - 5
      }
      if (nc$dim$TAX$units == "YRS AFTER 21K BP") {
        #time_sim <- ncvar_get(nc,"TAX")-21000
        time_sim <- -21000 + decs - 5
      }
      tas_sim <- sapply(decs, function(t) apply(tas_sim[,,t:t+9], c(1,2), mean)) %>% array(dim = c(dim(tas_sim)[c(1,2)],length(decs)))
    }
    sim_gmst[[i]] <- normalize(paleodata_windowing(zoo(apply(tas_sim,3,function(x) spatial_means(lon_sim,lat_sim,x)),order.by=time_sim),-22000,-6000),scale = FALSE)
    rm(tas_sim); nc_close(nc); gc()
  }
  
  # TraCE
  trace_ts <- sim_gmst[[1]]
  # LOVECLIM
  lclim_ts <- sim_gmst[[2]]
  
  # BRIDGE HadCM3
  nc <- nc_open("/stacydata/data/BRIDGE/bbc_all_triff_rec_dyn04/teii01.temp_mm_1_5m.monthly.nc")
  lat <- rev(ncvar_get(nc,"latitude"))
  lon <- ncvar_get(nc,"longitude")
  time <- ncvar_get(nc,"t")
  sim_names <- paste(read.csv("/stacydata/data/BRIDGE/BRIDGE.csv",header=TRUE)$bbc_all_triff_rev_dyn04)[1:62]
  bridge_time <- c(seq(120,80,by=-4),seq(78,22,by=-2),seq(21,0,by=-1))
  # Mean temperature
  bridge_ts <- array(0,dim=c(length(bridge_time),length(lon),length(lat)))
  for (i in 1:length(bridge_time)) {
    nc <- nc_open(paste("/stacydata/data/BRIDGE/bbc_all_triff_rec_dyn04/",sim_names[i],".temp_mm_1_5m.monthly.nc",sep=""))
    bridge_ts[i,,] <- apply(ncvar_get(nc,"temp_mm_1_5m")[,length(lat):1,9601:12000],c(1,2),mean)
  }
  bridge_lon <- lon
  bridge_lat <- lat
  bridge_ts <- paleodata_windowing(zoo(apply(bridge_ts,1,function(x) spatial_means(bridge_lon,bridge_lat,x)),order.by=-bridge_time*1000),-22000,-6000)
  
  ## load GMPRCP data from simulations
  # Trace
  nc <- nc_open(file.path(DIR_MODEL_DATA,"trace21ka/trace.01-36.22000BP.cam2.PRECT.22000BP_decavg_400BCE.nc"))
  pr_sim <- ncvar_get(nc,"PRECT")*60*60*24*30*12*1000
  lon_sim <- ncvar_get(nc,"lon")
  lat_sim <- ncvar_get(nc,"lat")
  time_sim <- ncvar_get(nc,"time")*1000
  trace_pr <- paleodata_windowing(zoo(apply(pr_sim,3,function(x) spatial_means(lon_sim,lat_sim,x)),order.by=time_sim),-22000,-6000)
  
  # LOVECLIM
  nc <- nc_open(file.path(DIR_MODEL_DATA,"LOVECLIM/deglaciation_ns/annual_total_precipitation.nc"))
  lon <- ncvar_get(nc,"LONN31_33")
  lat <- ncvar_get(nc,"LAT")
  time <- ncvar_get(nc,"TAX")-18000  #("YRS AFTER 18K BP")
  clim_var <- ncvar_get(nc,"PP")*10
  decs <- seq(1,dim(clim_var)[3],by=10)
  time <- -18000 + decs - 5
  lclim_pr <- sapply(decs, function(t) apply(clim_var[,,t:t+9], c(1,2), mean)) %>% array(dim = c(dim(clim_var)[c(1,2)],length(decs)))
  lclim_pr <- paleodata_windowing(zoo(apply(lclim_pr,3,function(x) spatial_means(lon,lat,x)),order.by=time),-22000,-6000)
  
  # BRIDGE HadCM3
  # Mean precip
  bridge_precip <- array(0,dim=c(length(bridge_time),length(bridge_lon),length(bridge_lat)))
  for (i in 1:length(bridge_time)) {
    nc <- nc_open(paste("/stacydata/data/BRIDGE/bbc_all_triff_rec_dyn04/",sim_names[i],".precip_mm_srf.monthly.nc",sep=""))
    bridge_precip[i,,] <- apply(ncvar_get(nc,"precip_mm_srf")[,length(bridge_lat):1,9601:12000],c(1,2),mean)
  }
  bridge_precip <- bridge_precip*60*60*24*30*12 # = Monthly sum in mm from kg/(s*m^2)
  bridge_pr <- paleodata_windowing(zoo(apply(bridge_precip,1,function(x) spatial_means(bridge_lon,bridge_lat,x)),order.by=-bridge_time*1000),-22000,-6000)
  
  ### compute orbital insolation forcing data with palinsol (Berger '78) ----
  # lat in degree
  # long: 1=June solstice, 3=December solstice
  insolation <- function(times, astrosol=palinsol::ber78,long=1,lat=65) {
    lat <- lat*pi/180
    long <- long*pi/2 
    return(sapply(times, function(tt) palinsol::Insol(orbit=astrosol(tt),long=long,lat=lat)))
  }
  
  insol_times <- seq(from=-22000,to=-6000,by=200)
  insol45N_J_forc <- insolation(times=insol_times, lat = 45)
  insol45S_D_forc <- insolation(times=insol_times, long=3, lat = -45)
  
  if (all(is.na(insol45N_J_forc))) {
    load(file.path(DIR_CACHE,'pre_processed_insol.RData'))
  }
} # end !FIG1DATA_CACHED

## cache data
if (!FIG1DATA_CACHED) {
  save(
    insol_times,insol45N_J_forc,insol45S_D_forc,
    lclim_pr,lclim_ts,
    trace_pr,trace_ts,
    bridge_pr,bridge_ts,bridge_time,
    CO2_forc,trace_lsm_lgm,
    file = FILE_FIG1DATA_CACHE
  )
}

## labels and styling for time series plot ----
all_data_tseries_labs <- c("bold(atop(NA,atop('\u0394'*S[0],'[W/'*m^{2}*']')))", "bold(atop(NA,atop('\u0394'*S[0],'[W/'*m^{2}*']')))",
                           "bold(atop(NA,atop('CO2','[ppm]')))", 
                           "bold(atop(NA,atop('\u0394GMST','[K]')))", "bold(atop(NA,atop('\u0394GMPR','[mm/a]')))", 
                           "bold(atop(NA,atop('\u0394GMST','[K]')))", "bold(atop(NA,atop('\u0394GMPR','[mm/a]')))", 
                           "bold(atop(NA,atop('\u0394GMST','[K]')))", "bold(atop(NA,atop('\u0394GMPR','[mm/a]')))")

all_data_tseries_grps <- c('45\u00B0N, Jun', '45\u00B0S, Dec', 'CO2', 'TraCE', 'TraCE', 'LOVECLIM', 'LOVECLIM', 'HadCM3', 'HadCM3')

all_data_tseries <- tibble(signal_name = all_data_tseries_labs, 
                           group = all_data_tseries_grps,
                           time = list(insol_times, insol_times, index(CO2_forc), index(trace_ts), index(trace_pr),
                                       index(lclim_ts), index(lclim_pr), index(bridge_ts), index(bridge_pr)),
                           signal = list(insol45N_J_forc-coredata(insol45N_J_forc)[length(coredata(insol45N_J_forc))], insol45S_D_forc-coredata(insol45S_D_forc)[length(coredata(insol45S_D_forc))], 
                                         coredata(CO2_forc), 
                                         coredata(trace_ts)-coredata(trace_ts)[length(coredata(trace_ts))], coredata(trace_pr)-coredata(trace_pr)[length(coredata(trace_pr))],
                                         coredata(lclim_ts)-coredata(lclim_ts)[length(coredata(lclim_ts))], coredata(lclim_pr)-coredata(lclim_pr)[length(coredata(lclim_pr))], 
                                         coredata(bridge_ts)-coredata(bridge_ts)[length(coredata(bridge_ts))], 
                                         coredata(bridge_pr)-coredata(bridge_pr)[length(coredata(bridge_pr))]))

all_data_tseries_u <- all_data_tseries %>% 
  unnest(cols = c(time, signal))
all_data_tseries_u$signal_name <- factor(all_data_tseries_u$signal_name, 
                                         levels = c(all_data_tseries_labs[1], all_data_tseries_labs[3], all_data_tseries_labs[4], all_data_tseries_labs[5]))

all_data_tseries_u$group <- factor(all_data_tseries_u$group,
                                   levels = c(all_data_tseries_grps[1], all_data_tseries_grps[2], all_data_tseries_grps[3], 
                                              all_data_tseries_grps[6], all_data_tseries_grps[4], all_data_tseries_grps[8]))

clrs <- c('#6e77a7', '#d7191c', '#e69f00', '#cc79a7', '#009e73', '#0072b3')
alphas <- c(1,1,1,1,1,1)


## Panel A plot: global mean time series & forcings plot ----
tseries_models_plt <- ggplot(all_data_tseries_u,
                             aes(color = group, x = -1*time/1000, y = signal, alpha = group)) + 
  geom_line() + 
  facet_wrap(. ~signal_name, scales = 'free_y', strip.position = 'left', ncol = 1, labeller = label_parsed, drop = F) + 
  scale_y_continuous(position = 'right', breaks = scales::pretty_breaks(n = 3)) +#sec.axis = sec_axis(~ .)) + 
  scale_x_continuous(labels = c(6,14,22), breaks = c(6,14,22), expand = c(0,0)) + 
  scale_alpha_manual(values = alphas, guide = FALSE) + 
  scale_color_manual(values = clrs, guide = guide_legend('Forcing/Simulation', nrow=2, override.aes = list(size = 4), title.position = 'top')) + 
  global_title_and_axis() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.y = element_blank(), axis.text = element_text(size = GLOBAL_FONT_SIZE + 2),
        strip.placement = 'outside', strip.background.y = element_blank(), strip.text = element_text(size = GLOBAL_FONT_SIZE + 8),
        panel.spacing.y = unit(0.02, 'npc'), legend.position = 'bottom', legend.direction = 'horizontal') + #, panel.background = element_rect(color = 'black', inherit.blank = F, size = 1)) + 
  labs(x = 'time [ka BP]')

tseries_models_grb <- ggplotGrob(tseries_models_plt)
tseries_models_grb$grobs[[2]]$children[[5]] <- segmentsGrob(c(1,0),c(0,0),c(1,0),c(1,1))
tseries_models_grb$grobs[[3]]$children[[5]] <- segmentsGrob(c(1,0),c(0,0),c(1,0),c(1,1))
tseries_models_grb$grobs[[4]]$children[[5]] <- segmentsGrob(c(1,0),c(0,0),c(1,0),c(1,1))
tseries_models_grb$grobs[[5]]$children[[5]] <- segmentsGrob(c(1,0,0),c(0,0,0),c(1,0,1),c(1,1,0))

## for saving model & forcing time series plot separately
#cmbplt_mseries <- gridExtra::arrangeGrob(tseries_models_grb)
#ggsave(plot = cmbplt_mseries, file.path(DIR_FIGURES, 'map_ts_modelseries_only_update21.pdf'), 
#       width = 38, height = 22, units = 'cm', dpi = 'print', limitsize = FALSE, device = cairo_pdf)


## combine time series with map ----
cmbplt_mseries <- plot_grid(tseries_models_grb, cmbplt, nrow = 1,
                            rel_widths = c(1,2.6), rel_heights = c(0.8,1),
                            labels = "AUTO", label_size = GLOBAL_FONT_SIZE + 5)

## save combined plot ----
ggsave(plot = cmbplt_mseries, file.path(DIR_FIGURES, 'map_ts_modelseries_update21_lr.pdf'), 
       width = 52, height = 22, units = 'cm', dpi = 'print', limitsize = FALSE, device = cairo_pdf)
ggsave(plot = cmbplt_mseries, file.path(DIR_FIGURES, 'map_ts_modelseries_update21_lr.png'), 
       width = 52, height = 22, units = 'cm', dpi = 'print', limitsize = FALSE)

