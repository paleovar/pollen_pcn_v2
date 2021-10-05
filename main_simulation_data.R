# Wrapper script for initalizing the climate simulation data
# as proxy-like time series in the same data structure as the ACER pollen data
# The wrapper script is `source()`ed at several occasions in the create_fig* scripts and in `main_pcns.R`

# NOTE:
# Per default model data that has been pre-processed is used in this repository if available.
# This is to allow re-creating the plots without the need to obtain external data. Still,
# the code is fully functional if this flag is set to `FALSE` and the data has been obtained
# as described in the tutorial `get_started.Rmd`.
# In any of these two cases the processed model data will not be re-loaded or computed if the 
# required objects are already present in the global environment. If you want to force re-loading
# or re-computed (depending on the flag state of `USE_PREPROCESSED_DATA`), set USE_MODEL_DATA_GLOBAL_ENV to `FALSE`
# CAUTION: Forcing re-intialisation by USE_MODEL_DATA_GLOBAL_ENV <- FALSE causes `main_pcns.R` to re-compute
# correlations
USE_PREPROCESSED_DATA <- TRUE
USE_MODEL_DATA_GLOBAL_ENV <- TRUE

FILE_MODEL_DATA_CACHE <- file.path(DIR_CACHE,'pre_processed_model_data.RData')


# script flag preventing script to be executed if this has been done in the session before
# overwrite STATUS_MODEL_SIMULATION_DATA_RUN <- FALSE in global environment to force re-execution
# Note that this does not influence if data is read from cache or re-produced from the raw data
if (!exists('STATUS_MODEL_SIMULATION_DATA_RUN')) STATUS_MODEL_SIMULATION_DATA_RUN <- FALSE

if (STATUS_MODEL_SIMULATION_DATA_RUN == FALSE) {
  message('main_simulation_data.R has not been run in session, executing script with the defined flags')
} else {
  message('main_simulation_data.R has been run in session, aborting executing the script.\nTo force execution and reload the data from file/cache, set\nSTATUS_MODEL_SIMULATION_DATA_RUN <- FALSE in global environment')
  invokeRestart('abort')
}

# if pre-processed simulation data should be used
if (USE_PREPROCESSED_DATA) {
  if (USE_MODEL_DATA_GLOBAL_ENV) {
    status_model_data_global_env <- all(sapply(c('ACERtrc', 'ACERtrcmil', 'ACERlc', 'ACERlcmil', 'ACERhcm'), exists))
  }
  if (status_model_data_global_env == FALSE & file.exists(FILE_MODEL_DATA_CACHE)) {
    message('Using pre-processed pseudo proxies from simulation data, therefore skipping aggregating them from scratch')
    load(FILE_MODEL_DATA_CACHE)
    STATUS_MODEL_SIMULATION_DATA_RUN <- TRUE
    ## this loads main data objects into global environment:
    ## `ACERtrc`, `ACERtrcmil` - TraCE21ka simulation
    ## `ACERlc`, `ACERlcmil` - LOVECLIM Last Deglaciation simulation
    ## `ACERhcm` - HadCM3 BRIDGE simulation of the Last Glacial
    invokeRestart('abort')
  } else if (status_model_data_global_env) {
    message('Pseudo proxies are in global environment already. ')
    invokeRestart('abort')
  } else {
    message('USE_PREPROCESSED_DATA == TRUE but no matching .RData file found, attempting to aggregate pseudo proxies from scratch')
  }
}

# otherwise, aggregate the pseudo proxies from scratch
# TRACE ----
ACERhere <- PCNdata() %>% 
  mix_sample_dating_and_err() %>% 
  harmonized_pollen_to_db() %>% 
  arboreal_pollen_to_db()

nc <- nc_open(file.path(DIR_MODEL_DATA,'trace21ka/processed/trace.01-36.22000BP.PFTFRAC_decavg_400BCE.nc'))
trace_pfts <- ncvar_get(nc, 'PFTFRAC')
trace_time <- ncvar_get(nc, 'time')
trace_lats <- ncvar_get(nc, 'lat')
trace_lons <- ncvar_get(nc, 'lon')
nc_close(nc)
trace_time <- trace_time * -1e3
#trace_ap <- apply(trace_pfts[,,,c(2:9)], 1:3, sum, na.rm = T)/apply(trace_pfts[,,,c(2:15)], 1:3, sum, na.rm = T)
#trace_sb <- apply(trace_pfts[,,,c(10:12)], 1:3, sum, na.rm = T)/apply(trace_pfts[,,,c(2:15)], 1:3, sum, na.rm = T)
trace_apsb <- apply(trace_pfts[,,,c(2:12)], 1:3, sum, na.rm = T)/apply(trace_pfts[,,,c(2:15)], 1:3, sum, na.rm = T)
rm(trace_pfts)
trace_apsb <- aperm(trace_apsb, c(2,3,1)) # lon,lat,time

nc <- nc_open(file.path(DIR_MODEL_DATA,'trace21ka/trace.01-36.22000BP.cam2.PRECT.22000BP_decavg_400BCE.nc'))
trace_pr <- ncvar_get(nc, 'PRECT') # has lon,lat,time
nc_close(nc)

nc <- nc_open(file.path(DIR_MODEL_DATA,'trace21ka/trace.01-36.22000BP.cam2.TREFHT.22000BP_decavg_400BCE.nc'))
trace_ts <- ncvar_get(nc, 'TREFHT') # has lon,lat,time
nc_close(nc)

# use only ACER grid cells and create time series the ACER way
p4sCRS <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
ACERlocs <- ACERhere$sites %>% 
  select(long, lat) %>% 
  SpatialPoints()
proj4string(ACERlocs) <- p4sCRS

TRACEgrid <- expand.grid(trace_lons, trace_lats) %>% 
  as_tibble() %>% 
  rename(long = Var1, lat = Var2) %>% 
  mutate(long = if_else(long > 180, long-360, long)) %>% 
  SpatialPoints()
proj4string(TRACEgrid) <- p4sCRS

trace_lons[trace_lons > 180] <- trace_lons[trace_lons > 180]-360

TRACE_apsb_r <- retrieve_nearest_tseries(field = trace_apsb, match_locs = ACERlocs, 
                                         match_sites = ACERhere$sites$site_id,
                                         grid = TRACEgrid, lons = trace_lons, lats = trace_lats, time = trace_time,
                                         search_depth = 4) # tested up to `search_depth = 4` needed for TraCE-21ka simulation
# use `View(TRACE_apsb_r$iterations)` to check the progress of retrieving
TRACE_apsb <- TRACE_apsb_r$tseries
TRACEijs <- TRACE_apsb_r$center_ixs_used
TRACE_ts <- retrieve_tseries(field = trace_ts, center_ixs = TRACEijs,
                             match_sites = ACERhere$sites$site_id,
                             time = trace_time)
TRACE_pr <- retrieve_tseries(field = trace_pr, center_ixs = TRACEijs,
                             match_sites = ACERhere$sites$site_id,
                             time = trace_time)

rm(ACERlocs,TRACEgrid,TRACE_apsb_r)

# irregular sampling
# aggregate to ACER time points, read null model sets
ACERtrc <- lapply(list(TRACE_apsb, TRACE_pr, TRACE_ts), function(tr) {
  tr <- tr %>% 
    rownames_to_column('depth') %>% # just as pseudo if ordering by depth somewhere
    gather(key = 'site_id', value = 'pcnt_arb_pollen', -mixed_age, -depth) %>% 
    mutate_at(vars(site_id), as.numeric)
  
  tr <- tr %>% 
    group_by(site_id) %>% 
    nest(.key = 'model') %>% 
    inner_join(ACERhere$sample_dating %>% 
                 inner_join(ACERhere$arb_pollen_data %>% select(sample_id)) %>% # make sure samples match (filters != 'COUN' and full zero samples, see `pollen_percentize`)
                 group_by(site_id) %>%
                 nest(.key = 'ACER')) %>% ungroup() %>% 
    mutate(model_agg = purrr::map2(.$model, .$ACER, function(M,A) {
      sample.dataset.box(M %>% rename(age = mixed_age) %>% select(age, pcnt_arb_pollen),
                         A %>% arrange(mixed_age) %>% .$mixed_age,
                         'pcnt_arb_pollen')
    })) %>% #print() %>% 
    unnest(model_agg) %>% 
    select(-bin) %>% 
    rename(mixed_age = age) %>% 
    rownames_to_column('arb_pollen_data_id') %>% 
    mutate_at(vars(arb_pollen_data_id), as.numeric) %>% 
    inner_join(ACERhere$sample_dating %>% 
                 group_by(site_id) %>% 
                 select(site_id, sample_id, mixed_age),
               by = c('site_id', 'mixed_age'))
  
  ACERloc <- PCNdata() %>% 
    mix_sample_dating_and_err()
  ACERloc$arb_pollen_data <- select(tr, site_id, sample_id, arb_pollen_data_id, pcnt_arb_pollen)
  ACERloc$pollen_data <- NULL
  return(ACERloc)
}) %>% 
  setNames(c('TRACE_apsb', 'TRACE_pr', 'TRACE_ts'))
rm(TRACE_apsb, TRACE_pr, TRACE_ts)

# save aggregated data to disk (without site filtering)
lapply(1:length(ACERtrc), function(i) {
  nms <- names(ACERtrc)
  dat <- ACERtrc[[i]]
  dir <- file.path(DIR_DATASETS_SUPPLEMENTARY, "model_data/6_22_TRACE/dat_updt3") # was /dat, /dat_updt only previously
  if (!dir.exists(dir)) dir.create(dir, recursive = T)
  file <- file.path(dir, paste0('6_22_TRACE_', nms[i], '.RData'))
  save(dat, file = file)
})

# save windowed & filtered data to disk (i.e. only sites that are used in the following)
ACERtrc_temp <- ACERtrc
ACERtrc_temp[['ACER_ap']] <- ACERhere

lapply(1:length(ACERtrc_temp), function(i) {
  nms <- names(ACERtrc_temp)
  dat <- filter_sites(ACERtrc_temp[[i]]$sites, ACERtrc_temp[[i]]$sample_dating, hres_only = T, 'all', 'all') %>% 
    inner_join(ACERtrc_temp[[i]]$arb_pollen_data) %>% 
    inner_join(ACERtrc_temp[[i]]$sample_dating) %>% 
    select(arb_pollen_data_id,site_id,sample_id,mixed_age,pcnt_arb_pollen) %>% 
    qtransform_data(data = ., type_data = 'arboreal_pollen', transform = 'identity', sd_one = F) %>% 
    window_and_detrend(data = ., windows = WINDOWS_TRACE, labels = WINDOW_LABELS_TRACE, type_data = 'arboreal_pollen', transform = 'identity', detrend = list(activate = FALSE), sd_one = F) %>% 
    unnest(cols = c(data_win))
  dir <- file.path(DIR_DATASETS_SUPPLEMENTARY, "model_data/6_22_TRACE/dat_updt3") # was /dat, /dat_updt only previously
  if (!dir.exists(dir)) dir.create(dir, recursive = T)
  file <- file.path(dir, paste0('6_22_TRACE_', nms[i], '_proc.RData'))
  save(dat, file = file)
})

# save unharmonized ACER data to disk
dat <- ACERhere$sites %>% select(site_id) %>% 
  inner_join(ACERhere$pollen_data) %>% 
  inner_join(ACERhere$sample_dating) %>% 
  select(pollen_data_id,site_id,sample_id,mixed_age,taxon_pcnt,taxon) %>% 
  window_data(data = ., type_data = 'pollen', transform = NULL, windows = WINDOWS_TRACE, labels = WINDOW_LABELS_TRACE)
dir <- file.path(DIR_DATASETS_SUPPLEMENTARY, "model_data/6_22_TRACE/dat_updt3") # was /dat, /dat_updt only previously
if (!dir.exists(dir)) dir.create(dir, recursive = T)
file <- file.path(dir, paste0('6_22_TRACE_ACER_orig.RData'))
save(dat, file = file)

rm(dat,dir,file,ACERtrc_temp,ACERhere)


# repeat initialisation for ACER AP + ice core d18O and TraCE pseudo proxies
# for millennial scales 
# Note that links to EPICA and NGRIP are filtered for all network measures in
# `main_pcns.R`. Including EPICA and NGRIP here saves some additional data handling.
# For LOVECLIM on millennial scales ACERlcmil is simply copied from the same initialisation instead. ----
ACERhere <- PCNdata() %>% 
  mix_sample_dating_and_err() %>% 
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
ACERhere$arb_pollen_data <- ACERhere$ACERarbpoll_NGECice_data %>% 
  select(-proxy_type) %>% 
  rename(arb_pollen_data_id = proxy_id, pcnt_arb_pollen = proxy)

# use only ACER grid cells and create time series the ACER way
p4sCRS <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
ACERlocs <- ACERhere$sites %>% 
  select(long, lat) %>% 
  SpatialPoints()
proj4string(ACERlocs) <- p4sCRS

TRACEgrid <- expand.grid(trace_lons, trace_lats) %>% 
  as_tibble() %>% 
  rename(long = Var1, lat = Var2) %>% 
  mutate(long = if_else(long > 180, long-360, long)) %>% 
  SpatialPoints()
proj4string(TRACEgrid) <- p4sCRS

trace_lons[trace_lons > 180] <- trace_lons[trace_lons > 180]-360

TRACE_apsb_r <- retrieve_nearest_tseries(field = trace_apsb, match_locs = ACERlocs, 
                                         match_sites = ACERhere$sites$site_id,
                                         grid = TRACEgrid, lons = trace_lons, lats = trace_lats, time = trace_time,
                                         search_depth = 4, # tested up to `search_depth = 4` needed for TraCE-21ka simulation (apart from ice core locations obviously)
                                         overwrite_sites = which(ACERhere$sites$site_id %in% c(101,102))) # EPICA, NGRIP return NA for apsb, but should remain at their site for TS, PR 
# use `View(TRACE_apsb_r$iterations)` to check the progress of retrieving
TRACE_apsb <- TRACE_apsb_r$tseries
TRACEijs <- TRACE_apsb_r$center_ixs_used
TRACE_ts <- retrieve_tseries(field = trace_ts, center_ixs = TRACEijs,
                             match_sites = ACERhere$sites$site_id,
                             time = trace_time)
TRACE_pr <- retrieve_tseries(field = trace_pr, center_ixs = TRACEijs,
                             match_sites = ACERhere$sites$site_id,
                             time = trace_time)

rm(trace_time, trace_lats, trace_lons, trace_apsb, trace_pr, trace_ts,
   TRACE_apsb_r, ACERlocs, TRACEgrid, TRACEijs)

# aggregate to ACER time points, read null model sets
ACERtrcmil <- lapply(list(TRACE_apsb, TRACE_pr, TRACE_ts), function(tr) {
  tr <- tr %>% 
    rownames_to_column('depth') %>% # just as pseudo if ordering by depth somewhere
    gather(key = 'site_id', value = 'pcnt_arb_pollen', -mixed_age, -depth) %>% 
    mutate_at(vars(site_id), as.numeric)
  
  tr <- tr %>% 
    group_by(site_id) %>% 
    nest(.key = 'model') %>% 
    inner_join(ACERhere$sample_dating %>% 
                 inner_join(ACERhere$arb_pollen_data %>% select(sample_id)) %>% # make sure samples match (filters != 'COUN' and full zero samples, see `pollen_percentize`)
                 group_by(site_id) %>%
                 nest(.key = 'ACER')) %>%
    ungroup() %>% 
    mutate(model_agg = purrr::map2(.$model, .$ACER, function(M,A) {
      sample.dataset.box(M %>% rename(age = mixed_age) %>% select(age, pcnt_arb_pollen),
                         A %>% arrange(mixed_age) %>% .$mixed_age,
                         'pcnt_arb_pollen')
    })) %>% 
    unnest(model_agg) %>% 
    select(-bin) %>% 
    rename(mixed_age = age) %>% 
    rownames_to_column('arb_pollen_data_id') %>% 
    mutate_at(vars(arb_pollen_data_id), as.numeric) %>% 
    inner_join(ACERhere$sample_dating %>% 
                 group_by(site_id) %>% 
                 select(site_id, sample_id, mixed_age),
               by = c('site_id', 'mixed_age'))
  
  ACERloc <- PCNdata() %>% 
    mix_sample_dating_and_err() %>% 
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
  ACERloc$arb_pollen_data <- select(tr, site_id, sample_id, arb_pollen_data_id, pcnt_arb_pollen)
  ACERloc$pollen_data <- NULL
  ACERloc$ACERarbpoll_NGECice_data <- NULL
  ACERloc$ACERdefaultsites <- NULL
  return(ACERloc)
}) %>% 
  setNames(c('TRACE_apsb', 'TRACE_pr', 'TRACE_ts'))
rm(TRACE_apsb,TRACE_pr,TRACE_ts,ACERhere)

# LOVECLIM ----
ACERhere <- PCNdata() %>% 
  mix_sample_dating_and_err() %>% 
  harmonized_pollen_to_db() %>% 
  arboreal_pollen_to_db()

nc <- nc_open(file.path(DIR_MODEL_DATA,'/LOVECLIM/deglaciation_ns/annual_surface_temperature_at_2_meter.nc'))
lc_ts <- ncvar_get(nc, 'T2M')
lc_time <- (ncvar_get(nc, 'TAX')-18000)*-1
lc_lats <- ncvar_get(nc, 'LAT')
lc_lons <- ncvar_get(nc, 'LONN31_33')
nc_close(nc)
nc <- nc_open(file.path(DIR_MODEL_DATA,'/LOVECLIM/deglaciation_ns/annual_total_precipitation.nc'))
lc_pr <- ncvar_get(nc, 'PP')
nc_close(nc)

# use only ACER grid cells and create time series the ACER way
p4sCRS <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
ACERlocs <- ACERhere$sites %>% 
  select(long, lat) %>% 
  SpatialPoints()
proj4string(ACERlocs) <- p4sCRS

LCgrid <- expand.grid(lc_lons, lc_lats) %>% 
  as_tibble() %>% 
  rename(long = Var1, lat = Var2) %>% 
  mutate(long = if_else(long > 180, long-360, long)) %>% 
  SpatialPoints()
proj4string(LCgrid) <- p4sCRS

lc_lons[lc_lons > 180] <- lc_lons[lc_lons > 180]-360

LC_ts_r <- retrieve_nearest_tseries(field = lc_ts, match_locs = ACERlocs, 
                                    match_sites = ACERhere$sites$site_id,
                                    grid = LCgrid, lons = lc_lons, lats = lc_lats, time = lc_time,
                                    search_depth = 1) # LOVECLIM ts is a full field, thus `search_depth = 1` is sufficient
# use `View(LC_ts_r$iterations)` to check the progress of retrieving
LC_ts <- LC_ts_r$tseries
LCijs <- LC_ts_r$center_ixs_used
LC_pr <- retrieve_tseries(field = lc_pr, center_ixs = LCijs,
                          match_sites = ACERhere$sites$site_id,
                          time = lc_time)

rm(lc_time, lc_lats, lc_lons, lc_pr, lc_ts,
   ACERlocs, LCgrid, LC_ts_r, LCijs)

# irregular sampling
ACERlc <- lapply(list(LC_ts, LC_pr), function(lc) {
  lc <- lc %>% 
    rownames_to_column('depth') %>% # just as pseudo if ordering by depth somewhere
    gather(key = 'site_id', value = 'pcnt_arb_pollen', -mixed_age, -depth) %>% 
    mutate_at(vars(site_id), as.numeric)
  
  lc <- lc %>% 
    group_by(site_id) %>% 
    nest(.key = 'model') %>% 
    inner_join(ACERhere$sample_dating %>% 
                 inner_join(ACERhere$arb_pollen_data %>% select(sample_id)) %>% # make sure samples match (filters != 'COUN' and full zero samples, see `pollen_percentize`)
                 group_by(site_id) %>%
                 #filter(mixed_age > 6000 & mixed_age < 18000) %>% 
                 nest(.key = 'ACER')) %>% ungroup() %>% 
    mutate(model_agg = purrr::map2(.$model, .$ACER, function(M,A) {
      sample.dataset.box(M %>% rename(age = mixed_age) %>% select(age, pcnt_arb_pollen),
                         A %>% arrange(mixed_age) %>% .$mixed_age,
                         'pcnt_arb_pollen')
    })) %>% 
    unnest(model_agg) %>% 
    select(-bin) %>% 
    rename(mixed_age = age) %>% 
    rownames_to_column('arb_pollen_data_id') %>% 
    mutate_at(vars(arb_pollen_data_id), as.numeric) %>% 
    inner_join(ACERhere$sample_dating %>% 
                 group_by(site_id) %>% 
                 select(site_id, sample_id, mixed_age),
               by = c('site_id', 'mixed_age'))
  
  ACERloc <- PCNdata() %>% 
    mix_sample_dating_and_err()
  ACERloc$arb_pollen_data <- select(lc, site_id, sample_id, arb_pollen_data_id, pcnt_arb_pollen)
  ACERloc$pollen_data <- NULL
  return(ACERloc)
}) %>% 
  setNames(c('LC_ts', 'LC_pr'))
rm(LC_ts,LC_pr)

# save aggregated data to disk
lapply(1:length(ACERlc), function(i) {
  nms <- names(ACERlc)
  dat <- ACERlc[[i]]
  dir <- file.path(DIR_DATASETS_SUPPLEMENTARY, "model_data/6.2_18_LC/dat_updt3") # was /dat, /dat_updt2 only previously
  if (!dir.exists(dir)) dir.create(dir, recursive = T)
  file <- file.path(dir, paste0('6.2_18_LC', nms[i], '.RData'))
  save(dat, file = file)
})

# save windowed data to disk
ACERlc_temp <- ACERlc
ACERlc_temp[['ACER_ap']] <- ACERhere

lapply(1:length(ACERlc), function(i) {
  nms <- names(ACERlc)
  dat <- filter_sites(ACERlc[[i]]$sites, ACERlc[[i]]$sample_dating, hres_only = T, 'all', 'all') %>% 
    inner_join(ACERlc[[i]]$arb_pollen_data) %>% 
    inner_join(ACERlc[[i]]$sample_dating) %>% 
    select(arb_pollen_data_id,site_id,sample_id,mixed_age,pcnt_arb_pollen) %>% 
    qtransform_data(data = ., type_data = 'arboreal_pollen', transform = 'identity', sd_one = F) %>% 
    window_and_detrend(data = ., windows = WINDOWS_LC, labels = WINDOW_LABELS_LC, type_data = 'arboreal_pollen', transform = 'identity', detrend = list(activate = FALSE), sd_one = F) %>% 
    unnest()
  dir <- file.path(DIR_DATASETS_SUPPLEMENTARY, "model_data/6.2_18_LC/dat_updt3") # was /dat, /dat_updt2 only previously
  if (!dir.exists(dir)) dir.create(dir, recursive = T)
  file <- file.path(dir, paste0('6.2_18_LC_', nms[i], '_proc.RData'))
  save(dat, file = file)
})

# save unharmonized ACER data to disk
dat <- ACERhere$sites %>% select(site_id) %>% 
  inner_join(ACERhere$pollen_data) %>% 
  inner_join(ACERhere$sample_dating) %>% 
  select(pollen_data_id,site_id,sample_id,mixed_age,taxon_pcnt,taxon) %>% 
  window_data(data = ., type_data = 'pollen', transform = NULL, windows = WINDOWS_LC, labels = WINDOW_LABELS_LC)
dir <- file.path(DIR_DATASETS_SUPPLEMENTARY, "model_data/6.2_18_LC/dat_updt3")
if (!dir.exists(dir)) dir.create(dir, recursive = T)
file <- file.path(dir, paste0('6._18_LC_ACER_orig.RData'))
save(dat, file = file)

rm(dat,dir,file,ACERlc_temp,ACERhere)

# copy for millennial scale processing
ACERlcmil <- ACERlc


# HadCM3 BRIDGE simulation ----
ACERhere <- PCNdata() %>% 
  mix_sample_dating_and_err() %>% 
  harmonized_pollen_to_db() %>% 
  arboreal_pollen_to_db()

BBC_sim_names <- paste(read.csv("/stacydata/data/BRIDGE/BRIDGE.csv",header=TRUE)$bbc_all_triff_rev_dyn04)[1:62]
BBC_sim_dates <- c(seq(120,80,by=-4),seq(78,22,by=-2),seq(21,0,by=-1))

nc <- nc_open('/stacywork/fredweasley/data/BRIDGE_work/21ka.PFTsumnorm.annual.nc')#'/stacywork/fredweasley/data/BRIDGE_work/21ka.PFTsum.annual.nc')
bbc_lons <- ncvar_get(nc, 'longitude')
bbc_lats <- ncvar_get(nc, 'latitude')
bbc_time <- 1:22000

bbc_ap <- ncvar_get(nc,'fracPFTs_mm_srf')
nc_close(nc)

nc <- nc_open('/stacywork/fredweasley/data/BRIDGE_work/21ka.precip_mm_srf.annual.nc')
bbc_pr <- ncvar_get(nc,'precip_mm_srf')
nc_close(nc)

nc <- nc_open('/stacywork/fredweasley/data/BRIDGE_work/21ka.temp_mm_1_5m.annual.nc') #temp_mm_srf.annual.nc') for surface temp
bbc_ts <- ncvar_get(nc,'temp_mm_1_5m')#'temp_mm_srf')
nc_close(nc)

# use only ACER grid cells and create time series the ACER way
p4sCRS <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
ACERlocs <- ACERhere$sites %>% 
  select(long, lat) %>% 
  SpatialPoints()
proj4string(ACERlocs) <- p4sCRS

BBCgrid <- expand.grid(bbc_lons, bbc_lats) %>% 
  as_tibble() %>% 
  rename(long = Var1, lat = Var2) %>% 
  mutate(long = if_else(long > 180, long-360, long)) %>% 
  SpatialPoints()
proj4string(BBCgrid) <- p4sCRS

bbc_lons[bbc_lons > 180] <- bbc_lons[bbc_lons > 180]-360
BBC_ap_r <- retrieve_nearest_tseries(field = bbc_ap, match_locs = ACERlocs, 
                                     match_sites = ACERhere$sites$site_id,
                                     grid = BBCgrid, lons = bbc_lons, lats = bbc_lats, time = bbc_time,
                                     search_depth = 9) # tested `search_depth = 9` needed for HadCM3 BRIDGE simulation (all but site 75 need `search_depth = 5` at most)
# use `View(BBC_ap_r$iterations)` to check the progress of retrieving
BBC_ap <- BBC_ap_r$tseries
BBCijs <- BBC_ap_r$center_ixs_used
BBC_ts <- retrieve_tseries(field = bbc_ts, center_ixs = BBCijs,
                           match_sites = ACERhere$sites$site_id,
                           time = bbc_time)
BBC_pr <- retrieve_tseries(field = bbc_pr, center_ixs = BBCijs,
                           match_sites = ACERhere$sites$site_id,
                           time = bbc_time)
rm(BBC_ap_r,BBCgrid,BBCijs)

# use only ACER grid cells and create time series the ACER way

# irregular sampling
ACERhcm <- lapply(list(BBC_ap, BBC_pr, BBC_ts), function(hcm) {
  hcm <- hcm %>% 
    mutate(mixed_age = mixed_age + 500) %>%  # for the windowing
    rownames_to_column('depth') %>% # just as pseudo if ordering by depth somewhere
    gather(key = 'site_id', value = 'pcnt_arb_pollen', -mixed_age, -depth) %>% 
    mutate_at(vars(site_id), as.numeric)
  
  hcm <- hcm %>% 
    group_by(site_id) %>% 
    nest(.key = 'model') %>% 
    inner_join(ACERhere$sample_dating %>% 
                 inner_join(ACERhere$arb_pollen_data %>% select(sample_id)) %>% # make sure samples match (filters != 'COUN' and full zero samples, see `pollen_percentize`)
                 group_by(site_id) %>%
                 #filter(mixed_age < 22000 & mixed_age > 0) %>% 
                 nest(.key = 'ACER')) %>%
    ungroup() %>% 
    mutate(model_agg = purrr::map2(.$model, .$ACER, function(M,A) {
      sample.dataset.box(M %>% rename(age = mixed_age) %>% select(age, pcnt_arb_pollen),
                         A %>% arrange(mixed_age) %>% .$mixed_age,
                         'pcnt_arb_pollen')
    })) %>% 
    unnest(model_agg) %>% 
    select(-bin) %>% 
    rename(mixed_age = age) %>% 
    rownames_to_column('arb_pollen_data_id') %>% 
    mutate_at(vars(arb_pollen_data_id), as.numeric) %>% 
    inner_join(ACERhere$sample_dating %>% 
                 group_by(site_id) %>% 
                 select(site_id, sample_id, mixed_age),
               by = c('site_id', 'mixed_age'))
  
  ACERloc <- PCNdata() %>% 
    mix_sample_dating_and_err()
  ACERloc$arb_pollen_data <- select(hcm, site_id, sample_id, arb_pollen_data_id, pcnt_arb_pollen)
  ACERloc$pollen_data <- NULL
  return(ACERloc)
}) %>% 
  setNames(c('BBC_ap', 'BBC_pr', 'BBC_ts'))
rm(BBC_ap, BBC_pr, BBC_ts,bbc_ap, bbc_pr, bbc_ts,
   bbc_time, bbc_lats, bbc_lons, ACERlocs)

# save aggregated data to disk
lapply(1:length(ACERhcm), function(i) {
  nms <- names(ACERhcm)
  dat <- ACERhcm[[i]]
  dir <- file.path(DIR_DATASETS_SUPPLEMENTARY, "model_data/6_22_BRIDGE/dat_updt3") # was /dat, /dat_updt only previously
  if (!dir.exists(dir)) dir.create(dir, recursive = T)
  file <- file.path(dir, paste0('6_22_BRIDGE_', nms[i], '.RData'))
  save(dat, file = file)
})

# save windowed data to disk
ACERhcm_temp <- ACERhcm
ACERhcm_temp[['ACER_ap']] <- ACERhere

lapply(1:length(ACERhcm), function(i) {
  nms <- names(ACERhcm)
  dat <- filter_sites(ACERhcm[[i]]$sites, ACERhcm[[i]]$sample_dating, hres_only = T, 'all', 'all') %>% 
    inner_join(ACERhcm[[i]]$arb_pollen_data) %>% 
    inner_join(ACERhcm[[i]]$sample_dating) %>% 
    select(arb_pollen_data_id,site_id,sample_id,mixed_age,pcnt_arb_pollen) %>% 
    qtransform_data(data = ., type_data = 'arboreal_pollen', transform = 'identity', sd_one = F) %>% 
    window_and_detrend(data = ., windows = WINDOWS_BBC, labels = WINDOW_LABELS_BBC, type_data = 'arboreal_pollen', transform = 'identity', detrend = list(activate = FALSE), sd_one = F) %>% 
    unnest()
  dir <- file.path(DIR_DATASETS_SUPPLEMENTARY, "model_data/6_22_BRIDGE/dat_updt3") # was /dat, /dat_updt2 only previously
  if (!dir.exists(dir)) dir.create(dir, recursive = T)
  file <- file.path(dir, paste0('6_22_BRIDGE_', nms[i], '_proc.RData'))
  save(dat, file = file)
})
# save unharmonized ACER data to disk
dat <- ACERhere$sites %>% select(site_id) %>% 
  inner_join(ACERhere$pollen_data) %>% 
  inner_join(ACERhere$sample_dating) %>% 
  select(pollen_data_id,site_id,sample_id,mixed_age,taxon_pcnt,taxon) %>% 
  window_data(data = ., type_data = 'pollen', transform = NULL, windows = WINDOWS_BBC, labels = WINDOW_LABELS_BBC)
dir <- file.path(DIR_DATASETS_SUPPLEMENTARY, "model_data/6_22_BRIDGE/dat_updt3")
if (!dir.exists(dir)) dir.create(dir, recursive = T)
file <- file.path(dir, paste0('6_22_BRIDGE_ACER_orig.RData'))
save(dat, file = file)

rm(dat,dir,file,ACERhcm_temp,ACERhere)

# garbage collect ----
gc()

# cache pseudo proxy data if permitted ----
if (USE_PREPROCESSED_DATA) {
  save(
    ACERtrc, ACERtrcmil, ACERlc, ACERlcmil, ACERhcm,
    file = FILE_MODEL_DATA_CACHE
  )
}

# modifying run status of this script
STATUS_MODEL_SIMULATION_DATA_RUN <- TRUE
