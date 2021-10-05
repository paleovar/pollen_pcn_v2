# Preliminary wrapper script to compute
# - load and prepare emulated data of TraCE and HadCM3
#   for computation of paleoclimate networks
# as preparation for `main_emulation_pcns.R`.
#
# For fitting of the generalized additive model (GAM) see `main_fit_emulators.R`.


# NOTE (as in `main_simulation_data.R`):
# Per default model data that has been pre-processed is used in this repository if available.
# This is to allow re-creating the plots without the need to obtain external data. Still,
# the code is fully functional if this flag is set to `FALSE` and the data has been obtained
# as described in the tutorial `get_started.Rmd`.
# In any of these two cases the processed model data will not be re-loaded or computed if the 
# required objects are already present in the global environment. If you want to force re-loading
# or re-computed (depending on the flag state of `USE_PREPROCESSED_EMUL_DATA`), set USE_EMUL_DATA_GLOBAL_ENV to `FALSE`
# CAUTION: Forcing re-intialisation by USE_EMUL_DATA_GLOBAL_ENV <- FALSE causes `main_pcns.R` to re-compute
# correlations
USE_PREPROCESSED_EMUL_DATA <- TRUE
USE_EMUL_DATA_GLOBAL_ENV <- TRUE

FILE_EMUL_DATA_CACHE <- file.path(DIR_CACHE,'pre_processed_emul_data.RData')


# script flag preventing script to be executed if this has been done in the session before
# overwrite STATUS_EMUL_DATA_RUN <- FALSE in global environment to force re-execution
if (!exists('STATUS_EMUL_DATA_RUN')) STATUS_EMUL_DATA_RUN <- FALSE

if (STATUS_EMUL_DATA_RUN == FALSE) {
  message('main_emulation_data.R has not been run in session, executing script with the defined flags')
} else {
  message('main_emulation_data.R has been run in session, aborting executing the script. To force execution, set STATUS_EMUL_DATA_RUN <- FALSE in global environment')
  invokeRestart('abort')
}


# if pre-processed simulation data should be used
if (USE_PREPROCESSED_EMUL_DATA) {
  if (USE_EMUL_DATA_GLOBAL_ENV) {
    status_emul_data_global_env <- all(sapply(c('ACERtrc_em', 'ACERhcm_em'), exists))
  }
  if (status_emul_data_global_env == FALSE & file.exists(FILE_EMUL_DATA_CACHE)) {
    message('Using pre-processed emulated surrogates, therefore skipping aggregating them from scratch')
    load(FILE_EMUL_DATA_CACHE)
    STATUS_EMUL_DATA_RUN <- TRUE
    ## this loads main data objects into global environment:
    ## `ACERtrc_em` [$forceall, $forcets, $forcepr, $forceco2] - Emulated LPJ-DGVM (TraCE21ka simulation)
    ## `ACERhcm_em` [$forceall, $forcets, $forcepr, $forceco2] - Emulated TRIFFID (HadCM3 BRIDGE simulation of the Last Glacial)
    invokeRestart('abort')
  } else if (status_emul_data_global_env) {
    message('Pseudo proxies are in global environment already. ')
  } else {
    message('USE_PREPROCESSED_EMUL_DATA == TRUE but no matching .RData file found, attempting to aggregate emulated surrogates from scratch')
  }
}


# Initialize pseudo-proxies from simulation data first
# to feed into fitted emulators
source('main_simulation_data.R')

# load the emulator fit & interpolated CO2 forcing
# - `trace_gam_all_te`, `hadcm3_gam_all_te`, `trace_co2_ts`, `hadcm3_co2_ts`
source('main_fit_emulators.R')

## wrapper for the prediction
em_exps <- list(
  'forceall' = rep(TRUE,3), # combined forcing
  'forcets' = c(TRUE,FALSE,FALSE), # single temperature forcing
  'forcepr' = c(FALSE,TRUE,FALSE), # single precip forcing
  'forceco2' = c(FALSE,FALSE,TRUE) # single CO2 forcing
)

gams_list <- list('TRACE' = trace_gam_all_te,
                  'HADCM3' = hadcm3_gam_all_te) 

em_predict <- function(ts_data,pr_data,co2_data,gam,cfg) {
  if (cfg[1] == FALSE) {
    ts_data <- rep(mean(coredata(ts_data),na.rm=TRUE),times=length(ts_data))
  }
  if (cfg[2] == TRUE) {
    pr_data <- coredata(pr_data)
  } else {
    pr_data <- rep(mean(coredata(pr_data),na.rm=TRUE),times=length(pr_data))
  }
  if (cfg[3] == TRUE) {
    co2_data <- coredata(co2_data)
  } else {
    co2_data <- rep(mean(coredata(co2_data),na.rm=TRUE),times=length(co2_data))
  }
  out_data <- c(predict(gam,
                        type="response",
                        newdata=data.frame(x1=ts_data,x2=pr_data,x3=co2_data)))
  return(out_data)
}

# Emulated TraCE ----
ACERhere <- PCNdata() %>% 
  mix_sample_dating_and_err() %>% 
  harmonized_pollen_to_db() %>% 
  arboreal_pollen_to_db()

## predict TSF for combined and single forcings
precip_fac_trc <- 3600*24*360*1000

ACERtrc_co2 <- ACERhere$sample_dating %>% 
  group_by(site_id) %>%
  nest(.key = 'ACER') %>%
  ungroup() %>% 
  mutate(co2_agg = purrr::map(.$ACER, function(A) {
    sample.dataset.box(tibble(co2 = coredata(trace_co2_ts),
                              age = -1*index(trace_co2_ts)),
                       A %>% arrange(mixed_age) %>% .$mixed_age,
                       'co2')
  })) %>% 
  select(-ACER) %>% 
  unnest(co2_agg) %>% 
  select(-bin) %>% 
  rename(mixed_age = age) %>% 
  rownames_to_column('arb_pollen_data_id') %>% 
  mutate_at(vars(arb_pollen_data_id), as.numeric) %>% 
  inner_join(ACERhere$sample_dating %>% 
               group_by(site_id) %>% 
               select(site_id, sample_id, mixed_age) %>% 
               arrange(mixed_age),
             by = c('site_id', 'mixed_age')) %>% 
  filter(mixed_age > WINDOWS_TRACE[[1]] & mixed_age < WINDOWS_TRACE[[2]])

# predict tree-and-shrub fraction
ACERtrcloc_dat <- plyr::join_all(list(
    ACERtrc$TRACE_ts$arb_pollen_data %>% 
      rename(ts = pcnt_arb_pollen) %>% 
      inner_join(ACERtrc$TRACE_ts$sample_dating %>% 
                   select(sample_id,mixed_age)) %>% 
      select(-arb_pollen_data_id),
    ACERtrc$TRACE_pr$arb_pollen_data %>% 
      rename(pr = pcnt_arb_pollen) %>% 
      mutate(pr = precip_fac_trc * pr) %>% 
      inner_join(ACERtrc$TRACE_pr$sample_dating %>% 
                   select(sample_id,mixed_age)) %>% 
      select(-arb_pollen_data_id),
    ACERtrc_co2 %>% 
      select(-arb_pollen_data_id)
    ),
  type = 'inner',
  by = c('site_id', 'sample_id', 'mixed_age')) %>% 
  nest_by(site_id, .key = 'sim')

ACERtrc_em <- lapply(1:length(names(em_exps)), function(i) {
  new_arb_pollen <- tibble(
    site_id = ACERtrcloc_dat$site_id,
    emul = lapply(ACERtrcloc_dat$sim, function(x) tibble(sample_id = x$sample_id,
                                                         pcnt_arb_pollen = em_predict(x$ts,x$pr,x$co2,gam=gams_list[['TRACE']],cfg=em_exps[[i]])))
  ) %>% 
    unnest(emul) %>% 
    inner_join(ACERtrc$TRACE_ts$arb_pollen_data %>% select(-pcnt_arb_pollen),
               by = c('site_id', 'sample_id'))
  # recycle existing object
  new_PCNdata <- ACERtrc$TRACE_ts
  new_PCNdata$corr_list_gaussian_kernel_hres_identity_arb_pollen_data <- NULL
  new_PCNdata$corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data <- NULL
  new_PCNdata$arb_pollen_data <- new_arb_pollen
  new_PCNdata$pollen_data_percentized <- TRUE
  return(new_PCNdata)
}) %>% 
  setNames(names(em_exps))

rm(ACERtrcloc_dat,ACERhere)


# Emulated TRIFFID (HadCM3) ----
ACERhere <- PCNdata() %>% 
  mix_sample_dating_and_err() %>% 
  harmonized_pollen_to_db() %>% 
  arboreal_pollen_to_db()

## predict TSF for combined and single forcings
precip_fac_hcm <- 3600*24*360

ACERhcm_co2 <- ACERhere$sample_dating %>% 
  group_by(site_id) %>%
  nest(.key = 'ACER') %>%
  ungroup() %>% 
  mutate(co2_agg = purrr::map(.$ACER, function(A) {
    sample.dataset.box(tibble(co2 = coredata(hadcm3_co2_ts),
                              age = -1*index(hadcm3_co2_ts)),
                       A %>% arrange(mixed_age) %>% .$mixed_age,
                       'co2')
  })) %>% 
  select(-ACER) %>% 
  unnest(co2_agg) %>% 
  select(-bin) %>% 
  rename(mixed_age = age) %>% 
  rownames_to_column('arb_pollen_data_id') %>% 
  mutate_at(vars(arb_pollen_data_id), as.numeric) %>% 
  inner_join(ACERhere$sample_dating %>% 
               group_by(site_id) %>% 
               select(site_id, sample_id, mixed_age) %>% 
               arrange(mixed_age),
             by = c('site_id', 'mixed_age')) %>% 
  filter(mixed_age > WINDOWS_BBC[[1]] & mixed_age < WINDOWS_BBC[[2]])

# predict tree-and-shrub fraction
ACERhcmloc_dat <- plyr::join_all(list(
  ACERhcm$BBC_ts$arb_pollen_data %>% 
    rename(ts = pcnt_arb_pollen) %>% 
    inner_join(ACERhcm$BBC_ts$sample_dating %>% 
                 select(sample_id,mixed_age)) %>% 
    select(-arb_pollen_data_id),
  ACERhcm$BBC_pr$arb_pollen_data %>% 
    rename(pr = pcnt_arb_pollen) %>% 
    mutate(pr = precip_fac_hcm * pr) %>% 
    inner_join(ACERhcm$BBC_pr$sample_dating %>% 
                 select(sample_id,mixed_age)) %>% 
    select(-arb_pollen_data_id),
  ACERhcm_co2 %>% 
    select(-arb_pollen_data_id)
  ),
  type = 'inner',
  by = c('site_id', 'sample_id', 'mixed_age')) %>% 
  nest_by(site_id, .key = 'sim')

ACERhcm_em <- lapply(1:length(names(em_exps)), function(i) {
  new_arb_pollen <- tibble(
    site_id = ACERhcmloc_dat$site_id,
    emul = lapply(ACERhcmloc_dat$sim, function(x) tibble(sample_id = x$sample_id,
                                                         pcnt_arb_pollen = em_predict(x$ts,x$pr,x$co2,gam=gams_list[['HADCM3']],cfg=em_exps[[i]])))
  ) %>% 
    unnest(emul) %>% 
    inner_join(ACERhcm$BBC_ts$arb_pollen_data %>% select(-pcnt_arb_pollen),
               by = c('site_id', 'sample_id'))
  # recycle existing object
  new_PCNdata <- ACERhcm$BBC_ts
  new_PCNdata$corr_list_gaussian_kernel_hres_identity_arb_pollen_data <- NULL
  new_PCNdata$corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data <- NULL
  new_PCNdata$arb_pollen_data <- new_arb_pollen
  new_PCNdata$pollen_data_percentized <- TRUE
  return(new_PCNdata)
}) %>% 
  setNames(names(em_exps))

rm(ACERhcmloc_dat,ACERhere)


# cache pseudo proxy data if permitted ----
if (USE_PREPROCESSED_EMUL_DATA) {
  save(
    ACERtrc_em,
    ACERhcm_em,
    file = FILE_EMUL_DATA_CACHE
  )
}

# modifying run status of this script
STATUS_EMUL_DATA_RUN <- TRUE

