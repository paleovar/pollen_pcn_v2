# Wrapper script for computing 
# - the similarity of ACER and simulation data
# Note that this script relies on the individual draws from the spatio-temporal statistical null
# model and, therefore, need a decent amount of computational resources and time to complete.
# For details see the `get_started.Rmd` tutorial and for a small test case see 
# `main_test_pcn.R`.

# NOTE: 
# Per default null model draws and all computed similarity values are cached to disk
# to speed up re-loading data. Set this flag to `FALSE` if this behavior is not wanted 
# (not recommended, see below).
DO_CACHEING_PCN <- TRUE

# CAUTION: Drawing the many null model draws and computing the reference
# similarity values for the null model necessitates relatively large
# computational resources. Therefore, running the script for the very first time
# when the respective data have not yet been cached can take several hours or
# days, depending on the machine used.
#
# In this repository we provide readily-computed and pre-cached data for PCNs of
# the ACER AP and of PCNS of simulation data tested against the reference model.
# Due to storage restrictions the individual null model draws cannot be archived
# in this repository. To enact using this data instead of re-computing null
# model draws and similarity measures, set this flag to TRUE (default). This
# flag will also lead to caches being created after PCNs have been computed in
# case no repository data has been found. Thus, deleting the cached data causes
# that a new cache will be created the next time this script is executed.
USE_REPO_DATA <- TRUE

FILE_PCN_CACHE <- file.path(DIR_CACHE,'pre_processed_pcn_data.RData')


### if applicable read network measures from disk instead of re-computing
#OVERWRITE_PCN_MEAS <- FALSE

## flag to write similarity measures to csv for external processing
SAVE_CORR_TO_CSV <- FALSE


# script flag preventing script to be executed if this has been done in the session before
# overwrite STATUS_MODEL_PCNS_RUN <- FALSE in global environment to force re-execution
if (!exists('STATUS_MODEL_PCNS_RUN')) STATUS_MODEL_PCNS_RUN <- FALSE

if (STATUS_MODEL_PCNS_RUN == FALSE) {
  message('main_pcns.R has not been run in session, executing script with the defined flags')
} else {
  message('main_pcns.R has been run in session, aborting executing the script. To force execution, set STATUS_MODEL_PCNS_RUN <- FALSE in global environment')
  invokeRestart('abort')
}


# if pre-processed simulation data should be used
if (USE_REPO_DATA) {
  #if (USE_MODEL_DATA_GLOBAL_ENV) {
  #  status_model_data_global_env <- all(sapply(c('ACERtrc', 'ACERtrcmil', 'ACERlc', 'ACERlcmil', 'ACERhcm'), exists))
  #}
  if (file.exists(FILE_PCN_CACHE)) {#(status_model_data_global_env == FALSE & file.exists(FILE_PCN_CACHE)) {
    message('Using pre-processed PCN data, therefore skipping computing them from scratch')
    load(FILE_PCN_CACHE)
    STATUS_MODEL_PCNS_RUN <- TRUE
    ## this loads main data objects into global environment:
    ## `ACERtrc` - TraCE21ka simulation
    ## `ACERlc` - LOVECLIM Last Deglaciation simulation
    ## `ACERhcm` - HadCM3 BRIDGE simulation of the Last Glacial
    invokeRestart('abort')
  #} else if (status_model_data_global_env) {
  #  message('Pseudo proxies are in global environment already. ')
  } else {
    message('USE_REPO_DATA == TRUE but no matching .RData file found, attempting to compute PCNs from scratch. This will take a while.')
  }
}


## Initialize all simulation and proxy data (per default, skipped automatically in `main_simulation_data.R` if this has been done before) ----
## also skipped if PCN data has been loaded from .RData cache because this cache contains the pseudo proxies as well.
## (Note that skipping is implicit by `invokeRestart('abort')` above)
## sourcing initializes main data objects for pseudo proxies from simulations:
## `ACERtrc`, `ACERtrcmil` - TraCE21ka simulation
## `ACERlc`, `ACERlcmil` - LOVECLIM Last Deglaciation simulation
## `ACERhcm` - HadCM3 BRIDGE simulation of the Last Glacial
source('main_simulation_data.R')


## Compute/load all PCNs ----

### Raw Trace21ka + ACER AP networks in 6-22 ka BP ----
ACERhere <- PCNdata() %>% 
  mix_sample_dating_and_err() %>% 
  harmonized_pollen_to_db() %>% 
  arboreal_pollen_to_db()

# generate hash for record set
hsh <- gen_site_hash(ACERtrc[[1]]$sites, ACERtrc[[1]]$sample_dating)

# generate/load null model and compute/load correlations
invisible(
  lapply(1:length(ACERtrc), function(i) {
    # load/generate null models draws
    ACERtrc[[i]] <<- ACERtrc[[i]] %>% 
      set_hybrid_icecores_to_db(n = 1000, write_trial = list(activate = DO_CACHEING_PCN, mode = 'slice100'), write_noise = list(activate = DO_CACHEING_PCN, mode = 'slice100'),
                                            use_cached_noise = TRUE, params = list(noise_params = list(snvr = 0,
                                                                                                       temporal = 'ar1',
                                                                                                       temporal_params = list(alpha = 0.8),
                                                                                                       spatial = 'matern',
                                                                                                       spatial_params = list(theta = 1000, smoothness = 1.5)),
                                                                                   sampling_params = list(sampling_mode = 'box',
                                                                                                          polamp = 4)))
    # compute distribution of similarity measure for the set of null model draws
    ACERtrc[[i]] <<- ACERtrc[[i]] %>% 
      corr_list_from_set_to_db(orig_data = paste0('hic_set_box-4-1000-', 0, '-ar1-0.8-matern-1000-1.5-slice100'),
                                           type_data = 'hybrid_ice_core', 
                                           transform = 'identity',
                                           detrend = list(activate = FALSE),
                                           hres_only = TRUE, interp = list(method = 'gaussian_kernel'), sd_one = FALSE,
                                           write_corr = list(activate = DO_CACHEING_PCN, dir = file.path(DIR_NM_DATA, 'correlations', '6_22')),
                                           windows = WINDOWS_TRACE, labels = WINDOW_LABELS_TRACE, ncores = NCORES_PCN) 
    
    # dump the null model draws because they just take up space
    # keep the similarity measurements
    ACERtrc[[i]]$`hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100` <<- NULL
  })
)

# compute/load correlations of simulation pseudo-proxies
# and test significance/load significance test results
invisible(
  lapply(1:length(ACERtrc), function(i) {
    # compute the similarity of the data
    ACERtrc[[i]] <<- ACERtrc[[i]] %>% 
      corr_list_to_db(orig_data = 'arb_pollen_data', 
                                  type_data = 'arboreal_pollen',
                                  transform = 'identity',
                                  detrend = list(activate = FALSE),
                                  hres_only = TRUE,
                                  interp = list(method = 'gaussian_kernel'),
                                  sd_one = FALSE, windows = WINDOWS_TRACE, labels = WINDOW_LABELS_TRACE,
                                  in_place = TRUE)
    # evaluate against the similarity measure distribution of the null model
    ACERtrc[[i]] <<- ACERtrc[[i]] %>% 
      corr_set_summary_to_db(orig_data = 'corr_list_gaussian_kernel_hres_identity_arb_pollen_data', 
                                         noise_data = 'corr_list_gaussian_kernel_hres_identity_hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100', test_option = 'ratio')
    # drop the similarity measurements of the null model draws
    # keep the statistics
    ACERtrc[[i]]$`corr_list_gaussian_kernel_hres_identity_hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100` <<- NULL
  })
)

# same steps for ACER proxies
ACERhere <- ACERhere %>% 
  set_hybrid_icecores_to_db(n = 1000, write_trial = list(activate = DO_CACHEING_PCN, mode = 'slice100'), write_noise = list(activate = DO_CACHEING_PCN, mode = 'slice100'),
                               use_cached_noise = TRUE, params = list(noise_params = list(snvr = 0,
                                                                                          temporal = 'ar1',
                                                                                          temporal_params = list(alpha = 0.8),
                                                                                          spatial = 'matern',
                                                                                          spatial_params = list(theta = 1000, smoothness = 1.5)),
                                                                      sampling_params = list(sampling_mode = 'box',
                                                                                             polamp = 4)))
ACERhere <- ACERhere %>% 
  corr_list_from_set_to_db(orig_data = paste0('hic_set_box-4-1000-', 0, '-ar1-0.8-matern-1000-1.5-slice100'),
                              type_data = 'hybrid_ice_core', 
                              transform = 'identity',
                              detrend = list(activate = FALSE),
                              hres_only = TRUE, interp = list(method = 'gaussian_kernel'), sd_one = FALSE,
                              write_corr = list(activate = DO_CACHEING_PCN, dir = file.path(DIR_NM_DATA, 'correlations', '6_22')),
                              windows = WINDOWS_TRACE, labels = WINDOW_LABELS_TRACE)
ACERhere$`hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100` <- NULL
ACERhere <- ACERhere %>% 
  corr_list_to_db(orig_data = 'arb_pollen_data', 
                     type_data = 'arboreal_pollen',
                     transform = 'probit',
                     detrend = list(activate = FALSE),
                     hres_only = TRUE,
                     interp = list(method = 'gaussian_kernel'),
                     sd_one = FALSE, windows = WINDOWS_TRACE, labels = WINDOW_LABELS_TRACE,
                     in_place = TRUE)
ACERhere <- ACERhere %>% 
  corr_set_summary_to_db(orig_data = 'corr_list_gaussian_kernel_hres_probit_arb_pollen_data', 
                            noise_data = 'corr_list_gaussian_kernel_hres_identity_hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100', test_option = 'ratio')
ACERhere$`corr_list_gaussian_kernel_hres_identity_hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100` <- NULL

ACERtrc[['ACER_ap']] <- ACERhere
rm(ACERhere)


### Millennial-scale Trace21ka PR, TS + ACER AP and icecore networks in 6-22 ka BP ----
### Note that for network measures links to the icecores are filtered out
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


# generate hash for record set
hsh <- gen_site_hash(ACERtrcmil[[1]]$sites, ACERtrcmil[[1]]$sample_dating)

# generate/load null model and compute/load correlations
# note that Gaussian filtering is enabled now for 2-8ka window
invisible(
  lapply(1:length(ACERtrcmil), function(i) {
    ACERtrcmil[[i]] <<- ACERtrcmil[[i]] %>% 
      set_hybrid_icecores_to_db(n = 1000, write_trial = list(activate = DO_CACHEING_PCN, mode = 'slice100'), write_noise = list(activate = DO_CACHEING_PCN, mode = 'slice100'),
                                            use_cached_noise = TRUE, params = list(noise_params = list(snvr = 0,
                                                                                                       temporal = 'ar1',
                                                                                                       temporal_params = list(alpha = 0.8),
                                                                                                       spatial = 'matern',
                                                                                                       spatial_params = list(theta = 190, smoothness = 1.5)),
                                                                                   sampling_params = list(sampling_mode = 'box',
                                                                                                          polamp = 4)), ncores = NCORES_PCN)
    ACERtrcmil[[i]] <<- ACERtrcmil[[i]] %>% 
      corr_list_from_set_to_db(orig_data = paste0('hic_set_box-4-1000-', 0, '-ar1-0.8-matern-190-1.5-slice100'),
                                           type_data = 'hybrid_ice_core', 
                                           transform = 'identity',
                                           detrend = list(activate = TRUE, method = 'gaussbp', args = list(per1 = 2000, per2 = 8000)),
                                           hres_only = TRUE, interp = list(method = 'gaussian_kernel'), sd_one = FALSE,
                                           write_corr = list(activate = DO_CACHEING_PCN, dir = file.path(DIR_NM_DATA, 'correlations', '6_22')),
                                           windows = WINDOWS_TRACE, labels = WINDOW_LABELS_TRACE, ncores = NCORES_PCN)
    
    ACERtrcmil[[i]]$`hic_set_box-4-1000-0-ar1-0.8-matern-190-1.5-slice100` <<- NULL
  })
)


# compute/load correlations of simulation pseudo-proxies
# and test significance/load significance test results
invisible(
  lapply(1:length(ACERtrcmil), function(i) {
    ACERtrcmil[[i]] <<- ACERtrcmil[[i]] %>% 
      corr_list_to_db(orig_data = 'arb_pollen_data', 
                                  type_data = 'arboreal_pollen',
                                  transform = 'identity',
                                  detrend = list(activate = TRUE, method = 'gaussbp', args = list(per1 = 2000, per2 = 8000)),
                                  hres_only = TRUE,
                                  interp = list(method = 'gaussian_kernel'),
                                  sd_one = FALSE, windows = WINDOWS_TRACE, labels = WINDOW_LABELS_TRACE,
                                  in_place = TRUE)
    ACERtrcmil[[i]] <<- ACERtrcmil[[i]] %>% 
      corr_set_summary_to_db(orig_data = 'corr_list_gaussian_kernel_hres_detrend_gaussbp2000-8000_identity_arb_pollen_data', 
                             noise_data = 'corr_list_gaussian_kernel_hres_detrend_gaussbp2000-8000_identity_hic_set_box-4-1000-0-ar1-0.8-matern-190-1.5-slice100', test_option = 'ratio')
    ACERtrcmil[[i]]$`corr_list_gaussian_kernel_hres_detrend_gaussbp2000-8000_identity_hic_set_box-4-1000-0-ar1-0.8-matern-190-1.5-slice100` <<- NULL
  })
)


# same steps for ACER proxies
ACERhere <- ACERhere %>% 
  set_hybrid_icecores_to_db(n = 1000, write_trial = list(activate = DO_CACHEING_PCN, mode = 'slice100'), write_noise = list(activate = DO_CACHEING_PCN, mode = 'slice100'),
                               use_cached_noise = TRUE, params = list(noise_params = list(snvr = 0,
                                                                                          temporal = 'ar1',
                                                                                          temporal_params = list(alpha = 0.8),
                                                                                          spatial = 'matern',
                                                                                          spatial_params = list(theta = 190, smoothness = 1.5)),
                                                                      sampling_params = list(sampling_mode = 'box',
                                                                                             polamp = 4)))
ACERhere <- ACERhere %>% 
  corr_list_from_set_to_db(orig_data = paste0('hic_set_box-4-1000-', 0, '-ar1-0.8-matern-190-1.5-slice100'),
                              type_data = 'hybrid_ice_core', 
                              transform = 'identity',
                              detrend = list(activate = TRUE, method = 'gaussbp', args = list(per1 = 2000, per2 = 8000)),
                              hres_only = TRUE, interp = list(method = 'gaussian_kernel'), sd_one = FALSE,
                              write_corr = list(activate = DO_CACHEING_PCN, dir = file.path(DIR_NM_DATA, 'correlations', '6_22')),
                              windows = WINDOWS_TRACE, labels = WINDOW_LABELS_TRACE)
ACERhere$`hic_set_box-4-1000-0-ar1-0.8-matern-190-1.5-slice100` <- NULL
ACERhere <- ACERhere %>% 
  corr_list_to_db(orig_data = 'arb_pollen_data', 
                     type_data = 'arboreal_pollen',
                     transform = 'probit',
                     detrend = list(activate = TRUE, method = 'gaussbp', args = list(per1 = 2000, per2 = 8000)),
                     hres_only = TRUE,
                     interp = list(method = 'gaussian_kernel'),
                     sd_one = FALSE, windows = WINDOWS_TRACE, labels = WINDOW_LABELS_TRACE,
                     in_place = TRUE)
ACERhere <- ACERhere %>% 
  corr_set_summary_to_db(orig_data = 'corr_list_gaussian_kernel_hres_detrend_gaussbp2000-8000_probit_arb_pollen_data', 
                            noise_data = 'corr_list_gaussian_kernel_hres_detrend_gaussbp2000-8000_identity_hic_set_box-4-1000-0-ar1-0.8-matern-190-1.5-slice100', test_option = 'ratio')
ACERhere$`corr_list_gaussian_kernel_hres_detrend_gaussbp2000-8000_identity_hic_set_box-4-1000-0-ar1-0.8-matern-190-1.5-slice100` <- NULL

ACERtrcmil[['ACER_ap']] <- ACERhere
rm(ACERhere)


### Raw LOVECLIM + ACER AP networks in 6.2-18 ka BP ----
ACERhere <- PCNdata() %>% 
  mix_sample_dating_and_err() %>% 
  harmonized_pollen_to_db() %>% 
  arboreal_pollen_to_db()

# generate hash for record set
hsh <- gen_site_hash(ACERlc[[1]]$sites, ACERlc[[1]]$sample_dating)

# generate/load null model and compute/load correlations
invisible(
  lapply(1:length(ACERlc), function(i) {
    ACERlc[[i]] <<- ACERlc[[i]] %>% 
      set_hybrid_icecores_to_db(n = 1000, write_trial = list(activate = DO_CACHEING_PCN, mode = 'slice100'), write_noise = list(activate = DO_CACHEING_PCN, mode = 'slice100'),
                                          use_cached_noise = TRUE, params = list(noise_params = list(snvr = 0,
                                                                                                     temporal = 'ar1',
                                                                                                     temporal_params = list(alpha = 0.8),
                                                                                                     spatial = 'matern',
                                                                                                     spatial_params = list(theta = 1000, smoothness = 1.5)),
                                                                                 sampling_params = list(sampling_mode = 'box',
                                                                                                        polamp = 4)), ncores = NCORES_PCN)
    ACERlc[[i]] <<- ACERlc[[i]] %>% 
      corr_list_from_set_to_db(orig_data = paste0('hic_set_box-4-1000-', 0, '-ar1-0.8-matern-1000-1.5-slice100'),
                                         type_data = 'hybrid_ice_core', 
                                         transform = 'identity',
                                         detrend = list(activate = FALSE),
                                         hres_only = TRUE, interp = list(method = 'gaussian_kernel'), sd_one = FALSE,
                                         write_corr = list(activate = DO_CACHEING_PCN, dir = file.path(DIR_NM_DATA, 'correlations', '6.2-18')),
                                         windows = WINDOWS_LC, labels = WINDOW_LABELS_LC, ncores = NCORES_PCN)
    
    ACERlc[[i]]$`hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100` <<- NULL
  })
)

# compute/load correlations of simulation pseudo-proxies
# and test significance/load significance test results
invisible(
  lapply(1:length(ACERlc), function(i) {
    ACERlc[[i]] <<- ACERlc[[i]] %>% 
      corr_list_to_db(orig_data = 'arb_pollen_data', 
                                type_data = 'arboreal_pollen',
                                transform = 'identity',
                                detrend = list(activate = FALSE),
                                hres_only = TRUE,
                                interp = list(method = 'gaussian_kernel'),
                                sd_one = FALSE, windows = WINDOWS_LC, labels = WINDOW_LABELS_LC,
                                in_place = TRUE)
    ACERlc[[i]] <<- ACERlc[[i]] %>% 
      corr_set_summary_to_db(orig_data = 'corr_list_gaussian_kernel_hres_identity_arb_pollen_data', 
                             noise_data = 'corr_list_gaussian_kernel_hres_identity_hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100', test_option = 'ratio')
    ACERlc[[i]]$`corr_list_gaussian_kernel_hres_identity_hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100` <<- NULL
  })
)

# same steps for ACER proxies
ACERhere <- ACERhere %>% 
  set_hybrid_icecores_to_db(n = 1000, write_trial = list(activate = DO_CACHEING_PCN, mode = 'slice100'), write_noise = list(activate = DO_CACHEING_PCN, mode = 'slice100'),
                               use_cached_noise = TRUE, params = list(noise_params = list(snvr = 0,
                                                                                          temporal = 'ar1',
                                                                                          temporal_params = list(alpha = 0.8),
                                                                                          spatial = 'matern',
                                                                                          spatial_params = list(theta = 1000, smoothness = 1.5)),
                                                                      sampling_params = list(sampling_mode = 'box',
                                                                                             polamp = 4)))
ACERhere <- ACERhere %>% 
  corr_list_from_set_to_db(orig_data = paste0('hic_set_box-4-1000-', 0, '-ar1-0.8-matern-1000-1.5-slice100'),
                              type_data = 'hybrid_ice_core', 
                              transform = 'identity',
                              detrend = list(activate = FALSE),
                              hres_only = TRUE, interp = list(method = 'gaussian_kernel'), sd_one = FALSE,
                              write_corr = list(activate = DO_CACHEING_PCN, dir = file.path(DIR_NM_DATA, 'correlations', '6.2-18')),
                              windows = WINDOWS_LC, labels = WINDOW_LABELS_LC)
ACERhere$`hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100` <- NULL
ACERhere <- ACERhere %>% 
  corr_list_to_db(orig_data = 'arb_pollen_data', 
                     type_data = 'arboreal_pollen',
                     transform = 'probit',
                     detrend = list(activate = FALSE),
                     hres_only = TRUE,
                     interp = list(method = 'gaussian_kernel'),
                     sd_one = FALSE, windows = WINDOWS_LC, labels = WINDOW_LABELS_LC,
                     in_place = TRUE)
ACERhere <- ACERhere %>% 
  corr_set_summary_to_db(orig_data = 'corr_list_gaussian_kernel_hres_probit_arb_pollen_data', 
                            noise_data = 'corr_list_gaussian_kernel_hres_identity_hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100', test_option = 'ratio')
ACERhere$`corr_list_gaussian_kernel_hres_identity_hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100` <- NULL

ACERlc[['ACER_ap']] <- ACERhere
rm(ACERhere)


### Millennial-scale LOVECLIM + ACER AP networks in 6.2-18 ka BP ----
ACERhere <- PCNdata() %>% 
  mix_sample_dating_and_err() %>% 
  harmonized_pollen_to_db() %>% 
  arboreal_pollen_to_db()

# generate hash for record set
hsh <- gen_site_hash(ACERlcmil[[1]]$sites, ACERlcmil[[1]]$sample_dating)

# generate/load null model and compute/load correlations
invisible(
  lapply(1:length(ACERlcmil), function(i) {
    ACERlcmil[[i]] <<- ACERlcmil[[i]] %>% 
      set_hybrid_icecores_to_db(n = 1000, write_trial = list(activate = DO_CACHEING_PCN, mode = 'slice100'), write_noise = list(activate = DO_CACHEING_PCN, mode = 'slice100'),
                                          use_cached_noise = TRUE, params = list(noise_params = list(snvr = 0,
                                                                                                     temporal = 'ar1',
                                                                                                     temporal_params = list(alpha = 0.8),
                                                                                                     spatial = 'matern',
                                                                                                     spatial_params = list(theta = 190, smoothness = 1.5)),
                                                                                 sampling_params = list(sampling_mode = 'box',
                                                                                                        polamp = 4)), ncores = NCORES_PCN)
    ACERlcmil[[i]] <<- ACERlcmil[[i]] %>% 
      corr_list_from_set_to_db(orig_data = paste0('hic_set_box-4-1000-', 0, '-ar1-0.8-matern-190-1.5-slice100'),
                                         type_data = 'hybrid_ice_core', 
                                         transform = 'identity',
                                         detrend = list(activate = TRUE, method = 'gaussbp', args = list(per1 = 2000, per2 = 8000)),
                                         hres_only = TRUE, interp = list(method = 'gaussian_kernel'), sd_one = FALSE,
                                         write_corr = list(activate = DO_CACHEING_PCN, dir = file.path(DIR_NM_DATA, 'correlations', '6.2-18')),
                                         windows = WINDOWS_LC, labels = WINDOW_LABELS_LC, ncores = NCORES_PCN)
    
    ACERlcmil[[i]]$`hic_set_box-4-1000-0-ar1-0.8-matern-190-1.5-slice100` <<- NULL
  })
)

# compute/load correlations of simulation pseudo-proxies
# and test significance/load significance test results
invisible(
  lapply(1:length(ACERlcmil), function(i) {
    ACERlcmil[[i]] <<- ACERlcmil[[i]] %>% 
      corr_list_to_db(orig_data = 'arb_pollen_data', 
                                type_data = 'arboreal_pollen',
                                transform = 'identity',
                                detrend = list(activate = TRUE, method = 'gaussbp', args = list(per1 = 2000, per2 = 8000)),
                                hres_only = TRUE,
                                interp = list(method = 'gaussian_kernel'),
                                sd_one = FALSE, windows = WINDOWS_LC, labels = WINDOW_LABELS_LC,
                                in_place = TRUE)
    ACERlcmil[[i]] <<- ACERlcmil[[i]] %>% 
      corr_set_summary_to_db(orig_data = 'corr_list_gaussian_kernel_hres_detrend_gaussbp2000-8000_identity_arb_pollen_data', 
                             noise_data = 'corr_list_gaussian_kernel_hres_detrend_gaussbp2000-8000_identity_hic_set_box-4-1000-0-ar1-0.8-matern-190-1.5-slice100', test_option = 'ratio')
    ACERlcmil[[i]]$`corr_list_gaussian_kernel_hres_detrend_gaussbp2000-8000_identity_hic_set_box-4-1000-0-ar1-0.8-matern-190-1.5-slice100` <<- NULL
  })
)

# same steps for ACER proxies
ACERhere <- ACERhere %>% 
  set_hybrid_icecores_to_db(n = 1000, write_trial = list(activate = DO_CACHEING_PCN, mode = 'slice100'), write_noise = list(activate = DO_CACHEING_PCN, mode = 'slice100'),
                               use_cached_noise = TRUE, params = list(noise_params = list(snvr = 0,
                                                                                          temporal = 'ar1',
                                                                                          temporal_params = list(alpha = 0.8),
                                                                                          spatial = 'matern',
                                                                                          spatial_params = list(theta = 190, smoothness = 1.5)),
                                                                      sampling_params = list(sampling_mode = 'box',
                                                                                             polamp = 4)))
ACERhere <- ACERhere %>% 
  corr_list_from_set_to_db(orig_data = paste0('hic_set_box-4-1000-', 0, '-ar1-0.8-matern-190-1.5-slice100'),
                              type_data = 'hybrid_ice_core', 
                              transform = 'identity',
                              detrend = list(activate = TRUE, method = 'gaussbp', args = list(per1 = 2000, per2 = 8000)),
                              hres_only = TRUE, interp = list(method = 'gaussian_kernel'), sd_one = FALSE,
                              write_corr = list(activate = DO_CACHEING_PCN, dir = file.path(DIR_NM_DATA, 'correlations', '6.2-18')),
                              windows = WINDOWS_LC, labels = WINDOW_LABELS_LC)
ACERhere$`hic_set_box-4-1000-0-ar1-0.8-matern-190-1.5-slice100` <- NULL
ACERhere <- ACERhere %>% 
  corr_list_to_db(orig_data = 'arb_pollen_data', 
                     type_data = 'arboreal_pollen',
                     transform = 'probit',
                     detrend = list(activate = TRUE, method = 'gaussbp', args = list(per1 = 2000, per2 = 8000)),
                     hres_only = TRUE,
                     interp = list(method = 'gaussian_kernel'),
                     sd_one = FALSE, windows = WINDOWS_LC, labels = WINDOW_LABELS_LC,
                     in_place = TRUE)
ACERhere <- ACERhere %>% 
  corr_set_summary_to_db(orig_data = 'corr_list_gaussian_kernel_hres_detrend_gaussbp2000-8000_probit_arb_pollen_data', 
                         noise_data = 'corr_list_gaussian_kernel_hres_detrend_gaussbp2000-8000_identity_hic_set_box-4-1000-0-ar1-0.8-matern-190-1.5-slice100', test_option = 'ratio')
ACERhere$`corr_list_gaussian_kernel_hres_detrend_gaussbp2000-8000_identity_hic_set_box-4-1000-0-ar1-0.8-matern-190-1.5-slice100` <- NULL

ACERlcmil[['ACER_ap']] <- ACERhere
rm(ACERhere)


### Raw HadCM3 + ACER AP networks in 6-22 ka BP ----
ACERhere <- PCNdata() %>% 
  mix_sample_dating_and_err() %>% 
  harmonized_pollen_to_db() %>% 
  arboreal_pollen_to_db()

# generate hash for record set
hsh <- gen_site_hash(ACERhcm[[1]]$sites, ACERhcm[[1]]$sample_dating)

# generate/load null model and compute/load correlations
invisible(
  lapply(1:length(ACERhcm), function(i) {
    ACERhcm[[i]] <<- ACERhcm[[i]] %>% 
      set_hybrid_icecores_to_db(n = 1000,
                                write_trial = list(activate = DO_CACHEING_PCN, mode = 'slice100'),
                                write_noise = list(activate = DO_CACHEING_PCN, mode = 'slice100'),
                                use_cached_noise = TRUE, params = list(noise_params = list(snvr = 0,
                                                                                           temporal = 'ar1',
                                                                                           temporal_params = list(alpha = 0.8),
                                                                                           spatial = 'matern',
                                                                                           spatial_params = list(theta = 1000, smoothness = 1.5)),
                                                                       sampling_params = list(sampling_mode = 'box',
                                                                                              polamp = 4)), ncores = NCORES_PCN)
    ACERhcm[[i]] <<- ACERhcm[[i]] %>% 
      corr_list_from_set_to_db(orig_data = paste0('hic_set_box-4-1000-', 0, '-ar1-0.8-matern-1000-1.5-slice100'),
                               type_data = 'hybrid_ice_core', 
                               transform = 'identity',
                               detrend = list(activate = FALSE),
                               hres_only = TRUE, interp = list(method = 'gaussian_kernel'), sd_one = FALSE,
                               write_corr = list(activate = DO_CACHEING_PCN, dir = file.path(DIR_NM_DATA, 'correlations', '6_22')),
                               windows = WINDOWS_BBC, labels = WINDOW_LABELS_BBC, ncores = NCORES_PCN)
    
    ACERhcm[[i]]$`hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100` <<- NULL
  })
)

# compute/load correlations of simulation pseudo-proxies
# and test significance/load significance test results
invisible(
  lapply(1:length(ACERhcm), function(i) {
    ACERhcm[[i]] <<- ACERhcm[[i]] %>% 
      corr_list_to_db(orig_data = 'arb_pollen_data', 
                                 type_data = 'arboreal_pollen',
                                 transform = 'identity',
                                 detrend = list(activate = FALSE),
                                 hres_only = TRUE,
                                 interp = list(method = 'gaussian_kernel'),
                                 sd_one = FALSE, windows = WINDOWS_BBC, labels = WINDOW_LABELS_BBC,
                                 in_place = TRUE)
    ACERhcm[[i]] <<- ACERhcm[[i]] %>% 
      corr_set_summary_to_db(orig_data = 'corr_list_gaussian_kernel_hres_identity_arb_pollen_data', 
                                        noise_data = 'corr_list_gaussian_kernel_hres_identity_hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100', test_option = 'ratio')
    ACERhcm[[i]]$`corr_list_gaussian_kernel_hres_identity_hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100` <<- NULL
  })
)

# same steps for ACER proxies
ACERhere <- ACERhere %>% 
  set_hybrid_icecores_to_db(n = 1000, write_trial = list(activate = DO_CACHEING_PCN, mode = 'slice100'), write_noise = list(activate = DO_CACHEING_PCN, mode = 'slice100'),
                                    use_cached_noise = TRUE, params = list(noise_params = list(snvr = 0,
                                                                                               temporal = 'ar1',
                                                                                               temporal_params = list(alpha = 0.8),
                                                                                               spatial = 'matern',
                                                                                               spatial_params = list(theta = 1000, smoothness = 1.5)),
                                                                           sampling_params = list(sampling_mode = 'box',
                                                                                                  polamp = 4)))
ACERhere <- ACERhere %>% 
  corr_list_from_set_to_db(orig_data = paste0('hic_set_box-4-1000-', 0, '-ar1-0.8-matern-1000-1.5-slice100'),
                                   type_data = 'hybrid_ice_core', 
                                   transform = 'identity',
                                   detrend = list(activate = FALSE),
                                   hres_only = TRUE, interp = list(method = 'gaussian_kernel'), sd_one = FALSE,
                                   write_corr = list(activate = DO_CACHEING_PCN, dir = file.path(DIR_NM_DATA, 'correlations', '6_22')),
                                   windows = WINDOWS_BBC, labels = WINDOW_LABELS_BBC)
ACERhere$`hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100` <- NULL
ACERhere <- ACERhere %>% 
  corr_list_to_db(orig_data = 'arb_pollen_data', 
                          type_data = 'arboreal_pollen',
                          transform = 'probit',
                          detrend = list(activate = FALSE),
                          hres_only = TRUE,
                          interp = list(method = 'gaussian_kernel'),
                          sd_one = FALSE, windows = WINDOWS_BBC, labels = WINDOW_LABELS_BBC,
                          in_place = TRUE)
ACERhere <- ACERhere %>% 
  corr_set_summary_to_db(orig_data = 'corr_list_gaussian_kernel_hres_probit_arb_pollen_data', 
                                 noise_data = 'corr_list_gaussian_kernel_hres_identity_hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100', test_option = 'ratio')
ACERhere$`corr_list_gaussian_kernel_hres_identity_hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100` <- NULL

ACERhcm[['ACER_ap']] <- ACERhere
rm(ACERhere)

# garbage collect for parallel processes/file loading
gc()

#### end PCN computation


# cache PCNs if permitted ----
if (USE_REPO_DATA) {
  save(
    ACERtrc, ACERtrcmil, ACERlc, ACERlcmil, ACERhcm,
    file = FILE_PCN_CACHE
  )
}

# modifying run status of this script
STATUS_MODEL_PCNS_RUN <- TRUE
