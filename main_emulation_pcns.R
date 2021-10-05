# Preliminary wrapper script to compute
# - paleoclimate networks from emulated TSF response of TraCE and HadCM3
#   * combined forcing of mean annual temperature, mean annual precipitation, CO2 level
#   * single forcing experiments with one of MAT, MAP, CO2 and others fixed to sample
#     mean of respective pseudo-record
# - link densities of these networks
# as preparation for `create_fig_emul.R` (preliminary name) and `main_emulation_pcn_meas.R`.


# NOTE (as in `main_pcns.R`): 
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
USE_REPO_DATA_EMUL <- TRUE

FILE_EMUL_PCN_CACHE <- file.path(DIR_CACHE,'pre_processed_emul_pcn_data.RData')


## flag to write similarity measures to csv for external processing
SAVE_CORR_TO_CSV <- FALSE


# script flag preventing script to be executed if this has been done in the session before
# overwrite STATUS_EMUL_PCNS_RUN <- FALSE in global environment to force re-execution
if (!exists('STATUS_EMUL_PCNS_RUN')) STATUS_EMUL_PCNS_RUN <- FALSE

if (STATUS_EMUL_PCNS_RUN == FALSE) {
  message('main_emulation_pcns.R has not been run in session, executing script with the defined flags')
} else {
  message('main_emulation_pcns.R has been run in session, aborting executing the script. To force execution, set STATUS_EMUL_PCNS_RUN <- FALSE in global environment')
  invokeRestart('abort')
}


# if pre-processed simulation data should be used
if (USE_REPO_DATA_EMUL) {
  if (file.exists(FILE_EMUL_PCN_CACHE)) {#(status_model_data_global_env == FALSE & file.exists(FILE_EMUL_PCN_CACHE)) {
    message('Using pre-processed emulated surrogates PCN data, therefore skipping computing them from scratch')
    load(FILE_EMUL_PCN_CACHE)
    STATUS_EMUL_PCNS_RUN <- TRUE
    ## this loads main data objects into global environment:
    ## `ACERtrc_em` [$forceall, $forcets, $forcepr, $forceco2] - Emulated TraCE21ka response
    ## `ACERhcm_em` [$forceall, $forcets, $forcepr, $forceco2] - Emulated HadCM3 BRIDGE simulation
    invokeRestart('abort')
  } else {
    message('USE_REPO_DATA_EMUL == TRUE but no matching .RData file found, attempting to compute emulated surrogate PCNs from scratch')
  }
}


## Initialize all emulated data (per default, skipped automatically in `main_emulation_data.R` if this has been done before) ----
## also skipped if PCN data has been loaded from .RData cache because this cache contains the pseudo proxies as well.
## (Note that skipping is implicit by `invokeRestart('abort')` above)
## sourcing initializes main data objects for pseudo proxies from simulations:
## `ACERtrc_em` [$forceall, $forcets, $forcepr, $forceco2] - Emulated TraCE21ka response
## `ACERhcm_em` [$forceall, $forcets, $forcepr, $forceco2] - Emulated HadCM3 BRIDGE simulation
source('main_emulation_data.R')

#### compute PCNs for emulated DGVM experiments ----

# Emulated TraCE response in 6-22 ka BP ----
# generate hash for record set
hsh <- gen_site_hash(ACERtrc_em[[1]]$sites, ACERtrc_em[[1]]$sample_dating)

# generate/load null model and compute/load correlations
invisible(
  lapply(1:length(ACERtrc_em), function(i) {
    # load/generate null models draws
    ACERtrc_em[[i]] <<- ACERtrc_em[[i]] %>% 
      set_hybrid_icecores_to_db(n = 1000, write_trial = list(activate = DO_CACHEING_PCN, mode = 'slice100'), write_noise = list(activate = DO_CACHEING_PCN, mode = 'slice100'),
                                use_cached_noise = TRUE, params = list(noise_params = list(snvr = 0,
                                                                                           temporal = 'ar1',
                                                                                           temporal_params = list(alpha = 0.8),
                                                                                           spatial = 'matern',
                                                                                           spatial_params = list(theta = 1000, smoothness = 1.5)),
                                                                       sampling_params = list(sampling_mode = 'box',
                                                                                              polamp = 4)))
    # compute distribution of similarity measure for the set of null model draws
    ACERtrc_em[[i]] <<- ACERtrc_em[[i]] %>% 
      corr_list_from_set_to_db(orig_data = paste0('hic_set_box-4-1000-', 0, '-ar1-0.8-matern-1000-1.5-slice100'),
                               type_data = 'hybrid_ice_core', 
                               transform = 'identity',
                               detrend = list(activate = FALSE),
                               hres_only = TRUE, interp = list(method = 'gaussian_kernel'), sd_one = FALSE,
                               write_corr = list(activate = DO_CACHEING_PCN, dir = file.path(DIR_NM_DATA, 'correlations', '6_22')),
                               windows = WINDOWS_TRACE, labels = WINDOW_LABELS_TRACE, ncores = NCORES_PCN) 
    
    # dump the null model draws because they just take up space
    # keep the similarity measurements
    ACERtrc_em[[i]]$`hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100` <<- NULL
  })
)

# compute/load correlations of simulation pseudo-proxies
# and test significance/load significance test results
invisible(
  lapply(1:length(ACERtrc_em), function(i) {
    # compute the similarity of the data
    ACERtrc_em[[i]] <<- ACERtrc_em[[i]] %>% 
      corr_list_to_db(orig_data = 'arb_pollen_data', 
                      type_data = 'arboreal_pollen',
                      transform = 'identity',
                      detrend = list(activate = FALSE),
                      hres_only = TRUE,
                      interp = list(method = 'gaussian_kernel'),
                      sd_one = FALSE, windows = WINDOWS_TRACE, labels = WINDOW_LABELS_TRACE,
                      in_place = TRUE)
    # evaluate against the similarity measure distribution of the null model
    ACERtrc_em[[i]] <<- ACERtrc_em[[i]] %>% 
      corr_set_summary_to_db(orig_data = 'corr_list_gaussian_kernel_hres_identity_arb_pollen_data', 
                             noise_data = 'corr_list_gaussian_kernel_hres_identity_hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100', test_option = 'ratio')
    # drop the similarity measurements of the null model draws
    # keep the statistics
    ACERtrc_em[[i]]$`corr_list_gaussian_kernel_hres_identity_hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100` <<- NULL
  })
)


# Emulated HadCM3 response ----
# generate hash for record set
hsh <- gen_site_hash(ACERhcm_em[[1]]$sites, ACERhcm_em[[1]]$sample_dating)

# generate/load null model and compute/load correlations
invisible(
  lapply(1:length(ACERhcm_em), function(i) {
    ACERhcm_em[[i]] <<- ACERhcm_em[[i]] %>% 
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
    ACERhcm_em[[i]] <<- ACERhcm_em[[i]] %>% 
      corr_list_from_set_to_db(orig_data = paste0('hic_set_box-4-1000-', 0, '-ar1-0.8-matern-1000-1.5-slice100'),
                               type_data = 'hybrid_ice_core', 
                               transform = 'identity',
                               detrend = list(activate = FALSE),
                               hres_only = TRUE, interp = list(method = 'gaussian_kernel'), sd_one = FALSE,
                               write_corr = list(activate = DO_CACHEING_PCN, dir = file.path(DIR_NM_DATA, 'correlations', '6_22')),
                               windows = WINDOWS_BBC, labels = WINDOW_LABELS_BBC, ncores = NCORES_PCN)
    
    ACERhcm_em[[i]]$`hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100` <<- NULL
  })
)

# compute/load correlations of simulation pseudo-proxies
# and test significance/load significance test results
invisible(
  lapply(1:length(ACERhcm_em), function(i) {
    ACERhcm_em[[i]] <<- ACERhcm_em[[i]] %>% 
      corr_list_to_db(orig_data = 'arb_pollen_data', 
                      type_data = 'arboreal_pollen',
                      transform = 'identity',
                      detrend = list(activate = FALSE),
                      hres_only = TRUE,
                      interp = list(method = 'gaussian_kernel'),
                      sd_one = FALSE, windows = WINDOWS_BBC, labels = WINDOW_LABELS_BBC,
                      in_place = TRUE)
    ACERhcm_em[[i]] <<- ACERhcm_em[[i]] %>% 
      corr_set_summary_to_db(orig_data = 'corr_list_gaussian_kernel_hres_identity_arb_pollen_data', 
                             noise_data = 'corr_list_gaussian_kernel_hres_identity_hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100', test_option = 'ratio')
    ACERhcm_em[[i]]$`corr_list_gaussian_kernel_hres_identity_hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100` <<- NULL
  })
)


# garbage collect for parallel processes/file loading
gc()

#### end PCN computation

# cache PCNs if permitted ----
if (USE_REPO_DATA_EMUL) {
  save(
    ACERtrc_em,
    ACERhcm_em,
    file = FILE_EMUL_PCN_CACHE
  )
}

# modifying run status of this script
STATUS_EMUL_PCNS_RUN <- TRUE
