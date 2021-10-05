# Wrapper script for a small demonstration of the routines used to generate the statistical null model
# and to compute the paleoclimate networks
# Note, that results have no significance given the small number of null model draws (see `get_started.Rmd` and the Gelman-Rubin criterion).

# remember to initialize packages and constants using `source('init_all.R)`

# initialize ACER data
ACERhere <- PCNdata() %>% 
  mix_sample_dating_and_err() %>% 
  harmonized_pollen_to_db() %>% 
  arboreal_pollen_to_db()

# generate hash for record set
hsh <- gen_site_hash(ACERhere$sites, ACERhere$sample_dating)

# generate/load noise drawn with the spatio-temporal properties of the null model
## make sure that you have created the directories
## dir.create(file.path(DIR_NM_DATA, 'hic'), recursive = TRUE)
## dir.create(file.path(DIR_NM_DATA, 'noise'), recursive = TRUE)
ACERhere <- ACERhere %>% 
  set_hybrid_icecores_to_db(n = 4, write_trial = list(activate = TRUE, mode = 'slice2'), write_noise = list(activate = TRUE, mode = 'slice2'),
                            use_cached_noise = TRUE, params = list(noise_params = list(snvr = 0,
                                                                                       temporal = 'ar1',
                                                                                       temporal_params = list(alpha = 0.8),
                                                                                       spatial = 'matern',
                                                                                       spatial_params = list(theta = 1000, smoothness = 1.5)),
                                                                   sampling_params = list(sampling_mode = 'box',
                                                                                          polamp = 4)))

# calculate correlation ensemble for the individual draws from the null model
## make sure that you have created the directories
## dir.create(file.path(DIR_NM_DATA, 'correlations'), recursive = TRUE)
ACERhere <- ACERhere %>% 
  corr_list_from_set_to_db(orig_data = paste0('hic_set_box-4-4-', 0, '-ar1-0.8-matern-1000-1.5-slice2'),
                           type_data = 'hybrid_ice_core', 
                           transform = 'identity',
                           detrend = list(activate = FALSE),
                           hres_only = TRUE, interp = list(method = 'gaussian_kernel'), sd_one = FALSE,
                           write_corr = list(activate = TRUE, dir = file.path(DIR_NM_DATA, 'correlations', '6_22')),
                           windows = WINDOWS_TRACE, labels = WINDOW_LABELS_TRACE)

# remove the time series drawn from the null model
# as similarity estimates are compared to the correlation ensemble
ACERhere$`hic_set_box-4-4-0-ar1-0.8-matern-1000-1.5-slice2` <- NULL

# compute the correlations of the pollen proxies
ACERhere <- ACERhere %>% 
  corr_list_to_db(orig_data = 'arb_pollen_data', 
                  type_data = 'arboreal_pollen',
                  transform = 'probit',
                  detrend = list(activate = FALSE),
                  hres_only = TRUE,
                  interp = list(method = 'gaussian_kernel'),
                  sd_one = FALSE, windows = WINDOWS_TRACE, labels = WINDOW_LABELS_TRACE,
                  in_place = TRUE)

# quantify the significance of similarity estimates
ACERhere <- ACERhere %>% 
  corr_set_summary_to_db(orig_data = 'corr_list_gaussian_kernel_hres_probit_arb_pollen_data', 
                         noise_data = 'corr_list_gaussian_kernel_hres_identity_hic_set_box-4-4-0-ar1-0.8-matern-1000-1.5-slice2', test_option = 'ratio')

# remove the correlation ensemble from memory because it is not needed anymore
ACERhere$`corr_list_gaussian_kernel_hres_identity_hic_set_box-4-4-0-ar1-0.8-matern-1000-1.5-slice2` <- NULL

# the similarity estimates and the results of the significance test are contained here
# (Note that the results are not meaningful, given the low number of null model draws)
ACERhere$corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data

# remove object in the end to clear up some space
rm(ACERhere)
