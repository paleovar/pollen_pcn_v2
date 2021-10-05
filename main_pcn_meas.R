# Wrapper script for computing 
# - network measures of the networks generated from ACER proxy data and pseudo proxies from
#   TraCE21ka, HadCM3, LOVECLIM

USE_REPO_DATA_PCN_MEAS <- TRUE

FILE_PCN_MEAS_CACHE <- file.path(DIR_CACHE,'pre_processed_pcn_meas_data.RData')


# script flag preventing script to be executed if this has been done in the session before
# overwrite STATUS_PCNS_MEAS_RUN <- FALSE in global environment to force re-execution
if (!exists('STATUS_PCNS_MEAS_RUN')) STATUS_PCNS_MEAS_RUN <- FALSE

if (STATUS_PCNS_MEAS_RUN == FALSE) {
  message('main_pcn_meas.R has not been run in session, executing script with the defined flags')
} else {
  message('main_pcn_meas.R has been run in session, aborting executing the script. To force execution, set STATUS_PCNS_MEAS_RUN <- FALSE in global environment')
  invokeRestart('abort')
}


# if pre-processed simulation data should be used
if (USE_REPO_DATA_PCN_MEAS) {
  if (file.exists(FILE_PCN_MEAS_CACHE)) {
    message('Using pre-processed PCN measures, therefore skipping computing them from scratch')
    load(FILE_PCN_MEAS_CACHE)
    STATUS_PCNS_MEAS_RUN <- TRUE
    ## this loads main data objects into global environment:
    ## lds_trc_raw, lds_trc_mil, lds_lc_raw, lds_lc_mil, lds_hcm_raw ~ link densities
    ## ACERtrc_clp_raw, ACERtrc_clp_mil, ACERlc_clp_raw, ACERlc_clp_mil, ACERhcm_clp_raw ~ cross-link fractions
    ## ACERtrc_ndg_raw, ACERtrc_ndg_mil, ACERlc_ndg_raw, ACERlc_ndg_mil, ACERhcm_ndg_raw ~ node degrees (site-wise)
    invokeRestart('abort')
  } else {
    message('USE_REPO_DATA_PCN_MEAS == TRUE but no matching .RData file found, attempting to compute PCN measures from scratch')
  }
}

## Initialize all simulation and proxy data & PCN data (per default, skipped automatically in `main_simulation_data.R` & `main_pcns.R` if this has been done before) ----
## also skipped if network measure data has been loaded from .RData.
## (Note that skipping is implicit by `invokeRestart('abort')` above)
source('main_pcns.R')


## Compute/load PCN measures ----

#REGIONS_EXT <- REGIONS %>% 
#  bind_rows(
#    tibble(site_id = c(10, 54, 91, 58, 74, 27,  7, 55, 42,  6, 90, 98, 90.98, 53, 47, 57, 59,  2, 39, 82, 41, 24, 96),
#           region = c(rep('NAm-W',20),rep('NAm-E',3)))
#  ) %>% 
#  bind_rows(
#    tibble(site_id = c( 30,  43,  62,  75,  60,  44,  45,  80,  85, 100,  85.1,  83,  84,  93,  49,  61,  81,  63,  69,  68,  66,  70),
#           region = c('SAm-T/S', rep('SAm-M',3),rep('SAm-T/S',3),rep('SAm-M',4),rep('SAm-T/S',2),rep('SAm-M',2),rep('SAm-T/S',2),rep('SAm-M',5)))
#  ) %>% 
#  bind_rows(
#    tibble(site_id = c( 17, 77, 12, 23, 21, 20, 38, 67, 71, 97, 72, 87, 72.97),
#           region = c('Afr-CW', rep('Afr-CSE',5),rep('Afr-CW',2),rep('Afr-CSE',5)))
#  )


### TraCE ----

# check if network measures have to be computed or if they are available from disk
#FILES_ACERTRACE_MEAS <- list(
#  'ACERtrc_ld_raw' = file.path(DIR_FIGURES, 'paperMIS1-2', 'ACERtrc_ld_updnorm_raw.csv'),
#  'ACERtrc_ld_mil' = file.path(DIR_FIGURES, 'paperMIS1-2', 'ACERtrc_ld_updnorm_mil.csv'),
#  'ACERtrc_clp_raw' = file.path(DIR_FIGURES, 'paperMIS1-2', 'ACERtrc_clp_updnorm_correctedNA_woanom_raw.csv'),
#  'ACERtrc_clp_mil' = file.path(DIR_FIGURES, 'paperMIS1-2', 'ACERtrc_clp_updnorm_correctedNA_woanom_mil.csv'),
#  'ACERtrc_ndg_raw' = file.path(DIR_FIGURES, 'paperMIS1-2', 'ACERtrc_ndg_updnorm_woanom_raw.csv'),
#  'ACERtrc_ndg_mil' = file.path(DIR_FIGURES, 'paperMIS1-2', 'ACERtrc_ndg_updnorm_woanom_mil.csv')
#)

ACERtrc_nw_meas <- function(ACERtrc,ACERtrcmil) {
  # LDs
  lds_trc_raw <- lapply(1:length(ACERtrc), function(i){
    nm <- names(ACERtrc)[i]
    print(nm)
    if (nm == 'ACER_ap') {arb <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data'} else {arb <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data'}
    out <- compute_nw_link_dens_pn(ACERtrc[[i]][[arb]]) %>% 
      mutate(data = 'TRACE', signal = nm, tscale = 'raw')
    out
  }) %>% 
    bind_rows()
  
  lds_trc_mil <- lapply(1:length(ACERtrcmil), function(i){
    nm <- names(ACERtrcmil)[i]
    print(nm)
    if (nm == 'ACER_ap') {arb <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_detrend_gaussbp2000-8000_probit_arb_pollen_data'} else {arb <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_detrend_gaussbp2000-8000_identity_arb_pollen_data'}
    out <- ACERtrcmil[[i]][[arb]] %>% 
      filter(from <= 100 & to <= 100) %>% # filter NGRIP, EPICA
      compute_nw_link_dens_pn(.) %>% 
      mutate(data = 'TRACE', signal = nm, tscale = 'mil')
    out
  }) %>% 
    bind_rows()
  
  # CLFs
  ACERtrc_clp_raw <- lapply(1:length(ACERtrc), function(i) {
    nms <- names(ACERtrc)
    if (nms[i] == 'ACER_ap') {
      dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data'
    } else {
      dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data'
    }
    out <- ACERtrc[[i]][[dat]] %>%  
      compute_nw_crosslink_prob_pn(., save_to_file = list(activate = F), regions = REGIONS) %>% 
      rename(!!sym(paste0('pcross_link.p.', nms[i])) := pclp,
             !!sym(paste0('pcross_link.n.', nms[i])) := pcln)
    return(out)
  }) %>% 
    plyr::join_all(., by = c('window', 'cross_reg1', 'cross_reg2'), type = 'inner') %>% 
    as_tibble()
  
  ACERtrc_clp_mil <- lapply(1:length(ACERtrcmil), function(i) {
    nms <- names(ACERtrcmil)
    if (nms[i] == 'ACER_ap') {
      dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_detrend_gaussbp2000-8000_probit_arb_pollen_data'
    } else {
      dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_detrend_gaussbp2000-8000_identity_arb_pollen_data'
    }
    out <- ACERtrcmil[[i]][[dat]] %>% 
      filter(from <= 100 & to <= 100) %>% # filter NGRIP, EPICA
      compute_nw_crosslink_prob_pn(., save_to_file = list(activate = F), regions = REGIONS) %>% 
      rename(!!sym(paste0('pcross_link.p.', nms[i])) := pclp,
             !!sym(paste0('pcross_link.n.', nms[i])) := pcln)
    return(out)
  }) %>% 
    plyr::join_all(., by = c('window', 'cross_reg1', 'cross_reg2'), type = 'inner') %>% 
    as_tibble()
  
  # NDGs
  ACERtrc_ndg_raw <- lapply(1:length(ACERtrc), function(i) {
    nms <- names(ACERtrc)
    if (nms[i] == 'ACER_ap') {
      dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data'
    } else {
      dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data'
    }
    out <- ACERtrc[[i]][[dat]] %>%  
      compute_nw_node_degree_pn(., save_to_file = list(activate = F)) %>% 
      rename(!!sym(paste0('node_degree.p.', nms[i])) := ndgp,
             !!sym(paste0('node_degree.n.', nms[i])) := ndgn,
             !!sym(paste0('node_degree.', nms[i])) := ndg)
    return(out)
  }) %>% 
    plyr::join_all(., by = c('window', 'site_id'), type = 'full') %>% 
    as_tibble()
  
  ACERtrc_ndg_mil <- lapply(1:length(ACERtrcmil), function(i) {
    nms <- names(ACERtrcmil)
    if (nms[i] == 'ACER_ap') {
      dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_detrend_gaussbp2000-8000_probit_arb_pollen_data'
    } else {
      dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_detrend_gaussbp2000-8000_identity_arb_pollen_data'
    }
    out <- ACERtrcmil[[i]][[dat]] %>% 
      filter(from <= 100 & to <= 100) %>% # filter NGRIP, EPICA
      compute_nw_node_degree_pn(., save_to_file = list(activate = F)) %>% 
      rename(!!sym(paste0('node_degree.p.', nms[i])) := ndgp,
             !!sym(paste0('node_degree.n.', nms[i])) := ndgn,
             !!sym(paste0('node_degree.', nms[i])) := ndg)
    return(out)
  }) %>% 
    plyr::join_all(., by = c('window', 'site_id'), type = 'full') %>% 
    as_tibble()
  
  return(list('lds_trc_raw'=lds_trc_raw,'ACERtrc_clp_raw'=ACERtrc_clp_raw,'ACERtrc_ndg_raw'=ACERtrc_ndg_raw,
              'lds_trc_mil'=lds_trc_mil,'ACERtrc_clp_mil'=ACERtrc_clp_mil,'ACERtrc_ndg_mil'=ACERtrc_ndg_mil))
}

ACERtrc_nw_meas_dat <- ACERtrc_nw_meas(ACERtrc,ACERtrcmil)
lds_trc_raw <- ACERtrc_nw_meas_dat[['lds_trc_raw']]
lds_trc_mil <- ACERtrc_nw_meas_dat[['lds_trc_mil']]
ACERtrc_clp_raw <- ACERtrc_nw_meas_dat[['ACERtrc_clp_raw']]
ACERtrc_clp_mil <- ACERtrc_nw_meas_dat[['ACERtrc_clp_mil']]
ACERtrc_ndg_raw <- ACERtrc_nw_meas_dat[['ACERtrc_ndg_raw']]
ACERtrc_ndg_mil <- ACERtrc_nw_meas_dat[['ACERtrc_ndg_mil']]
rm(ACERtrc_nw_meas_dat)

#rm(ACERtrc,ACERtrcmil)


### LOVECLIM ----

# check if network measures have to be computed or if they are available from disk
#FILES_ACERLC_MEAS <- list(
#  'ACERlc_ld_raw' = file.path(DIR_FIGURES, 'paperMIS1-2', 'ACERlc_ld_updnorm_raw.csv'),
#  'ACERlc_ld_mil' = file.path(DIR_FIGURES, 'paperMIS1-2', 'ACERlc_ld_updnorm_mil.csv'),
#  'ACERlc_clp_raw' = file.path(DIR_FIGURES, 'paperMIS1-2', 'ACERlc_clp_updnorm_correctedNA_woanom_raw.csv'),
#  'ACERlc_clp_mil' = file.path(DIR_FIGURES, 'paperMIS1-2', 'ACERlc_clp_updnorm_correctedNA_woanom_mil.csv'),
#  'ACERlc_ndg_raw' = file.path(DIR_FIGURES, 'paperMIS1-2', 'ACERlc_ndg_updnorm_woanom_raw.csv'),
#  'ACERlc_ndg_mil' = file.path(DIR_FIGURES, 'paperMIS1-2', 'ACERlc_ndg_updnorm_woanom_mil.csv')
#)

ACERlc_nw_meas <- function(ACERlc,ACERlcmil) {
  lds_lc_raw <- lapply(1:length(ACERlc), function(i){
    nm <- names(ACERlc)[i]
    print(nm)
    if (nm == 'ACER_ap') {arb <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data'} else {arb <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data'}
    out <- compute_nw_link_dens_pn(ACERlc[[i]][[arb]]) %>% 
      mutate(data = 'LCLIM', signal = nm, tscale = 'raw')
    out
  }) %>% 
    bind_rows()
  
  lds_lc_mil <- lapply(1:length(ACERlcmil), function(i){
    nm <- names(ACERlcmil)[i]
    print(nm)
    if (nm == 'ACER_ap') {arb <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_detrend_gaussbp2000-8000_probit_arb_pollen_data'} else {arb <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_detrend_gaussbp2000-8000_identity_arb_pollen_data'}
    out <- ACERlcmil[[i]][[arb]] %>% 
      filter(from <= 100 & to <= 100) %>% # filter NGRIP, EPICA
      compute_nw_link_dens_pn(.) %>% 
      mutate(data = 'LCLIM', signal = nm, tscale = 'mil')
    out
  }) %>% 
    bind_rows()
  
  ACERlc_clp_raw <- lapply(1:length(ACERlc), function(i) {
    nms <- names(ACERlc)
    if (nms[i] == 'ACER_ap') {
      dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data'
    } else {
      dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data'
    }
    out <- ACERlc[[i]][[dat]] %>%  
      compute_nw_crosslink_prob_pn(., save_to_file = list(activate = F)) %>% 
      rename(!!sym(paste0('pcross_link.p.', nms[i])) := pclp,
             !!sym(paste0('pcross_link.n.', nms[i])) := pcln)
    return(out)
  }) %>% 
    plyr::join_all(., by = c('window', 'cross_reg1', 'cross_reg2'), type = 'inner') %>% 
    as_tibble()
  
  ACERlc_clp_mil <- lapply(1:length(ACERlcmil), function(i) {
    nms <- names(ACERlcmil)
    if (nms[i] %in% c('ACER_ap')) {#, 'LC_vc')) {
     dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_detrend_gaussbp2000-8000_probit_arb_pollen_data'
    } else {
     dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_detrend_gaussbp2000-8000_identity_arb_pollen_data'
    }
    out <- ACERlcmil[[i]][[dat]] %>%  
     compute_nw_crosslink_prob_pn(., save_to_file = list(activate = F)) %>% 
      rename(!!sym(paste0('pcross_link.p.', nms[i])) := pclp,
             !!sym(paste0('pcross_link.n.', nms[i])) := pcln)
    return(out)
  }) %>% 
    plyr::join_all(., by = c('window', 'cross_reg1', 'cross_reg2'), type = 'inner') %>% 
    as_tibble()
  
  ACERlc_ndg_raw <- lapply(1:length(ACERlc), function(i) {
    nms <- names(ACERlc)
    if (nms[i] == 'ACER_ap') {
      dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data'
    } else {
      dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data'
    }
    out <- ACERlc[[i]][[dat]] %>%  
      compute_nw_node_degree_pn(., save_to_file = list(activate = F)) %>% 
      rename(!!sym(paste0('node_degree.p.', nms[i])) := ndgp,
             !!sym(paste0('node_degree.n.', nms[i])) := ndgn,
             !!sym(paste0('node_degree.', nms[i])) := ndg)
    return(out)
  }) %>% 
    plyr::join_all(., by = c('window', 'site_id'), type = 'full') %>% 
    as_tibble()
  
  ACERlc_ndg_mil <- lapply(1:length(ACERlcmil), function(i) {
    nms <- names(ACERlcmil)
    if (nms[i] %in% c('ACER_ap')) {#, 'LC_vc')) {
      dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_detrend_gaussbp2000-8000_probit_arb_pollen_data'
    } else {
      dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_detrend_gaussbp2000-8000_identity_arb_pollen_data'
    }
    out <- ACERlcmil[[i]][[dat]] %>%  
      compute_nw_node_degree_pn(., save_to_file = list(activate = F)) %>% 
      rename(!!sym(paste0('node_degree.p.', nms[i])) := ndgp,
             !!sym(paste0('node_degree.n.', nms[i])) := ndgn,
             !!sym(paste0('node_degree.', nms[i])) := ndg)
    return(out)
  }) %>% 
    plyr::join_all(., by = c('window', 'site_id'), type = 'full') %>% 
    as_tibble() 
  return(list('lds_lc_raw'=lds_lc_raw,'ACERlc_clp_raw'=ACERlc_clp_raw,'ACERlc_ndg_raw'=ACERlc_ndg_raw,
              'lds_lc_mil'=lds_lc_mil,'ACERlc_clp_mil'=ACERlc_clp_mil,'ACERlc_ndg_mil'=ACERlc_ndg_mil))
}

ACERlc_nw_meas_dat <- ACERlc_nw_meas(ACERlc,ACERlcmil)
lds_lc_raw <- ACERlc_nw_meas_dat[['lds_lc_raw']]
lds_lc_mil <- ACERlc_nw_meas_dat[['lds_lc_mil']]
ACERlc_clp_raw <- ACERlc_nw_meas_dat[['ACERlc_clp_raw']]
ACERlc_clp_mil <- ACERlc_nw_meas_dat[['ACERlc_clp_mil']]
ACERlc_ndg_raw <- ACERlc_nw_meas_dat[['ACERlc_ndg_raw']]
ACERlc_ndg_mil <- ACERlc_nw_meas_dat[['ACERlc_ndg_mil']]
rm(ACERlc_nw_meas_dat)

#rm(ACERlc,ACERlcmil)


### HadCM3 ----

# check if network measures have to be computed or if they are available from disk
#FILES_ACERHCM_MEAS <- list(
#  'ACERhcm_ld_raw' = file.path(DIR_FIGURES, 'paperMIS1-2', 'ACERhcm_ld_updnorm_raw.csv'),
#  'ACERhcm_clp_raw' = file.path(DIR_FIGURES, 'paperMIS1-2', 'ACERhcm_clp_updnorm_correctedNA_woanom_raw.csv'),
#  'ACERhcm_ndg_raw' = file.path(DIR_FIGURES, 'paperMIS1-2', 'ACERhcm_ndg_updnorm_woanom_raw.csv')
#)

ACERhcm_nw_meas <- function(ACERhcm) {
  lds_hcm_raw <- lapply(1:length(ACERhcm), function(i){
    nm <- names(ACERhcm)[i]
    print(nm)
    if (nm == 'ACER_ap') {arb <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data'} else {arb <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data'}
    out <- compute_nw_link_dens_pn(ACERhcm[[i]][[arb]]) %>% 
     mutate(data = 'BBC', signal = nm, tscale = 'raw')
    out
  }) %>% 
    bind_rows()
  
  ACERhcm_clp_raw <- lapply(1:length(ACERhcm), function(i) {
    nms <- names(ACERhcm)
    if (nms[i] == 'ACER_ap') {
      dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data'
    } else {
      dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data'
    }
    out <- ACERhcm[[i]][[dat]] %>%  
      compute_nw_crosslink_prob_pn(., save_to_file = list(activate = F)) %>% 
      rename(!!sym(paste0('pcross_link.p.', nms[i])) := pclp,
             !!sym(paste0('pcross_link.n.', nms[i])) := pcln)
    return(out)
  }) %>% 
    plyr::join_all(., by = c('window', 'cross_reg1', 'cross_reg2'), type = 'inner') %>% 
    as_tibble()
  
  ACERhcm_ndg_raw <- lapply(1:length(ACERhcm), function(i) {
    nms <- names(ACERhcm)
    if (nms[i] == 'ACER_ap') {
      dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data'
    } else {
      dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data'
    }
    out <- ACERhcm[[i]][[dat]] %>%  
      compute_nw_node_degree_pn(., save_to_file = list(activate = F)) %>% 
      rename(!!sym(paste0('node_degree.p.', nms[i])) := ndgp,
             !!sym(paste0('node_degree.n.', nms[i])) := ndgn,
             !!sym(paste0('node_degree.', nms[i])) := ndg)
    return(out)
  }) %>% 
    plyr::join_all(., by = c('window', 'site_id'), type = 'full') %>% 
    as_tibble()
  
  return(list('lds_hcm_raw'=lds_hcm_raw,'ACERhcm_clp_raw'=ACERhcm_clp_raw,'ACERhcm_ndg_raw'=ACERhcm_ndg_raw))
}

ACERhcm_nw_meas_dat <- ACERhcm_nw_meas(ACERhcm)
lds_hcm_raw <- ACERhcm_nw_meas_dat[['lds_hcm_raw']]
ACERhcm_clp_raw <- ACERhcm_nw_meas_dat[['ACERhcm_clp_raw']]
ACERhcm_ndg_raw <- ACERhcm_nw_meas_dat[['ACERhcm_ndg_raw']]
rm(ACERhcm_nw_meas_dat)

#rm(ACERhcm)


#### end network measures

# cache network measures if permitted ----
if (USE_REPO_DATA_PCN_MEAS) {
  save(
    lds_trc_raw, lds_trc_mil, lds_lc_raw, lds_lc_mil, lds_hcm_raw,
    ACERtrc_clp_raw, ACERtrc_clp_mil, ACERlc_clp_raw, ACERlc_clp_mil, ACERhcm_clp_raw,
    ACERtrc_ndg_raw, ACERtrc_ndg_mil, ACERlc_ndg_raw, ACERlc_ndg_mil, ACERhcm_ndg_raw,
    file = FILE_PCN_MEAS_CACHE
  )
}

# garbage collect
gc()

# modifying run status of this script
STATUS_PCNS_MEAS_RUN <- TRUE
