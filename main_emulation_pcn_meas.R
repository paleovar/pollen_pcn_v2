# Preliminary wrapper script to compute
# - link densities of these networks of paleoclimate networks from emulated TSF response of TraCE and HadCM3
# as preparation for `create_fig_emul.R` (preliminary name).


# Same flags for network measures (as in `main_pcn_meas.R`)
USE_REPO_DATA_EMUL_PCN_MEAS <- TRUE

FILE_EMUL_PCN_MEAS_CACHE <- file.path(DIR_CACHE,'pre_processed_emul_pcn_meas_data.RData')


# script flag preventing script to be executed if this has been done in the session before
# overwrite STATUS_EMUL_PCNS_MEAS_RUN <- FALSE in global environment to force re-execution
if (!exists('STATUS_EMUL_PCNS_MEAS_RUN')) STATUS_EMUL_PCNS_MEAS_RUN <- FALSE

if (STATUS_EMUL_PCNS_MEAS_RUN == FALSE) {
  message('main_emul_pcn_meas.R has not been run in session, executing script with the defined flags')
} else {
  message('main_emul_pcn_meas.R has been run in session, aborting executing the script. To force execution, set STATUS_EMUL_PCNS_MEAS_RUN <- FALSE in global environment')
  invokeRestart('abort')
}


# if pre-processed simulation data should be used
if (USE_REPO_DATA_EMUL_PCN_MEAS) {
  if (file.exists(FILE_EMUL_PCN_MEAS_CACHE)) {
    message('Using pre-processed emulated surrogates, therefore skipping computing them from scratch')
    load(FILE_EMUL_PCN_MEAS_CACHE)
    STATUS_EMUL_PCNS_MEAS_RUN <- TRUE
    ## this loads main data objects into global environment:
    ## lds_trc_em, lds_hcm_em ~ link densities
    ## ACERtrc_em_clp, ACERhcm_em_clp ~ cross-link fractions
    ## ACERtrc_em_ndg, ACERhcm_em_ndg ~ node degrees (site-wise)
    invokeRestart('abort')
  } else {
    message('USE_REPO_DATA_EMUL_PCN_MEAS == TRUE but no matching .RData file found, attempting to compute PCN measure data from scratch')
  }
}

## Initialize all simulation and proxy data & PCN data (per default, 
## skipped automatically in `main_emulation_data.R` & `main_emulation_pcns.R` if this has been done before) ----
## also skipped if network measure data has been loaded from .RData.
## (Note that skipping is implicit by `invokeRestart('abort')` above)
source('main_emulation_pcns.R')

## Compute/load PCN measures ----
### Emulated TraCE response ----
ACERtrc_em_nw_meas <- function(ACERtrc_em) {
  # LDs
  lds_trc_em <- lapply(1:length(ACERtrc_em), function(i){
    nm <- names(ACERtrc_em)[i]
    print(nm)
    if (nm == 'ACER_ap') {arb <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data'} else {arb <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data'}
    out <- compute_nw_link_dens_pn(ACERtrc_em[[i]][[arb]]) %>% 
      mutate(data = 'TRACE', signal = nm, tscale = 'raw_em')
    out
  }) %>% 
    bind_rows()
  
  # CLFs
  ## continents
  ACERtrc_em_clp <- lapply(1:length(ACERtrc_em), function(i) {
    nms <- names(ACERtrc_em)
    if (nms[i] == 'ACER_ap') {
      dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data'
    } else {
      dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data'
    }
    out <- ACERtrc_em[[i]][[dat]] %>%  
      compute_nw_crosslink_prob_pn(., save_to_file = list(activate = F), regions = REGIONS) %>% 
      rename(!!sym(paste0('pcross_link.p.', nms[i])) := pclp,
             !!sym(paste0('pcross_link.n.', nms[i])) := pcln)
    return(out)
  }) %>% 
    plyr::join_all(., by = c('window', 'cross_reg1', 'cross_reg2'), type = 'inner') %>% 
    as_tibble()
  
  # NDGs
  ACERtrc_em_ndg<- lapply(1:length(ACERtrc_em), function(i) {
    nms <- names(ACERtrc_em)
    if (nms[i] == 'ACER_ap') {
      dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data'
    } else {
      dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data'
    }
    out <- ACERtrc_em[[i]][[dat]] %>%  
      compute_nw_node_degree_pn(., save_to_file = list(activate = F)) %>% 
      rename(!!sym(paste0('node_degree.p.', nms[i])) := ndgp,
             !!sym(paste0('node_degree.n.', nms[i])) := ndgn,
             !!sym(paste0('node_degree.', nms[i])) := ndg)
    return(out)
  }) %>% 
    plyr::join_all(., by = c('window', 'site_id'), type = 'full') %>% 
    as_tibble()
  
  return(list('lds_trc_em'=lds_trc_em,'ACERtrc_em_clp'=ACERtrc_em_clp,'ACERtrc_em_ndg'=ACERtrc_em_ndg))
}

ACERtrc_em_nw_meas_dat <- ACERtrc_em_nw_meas(ACERtrc_em)
lds_trc_em <- ACERtrc_em_nw_meas_dat[['lds_trc_em']]
ACERtrc_em_clp <- ACERtrc_em_nw_meas_dat[['ACERtrc_em_clp']]
ACERtrc_em_ndg <- ACERtrc_em_nw_meas_dat[['ACERtrc_em_ndg']]
rm(ACERtrc_em_nw_meas_dat)


### emulated HadCM3 response ----
ACERhcm_em_nw_meas <- function(ACERhcm_em) {
  lds_hcm_em <- lapply(1:length(ACERhcm_em), function(i){
    nm <- names(ACERhcm_em)[i]
    print(nm)
    if (nm == 'ACER_ap') {arb <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data'} else {arb <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data'}
    out <- compute_nw_link_dens_pn(ACERhcm_em[[i]][[arb]]) %>% 
      mutate(data = 'BBC', signal = nm, tscale = 'raw_em')
    out
  }) %>% 
    bind_rows()
  
  ACERhcm_em_clp <- lapply(1:length(ACERhcm_em), function(i) {
    nms <- names(ACERhcm_em)
    if (nms[i] == 'ACER_ap') {
      dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data'
    } else {
      dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data'
    }
    out <- ACERhcm_em[[i]][[dat]] %>%  
      compute_nw_crosslink_prob_pn(., save_to_file = list(activate = F)) %>% 
      rename(!!sym(paste0('pcross_link.p.', nms[i])) := pclp,
             !!sym(paste0('pcross_link.n.', nms[i])) := pcln)
    return(out)
  }) %>% 
    plyr::join_all(., by = c('window', 'cross_reg1', 'cross_reg2'), type = 'inner') %>% 
    as_tibble()
  
  ACERhcm_em_ndg <- lapply(1:length(ACERhcm_em), function(i) {
    nms <- names(ACERhcm_em)
    if (nms[i] == 'ACER_ap') {
      dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data'
    } else {
      dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data'
    }
    out <- ACERhcm_em[[i]][[dat]] %>%  
      compute_nw_node_degree_pn(., save_to_file = list(activate = F)) %>% 
      rename(!!sym(paste0('node_degree.p.', nms[i])) := ndgp,
             !!sym(paste0('node_degree.n.', nms[i])) := ndgn,
             !!sym(paste0('node_degree.', nms[i])) := ndg)
    return(out)
  }) %>% 
    plyr::join_all(., by = c('window', 'site_id'), type = 'full') %>% 
    as_tibble()
  
  return(list('lds_hcm_em'=lds_hcm_em,'ACERhcm_em_clp'=ACERhcm_em_clp,'ACERhcm_em_ndg'=ACERhcm_em_ndg))
}

ACERhcm_em_nw_meas_dat <- ACERhcm_em_nw_meas(ACERhcm_em)
lds_hcm_em <- ACERhcm_em_nw_meas_dat[['lds_hcm_em']]
ACERhcm_em_clp <- ACERhcm_em_nw_meas_dat[['ACERhcm_em_clp']]
ACERhcm_em_ndg <- ACERhcm_em_nw_meas_dat[['ACERhcm_em_ndg']]
rm(ACERhcm_em_nw_meas_dat)

#### end network measures

# cache network measures if permitted ----
if (USE_REPO_DATA_EMUL_PCN_MEAS) {
  save(
    lds_trc_em, ACERtrc_em_clp, ACERtrc_em_ndg, 
    lds_hcm_em, ACERhcm_em_clp, ACERhcm_em_ndg,
    file = FILE_EMUL_PCN_MEAS_CACHE
  )
}

# garbage collect
gc()

# modifying run status of this script
STATUS_EMUL_PCNS_MEAS_RUN <- TRUE
