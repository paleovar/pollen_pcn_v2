# Wrapper script for testing ACER AP PCN sensitivity to
# - marine vs terrestrial sites
# - square root transform instead of probit transformation (generally only done for ACER AP)
# Computations of PCNs are done in to `main_pcns.R`.
# Computations of network measures are analogous to `main_pcn_meas.R`.
# Used in `create_fig_sensitivity.R`

USE_REPO_SENS <- TRUE

FILE_SENS_CACHE <- file.path(DIR_CACHE,'pre_processed_sensitivity_data.RData')

# script flag preventing script to be executed if this has been done in the session before
# overwrite STATUS_SENS_RUN <- FALSE in global environment to force re-execution
if (!exists('STATUS_SENS_RUN')) STATUS_SENS_RUN <- FALSE

if (STATUS_SENS_RUN == FALSE) {
  message('main_pcn_sensitivity.R has not been run in session, executing script with the defined flags')
} else {
  message('main_pcn_sensitivity.R has been run in session, aborting executing the script. To force execution, set STATUS_SENS_RUN <- FALSE in global environment')
  invokeRestart('abort')
}


# if pre-processed simulation data should be used
if (USE_REPO_SENS) {
  if (file.exists(FILE_SENS_CACHE)) {
    message('Using pre-processed PCN data, therefore skipping computing them from scratch')
    load(FILE_SENS_CACHE)
    STATUS_SENS_RUN <- TRUE
    ## this loads main data objects into global environment:
    ## `ACERsens_stype_err_sgn`, `ld_rank`, `ACERsens_clp_raw`, `ACERsens_ndg_raw`, `ACERtrf`, `ACERtrf_clp_raw`, `ACERtrf_ndg_raw`, `lds_trf_raw`
    invokeRestart('abort')
  } else {
    message('USE_REPO_SENS == TRUE but no matching .RData file found, attempting to compute PCNs from scratch')
  }
}


# load network data
source('main_pcns.R')
source('main_pcn_meas.R')

# sensitivity of ACER AP network to terrestrial vs marine sites
## filter data depending for sensitivity cases ----
ACERsens <- ACERtrc$ACER_ap

site_types <- c('TERR', 'MARI')
summary_name <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data'

ACERsens_stype <- lapply(site_types, function(t) {
  a <- ACERsens
  a$sites <- filter(a$sites, site_type == t)
  a[[summary_name]] <- a[[summary_name]] %>% 
    filter(from %in% a$sites$site_id & to %in% a$sites$site_id)
  # don't filter the other data here because they are not needed for LD analysis, contingency matrices, and PCN plots
  # just setting them to zero here to avoid confusion
  a$corr_list_gaussian_kernel_hres_probit_arb_pollen_data <- NULL
  a$arb_pollen_data <- NULL
  a$harm_pollen_data <- NULL
  return(a)
}) %>% 
  setNames(site_types)

## compute link densities ----
lds_stype <- lapply(1:length(ACERsens_stype), function(i){
  nm <- names(ACERsens_stype)[i]
  print(nm)
  out <- compute_nw_link_dens_pn(ACERsens_stype[[i]][[summary_name]]) %>% 
    mutate(data = 'ACER', signal = nm, tscale = 'raw')
  out
}) %>% 
  bind_rows()


## compare link density against randomly selected sets of links ----
n_sets <- 100

# marine sites contained in the ACER AP link set
sites_contained <- ACERsens[[summary_name]] %>% 
  select(from, to) %>% 
  gather(value = 'site_id') %>% 
  distinct(site_id) %>% 
  inner_join(ACERsens$sites %>% select(site_id,site_type))
n_site_type <- sites_contained %>% group_by(site_type) %>% summarise(n = n())

# random assignment of sites to surrogate categories "MARI" and "TERR"
site_sets <- lapply(1:n_sets, function(i) {
  a <- sample(1:(length(sites_contained$site_id)), n_site_type %>% filter(site_type == 'MARI') %>% .$n, replace = FALSE)
  set_m <- sites_contained %>% select(-site_type) %>% .[a,] %>% 
    mutate(site_type = 'MARI')
  set_t <- sites_contained %>% select(-site_type) %>% filter(!(site_id %in% set_m$site_id)) %>% 
    mutate(site_type = 'TERR')
  return(bind_rows(set_m,set_t))
})

# compute link densities for surrogate assignment to site types
ACERsens_stype_err <- lapply(1:n_sets, function(n) {
  ACERsens_lcl <- lapply(site_types, function(t) {
    a <- ACERsens
    a$sites <- a$sites %>% select(-site_type) %>% inner_join(site_sets[[n]], by = 'site_id') %>% filter(site_type == t)
    a[[summary_name]] <- a[[summary_name]] %>% 
      filter(from %in% a$sites$site_id & to %in% a$sites$site_id)
    a$corr_list_gaussian_kernel_hres_probit_arb_pollen_data <- NULL
    a$arb_pollen_data <- NULL
    a$harm_pollen_data <- NULL
    return(a)
  }) %>% 
    setNames(site_types)
  lds_lcl <- lapply(1:length(ACERsens_lcl), function(i){
    nm <- names(ACERsens_lcl)[i]
    out <- compute_nw_link_dens_pn(ACERsens_lcl[[i]][[summary_name]]) %>% 
      mutate(data = 'ACER', signal = nm, tscale = 'raw')
    out
  }) %>% 
    bind_rows() %>% 
    mutate(n_set = n)
  return(lds_lcl)
}) %>% 
  bind_rows()

# determine rank of actual link density compared to randomly assigned sites
ld_rank <- lapply(1:nrow(lds_stype), function(i) {
  x <- ACERsens_stype_err %>% select(link_dens, ld_sgn, signal, n_set) %>% filter(ld_sgn == lds_stype[i,]$ld_sgn, signal == lds_stype[i,]$signal) %>% 
    arrange(link_dens) %>% .$link_dens
  min(which(x > lds_stype[i,]$link_dens))
}) %>% unlist()
ld_rank <- lds_stype %>% select(link_dens, ld_sgn, signal) %>% 
  inner_join(tibble(ld_sgn = lds_stype$ld_sgn, signal = lds_stype$signal, 
                    ld_rank = ld_rank)) %>% 
  bind_rows(lds_trc_raw %>% filter(signal == 'ACER_ap') %>% select(link_dens,ld_sgn) %>% mutate(signal = 'ALL'))

ACERsens_stype_err_sgn <- ACERsens_stype_err %>% mutate(link_dens = if_else(ld_sgn == 'ldp', link_dens, -1*link_dens), ld_sgn = if_else(ld_sgn == 'ldp', '+', '-'))
ACERsens_stype_err_sgn$ld_sgn <- factor(ACERsens_stype_err_sgn$ld_sgn, levels = c('+','-'))
ld_rank <- ld_rank %>% mutate(link_dens = if_else(ld_sgn == 'ldp', link_dens, -1*link_dens), 
                              ld_sgn = if_else(ld_sgn == 'ldp', '+', '-')) %>% 
  mutate(link_dens_lab = if_else(ld_sgn == '+', link_dens+0.012, link_dens-0.012))
ld_rank$ld_sgn <- factor(ld_rank$ld_sgn, levels = c('+','-'))


## compare cross-link fractions and node degrees ----
ACERsens_clp_raw <- lapply(1:length(ACERsens_stype), function(i) {
  nms <- names(ACERsens_stype)
  out <- ACERsens_stype[[i]][[summary_name]] %>%  
    compute_nw_crosslink_prob_pn(., save_to_file = list(activate = F), regions = REGIONS) %>% 
    rename(!!sym(paste0('pcross_link.p.', nms[i])) := pclp,
           !!sym(paste0('pcross_link.n.', nms[i])) := pcln)
  return(out)
}) %>% 
  plyr::join_all(., by = c('window', 'cross_reg1', 'cross_reg2'), type = 'inner') %>% 
  as_tibble()

ACERsens_ndg_raw <- lapply(1:length(ACERsens_stype), function(i) {
  nms <- names(ACERsens_stype)
  out <- ACERsens_stype[[i]][[summary_name]] %>%  
    compute_nw_node_degree_pn(., save_to_file = list(activate = F)) %>% 
    rename(!!sym(paste0('node_degree.p.', nms[i])) := ndgp,
           !!sym(paste0('node_degree.n.', nms[i])) := ndgn,
           !!sym(paste0('node_degree.', nms[i])) := ndg)
  return(out)
}) %>% 
  plyr::join_all(., by = c('window', 'site_id'), type = 'full') %>% 
  as_tibble()

## some summarizing stats
ld_rank %>% group_by(signal) %>% summarise(sum = sum(abs(link_dens))) # total LD
ld_rank %>% select(-link_dens_lab,-ld_rank) %>% group_by(signal) %>% spread(ld_sgn, link_dens) %>% mutate(frac_n = abs(`-`)/(`+`-`-`)) # frac of pos & neg links


# sensititvity of ACER AP network to sqrt transform/no transformation
## prepare data ----
ACERhere <- PCNdata() %>% 
  mix_sample_dating_and_err() %>% 
  harmonized_pollen_to_db() %>% 
  arboreal_pollen_to_db()
ACERtrf <- list('identity' = ACERhere, 'sqrt' = ACERhere)
rm(ACERhere)

## compute networks analogous to `main_pcns.R` ----
# generate hash for record set
hsh <- gen_site_hash(ACERtrf[[1]]$sites, ACERtrf[[1]]$sample_dating)

# generate/load null model and compute/load correlations
invisible(
  lapply(1:length(ACERtrf), function(i) {
    # load/generate null models draws
    ACERtrf[[i]] <<- ACERtrf[[i]] %>% 
      set_hybrid_icecores_to_db(n = 1000, write_trial = list(activate = DO_CACHEING_PCN, mode = 'slice100'), write_noise = list(activate = DO_CACHEING_PCN, mode = 'slice100'),
                                use_cached_noise = TRUE, params = list(noise_params = list(snvr = 0,
                                                                                           temporal = 'ar1',
                                                                                           temporal_params = list(alpha = 0.8),
                                                                                           spatial = 'matern',
                                                                                           spatial_params = list(theta = 1000, smoothness = 1.5)),
                                                                       sampling_params = list(sampling_mode = 'box',
                                                                                              polamp = 4)))
    # compute distribution of similarity measure for the set of null model draws
    ACERtrf[[i]] <<- ACERtrf[[i]] %>% 
      corr_list_from_set_to_db(orig_data = paste0('hic_set_box-4-1000-', 0, '-ar1-0.8-matern-1000-1.5-slice100'),
                               type_data = 'hybrid_ice_core', 
                               transform = 'identity',
                               detrend = list(activate = FALSE),
                               hres_only = TRUE, interp = list(method = 'gaussian_kernel'), sd_one = FALSE,
                               write_corr = list(activate = DO_CACHEING_PCN, dir = file.path(DIR_NM_DATA, 'correlations', '6_22')),
                               windows = WINDOWS_TRACE, labels = WINDOW_LABELS_TRACE, ncores = NCORES_PCN) 
    
    # dump the null model draws because they just take up space
    # keep the similarity measurements
    ACERtrf[[i]]$`hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100` <<- NULL
  })
)

# compute/load correlations of simulation pseudo-proxies
# and test significance/load significance test results
invisible(
  lapply(1:length(ACERtrf), function(i) {
    # compute the similarity of the data
    ACERtrf[[i]] <<- ACERtrf[[i]] %>% 
      corr_list_to_db(orig_data = 'arb_pollen_data', 
                      type_data = 'arboreal_pollen',
                      transform = names(ACERtrf)[i],
                      detrend = list(activate = FALSE),
                      hres_only = TRUE,
                      interp = list(method = 'gaussian_kernel'),
                      sd_one = FALSE, windows = WINDOWS_TRACE, labels = WINDOW_LABELS_TRACE,
                      in_place = TRUE)
    # evaluate against the similarity measure distribution of the null model
    ACERtrf[[i]] <<- ACERtrf[[i]] %>% 
      corr_set_summary_to_db(orig_data = paste0('corr_list_gaussian_kernel_hres_', names(ACERtrf)[i],'_arb_pollen_data'), 
                             noise_data = 'corr_list_gaussian_kernel_hres_identity_hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100', test_option = 'ratio')
    # drop the similarity measurements of the null model draws
    # keep the statistics
    ACERtrf[[i]][[paste0('corr_list_gaussian_kernel_hres_', names(ACERtrf)[i],'_hic_set_box-4-1000-0-ar1-0.8-matern-1000-1.5-slice100')]] <<- NULL
  })
)

# compute network measures analogous to `main_pcn_meas.R` ----
## compute link densities and add ACER probit computed in `main_pcn_meas.R` ----
lds_trf_raw <- lapply(1:length(ACERtrf), function(i){
  nm <- names(ACERtrf)[i]
  print(nm)
  sum_name <- paste0('corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_', nm,'_arb_pollen_data')
  out <- compute_nw_link_dens_pn(ACERtrf[[i]][[sum_name]]) %>% 
    mutate(data = 'ACER', signal = nm, tscale = 'raw')
  out
}) %>% 
  bind_rows() %>% 
  bind_rows(lds_trc_raw %>% filter(signal == 'ACER_ap') %>% mutate(data = 'ACER', signal = 'probit'))

## compare cross-link fractions and node degrees ----
ACERtrf_clp_raw <- lapply(1:length(ACERtrf), function(i) {
  nm <- names(ACERtrf)[i]
  sum_name <- paste0('corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_', nm,'_arb_pollen_data')
  out <- ACERtrf[[i]][[sum_name]] %>%  
    compute_nw_crosslink_prob_pn(., save_to_file = list(activate = F), regions = REGIONS) %>% 
    rename(!!sym(paste0('pcross_link.p.', nm)) := pclp,
           !!sym(paste0('pcross_link.n.', nm)) := pcln)
  return(out)
}) %>% 
  plyr::join_all(., by = c('window', 'cross_reg1', 'cross_reg2'), type = 'inner') %>% 
  as_tibble()

ACERtrf_ndg_raw <- lapply(1:length(ACERtrf), function(i) {
  nm <- names(ACERtrf)[i]
  sum_name <- paste0('corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_', nm,'_arb_pollen_data')
  out <- ACERtrf[[i]][[sum_name]] %>%  
    compute_nw_node_degree_pn(., save_to_file = list(activate = F)) %>% 
    rename(!!sym(paste0('node_degree.p.', nm)) := ndgp,
           !!sym(paste0('node_degree.n.', nm)) := ndgn,
           !!sym(paste0('node_degree.', nm)) := ndg)
  return(out)
}) %>% 
  plyr::join_all(., by = c('window', 'site_id'), type = 'full') %>% 
  as_tibble()


# garbage collect for parallel processes/file loading
gc()

# cache data if permitted ----
if (USE_REPO_SENS) {
  save(
    ACERsens_stype_err_sgn, ld_rank, ACERsens_clp_raw, ACERsens_ndg_raw, ACERtrf, ACERtrf_clp_raw, ACERtrf_ndg_raw, lds_trf_raw,
    file = FILE_SENS_CACHE
  )
}

# modifying run status of this script
STATUS_SENS_RUN <- TRUE

