# helpers for the null model (contained in util_networkres.R in v1)


# hic_ratio = 0 -> only epica, hic_ratio = 1 -> only greenland
snhem_hybr_sign <- function(min_age_ky = 0, max_age_ky = 115) {
  if(exists('epica_ngrip_cores')) {cores <- epica_ngrip_cores}
  else {
    if(exists('epica_raw')) {epica <- epica_raw}else {epica <- load_data_icecore('epica'); epica_raw <<- epica}
    if(exists('ngrip_raw')) {ngrip <- ngrip_raw}else {ngrip <- load_data_icecore('ngrip'); ngrip_raw <<- ngrip}
    epica <- epica %>% mutate(core = 'epica', d18O_corr2 = d18O_corr - mean(d18O_corr)) %>% filter(age <= max_age_ky & age >= min_age_ky) %>% mutate(diff = lead(age) - age)
    ngrip <- ngrip %>% mutate(core = 'ngrip', d18O_corr2 = d18O_corr - mean(d18O_corr)) %>% filter(age <= max_age_ky & age >= min_age_ky) %>% mutate(diff = lead(age) - age)
    
    dt <- max(mean(epica$diff, na.rm = TRUE), mean(ngrip$diff, na.rm = TRUE))
    tst <- seq(min_age_ky, max_age_ky, by = dt)
    
    cores <- tibble(age = tst,
                    epica = MakeEquidistant(t.x = epica$age, t.y = epica$d18O_corr2, time.target = tst), 
                    ngrip = MakeEquidistant(t.x = ngrip$age, t.y = ngrip$d18O_corr2, time.target = tst)) %>% 
      gather(key = 'core', value = 'd18O', -age) %>% filter(!is.na(d18O)) %>% 
      spread(core, d18O) %>% 
      filter(!is.na(epica) & !is.na(ngrip))
    epica_ngrip_cores <<- cores
  }
  return(cores)
}

cat('pre-loading ice core data')
invisible(snhem_hybr_sign())


snhem_hybr_sign_stage <- function(min_age = 0, max_age = 115000, tout = NULL, hic_ratio = list(mode = 'lat', mode_args = list(eqt_ratio = 0.5, npole_lat = 90, spole_lat = -90)), lat) {
  if (!(hic_ratio$mode %in% c('lat', 'fixed'))) {stop('unknown hic_ratio mode given to sample_hic_sign')}
  if (hic_ratio$mode == 'lat') {if (lat > hic_ratio$mode_args$npole_lat) {hic_ratio_loc <- 1}
    else if (lat < hic_ratio$mode_args$spole_lat) {hic_ratio_loc <- 0}
    else {hic_ratio_loc <- lat / (hic_ratio$mode_args$npole_lat - hic_ratio$mode_args$spole_lat) + 0.5}}
  else {hic_ratio_loc <- hic_ratio$mode_args$eqt_ratio}
  
  cores <- snhem_hybr_sign(min_age_ky = min_age / 1000, max_age_ky = max_age / 1000) %>% 
    mutate(age = age * 1000) %>% 
    mutate(hic.hres = (1 - hic_ratio_loc) * epica + hic_ratio_loc * ngrip) %>% 
    mutate(hic.hres = (hic.hres - mean(hic.hres, na.rm = TRUE)) / sd(hic.hres, na.rm = TRUE)) %>% 
    #mutate(hic.hres = hic.hres / sum(abs(hic.hres), na.rm = TRUE)) %>% 
    select(age, hic.hres) %>% 
    filter(!is.na(hic.hres))
  return(cores)
}


# sample hybrid ice core for one site of given latitude and with time steps
sample_hic_sign <- function(data = NULL, sites = NULL, params = NULL, lat = NULL) {
  if (!is.null(params) & !is.null(sites) & !is.null(data)) {
    data <- do(data, sample_hic_sign_site_(tout = .$mixed_age, lat = unique(.$lat), sampling_params = params$sampling_params))
  }
  else if (!is.null(lat)) {
    data <- sample_hic_sign_site_(lat = lat)
  }
  return(data)
}


# ext call = TRUE if all distributions should be returned instead of only the final signal
sample_hic_sign_withnoise <- function(data, sites, params, ext_call = FALSE) {
  def_noise_params <- list(snvr = 9)
  if ('noise_params' %in% names(params)) {noise_params <- params$noise_params}else {noise_params <- list()}
  noise_params <- merge.list(noise_params, def_noise_params)
  if ('hic_noise' %in% names(params)) {
    hic_noise <- params$hic_noise
  } else if (params$noise_params$temporal == 'ar1' & params$noise_params$spatial == 'matern') {
    hic_noise <- ar1matern_noise_for_ACERsites(sites = sites, spatial_params = params$noise_params$spatial_params, temporal_params = params$noise_params$temporal_params)
  } else if (params$noise_params$temporal == 'ar1' & params$noise_params$spatial == 'white') {
    hic_noise <- ar1white_noise_for_ACERsites(sites = sites, spatial_params = params$noise_params$spatial_params, temporal_params = params$noise_params$temporal_params)
  } else {
    stop('unknown noise params given to sample_hic_sign_withnoise - known: temporal: ar1, spatial: matern, white')
  }
  
  hic_data <- select(sites, site_id, lat) %>% 
    group_by(site_id) %>% 
    do(hic_data = snhem_hybr_sign_stage(lat = unique(.$lat))) %>% 
    ungroup() %>% 
    inner_join(hic_noise, by = 'site_id') %>% 
    group_by(site_id) %>% 
    mutate(nhic_data = map_dbl(hic_data, nrow), nnoise_data = map_dbl(noise_data, nrow)) %>% 
    mutate(n_match = if_else(nhic_data == nnoise_data, TRUE, FALSE))
  
  if (!all(hic_data$n_match)) {stop('lengths of hic_data and corresponding noise do not match, check NA in snhem_hybr_sign -> check age limits -> rm(epica_ngrip_cores)')}
  hic_data <- hic_data %>% 
    select(-n_match, -nhic_data, -nnoise_data) %>% 
    mutate(hicnoise_data = purrr::map2(hic_data, noise_data, san, snvr = noise_params$snvr))
  
  if(ext_call) {return(hic_data)}
  else {
    data <- group_by(data, site_id) %>% 
      nest(.key = 'sample_data') %>% 
      inner_join(select(hic_data, site_id, hicnoise_data), by = 'site_id') %>% 
      mutate(hicnoise_data_sampled = purrr::map2(sample_data, hicnoise_data, sample_hic_sign_withnoise_site_, sampling_params = params$sampling_params)) %>% 
      select(site_id, hicnoise_data_sampled) %>% 
      unnest()
    return(data)
  }
}


# mix signal and noise tibbles
# hic signal and noise need to have zero mean and unite variance (and have to be normalised ?)
# snvr (signal to noise variance ratio) -> unite variance maintained
san <- function(s, n, snvr) {
  sn <- inner_join(s, n, by = 'age') %>% 
    mutate(hic.hres = sqrt(snvr / (1 + snvr)) * hic.hres + sqrt(1 / (1 + snvr)) * hic_noise) %>% 
    select(-hic_noise)
  return(sn)
}


sample_hic_sign_site_ <- function(min_age = 0, max_age = 115000, tout = NULL, hic_ratio = list(mode = 'lat', mode_args = list(eqt_ratio = 0.5, npole_lat = 90, spole_lat = -90)),
                                  lat, sampling_params = list(sampling_mode = 'box', polamp = 4)) {
  def_sampling_params <- list(sampling_mode = 'box', polamp = 4)
  sampling_params <- modifyList(def_sampling_params, sampling_params)
  if (is.null(tout)) {warning('no tout given to sample_hic_sign using default resolution'); sampling_params$sampling_mode <- 'default'}
  if (!(hic_ratio$mode %in% c('lat', 'fixed'))) {stop('unknown hic_ratio mode given to sample_hic_sign')}
  if (!(sampling_params$sampling_mode %in% c('interpolation', 'box', 'point', 'default'))) {stop('unknown sampling_mode given to sample_hic_sign_site_')}
  
  if (hic_ratio$mode == 'lat') {if (lat > hic_ratio$mode_args$npole_lat) {hic_ratio_loc <- 1}
    else if (lat < hic_ratio$mode_args$spole_lat) {hic_ratio_loc <- 0}
    else {hic_ratio_loc <- lat / (hic_ratio$mode_args$npole_lat - hic_ratio$mode_args$spole_lat) + 0.5}}
  else {hic_ratio_loc <- hic_ratio$mode_args$eqt_ratio}
  
  hic_sign <- snhem_hybr_sign(min_age_ky = min_age / 1000, max_age_ky = max_age / 1000) %>% 
    mutate(age = age * 1000) %>% 
    mutate(hic.hres = (1 - hic_ratio_loc) * epica + hic_ratio_loc * ngrip) %>% 
    mutate(hic.hres = hic.hres / sum(abs(hic.hres), na.rm = TRUE)) %>% 
    mutate(hic.hres = hic.hres / sd(hic.hres, na.rm = TRUE)) %>% 
    select(age, hic.hres) %>% 
    mutate(hic.hres = ((abs(lat) - 90) / 90 * (1 / sampling_params$polamp) + 1) * hic.hres) %>% 
    {if (sampling_params$sampling_mode == 'interpolation') {interp.dataset(y = ., x = .$age, xout = tout, method = 'linear', rep.negt = FALSE) %>% as_tibble(.)}
      else if (sampling_params$sampling_mode == 'box') {sample.dataset.box(y = ., tout = tout)}
      else if (sampling_params$sampling_mode == 'point') {sample.dataset.point(y = ., tout = tout)}
      else {.}} %>% 
    rename(hic = hic.hres) %>% 
    drop_na()
  
  return(hic_sign)
}


sample_hic_sign_withnoise_site_ <- function(sample_data, hicnoise_data, sampling_params = list(sampling_mode = 'box', polamp = 4), ...) {
  def_sampling_params <- list(sampling_mode = 'box', polamp = 4)
  sampling_params <- merge.list(sampling_params, def_sampling_params)
  if (!(sampling_params$sampling_mode %in% c('interpolation', 'box', 'point', 'default'))) {stop('unknown sampling_mode given to sample_hic_sign_withnoise_site_')}
  lat <- sample_data$lat %>% unique(.)
  tout <- arrange(sample_data, mixed_age) %>% .$mixed_age
  
  hic_sign <- hicnoise_data %>% 
    mutate(hic.hres = ((abs(lat) - 90) / 90 * (1 / sampling_params$polamp) + 1) * hic.hres) %>% 
    {if (sampling_params$sampling_mode == 'interpolation') {interp.dataset(y = ., x = .$age, xout = tout, method = 'linear', rep.negt = FALSE) %>% as_tibble(.)}
      else if (sampling_params$sampling_mode == 'box') {sample.dataset.box(y = ., tout = tout)}
      else if (sampling_params$sampling_mode == 'point') {sample.dataset.point(y = ., tout = tout)}
      else {.}} %>% 
    rename(hic = hic.hres) %>% 
    drop_na()
  
  return(hic_sign)
}


make_single_set_hybrid_icecores_for_ACERsites_ <- function(use_cached_noise, hic_noise_iscached, n, sites, samples, withnoise, write_noise, file_test_noise, params, run = 1, ncores = NULL) {
  export <- c('make_hybrid_icecores_for_ACERsites', 'sample_hic_sign_withnoise', 'sample_hic_sign', 'snhem_hybr_sign_stage', 'snhem_hybr_sign',
              'load_data_icecore', 'ar1matern_noise_for_ACERsites', 'san', 'sample_hic_sign_withnoise_site_', 'sample.dataset.box', 'sample.dataset.point', 
              'make.bins.aggregate', 'which.min.vec', 'bin.ids')
  packages <- c('mvnfast', 'fields', 'dplyr', 'RCurl', 'readr', 'PaleoSpec', 'tibble', 'tidyr', 'stringr', 'purrr')
  
  if (!use_cached_noise) {
    cat(cr('calculating new noise sample\n'))
    args <- list(sites = sites, samples = samples, withnoise = withnoise, params = params)
    cl <- detectCores(); if (!is.null(ncores)) {if (ncores > cl) {stop('ncores provided to make_single_set_hybrid_icecores_for_ACERsites_ is larger than detected number of cores')} else {cl <- ncores}}
    registerDoParallel(cl)
    print(cl)
    data <- foreach (i = (1 + (run - 1) * n):(n * run), .combine = 'list', .inorder = FALSE, .export = export, .packages = packages) %dopar% {do.call(make_hybrid_icecores_for_ACERsites, args = args)}
    stopImplicitCluster()
  }
  else {
    if (!hic_noise_iscached) {
      cat(cr('no noise of given length cached, calculating new noise sample\n'))
      #hic_noise_cached <- make_set_hic_noise(n = n, sites = sites, matern_params = params$matern_params, ar1_params = params$ar1_params) %>% 
      hic_noise_cached <- make_set_hic_noise(n = n, sites = sites, params = params, ncores = ncores) %>% 
        mutate(noise_sample_id = noise_sample_id + ((run - 1) * n))
    }
    else {cat(cr('reading noise from cached file\n')); hic_noise_cached <- read_set_noise_hybrid_icecores_for_ACERsites(mode = 'single_slice', dir = file_test_noise)}
    cl <- detectCores(); if (!is.null(ncores)) {if (ncores > cl) {stop('ncores provided to make_single_set_hybrid_icecores_for_ACERsites_ is larger than detected number of cores')} else {cl <- ncores}}
    registerDoParallel(cl)
    data <- foreach (i = (1 + (run - 1) * n):(n * run), .combine = 'list', .inorder = FALSE, .export = export, .packages = packages) %dopar% {
      hic_noise_loc <- filter(hic_noise_cached, noise_sample_id == i) %>% select(-noise_sample_id)
      args <- list(sites = sites, samples = samples, withnoise = withnoise,
                   params = merge.list(params, list(hic_noise = hic_noise_loc)))
      do.call(make_hybrid_icecores_for_ACERsites, args = args)}
    stopImplicitCluster()
  }
  
  renquote <- function(l) if (class(l)[1] == 'list') lapply(l, renquote) else enquote(l)
  
  if (all(class(data) == 'list')) {
    data <- lapply(unlist(renquote(data)), eval)
    data <- bind_rows(data, .id = 'noise_sample_id') %>% 
      mutate_at(vars(noise_sample_id), as.numeric) %>% 
      mutate(noise_sample_id = noise_sample_id + ((run - 1) * n))
  } else {
    data <- mutate(data, noise_sample_id = ((run - 1) * n))
  }
  
  if (use_cached_noise) {
    return(list(trial = data, hic_noise_cached = hic_noise_cached))
  }
  else {
    return(data)
  }
}


# main call
make_hybrid_icecores_for_ACERsites <- function(sites, samples, withnoise = FALSE, params) {
  # for correlated noise in space and time: sample one climate field: fun(sites, range, theta, smoothness) -> AR(1) process per site on EP/NG common high res time scale: fun(climate_field, core_time_scale, ret:nested series)
  # -> nest noise time series by site id: sub() -> join to nested version of ACER site data (nesting of mixed_age by site_id and lat) : sub() -> sample hic_sign with error fct dealing with the two nested series: fun() analog to sample_hic_sign
  # -> within: add/multiply noise to high res individual hic core -> within: call sampling stage -> return as unnested, sampled time series to main call
  calls <- list('TRUE' = sample_hic_sign_withnoise, 'FALSE' = sample_hic_sign)
  data <- inner_join(select(sites, site_id, lat), select(samples, site_id, sample_id, mixed_age), by = 'site_id') %>% 
    group_by(site_id) %>% 
    arrange(mixed_age, .by_group = TRUE) %>% 
    calls[[as.character(withnoise)]](data = ., sites = sites, params = params) %>% 
    rename(mixed_age = age) %>% 
    rownames_to_column(var = 'hic_id') %>% 
    inner_join(select(samples, mixed_age, sample_id, site_id), by = c('site_id', 'mixed_age')) %>% 
    select(hic_id, sample_id, site_id, hic)
  data$hic_id <- as.numeric(data$hic_id)
  return(data)
}


# main call (set)
# cached trials are always recycled automatically, cached noise is only recycled if permitted to do so
make_set_hybrid_icecores_for_ACERsites <- function(n, sites, samples, withnoise = TRUE,
                                                   write_trial = list(activate = TRUE, mode = 'single', dir = paste(DIR_NM_DATA, 'hic', sep = '/')),
                                                   write_noise = list(activate = TRUE, mode = 'single', dir = paste(DIR_NM_DATA, 'noise', sep = '/')),
                                                   use_cached_noise = FALSE,
                                                   params = list(noise_params = list(snvr = 9,
                                                                                     temporal = 'ar1',
                                                                                     temporal_params = list(alpha = 0.8),
                                                                                     spatial = 'matern',
                                                                                     spatial_params = list(theta = 1000, smoothness = 1.5)),
                                                                 sampling_params = list(sampling_mode = 'box',
                                                                                        polamp = 4)), 
                                                   ncores = NULL) {
  #params = list(matern_params = list(theta = 1000, smoothness = 1.5),
  #              ar1_params = list(alpha = 0.8),
  #              noise_params = list(snvr = 9),
  #              sampling_params = list(sampling_mode = 'box', polamp = 4))) {
  # if more modes implemented: add mode name to vector in any
  if (write_trial$activate & !any(str_detect(write_trial$mode, c('single', 'slice')))) {stop('unknown write mode given to make_set_hybrid_icecores_for_ACERsites')}
  write_trial <- modifyList(list(mode = 'single', dir = paste(DIR_NM_DATA, 'hic', sep = '/')), write_trial)
  if (write_noise$activate & !any(str_detect(write_noise$mode, c('single', 'slice')))) {stop('unknown write mode given to make_set_hybrid_icecores_for_ACERsites')}
  write_noise <- modifyList(list(mode = 'single', dir = paste(DIR_NM_DATA, 'noise', sep = '/')), write_noise)
  
  #def... <- list(matern_params = list(theta = 1000, smoothness = 1.5), ar1_params = list(alpha = 0.8), noise_params = list(snvr = 9), sampling_params = list(sampling_mode = 'box', polamp = 4))
  if (params$noise_params$spatial == 'matern') {
    def... <- list(noise_params = list(snvr = 9, temporal = 'ar1', temporal_params = list(alpha = 0.8), spatial = 'matern', spatial_params = list(theta = 1000, smoothness = 1.5)), sampling_params = list(sampling_mode = 'box', polamp = 4))
  } else if (params$noise_params$spatial == 'white') {
    def... <- list(noise_params = list(snvr = 9, temporal = 'ar1', temporal_params = list(alpha = 0.8), spatial = 'white', spatial_params = list()), sampling_params = list(sampling_mode = 'box', polamp = 4))
  } else {
    stop('unknown spatial noise mode given to make_set_hybrid_icecores_for_ACERsites; known: matern, white')
  }
  if (!is.null(params) & any(!(names(params) %in% names(def...)))) {stop('unsupported parameter given to make_set_hybrid_icecores_for_ACERsites')}
  params <- modifyList(def..., params)
  
  cat('using params:\n')
  print(hash::hash(modifyList(list(n = n, nwrite_mode = write_trial$mode), params)))
  
  if (write_trial$mode == 'single' & write_noise$mode == 'single') {
    cat(cr('set_hybrid_icecores in mode single\n'))
    # check whether trial file is cached
    dir_test_trial <- paste(write_trial$dir, paste(as.character(unlist(as.list(hash::hash(modifyList(list(n = n, nwrite_mode = write_trial$mode), params))))), collapse = '-'), sep = '/')
    file_test_trial <- paste(paste(dir_test_trial, paste(as.character(unlist(as.list(hash::hash(modifyList(list(n = n, nwrite_mode = write_trial$mode), params))))), collapse = '-'), sep = '/'), DATAFILES_TYPE, sep = '.')
    if (dir.exists(dir_test_trial) & file.exists(file_test_trial)){file_cached <- TRUE; read_mode <- write_trial$mode}
    #else if (any(str_detect(list.dirs(write$dir), paste0(n, '-slice')))) {file_cached <- TRUE; read_mode <- paste(n, '-slice')}
    else {file_cached <- FALSE}
    
    if (!file_cached) {
      cat(cr('trial not yet cached, calculating noisy hybrid ice cores\n'))
      # check eventually if noise of given length is cached
      if (use_cached_noise) {
        noise_params_test <- params$noise_params; noise_params_test$snvr <- NULL
        dir_test_noise <- paste(write_noise$dir, paste('hic_noise', paste(as.character(unlist(as.list(hash::hash(list(n = n, nwrite_mode = write_noise$mode, noise_params = noise_params_test))))), collapse = '-'), sep = '_'), sep = '/')
        file_test_noise <- paste(paste(dir_test_noise, paste('hic_noise', paste(as.character(unlist(as.list(hash::hash(list(n = n, nwrite_mode = write_noise$mode, noise_params = noise_params_test))))), collapse = '-'), sep = '_'), sep = '/'), DATAFILES_TYPE, sep = '.')
        if (dir.exists(dir_test_noise) & file.exists(file_test_noise)) {hic_noise_iscached <- TRUE} 
        else {hic_noise_iscached <- FALSE}
      }
      else {hic_noise_iscached <- FALSE}
      
      data <- make_single_set_hybrid_icecores_for_ACERsites_(use_cached_noise = use_cached_noise, hic_noise_iscached = hic_noise_iscached, n = n,
                                                             sites = sites, samples = samples, withnoise = withnoise, write_noise = write_noise, file_test_noise = file_test_noise,
                                                             params = params, ncores = ncores)
      if (use_cached_noise) {
        hic_noise_cached <- data[['hic_noise_cached']]
        data <- data[['trial']]
      }
      
      if (write_trial$activate) {
        cat(cr('writing trial to disk\n'))
        name <- paste(as.character(unlist(as.list(hash::hash(modifyList(list(n = n, nwrite_mode = write_trial$mode), params))))), collapse = '-')
        dir <- paste(write_trial$dir, name, sep = '/')
        if (!dir.exists(dir)) {dir.create(dir)}
        write_csv(x = data, path = paste(paste(dir, name, sep = '/'), DATAFILES_TYPE, sep = '.'), col_names = TRUE)
      }
      if (write_noise$activate & !hic_noise_iscached & use_cached_noise) {
        cat(cr('writing noise to disk\n'))
        hic_noise_cached <- unnest(hic_noise_cached) %>% 
          inner_join(select(snhem_hybr_sign(), age) %>% arrange(age) %>% rownames_to_column(var = 'age_id') %>%
                       mutate_at(vars(age_id), as.numeric) %>% mutate(age = age * 1000),
                     by = 'age') %>% 
          select(noise_sample_id, site_id, age_id, hic_noise)
        if (!dir.exists(dir_test_noise)) {dir.create(dir_test_noise, recursive = TRUE)}
        write_csv(x = hic_noise_cached, path = file_test_noise, col_names = TRUE)
      }
    }
    else {
      cat(cr('reading trial from cached file\n'))
      data <- read_set_hybrid_icecores_for_ACERsites(mode = read_mode, dir = file_test_trial)
      addsites <- dplyr::union(data$site_id, sites$site_id ) %>% dplyr::setdiff(sites$site_id)
      if (length(addsites) != 0) {
        
      }
    }
  }
  else if (str_detect(write_trial$mode, 'slice') & str_detect(write_noise$mode, 'slice')) {
    slice <- as.integer(str_replace(write_trial$mode, 'slice', ''))
    if (slice != as.integer(str_replace(write_noise$mode, 'slice', ''))) {stop('write_noise$mode has to equal write_trial$mode')}
    if (n %% slice != 0) {stop('n %% slice != 0 in make_set_hybrid_icecores_for_ACERsites')}
    else {slices <- n / slice}
    data <- list()
    for (r in 1:slices) {
      cat(cr(paste0('slice ', r, '\n')))
      
      # check whether trial file is cached
      sites <- distinct(sites) %>% arrange(site_id)
      hsh <- paste0(paste0(sites$site_id, collapse = '.'), paste0(sites$site_name, collapse = '.')) %>%  
        digest::digest(., serialize = F, algo = 'crc32')
      dir_test_trial <- paste(as.character(unlist(as.list(hash::hash(modifyList(list(n = n, nwrite_mode = write_trial$mode), params))))), collapse = '-') %>% 
        paste0('-prx', hsh)
      file_test_trial <- paste(as.character(unlist(as.list(hash::hash(modifyList(list(n = n, nwrite_mode = paste0(write_trial$mode, 'run', r)), params))))), collapse = '-') #%>% 
      #paste0('-prx', hsh)
      dir_test_trial <- paste(write_trial$dir, dir_test_trial, sep = '/')
      file_test_trial <- paste(paste(dir_test_trial, file_test_trial, sep = '/'), DATAFILES_TYPE, sep = '.')
      if (dir.exists(dir_test_trial) & file.exists(file_test_trial)){file_cached <- TRUE; read_mode <- 'single_slice'}
      else {file_cached <- FALSE}
      
      if (!file_cached) {
        cat(cr('trial not yet cached, calculating noisy hybrid ice cores\n'))
        # check eventually if noise of given length is cached
        if (use_cached_noise) {
          noise_params_test <- params$noise_params; noise_params_test$snvr <- NULL
          dir_test_noise <- paste('hic_noise', paste(as.character(unlist(as.list(hash::hash(list(n = n, nwrite_mode = write_noise$mode, noise_params = noise_params_test))))), collapse = '-'), sep = '_') %>% 
            paste0('-prx', hsh)
          file_test_noise <- paste('hic_noise', paste(as.character(unlist(as.list(hash::hash(list(n = n, nwrite_mode = paste0(write_noise$mode, 'run', r), noise_params = noise_params_test))))), collapse = '-'), sep = '_') #%>% 
          #paste0('-prx', hsh)
          dir_test_noise <- paste(write_noise$dir, dir_test_noise, sep = '/')
          file_test_noise <- paste(paste(dir_test_noise, file_test_noise, sep = '/'), DATAFILES_TYPE, sep = '.')
          if (dir.exists(dir_test_noise) & file.exists(file_test_noise)) {hic_noise_iscached <- TRUE} 
          else {hic_noise_iscached <- FALSE}
        }
        else {hic_noise_iscached <- FALSE}
        
        trial <- make_single_set_hybrid_icecores_for_ACERsites_(use_cached_noise = use_cached_noise, hic_noise_iscached = hic_noise_iscached, n = slice,
                                                                sites = sites, samples = samples, withnoise = withnoise, write_noise = write_noise,
                                                                file_test_noise = file_test_noise, params = params, run = r, ncores = ncores)
        if (use_cached_noise) {
          hic_noise_cached <- trial[['hic_noise_cached']]
          trial <- trial[['trial']]
        }
        
        if (write_trial$activate) {
          cat(cr('writing trial to disk\n'))
          #name <- paste(as.character(unlist(as.list(hash::hash(modifyList(list(n = n, nwrite_mode = write_trial$mode), params))))), collapse = '-')
          #dir <- paste(write_trial$dir, name, sep = '/')
          if (!dir.exists(dir_test_trial)) {dir.create(dir_test_trial)}
          write_csv(x = trial, path = file_test_trial, col_names = TRUE)
        }
        if (write_noise$activate & !hic_noise_iscached & use_cached_noise) {
          cat(cr('writing noise to disk\n'))
          hic_noise_cached <- unnest(hic_noise_cached) %>% 
            inner_join(select(snhem_hybr_sign(), age) %>% arrange(age) %>% rownames_to_column(var = 'age_id') %>%
                         mutate_at(vars(age_id), as.numeric) %>% mutate(age = age * 1000),
                       by = 'age') %>% 
            select(noise_sample_id, site_id, age_id, hic_noise)
          if (!dir.exists(dir_test_noise)) {dir.create(dir_test_noise, recursive = TRUE)}
          write_csv(x = hic_noise_cached, path = file_test_noise, col_names = TRUE)
        }
        data[[as.character(r)]] <- trial
      }
      else {
        cat(cr('reading trial from cached file\n'))
        data[[as.character(r)]] <- read_set_hybrid_icecores_for_ACERsites(mode = read_mode, dir = file_test_trial)
      }
    }
  }
  else {
    stop('unknown write mode')
    #stop('not used code part')
    #slice <- as.integer(str_replace(write_trial$mode, 'slice', ''))
    #if (n %% slice != 0) {stop('n %% slice != 0 in make_set_hybrid_icecores_for_ACERsites')}
    ##data <- vector(mode = 'list', length = slice)
    ##tasks <- seq(1, n, 1)
    ##for (i in 1:n) {
    ##  data[[i]] <- do.call(make_hybrid_icecores_for_ACERsites, args = args)
    ##  if (i %% slice == 0) {data <- bind_rows(data, .id = 'noise_sample_id')
    ##    write_csv(x = data, path = paste(paste(dir, paste(name, paste0('slice', i), sep = '_'), sep = '/'), DATAFILES_TYPE, sep = '.'), col_names = TRUE)
    ##    data <- vector(mode = 'list', length = slice)
    ##  }
    ##}
    ##data <- mclapply(tasks, function(x, args) {do.call(make_hybrid_icecores_for_ACERsites, args = args)}, args = args)
    #cl <- detectCores()
    #registerDoParallel(cl)
    #data <- foreach (i = 1:n, .combine = 'c', .inorder = FALSE) %:% 
    #  foreach(j = 1:slice, .combine = 'c', .inorder = FALSE) %do% {do.call(make_hybrid_icecores_for_ACERsites, args = args)}
    ##loc_data <- bind_rows(loc_data, .id = 'noise_sample_id1')
    ##write_csv(x = loc_data, path = paste(paste(dir, paste(name, paste0('slice', i), sep = '_'), sep = '/'), DATAFILES_TYPE, sep = '.'), col_names = TRUE)
    ##loc_data
    #stopCluster(cl)
    #
    #data <- bind_rows(data, .id = 'noise_sample_id') #2') %>% 
    ##group_by(noise_sample_id1, noise_sample_id2) %>% 
    ##nest() %>% 
    ##rownames_to_column(var = 'noise_sample_id') %>% 
    ##ungroup() %>% 
    ##select(noise_sample_id, data) %>% 
    ##unnest()
  }
  return(data)
}


read_set_hybrid_icecores_for_ACERsites <- function(mode, dir) {
  if (mode == 'single' | mode == 'single_slice') {data <- read_csv(file = dir, col_names = TRUE, col_types = list(col_double(), col_double(), col_double(), col_double(), col_double()))}
  else {stop('unimplemented read mode in read_set_gybrid_icecores')}
  return(data)
}


read_set_noise_hybrid_icecores_for_ACERsites <- function(mode, dir) {
  if (mode == 'single' | mode == 'single_slice') {
    data <- read_csv(file = dir, col_names = TRUE, col_types = list(col_double(), col_double(), col_double(), col_double())) %>% 
      inner_join(select(snhem_hybr_sign(), age) %>% arrange(age) %>% rownames_to_column(var = 'age_id') %>% 
                   mutate_at(vars(age_id), as.numeric) %>% mutate(age = age * 1000),
                 by = 'age_id') %>% 
      select(noise_sample_id, site_id, age, hic_noise) %>% 
      group_by(noise_sample_id, site_id) %>% 
      nest(.key = 'noise_data')
  }
  else {stop('unimplemented read mode in read_set_gybrid_icecores')}
  return(data)
}


sample.dataset.box <- function(y, tout, signal = 'hic.hres') {
  y <- arrange(y, age) %>% 
    rownames_to_column(var = 'id') %>% 
    mutate(id = as.integer(id)) %>% 
    make.bins.aggregate(data = ., tout = tout, signal = signal)
  #stop()
  #print(y)
  return(y)
}


sample.dataset.point <- function(y, tout) {
  y <- rownames_to_column(y, var = 'id') %>% 
    mutate(id = as.integer(id)) %>% 
    slice(which.min.vec(tin = .$age, tout = tout)) %>% 
    bind_cols(tibble(tout_age = tout)) %>% 
    select(-age) %>% 
    rename(age = tout_age)
  return(y)
}


# find closest data points (in time) in hybrid icecore for given time series of a record
which.min.vec <- function(tin, tout) {
  mins <- vector()
  i <- 1
  for (t in tout) {
    mins[[i]] <- which.min(abs(tin - t))
    i <- i + 1
  }
  return(mins)
}


# bin hybrid ice core data and box sample to given time series record of a record
make.bins.aggregate <- function(data, tout, signal = 'hic.hres') {
  tout <- tibble(tout = tout) %>% filter(tout > 0 & tout < max(data$age)) %>% .$tout
  orig_data <- data
  id_data <- slice(data, which.min.vec(tin = data$age, tout = tout)) %>% 
    select(id) %>% 
    arrange(id) %>% 
    mutate(id_diff_ahead = lead(id) - id, id_diff_aft = id - lag(id)) %>% 
    mutate(id_diff_aft = if_else(is.na(id_diff_aft), id_diff_ahead, id_diff_aft), id_diff_ahead = if_else(is.na(id_diff_ahead), id_diff_aft, id_diff_ahead)) %>%
    mutate(id_aft = if_else(id - floor(id_diff_aft / 2) < 1, 1, id - floor(id_diff_aft / 2)), id_cen = id, id_ahead = id + floor(id_diff_ahead / 2)) %>% # find some way to deal with overlap!!
    select(-id_diff_ahead, -id_diff_aft, -id) %>% 
    gather(key = 'id_kind', value = 'id')
  #id_data$id_kind <- ordered(id_data$id_kind, levels = c('id_aft', 'id_cen', 'id_ahead'))
  data <- left_join(arrange(data, id), id_data, by = 'id')
  data$bin <- bin.ids(data$id_kind)
  data_id <- data
  data <- filter(data, !is.na(bin)) %>% 
    group_by(bin) %>% 
    summarise(!!sym(signal) := mean(!!sym(signal)), .groups = 'drop')
  
  if ((length(data$bin) != length(tout)) & (nrow(data) != 0)) {
    stop('data might be resolved higher than input data, check make.bins.aggregate')
    #warning('data might be resolved higher than input data, returning point sampled values instead of box sampled values')
    #prob[[prob_count]] <<- list(data_id, id_data, data, tout); prob_count <<- prob_count + 1
    #data <- orig_data %>% 
    #  slice(which.min.vec(tin = .$age, tout = tout)) %>% 
    #  bind_cols(tibble(tout_age = tout)) %>% 
    #  select(-age) %>% 
    #  rename(age = tout_age)
  }
  else if (nrow(data) > 0) {
    data <- bind_cols(data, tibble(age = tout))
  } else {
    data <- tibble(bin = numeric(0), !!sym(signal) := numeric(0))
  }
  
  return(data)
}

make_set_hic_noise <- function(n = n, sites = sites, params, ncores = NULL) {
  if (!params$noise_params$temporal %in% c('ar1') | !params$noise_params$spatial %in% c('matern', 'white')) {stop('unknown noise params given to make_set_hic_noise - known: temporal: ar1, spatial: matern, white')}
  cl <- detectCores(); if (!is.null(ncores)) {if (ncores > cl) {stop('ncores provided to make_single_set_hybrid_icecores_for_ACERsites_ is larger than detected number of cores')} else {cl <- ncores}}
  registerDoParallel(cl)
  export <- c('ar1matern_noise_for_ACERsites', 'ar1white_noise_for_ACERsites', 'snhem_hybr_sign', 'load_data_icecore')
  packages <- c('dplyr', 'mvnfast', 'tibble', 'tidyr', 'fields', 'RCurl', 'readr', 'PaleoSpec')
  args <- list(sites = sites, spatial_params = params$noise_params$spatial_params, temporal_params = params$noise_params$temporal_params)
  #args <- list(sites = sites, matern_params = matern_params, ar1_params = ar1_params)
  
  if (params$noise_params$temporal == 'ar1' & params$noise_params$spatial == 'matern') {
    data <- foreach (i = 1:n, .combine = 'list', .inorder = FALSE, .export = export, .packages = packages) %dopar% {
      j <- i
      do.call(ar1matern_noise_for_ACERsites, args = args)
    }
    stopImplicitCluster()
  } else if (params$noise_params$temporal == 'ar1' & params$noise_params$spatial == 'white') {
    data <- foreach (i = 1:n, .combine = 'list', .inorder = FALSE, .export = export, .packages = packages) %dopar% {
      j <- i
      do.call(ar1white_noise_for_ACERsites, args = args)
    }
    stopImplicitCluster()
  } else {
    stop('unknown noise params given to make_set_hic_noise - known: temporal: ar1, spatial: matern, white')
  }
  
  renquote <- function(l) if (class(l)[1] == 'list') lapply(l, renquote) else enquote(l)
  if (all(class(data) == 'list')) {
    data <- lapply(unlist(renquote(data)), eval)
    data <- bind_rows(data, .id = 'noise_sample_id') %>% 
      mutate_at(vars(noise_sample_id), as.numeric)
  } else {
    data <- mutate(data, noise_sample_id = 1)
  }
  
  return(data)
}


# helper for binning
bin.ids <- function(data) {
  group <- 0
  status <- FALSE
  new_data <- vector()
  for (i in seq(1, length(data), 1)) {
    if(!is.na(data[[i]]) & !str_detect(data[[i]], 'id_cen')) {if (str_detect(data[[i]], 'id_aft') == TRUE) {group <- group + 1; status <- TRUE; new_data[[i]] <- group}
      else if(str_detect(data[[i]], 'id_ahead') == TRUE) {status <- FALSE; new_data[[i]] <- group}}
    else if (status) {new_data[[i]] <- group}else {new_data[[i]] <- NA}
  }
  return(new_data)
}


# helper for retrieving window list from label list
window.char.to.num <- function(label_list) {
  window_list <- vector()
  for (i in label_list) {
    w <- paste0(as.vector(unlist(str_split(i, '-'))), '000') %>% as.numeric(.)
    window_list <- c(window_list, w)
  }
  return(as.list(sort(unique(window_list))))
}


# noising & noise sampling ------------
# noise routines are separated on higher level to accelerate parallelisation

# matern covariance routines (adapted from N. Weitzel) ---------
ar1matern_noise_for_ACERsites <- function(sites, spatial_params = list(theta = 1000, smoothness = 1.5), temporal_params = list(alpha = 0.8), 
                                          tout_wrapper = 'snhem_hybr_sign') {
  def_matern_params <- list(theta = 1000, smoothness = 1.5); def_ar1_params <- list(alpha = 0.8)
  spatial_params <- merge.list(spatial_params, def_matern_params); temporal_params <- merge.list(temporal_params, def_ar1_params)
  
  sites <- arrange(sites, site_id)
  tout <-  do.call(tout_wrapper, args = list())$age# snhem_hybr_sign()$age
  dim <- sites$site_id %>% length(.)
  acer_coord <- array(c(sites$long, sites$lat), dim = c(dim, 2))
  # Matern covariance "climate field"
  matern_cov <- stationary.cov(acer_coord, Covariance="Matern", Distance="rdist.earth", theta = spatial_params$theta, smoothness = spatial_params$smoothness)
  #Sample records from field (initial rd sample from field, then AR1 using Cholesky decomposition)
  time_steps <- length(tout)
  alpha <- temporal_params$alpha
  sim_ar1 <- array(0, dim = c(time_steps, dim))
  #sim_ar1[1,] <- as.numeric(rmvnorm(1, sigma = matern_cov))
  mu <- rep(0, dim)
  sim_ar1[1,] <- as.numeric(rmvn(1, mu = mu, sigma = matern_cov, kpnames = TRUE))
  for (i in 2:time_steps) {
    sim_ar1[i,] <- alpha * sim_ar1[i-1,] + as.numeric(rmvn(1, mu = mu, sigma = matern_cov, kpnames = TRUE)) #as.numeric(rmvnorm(1, sigma = matern_cov))
  }
  
  sim_ar1 <- as.tibble(sim_ar1)
  colnames(sim_ar1) <- sites$site_id
  sim_ar1$age <- tout
  sim_ar1 <- gather(sim_ar1, key = 'site_id', value = 'hic_noise', -age) %>% 
    select(site_id, age, hic_noise) %>% 
    group_by(site_id) %>% 
    mutate(age = age * 1000, hic_noise = (hic_noise - mean(hic_noise))/ sd(hic_noise)) %>% 
    #mutate(hic_noise = hic_noise / sum(abs(hic_noise))) %>% 
    nest(.key = 'noise_data') %>% 
    ungroup() %>% 
    mutate_at(vars(site_id), as.numeric)
  return(sim_ar1)
}
