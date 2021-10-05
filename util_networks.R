# helpers to compute, handle, and evaluate PCNs

#' @rdname acer.networks
#' @param sites tibble of ACER sites
#' @param p overall probability for a link occurence
#'
#' @return tidygraph object
#' @export
create_network_ACER_random <- function(sites, p = 0.1) {
  sites <- distinct(sites, lat, long, .keep_all = TRUE) %>% 
    rownames_to_column(var = 'node_id')
  sites$node_id <- as.numeric(sites$node_id)
  
  graph <- as_tbl_graph(x = sample_gnp(n = count(sites)$n, p = p)) %>% 
    activate(nodes) %>% 
    mutate(node_id = seq(1, count(sites)$n, 1)) %>% 
    inner_join(sites, by = 'node_id')
  
  return(graph)
}


#' @rdname acer.networks
#' @param corr_list 
#' @param sites 
#' @param frac_abs_strongest is a two sided fraction!
#' @param region_tree logical, TRUE will return graph object for region hierarchy and corr_list with corresponding nodes; FALSE will return graph object based on corr_list (both can be used within \code{\link{plot_network_spatial}})
#' @param split_window logical, should the returned graph contain all windows occuring in corr_list at once or should a list of multiple graphs, which cover different time spans be returned?
#' @param split_regions_by_sign logical, should different refions be assigned for positive (incl zero) and negative (excl zero) links? Might enhance readability when plotted with plot_network_spatial()
#' @param get_node_occur logical, should the occurence of a node in each specific time window be returned as well
#' @param orig_corr_bundle logical. should the original corr_bundle be returned for a bundled network? should normally be enabled to retain original stats for network plotting (occurences and node degree)
#' @param filter_sites vector or named list. vector of site ids to filter for all time windows, list named list to filter in specific windows with window labels as in corr_list, e.g. list('0-25' = c(1,2))
#'
#' @return list, list(graph = tidygraph object or list of tidygraph objects, occurs = tibble of node occurences)
#' @export
#'
#' @examples
create_network_from_corr_list <- function(corr_list, sites, frac_abs_strongest = 1, region_tree = TRUE, split_window = TRUE, split_regions_by_sign = F,
                                          get_node_occur = FALSE, orig_corr_bundle = FALSE, filter_sites = NULL) {
  if (orig_corr_bundle) {
    orig_corr_bdl <- combo_site_node_id(corr_list = corr_list, sites = sites, region_tree = region_tree, regions = REGIONS, region_locs = REGION_LOCS) %>% 
      .$corr_list %>% {if ('diff_sig' %in% names(corr_list)) {filter(., diff_sig == TRUE)}else {.}}
  }
  if (!is.null(filter_sites)) {
    if (is.null(names(filter_sites))) {
      corr_list <- filter(corr_list, from %in% filter_sites & to %in% filter_sites)
    }
    else {
      if (!all(unique(sort(corr_list$window)) == sort(names(filter_sites)))) {stop('for named list as filter_sites: windows in list and corr_list do not match')}
      mxl <- lapply(filter_sites, function(x) length(x)) %>% unlist(., use.names = FALSE) %>% max(.)
      filter_sites <- filter_sites %>% lapply(., function(x) return(c(x, rep(NA_real_, mxl - length(x))))) %>% as_tibble() %>% 
        gather(key = 'window', value = 'site_id') %>% filter(!is.na(site_id)) %>% group_by(window) %>% nest(.key = 'filter')
      corr_list <- corr_list %>% group_by(window) %>% nest() %>% inner_join(filter_sites, by = 'window') %>% 
        mutate(data_filt = purrr::map2(data, filter, function(x, y) return(filter(x, from %in% y$site_id | to %in% y$site_id)))) %>% 
        unnest(data_filt)
    }
  }
  if (split_regions_by_sign) {
    combo <- combo_site_node_id(corr_list = corr_list, sites = sites, region_tree = region_tree)
  } else {
    combo <- combo_site_node_id(corr_list = corr_list, sites = sites, region_tree = region_tree, region_locs = REGION_LOCS)
  }
  sites <- combo$sites
  corr_list <- combo$corr_list %>% {if ('diff_sig' %in% names(corr_list)) {filter(., diff_sig == TRUE)}else {.}} #%>% select(-diff_sig)
  if (nrow(corr_list) == 0){warning('corr_list contains no significant values'); return(NA)}
  if (get_node_occur & orig_corr_bundle) {corrs <- orig_corr_bdl} else if (get_node_occur) {corrs <- corr_list}
  
  corr_list <- lapply(frac_abs_strongest, function(x, corr_list) return(mutate(corr_list, abs_corr = abs(corr)) %>% 
                                                                          group_by(window) %>% arrange(desc(abs_corr)) %>% 
                                                                          slice(1:ceiling(n() * x)) %>% ungroup() %>% 
                                                                          mutate(frac_abs_strongest = x)), corr_list = corr_list) %>% bind_rows(.)

  #if (length(frac_abs_strongest) == 1) {corr_list <- select(corr_list, -frac_abs_strongest)}
  if (region_tree) {
    if (split_window){stop('split windows supported only for ordinary non region tree')}
    else {
      graph <- list(graph = tbl_graph(nodes = sites, edges = combo$nw), corr_bundle = corr_list)
    }
  }
  else {
    if (split_window) {
      graph <- group_by(corr_list, window) %>% 
        nest(.key = 'corr_list') %>% 
        mutate(tbl_graph = purrr::map(corr_list, 
                                      function(x, sites) {x <- as_tbl_graph(x, directed = F) %>% activate(nodes) %>% rename(node_id = name) %>% mutate_at(vars(node_id), as.numeric) %>% full_join(sites, by = 'node_id') %>% arrange(node_id); return(x)},
                                      sites = sites))
    }
    else {
      graph <- corr_list %>% 
        as_tbl_graph(., directed = F) %>% activate(nodes) %>% rename(node_id = name) %>% mutate_at(vars(node_id), as.numeric) %>% full_join(sites, by = 'node_id') %>% arrange(node_id)
    }
    graph <- list(graph = graph)
  }
  
  if (get_node_occur) {
    windows <- corrs$window %>% unique()
    occurs <- bind_rows(corrs %>% group_by(window) %>% distinct(from, .keep_all = T) %>% select(from, window) %>% rename(node_id = from), 
                        corrs %>% group_by(window) %>% distinct(to, .keep_all = T) %>% select(to, window) %>% rename(node_id = to)) %>% 
      mutate(occurs = TRUE) %>% 
      distinct(window, node_id, .keep_all = T)
    occurs <- bind_rows(occurs,
                        {if (region_tree) {filter(sites, !hierarchial)}else {sites}} %>% select(node_id) %>%
                          do(tibble(node_id = rep(.$node_id, length(windows)), window = as.vector(sapply(windows, function(x, l) rep(x, l), l = length(.$node_id))))) %>% 
                          anti_join(occurs, by = c('node_id', 'window')) %>% 
                          mutate(occurs = FALSE)) %>% distinct() %>% ungroup()
    occurs$occurs <- factor(occurs$occurs, levels = c(TRUE, FALSE))
    graph[['occurs']] <- occurs
  }
  if (orig_corr_bundle) {
    graph$orig_corr_bundle <- orig_corr_bdl
  }
  return(graph)
}



create_bubble_network_from_corr_list <- function(corr_list, regions = REGIONS, region_locs = REGION_LOCS_BUBBLE,
                                                 clf = NULL) {
  
  # summarise by regions
  if (is.null(clf)) { # no clf data: mode nlinks
    rgs <- distinct(regions, region) %>% 
      add_column(reg_n = 1:length(unique(regions$region))) %>% 
      inner_join(regions, by = 'region') %>% 
      inner_join(region_locs, by = 'region')
    clist <- corr_list %>% 
      adj_list_from_corr_list_pn() %>% 
      select(from,to,corr,adj) %>% 
      inner_join(rgs %>% select(-region,-lat,-long) %>% rename(from = site_id, from_reg_n = reg_n), by = 'from') %>% 
      inner_join(rgs %>% select(-region,-lat,-long) %>% rename(to = site_id, to_reg_n = reg_n), by = 'to') %>% 
      select(-from,-to) %>% 
      rename(from = from_reg_n, to = to_reg_n) %>% 
      group_by(from,to,adj) %>% 
      summarise(nlinks = n(), .groups = 'drop') %>% 
      filter(adj != 0) # filter non-significant
    clist <- bind_rows(clist,
                       clist %>% rename(to = from, from = to)) %>% 
      group_by(from,to,adj) %>% 
      summarise(nlinks = sum(nlinks)) %>% 
      group_by(adj) %>% 
      filter(!duplicated(paste0(pmax(from, to), pmin(from, to))))
    clist$adj <- factor(clist$adj, levels = c(1,-1))
    
    rgs <- rgs %>% 
      filter(reg_n %in% c(clist$from,clist$to))
    
    graph <- tbl_graph(nodes = rgs %>% distinct(region,reg_n,lat,long) %>% arrange(reg_n) %>% rename(node_id = reg_n),
                       edges = clist)
    
  } else { # clf data provided: mode clf
    stop()
  }
  
  return(graph)
}




# correlation determination routines and helpers
# transform first to avoid insertion of new trends after alternative first detrending and to symmetrize distribution for detrending
# linear detrend in windows; linear detrend does not smooth variability as sampling distance is fairly large wrt window

# main calls for similarity measurement ----------
#' @rdname acer.networks
#' @param data 
#' @param corr_signal 
#' @param interp list of shape list(activate = logical, method = c('linear', 'gaussian_kernel')) to indicate interpolation method used for time series regularisation;
#' can easily be extended to additional interpolation methods from \code{\link[rioja]{rioja-package}}.
#' @param gaussbp not used at the moment; would be for additional gaussian filter
#' @param dopar logical, should parallelisation be enabled?
#'
#' @return a tibble (corr_list) with columns indicating window, site-site pair and corresponding correlation
#' @export 
#'
#' @examples
make_corr_list <- function(data, corr_signal, interp, gaussbp, dopar = TRUE) {
  # determine correlation in one window for each site combination
  data <- arrange(data, site_id)
  corr_list <- tibble(from = numeric(0), to = numeric(0), corr = numeric(0))
  
  dim <- length(data$site_id)
  if(!dopar) {
    for (i in seq(1, length(data$site_id) - 1, by = 1)) {
      for (j in seq(i, length(data$site_id), by = 1)) {
        
        x <- unnest(data[i, 2])
        y <- unnest(data[j, 2])
        
        xid <- distinct(x, site_id)$site_id
        yid <- distinct(y, site_id)$site_id
        
        corr <- time_series_correlation(select(x, site_id, mixed_age, !!corr_signal),
                                        select(y, site_id, mixed_age, !!corr_signal),
                                        corr_signal = corr_signal, 
                                        interp = interp, 
                                        gaussbp = gaussbp)
        
        corr_list <- add_row(corr_list, from = xid, to = yid, corr = corr[1])
      }
    }
  }
  else {
    export <- c('time_series_correlation', 'GLOBAL_SERIES_SIGNALS')
    packages <- c('rioja', 'dplyr', 'tidyr', 'nest', 'zoo', 'RCurl')
    
    cl <- detectCores()
    registerDoParallel(cl)
    
    corr_list <- foreach(i = 1:(dim - 1), .inorder = FALSE, .export = export, .packages = packages) %:%
      foreach(j = i:dim) %dopar% {
        x <- unnest(data[i, 2])
        y <- unnest(data[j, 2])
        
        xid <- distinct(x, site_id)$site_id
        yid <- distinct(y, site_id)$site_id
        
        corr <- time_series_correlation(select(x, site_id, mixed_age, !!corr_signal),
                                        select(y, site_id, mixed_age, !!corr_signal),
                                        corr_signal = corr_signal, 
                                        interp = interp, 
                                        gaussbp = gaussbp)
        tibble(from = xid, to = yid, corr = corr[1])
      }
    #stopCluster(cl)
    stopImplicitCluster()
    
    renquote <- function(l) if (class(l)[1] == 'list') lapply(l, renquote) else enquote(l)
    corr_list <- lapply(unlist(renquote(corr_list)), eval)
    corr_list <- bind_rows(corr_list)
  }
  return(corr_list)
}


gen_site_hash <- function(sites, dating,  hres_only = T, site_ids = 'all', site_lats = 'all') {
  sites <- filter_sites(sites = sites, dating = dating, hres_only = hres_only, site_ids = site_ids, site_lats = site_lats)
  sites <- distinct(sites) %>% arrange(site_id)
  hsh <- paste0(paste0(sites$site_id, collapse = '.')) %>%  
    digest::digest(., serialize = F, algo = 'crc32')
  return(hsh)
}

#' @rdname acer.networks
#' @param data 
#' @param dating tibble, ACER sample dating
#' @param sites_loc tibble, ACER site locations
#' @param orig_data character, name of the original signal tibble from ACER to be incorporated into cached filenames
#' @param type_data character, has to be a supported type data from ...
#' @param hres_only logical, should only hres sites be considered?
#' @param site_ids numeric vector, only specified sites to consider
#' @param site_lats ?
#' @param sd_one logical, should 
#' @param windows numeric list, gives the different window cut points in years b.p., e.g. list(0, 25000, 64000)
#' @param labels character list, gives the corresponding labeling, e.g. list('0-25', '25-64')
#' @param transform character, c('identity', 'logit', 'probit') indicating initial signal transformation
#' @param detrend list of shape list(activate = logical, method = c('linear', 'gaussbp'), args = list(per1 = numeric, per2 = numeric)) with args considered for gaussbp option only according to \code{\link[nest]{gaussdetr}}
#' @param interp list of shape list(activate = logical, method = c('linear', 'gaussian_kernel')) to indicate interpolation method used for time series regularisation;
#' can easily be extended to additional interpolation methods from \code{\link[rioja]{rioja-package}}.
#' @param write_corr list of shape list(activate = logical - should correlation files be saved to disk?, dir = character - path/to/cached/files)). 
#' Cacheing is strongly recommended for sets larger than around 50 sites. 
#' @param use_cached_corr logical - should cached correlation list be used if present?
#'
#' @return list, list of correlation lists ? 
#' @export
#'
#' @examples
make_set_corr_list_from_dataset <- function(data, dating, sites_loc, orig_data, type_data, hres_only = TRUE, site_ids = 'all', site_lats = 'all',
                                            sd_one = FALSE, windows = GLOBAL_SERIES_WINDOWS, labels = GLOBAL_SERIES_WINDOWS_LABELS, transform = 'probit',
                                            detrend = list(activate = TRUE, method = 'linear', args = NULL),
                                            interp = list(method = 'linear'), write_corr = list(activate = TRUE, dir = DIR_NM_DATA), use_cached_corr = TRUE, ncores = NULL) {
  if (use_cached_corr) {
    if(detrend$activate) {name <- paste('detrend', paste0(detrend$method, paste(detrend$args, collapse = '-')), transform, sep = '_')}else {name <- paste(transform, sep = '_')}
    if(hres_only) {name <- paste('hres', name, sep = '_')}
    if(interp$method != 'linear') {name <- paste(name, interp$method, sep = '_')}
  }
  print(name)
  if (all(class(data) == 'list')) {
    r <- 1
    all_corr_list <- list()
    for (trial in data) {
      cat(cr(paste('run', r, '\n')))
      if (use_cached_corr) {
        sites_loc <- filter_sites(sites = sites_loc, dating = dating, hres_only = hres_only, site_ids = site_ids, site_lats = site_lats)
        sites_loc <- distinct(sites_loc) %>% arrange(site_id)
        hsh <- paste0(paste0(sites_loc$site_id, collapse = '.')) %>% #, paste0(sites_loc$site_name, collapse = '.')) %>%  
          digest::digest(., serialize = F, algo = 'crc32')
        
        dir_test <- paste(write_corr$dir, 'correlations', 
                          paste0(unlist(windows) / 1000, collapse = '_'), 
                          paste0('corr_', str_replace(orig_data, 'hic_set_', ''), '-prx', hsh),
                          name, sep = '/')
        file_test <- paste(paste(dir_test, paste0('corr_', str_replace(orig_data, 'hic_set_', ''), 'run', r), sep = '/'), DATAFILES_TYPE, sep = '.')
        if (dir.exists(dir_test) & file.exists(file_test)){file_cached <- TRUE}
        else {file_cached <- FALSE}
        #print(file_test)
        #stop()
      }
      if (!file_cached) {
        cat(cr('correlation file not cached, calculating correlations\n'))
        n <- distinct(trial, noise_sample_id) %>% .$noise_sample_id %>% length(.)
        one <- distinct(trial, noise_sample_id) %>% .$noise_sample_id %>% min(.)
        fin <- distinct(trial, noise_sample_id) %>% .$noise_sample_id %>% max(.)
        cl <- detectCores(); if (!is.null(ncores)) {if (ncores > cl) {stop('ncores provided to make_set_corr_list_from_dataset is larger than detected number of cores')} else {cl <- ncores}}
        registerDoParallel(cl)
        
        export <- c('make_corr_list_from_dataset', 'windowed_time_series_correlation', 'qtransform_data', 'window_and_detrend', 'make_corr_list', 'GLOBAL_SERIES_SIGNALS', 
                    'time_series_correlation', 'filter_sites', 'categorize_age_interval_stats', 'compute_age_interval_stats', 'compute_age_intervals', 
                    'logit', 'id_trf', 'cr', 'window_data', 'detrend_data_window', 'detrend_series_linear', 'detrend_series_gaussbp')
        packages <- c('dplyr', 'tibble', 'nest', 'rioja', 'purrr', 'tidyr', 'zoo', 'RCurl', 'broom')
        
        corr_list <- foreach (i = one:fin, .combine = 'list', .inorder = FALSE, .export = export, .packages = packages) %dopar% {
          data_loc <- filter(trial, noise_sample_id == i) %>% select(-noise_sample_id)
          make_corr_list_from_dataset(data = data_loc, dating = dating, sites_loc = sites_loc, type_data = type_data, hres_only = hres_only, site_ids = site_ids, site_lats = site_lats,
                                      sd_one = sd_one, windows = windows, labels = labels, transform = transform,
                                      detrend = detrend, interp = interp, gaussbp = list(activate = FALSE), dopar = FALSE, verbose = FALSE)
        }
        stopImplicitCluster()
        renquote <- function(l) if (class(l)[1] == 'list') lapply(l, renquote) else enquote(l)
        corr_list <- lapply(unlist(renquote(corr_list)), eval)
        corr_list <- bind_rows(corr_list, .id = 'corr_set_id') %>% 
          mutate_at(vars(corr_set_id), as.numeric) %>% 
          mutate(corr_set_id = corr_set_id + (r - 1) * n)
      }
      else {
        cat(cr('reading correlation list from file\n'))
        corr_list <- read_csv(file = file_test, col_names = TRUE, col_types = list(col_double(), col_factor(), col_double(), col_double(), col_double()))
      }
      
      if (write_corr$activate & !file_cached) {
        cat(cr('writing correlation list to disk\n'))
        if (!dir.exists(dir_test)) {dir.create(dir_test, recursive = TRUE)}
        write_csv(x = corr_list, path = file_test, col_names = TRUE)
      }
      all_corr_list[[as.character(r)]] <- corr_list
      r <- r + 1
    }
    return(all_corr_list)
  }
  else {
    if (use_cached_corr) {
      sites_loc <- filter_sites(sites = sites_loc, dating = dating, hres_only = hres_only, site_ids = site_ids, site_lats = site_lats)
      sites_loc <- distinct(sites_loc) %>% arrange(site_id)
      hsh <- paste0(paste0(sites_loc$site_id, collapse = '.'), paste0(sites_loc$site_name, collapse = '.')) %>%  
        digest::digest(., serialize = F, algo = 'crc32')
      
      dir_test <- paste(write_corr$dir, 'correlations', 
                        paste0(unlist(windows) / 1000, collapse = '_'), 
                        paste0('corr_', str_replace(orig_data, 'hic_set_', ''), '-prx', hsh),
                        name, sep = '/')
      file_test <- paste(paste(dir_test, paste('corr', str_replace(orig_data, 'hic_set_', ''), sep = '_'), sep = '/'), DATAFILES_TYPE, sep = '.')
      if (dir.exists(dir_test) & file.exists(file_test)){file_cached <- TRUE}
      else {file_cached <- FALSE}
    }
    
    if (!file_cached) {
      cat(cr('correlation file not cached, calculating correlations\n'))
      n <- distinct(data, noise_sample_id) %>% .$noise_sample_id %>% length(.)
      cl <- detectCores(); if (!is.null(ncores)) {if (ncores > cl) {stop('ncores provided to make_set_corr_list_from_dataset is larger than detected number of cores')} else {cl <- ncores}}
      registerDoParallel(cl)
      
      export <- c('make_corr_list_from_dataset', 'windowed_time_series_correlation', 'qtransform_data', 'window_and_detrend', 'make_corr_list', 'GLOBAL_SERIES_SIGNALS', 
                  'time_series_correlation', 'filter_sites', 'categorize_age_interval_stats', 'compute_age_interval_stats', 'compute_age_intervals', 
                  'logit', 'id_trf', 'cr', 'window_data', 'detrend_data_window', 'detrend_series_linear', 'detrend_series_gaussbp')
      packages <- c('dplyr', 'tibble', 'nest', 'rioja', 'purrr', 'tidyr', 'zoo', 'RCurl', 'broom')
      
      corr_list <- foreach (i = 1:n, .combine = 'list', .inorder = FALSE, .export = export, .packages = packages) %dopar% {
        data_loc <- filter(data, noise_sample_id == i) %>% select(-noise_sample_id)
        make_corr_list_from_dataset(data = data_loc, dating = dating, sites_loc = sites_loc, type_data = type_data, hres_only = hres_only, site_ids = site_ids, site_lats = site_lats,
                                    sd_one = sd_one, windows = windows, labels = labels, transform = transform,
                                    detrend = detrend, interp = interp, gaussbp = list(activate = FALSE), dopar = FALSE, verbose = FALSE)
      }
      stopImplicitCluster()
      renquote <- function(l) if (class(l)[1] == 'list') lapply(l, renquote) else enquote(l)
      corr_list <- lapply(unlist(renquote(corr_list)), eval)
      corr_list <- bind_rows(corr_list, .id = 'corr_set_id') %>% 
        mutate_at(vars(corr_set_id), as.numeric)
    }
    else {
      cat(cr('reading correlation list from file\n'))
      corr_list <- read_csv(file = file_test, col_names = TRUE, col_types = list(col_double(), col_factor(), col_double(), col_double(), col_double()))
    }
    
    if (write_corr$activate & !file_cached) {
      cat(cr('writing correlation list to disk\n'))
      if (!dir.exists(dir_test)) {dir.create(dir_test, recursive = TRUE)}
      write_csv(x = corr_list, path = file_test, col_names = TRUE)
    }
    return(corr_list)
  }
}

# wrappers and helpers -----------
windowed_time_series_correlation <- function(data, windows = GLOBAL_SERIES_WINDOWS, lables = GLOBAL_SERIES_WINDOWS_LABELS, type_data = 'arboreal_pollen',
                                             transform = 'logit', detrend = TRUE, interp = list(method = 'linear'), sd_one = TRUE, gaussbp = list(activate = FALSE), 
                                             dopar = TRUE, verbose = TRUE) {
  # to be used on an entire dataset of proxy data
  # data should have mixed_age
  if (class(transform) != 'list') {
    data <- qtransform_data(data = data, type_data = type_data, transform = transform, sd_one = if_else(detrend$activate == TRUE, true = FALSE, false = sd_one))
    if (verbose) {cat(cr('transformed data if applicable\n'))}
    #print(data)
    data <- window_and_detrend(data = data, windows = windows, labels = lables, type_data = type_data, transform = transform, detrend = detrend, sd_one = sd_one)
    if (verbose) {cat(cr('detrended data\n'))}
    #print(data)
    # correlations
    if (verbose) {cat(cr('calcuating correlation list\n'))}
    data <- data %>% 
      ungroup() %>% 
      #filter(window == '50-100') %>% 
      #unnest(data_win) %>% 
      #select(-window) %>% 
      #make_corr_list(., corr_signal = {if(detrend$activate) {paste('detrend', transform, GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name, sep = '_')}
      #                                  else{paste(transform, GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name, sep = '_')}},
      #               interp = interp, gaussbp = gaussbp, dopar = dopar) %>% print()
      mutate(corr_list = purrr::map(.$data_win,
                                    make_corr_list,
                                    corr_signal = {if(detrend$activate) {paste('detrend', transform, GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name, sep = '_')}
                                      else{paste(transform, GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name, sep = '_')}}, 
                                    #pb = pb, 
                                    interp = interp, 
                                    gaussbp = gaussbp, 
                                    dopar = dopar))
    #stop()
    if (verbose) {cat(cr('calculated correlation list\n'))}
    # easy access: data %>% group_by(window) %>% unnest(corr_list)
  } else if (class(transform) == 'list') {
    if (!all(sort(names(transform)) == sort(unique(data[[GLOBAL_SERIES_SIGNALS[[type_data]]$signal_grouper]])))) {stop('transform and grouper have to meet in windowed_time_series_correlation')}
    data <- group_by(data, !!sym(GLOBAL_SERIES_SIGNALS[[type_data]]$signal_grouper)) %>% 
      nest() %>% 
      mutate(data_trf = purrr::map2(data, !!sym(GLOBAL_SERIES_SIGNALS[[type_data]]$signal_grouper),
                                    function(x, y) qtransform_data(data = x, type_data = type_data, transform = transform[[y]], sd_one = if_else(detrend$activate == TRUE, true = FALSE, false = sd_one)))) %>% 
      select(-data)
    if (verbose) {cat(cr('transformed data\n'))}
    
    data <- data %>% 
      mutate(data_detr = purrr::map2(data_trf, !!sym(GLOBAL_SERIES_SIGNALS[[type_data]]$signal_grouper),
                                     function(x, y) window_and_detrend(data = x, windows = windows, labels = lables, type_data = type_data, transform = transform[[y]], detrend = detrend, sd_one = sd_one))) %>% 
      mutate(data_detr_ren = purrr::map2(data_detr, !!sym(GLOBAL_SERIES_SIGNALS[[type_data]]$signal_grouper), 
                                         function(x, y) {
                                           if(detrend$activate) {
                                             nm <- paste('detrend', transform[[y]], GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name, sep = '_')
                                           } else{
                                             nm <- paste(transform[[y]], GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name, sep = '_')
                                           }
                                           
                                           x <- unnest(x) %>% mutate(site_datan = purrr::map(site_data, 
                                                                                             function(x) rename(x, corr_sign = nm))) %>% 
                                             select(-site_data) %>% rename(site_data = site_datan)
                                           #%>% unnest() %>% rename(corr_sign = nm)
                                           return(x)
                                         }
      )) %>% 
      unnest(data_detr_ren) %>% 
      select(-!!sym(GLOBAL_SERIES_SIGNALS[[type_data]]$signal_grouper)) %>% 
      group_by(window) %>% nest(.key = 'data_win')
    if (verbose) {cat(cr('detrended data\n'))}
    if (verbose) {cat(cr('calcuating correlation list\n'))}
    data <- data %>% 
      ungroup() %>% 
      mutate(corr_list = purrr::map(.$data_win,
                                    make_corr_list,
                                    corr_signal = 'corr_sign', 
                                    #pb = pb, 
                                    interp = interp, 
                                    gaussbp = gaussbp, 
                                    dopar = dopar))
    if (verbose) {cat(cr('calculated correlation list\n'))}
  } else {
    stop('unhandable class of transform given to windowed_time_series_correlation')
  }
  
  
  return(data)
}


# for rioja methods (NOT USED IN PUBLICATION): interpolates to common time axis with larger of the two mean intervals, then cuts series to same window and finally determines ccf
# if no second series is given: acf is calculated
# x and y have to provide mixed age (thus already have to be join of sample_dating with data)
time_series_correlation <- function(x, y, lagging = list(lag_max = 0, xplus = tibble(), xminus = tibble(), yplus = tibble(), yminus = tibble()),
                                    interp = list(method = 'linear', scale = 'sparse'),
                                    type_data = 'arboreal_pollen',
                                    corr_signal = GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name, 
                                    min_smpl_overl = 10, 
                                    min_time_overl = 5000,
                                    gaussbp = list(activate = FALSE, per1 = 1250, per2 = 3750)) {
  # returns corr for a pair of time series
  rioja_interp <- list('linear', 'loess', 'sspline')
  nest_interp <- list('gaussian_kernel')
  scales <- list('sparse', 'dense', 'global') # different scales have no effect with gaussian kernel, gaussian kernel interpolation is 'sparse' always
  dt_global <- 150
  interp_def <- list(method = 'linear', scale = 'sparse')
  interp <- merge.list(interp, interp_def)
  
  if (!(type_data %in% names(GLOBAL_SERIES_SIGNALS))) {stop('unknown type_data provided to homogenize_time_series_time')}
  if (!('mixed_age' %in% names(x))) {stop('time series provided to homogenize_time_series_time does not contain column mixed_age')}
  if (!(interp$method %in% rioja_interp | interp$method %in% nest_interp)) {stop('unknown time interpolation method provided to time_series_correlation')}
  if (!(interp$scale %in% scales) & !str_detect(interp$scale, 'global')) {stop('unknown time interpolation scale provided to time_series_correlation')}
  
  if (str_detect(interp$scale, 'global') & !(interp$scale == 'global')) {
    dt_global <- str_extract_all(interp$scale, '[:digit:]') %>% unlist() %>% paste0(collapse = '') %>% as.numeric()
  }
  
  xs <- arrange(x, mixed_age) %>% 
    select(mixed_age, corr_signal)
  ys <- arrange(y, mixed_age) %>% 
    select(mixed_age, corr_signal)
  
  if((length(distinct(xs, UQ(sym(corr_signal)))[[corr_signal]] == 1) & is.na(xs[[corr_signal]][1])) | (length(distinct(ys, UQ(sym(corr_signal)))[[corr_signal]] == 1) & is.na(ys[[corr_signal]][1]))) {return(NA)}
  
  # test overlap of series, cut to overlap window and determine series to interpolate
  # this is contained in nexcf as well despite min_time_overl
  corr <- NA
  
  t0 <- max(filter(xs, row_number() == 1)$mixed_age, filter(ys, row_number() == 1)$mixed_age)
  t1 <- min(filter(xs, row_number() == n())$mixed_age, filter(ys, row_number() == n())$mixed_age)
  
  if (t1 - t0 <= 0) {
    warning(paste('time series', (distinct(x, site_id))$site_id, 'and', (distinct(y, site_id))$site_id, 'do not overlap, returning NA'))
    return(corr)
  }
  if (t1 - t0 < min_time_overl){
    warning(paste('time series', (distinct(x, site_id))$site_id, 'and', (distinct(y, site_id))$site_id,
                  'do overlap at less than', min_time_overl, 'years, returning NA'))
    return(corr)
  }
  else if (length(filter(xs, mixed_age >= t0 & mixed_age <= t1)$mixed_age) < min_smpl_overl | length(filter(ys, mixed_age >= t0 & mixed_age <= t1)$mixed_age) < min_smpl_overl) {
    warning(paste('time series', (distinct(x, site_id))$site_id, 'and', (distinct(y, site_id))$site_id,
                  'do overlap at less than', min_smpl_overl, 'data points, returning NA'))
    return(corr)
  }
  
  if (interp$method %in% rioja_interp) {
    #xs <- filter(xs, mixed_age >= t0 & mixed_age <= t1)
    #ys <- filter(ys, mixed_age >= t0 & mixed_age <= t1) # should not be used to avoid wrong overlaps
    
    if (interp$scale == 'sparse' | interp$scale == 'dense') {
      xres <- (mutate(xs, diff = lead(mixed_age) - mixed_age) %>% summarise(mean_diff = mean(diff, na.rm = TRUE)))$mean_diff
      yres <- (mutate(ys, diff = lead(mixed_age) - mixed_age) %>% summarise(mean_diff = mean(diff, na.rm = TRUE)))$mean_diff
    }
    
    # interpolate to common time axis with larger/smaller mean time interval (in overlap window) as interval
    if (interp$scale == 'sparse') {dt <- max(xres, yres)}
    else if (interp$scale == 'dense') {dt <- min(xres, yres)}
    else {dt <- dt_global}
    
    # eventual gaussian bandpass
    if(gaussbp$activate) {
      xsc <- gaussbandpass(as.zoo(as.data.frame(xs[[corr_signal]]), order.by = xs$mixed_age), per1 = gaussbp$per1, per2 = gaussbp$per2)$filt
      ysc <- gaussbandpass(as.zoo(as.data.frame(ys[[corr_signal]]), order.by = ys$mixed_age), per1 = gaussbp$per1, per2 = gaussbp$per2)$filt
      if (length(xsc) != length(xs[[corr_signal]]) | length(ysc) != length(ys[[corr_signal]])) {warning('no values from gaussbandpass, returning NA'); return(NA)}
      else{
        xs[[corr_signal]] <- xsc
        ys[[corr_signal]] <- ysc
      }
    }
    
    # normalize time scale
    xs$mixed_age <- xs$mixed_age / dt
    ys$mixed_age <- ys$mixed_age / dt
    
    tint <- seq(from = t0, to = t1, by = dt) / dt
    xs <- as_tibble(interp.dataset(y = xs,
                                   x = xs$mixed_age,
                                   xout = tint,
                                   method = interp$method, 
                                   rep.negt = FALSE)) %>% drop_na(UQ(sym(corr_signal)))
    
    ys <- as_tibble(interp.dataset(y = ys,
                                   x = ys$mixed_age,
                                   xout = tint,
                                   method = interp$method, 
                                   rep.negt = FALSE)) %>% drop_na(UQ(sym(corr_signal)))
    # window again?
    #if(is.null(xs) | is.null(ys)) {return(corr)}
    
    xs <- as.ts(as.zoo(as.data.frame(xs[[corr_signal]]), order.by = tint))
    ys <- as.ts(as.zoo(as.data.frame(ys[[corr_signal]]), order.by = tint))
    
    # to univariate series
    xs <- as.numeric(xs)
    ys <- as.numeric(ys)
    corr <- ccf(x = xs, y = ys, type = 'correlation', plot = FALSE, lag.max = lagging$lag_max)$acf
    #plot(ccf(x = xs, y = ys, type = 'correlation', plot = FALSE))
  }
  else {
    if(gaussbp$activate) {
      xs <- as.zoo(gaussbandpass(as.zoo(as.data.frame(xs[[corr_signal]]), order.by = xs$mixed_age), per1 = gaussbp$per1, per2 = gaussbp$per2)$filt, order.by = xs$mixed_age)
      ys <- as.zoo(gaussbandpass(as.zoo(as.data.frame(ys[[corr_signal]]), order.by = ys$mixed_age), per1 = gaussbp$per1, per2 = gaussbp$per2)$filt, order.by = ys$mixed_age)
    }
    else {
      #xs <- as.numeric(as.ts(as.zoo(as.data.frame(xs[[corr_signal]]), order.by = as.double(xs$mixed_age))))
      #ys <- as.numeric(as.ts(as.zoo(as.data.frame(ys[[corr_signal]]), order.by = as.double(ys$mixed_age)))) #-> strange conversion for all digit zero ages!
      #xs <- as.numeric(arrange(xs, mixed_age)[[corr_signal]]); ys <- as.numeric(arrange(ys, mixed_age)[[corr_signal]]) #-> just wrong object for nexcf!!
      xs <- zoo(x = arrange(xs, mixed_age)[[corr_signal]], order.by = arrange(xs, mixed_age)[['mixed_age']])
      ys <- zoo(x = arrange(y, mixed_age)[[corr_signal]], order.by = arrange(ys, mixed_age)[['mixed_age']])
    }
    corr <- nexcf(x = xs, y = ys, lag = lagging$lag_max, h = 0.25)
  }
  if (lagging$lag_max != 0) {corr <- max(corr)}
  return(corr)
}


make_corr_list_from_dataset <- function(data, dating, sites_loc, type_data, hres_only = TRUE, site_ids = 'all', site_lats = 'all',
                                        sd_one = FALSE, windows = GLOBAL_SERIES_WINDOWS, labels = GLOBAL_SERIES_WINDOWS_LABELS, transform = 'probit', detrend = list(activate = TRUE, method = 'linear', args = NULL),
                                        interp = list(method = 'linear'), gaussbp = FALSE, dopar = TRUE, verbose = TRUE) {
  sites <- filter_sites(sites = sites_loc, dating = dating, hres_only = hres_only, site_ids = site_ids, site_lats = site_lats)
  
  data <- data %>% 
    inner_join(sites, by = 'site_id') %>% 
    inner_join(dating, by = c('site_id', 'sample_id')) %>% 
    {if (is.null(GLOBAL_SERIES_SIGNALS[[type_data]]$signal_grouper)) {
      select(., GLOBAL_SERIES_SIGNALS[[type_data]]$signal_id, 'site_id', 'sample_id', GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name, 'mixed_age')
    } else if (GLOBAL_SERIES_SIGNALS[[type_data]]$signal_grouper %in% names(.)) {
      select(., GLOBAL_SERIES_SIGNALS[[type_data]]$signal_id, 'site_id', 'sample_id', GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name, 'mixed_age', GLOBAL_SERIES_SIGNALS[[type_data]]$signal_grouper)
    } else {
      select(., GLOBAL_SERIES_SIGNALS[[type_data]]$signal_id, 'site_id', 'sample_id', GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name, 'mixed_age')
    }} %>% 
    windowed_time_series_correlation(data = ., windows = windows, lables = labels, type_data = type_data, transform = transform, detrend = detrend, interp = interp,
                                     sd_one = sd_one, gaussbp = gaussbp, dopar = dopar, verbose = verbose) %>% 
    select(window, corr_list) %>% 
    unnest() %>% 
    filter(!is.na(corr))
  return(data)
}


filter_sites <- function(sites, dating, hres_only, site_ids, site_lats) {
  sites <- dating %>% 
    inner_join(sites, by = 'site_id') %>% 
    {if (!site_ids == 'all') {filter(., site_id %in% site_ids)}else {.}} %>% 
    {if (!site_lats == 'all') {filter(., lat > site_lats$lower & lat < site_lats$upper)}else {.}} %>% 
    inner_join(categorize_age_interval_stats(compute_age_interval_stats(dating, ranges = list(list(start = 6000, stop = 22000)))), by = 'site_id') %>%
    #inner_join(categorize_age_interval_stats(compute_age_interval_stats(dating, ranges = list(list(start = 0, stop = 65000)))), by = 'site_id') %>%
    #{if (hres_only) {filter(., status_med_smp_res %in% c('0 - 200', '200 - 300'))}else {.}} %>%
    {if (hres_only) {filter(., status_med_smp_res %in% c('0 - 200', '200 - 300', '300 - 500'))}else {.}} %>% 
    distinct(site_id)
  return(sites)
}


corr_list_to_matrix <- function(corr_list) {
  return(corr_list)
}

# summarising and testing routines
# summary (mean corr, etc.) -----------
# if noise corr_list provided a t_test or Wilcox rank or "ratio" test will be performed for each window-site-site combination
#' @rdname acer.networks
#' @param set_corr_list set of corr_lists
#' @param write 
#' @param set_corr_list_noise set of null model corr_lists
#' @param test_option 'rank' for Wilcoxon rank test and 'ratio' for frequentistic test
#' @param summary_option character. one of 'single', 'set', indicating if one random corr_list from the set or the set mean should be tested against the null model
#'
#' @return
#' @export
calculate_summary_set_corr_list <- function(set_corr_list, write = list(activate = FALSE, dir = DIR_NM_DATA), set_corr_list_noise = NULL, test_option = 'rank', summary_option = 'single') {
  if (write$activate) {stop('writing has yet to be implemented for calculate_summary_set_corr_list')}
  if (!(summary_option %in% c('single', 'set'))) {stop('unknown summary option given to calculate_summary_set_corr_list')}
  tests <- list(t = t_test_corr_distributions, rank = rank_test_corr_distributions, ratio = ratio_test_corr_distributions)
  if (!is.null(set_corr_list_noise) & (!all((test_option %in% c('ratio', 'rank', 't'))) | is.null(test_option))) {stop(paste('no or unsupported test option given to calculate_summary_set_corr_list; supported: ratio, rank, t'))}
  if (all(class(set_corr_list) == 'list')) {set_corr_list <- bind_rows(set_corr_list)}
  if (all(class(set_corr_list_noise) == 'list')) {set_corr_list_noise <- bind_rows(set_corr_list_noise)}
  if (!is.null(set_corr_list_noise) & !(length(set_corr_list$corr) == length(set_corr_list_noise$corr))) {warning('noise and set do not have the same length')}
  if (summary_option == 'set') {
    set_corr_list <- set_corr_list %>% 
      group_by(window, from, to) %>% 
      summarise(mean_corr = mean(corr, na.rm = T), sd_corr = sd(corr, na.rm = T)) %>% 
      rename(corr = mean_corr) %>% 
      ungroup() %>% 
      {if (!is.null(set_corr_list_noise)) {inner_join(., tests[[test_option]](x = set_corr_list, y = set_corr_list_noise), by = c('window', 'from', 'to'))}else {.}} %>% 
      ungroup()
  } else if (summary_option == 'single') {
    tests$ratio <- ratio_test_ap_corr_distributions
    set_corr_list <- set_corr_list %>% filter(corr_set_id == sample(1:length(unique(set_corr_list$corr_set_id)), size = 1)) %>% 
      select(-corr_set_id) %>% 
      {if (!is.null(set_corr_list_noise)) {inner_join(., tests[[test_option]](x = mutate(., corr = abs(corr)), y = set_corr_list_noise), by = c('window', 'from', 'to'))}else {.}} %>% 
      ungroup()
  }
  return(set_corr_list)
}

# ratio test AP against matern noise
#' @rdname acer.networks
#' @param set_corr_list set of corr_lists
#' @param set_corr_list_noise set of null model corr_lists
#' @param write 
#' @param test_option 'rank' for Wilcoxon rank test and 'ratio' for frequentistic test
#' @param ... ignored
#'
#' @return
#' @export
calculate_summary_ap_corr_list <- function(set_corr_list, set_corr_list_noise, write = list(activate = FALSE, dir = DIR_NM_DATA), test_option = 'ratio', ...) {
  if (write$activate) {stop('writing has yet to be implemented for calculate_summary_ap_corr_listS')}
  tests <- list(t = t_test_corr_distributions, rank = rank_test_corr_distributions, ratio = ratio_test_ap_corr_distributions)
  if (all(class(set_corr_list_noise) == 'list')) {set_corr_list_noise <- bind_rows(set_corr_list_noise)}
  set_corr_list <- set_corr_list %>% 
    group_by(window, from, to) %>% 
    inner_join(., tests[[test_option]](x = mutate(set_corr_list, corr = abs(corr)), y = set_corr_list_noise), by = c('window', 'from', 'to')) %>% 
    ungroup()
  return(set_corr_list)
}


# network measures applied to corr_list ----------
#' @name acer.network.measures
#' @export
NULL

#' @rdname acer.network.measures
#' @param corr_list tibble or list of tibble. each a corr_list summary containing a column diff_sig indicating significant links, corr_list's are assumed to be of implicit upper triangular shape, thus links are uniquely defined
#'
#' @return tibble. adjacency list
#' @export
adj_list_from_corr_list <- function(corr_list) {
  func <- function(c) {
    if (!'diff_sig' %in% names(c)) {stop('no column diff_sig in given corr_list summary, call a test and / or summary method first')}
    if (!('adj' %in% names(c))) {
      if ('diff_sgn' %in% names(c)) {
        c <- mutate(c, adj = if_else(diff_sig == TRUE & diff_sgn == 1, true = 1, false = 0))
      } else {
        c <- mutate(c, adj = if_else(diff_sig, true = 1, false = 0))
      }
    } else {} # pass if adjacency is already computed for list
    return(c)
  }
  
  if (all(class(corr_list) == 'list')) {
    corr_list <-  lapply(corr_list, func)
  } else {
    corr_list <- func(corr_list)
  }
  
  return(corr_list)
}


adj_list_from_corr_list_pn <- function(corr_list) {
  func <- function(c) {
    if (!'diff_sig' %in% names(c)) {stop('no column diff_sig in given corr_list summary, call a test and / or summary method first')}
    if (!('adj' %in% names(c))) {
      if ('diff_sgn' %in% names(c)) {
        c <- mutate(c, adj = case_when(corr >= 0 & diff_sig & diff_sgn == 1 ~ 1, corr < 0 & diff_sig & diff_sgn == 1 ~ -1, TRUE ~ 0))
      } else {
        c <- mutate(c, adj = case_when(corr >= 0 & diff_sig ~ 1, corr < 0 & diff_sig ~ -1, TRUE ~ 0))
      }
    } else {} # pass if adjacency is already computed for list
    return(c)
  }
  
  if (all(class(corr_list) == 'list')) {
    corr_list <-  lapply(corr_list, func)
  } else {
    corr_list <- func(corr_list)
  }
  
  return(corr_list)
}


#' @rdname acer.network.measures
#' @param corr_list tibble or list of tibble. each a corr_list summary containing a column diff_sig indicating significant links, corr_list's are assumed to be of implicit upper triangular shape, thus links are uniquely defined
#' @param id_type character. one of c('site_id', 'node_id') to indicate whether 'from' and 'to' column in corr_list are original site_id's or modified node_id's from a tidygraph object
#' @param save_to_file list(activate = logical, filename = character) option to save table as .tex table to file
#'
#' @return tibble. gives link density by window in an individual column for each corr_list in input
#' @export
compute_nw_link_dens <- function(corr_list, id_type = 'site_id', save_to_file = list(activate = FALSE, filename = file.path(DIR_EVALUATIONS, 'nw_link_dens'))) {
  if (!(id_type %in% c('site_id', 'node_id'))) {stop('unknown id_type given to calculate_nw_link_dens')}
  if (all(class(corr_list) == 'list') & is.null(names(corr_list))) {stop('corr_list must be named list for multiple corr_lists')}
  func <- function(c) {
    c <- c %>% filter(to != from) %>% group_by(window) %>% summarise(sum_adj = sum(adj), .groups = 'drop') %>% 
      inner_join(bind_rows(c %>% group_by(window) %>% distinct(from, .keep_all = T) %>% select(from, window) %>% rename(UQ(sym(id_type)) := from), 
                           c %>% group_by(window) %>% distinct(to, .keep_all = T) %>% select(to, window) %>% rename(UQ(sym(id_type)) := to)) %>% 
                   group_by(window) %>% distinct(UQ(sym(id_type)), .keep_all = TRUE) %>% summarise(nnode = n(), .groups = 'drop'), 
                 by = 'window') %>% 
      inner_join(bind_rows(c %>% filter(diff_sig) %>% group_by(window) %>% distinct(from, .keep_all = T) %>% select(from, window) %>% rename(UQ(sym(id_type)) := from), 
                           c %>% filter(diff_sig) %>% group_by(window) %>% distinct(to, .keep_all = T) %>% select(to, window) %>% rename(UQ(sym(id_type)) := to)) %>% 
                   group_by(window) %>% distinct(UQ(sym(id_type)), .keep_all = TRUE) %>% summarise(nnode_1sign = n(), .groups = 'drop'), 
                 by = 'window') %>% 
      mutate(link_dens = 2 * sum_adj / (nnode * (nnode - 1)))# %>% ungroup()
    return(c)
  }
  
  if (all(class(corr_list) == 'list')) {
    corr_list <- adj_list_from_corr_list(corr_list) %>% lapply(., func)
    corr_list <- lapply(1:length(corr_list), function(i, c, n) {rename(c[[i]],
                                                                       UQ(sym(paste0('link_dens.', n[[i]]))) := link_dens,
                                                                       UQ(sym(paste0('sum_adj.', n[[i]]))) := sum_adj, 
                                                                       UQ(sym(paste0('nnode.', n[[i]]))) := nnode)}, c = corr_list, n = names(corr_list)) %>% 
      plyr::join_all(., by = 'window', type = 'full') %>% as_tibble()
  } else {
    corr_list <- adj_list_from_corr_list(corr_list) %>% func()
  }
  if (save_to_file$activate) {
    stop()
    #save_to_file <- merge.list(save_to_file, list(filename = file.path(DIR_EVALUATIONS, 'nw_link_dens')))
    #print(xtable(corr_list), file = paste0(save_to_file$filename, '.tex'))
  }
  return(corr_list)
}

#' @rdname acer.network.measures
#' @param corr_list tibble or list of tibble. each a corr_list summary containing a column diff_sig indicating significant links, corr_list's are assumed to be of implicit upper triangular shape, thus links are uniquely defined
#' @param id_type character. one of c('site_id', 'node_id') to indicate whether 'from' and 'to' column in corr_list are original site_id's or modified node_id's from a tidygraph object
#' @param save_to_file list(activate = logical, filename = character) option to save table as .tex table to file
#'
#' @return tibble. gives link density by window in an individual column for each corr_list in input
#' @export
compute_nw_link_dens_pn <- function(corr_list, id_type = 'site_id', save_to_file = list(activate = FALSE, filename = file.path(DIR_EVALUATIONS, 'nw_link_dens_pn'))) {
  corr_list <- bind_rows(corr_list %>% 
                           filter(corr >= 0) %>% 
                           compute_nw_link_dens(., id_type) %>% 
                           mutate(ld_sgn = 'ldp'),
                         corr_list %>% 
                           filter(corr < 0) %>% 
                           compute_nw_link_dens(., id_type) %>% 
                           mutate(ld_sgn = 'ldn'))
  # make sure LDs reference to same fully connected network
  corr_list <- corr_list %>% 
    mutate(nnode = max(corr_list$nnode)) %>% 
    mutate(link_dens = 2*sum_adj/(nnode*(nnode-1)))
  if (save_to_file$activate) {
    stop()
    #save_to_file <- merge.list(save_to_file, list(filename = file.path(DIR_EVALUATIONS, 'nw_link_dens')))
    #print(xtable(corr_list), file = paste0(save_to_file$filename, '.tex'))
  }
  return(corr_list)
}


#' @rdname acer.network.measures
#' @param corr_list tibble or list of tibble. each a corr_list summary containing a column diff_sig indicating significant links, corr_list's are assumed to be of implicit upper triangular shape, thus links are uniquely defined
#' @param sites_match tibble. only necessary for id_type = 'node_id' and regions tibble based on 'site_id'; giving the assignment of sites to nodes
#' @param regions tibble. giving the assignment of site_id or node_id to regions, both as one column in a tibble
#' @param id_type character. one of c('site_id', 'node_id') to indicate whether 'from' and 'to' column in corr_list are original site_id's or modified node_id's from a tidygraph object
#' @param self logical. compute intra-region links? 
#' @param save_to_file list(activate = logical, filename = character) option to save table as .tex table to file
#'
#' @return tibble. giving cross link ratio, number of cross links and cross link probability for each region pair and time window contained in corr_list(s)
#' @export
compute_nw_crosslink_ratio <- function(corr_list, 
                                       sites_match = NULL, 
                                       regions = REGIONS, 
                                       id_type = 'site_id', 
                                       self = FALSE,
                                       #split_pn = FALSE,
                                       save_to_file = list(activate = FALSE, filename = file.path(DIR_EVALUATIONS, 'nw_cl_ratio'))) {
  if (!(id_type %in% c('site_id', 'node_id'))) {stop('unknown id_type given to calculate_nw_crosslink_ratio')}
  if (all(class(corr_list) == 'list') & is.null(names(corr_list))) {stop('corr_list must be named list for multiple corr_lists')}
  if (id_type == 'node_id' & !('node_id' %in% names(regions)) & is.null(sites_match)) {stop('need to either provide node_id to site_id match as sites_match or regions with node_id column')}
  if (id_type == 'node_id' & !('node_id' %in% names(regions))) {
    regions <- corr_list %>% inner_join(sites_match, by = 'site_id') %>% select(node_id, region)
  }
  
  func <- function(c) {
    ldens <- compute_nw_link_dens(c, id_type = id_type) %>% 
      select(window, link_dens)
    c <- adj_list_from_corr_list(c) %>% 
      filter(to != from) %>% 
      inner_join(rename(regions, from = UQ(sym(id_type)), from_reg = region), by = 'from') %>% 
      inner_join(rename(regions, to = UQ(sym(id_type)), to_reg = region), by = 'to') %>% 
      rowwise() %>% 
      mutate(cross_reg1 = if_else(from_reg <= to_reg, from_reg, to_reg)) %>% 
      mutate(cross_reg2 = if_else(from_reg <= to_reg, to_reg, from_reg)) %>% 
      ungroup() %>% 
      #{if (!self) {filter(.,cross_reg1 != cross_reg2)} else {.}} %>% # need equal case for correct normalisation
      select(window, from, to, adj, from_reg, to_reg, cross_reg1, cross_reg2) # corr (for different normalisation below)
    nnodes <- bind_rows(c %>% 
                          group_by(window, from_reg) %>% 
                          distinct(from, .keep_all = T) %>% 
                          select(from, from_reg, window) %>% 
                          rename(UQ(sym(id_type)) := from, reg = from_reg),
                        c %>% 
                          group_by(window, to_reg) %>% 
                          distinct(to, .keep_all = T) %>% 
                          select(to, to_reg, window) %>% 
                          rename(UQ(sym(id_type)) := to, reg = to_reg)) %>% 
      group_by(window, reg) %>% 
      distinct(UQ(sym(id_type)), .keep_all = TRUE) %>% 
      summarise(nnode = n())
    fct <- function(c) {
      c <- c %>% select(window, cross_reg1, cross_reg2, adj) %>% 
        {if (!self) {filter(.,cross_reg1 != cross_reg2)} else {.}} %>% 
        group_by(window, cross_reg1, cross_reg2) %>% 
        summarise(ncross_links = sum(adj)) %>% 
        inner_join(rename(nnodes, cross_reg1 = reg, nnode1 = nnode), by = c('window', 'cross_reg1')) %>% 
        inner_join(rename(nnodes, cross_reg2 = reg, nnode2 = nnode), by = c('window', 'cross_reg2')) %>% 
        {if (!self) {
          mutate(., prcross_link = ncross_links / (nnode1 * nnode2))
        } else {
          mutate(., prcross_link = if_else(cross_reg1 == cross_reg2, ncross_links / (nnode1 * (nnode2-1)), ncross_links / (nnode1 * nnode2))) }
        } %>% # self-links excluded above
        inner_join(ldens, by = 'window') %>% 
        mutate(crosslink_ratio = prcross_link / link_dens) %>% 
        select(-nnode1, -nnode2, -link_dens) %>% 
        ungroup()
      return(c)
    }
    # This would normalize on all links (for manuscript we subdivide adjacency matrix)
    #if (split_pn) {
    #  cc <- lapply(c(-1,1), function(s) {
    #    if (s==-1) {
    #      c %>% 
    #        filter(corr < 0) %>% 
    #        fct()
    #    } else {
    #      c %>% 
    #        filter(corr >= 0) %>% 
    #        fct()
    #    }
    #    
    #  })
    #  
    #  if (nrow(cc[[1]]) > 0 & nrow(cc[[2]]) > 0) {
    #    c <- full_join(cc[[2]] %>% rename(pclp = prcross_link) %>% select(-ncross_links, -crosslink_ratio),
    #                   cc[[1]] %>% rename(pcln = prcross_link) %>% select(-ncross_links, -crosslink_ratio),
    #                   by = c('window', 'cross_reg1', 'cross_reg2')) %>% 
    #      mutate(pclp = replace_na(pclp,0.),pcln = replace_na(pcln,0.))
    #  } else if (nrow(cc[[2]]) > 0) {
    #    c <- cc[[2]] %>% 
    #      select(-ncross_links, -crosslink_ratio) %>% 
    #      rename(pclp = prcross_link) %>% 
    #      mutate(pcln = 0)
    #  } else {
    #    c <- cc[[1]] %>% 
    #      select(-ncross_links, -crosslink_ratio) %>% 
    #      rename(pcln = prcross_link) %>% 
    #      mutate(pclp = 0)
    #  }
    #} else {
    c <- fct(c)
    #}
    
    return(c)
  }
  
  if (all(class(corr_list) == 'list')) {
    corr_list <- lapply(corr_list, func)
    corr_list <- lapply(1:length(corr_list), function(i, c, n) {rename(c[[i]],
                                                                       UQ(sym(paste0('ncross_links.', n[[i]]))) := ncross_links,
                                                                       UQ(sym(paste0('prcross_link.', n[[i]]))) := prcross_link, 
                                                                       UQ(sym(paste0('crosslink_ratio.', n[[i]]))) := crosslink_ratio)}, c = corr_list, n = names(corr_list)) %>% 
      plyr::join_all(., by = c('window', 'cross_reg1', 'cross_reg2'), type = 'full') %>% as_tibble()
  } else {
    corr_list <- adj_list_from_corr_list(corr_list) %>% func()
  }
  if (save_to_file$activate) {
    stop()
    #save_to_file <- merge.list(save_to_file, list(filename = file.path(DIR_EVALUATIONS, 'nw_cl_ratio')))
    #print(xtable(corr_list), file = paste0(save_to_file$filename, '.tex'))
  }
  return(corr_list)
}

#' @rdname acer.network.measures
#' @param corr_list tibble or list of tibble. each a corr_list summary containing a column diff_sig indicating significant links, corr_list's are assumed to be of implicit upper triangular shape, thus links are uniquely defined
#' @param sites_match tibble. only necessary for id_type = 'node_id' and regions tibble based on 'site_id'; giving the assignment of sites to nodes
#' @param regions tibble. giving the assignment of site_id or node_id to regions, both as one column in a tibble
#' @param id_type character. one of c('site_id', 'node_id') to indicate whether 'from' and 'to' column in corr_list are original site_id's or modified node_id's from a tidygraph object
#' @param self logical. compute intra-region links? 
#' @param save_to_file list(activate = logical, filename = character) option to save table as .tex table to file
#'
#' @return tibble. giving cross link ratio, number of cross links and cross link probability for each region pair and time window contained in corr_list(s)
#' @export
compute_nw_crosslink_prob_pn <- function(corr_list, sites_match = NULL, regions = REGIONS, id_type = 'site_id', self = FALSE, save_to_file = list(activate = FALSE, filename = file.path(DIR_EVALUATIONS, 'nw_cl_ratio'))) {
  if (!(id_type %in% c('site_id', 'node_id'))) {stop('unknown id_type given to calculate_nw_crosslink_ratio_rel')}
  if (all(class(corr_list) == 'list') & is.null(names(corr_list))) {stop('corr_list must be named list for multiple corr_lists')}
  if (id_type == 'node_id' & !('node_id' %in% names(regions)) & is.null(sites_match)) {stop('need to either provide node_id to site_id match as sites_match or regions with node_id column')}
  if (id_type == 'node_id' & !('node_id' %in% names(regions))) {
    regions <- corr_list %>% inner_join(sites_match, by = 'site_id') %>% select(node_id, region)
  }
  func <- function(c) {
    cp <- c %>%
      filter(corr >= 0) %>% 
      compute_nw_crosslink_ratio(corr_list = ., sites_match = sites_match, regions = regions, self = self, id_type = id_type,
                                 save_to_file = list(activate = FALSE))
    cn <- c %>% 
      filter(corr < 0) %>% 
      compute_nw_crosslink_ratio(corr_list = ., sites_match = sites_match, regions = regions, self = self, id_type = id_type,
                                 save_to_file = list(activate = FALSE))
    
    if (nrow(cp) > 0 & nrow(cn) > 0) {
      c <- full_join(cp %>% rename(pclp = prcross_link) %>% select(-ncross_links, -crosslink_ratio),
                     cn %>% rename(pcln = prcross_link) %>% select(-ncross_links, -crosslink_ratio),
                     by = c('window', 'cross_reg1', 'cross_reg2')) %>% 
        mutate(pclp = replace_na(pclp,0.),pcln = replace_na(pcln,0.))
    } else if (nrow(cp) > 0) {
      c <- cp %>% 
        select(-ncross_links, -crosslink_ratio) %>% 
        rename(pclp = prcross_link) %>% 
        mutate(pcln = 0)
    } else {
      c <- cn %>% 
        select(-ncross_links, -crosslink_ratio) %>% 
        rename(pcln = prcross_link) %>% 
        mutate(pclp = 0)
    }
    # for normalisation on entire adjacency matrix comment above part and uncomment this (see compute_nw_crosslink_ratio; not used in manuscript)
    #c <- compute_nw_crosslink_ratio(corr_list = c, sites_match = sites_match, regions = regions, self = self, id_type = id_type,
    #                                split_pn = TRUE,
    #                                save_to_file = list(activate = FALSE))
    return(c)
  }
  
  if (all(class(corr_list) == 'list')) {
    corr_list <- lapply(corr_list, func)
  } else {
    corr_list <- corr_list %>% func()
  }
  if (save_to_file$activate) {
    stop()
    #save_to_file <- merge.list(save_to_file, list(filename = file.path(DIR_EVALUATIONS, 'nw_cl_prob')))
    #print(xtable(corr_list), file = paste0(save_to_file$filename, '.tex'))
  }
  return(corr_list)
}


#' @rdname acer.network.measures
#' @param corr_list tibble or list of tibble. each a corr_list summary containing a column diff_sig indicating significant links, corr_list's are assumed to be of implicit upper triangular shape, thus links are uniquely defined
#' @param id_type character. one of c('site_id', 'node_id') to indicate whether 'from' and 'to' column in corr_list are original site_id's or modified node_id's from a tidygraph object
#' @param save_to_file list(activate = logical, filename = character) option to save table as .tex table to file
#'
#' @return tibble. gives node degree by window and node in an individual column for each corr_list in input
#' @export
compute_nw_node_degree <- function(corr_list, id_type = 'site_id', save_to_file = list(activate = FALSE, filename = file.path(DIR_EVALUATIONS, 'nw_nd_degr'))) {
  if (!(id_type %in% c('site_id', 'node_id'))) {stop('unknown id_type given to calculate_nw_node_degree')}
  if (all(class(corr_list) == 'list') & is.null(names(corr_list))) {stop('corr_list must be named list for multiple corr_lists')}
  func <- function(c) {
    c <- adj_list_from_corr_list(c) %>% filter(to != from) %>% 
      select(window, from, to, adj) %>% 
      bind_rows(., rename(., from_ = to, to_ = from) %>% rename(from = from_, to = to_)) %>% 
      group_by(window) %>% distinct(from, to, .keep_all = TRUE) %>% filter(from != to) %>% 
      select(-to) %>% rename(UQ(sym(id_type)) := from) %>% 
      group_by(window, UQ(sym(id_type))) %>% summarise(node_degree = sum(adj)) %>% ungroup()
    return(c)
  }
  
  if (all(class(corr_list) == 'list')) {
    corr_list <-  adj_list_from_corr_list(corr_list) %>% lapply(., func)
    corr_list <- lapply(1:length(corr_list), function(i, c, n) {rename(c[[i]], UQ(sym(paste0('node_degree.', n[[i]]))) := node_degree)}, c = corr_list, n = names(corr_list)) %>% 
      plyr::join_all(., by = c('window', id_type), type = 'full') %>% as_tibble()
  } else {
    corr_list <- adj_list_from_corr_list(corr_list) %>% func()
  }
  if (save_to_file$activate) {
    stop()
    #save_to_file <- merge.list(save_to_file, list(filename = file.path(DIR_EVALUATIONS, 'nw_nd_degr')))
    #print(xtable(corr_list), file = paste0(save_to_file$filename, '.tex'))
  }
  return(corr_list)
}

#' @rdname acer.network.measures
#' @param corr_list tibble or list of tibble. each a corr_list summary containing a column diff_sig indicating significant links, corr_list's are assumed to be of implicit upper triangular shape, thus links are uniquely defined
#' @param id_type character. one of c('site_id', 'node_id') to indicate whether 'from' and 'to' column in corr_list are original site_id's or modified node_id's from a tidygraph object
#' @param save_to_file list(activate = logical, filename = character) option to save table as .tex table to file
#'
#' @return tibble. gives node degree by window and node in an individual column for each corr_list in input
#' @export
compute_nw_node_degree_pn <- function(corr_list, id_type = 'site_id', save_to_file = list(activate = FALSE, filename = file.path(DIR_EVALUATIONS, 'nw_nd_degr'))) {
  if (!(id_type %in% c('site_id', 'node_id'))) {stop('unknown id_type given to calculate_nw_node_degree')}
  if (all(class(corr_list) == 'list') & is.null(names(corr_list))) {stop('corr_list must be named list for multiple corr_lists')}
  func <- function(c) {
    cp <- c %>%
      filter(corr >= 0) %>% 
      compute_nw_node_degree(corr_list = ., id_type = id_type, save_to_file = list(activate = FALSE))
    cn <- c %>% 
      filter(corr < 0) %>% 
      compute_nw_node_degree(corr_list = ., id_type = id_type, save_to_file = list(activate = FALSE))
    c <- full_join(cp %>% rename(ndgp = node_degree),
                   cn %>% rename(ndgn = node_degree),
                   by = c('window', 'site_id'))  %>% 
      mutate(ndgp = if_else(is.na(ndgp), 0, ndgp),
             ndgn = if_else(is.na(ndgn), 0, ndgn),
             ndg = ndgp + ndgn)
    return(c)
  }
  
  if (all(class(corr_list) == 'list')) {
    corr_list <- lapply(corr_list, func)
  } else {
    corr_list <- corr_list %>% func()
  }
  if (save_to_file$activate) {
    stop()
    #save_to_file <- merge.list(save_to_file, list(filename = file.path(DIR_EVALUATIONS, 'nw_nd_degr_pn')))
    #print(xtable(corr_list), file = paste0(save_to_file$filename, '.tex'))
  }
  return(corr_list)
}


# contingency table of two networks
compute_nw_contingency_table <- function(a,b,mode='global',normalize=FALSE,adj_vals = c(-1,0,1)) {
  modes <- c('global', 'sitewise')
  if (!(mode %in% modes)) stop(paste0('mode ',mode,' not supported'))
  ct <- full_join(a %>% 
                    adj_list_from_corr_list_pn() %>% 
                    select(-corr_trial_g_noise,-corr_trial_s_noise,-corr_ratio_trial_greater_noise,-corr,-diff_sig) %>% 
                    #mutate(name = a_name) %>% 
                    rename(adj.a = adj),
                  b %>% 
                    adj_list_from_corr_list_pn() %>% 
                    select(-corr_trial_g_noise,-corr_trial_s_noise,-corr_ratio_trial_greater_noise,-corr,-diff_sig) %>% 
                    #mutate(name = b_name) %>% 
                    rename(adj.b = adj),
                  by = c('window', 'from', 'to')) %>% 
    group_by(window)
  
  fct <- function(a,b) {
    # this order gives table with rows ~ a, cols ~ b
    contingency_table <- sapply(adj_vals,function(x) sapply(adj_vals, function(y) length(which(b[which(b %in% adj_vals & a %in% adj_vals)] == x & a[which(b %in% adj_vals & a %in% adj_vals)] == y)))    )
    if (normalize==TRUE) {
      contingency_table <- contingency_table / length(which(a %in% adj_vals & b %in% adj_vals))
    }
    return(contingency_table)
  }
  
  if (mode == 'global') {
    ct <- ct %>%
      summarise(cont_table = list(fct(.$adj.a,.$adj.b)), .groups = 'drop')
  }
  
  return(ct)
}


# Cohen kappa statistics on the contingency table of two networks
cohens_kappa <- function(contingency_table) {
  p_0 <- sum(diag(contingency_table))/sum(c(contingency_table))
  p_c <- sum(rowSums(contingency_table)*colSums(contingency_table))/sum(c(contingency_table))^2
  se <- sqrt(p_0*(1-p_0)/(sum(c(contingency_table))*(1-p_c)^2))
  kappa_coeff <- (p_0-p_c)/(1-p_c)
  return(c(kappa_coeff,kappa_coeff-1.965*se,kappa_coeff+1.965*se))
}

# tests ----------
# t test of two correlation lists; summary will be returned with added diff_sig column
#' Significance tests links based on correlation scores
#'
#' @name acer.networks.test
#' @export
NULL

t_test_corr_distributions <- function(x, y, conf.level = 0.95) {
  x <- rename(y, corr.y = corr) %>% 
    inner_join(x, by = c('corr_set_id', 'window', 'from', 'to')) %>% 
    select(-corr_set_id) %>% 
    mutate_at(vars(corr, corr.y), as.numeric) %>% 
    filter(from != to) %>% 
    group_by(window, from, to) %>% 
    nest(.key = 'corr_vals') %>% 
    mutate(t_test_ = purrr::map(corr_vals, ~ t.test(x = .x$corr, y = .x$corr.y, conf.level = conf.level, alternative = 'two.sided')),
           t_test = purrr::map(t_test_, ~ .$p.value %>% as_tibble())) %>% 
    select(-t_test_, -corr_vals) %>%
    unnest() %>% 
    mutate(diff_sig = if_else(value < (1 - conf.level), TRUE, FALSE)) %>% 
    select(-value)
  return(x)
}


rank_test_corr_distributions <- function(x, y, conf.level = 0.95) {
  if ('corr_set_id' %in% names(x)) {bys <- c('corr_set_id', 'window', 'from', 'to')}else {bys <- c('window', 'from', 'to')}
  x <- rename(y, corr.y = corr) %>% 
    inner_join(x, by = bys) %>% 
    select(-corr_set_id) %>% 
    mutate_at(vars(corr, corr.y), as.numeric) %>% 
    filter(from != to) %>% 
    group_by(window, from, to) %>% 
    nest(.key = 'corr_vals') %>% 
    mutate(rank_test_ = purrr::map(corr_vals, ~ wilcox.test(x = .x$corr, y = .x$corr.y, conf.level = conf.level, alternative = 'two.sided', paired = FALSE)),
           rank_test = purrr::map(rank_test_, ~ .$p.value %>% as_tibble())) %>% 
    select(-rank_test_, -corr_vals) %>%
    unnest() %>% 
    mutate(diff_sig = if_else(value < (1 - conf.level), TRUE, FALSE)) %>% 
    select(-value)
  return(x)
}


#' test based on frequentistic ratio
#'
#' @rdname acer.networks.test
#' @param x tibble. corr_list to be tested
#' @param y tibble. corr_list set from null model
#' @param conf.level double. confidence level
#' @param neg_sgn_as_sig logical. should correlation scores deviating significantly from null model with negative signum be considered as significant
#'
#' @return tibble. original corr_list with added column indicating significant links
ratio_test_corr_distributions <- function(x, y, conf.level = 0.95, neg_sgn_as_sig = FALSE) {
  x <- rename(y, corr.y = corr) %>% 
    inner_join(x, by = c('corr_set_id', 'window', 'from', 'to')) %>% 
    select(-corr_set_id) %>% 
    mutate_at(vars(corr, corr.y), as.numeric) %>% 
    filter(from != to) %>% 
    mutate(corr_trial_g_noise_ = if_else(corr.y < corr, TRUE, FALSE)) %>% 
    select(-corr.y, -corr) %>% 
    group_by(window, from, to)
  x <- full_join(filter(x, corr_trial_g_noise_ == TRUE) %>% count() %>% rename(corr_trial_g_noise_ = n), 
                 filter(x, corr_trial_g_noise_ == FALSE) %>% count() %>% rename(corr_trial_s_noise_ = n), 
                 by = c('window', 'from', 'to')) %>% 
    ungroup() %>%
    group_by(1:n()) %>% 
    mutate_at(vars(corr_trial_g_noise_, corr_trial_s_noise_), function(x) {if(is.na(x)) {0}else{x}}) %>% 
    mutate(corr_trial_g_noise = corr_trial_g_noise_ / (corr_trial_g_noise_ + corr_trial_s_noise_), corr_trial_s_noise = corr_trial_s_noise_ / (corr_trial_g_noise_ + corr_trial_s_noise_)) %>% 
    {if (neg_sgn_as_sig) {
      mutate(., diff_sig = if_else(corr_trial_g_noise >= conf.level | corr_trial_s_noise >= conf.level, TRUE, FALSE), diff_sgn = if_else(diff_sig, if_else(corr_trial_g_noise >= conf.level, 1, -1), NA_real_))
    } else {
      mutate(., diff_sig = if_else(corr_trial_g_noise >= conf.level, TRUE, FALSE))
    }} %>%  
    select(-corr_trial_g_noise_, -corr_trial_s_noise_) %>% 
    ungroup() %>% 
    mutate(corr_ratio_trial_greater_noise = corr_trial_g_noise / corr_trial_s_noise) %>% 
    select(-UQ(sym('1:n()')))
  return(x)
}


#' test based on frequentistic ratio
#'
#' @rdname acer.networks.test
#' @param x tibble. corr_list to be tested
#' @param y tibble. corr_list set from null model
#' @param conf.level double. confidence level
#' @param neg_sgn_as_sig logical. should correlation scores deviating significantly from null model with negative signum be considered as significant
#'
#' @return tibble. original corr_list with added column indicating significant links
ratio_test_ap_corr_distributions <- function(x, y, conf.level = 0.95, neg_sgn_as_sig = FALSE) {
  x <- rename(y, corr.y = corr) %>% 
    inner_join(x, by = c('window', 'from', 'to')) %>% 
    mutate_at(vars(corr, corr.y), as.numeric) %>% 
    filter(from != to) %>% 
    mutate(corr_trial_g_noise_ = if_else(corr.y < corr, TRUE, FALSE)) %>% 
    select(-corr.y, -corr) %>% 
    group_by(window, from, to)
  x <- full_join(filter(x, corr_trial_g_noise_ == TRUE) %>% count() %>% rename(corr_trial_g_noise_ = n), 
                 filter(x, corr_trial_g_noise_ == FALSE) %>% count() %>% rename(corr_trial_s_noise_ = n), 
                 by = c('window', 'from', 'to')) %>% 
    ungroup() %>%
    group_by(1:n()) %>% 
    mutate_at(vars(corr_trial_g_noise_, corr_trial_s_noise_), function(x) {if(is.na(x)) {0}else{x}}) %>% 
    mutate(corr_trial_g_noise = corr_trial_g_noise_ / (corr_trial_g_noise_ + corr_trial_s_noise_), corr_trial_s_noise = corr_trial_s_noise_ / (corr_trial_g_noise_ + corr_trial_s_noise_)) %>% 
    {if (neg_sgn_as_sig) {
      mutate(., diff_sig = if_else(corr_trial_g_noise >= conf.level | corr_trial_s_noise >= conf.level, TRUE, FALSE), diff_sgn = if_else(diff_sig, if_else(corr_trial_g_noise >= conf.level, 1, -1), NA_real_))
    } else {
      mutate(., diff_sig = if_else(corr_trial_g_noise >= conf.level, TRUE, FALSE))
    }} %>%  
    select(-corr_trial_g_noise_, -corr_trial_s_noise_) %>% 
    ungroup() %>% 
    mutate(corr_ratio_trial_greater_noise = corr_trial_g_noise / corr_trial_s_noise) %>% 
    select(-UQ(sym('1:n()')))
  return(x)
}


# network measure tests
#' test network measure based on randomisation and frequentistic ratio
#'
#' @rdname acer.networks.test
#' @param corr_list tibble or list of tibble to be tested. each a corr_list summary containing a column diff_sig indicating significant links, corr_list's are assumed to be of implicit upper triangular shape, thus links are uniquely defined
#' @param nw_meas character. one of c('link_dens', 'node_degree', 'crosslink_ratio')
#' @param conf_interv integer between 0 and 1
#' @param save_to_file list(activate = logical, filename = character) option to save table as .tex table to file
#' @param ... passed to correponding network measure
#'
#' @return tibble. 
#' @export
compute_nw_measure_test <- function(corr_list, nw_meas, n = 10, conf_interv = 0.95, save_to_file = list(activate = FALSE, filename = file.path(DIR_EVALUATIONS, paste0('test_', nw_meas))), ...) {
  if (!(nw_meas %in% c('link_dens', 'node_degree', 'crosslink_ratio', 'crosslink_ratio_rel', 'crosslink_ratio_anom', 'crosslink_prob_anom', 'node_degree_rel', 'node_degree_anom'))) {
    stop('unknown network measure given to compute_nw_measure_test')
  }
  if (all(class(corr_list) == 'list') & is.null(names(corr_list))) {stop('corr_list must be named list for multiple corr_lists')}
  nw_measm <- paste0('compute_nw_', nw_meas)
  
  func <- function(cl) {
    clt <- do.call(nw_measm, args = list(corr_list = cl, ...))
    nullm <- tibble(perm_id = 1:n,
                    perm = lapply(1:n, function(i) return(cl %>% group_by(window) %>% nest() %>% mutate(perm = purrr::map(.$data, permute_corr_list)) %>% unnest(perm)))) %>% 
      mutate(meas = purrr::map(.$perm, function(x)return(do.call(nw_measm, args = list(corr_list = x, ...))))) %>% unnest(meas) %>% {if (nw_meas == 'crosslink_ratio') {select(., -ncross_links, -prcross_link)} else {.}} %>% 
      rename(!!sym(paste0(nw_meas, '.p')) := !!sym(nw_meas)) %>% inner_join(clt, by = names(clt)[!(names(clt) %in% c(nw_meas, 'ncross_links', 'prcross_link'))]) %>% 
      mutate(dev = case_when(!!sym(paste0(nw_meas, '.p')) > !!sym(nw_meas) ~ -1, !!sym(paste0(nw_meas, '.p')) < !!sym(nw_meas) ~ 1, !!sym(paste0(nw_meas, '.p')) == !!sym(nw_meas) ~ 0)) %>% 
      group_by(!!!syms(names(.)[!(names(.) %in% c('perm_id', paste0(nw_meas, '.p'), 'dev'))])) %>% nest() %>% 
      mutate(!!sym(paste0(nw_meas, '.p')) := purrr::map(.$data, function(x) return(mean(x[[paste0(nw_meas, '.p')]]))),
             devp = purrr::map(.$data, function(x) return(filter(x, dev == 1) %>% .$dev %>% sum())), 
             devn = purrr::map(.$data, function(x) return(filter(x, dev == -1) %>% .$dev %>% abs() %>% sum())), 
             deve = purrr::map(.$data, function(x) return(filter(x, dev == 0) %>% .$dev %>% length()))) %>% select(-data) %>% unnest() %>%
      rowwise() %>% mutate(diff_sig = if_else(devp / (devp + devn + deve) >= conf_interv | devn / (devp + devn + deve) >= conf_interv, TRUE, FALSE)) %>% ungroup()
    #if (any(str_detect(names(nullm), '.x'))) {print('#'); nullm <- select_at(nullm, vars(-contains('.x'))) %>% rename_if(contains('.y'), ~ str_remove(., '.y'))}
    # compare distrs
    return(nullm)
  }
  
  if (all(class(corr_list) == 'list')) {
    corr_list <- lapply(corr_list, func)
    corr_list <- lapply(1:length(corr_list), function(i, c, n) {rename_at(c[[i]], vars(!!!syms(names(c[[i]])[!(names(c[[i]]) %in% c('window', 'site_id', 'cross_reg1', 'cross_reg2'))])),
                                                                          function(nm){return(paste0(nm, '.', n[[i]]))})},
                        c = corr_list, n = names(corr_list)) %>% 
      plyr::join_all(., by = {if (nw_meas == 'crosslink_ratio') {c('window', 'cross_reg1', 'cross_reg2')} else if (nw_meas == 'node_degree') {c('window', 'site_id')} else {'window'}}, type = 'full') %>% as_tibble()
  } else {
    corr_list <- func(corr_list)
  }
  if (save_to_file$activate) {
    stop()
    #save_to_file <- merge.list(save_to_file, list(filename = file.path(DIR_EVALUATIONS, paste0('nw_', nw_meas, '_test'))))
    #print(xtable(corr_list), file = paste0(save_to_file$filename, '.tex'))
    #write_csv(x = corr_list, path = paste0(save_to_file$filename, '.csv'), col_names = TRUE)
  }
  return(corr_list)
}


permute_corr_list <- function(corr_list) {
  smpl <- NULL
  nms <- select(corr_list, -from, -to) %>% names() 
  all <- bind_rows(distinct(corr_list, to) %>% rename(from = to), distinct(corr_list, from)) %>% distinct(from) %>% arrange(from)
  corr_list <- bind_rows(corr_list, tibble(from = all$from, to = all$from)) %>% arrange(from)
  cls <- lapply(nms, function(nm) {
    cl <- select(corr_list, from, to, UQ(sym(nm)))
    #rnone <- cl$from[1]
    cl <- cl %>% tidyr::complete(from, to) %>% spread(to, UQ(sym(nm)), fill = 10000)
    rn <- cl$from
    cl <- cl %>% select(-from) %>% as.matrix()#; rownames(cl) <- rn
    ut <- upper.tri(cl)
    #print(colnames(cl)); print(rn)
    #print(smpl)
    if (is.null(smpl)) {smpl <<- sample.int(n = length(cl[ut]))}; cl[ut] <- cl[ut] %>% .[smpl]
    cl <- cl %>% as_tibble(rownames = NULL) %>% mutate(from = rn) %>% gather(key = 'to', value = !!nm, -from) %>% mutate_at(vars(to), as.double) %>% arrange(from) %>% filter(UQ(sym(nm)) != 10000)
  }) %>% 
    plyr::join_all(by = c('from', 'to'), type = 'inner') %>% mutate(diff_sig = if_else(diff_sig == 1, TRUE, FALSE)) %>% as_tibble()
  return(cls)
}


# translate between site_id (not numbered continuously) and node_id (numbered continuously for ggraph)
#' @rdname acer.networks
#' @param corr_list tibble. from and to columns based on site_id identifier
#' @param sites 
#' @param region_tree 
#' @param region_locs
#'
#' @return list(corr_list = tibble, sites = tibble). corr_list with from and to columns based on dense node_id identifier for use in a tbl_graph object, sites corresponding tibble to match site_id to node_id
#' @export
combo_site_node_id <- function(corr_list, sites, region_tree = FALSE, regions = REGIONS, region_locs = REGION_LOCS_pn) {
  regions_in <- regions
  sites <- distinct(sites, site_id, lat, long) %>% 
    arrange(site_id)
  if (all(str_detect(region_locs$region, '.p') | str_detect(region_locs$region, '.n')) & region_tree) {
    all(str_detect(region_locs$region, '.p') | str_detect(region_locs$region, '.n'))
    split_reg_by_sign <- T
    sites_offset <- ceiling(max(sites$site_id)) + 1
    sites <- sites %>% 
      bind_rows(sites %>% mutate(site_id = sites_offset + site_id))
    regions_md <- bind_rows(
      regions_in %>% mutate(region = paste0(region, '.p')),
      regions_in %>% mutate(region = paste0(region, '.n'),site_id = site_id + sites_offset)
    )
    sites <- bind_rows(sites, seq(from = ceiling(max(sites$site_id)) + 1, to = ceiling(max(sites$site_id)) + (distinct(regions_md, region) %>% add_row(region = 'earth.p') %>% add_row(region = 'earth.n') %>% .$region %>% length()), by = 1)
                       %>% tibble(site_id = .)) %>% 
      rownames_to_column(var = 'node_id') %>% 
      mutate_at(vars(node_id), as.numeric)
  } else {
    split_reg_by_sign <- F
    sites <- {if (region_tree) {bind_rows(sites, seq(from = ceiling(max(sites$site_id)) + 1, to = ceiling(max(sites$site_id)) + (distinct(regions_in, region) %>% add_row(region = 'earth') %>% .$region %>% length()), by = 1)
                                          %>% tibble(site_id = .))}else{sites}} %>% 
      rownames_to_column(var = 'node_id') %>% 
      mutate_at(vars(node_id), as.numeric)
  }
  
  if (split_reg_by_sign) {
    corr_list <- corr_list %>% 
      mutate(from = if_else(corr >= 0, from, sites_offset + from),
             to = if_else(corr >= 0, to, sites_offset + to)) %>% 
      inner_join(select(sites, site_id, node_id) %>% rename(from = site_id, from_ = node_id), by = 'from') %>% 
      inner_join(select(sites, site_id, node_id) %>% rename(to = site_id, to_ = node_id), by = 'to') %>% 
      select(-from, -to) %>% 
      rename(from = from_, to = to_)
  } else {
    corr_list <- corr_list %>% 
      inner_join(select(sites, site_id, node_id) %>% rename(from = site_id, from_ = node_id), by = 'from') %>% 
      inner_join(select(sites, site_id, node_id) %>% rename(to = site_id, to_ = node_id), by = 'to') %>% 
      select(-from, -to) %>% 
      rename(from = from_, to = to_)
  }
  
  if (region_tree) {
    if (split_reg_by_sign) {
      regions <- filter(sites, is.na(lat) & is.na(long)) %>% select(node_id, site_id) %>% bind_cols(distinct(regions_md, region) %>% add_row(region = 'earth.p') %>% add_row(region = 'earth.n'))
      top_level <- tibble(from = filter(regions, str_detect(region,'earth')) %>% .$node_id %>% rep(., each = length(filter(regions, !str_detect(region,'earth')) %>% .$node_id)/length(filter(regions, str_detect(region,'earth')) %>% .$node_id)), 
                          to = filter(regions, !str_detect(region,'earth')) %>% .$node_id) %>% 
        lapply(unique(corr_list$window), function(x, t) return(mutate(t, window = x)), t = .) %>% bind_rows()
      bottom_level <- filter(regions, !str_detect(region,'earth')) %>% rename(from = node_id) %>% select(-site_id) %>% inner_join(regions_md, by = 'region') %>% 
        inner_join(filter(sites, !is.na(lat) & !is.na(long)) %>% select(node_id, site_id) %>% rename(to = node_id), by = 'site_id') %>% select(from, to)%>% lapply(unique(corr_list$window), function(x, t) return(mutate(t, window = x)), t = .) %>% bind_rows()
    } else {
      regions <- filter(sites, is.na(lat) & is.na(long)) %>% select(node_id, site_id) %>% bind_cols(distinct(regions_in, region) %>% add_row(region = 'earth'))
      top_level <- tibble(from = filter(regions, region == 'earth') %>% .$node_id, to = filter(regions, region != 'earth') %>% .$node_id) %>% lapply(unique(corr_list$window), function(x, t) return(mutate(t, window = x)), t = .) %>% bind_rows()
      bottom_level <- filter(regions, region != 'earth') %>% rename(from = node_id) %>% select(-site_id) %>% inner_join(regions_in, by = 'region') %>% 
        inner_join(filter(sites, !is.na(lat) & !is.na(long)) %>% select(node_id, site_id) %>% rename(to = node_id), by = 'site_id') %>% select(from, to)%>% lapply(unique(corr_list$window), function(x, t) return(mutate(t, window = x)), t = .) %>% bind_rows()
    }
    sites <- inner_join(regions, region_locs, by = 'region') %>% select(-region) %>% mutate(hierarchial = TRUE) %>% bind_rows(filter(sites, !is.na(lat) & !is.na(long)) %>% mutate(hierarchial = FALSE)) %>% arrange(node_id)
    #View(top_level); View(bottom_level); View(sites); View(corr_list)
    combo <- list(nw = bind_rows(list(top_level, bottom_level)), sites = sites, corr_list = corr_list)
  }
  else {
    combo <- list(corr_list = corr_list, sites = sites)
  }
  return(combo)
}
