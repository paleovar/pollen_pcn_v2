# PCNdata S3 class which replaces ACER_dataset RC class from pollen_pcn (v1)
# for easier handling
# and methods of the PCNdata class

# PCNdata class ----
PCNdata <- function() {
  p <- structure(
    list(
    'sites' = load_db_flat('site'), 
    'sample_dating' = load_db_flat('sample_ct_corrected'), 
    'pollen_data' = load_db_flat('pollen_data'),
    pollen_data_percentized = FALSE,
    cores_merged = FALSE
    ), class = "PCNdata")
  return(p)
}

is.PCNdata <- function(x) inherits(x, "PCNdata")


# PCNdata class methods for handling of chronologies ----
## merge the data from different publications belonging to the sediment core
merge_cores_in_db <- function(x, ...) UseMethod('merge_cores_in_db')

merge_cores_in_db.PCNdata <- function(db, in_place = FALSE, type_data_not_in_place = NULL) {
  # merges default cores in database
  types <- list('pollen_data', 'sample_dating') #'biome_percentages', 'charcoal_data',
  types_not_in_place <- list('pollen' = 'pollen_data') #, 'biomes' = 'biome_percentages', 'charcoal' = 'charcoal_data')
  if (in_place) {
    if (db$cores_merged) {
      cat('cores already merged, skipping merge\n')
      return(db)
    }
    for (type_data in types) {
      merge <- merge_cores_in_data(data = db[[type_data]], sites = db$sites)
      db[[type_data]] <- merge$data %>% ungroup()
    }
    db$sites <- merge_cores_in_data(data = db[['pollen_data']], sites = db$sites)$new_sites
    db$cores_merged <- TRUE
    cat('cores in database are merged now\n')
    return(db)
  }
  else {
    stop('removed!')
    #if(is.null(type_data_not_in_place)){stop('no type data provided to merge_cores_in_db')}
    #return(merge_cores_in_data(data = db[[types_not_in_place[[type_data_not_in_place]]]], sites = db$sites)$data %>% ungroup())
  }
}


## concatenate the age models (w/o error where available)
mix_sample_dating <- function(x, ...) UseMethod('mix_sample_dating')

mix_sample_dating.PCNdata <- function(db, in_place = TRUE, orig_age_only = FALSE){
  db <- merge_cores_in_db(db, in_place = TRUE)
  sample_dating <- mix_dating(db$sample_dating, orig_age_only = orig_age_only)
  if (in_place) {
    db$sample_dating <- sample_dating
    return(db)
  }
  else {
    return(sample_dating)
  }
}

## concatenate the age models (w error where available)
mix_sample_dating_and_err <- function(x, ...) UseMethod('mix_sample_dating_and_err')

mix_sample_dating_and_err.PCNdata <- function(db, in_place = TRUE, orig_age_only = FALSE){
  db <- merge_cores_in_db(db, in_place = TRUE)
  sample_dating <- mix_dating_and_err(db$sample_dating, orig_age_only = orig_age_only)
  if (in_place) {
    db$sample_dating <- sample_dating
    return(db)
  }
  else {
    return(sample_dating)
  }
}


# PCNdata class for percentizing
pollen_percentize <- function(x, ...) UseMethod('pollen_percentize')

pollen_percentize.PCNdata <- function(db, overwrite_original_col = FALSE, in_place = TRUE){
  if(!db$pollen_data_percentized){
    pollen <- inner_join(db$pollen_data, select(db$sample_dating, sample_id, count_type), by = 'sample_id')
    
    coun <- pollen[pollen$count_type == 'COUN',] %>% # filters out charcoal samples from sample_dating implicitly
      group_by(sample_id) %>% 
      filter(!is.na(taxon_count)) %>% # filters some "Unknown"/"Indeterminable" counts that are NA in the data
      summarise(total_taxon_count = sum(taxon_count)) %>% 
      full_join(., pollen[pollen$count_type == 'COUN',], by = 'sample_id') %>% 
      filter(total_taxon_count != 0.) %>% # filters samples which contain only fully-zero counts
      {if(overwrite_original_col){
        mutate(., taxon_count = taxon_count / total_taxon_count * 100.0)
      }
        else{
          mutate(., taxon_pcnt = taxon_count / total_taxon_count * 100.0)
        }} %>% select(-total_taxon_count)
    pcnt <- pollen[pollen$count_type %in% list('DIGI', 'PCNT'),] %>% 
      {if(overwrite_original_col){.}
        else{
          mutate(., taxon_pcnt = taxon_count)
        }}
    
    if(overwrite_original_col){
      stop('removed!')
    }
    else{
      pollen <- full_join(coun, pcnt, by = c('pollen_data_id', 'sample_id', 'site_id', 'taxon', 'taxon_count', 'taxon_pcnt', 'count_type')) %>%
        select(pollen_data_id, sample_id, site_id, taxon, taxon_count, taxon_pcnt)
      if(in_place){
        db$pollen_data <- select(pollen, pollen_data_id, sample_id, site_id, taxon, taxon_count, taxon_pcnt)
        db$pollen_data_percentized <- TRUE
        cat('wrote taxon percentages to db$pollen_data\n')
      }
      else{
        stop('removed!')
      }
    }
  }
  else if(in_place){
    cat('taxon percentages already in db$pollen_data, skipping percentize')#'wrote taxon percentages to db$pollen_data\n')
  }
  return(db)
}

# PCNdata class methods for data set extension (e.g. adding icecore data to the same data table & convention as ACER pollen records) ----
extend_dataset_in_db <- function(x, ...) UseMethod('extend_dataset_in_db')

extend_dataset_in_db.PCNdata <- function(db, ext_data, ext_proxy,data_to_ext, 
                                         from_type_data,from_proxy,
                                         to_type_data,to_data_name,to_signal_name, to_signal_id, to_signal_grouper, 
                                         add_site_info = NULL) {
  db <- mix_sample_dating_and_err(db)
  
  sample_id_start <- floor(max(db$sample_dating$sample_id) + 1)
  signal_id_start <- floor(max(db[[data_to_ext]][[GLOBAL_SERIES_SIGNALS[[from_type_data]]$signal_id]]) %>% as.numeric() + 1)
  site_id_start <- floor(max(db$sample_dating$site_id) + 1)
  
  ext_data <- ext_data %>% 
    rename(age = 1) %>% 
    gather(key = 'site_name', value = to_name, -age) %>% rename(!!sym(to_signal_name) := to_name) %>% 
    group_by(site_name) %>% nest() %>% rownames_to_column(var = 'site_id') %>% mutate_at(vars(site_id), as.numeric) %>% 
    mutate(site_id = site_id + site_id_start) %>% unnest() %>% 
    rownames_to_column(var = 'sample_id') %>% 
    mutate_at(vars(sample_id), as.numeric) %>%
    mutate(!!sym(to_signal_id) := sample_id + signal_id_start) %>% 
    mutate(sample_id = sample_id + sample_id_start) %>% 
    {
      if (to_signal_grouper %in% colnames(add_site_info)) {
        inner_join(., select(add_site_info, site_name, to_signal_grouper), by = 'site_name')
      } else {
        mutate(., !!sym(to_signal_grouper) := ext_proxy)
      }
    } %>% 
    ungroup()
  
  new_data_db <- db[[data_to_ext]] %>%
    rename(!!sym(to_signal_id) := !!sym(GLOBAL_SERIES_SIGNALS[[from_type_data]]$signal_id), 
           !!sym(to_signal_name) := !!sym(GLOBAL_SERIES_SIGNALS[[from_type_data]]$signal_name)) %>% 
    mutate_at(vars(!!sym(to_signal_id)), as.numeric) %>% 
    {
      if (to_signal_grouper %in% colnames(db[[data_to_ext]])) {
        .
      } else {
        mutate(., !!sym(to_signal_grouper) := from_proxy)
      }
    } %>% 
    bind_rows(select(.data = ext_data, -age, -site_name))
  
  if (is.null(add_site_info)) {
    add_site_info <- select(ext_data, site_id, site_name) %>% distinct()
  } else {
    add_site_info <- add_site_info %>% inner_join(select(ext_data, site_id, site_name) %>% distinct(), by = 'site_name')
  }
  overlap <- intersect(colnames(add_site_info), colnames(db$sites))
  
  db$ACERdefaultsites <- db$sites
  db$sites <- db$sites %>% 
    full_join(select(add_site_info, overlap), by = overlap) %>% arrange(site_id)
  
  db$sample_dating <- db$sample_dating %>% 
    full_join(select(ext_data, site_id, sample_id, age) %>% rename(mixed_age = age), by = c('site_id', 'sample_id', 'mixed_age'))
  
  db[[to_data_name]] <- new_data_db
  
  GLOBAL_SERIES_SIGNALS[[to_type_data]] <- merge.list(list(signal_name = to_signal_name, signal_id = to_signal_id, signal_grouper = to_signal_grouper), 
                                                       GLOBAL_SERIES_SIGNALS[[to_type_data]])
  return(db)
}


# PCNdata class methods for harmonizing taxa and computing arboreal pollen percentages ----
harmonized_pollen_to_db <- function(x, ...) UseMethod('harmonized_pollen_to_db')

harmonized_pollen_to_db.PCNdata <- function(db, harmonization_level = 2, in_place = TRUE, write_to_flat_file = FALSE) {
  # core-merged data only for harmonized pollen_data
  db <- pollen_percentize(db) %>% 
    merge_cores_in_db(in_place = TRUE)
  harm_pollen_data <- harmonize_taxa_ax(pollen_data = db$pollen_data,
                                        harmonization_level = harmonization_level) %>% 
    ungroup()
  
  if (write_to_flat_file) {
    write_csv(x = harm_pollen_data,
              path = paste(DIR_TAXA, paste(paste('harmonized_pollen_data', harmonization_level, 'harm_level', sep = '_'), DATAFILES_TYPE, sep = '.'), sep = '/'))
  }
  
  if (in_place) {
    db$harm_pollen_data <- harm_pollen_data
    return(db)
  }
  else {
    return(harm_pollen_data)
  }
}


arboreal_pollen_to_db <- function(x, ...) UseMethod('arboreal_pollen_to_db')

arboreal_pollen_to_db.PCNdata <- function(db, min_occurence_fraction = 0.01,
                                          in_place = TRUE, write_to_flat_file = FALSE) {
  if (!('harm_pollen_data' %in% names(db))) {
    stop('call harmonized_pollen_to_db before calculating arboreal pollen data match')
  }else {
    # percentized and core-merged data only for arboreal pollen_data (both called inside harmonized pollen_data)
    arb_pollen_data <- match_top_taxa_to_harmonized_arb_pollen_by_min_occ_frac(harm_pollen_data = db$harm_pollen_data,
                                                                               harmonization_level = 2,
                                                                               min_occurence_fraction = min_occurence_fraction) %>%
      ungroup()
    
    if (write_to_flat_file) {
      write_csv(x = arb_pollen_data,
                path = paste(DIR_ARB_POLLEN, paste(paste('arboreal_pollen_data', min_occurence_fraction, 'min_occ_frac', sep = '_'), DATAFILES_TYPE, sep = '.'), sep = '/'))
    }
    
    if (in_place) {
      db$arb_pollen_data <- arb_pollen_data
      return(db)
    }
    else {
      return(arb_pollen_data)
    }
  }
}


# PCNdata class methods for null model draws (called "hybrid ice cores / HIC", i.e. noised ice core data in v1) ----
hybrid_icecores_to_db <- function(x, ...) UseMethod('hybrid_icecores_to_db')

hybrid_icecores_to_db = function(db, in_place = TRUE, write_to_flat_file = FALSE, params) {
  mix_sample_dating_and_err()
  hic_data <- make_hybrid_icecores_for_ACERsites(db$sites, db$sample_dating, withnoise = FALSE, params = params) %>%
    ungroup()
  
  if (write_to_flat_file) {
    if ('sampling_mode' %in% names(params$sampling_params)) {name <- paste('hic_data', params$sampling_params$sampling_mode, sep = '_')}
    else {name <- 'hic_data'}
    write_csv(x = hic_data,
              path = paste(DIR_NM_DATA, paste(name, DATAFILES_TYPE, sep = '.'), sep = '/'))
  }
  
  if (in_place) {
    if ('sampling_mode' %in% names(params$sampling_params)) {name <- paste('hic_data', params$sampling_params$sampling_mode, sep = '_')}
    else {name <- 'hic_data'}
    db[[name]] <- hic_data
    return(db)
  }
  else {
    return(hic_data)
  }
  
}


set_hybrid_icecores_to_db <- function(x, ...) UseMethod('set_hybrid_icecores_to_db')

set_hybrid_icecores_to_db.PCNdata <- function(db, n = 100, write_trial = list(activate = TRUE, mode = 'single'), 
                                              write_noise = list(activate = TRUE, mode = 'single'), use_cached_noise = TRUE, 
                                     params = list(noise_params = list(snvr = 9,
                                                                       temporal = 'ar1',
                                                                       temporal_params = list(alpha = 0.8),
                                                                       spatial = 'matern',
                                                                       spatial_params = list(theta = 1000, smoothness = 1.5)),
                                                   sampling_params = list(sampling_mode = 'box',
                                                                          polamp = 4)), 
                                     unnest_data = FALSE, ncores = NULL) {
  if (!('mixed_age' %in% names(db$sample_dating))) {mix_dating_and_err()}
  data <- make_set_hybrid_icecores_for_ACERsites(n = n, samples = db$sample_dating, sites = db$sites, withnoise = TRUE,
                                                 write_trial = write_trial, write_noise = write_noise,
                                                 use_cached_noise = use_cached_noise, params = params, ncores = ncores)
  name <- paste('hic_set', paste(as.character(unlist(as.list(hash::hash(modifyList(list(n = n, nwrite_mode = write_trial$mode), params))))), collapse = '-'), sep = '_')
  db[[name]] <- data
  return(db)
}


# PCNdata class methods for transforming data (probit transform, detrending with a Gaussian filter) ----
corrected_data_to_db <- function(x, ...) UseMethod('corrected_data_to_db')

corrected_data_to_db.PCNdata <- function(db, orig_data = 'arb_pollen_data', transform = 'probit', detrend = list(activate = TRUE, method = 'linear', args = NULL), sd_one = FALSE,
                                unnest_data = FALSE, in_place, windows = GLOBAL_SERIES_WINDOWS, labels = GLOBAL_SERIES_WINDOWS_LABELS) {
  if(!(orig_data %in% names(db))) {stop('unknown orig_data passed to ACER_dataset$corrected_data_to_db, maybe need to call a ..._to_db member fct first')}
  if(!('mixed_age' %in% names(db$sample_dating))) {stop('need to call ACER_dataset$mix_dating before calling ACER_dataset$corrected_data_to_db')}
  def_detr <- list(activate = TRUE, method = 'linear', args = NULL)
  detrend <- merge.list(detrend, def_detr)
  
  orig_datas <- list(arb_pollen_data = 'arboreal_pollen', pollen_data = 'pollen', pca_harm_pollen_10_unscaled = 'pca_pollen', hic_data_box = 'hybrid_ice_core',
                     hic_data_point = 'hybrid_ice_core')
  
  # requires mixed_age
  data <- inner_join(db[[orig_data]], db$sample_dating, by = c('site_id', 'sample_id')) %>% 
    correct_data_stage(data = ., orig_data = orig_datas[[orig_data]],
                       transform = transform, detrend = detrend, sd_one = sd_one,
                       unnest_data = unnest_data, windows = windows, labels = labels,
                       signal_id = GLOBAL_SERIES_SIGNALS[[orig_datas[[orig_data]]]]$signal_id, signal_name = GLOBAL_SERIES_SIGNALS[[orig_datas[[orig_data]]]]$signal_name, 
                       write_to_flat_file = FALSE)
  
  if (in_place) {
    if(detrend$activate) {db[[paste('detrend', detrend$method, transform, orig_data, sep = '_')]] <- data
    } else {
      db[[paste(transform, orig_data, sep = '_')]] <- data
    }
    return(db)
    }
  else {return(data)}
  
}


# small wrapper, requires mixed_age
correct_data_stage <- function(data, orig_data, transform, detrend, sd_one, unnest_data, signal_id, signal_name, windows, labels, write_to_flat_file = FALSE, remove_age = TRUE) {
  data <- data %>% 
    filter(mixed_age >= 0) %>% 
    filter(!is.na(mixed_age)) %>% 
    select(UQ(sym(signal_id)), site_id, sample_id, mixed_age, UQ(sym(signal_name))) %>% 
    #filter(UQ(sym(signal_name)) <= 100) %>% 
    qtransform_data(data = ., type_data = orig_data, transform = transform, sd_one = if_else(detrend$activate == TRUE, true = FALSE, false = sd_one)) %>% 
    {if((detrend$activate == TRUE & unnest_data == TRUE) | write_to_flat_file) {window_and_detrend(data = ., type_data = orig_data, transform = transform, detrend = detrend, sd_one = sd_one, windows = windows, labels = labels) %>% 
        unnest(.) %>% select(-site_id) %>% unnest(.)}
      else if(detrend$activate == TRUE & unnest_data == FALSE) {window_and_detrend(data = ., type_data = orig_data, transform = transform, detrend = detrend, sd_one = sd_one, windows = windows, labels = labels)}
      else {.}} %>% 
    {if((detrend$activate == TRUE & unnest_data == TRUE | detrend$activate == FALSE) & remove_age) {select(., -mixed_age)}else {.}}
  return(data)
}


# PCNdata class methods for creating adjacency matrices based on the Gaussian kernel similarity measure and the null model draws as a reference for estimating significance ----
corr_list_from_set_to_db <- function(x, ...) UseMethod('corr_list_from_set_to_db')

corr_list_from_set_to_db.PCNdata <- function(db, orig_data = 'hic_set_0.8-1000-1.5-box-4-2-9-single', type_data = 'hybrid_ice_core',
                                    transform = 'probit', detrend = list(activate = TRUE, method = 'linear', args = NULL),
                                    interp = list(method = 'linear'), sd_one = TRUE, hres_only = TRUE,
                                    windows = GLOBAL_SERIES_WINDOWS, labels = GLOBAL_SERIES_WINDOWS_LABELS,
                                    site_ids = 'all', site_lats = 'all', in_place = TRUE, write_corr = list(activate = TRUE, dir = DIR_NM_DATA), ncores = NULL) {
  if(!(orig_data %in% names(db))) {stop('unknown orig_data passed to PCNdata$corrected_data_to_db, maybe need to call a ..._to_db member fct first')}
  if(!('mixed_age' %in% names(db$sample_dating))) {stop('need to call PCNdata$mix_dating before calling PCNdata$corrected_data_to_db')}
  
  data <- make_set_corr_list_from_dataset(data = db[[orig_data]], type_data = type_data, dating = db$sample_dating, sites_loc = db$sites, orig_data = orig_data,
                                          windows = windows, labels = labels, hres_only = hres_only, site_ids = site_ids, site_lats = site_lats, sd_one = sd_one, 
                                          transform = transform, detrend = detrend, interp = interp, write_corr = write_corr, ncores = ncores)
  if (in_place) {
    if(detrend$activate) {name <- paste('detrend', paste0(detrend$method, paste(detrend$args, collapse = '-')), transform, orig_data, sep = '_')}else {name <- paste(transform, orig_data, sep = '_')}
    if(hres_only) {name <- paste('hres', name, sep = '_')}
    if(interp$method != 'linear') {
      name <- paste(interp$method, name, sep = '_')
    } else if (str_detect(interp$scale, 'global')) {
      name <- paste(paste0(interp$method, interp$scale), name, sep = '_')
    }
    name <- paste('corr_list', name, sep = '_')
    db[[name]] <- data
    return(db)
    }
  else {return(data)}
}


corr_list_to_db <- function(x, ...) UseMethod('corr_list_to_db')

corr_list_to_db.PCNdata <- function(db, orig_data = 'arb_pollen_data', type_data = 'arboreal_pollen', transform = 'probit', detrend = list(activate = TRUE, method = 'linear', args = NULL),
                           interp = 'linear', gaussbp = list(activate = FALSE, per1 = 1250, per2 = 3750), sd_one = TRUE, hres_only = TRUE,
                           windows = GLOBAL_SERIES_WINDOWS, labels = GLOBAL_SERIES_WINDOWS_LABELS, dopar = TRUE,
                           site_ids = 'all', site_lats = 'all', in_place = TRUE, write_to_flat_file = FALSE) {
  if(!('mixed_age' %in% names(db$sample_dating))) {stop('need to call PCNdata$mix_dating before calling PCNdata$corrected_data_to_db')}
  if(str_detect(orig_data, 'pca') & type_data == 'pca_pollen') {
    if (str_detect(orig_data, '-pc\\d')) {pc_use <- paste0('pc', (str_replace(str_extract(orig_data, '-pc\\d'), '-pc', ''))); cat(paste('using pc', pc_use, '\n'))}else {stop('no pc component specified for correlation analysis')}
    if(!(str_replace(orig_data, '-pc\\d', '') %in% names(db))) {stop('unknown orig_data passed to PCNdata$corrected_data_to_db, maybe need to call a ..._to_db member fct first')}
    data <- db[[str_replace(orig_data, '-pc\\d', '')]] %>% filter(pc == pc_use) %>% select(-pc)
  }
  else {if(!(orig_data %in% names(db))) {stop('unknown orig_data passed to PCNdata$corrected_data_to_db, maybe need to call a ..._to_db member fct first')}; data <- db[[orig_data]]}
  
  data <- make_corr_list_from_dataset(data = data, type_data = type_data, dating = db$sample_dating, sites_loc = db$sites, windows = windows, labels = labels,
                                      hres_only = hres_only, site_ids = site_ids, site_lats = site_lats, sd_one = sd_one, transform = transform, detrend = detrend,
                                      interp = interp, gaussbp = gaussbp, dopar = dopar)
  
  if (write_to_flat_file) {
    if(detrend$activate) {
      name <- paste('detrend', paste0(detrend$method, paste(detrend$args, collapse = '-')), 
                    transform %>% unlist(use.names = T) %>% paste0(names(.), ., collapse = '-'), orig_data, sep = '_')
    } else {
      name <- paste(transform, orig_data, sep = '_')
    }
    if(hres_only) {name <- paste('hres', name, sep = '_')}
    name <- paste('corr_list', name, sep = '_')
    write_csv(x = data,
              path = paste(DIR_CORRECTED_DATA, name, DATAFILES_TYPE, sep = '.'), sep = '/')
  }
  
  if (in_place) {
    if(detrend$activate) {
      name <- paste('detrend', paste0(detrend$method, paste(detrend$args, collapse = '-')), 
                    transform %>% unlist(use.names = T) %>% paste0(names(.), ., collapse = '-'), orig_data, sep = '_')
    } else {
      name <- paste(transform %>% unlist(use.names = T) %>% paste0(names(.), ., collapse = '-'), orig_data, sep = '_')
    }
    if(hres_only) {name <- paste('hres', name, sep = '_')}
    if(interp$method != 'linear') {
      name <- paste(interp$method, name, sep = '_')
    } else if (str_detect(interp$scale, 'global')) {
      name <- paste(paste0(interp$method, interp$scale), name, sep = '_')
    }
    name <- paste('corr_list', name, sep = '_')
    db[[name]] <- data
    return(db)
    }
  else {return(data)}
}


corr_set_summary_to_db <- function(x, ...) UseMethod('corr_set_summary_to_db')

corr_set_summary_to_db.PCNdata <- function(db, orig_data = 'corr_list_hres_detrend_gaussbp_identity_hic_set_0.8-1000-1.5-box-4-2-9-single', 
                                  noise_data = NULL, test_option = 'rank', summary_option = 'single',
                                  write = list(activate = FALSE, dir = DIR_NM_DATA), in_place = TRUE) {
  if(!(orig_data %in% names(db))) {stop('unknown orig_data passed to PCNdata$corrected_data_to_db, maybe need to call a ..._to_db member fct first')}
  summs <- list(arb = calculate_summary_ap_corr_list, def = calculate_summary_set_corr_list)
  #if (str_detect(orig_data, 'arb_pollen')) {summ <- 'arb'}else if (str_detect(orig_data, 'pca')) {summ <- 'arb'}else {summ <- 'def'}
  if (all(class(db[[orig_data]]) == 'list')) {summ <- 'def'}else {summ <- 'arb'}
  if (!is.null(noise_data)) {
    data <- summs[[summ]](set_corr_list = db[[orig_data]], set_corr_list_noise = db[[noise_data]], test_option = test_option, write = write, summary_option = summary_option)
    orig_data <- paste0(test_option, 'test_', orig_data)
  }else {
    data <- summs[[summ]](set_corr_list = db[[orig_data]], write = write, summary_option = summary_option)
  }
  if (in_place) {
    if (summ != 'arb' & summary_option == 'single') {name <- paste('corr_set_summary-single', orig_data, sep = '_')} else {name <- paste('corr_set_summary', orig_data, sep = '_')}
    db[[name]] <- data
    return(db)
    }
  else {return(data)}
}

# PCNdata class methods for network plotting ----
run_network_plot <- function(x, ...) UseMethod('run_network_plot')

run_network_plot.PCNdata <- function(db, orig_data, projection = 'robinson', ndg_lims = NULL, #frac_abs_strongest
                                     black = FALSE, straight = FALSE, split_regions_by_sign = TRUE, leg_opt = 'single',
                                     bundled = TRUE, filter_windows = FALSE, node_degree = FALSE, filter_sites = NULL, bg = NULL,
                                     zoom = NULL, link_width = 0.4, color_strength = FALSE, save_plot = list(activate = FALSE), return = FALSE, return_nw = FALSE) {
  #for (f in frac_abs_strongest) {
    #for (b in bundled) {
  nw <- {if (!is.logical(filter_windows)) {corr_list <- filter(db[[orig_data]], window %in% filter_windows)} else {corr_list <- db[[orig_data]]}} %>% 
    create_network_from_corr_list(corr_list = .,
                                  sites = db$sites,
                                  filter_sites = filter_sites, split_regions_by_sign = split_regions_by_sign,
                                  split_window = F, get_node_occur = TRUE, region_tree = bundled, orig_corr_bundle = TRUE)
  #if (is.na(nw)) {break()}
  plot_corr_nw <- plot_network_spatial(nw, projection = projection, black = black, straight = straight, as_list = FALSE, bg = bg,
                                       split_regions_by_sign = split_regions_by_sign, ndg_lims = ndg_lims, leg_opt = leg_opt,
                                       zoom = zoom, bathy = FALSE, node_degree = node_degree, link_width = link_width, color_strength = color_strength)
  
  #plot_list <- list(corr_nw = plot_corr_nw)
    #}
  #}
  if (return & return_nw) {return(list('corr_nw' = plot_corr_nw, 'nw' = nw))} else if (return) {return(list('corr_nw' = plot_corr_nw))}
}
