# helper functions for file processing

load_db_flat <- function(filename, dir = DIR_DATASETS){
  path <- paste(paste(dir, filename, sep = '/'), DATAFILES_TYPE, sep = '.')
  print(paste('loading flat database file:', filename))
  return(as_tibble(read.csv(file = path, header = TRUE, sep = ',', encoding = 'UTF-8')))
}


print_db_flat_dir <- function(){
  print(list.files(path = DIR_DATASETS))
}


save_xtable <- function(x, filepath) {
  print.xtable(xtable(x), file = filepath, include.rownames = FALSE)
}


get_most_common_taxa <- function(pollen_dated, top_n = 10, taxon_signal){
  use_taxon <- pollen_dated %>% group_by(UQ(sym(taxon_signal))) %>% summarise(mean_taxon_pcnt = mean(taxon_pcnt, na.rm = TRUE)) %>% 
    .[.$mean_taxon_pcnt > 0,] %>% top_n(., n = top_n)
  use_taxon <- as.vector(use_taxon[[taxon_signal]], mode = 'character')
  
  print(paste('displaying the most common taxa which are'))
  print(use_taxon)
  
  return(list('pollen' = pollen_dated, 'use_taxon' = use_taxon))
}


load_forcing_dataset <- function() {
  if (!exists('forcings')) load(file.path(DIR_CACHE,'pre_processed_ts.RData'), envir = .GlobalEnv)
  return(forcings)
}


load_data_icecore <- function(core) { #, dir = file.path(DIR_DATASETS_SUPPLEMENTARY, 'icecores')) {
  cores <- c('epica', 'ngrip')
  if (!(core %in% cores)) {stop('unknown core given to load_data_icecore')}
  if (!exists('icecores')) load(file.path(DIR_CACHE,'pre_processed_ts.RData'), envir = .GlobalEnv)
  return(icecores[[core]])
}


get_core_existance <- function(data, core_types){
  l <- core_types
  for(x in names(l)){
    i <- l[[x]] %>% 
      distinct(site_id) %>%
      pull(site_id)
    
    u <- tibble(id = filter(data, site_id %in% i) %>% pull(site_id),
                ex = rep(TRUE, length(filter(data, site_id %in% i) %>% pull(site_id))))
    v <- tibble(id = filter(data, !site_id %in% i) %>% pull(site_id),
                ex = rep(FALSE, length(filter(data, !site_id %in% i) %>% pull(site_id))))
    
    w <- full_join(u, v, by = c('id', 'ex'))
    colnames(w) <- c('site_id', paste(x, '_exists', sep = ''))
    data <- left_join(data, w, by = 'site_id')
  }
  return(data)
}


write_dating_incoherences <- function(samples){
  break_indices <- vector('logical', nrow(samples))
  dating <- vector('logical', nrow(samples))
  
  for(i in 1:(nrow(samples)-1)){
    if(samples[i, 'CLAM_age_model_type'] != samples[i+1, 'CLAM_age_model_type'] 
       & samples[i, 'site_id'] == samples[i+1, 'site_id']){
      break_indices[[i]] <- TRUE
    }
    else{
      break_indices[[i]] <- FALSE
    }
    if(as.character(samples[i, 'CLAM_age_model_type']) == '-9999'){
      dating[[i]] <- FALSE
    }
    else{
      dating[[i]] <- TRUE
    }
  }
  samples$is_dating_break <- break_indices
  samples$uses_CLAM_dating <- dating
  return(samples)
}


time_restrict_data <- function(data, time_restrict){
  if(!is.null(time_restrict)){
    data <- data[data$mixed_age <= time_restrict$lower & data$mixed_age > time_restrict$upper,]
  }
  return(data)
}


sites_by_lat <- function(sites){
  dist_site_names <- distinct(sites, site_name, .keep_all = TRUE)
  sites$site_name <- factor(sites$site_name, levels = dist_site_names$site_name[order(-dist_site_names$lat)], ordered = TRUE)
  #sites$site_id <- factor(sites$site_id, levels = sites$site_id[order(-sites$lat)])
  return(sites)
}


filter_pollen_for_taxa_sites <- function(pollen, site_names, use_taxon, most_common_taxa, taxon_signal){
  if(site_names == 'all'){
    if(use_taxon == 'all' & most_common_taxa == 'all'){
      print('plotting all taxa from all sites')}
    else{
      pollen <- mutate_if(pollen, is.factor, as.character) %>% filter(UQ(sym(taxon_signal)) %in% use_taxon)}}
  else{
    site_names <- as.character(site_names)
    if(use_taxon == 'all'){
      pollen <- pollen[pollen$site_name %in% site_names,]}
    else{
      pollen <- mutate_if(pollen, is.factor, as.character) %>% filter(UQ(sym(taxon_signal)) %in% use_taxon)  %>% filter(as.character(site_name) %in% site_names)}
  }
  return(pollen)
}


filter_biomes_for_sites <- function(biomes, site_names){
  if(site_names != 'all'){
    biomes <- biomes[biomes$site_name %in% site_names,]
  }
  return(biomes)
}


lat_restrict_sites <- function(sites, lat = list('max' = 20, 'min' = -20)){
  sites <- sites[sites$lat <= lat$max & sites$lat >= lat$min,]
  return(as.character(sites$site_name))
}


signal.orig.ice.cores <- function() {
  # 101 NGRIP, 102 EPICA
  data_oc <- bind_rows(mutate(sample_hic_sign(lat = 90), site_id = 101), mutate(sample_hic_sign(lat = -90), site_id = 102)) %>% 
    rename(mixed_age = age) %>% 
    rownames_to_column(var = 'sample_id')
  return(data_oc)
}
