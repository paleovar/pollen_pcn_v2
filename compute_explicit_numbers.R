# Wrapper script for computing/extracting statistics/numbers which are referred to in the 
# text of the manuscript itself.
# - explicit network numbers for results & discussion
# - site-wise stats for table 1 in Supplement
# - site-wise correlation between emulated TSF of HadCM3 & TraCE and modeled TSF from the two models


# load/compute data
source('main_pcns.R')
source('main_pcn_meas.R')
source('main_pca.R')
source('main_emulation_data.R')

# additionally need to initialize some objects needed here in `create_fig4.R`
## `lds`
## `trc_suff`, `lc_suff`, `bbctrc_comb_suff`

# link density overview ----
lds_sum <- lds %>% group_by(tscale, signal, data) %>% summarise(ld_sum = sum(abs(ld)))
inner_join(lds_sum, lds) %>% mutate(ld_sign_perc = ld/ld_sum) %>% View()

lds_sum_em <- lds_em %>% group_by(tscale, signal, data) %>% summarise(ld_sum = sum(abs(ld)))
inner_join(lds_sum_em, lds_em) %>% mutate(ld_sign_perc = ld/ld_sum) %>% View()

# raw ACER node degrees ----
## bit ugly data extraction helper for CLF & NDG derived from plot function for bubble matrices
get_data_nw_clf_and_ndg_matr2x2_wosig <- function(datclf_raw, datclf_mil, datndg_raw, datndg_mil, suffix_pairs,
                                                  nw_meas_clf = 'pcross_link',
                                                  nw_meas_ndg = 'node_degree',
                                                  tscales = c('raw', 'mil'),
                                                  regions = REGIONS) {
  
  plts <- lapply(tscales, function(tscale) {
    datclf <- get(paste0('datclf_',tscale))
    datndg <- get(paste0('datndg_',tscale))
    plts <- lapply(1:length(suffix_pairs[[tscale]]), function(i){
      suffixes <- suffix_pairs[[tscale]][[i]]
      print(suffixes)
      # prepare clf data
      datclf <- datclf %>% 
        filter(cross_reg1 %in% unique(regions$region) & cross_reg2 %in% unique(regions$region)) %>% 
        mutate(cross_reg1 = substr(cross_reg1,1,2),cross_reg2 = substr(cross_reg2,1,2))
      datclf <- datclf %>% select(window, cross_reg1, cross_reg2, 
                                  !!!syms(mapply(function(nm, pr) return(paste0(pr, nm)), pr = c(paste0(nw_meas_clf, c('.p','.n'), '.')), nm = sort(rep(x = suffixes, 2)), USE.NAMES = FALSE)))
      datclf <- lapply(c(paste0(nw_meas_clf, c('.p','.n'))), function(c) {
        return(datclf %>% select(window, cross_reg1, cross_reg2, one_of(paste0(c, '.', suffixes))) %>% gather(key = 'type', value = !!sym(paste0(c, '_val')), -window, -cross_reg1, -cross_reg2) %>%
                 rowwise() %>% mutate(type = str_extract(type, suffixes) %>% .[!is.na(.)]) %>% ungroup())
      }) %>% plyr::join_all(., by = c('window', 'cross_reg1', 'cross_reg2', 'type'), type = 'inner') %>% as_tibble() %>% mutate_at(vars(window), as.character) %>% 
        rowwise %>% 
        rename(p = paste0(nw_meas_clf, '.p_val'), n = paste0(nw_meas_clf, '.n_val')) %>% 
        mutate(type = str_replace(type, '\\.', ' ')) %>% ungroup()
      datclf$type <- factor(datclf$type, levels = str_replace(suffixes, '\\.', ' '))
      
      regs_abv <- tibble(reg = distinct(bind_rows(select(datclf,cross_reg1),select(datclf,cross_reg2) %>% rename(cross_reg1 = cross_reg2)))$cross_reg1) %>% 
        rownames_to_column(var = 'reg_n') %>% 
        mutate_at(vars(reg_n), as.numeric)
      
      scaler_clf <- 0.5
      diag_start <- 0.5
      diag_end <- length(unique(c(as.vector(datclf$cross_reg1), as.vector(datclf$cross_reg2)))) + 0.5
      
      datclf <- datclf %>% 
        mutate_at(vars(cross_reg1,cross_reg2), as.character) %>% 
        inner_join(regs_abv %>% rename(cross_reg1 = reg, reg1_n = reg_n)) %>% 
        inner_join(regs_abv %>% rename(cross_reg2 = reg, reg2_n = reg_n)) %>% 
        mutate(sum = p+n, sum_scaled = sum*scaler_clf) %>% # mapping to srqt(sum) done in stat = 'pie' (consider legend below)
        filter(sum != 0) %>% 
        select(type,reg1_n,reg2_n,p,n,sum,sum_scaled)
      
      datclf <- datclf %>% 
        group_by(type) %>% 
        mutate(reg1_n_sort = case_when(type == suffixes[1] & reg1_n > reg2_n ~ reg2_n,
                                       type == suffixes[1] & reg1_n < reg2_n ~ reg1_n,
                                       type == suffixes[2] & reg1_n > reg2_n ~ reg1_n,
                                       type == suffixes[2] & reg1_n < reg2_n ~ reg2_n,
                                       TRUE ~ reg1_n),
               reg2_n_sort = case_when(type == suffixes[1] & reg1_n > reg2_n ~ reg1_n,
                                       type == suffixes[1] & reg1_n < reg2_n ~ reg2_n,
                                       type == suffixes[2] & reg1_n > reg2_n ~ reg2_n,
                                       type == suffixes[2] & reg1_n < reg2_n ~ reg1_n,
                                       TRUE ~ reg2_n)) %>% 
        rename('+' = p, '-' = n) %>% 
        inner_join(regs_abv %>% rename(region1 = reg, reg1_n = reg_n), by = 'reg1_n') %>%
        inner_join(regs_abv %>% rename(region2 = reg, reg2_n = reg_n), by = 'reg2_n')
      
      # prepare ndg data
      datndg <- datndg %>% select(window, site_id,
                                  !!!syms(mapply(function(nm, pr) return(paste0(pr, nm)), pr = c(paste0(nw_meas_ndg, c('','.p','.n'), '.')), nm = sort(rep(x = suffixes, 3)), USE.NAMES = FALSE)))
      datndg <- lapply(c(paste0(nw_meas_ndg, c('', '.p','.n'))), function(c) {
        return(datndg %>% select(window, site_id, one_of(paste0(c, '.', suffixes))) %>% gather(key = 'type', value = !!sym(paste0(c)), -window, -site_id) %>%
                 rowwise() %>% mutate(type = str_extract(type, suffixes) %>% .[!is.na(.)]) %>% ungroup())
      }) %>% plyr::join_all(., by = c('window', 'site_id', 'type'), type = 'inner') %>% as_tibble() %>% mutate_at(vars(window), as.character) %>% 
        rowwise %>%
        gather(key = 'ndg_type', value = 'ndg_val', !!!syms(paste0(paste0(nw_meas_ndg, c('', '.p','.n'))))) %>% 
        rowwise %>%
        mutate(ndg_type = case_when(str_detect(ndg_type, '.p') ~ 'p', str_detect(ndg_type, '.n') ~ 'n', TRUE ~ 'sum')) %>% 
        #rename(p = paste0(nw_meas_ndg, '.p_val'), n = paste0(nw_meas_ndg, '.n_val')) %>% 
        filter(type %in% suffixes)
      
      scaler_ndg <- 0.05
      
      datndg <- datndg %>% 
        inner_join(regions, by = 'site_id') %>% 
        mutate(region = substr(region,1,2)) %>% 
        spread(ndg_type,ndg_val) %>% 
        group_by(region, type, window) %>% 
        summarise(sum = mean(sum, na.rm = T), p = mean(p, na.rm = T), n = mean(n, na.rm = T)) %>% 
        ungroup() %>% 
        inner_join(regs_abv %>% rename(region = reg), by = 'region') %>% 
        mutate(sum_scaled = sqrt(sum)*scaler_ndg) 
      numndg <- datndg %>% 
        #filter(ndg_type == 'node_degree') %>% 
        group_by(region) %>% 
        count() %>% mutate_at(vars(n), as.character)
      datndg <- datndg %>% 
        mutate(frac_p = case_when(p != 0 & n != 0 ~ p/sum,
                                  p == 0 & n != 0 ~ 0,
                                  p != 0 & n == 0 ~ 1,
                                  T ~ 0.5),
               frac_n = case_when(p != 0 & n != 0 ~ n/sum,
                                  p == 0 & n != 0 ~ 1,
                                  p != 0 & n == 0 ~ 0,
                                  T ~ 0.5)) %>% 
        gather(key = 'sgn', value = 'single',p,n) %>% 
        mutate(start = case_when(type == suffixes[1] & sgn == 'p' ~ 1/4*pi - frac_p*pi,
                                 type == suffixes[2] & sgn == 'p' ~ 1/4*pi,
                                 type == suffixes[1] & sgn == 'n' ~ -3/4*pi,
                                 type == suffixes[2] & sgn == 'n' ~ 5/4*pi - frac_n*pi),
               end = case_when(type == suffixes[1] & sgn == 'p' ~ 1/4*pi,
                               type == suffixes[2] & sgn == 'p' ~ 1/4*pi + frac_p*pi,
                               type == suffixes[1] & sgn == 'n' ~ -3/4*pi + frac_n*pi,
                               type == suffixes[2] & sgn == 'n' ~ 5/4*pi)) %>% 
        mutate(type = str_replace_all(type, '_', ' ')) %>% ungroup() %>% #mutate(shapetype = str_replace(str_extract(type, '\\w+\\s'), ' ', '')) %>% ungroup()
        mutate(sgn = case_when(sgn == 'n' ~ '-',
                               sgn == 'p' ~ '+', 
                               TRUE ~ sgn))
      
      datndg$type <- factor(datndg$type, levels = str_replace_all(suffixes, '_', ' '))
      
      return(list(datndg=datndg,datclf=datclf))
    })
    out <- list()
    out[[tscale]] <- plts
    return(out)
  })
  return(plts)
}


datACERtrc <- get_data_nw_clf_and_ndg_matr2x2_wosig(ACERtrc_clp_raw %>% rename_if(., str_detect(names(.), 'TRACE_apsb'), ~ str_replace(., 'TRACE_apsb', 'TRACE_tsf')) %>% 
                                                      rename_if(., str_detect(names(.), 'TRACE'), ~ str_replace(., 'TRACE', 'TraCE')),
                                                    ACERtrc_clp_mil %>% rename_if(., str_detect(names(.), 'TRACE_apsb'), ~ str_replace(., 'TRACE_apsb', 'TRACE_tsf')) %>% 
                                                      rename_if(., str_detect(names(.), 'TRACE'), ~ str_replace(., 'TRACE', 'TraCE')),
                                                    ACERtrc_ndg_raw %>% rename_if(., str_detect(names(.), 'TRACE_apsb'), ~ str_replace(., 'TRACE_apsb', 'TRACE_tsf')) %>% 
                                                      rename_if(., str_detect(names(.), 'TRACE'), ~ str_replace(., 'TRACE', 'TraCE')),
                                                    ACERtrc_ndg_mil %>% filter(site_id <= 100) %>% rename_if(., str_detect(names(.), 'TRACE_apsb'), ~ str_replace(., 'TRACE_apsb', 'TRACE_tsf')) %>% 
                                                      rename_if(., str_detect(names(.), 'TRACE'), ~ str_replace(., 'TRACE', 'TraCE')),
                                                    suffix_pairs = trc_suff[c(1,2)])
datACERlc <- get_data_nw_clf_and_ndg_matr2x2_wosig(ACERlc_clp_raw %>% rename_if(., str_detect(names(.), 'LC_'), ~ str_replace(., 'LC_', 'LOVECLIM_')),
                                                    ACERlc_clp_mil %>% rename_if(., str_detect(names(.), 'LC_'), ~ str_replace(., 'LC_', 'LOVECLIM_')),
                                                    ACERlc_ndg_raw %>% rename_if(., str_detect(names(.), 'LC_'), ~ str_replace(., 'LC_', 'LOVECLIM_')),
                                                    ACERlc_ndg_mil %>% filter(site_id <= 100) %>% rename_if(., str_detect(names(.), 'LC_'), ~ str_replace(., 'LC_', 'LOVECLIM_')),
                                                    suffix_pairs = lc_suff)
datACERhcm <- get_data_nw_clf_and_ndg_matr2x2_wosig(ACERhcm_clp_raw %>% 
                                                     rename_if(., str_detect(names(.), 'BBC_ap'), ~ str_replace(., 'BBC_ap', 'BBC_tsf')) %>% 
                                                     rename_if(., str_detect(names(.), 'BBC'), ~ str_replace(., 'BBC', 'HadCM3')),
                                                   NULL,
                                                   ACERhcm_ndg_raw %>% 
                                                     rename_if(., str_detect(names(.), 'BBC_ap'), ~ str_replace(., 'BBC_ap', 'BBC_tsf')) %>% 
                                                     rename_if(., str_detect(names(.), 'BBC'), ~ str_replace(., 'BBC', 'HadCM3')),
                                                   NULL,
                                                   tscales = 'raw',
                                                   suffix_pairs = bbctrc_comb_suff)

## highest lowest (overall)
ACERtrc_ndg_raw %>% select(contains('ACER_ap'), site_id) %>% summarise(min=min(node_degree.ACER_ap, na.rm = T),max=max(node_degree.ACER_ap, na.rm = T))

## summary ACER
datACERtrc[[1]]$raw[[1]]$datndg %>% filter(type == 'ACER ap')

## average number a record is linked to (N-1 = 57 records theoretically available after filtering)
datACERtrc[[1]]$raw[[1]]$datndg %>% filter(type == 'ACER ap') %>% select(region,sum) %>% mutate(ratio = sum/57)
## overall average excluding Asia
datACERtrc[[1]]$raw[[1]]$datndg %>% filter(type == 'ACER ap') %>% select(region,sum) %>% filter(region!="As") %>% mutate(ratio = sum/57) %>% summarise(all = mean(ratio))


## raw ACER cross-link fractions
datACERtrc[[1]]$raw[[1]]$datclf %>% filter(type == 'ACER_ap') %>% arrange(sum) %>% mutate(negative_perc = `-`/sum) %>% select(-reg1_n_sort, -reg2_n_sort)
## without Asia involved
datACERtrc[[1]]$raw[[1]]$datclf %>% filter(type == 'ACER_ap') %>% arrange(sum) %>% mutate(negative_perc = `-`/sum) %>% select(-reg1_n_sort, -reg2_n_sort) %>% 
  filter(region1 != 'As' & region2 != 'As')

## ndg difference ACER AP TRACE TSF, PR, TS raw
bind_rows(datACERtrc[[1]]$raw[[1]]$datndg %>% group_by(region) %>% select(region,type,sum,frac_p,frac_n) %>% distinct(),
          datACERtrc[[1]]$raw[[2]]$datndg %>% group_by(region) %>% select(region,type,sum,frac_p,frac_n) %>% distinct()) %>% 
  arrange(region) %>% View()
bind_rows(datACERtrc[[1]]$raw[[1]]$datndg %>% group_by(region) %>% select(region,type,sum,frac_p,frac_n) %>% distinct(),
          datACERtrc[[1]]$raw[[2]]$datndg %>% group_by(region) %>% select(region,type,sum,frac_p,frac_n) %>% distinct()) %>% 
  group_by(type) %>% 
  summarise(mean = mean(sum))

## clf statistics
## mean
bind_rows(datACERtrc[[1]]$raw[[1]]$datclf,
          datACERtrc[[1]]$raw[[2]]$datclf,
          datACERlc[[1]]$raw[[1]]$datclf %>% filter(type != 'ACER_ap'),
          datACERlc[[1]]$raw[[2]]$datclf,
          datACERhcm[[1]]$raw[[1]]$datclf,
          datACERhcm[[1]]$raw[[2]]$datclf) %>% 
  group_by(type) %>% select(type,region1,region2,sum) %>% distinct() %>% 
  summarise(mean = mean(sum))

## additional stats
bind_rows(datACERtrc[[1]]$raw[[1]]$datclf,
          datACERtrc[[1]]$raw[[2]]$datclf,
          datACERlc[[1]]$raw[[1]]$datclf %>% filter(type != 'ACER_ap'),
          datACERlc[[1]]$raw[[2]]$datclf,
          datACERhcm[[1]]$raw[[1]]$datclf,
          datACERhcm[[1]]$raw[[2]]$datclf) %>% 
  group_by(type) %>% select(type,region1,region2,sum,`+`,`-`) %>% distinct() %>% 
  mutate(negfrac = `-`/sum, posfrac = `+`/sum) %>% 
  summarise(mean_posfrac = mean(posfrac, na.rm = T), min_posfrac = min(posfrac, na.rm = T), maxposfrac = max(posfrac, na.rm = T), #median_posfrac = median(posfrac, na.rm = T),
            mean_negfrac = mean(negfrac, na.rm = T), min_negfrac = min(negfrac, na.rm = T), maxnegfrac = max(negfrac, na.rm = T), #median_negfrac = median(negfrac, na.rm = T),
            n_geq_0.33_negfrac = negfrac[which(negfrac >= 1/3)] %>% length(),
            n_geq_0.66_negfrac = negfrac[which(negfrac >= 2/3)] %>% length(),
            n_geq_0.33_posfrac = negfrac[which(posfrac >= 1/3)] %>% length(),
            n_geq_0.66_posfrac = negfrac[which(posfrac >= 2/3)] %>% length())

# average ndg pos vs neg TRACE TS raw 
datACERtrc[[1]]$raw[[2]]$datndg %>% filter(type == 'TraCE ts') %>% select(region,type,sum, frac_p, frac_n) %>% distinct()# %>% spread(type,sum)#%>% arrange(sum) %>% mutate(negative_perc = `-`/sum) %>% select(-reg1_n_sort, -reg2_n_sort)
# and clf
datACERtrc[[1]]$raw[[2]]$datclf %>% filter(type == 'TraCE_ts') %>% arrange(sum) %>% 
  mutate(negative_perc = `-`/sum, pos_perc = `+`/sum) %>% select(-reg1_n_sort, -reg2_n_sort)

# stats of the up to 63 used records as basis for Table 1 in Supplement ----
ACERhere <- PCNdata() %>% 
  mix_sample_dating_and_err() %>% 
  harmonized_pollen_to_db() %>% 
  arboreal_pollen_to_db()

ACERstats <- compute_age_interval_stats(samples = ACERhere$sample_dating, ranges = list(list(start = 6000, stop = 22000))) %>% 
  select(site_id, med_smp_res, oldest_dp_in_interv, youngest_dp_in_interv) %>% 
  rename(ISR = med_smp_res, start_age = oldest_dp_in_interv, end_age = youngest_dp_in_interv) %>% 
  rowwise() %>% 
  mutate(start_age = round(start_age*1e-3, 1),
         end_age = round(end_age*1e-3, 1)) %>% 
  inner_join(ACERhere$sites %>% select(site_id, site_name, lat, long))

sites_fltr_stats <- filter_sites(ACERhere$sites, ACERhere$sample_dating, hres_only = T, 'all', 'all') %>% 
  mutate(hres = TRUE) %>% 
  filter(site_id != 36) # filter Megali which has only one sample in interval and is filtered out by the correlation routine 

publications <- read_csv(file = file.path(DIR_DATASETS, 'site_publication.csv')) %>% 
  inner_join(., read_csv(file = file.path(DIR_DATASETS, 'publication.csv'))) %>% 
  mutate(site_id = case_when(site_id %in% c(85.,100.) ~ 85.10,
                             site_id %in% c(90.,98.) ~ 90.98,
                             site_id %in% c(72.,97.) ~ 72.97,
                             TRUE ~ site_id))

publications_t <- publications %>% 
  rowwise() %>% 
  mutate(year = str_extract_all(CITATION, '\\([0-9].*\\)') %>% str_extract(., '[0-9]{4}'),
         unpublished = str_detect(CITATION, 'unpublished'),
         author1_m1 = str_extract(CITATION, '^(.+?)\\(') %>% str_split('\\(') %>% .[[1]] %>% .[1] %>% str_split(';') %>% .[[1]] %>% .[1], #%>% str_extract('^(.+?);') %>% 
         author1_m2 = str_extract(CITATION, '^(.+?)unpublished') %>% str_split(';') %>% .[[1]] %>% .[1],#%>% str_extract('^(.+?):') %>% str_split(':') %>% .[[1]] %>% .[1],
         https = str_extract_all(CITATION, 'https.*'),# %>% str_split('https://') %>% .[[1]] %>% .[2],
         doi = str_extract_all(CITATION, 'doi.*')) %>%  #%>% str_split('doi:') %>% .[[1]] %>% .[2]) %>% 
  unnest(c(doi,https,author1_m1, author1_m2, unpublished), keep_empty = T) %>% 
  rowwise() %>% 
  mutate(author1 = case_when(is.na(author1_m1) ~ author1_m2,
                             TRUE ~ author1_m1),
         year = if_else(unpublished, "unpublished", year)) %>% 
  mutate(m2_is_single = if_else(str_extract(CITATION, '^(.+?)unpublished') %>% str_split(';') %>% .[[1]] %>% length() == 2, T, F),
         m1_is_single = str_extract(CITATION, '^(.+?)\\(') %>% str_split('\\(') %>% .[[1]] %>% .[1] %>% str_split(';') %>% .[[1]] %>% length() == 1,
         author_is_single = if_else(m2_is_single | m1_is_single, T, F)) %>% 
  mutate(author1 = if_else((str_split(author1, ', ')[[1]][2] %>% str_split(' ') %>% .[[1]])[str_split(author1, ', ')[[1]][2] %>% str_split(' ') %>% .[[1]] != ""] %>% length() > 1,
                           paste0(str_split(author1, ',')[[1]][1], ', ', str_split(author1, ', ')[[1]][2] %>% str_split(' ') %>% .[[1]] %>% .[1] %>% str_sub(end = 1), '. ', str_split(author1, ', ')[[1]][2] %>% str_split(' ') %>% .[[1]] %>% .[2] %>% str_sub(end = 1), '.'),
                           paste0(str_split(author1, ',')[[1]][1], ', ', str_split(author1, ', ')[[1]][2] %>% str_sub(end = 1), '.'))) %>% 
  mutate(citation = if_else(author_is_single, paste0(author1, ' (', year, ')'), paste0(author1, ' et al. (', year, ')')),
         web = if_else(is.na(doi), https, doi)) %>% 
  mutate(citation = if_else(is.na(web), citation, paste0(citation, ', ', web))) %>% 
  select(site_id, pub_id, citation) %>% 
  group_by(site_id) %>% 
  nest(.key = 'citation') %>% 
  mutate(citation = purrr::map(citation, function(x) tibble(citation = paste0(x$citation, collapse = '; ')))) %>% 
  unnest(citation)

ACER_stats_out <- full_join(ACERstats, publications_t) %>% 
  inner_join(sites_fltr_stats) %>% 
  rowwise() %>% 
  mutate(coverage = paste0(start_age, "-", end_age)) %>% 
  ungroup() %>% 
  inner_join(REGIONS %>% inner_join(tibble(region = c("Europe", "Africa", "Asia", "SAmerica", "NAmerica", "Australia"), 
                                           region_short = c("EU", "AF", "AS", "SA", "NA", "OC"))) %>% 
               select(site_id, region_short)) %>% 
  select(site_id, site_name, region_short, lat, long, coverage, ISR, citation) %>% 
  rename(ID = site_id, Name = site_name, Continent = region_short, !!sym("Lon [$\\degree$]") := long, !!sym("Lat [$\\degree$]") := lat,
         !!sym("Coverage in 6-22 ka interval [ka]") := coverage, References = citation)

print(xtable(ACER_stats_out, type = "latex"), file = file.path(DIR_FIGURES, "ACER_stats.tex"))

# nmb terrestrial and marine records
ACERhere$sites %>% inner_join(ACER_stats_out %>% rename(site_id = ID)) %>% group_by(site_type) %>% summarise(n = n())

# mean of the median ISR
ACER_stats_out %>% summarise(mean_of_median_ISR = mean(ISR))

rm(ACERhere)

# site-wise correlation between emulated and modeled TSF ----
## TraCE
site_corr_em_mod <- bind_rows(ACERtrc_em$forceall$arb_pollen_data %>% 
            inner_join(ACERtrc_em$forceall$sample_dating) %>% 
            mutate(type = 'EM_TSF'),
          ACERtrc$TRACE_apsb$arb_pollen_data %>% 
            inner_join(ACERtrc$TRACE_apsb$sample_dating) %>% 
            filter(mixed_age > WINDOWS_TRACE[[1]] & mixed_age < WINDOWS_TRACE[[2]]) %>% # emulated forcing is filtered already
            mutate(type = 'TSF')) %>% 
  select(type, site_id, sample_id, mixed_age, pcnt_arb_pollen) %>% 
  filter(site_id %in% sites_in_nws) %>% 
  group_by(type,site_id) %>% 
  nest() %>% 
  mutate(data_zoo = purrr::map(data, function(x) {
    x <- x %>% arrange(mixed_age)
    return(zoo(x$pcnt_arb_pollen,order.by = x$mixed_age))
      })) %>% 
  select(-data) %>% 
  spread(type,data_zoo)
site_corr_em_mod_c <- site_corr_em_mod %>% 
  mutate(corr = purrr::map2(EM_TSF, TSF, function(x,y) cor(coredata(x),coredata(y))),
         corr_gkc = purrr::map2(EM_TSF, TSF, function(x,y) nexcf(x,y)[1])) %>% 
  unnest(c(corr,corr_gkc)) %>% 
  ungroup() %>% 
  summarise(mean_corr = mean(corr),mean_corr_gkc = mean(corr_gkc)) %>% 
  print()

## HadCM3
site_corr_em_mod <- bind_rows(ACERhcm_em$forceall$arb_pollen_data %>% 
                                inner_join(ACERhcm_em$forceall$sample_dating) %>% 
                                mutate(type = 'EM_TSF'),
                              ACERhcm$BBC_ap$arb_pollen_data %>% 
                                inner_join(ACERhcm$BBC_ap$sample_dating) %>% 
                                filter(mixed_age > WINDOWS_BBC[[1]] & mixed_age < WINDOWS_BBC[[2]]) %>% # emulated forcing is filtered already
                                mutate(type = 'TSF')) %>% 
  select(type, site_id, sample_id, mixed_age, pcnt_arb_pollen) %>% 
  filter(site_id %in% sites_in_nws) %>% 
  group_by(type,site_id) %>% 
  nest() %>% 
  mutate(data_zoo = purrr::map(data, function(x) {
    x <- x %>% arrange(mixed_age)
    return(zoo(x$pcnt_arb_pollen,order.by = x$mixed_age))
  })) %>% 
  select(-data) %>% 
  spread(type,data_zoo)
site_corr_em_mod_c <- site_corr_em_mod %>% 
  mutate(corr = purrr::map2(EM_TSF, TSF, function(x,y) cor(coredata(x),coredata(y))),
         corr_gkc = purrr::map2(EM_TSF, TSF, function(x,y) nexcf(x,y)[1])) %>% 
  unnest(c(corr,corr_gkc)) %>% 
  ungroup() %>% 
  summarise(mean_corr = mean(corr),mean_corr_gkc = mean(corr_gkc)) %>% 
  print()

