# helpers to handle chronologies

mix_dating <- function(samples, orig_age_only = FALSE) {
  if (!orig_age_only) {
    samples <- mutate(samples, mixed_age = if_else(CLAM_best %in% c(-9999, -1999, NA), original_age, CLAM_best)) %>% 
      filter(!is.na(mixed_age))
  } else {
    samples <- mutate(samples, mixed_age = original_age) %>% 
      filter(!is.na(mixed_age))
  }
  return(samples)
}


mix_dating_and_err <- function(samples, orig_age_only = FALSE) {
  if('CLAM_best' %in% names(samples) & !orig_age_only) {
    samples <- mutate(samples,
                      mixed_age = if_else(CLAM_best %in% c(-9999, -1999, NA), original_age, CLAM_best),
                      mixed_age_min = CLAM_min95,
                      mixed_age_max = CLAM_max95) %>% 
      filter(!is.na(mixed_age)) %>% 
      mutate_at(vars(mixed_age), as.double)#%>% 
    #select(-original_age)#, -starts_with('CLAM'))
  } else if(orig_age_only) {
    samples <- mutate(samples,
                      mixed_age = original_age,
                      mixed_age_min = original_age,
                      mixed_age_max = original_age) %>% 
      filter(!is.na(mixed_age)) %>% 
      mutate_at(vars(mixed_age), as.double)
  }
  else {cat('dating already mixed')}
  return(samples)
}


merge_cores_in_data <- function(data, sites) {
  # data has to contain site_id, site_name, lat, long
  sites <- arrange(sites, lat, long) %>% 
    mutate(to_be_merged = if_else(((lead(lat) == lat & lead(long) == long) | (lag(lat) == lat & lag(long) == long)), TRUE, FALSE, missing = FALSE))
  sites_not_merge <- filter(sites, to_be_merged == FALSE) %>% 
    mutate(new_site_name = site_name,
           new_site_id = site_id)
  sites_merge <- filter(sites, to_be_merged == TRUE) %>% 
    group_by(lat, long) %>% 
    mutate(new_site_name = if_else(lead(lat) == lat & lead(long) == long, paste(site_name, lead(site_name), sep = ' / '), paste(lag(site_name), site_name, sep = ' / '),
                                   missing = paste(lag(site_name), site_name, sep = ' / ')), 
           new_site_id = as.double(if_else(lead(lat) == lat & lead(long) == long, paste(site_id, lead(site_id), sep = '.'), paste(lag(site_id), site_id, sep = '.'),
                                           missing = paste(lag(site_id), site_id, sep = '.'))))
  sites <- bind_rows(sites_merge, sites_not_merge) %>% 
    ungroup() %>% 
    arrange(site_id)
  
  data <- inner_join(select(sites, site_id, new_site_id, new_site_name), data, by = 'site_id') %>% 
    {if('site_name' %in% names(data)){select(., -site_id) %>% rename(site_id = new_site_id, site_name = new_site_name)}
      else{select(., -site_id, -new_site_name) %>% rename(site_id = new_site_id)}}
  
  sites <- select(sites, -site_id, -site_name, -to_be_merged) %>% 
    rename(site_id = new_site_id, site_name = new_site_name) %>% 
    select(site_id, site_name, lat:catch_size) %>% 
    distinct()
  
  return(list(data = data, new_sites = sites))
}


compute_age_intervals <- function(samples) {
  if (!('mixed_age' %in% names(samples))) {
    samples <- mix_dating(samples = samples) %>% 
      filter(mixed_age > 0)
  }
  
  samples <- group_by(samples, site_id) 
  
  if ('corrected_depth' %in% colnames(samples)) {
    samples <- arrange(samples, corrected_depth, .by_group = FALSE)
  } else {
    samples <- arrange(samples, mixed_age, .by_group = FALSE)
  }
  
  samples <- samples %>% 
    mutate(age_interv_lead = lead(mixed_age) - mixed_age)
  
  return(samples)
}


compute_age_interval_stats <- function(samples, ranges = SAMPLING_RATE_INTERVALS_DEFAULT) {
  samples <- ungroup(compute_age_intervals(samples)) %>% 
    select(sample_id, site_id, mixed_age, age_interv_lead) %>% 
    filter(!(is.na(age_interv_lead) | age_interv_lead < 0)) %>%
    arrange(site_id) %>% 
    group_by(site_id)
  
  stats <- tibble()
  for (range in ranges) {
    stats <- filter(samples, mixed_age >= range$start & mixed_age <= range$stop) %>% 
      summarise(interv = paste(range$start, range$stop, sep = '-'), 
                smp_in_interv = n(), 
                smp_per_10k_in_interv = n() * 10000 / (max(mixed_age) - min(mixed_age)), 
                mean_smp_res = mean(age_interv_lead),
                sd_smp_res = sd(age_interv_lead), 
                med_smp_res = median(age_interv_lead),
                mad_smp_res = mad(age_interv_lead), 
                ev_res_in_interv = (vegan::diversity(age_interv_lead) / log(vegan::specnumber(age_interv_lead))), 
                oldest_dp_in_interv = max(mixed_age), 
                youngest_dp_in_interv = min(mixed_age)) %>% 
      {if ('interv' %in% names(stats)) {
        bind_rows(., stats)
      } else {
        stats <- .
      }}
  }
  
  return(arrange(stats, site_id))
}


categorize_age_interval_stats <- function(age_interval_stats, 
                                          cats = list(
                                            'smp_per_10k_in_interv' = c(0, 25, 28, Inf),
                                            'mean_smp_res' = c(0, 200, 300, 500, Inf),
                                            #'mean_smp_res' = c(0, 200, 300, Inf),
                                            'med_smp_res' = c(0, 200, 300, 500, Inf),
                                            #'med_smp_res' = c(0, 200, 300, Inf),
                                            'ev_res_in_interv' = c(0, 0.6, 0.8, 0.9, Inf)
                                          )) {
  catslabs <- lapply(cats, 
                     function(x) {
                       x <- paste(x, c(x[2:length(x)], 'tab'), sep = ' - ') %>% .[!str_detect(., ' - tab')] %>% str_replace(' - Inf', '')
                       x[length(x)] <- paste0('> ', x[length(x)])
                       return(x)
                     })
  
  age_interval_stats <- mutate(age_interval_stats, 
                               status_smp_per_10k = cut(smp_per_10k_in_interv,
                                                        breaks = cats[['smp_per_10k_in_interv']],
                                                        labels = catslabs[['smp_per_10k_in_interv']],
                                                        ordered_result = TRUE), 
                               status_mean_smp_res = cut(mean_smp_res,
                                                         breaks = cats[['mean_smp_res']],
                                                         labels = catslabs[['mean_smp_res']],
                                                         ordered_result = TRUE), 
                               status_med_smp_res = cut(med_smp_res,
                                                        breaks = cats[['med_smp_res']],
                                                        labels = catslabs[['med_smp_res']],
                                                        ordered_result = TRUE),
                               status_ev_res_in_interv = cut(ev_res_in_interv,
                                                             breaks = cats[['ev_res_in_interv']],
                                                             labels = catslabs[['ev_res_in_interv']],
                                                             ordered_result = TRUE)) %>%  
    select(site_id, interv, oldest_dp_in_interv, youngest_dp_in_interv, status_smp_per_10k, status_mean_smp_res, status_med_smp_res, status_ev_res_in_interv)
  return(age_interval_stats)
}


compute_moving_average <- function(samples, window = 28, limit = TRUE) {
  samples <- ungroup(compute_age_intervals(samples)) %>% 
    select(sample_id, site_id, corrected_depth, mixed_age, age_interv_lead) %>% 
    filter(!is.na(age_interv_lead)) %>% 
    group_by(site_id) %>% 
    arrange(corrected_depth, .by_group = TRUE)
  
  stats <- mutate(samples %>% select(-corrected_depth), 
                  status_mean_smp_res = rollapply(age_interv_lead, FUN = mean, width = window, align = 'center', fill = 'extend'), 
                  status_med_smp_res = rollapply(age_interv_lead, FUN = median, width = window, align = 'center', fill = 'extend'), 
                  status_smp_per_10k = rollapply(mixed_age, FUN = function(z) length(z) * 10000 / abs(max(z) - min(z)), width = window, align = 'center', fill = 'extend'))
  if (is.logical(limit)) {
    if (limit) {
      stats <- mutate(stats, status_mean_smp_res = replace(status_mean_smp_res, status_mean_smp_res > 400, 400),
                      status_med_smp_res = replace(status_med_smp_res, status_med_smp_res > 400, 400), 
                      status_smp_per_10k = replace(status_smp_per_10k, status_smp_per_10k > 40, 40))
    }
  } else if (is.numeric(limit)) {
    stats <- mutate(stats, status_mean_smp_res = replace(status_mean_smp_res, status_mean_smp_res > limit, limit),
                    status_med_smp_res = replace(status_med_smp_res, status_med_smp_res > limit, limit), 
                    status_smp_per_10k = replace(status_smp_per_10k, status_smp_per_10k > limit/10, limit/10))
  }
  return(stats)
}

