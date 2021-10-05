# helper functions for signal processing

# transformations (publication only uses probit transform for ACER AP)
logit <- function(p) {
  return(log(p / (1 - p)))
}


id_trf <- function(x) {
  return(x)
}


qtransform_data <-
  function(data, type_data, transform, sd_one = FALSE) {
    # can be used for both single site or multi-site tbl
    # transforms to zero mean always but to sd = 1 only if asked to; latter option should only be chosen if no further detrending is performed, otherwise chose option after detrending!
    # removes Inf values which are eventually introduced
    transforms = list(probit = qnorm,
                      logit = logit,
                      identity = id_trf,
                      sqrt = sqrt)
    if (!(transform %in% names(transforms))) {
      stop('unknown transform given to qtransform_data')
    }
    data <- group_by(data, site_id) %>%
      {
        if (transform == 'identity') {
          mutate(., UQ(sym(
            paste(transform, GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name, sep = '_')
          )) := transforms[[transform]](UQ(
            sym(GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name)
          )))
        }
        else {
          filter(., UQ(sym(
            GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name
          )) <= 100) %>%
            mutate(., UQ(sym(
              paste(transform, GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name, sep = '_')
            )) := transforms[[transform]](UQ(
              sym(GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name)
            ) / 100))
        }
      } %>%
      filter(!is.infinite(UQ(sym(
        paste(transform, GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name, sep = '_')
      )))) %>%
      filter(!is.na(UQ(sym(
        paste(transform, GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name, sep = '_')
      )))) %>%
      {
        if (sd_one) {
          mutate(., normalized = UQ(sym(
            paste(transform, GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name, sep = '_')
          )) / sd(UQ(sym(
            paste(transform, GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name, sep = '_')
          )))) %>%
            select(-UQ(sym(
              paste(transform, GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name, sep = '_')
            ))) %>%
            rename(UQ(sym(
              paste(transform, GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name, sep = '_')
            )) := normalized)
        }
        else {
          .
        }
      }
    if (sd_one) {
      warning('mapping to sd = 1 discarded all Inf values')
      warning('use sd_one = TRUE in qtransfrom_data only if data will not be detrended afterwards')
    }
    return(data)
  }



detrend_data_window <-
  function(data_windowed,
           type_data,
           transform,
           sd_one,
           method = 'linear',
           method_args = NULL) {
    # to be used on a single window of a single site
    methods <-
      list(linear = detrend_series_linear, gaussbp = detrend_series_gaussbp)
    if (!(method %in% names(methods))) {
      stop('unknown method provided to detrend_data_window or above layered function')
    }
    data_windowed <- mutate(
      data_windowed,
      UQ(sym(
        paste(
          'detrend',
          transform,
          GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name,
          sep = '_'
        )
      )) :=
        methods[[method]](
          series = select(data_windowed, mixed_age, UQ(sym(
            paste(transform, GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name, sep = '_')
          ))),
          type_data = type_data,
          transform = transform,
          sd_one = sd_one,
          args = method_args
        )
    )
    return(data_windowed)
  }


# detrending: in publication only Gaussian filtering is used for temporal detrending
# linear detrending is supported by the code whatsoever, without any warranty
detrend_series_linear <- function(series, type_data, transform, sd_one, args) {
  # to be used for on single time series (windowed)
  
  fit <- series %>% 
    arrange(mixed_age) %>% 
    do(regr = lm(UQ(sym(paste(transform, GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name, sep = '_'))) ~ mixed_age, data = .)) %>% 
    tidy(x = ., regr)
  
  series <- mutate(series,
                   detrend = UQ(sym(paste(transform, GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name, sep = '_'))) - (fit$estimate[2] * mixed_age + fit$estimate[1])) %>% 
    {if(sd_one) {mutate(., detrend_one = detrend / sd(detrend)) %>% select(-detrend) %>% rename(detrend = detrend_one)}else {.}}
  
  return(series$detrend)
  #return(series[[paste(transform, GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name, sep = '_')]])
}


detrend_series_gaussbp <-
  function(series, type_data, transform, sd_one, args) {
    if (is.null(args)) {
      stop('no valid method_args provided to function detrend_series_gaussbp')
    }
    series <- series %>%
      arrange(mixed_age) %>%
      filter(!is.infinite(UQ(sym(
        paste(transform, GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name, sep = '_')
      ))))
    
    xsc <-
      gaussbandpass(
        zoo(as.data.frame(series[[paste(transform, GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name, sep = '_')]]), order.by = series$mixed_age),
        per1 = args$per1,
        per2 = args$per2
      )$filt
    if (length(xsc) != length(series[[paste(transform, GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name, sep = '_')]])) {
      stop(paste('no values from gaussbandpass'))
    }
    else{
      series$detrend <- xsc
    }
    series <-
      {
        if (sd_one) {
          mutate(series, detrend_one = detrend / sd(detrend)) %>% select(-detrend) %>% rename(detrend = detrend_one)
        } else {
          series
        }
      }
    return(series$detrend)
  }


window_data <-
  function(data,
           type_data,
           transform = NULL,
           windows = GLOBAL_SERIES_WINDOWS,
           labels = GLOBAL_SERIES_WINDOWS_LABELS) {
    data <- filter(data, mixed_age >= 0) %>%
      {
        if (!is.null(transform)) {
          filter(.,!is.infinite(UQ(sym(
            paste(transform, GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name, sep = '_')
          ))))
        }
        else {
          filter(.,!is.infinite(UQ(
            sym(GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name)
          )))
        }
      } %>%
      group_by(site_id) %>%
      arrange(mixed_age) %>%
      mutate(window = cut(mixed_age, breaks = windows, labels = labels)) %>%
      filter(!is.na(window)) %>%
      group_by(window)
    return(data)
  }


window_and_detrend <-
  function(data,
           windows,
           labels,
           type_data,
           transform,
           detrend,
           sd_one = FALSE) {
    # cut time windows and eventually detrend
    # linear detrending will be carried out on cut windows, gaussbp detrending will be carried out on entire time series with window cutting afterwards
    def_detrend <-
      list(activate = TRUE,
           method = 'linear',
           args = NULL)
    detrend <-  RCurl::merge.list(detrend, def_detrend)
    
    if (detrend$activate) {
      if (detrend$method == 'linear') {
        data <-
          window_data(
            data = data,
            type_data = type_data,
            transform = transform,
            windows = windows,
            labels = labels
          )
        data <- data %>%
          group_by(window) %>%
          do(
            data_win = group_by(., site_id) %>%
              select(-window) %>%
              do(
                site_data := detrend_data_window(
                  data_windowed = .,
                  type_data = type_data,
                  transform = transform,
                  method = detrend$method,
                  method_args = detrend$args,
                  sd_one = sd_one
                )
              )
          )
      }
      else if (detrend$method == 'gaussbp') {
        data <- data %>%
          group_by(site_id) %>%
          do(
            site_data := detrend_data_window(
              data_windowed = select(.,-site_id),
              type_data = type_data,
              transform = transform,
              method = detrend$method,
              method_args = detrend$args,
              sd_one = sd_one
            )
          ) %>%
          ungroup() %>%
          unnest() %>%
          group_by(site_id) %>%
          window_data(
            data = .,
            type_data = type_data,
            transform = transform,
            windows = windows,
            labels = labels
          ) %>%
          mutate(site_id2 = site_id) %>%
          group_by(window) %>%
          nest(.key = 'data_win') %>%
          mutate(data_win = purrr::map(
            data_win,
            ~ .x %>% group_by(site_id2) %>% nest(.key = 'site_data') %>% rename(site_id = site_id2)
          ))
      }
    }
    else{
      data <-
        window_data(
          data = data,
          type_data = type_data,
          transform = transform,
          windows = windows,
          labels = labels
        )
      data <- data %>%
        group_by(window) %>%
        do(data_win = group_by(., site_id) %>%
             select(-window) %>%
             do(site_data := .)) # eventually select(., -site_id)
    }
    return(data)
  }


# retrieval of time series from model data

#' Title
#'
#' @param field 
#' @param center_ixs 
#' @param match_sites 
#' @param time 
#'
#' @return
#' @export
#'
#' @examples
retrieve_tseries <- function(field, center_ixs, 
                             match_sites,
                             time) {
  tseries <- lapply(1:length(match_sites), function(i) {
    field[center_ixs[i,1],center_ixs[i,2],] 
  }) %>%
    setNames(match_sites) %>% 
    bind_cols(tibble(mixed_age = time))
  return(tseries)
}


#' Retrieve pseudo proxy time series from field
#'
#' @param field matrix. expects dim(lon,lat,time)
#' @param match_locs 
#' @param match_sites 
#' @param grid 
#' @param lons 
#' @param lats 
#' @param time 
#' @param search_depth integer. `search_depth = 5` searches up to the von Neumann neighborhood, `search_depth = 9` up to the Moore neighborhood of the closest matching cell
#' @param overwrite_sites integer. indices of match_sites which should be forced to the first terieved time series
#'
#' @return list('tseries' = set of retrieved time series, 'center_ixs_used' = used indices of grid cells, 'iterations' = summary of the retrieval process)
#' @export
#'
#' @examples
retrieve_nearest_tseries <- function(field, match_locs, 
                                     match_sites,
                                     grid, lons, 
                                     lats, time, 
                                     #avg_surrounding = FALSE, 
                                     search_depth = 5,
                                     overwrite_sites = NULL) {
  #if (avg_surrounding == TRUE & search_depth > 1) stop('search_depth > 5 and avg_surrounding == TRUE; decide for one retrieval scheme!')
  
  dist2field <- raster::pointDistance(match_locs,grid,lonlat = T)
  ixs <- apply(dist2field, 1, function(x) sort(x,index.return = TRUE)$ix[1:search_depth])
  #print(dim(dist2field)); print(dim(ixs)); dim(field)
  if (search_depth == 1) {ixs <- matrix(ixs, nrow = 1)}
  tmp_list <- lapply(1:length(match_sites), function(i) {
    coordsl <- coordinates(grid[ixs[,i]])
    #print(i)
    #print(ixs[,i])
    #print(coordsl)
    #ijsl <- tibble(lon = match(coordsl[,1], lons), lat = match(coordsl[,2], lats))
    ijsl <- matrix(c(match(coordsl[,1], lons), match(coordsl[,2], lats)),ncol = 2)
    #print(ijsl)
    tseriesl <- lapply(1:search_depth, function(l) {
      field[ijsl[l,1],ijsl[l,2],]
    })
    #print(tseriesl)
    selix <- lapply(tseriesl, function(x) if_else(i %in% overwrite_sites, TRUE, all(!is.na(x)))) %>% # any(!is.na(x)))) %>% 
      unlist() %>% 
      which(. == TRUE) %>% 
      .[1]
    #if (is.na(selix)) stop(paste0('search_depth = ',search_depth,' insufficient'))
    return(list('ijs' = ijsl[selix,], 'tseries' = tseriesl[[selix]], 'selix' = selix))
  })
  
  tseries <- lapply(tmp_list, function(x) return(x$tseries))%>% # print()  %>% 
    setNames(match_sites) %>% 
    bind_cols(tibble(mixed_age = time)) 
  ixs_used <- lapply(tmp_list, function(x) return(x$ijs)) %>% 
    unlist() %>% 
    matrix(., ncol = 2, byrow = TRUE)
  iterations <- lapply(tmp_list, function(x) return(x$selix)) %>% 
    unlist()
  iterations <- tibble(site_id = match_sites, n_iterations = iterations)
  
  
  return(list('tseries' = tseries, 'center_ixs_used' = ixs_used, 'iterations' = iterations))
}


### helper functions ----
spatial_means <- function(lon,lat,clim_field) {
  spatial_weights <- array(0,dim=c(length(lon),length(lat)))
  for (i in 1:length(lat)) {spatial_weights[,i] <- (cos((lat[i]+mean(diff(lat))/2)*pi/180)+cos((lat[i]-mean(diff(lat))/2)*pi/180))/2}
  spatial_weights[which(is.na(clim_field[1:(length(lon)*length(lat))]))] <- NA
  spatial_weights <- spatial_weights/sum(spatial_weights,na.rm=T)
  return(sum(spatial_weights*clim_field,na.rm=T))
}

paleodata_windowing <- function(data,start_date,end_date) {
  if (class(coredata(data)) != "matrix") {
    return(zoo(data[which(index(data) >= start_date & index(data) <= end_date)],order.by=index(data)[which(index(data) >= start_date & index(data) <= end_date)]))
  } else {
    return(zoo(data[which(index(data) >= start_date & index(data) <= end_date),],order.by=index(data)[which(index(data) >= start_date & index(data) <= end_date)]))
  }
}

normalize <- function(x,center=TRUE,scale=TRUE) {
  if (center == TRUE) {
    x <- x - mean(x,na.rm=TRUE)
  }
  if (scale == TRUE) {
    x <- x/std(x)
  }
  return(x)
}

paleodata_interpolation <- function(data,interpolation_type,interpolation_dates) {
  if (interpolation_type=="spectral") {
    if (class(coredata(data)) != "matrix") {
      return(zoo(PaleoSpec::MakeEquidistant(index(data),data,time.target = interpolation_dates),order.by=interpolation_dates))
    } else {
      return(apply_zoo(data,function(x) PaleoSpec::MakeEquidistant(index(data),x,time.target = interpolation_dates),out_dates = interpolation_dates))
    }
  }
  if (interpolation_type=="spline") {
    if (class(coredata(data)) != "matrix") {
      return(zoo(spline(index(data),data,xout=interpolation_dates)$y,order.by=interpolation_dates))
    } else {
      return(apply_zoo(data,function(x) spline(index(data),x,xout=interpolation_dates)$y,out_dates = interpolation_dates))
    }
  }
}

rev_time_axis <- function(data) {
  data <- rev.zoo(data)
  index(data) <- -rev(index(data))
  return(data)
}

#t_test <- function(series, compare_windows, type_data, alternative) {
#  if (length(filter(series, window == compare_windows[1])$site_id) <= 1 | length(filter(series, window == compare_windows[2])$site_id) <= 1) {return(NA)}
#  else {
#    x <- filter(series, window == compare_windows[1])[[GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name]]
#    y <- filter(series, window == compare_windows[2])[[GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name]]
#    
#    t <- t.test(x = x,
#                y = y, 
#                alternative = alternative)$p.value
#    return(t)
#  }
#}


#var_test <- function(series, compare_windows, type_data, alternative) {
#  if (length(filter(series, window == compare_windows[1])$site_id) <= 1 | length(filter(series, window == compare_windows[2])$site_id) <= 1) {return(NA)}
#  else {
#    x <- filter(series, window == compare_windows[1])[[GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name]]
#    y <- filter(series, window == compare_windows[2])[[GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name]]
#    
#    var <- var.test(x = x,
#                    y = y, 
#                    alternative = alternative, 
#                    na.action = na.exclude)$p.value
#    return(var)
#  }
#}


#test_windowed_signal_chars <- function(data, dating, sites_loc, type_data, hres_only = TRUE, site_ids = 'all',
#                                       site_lats = 'all', windows = GLOBAL_SERIES_WINDOWS, labels = GLOBAL_SERIES_WINDOWS_LABELS, 
#                                       compare_windows = list('0-10', '20-30')) {
#  if(!(all(compare_windows %in% labels))) {stop('given compare_windows are not given as labels')}
#  if(length(compare_windows) != 2) {stop('need to provide exactly two compare windows')}
#  sites <- filter_sites(sites = sites_loc, dating = dating, hres_only = hres_only, site_ids = site_ids, site_lats = site_lats)
#  
#  data <- data %>%
#    inner_join(sites, by = 'site_id') %>% 
#    inner_join(dating, by = c('site_id', 'sample_id')) %>% 
#    window_data(data = ., type_data = type_data, windows = windows, labels = labels) %>% 
#    filter(window %in% compare_windows) %>% 
#    group_by(site_id) %>% 
#    filter(!is.na(UQ(sym(GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name)))) %>% 
#    filter(UQ(sym(GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name)) <= 100) %>% 
#    do(t_twosided = t_test(series = ., compare_windows = compare_windows, type_data = type_data, alternative = 'two.sided'), 
#       t_greater = t_test(series = ., compare_windows = compare_windows, type_data = type_data, alternative = 'greater'), 
#       t_less = t_test(series = ., compare_windows = compare_windows, type_data = type_data, alternative = 'less'), 
#       f_twosided = var_test(series = ., compare_windows = compare_windows, type_data = type_data, alternative = 'two.sided'), 
#       f_greater = var_test(series = ., compare_windows = compare_windows, type_data = type_data, alternative = 'greater'), 
#       f_less = var_test(series = ., compare_windows = compare_windows, type_data = type_data, alternative = 'less')) %>% 
#    unnest() %>% 
#    full_join(select(sites, site_id), by = 'site_id')
#  #  inner_join(select(sites_loc, site_id, long, lat), by = 'site_id')
#  return(data)
#}


#compare_keyvals <- function(data, dating, sites_loc, type_data, hres_only = TRUE, site_ids = 'all',
#                            site_lats = 'all', windows = GLOBAL_SERIES_WINDOWS, labels = GLOBAL_SERIES_WINDOWS_LABELS, 
#                            compare_windows = list('0-10', '20-30')) {
#  test_results <- test_windowed_signal_chars(data = data, dating = dating, sites_loc = sites_loc,
#                                             type_data = type_data, hres_only = hres_only, windows = windows, labels = labels,
#                                             compare_windows = compare_windows) %>% categorise_test(., signiveau = 0.05, test = c('t', 'f')) %>% 
#    select(site_id, t_twosided, f_twosided)
#  
#  if(!(all(compare_windows %in% labels))) {stop('given compare_windows are not given as labels')}
#  if(length(compare_windows) != 2) {stop('need to provide exactly two compare windows')}
#  sites <- filter_sites(sites = sites_loc, dating = dating, hres_only = hres_only, site_ids = site_ids, site_lats = site_lats)
#  
#  data <- data %>%
#    inner_join(sites, by = 'site_id') %>% 
#    inner_join(dating, by = c('site_id', 'sample_id')) %>% 
#    window_data(data = ., type_data = type_data, windows = windows, labels = labels) %>% 
#    filter(window %in% compare_windows) %>% 
#    group_by(site_id, window) %>% 
#    filter(!is.na(UQ(sym(GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name)))) %>% 
#    #filter(UQ(sym(GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name)) <= 100) %>% 
#    summarise(mean = mean(!!sym(GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name)), 
#              var = var(!!sym(GLOBAL_SERIES_SIGNALS[[type_data]]$signal_name))) %>% 
#    full_join(select(sites, site_id), by = 'site_id') %>% group_by(site_id)
#  
#  data <- data %>% inner_join(data %>% group_by(site_id) %>% count() %>% filter(n > 1) %>%
#                                select(site_id), by = 'site_id') %>% View()
#  inner_join(., rename(test_results, mean_diff_sig = t_twosided, var_diff_sig = f_twosided), by = 'site_id') %>% 
#    arrange(window, .by_group = T) %>% 
#    spread_n(key = 'window', value = c('mean', 'var')) %>% rowwise() %>% 
#    mutate(mean_diff = abs(!!sym(paste0(compare_windows[2], '_mean'))) - abs(!!sym(paste0(compare_windows[1], '_mean'))),
#           var_diff = abs(!!sym(paste0(compare_windows[2], '_var'))) - abs(!!sym(paste0(compare_windows[1], '_mean')))) %>% 
#    ungroup()
#  data <- bind_cols(filter(data, mean_diff >= 0 & mean_diff_sig) %>% count(name = 'mean_decrease'), 
#                    filter(data, mean_diff < 0 & mean_diff_sig) %>% count(name = 'mean_increase'), 
#                    filter(data, var_diff >= 0 & var_diff_sig) %>% count(name = 'var_decrease'), 
#                    filter(data, var_diff < 0 & var_diff_sig) %>% count(name = 'var_increase'))
#  #  inner_join(select(sites_loc, site_id, long, lat), by = 'site_id')
#  return(data)
#}


#spread_n <- function(df, key, value) {
#  # quote key
#  keyq <- rlang::enquo(key)
#  # break value vector into quotes
#  valueq <- rlang::enquo(value)
#  s <- rlang::quos(!!valueq)
#  df %>% gather(variable, value, !!!s) %>%
#    unite(temp, !!keyq, variable) %>%
#    spread(temp, value)
#}


#categorise_test <- function(data, signiveau, test) {
#  for (t in test) {
#    data <- data %>% 
#      mutate_at(vars(starts_with(as.character(t))), list(~if_else(. > signiveau, FALSE, TRUE)))
#  }
#  return(data)
#}
#
#
#interp_scl_sign <- function(series, interp) {
#  if (length(series$mixed_age) < 3) {warning('series shorter has less then 3 observations returning original series'); return(series)}
#  if (interp$activate_global) {
#    stats <- summarise(series, start = min(mixed_age), stop = max(mixed_age))
#    xout <- seq(stats$start, stats$stop, by = interp$global_scale)
#  }
#  else {
#    stats <- mutate(series, diff = lead(mixed_age) - mixed_age) %>% summarise(res = mean(diff, na.rm = TRUE), start = min(mixed_age), stop = max(mixed_age))
#    xout <- seq(stats$start, stats$stop, by = stats$res)
#  }
#  series <- interp.dataset(y = select(series, -mixed_age, -site_id), x = series$mixed_age, xout = xout, rep.negt = FALSE, method = interp$method) %>% 
#    as_tibble(.)
#  series$mixed_age <- xout
#  return(series)
#}


#interp_scl_all_sign <- function(data, interp) {
#  data <- group_by(data, site_id) %>% 
#    arrange(mixed_age, .by_group = TRUE) %>% 
#    do(data_interp = interp_scl_sign(series = ., interp = interp)) %>% 
#    unnest(data_interp)
#  return(data)
#}
#
#
#calc_acorr_sign <- function(series, signals) {
#  lag_step <- (select(series, mixed_age) %>% mutate(diff = lead(mixed_age) - mixed_age) %>% summarise(lag_step = mean(diff)))$lag_step
#  
#  if (nrow(series) < 3) {warning('series shorter than three observations, returning empty tibble'); return(tibble())}
#  else {
#    series <- lapply(1:length(signals), function(n) {
#      s <- select(series, site_id, !!!sym(signals[n]), mixed_age)
#      return(as_tibble(acf(s[[signals[n]]], type = 'correlation', plot = FALSE, lag.max = length(s$mixed_age) - 1)$acf) %>% rename(UQ(sym(signals[n])) := V1))
#    }) %>% bind_cols() %>% rownames_to_column(var = 'lag') %>% 
#      #series <- bind_cols(as_tibble(acf(series$original, type = 'correlation', plot = FALSE, lag.max = length(series$mixed_age) - 1)$acf) %>% rename(original = V1), 
#      #                    as_tibble(acf(series[['logit & detrend']], type = 'correlation', plot = FALSE, lag.max = length(series$mixed_age) - 1)$acf) %>% rename(UQ(sym('logit & detrend')) := V1), 
#      #                    as_tibble(acf(series[['probit & detrend']], type = 'correlation', plot = FALSE, lag.max = length(series$mixed_age) - 1)$acf) %>% rename(UQ(sym('probit & detrend')) := V1)) %>% 
#      #  rownames_to_column(var = 'lag') %>% 
#      inner_join((select(series, mixed_age) %>% transmute(lag_step = mixed_age - min(mixed_age)) %>% rownames_to_column(var = 'lag')), by = 'lag')
#  }
#  
#  return(series)
#}
#
#
#calc_acorr_all_sign <- function(data, signals) {
#  data <- group_by(data, site_id) %>% 
#    arrange(mixed_age, .by_group = TRUE) %>% 
#    do(data_acorr = calc_acorr_sign(series = ., signals = signals)) %>% 
#    #filter(!is.na(data_acorr)) %>% print() %>% 
#    unnest(data_acorr)
#  return(data)
#}
