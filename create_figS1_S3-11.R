# Script to compute supplementary Figures S1, and S3 to S11
# - Time series of 
#   S1: ACER AP for all records in the database which cover at least parts of 22-6 ka BP window
#       and of EPICA and NGRIP d18O isotope ratios colored with rolling median inter-sample range
#   S3: TraCE TS and PR time series sampled at the 63 ACER sites used for further analysis
#   S4: As S2 for TraCE TSF
#   S5: As S2 for emulated TraCE TSF response to TS and CO2
#   S6: As S2 for emulated TraCE TSF response to PR and ALL
#   S7: As S2 for HadCM3 pseudo-transient simulation of TS and PR
#   S8: As S2 for HadCM3 pseudo-transient simulation of TSF
#   S9: As S2 for emulated HadCM3 TSF response to TS and CO2
#   S10: As S2 for emulated HadCM3 TSF response to PR and ALL
#   S11: As S2 for LOVECLIM simulation of TS and PR in the period 18-6.2 ka BP

#DONE & TESTED

# load data ----
source('main_simulation_data.R')

# time series plot of ACER AP with coloured median ISR ----
## plot function
ACERp <- PCNdata() %>% 
  mix_sample_dating_and_err() %>% 
  harmonized_pollen_to_db() %>% 
  arboreal_pollen_to_db() %>% 
  extend_dataset_in_db(epica_ngrip_cores %>% mutate(age = age * 1e3), 
                           'd18O',
                           'arb_pollen_data', 
                           'arboreal_pollen', 
                           'AP',
                           'ACERarbpoll_NGECice',
                           'ACERarbpoll_NGECice_data',
                           'proxy', 
                           'proxy_id', 
                           'proxy_type', 
                           tibble(site_name = c('ngrip', 'epica'), 
                                  lat = c(75.1, -75.0), 
                                  long = c(-42.32, 125.0)))

plot_fct <- function(samples, sites, mav_dat = NULL, mav_spec = list(col='status_med_smp_res',title='Median ISR [yrs]'),
                     arrange_by = 'lat', window = NULL, window_tolerance = c(0,0), labs = FALSE, dating_pts = FALSE, rescale = FALSE) {
  if(!('mixed_age' %in% names(samples))) {stop('call mix_sample_and_dating function for ACER_dataset to merge cores')}
  palaeo_spans <- tibble(label = ordered(c('Hol', 'Pleistocene', 'LGM', 'Eemian'), levels = c('Hol', 'Pleistocene', 'LGM', 'Eemian')),
                         start = c(12, 115, 27, 130), end = c(0, 12, 19, 115), anchor = c(6, 51.5, 23, 122.5)) %>% 
    filter(label %in% c('Hol', 'LGM'))
  age_model_labels <- list('original_age' = 'original chronology', 'updated_age_model' = 'harmonised chronology')
  
  samples <- samples %>% 
    {if(!is.null(window)) {filter(., mixed_age < (window[2]+window_tolerance[2]) & mixed_age > (window[1]-window_tolerance[1]))}else {.}} %>% 
    inner_join(sites, by = 'site_id') %>% 
    ungroup() 
  
  if (rescale) {
    samples <- samples %>% 
      group_by(site_id) %>% 
      nest() %>% 
      mutate(data_n = purrr::map(data, function(x) {x$proxy <- scales::rescale(x$proxy, to = c(-1,1)); return(x)})) %>% 
      unnest(data_n)
  }
  
  if (!is.null(mav_dat)) {
    mav_dat <- mav_dat %>% 
      select(sample_id, site_id, mixed_age, !!sym(mav_spec$col))
    samples <- inner_join(samples, mav_dat, by = c('sample_id', 'site_id', 'mixed_age'))
  }
  
  if (dating_pts == TRUE) {
    merged_cores <- list('98' = 90, '97' = 72, '100' = 85)
    tie_pts <- read_csv(paste(DIR_DATASETS, paste('dating_info_ACER', DATAFILES_TYPE, sep = '.'), sep = '/')) %>% select(site_id, age) %>% 
      filter(age != 'HIATUS') %>% rename(tie_pt_age = age) %>% mutate(site_id = purrr::map(site_id, function(x) {if(as.character(x) %in% names(merged_cores)){merged_cores[[as.character(x)]]}else{x}})) %>% 
      unnest() %>% 
      filter(tie_pt_age < (window[2]+window_tolerance[2]) & tie_pt_age > (window[1]-window_tolerance[1]))
    samples <- left_join(samples, tie_pts, by = 'site_id')
  }
  
  samples <- samples %>% 
    mutate(site_id = case_when(site_id == '102' ~ 'NGP', site_id == '101' ~ 'EDC', TRUE ~ as.character(site_id))) %>%  
    mutate(site_id = str_replace_all(site_id, '\\.', '/')) %>% 
    mutate(site_id = if_else(site_id == '85/1', '85/100', site_id)) %>% 
    mutate(label = if_else(hres, str_pad(paste0(site_id, '*'),5,'left'), str_pad(str_pad(as.character(site_id),4,'left'),5,'right'))) 
  if (arrange_by == 'lat') {
    samples$site_id <- factor(samples$site_id, levels = unique(arrange(samples, desc(lat))$site_id), ordered = T)
  }else if(!(arrange_by == 'id')) {
    stop('unknown arranging parameter provided to plot_dating_by_site; allowed: lat, id')
  }
  
  samples <- samples %>% 
    mutate(mixed_age = mixed_age * -1)
  palaeo_spans <- palaeo_spans %>% 
    mutate(end = end * -1, start = start * -1)
  
  plot <- ggplot(data = arrange(samples, desc(lat)),
                 mapping = aes(x = mixed_age / 1000,
                               y = proxy))
  if (rescale) {
    plot <- plot + 
      geom_rect(data = palaeo_spans, inherit.aes = FALSE, mapping = aes(xmin = end, xmax = start, ymin = -1.05, ymax = 1.05), fill = GLOBAL_GREY_MEDIUM, alpha = I(0.5))
  }
  
  if (!is.null(mav_dat)) {
    plot <- plot + 
      geom_line(mapping = aes(color = !!sym(mav_spec$col))) + 
      scale_color_viridis_c(guide = guide_colorbar(title = mav_spec$title, ticks.linewidth = 2, barwidth = 15,
                                                   direction = 'horizontal', barheight = 1, title.hjust = 0.5, label.position = 'top'),
                            end = 0.8,
                            breaks = c(100,300,500), labels = c('100','300',expression(''>=500)))
  } else {
    plot <- plot + 
      geom_line()
  }
  
  if (dating_pts == TRUE) {
    plot <- plot + 
      geom_point(mapping = aes(x = tie_pt_age / 1000, y = 0), shape = I(23), color = 'black', fill = GLOBAL_GREY_MEDIUM, alpha = I(0.75))
  }
  
  if (rescale) {
    plot <- plot + 
      geom_label(mapping = aes(label = label),
                 x = -1*window[2]*1e-3 - 0.9, y = 0, label.size = 0.5, size = GLOBAL_FONT_SIZE - 11, fontface = GLOBAL_FONT_FACE_TITLE) + 
      facet_wrap(~ site_id, scales = 'fixed', ncol = 3, dir = 'v')
  } else {
    plot <- plot + 
      facet_wrap(~ site_id, scales = 'free_y', ncol = 3, dir = 'v')
  }
  plot <- plot + 
    scale_x_continuous(breaks = -1*c(window[1], (window[2]+window[1])/2, window[2]) * 1e-3, labels = c(window[1], (window[2]+window[1])/2, window[2]) * 1e-3) + 
    labs(y = 'proxy [a.u.]', x = 'age [ka BP]') + 
    global_title_and_axis() + 
    guides(fill = FALSE) 
  if (!is.null(window)) {
    plot <- plot + 
      coord_cartesian(xlim = -1*rev(c((window[1] - window_tolerance[1])*1e-3, (window[2] + window_tolerance[2])*1e-3 + 1.5)))
  }
  plot <- plot +
    theme(legend.position = 'bottom', legend.justification = 'bottom', legend.box = 'horizontal', 
          legend.text = element_text(size = GLOBAL_FONT_SIZE + 3), title = element_text(size = GLOBAL_FONT_SIZE + 3), axis.text = element_text(size = GLOBAL_FONT_SIZE + 3), 
          panel.grid.minor = element_blank(), panel.border = element_blank(),# axis.ticks.length = unit(15, 'pt'),
          legend.background = element_blank(),
          strip.text.x = element_blank(), strip.text.y = element_blank(), strip.background = element_blank(), legend.box.spacing = unit(3, 'mm'),
          axis.title.x = element_text(margin = margin(0,0,0,0,'cm')), axis.title.y = element_text(margin = margin(0,0,0,0,'cm'))) #+ 
  if (rescale) {
    plot <- plot + 
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), panel.grid.major = element_blank(), 
            panel.spacing.y = unit(0, units = 'lines'), axis.line.x = element_line(colour = 'black', size = 0.5))
  } else {
    plot <- plot + 
      theme(axis.line = element_line(colour = 'black', size = 0.5))
  }
  
  plot <- ggdraw(plot) + 
    annotate('text', x = c(0.275,0.6,0.925), y = c(0.0992, 0.0992, 0.158), label = 'HOL', size = GLOBAL_FONT_SIZE - 10, fontface = GLOBAL_FONT_FACE_TITLE, color = GLOBAL_GREY_DARK) + 
    annotate('text', x = c(0.087,0.413,0.74), y = c(0.0992, 0.0992, 0.158), label = 'LGM', size = GLOBAL_FONT_SIZE - 10, fontface = GLOBAL_FONT_FACE_TITLE, color = GLOBAL_GREY_DARK)
  
  return(plot)
}

acer_mav <- compute_moving_average(ACERp$sample_dating, window = 28, limit = 500) # moving average ISR
sites_fltr <- filter_sites(ACERp$sites, ACERp$sample_dating, hres_only = T, 'all', 'all') %>% 
  mutate(hres = TRUE)
plt <- plot_fct(inner_join(ACERp$sample_dating, ACERp$ACERarbpoll_NGECice_data) %>% select(site_id, sample_id, proxy_id, mixed_age, proxy) %>% 
                  full_join(sites_fltr, by = 'site_id') %>% mutate(hres = if_else(is.na(hres), FALSE, hres)),
                ACERp$sites, window = c(6000,22000), dating_pts = F, rescale = T, mav_dat = acer_mav)
global_save_plot_cairo(plt, list(filename = file.path(DIR_FIGURES,'all_ts_medISR'), width = 48, height = 35, png_and_pdf = F))

# pseudo proxies from simulation data ----
## plot function
plot_fct2 <- function(samples, sites, arrange_by = 'lat', window = NULL, window_tolerance = c(0,0), labs = FALSE, rescale = FALSE,
                      override.paleoy = NULL) {
  if(!('mixed_age' %in% names(samples))) {stop('call mix_sample_and_dating function for ACER_dataset to merge cores')}
  palaeo_spans <- tibble(label = ordered(c('Hol', 'Pleistocene', 'LGM', 'Eemian'), levels = c('Hol', 'Pleistocene', 'LGM', 'Eemian')),
                         start = c(12, 115, 27, 130), end = c(0, 12, 19, 115), anchor = c(6, 51.5, 23, 122.5)) %>% 
    filter(label %in% c('Hol', 'LGM'))
  age_model_labels <- list('original_age' = 'original chronology', 'updated_age_model' = 'harmonised chronology')
  
  samples <- samples %>% 
    {if(!is.null(window)) {filter(., mixed_age < (window[2]+window_tolerance[2]) & mixed_age > (window[1]-window_tolerance[1]))}else {.}} %>% 
    inner_join(sites, by = 'site_id') %>% 
    ungroup() #%>% 
  #mutate(site_id = floor(site_id) %>% as.integer(.))
  
  if (rescale) {
    samples <- samples %>% 
      group_by(site_id,model) %>% 
      nest() %>% 
      mutate(data_n = purrr::map(data, function(x) {x$proxy <- scales::rescale(x$proxy, to = c(-1,1)); return(x)})) %>% 
      unnest(data_n)
  }
  
  samples <- samples %>% 
    mutate(site_id = case_when(site_id == '102' ~ 'NGP', site_id == '101' ~ 'EDC', TRUE ~ as.character(site_id))) %>%  
    mutate(site_id = str_replace_all(site_id, '\\.', '/')) %>% 
    mutate(site_id = if_else(site_id == '85/1', '85/100', site_id)) %>% 
    mutate(label = if_else(hres, str_pad(paste0(site_id, '*'),5,'left'), str_pad(str_pad(as.character(site_id),4,'left'),5,'right'))) 
  
  if (arrange_by == 'lat') {
    samples$site_id <- factor(samples$site_id, levels = unique(arrange(samples, desc(lat))$site_id), ordered = T)
  }else if(!(arrange_by == 'id')) {
    stop('unknown arranging parameter provided to plot_dating_by_site; allowed: lat, id')
  }
  
  samples <- samples %>% 
    mutate(mixed_age = mixed_age * -1)
  palaeo_spans <- palaeo_spans %>% 
    mutate(end = end * -1, start = start * -1)
  
  plot <- ggplot(data = arrange(samples, desc(lat)),
                 mapping = aes(x = mixed_age / 1000,
                               y = proxy))
  if (rescale) {
    plot <- plot + 
      geom_rect(data = palaeo_spans, inherit.aes = FALSE, mapping = aes(xmin = end, xmax = start, ymin = -1.05, ymax = 1.05), fill = GLOBAL_GREY_MEDIUM, alpha = I(0.5)) 
  }
  
  plot <- plot + 
    geom_line(mapping = aes(color = model)) + 
    scale_color_viridis_d(end = 0.75, begin = 0.3,
                          guide = guide_legend('Model Data', override.aes = list(size = 5),
                                               direction = 'horizontal'))
  
  if (rescale) {
    plot <- plot + 
      geom_label(mapping = aes(label = label),
                 x = -1*window[2]*1e-3 - 0.9, y = 0, label.size = 0.5, size = GLOBAL_FONT_SIZE - 11, fontface = GLOBAL_FONT_FACE_TITLE) + 
      facet_wrap(~ site_id, scales = 'fixed', ncol = 3, dir = 'v')
    #facet_grid(scales = 'fixed', rows = vars(site_id))
  } else {
    plot <- plot + 
      facet_wrap(~ site_id, scales = 'free_y', ncol = 3, dir = 'v')
    #facet_grid(scales = 'free_y', rows = vars(site_id))
  }
  plot <- plot + 
    scale_x_continuous(breaks = -1*c(window[1], (window[2]+window[1])/2, window[2]) * 1e-3, labels = c(window[1], (window[2]+window[1])/2, window[2]) * 1e-3) + 
    labs(y = 'proxy [a.u.]', x = 'age [ka BP]') + 
    global_title_and_axis() + 
    guides(fill = FALSE) 
  if (!is.null(window)) {
    plot <- plot + 
      coord_cartesian(xlim = -1*rev(c((window[1] - window_tolerance[1])*1e-3 - 0.2, (window[2] + window_tolerance[2])*1e-3 + 1.5)))
  }
  plot <- plot +
    theme(legend.position = 'bottom', legend.justification = 'bottom', legend.box = 'horizontal', legend.direction = 'horizontal',
          legend.text = element_text(size = GLOBAL_FONT_SIZE + 3), title = element_text(size = GLOBAL_FONT_SIZE + 3), axis.text = element_text(size = GLOBAL_FONT_SIZE + 3), 
          panel.grid.minor = element_blank(), panel.border = element_blank(),# axis.ticks.length = unit(15, 'pt'),
          legend.background = element_blank(),
          strip.text.x = element_blank(), strip.text.y = element_blank(), strip.background = element_blank(), legend.box.spacing = unit(3, 'mm'),
          axis.title.x = element_text(margin = margin(0,0,0,0,'cm')), axis.title.y = element_text(margin = margin(0,0,0,0,'cm'))) #+ 
  if (rescale) {
    plot <- plot + 
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), panel.grid.major = element_blank(), 
            panel.spacing.y = unit(0, units = 'lines'), axis.line.x = element_line(colour = 'black', size = 0.5))
  } else {
    plot <- plot + 
      theme(axis.line = element_line(colour = 'black', size = 0.5))
  }
  
  plot <- ggdraw(plot) + 
    annotate('text', x = c(0.275,0.6,0.925), y = {if (!is.null(override.paleoy)) {override.paleoy} else {c(0.0972, 0.0972, 0.0972)}}, label = 'HOL', size = GLOBAL_FONT_SIZE - 10, fontface = GLOBAL_FONT_FACE_TITLE, color = GLOBAL_GREY_DARK) + 
    annotate('text', x = c(0.087,0.413,0.74), y = {if (!is.null(override.paleoy)) {override.paleoy} else {c(0.0972, 0.0972, 0.0972)}}, label = 'LGM', size = GLOBAL_FONT_SIZE - 10, fontface = GLOBAL_FONT_FACE_TITLE, color = GLOBAL_GREY_DARK)
  
  return(plot)
}

## TRACE TSF, PR, TS time series for used sites ----
plt_dat <- bind_rows(inner_join(ACERtrc$TRACE_apsb$sample_dating, ACERtrc$TRACE_apsb$arb_pollen_data) %>% mutate(model = 'TraCE TSF')) %>% #,
  rename(proxy_id = arb_pollen_data_id, proxy = pcnt_arb_pollen) %>% 
  select(site_id, sample_id, model, proxy_id, mixed_age, proxy) %>% 
  inner_join(sites_fltr, by = 'site_id')

plt <- plot_fct2(plt_dat,
                 ACERp$sites, window = c(6000,22000), rescale = T, override.paleoy = c(0.1063, 0.1063, 0.1063))
global_save_plot_cairo(plt, list(filename = file.path(DIR_FIGURES,'filtered_ts_TRACEtsf_updnorm'), width = 48, height = 24.5, png_and_pdf = F))

plt_dat <- bind_rows(inner_join(ACERtrc$TRACE_pr$sample_dating, ACERtrc$TRACE_pr$arb_pollen_data) %>% mutate(model = 'TraCE PR'),
                     inner_join(ACERtrc$TRACE_ts$sample_dating, ACERtrc$TRACE_ts$arb_pollen_data) %>% mutate(model = 'TraCE TS')) %>% 
  rename(proxy_id = arb_pollen_data_id, proxy = pcnt_arb_pollen) %>% 
  select(site_id, sample_id, model, proxy_id, mixed_age, proxy) %>% 
  inner_join(sites_fltr, by = 'site_id')

plt <- plot_fct2(plt_dat,
                 ACERp$sites, window = c(6000,22000), rescale = T)
global_save_plot_cairo(plt, list(filename = file.path(DIR_FIGURES,'filtered_ts_TRACEprts'), width = 48, height = 27, png_and_pdf = F))

## LCLIM TS, PR time series for used sites ----
plt_dat <- bind_rows(inner_join(ACERlc$LC_pr$sample_dating, ACERlc$LC_pr$arb_pollen_data) %>% mutate(model = 'LOVECLIM PR'),
                     inner_join(ACERlc$LC_ts$sample_dating, ACERlc$LC_ts$arb_pollen_data) %>% mutate(model = 'LOVECLIM TS')) %>% 
  rename(proxy_id = arb_pollen_data_id, proxy = pcnt_arb_pollen) %>% 
  select(site_id, sample_id, model, proxy_id, mixed_age, proxy) %>% 
  inner_join(sites_fltr, by = 'site_id')

plt <- plot_fct2(plt_dat,
                 ACERp$sites, window = c(6000,22000), rescale = T, override.paleoy = c(0.1063, 0.1063, 0.1063))
global_save_plot_cairo(plt, list(filename = file.path(DIR_FIGURES,'filtered_ts_LCLIMprts'), width = 48, height = 24.5, png_and_pdf = F))

## HadCM3 TSF, TS, PR time series for used sites ----
plt_dat <- bind_rows(inner_join(ACERhcm$BBC_pr$sample_dating, ACERhcm$BBC_pr$arb_pollen_data) %>% mutate(model = 'HadCM3 PR'),
                     inner_join(ACERhcm$BBC_ts$sample_dating, ACERhcm$BBC_ts$arb_pollen_data) %>% mutate(model = 'HadCM3 TS')) %>% 
  rename(proxy_id = arb_pollen_data_id, proxy = pcnt_arb_pollen) %>% 
  select(site_id, sample_id, model, proxy_id, mixed_age, proxy) %>% 
  inner_join(sites_fltr, by = 'site_id')

plt <- plot_fct2(plt_dat,
                 ACERp$sites, window = c(6000,22000), rescale = T, override.paleoy = c(0.1005, 0.1005, 0.142))
global_save_plot_cairo(plt, list(filename = file.path(DIR_FIGURES,'filtered_ts_BRIDGEprts'), width = 48, height = 26, png_and_pdf = F))

plt_dat <- bind_rows(inner_join(ACERhcm$BBC_ap$sample_dating, ACERhcm$BBC_ap$arb_pollen_data) %>% mutate(model = 'HadCM3 TSF')) %>% #,
  rename(proxy_id = arb_pollen_data_id, proxy = pcnt_arb_pollen) %>% 
  select(site_id, sample_id, model, proxy_id, mixed_age, proxy) %>% 
  inner_join(sites_fltr, by = 'site_id')

plt <- plot_fct2(plt_dat,
                 ACERp$sites, window = c(6000,22000), rescale = T, override.paleoy = c(0.1005, 0.1005, 0.142))
global_save_plot_cairo(plt, list(filename = file.path(DIR_FIGURES,'filtered_ts_BRIDGEtsf_updnorm'), width = 48, height = 26, png_and_pdf = F))

## emulated TRACE TSF response to TS, PR, CO2, and combined forcing for used sites ----
plt_dat <- bind_rows(inner_join(ACERtrc_em$forcets$sample_dating, ACERtrc_em$forcets$arb_pollen_data) %>% mutate(model = 'EM TSF TS'),
                     inner_join(ACERtrc_em$forceco2$sample_dating, ACERtrc_em$forceco2$arb_pollen_data) %>% mutate(model = 'EM TSF CO2')) %>% 
  rename(proxy_id = arb_pollen_data_id, proxy = pcnt_arb_pollen) %>% 
  select(site_id, sample_id, model, proxy_id, mixed_age, proxy) %>% 
  inner_join(sites_fltr, by = 'site_id')

plt <- plot_fct2(plt_dat,
                 ACERp$sites, window = c(6000,22000), rescale = T, override.paleoy = c(0.1063, 0.1063, 0.1063))
global_save_plot_cairo(plt, list(filename = file.path(DIR_FIGURES,'filtered_ts_TRACE_em_tsco2_updnorm'), width = 48, height = 24.5, png_and_pdf = F))

plt_dat <- bind_rows(inner_join(ACERtrc_em$forcepr$sample_dating, ACERtrc_em$forcepr$arb_pollen_data) %>% mutate(model = 'EM TSF PR'),
                     inner_join(ACERtrc_em$forceall$sample_dating, ACERtrc_em$forceall$arb_pollen_data) %>% mutate(model = 'EM TSF ALL')) %>% 
  rename(proxy_id = arb_pollen_data_id, proxy = pcnt_arb_pollen) %>% 
  select(site_id, sample_id, model, proxy_id, mixed_age, proxy) %>% 
  inner_join(sites_fltr, by = 'site_id')

plt <- plot_fct2(plt_dat,
                 ACERp$sites, window = c(6000,22000), rescale = T)
global_save_plot_cairo(plt, list(filename = file.path(DIR_FIGURES,'filtered_ts_TRACE_em_prall_prts'), width = 48, height = 27, png_and_pdf = F))

## emulated HadCM3 TSF response to TS, PR, CO2, and combined forcing for used sites ----
plt_dat <- bind_rows(inner_join(ACERhcm_em$forcets$sample_dating, ACERhcm_em$forcets$arb_pollen_data) %>% mutate(model = 'EM TSF TS'),
                     inner_join(ACERhcm_em$forceco2$sample_dating, ACERhcm_em$forceco2$arb_pollen_data) %>% mutate(model = 'EM TSF CO2')) %>% 
  rename(proxy_id = arb_pollen_data_id, proxy = pcnt_arb_pollen) %>% 
  select(site_id, sample_id, model, proxy_id, mixed_age, proxy) %>% 
  inner_join(sites_fltr, by = 'site_id')

plt <- plot_fct2(plt_dat,
                 ACERp$sites, window = c(6000,22000), rescale = T, override.paleoy = c(0.1063, 0.1063, 0.1063))
global_save_plot_cairo(plt, list(filename = file.path(DIR_FIGURES,'filtered_ts_HadCM3_em_tsco2_updnorm'), width = 48, height = 24.5, png_and_pdf = F))

plt_dat <- bind_rows(inner_join(ACERhcm_em$forcepr$sample_dating, ACERhcm_em$forcepr$arb_pollen_data) %>% mutate(model = 'EM TSF PR'),
                     inner_join(ACERhcm_em$forceall$sample_dating, ACERhcm_em$forceall$arb_pollen_data) %>% mutate(model = 'EM TSF ALL')) %>% 
  rename(proxy_id = arb_pollen_data_id, proxy = pcnt_arb_pollen) %>% 
  select(site_id, sample_id, model, proxy_id, mixed_age, proxy) %>% 
  inner_join(sites_fltr, by = 'site_id')

plt <- plot_fct2(plt_dat,
                 ACERp$sites, window = c(6000,22000), rescale = T)
global_save_plot_cairo(plt, list(filename = file.path(DIR_FIGURES,'filtered_ts_HadCM3_em_prall_prts'), width = 48, height = 27, png_and_pdf = F))

# garbage collect
gc()
