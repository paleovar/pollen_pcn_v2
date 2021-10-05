# Script to compute Figure 6 from the manuscript
# - Visualisation of millennial-scale paleoclimate network links between ACER AP sites and 
#   EPICA & NGRIP d18O isotopic ratio and PR, TS at ACER and icecore sites on a map, 
#   together with a panel of exemplary time series from model and data which are filtered for
#   millennial time scales of variability.
#   Code is similar to Figures 3, S15 to S19 (see `create_fig3.R` and `create_figS15-19.R`).

#DONE & TESTED

# load/compute data if not done previously
source('main_pcns.R')

# prepare network panels for links from/to NGRIP and EPICA DC only (site ids 101, 102) ----
## use ndg_lims from `create_fig3` to maintain same scaling throughout the manuscript
## to use individual scaling, uncomment these lines
#ndg_list <- lapply(names(ACERtrcmil), function(nm) ACERtrc_ndg_mil %>% rename(ndg = !!sym(paste0('node_degree.',nm))) %>% select(window,site_id,ndg)) %>% 
#  setNames(names(ACERtrcmil))
#
#ndg_lims <- lapply(ndg_list, function(n) n %>% .$ndg) %>% unlist(use.names = FALSE)
#ndg_lims <- c(min(ndg_lims),max(ndg_lims))

plts <- lapply(c(2,3,4), function(i) { # 1 would be TraCE TSF
  nm <- names(ACERtrcmil)[i]
  print(nm)
  if (nm == 'ACER_ap') dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_detrend_gaussbp2000-8000_probit_arb_pollen_data' else dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_detrend_gaussbp2000-8000_identity_arb_pollen_data'
  pltloc <- ACERtrcmil[[i]] %>% 
    run_network_plot(orig_data = dat, link_width = 0.85, ndg_lims = ndg_lims,
                     zoom = c(-180, -90, 180, 90), filter_sites = list('6-22' = c(101,102)),
                     bundled = T, save_plot = list(activate = F), color_strength = F, split_regions_by_sign = F,
                     filter_windows = c('6-22'), node_degree = T, return = T, bg = NULL)
  pltloc$corr_nw
}) %>% 
  setNames(names(ACERtrcmil)[c(2,3,4)])

# for saving the network panels separately
#anot <- annotate('label', label = c('NGP', 'EDC'), x = c(0.4425, 0.695), y = c(0.88, 0.12), hjust = c(0,1),
#                 fill = 'white', fontface = GLOBAL_FONT_FACE_TITLE, size = GLOBAL_FONT_SIZE - 11.5, label.size = 0.5)
#pltst1 <- plot_grid(ggdraw(plts$ACER_ap + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank())) + anot,
#                    ggdraw(plts$TRACE_pr + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank())) + anot,
#                    ggdraw(plts$TRACE_ts + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank())) + anot,
#                    align = "hv", axis = "tblr", nrow = 2, 
#                    labels = c('A | ACER AP', 'B | TRACE PR', 'C | TRACE TS'),
#                    label_size = 24, label_fontface = 'bold', label_y = 1.0)
#
#plt <- pltst1 + 
#  draw_grob(get_legend(plts$ACER_ap + theme(legend.justification = 'center', legend.box = 'vertical', legend.box.just = 'left', 
#                                            legend.text = element_text(size = GLOBAL_FONT_SIZE + 4), legend.title = element_text(size = GLOBAL_FONT_SIZE + 8),
#                                            legend.background = element_blank(), legend.box.background = element_rect(color = 'black', size = 0.8))),
#            x = 0.125, y = -0.25)
#
#ggsave(plot = plt, width = 61, height = 33, units = 'cm', dpi = 200,
#       filename = file.path(DIR_FIGURES, 'paperMIS1-2', 'PCN_ACERapTRACEtsfprts_mil2-8ka_6-22k_updnorm_EDC_NGP_theta1000.pdf'))
#ggsave(plot = plt, width = 61, height = 33, units = 'cm', dpi = 100,
#       filename = file.path(DIR_FIGURES, 'paperMIS1-2', 'PCN_ACERapTRACEtsfprts_mil2-8ka_6-22k_updnorm_EDC_NGP_theta1000.png'))


# prepare filtered time series panel ----
## data for the time series panel
ACERtrcdatmil <- lapply(1:length(ACERtrcmil), function(i){
  nm <- names(ACERtrcmil)[i]
  dat <- ACERtrcmil[[nm]]$arb_pollen_data %>% 
    rename(!!sym(paste0(nm, '_pcnt_arb_pollen')) := pcnt_arb_pollen) %>% 
    select(-arb_pollen_data_id)
  dat
})%>% 
  plyr::join_all(., type = 'inner', by = c('site_id', 'sample_id')) %>% 
  as_tibble() %>% 
  gather(key = 'data', value = 'pcnt_arb_pollen', paste0(names(ACERtrcmil), '_pcnt_arb_pollen')) %>% 
  mutate(data = str_replace(data,'_pcnt_arb_pollen','')) %>% 
  group_by(site_id, data) %>% 
  #mutate(pcnt_arb_pollen = (pcnt_arb_pollen-mean(pcnt_arb_pollen))/sd(pcnt_arb_pollen)) %>% 
  inner_join(ACERtrcmil$ACER_ap$sample_dating %>% select(mixed_age, site_id, sample_id), by = c('site_id', 'sample_id')) %>% 
  ungroup() %>% 
  mutate(data = str_replace(data, '_', ' ')) %>% 
  filter(mixed_age >= 6000 & mixed_age <= 22000) %>% 
  filter(data != 'TRACE vc') %>% 
  rowwise() %>% 
  mutate(data = if_else(data == 'TRACE apsb', 'TRACE tsf', data)) %>% 
  mutate(data = toupper(data)) %>% 
  ungroup() %>% 
  mutate(data = str_replace(data, 'TRACE', 'TraCE')) %>% 
  filter(site_id > 100) %>% 
  inner_join(ACERtrcmil$ACER_ap$sites %>% select(site_id,site_name), by = 'site_id') %>% 
  filter(data != "TraCE TSF") # filters TraCE PR, TS and ACER AP (contains d18O for NGP & EDC)

## filter the data to 2-8 ka window
ACERtrcdatmilflt <- ACERtrcdatmil %>% 
  group_by(site_id, site_name, data) %>% 
  nest(.key = 'tseries') %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(tseries_zoo = list(zoo::zoo(x = tseries$pcnt_arb_pollen, order.by = tseries$mixed_age))) %>% 
  mutate(tseries_flt = list(gaussbandpass(tseries_zoo, per1 = 2000, per2 = 8000)$filt)) %>% 
  select(site_name, data, tseries_flt) %>% 
  mutate(tseries_flt = list(tibble(time = index(tseries_flt), 
                                   signal = case_when(data == 'TraCE PR' ~ (coredata(tseries_flt))*60*60*24*30*12*1000,
                                                      TRUE ~ coredata(tseries_flt)))))

## labels and styling for the time series plot
all_data_tseriesmil_labs <- c("bold(atop(NA,atop('\u0394PR','[mm/a]')))",
                              "bold(atop(NA,atop('\u0394PR','[mm/a]')))", 
                              "bold(atop(NA,atop('\u0394TS','[K]')))", 
                              "bold(atop(NA,atop('\u0394TS','[K]')))",
                              "bold(atop(NA,atop(\u0394\u03B4^{18}*'O','[\u2030]')))", 
                              "bold(atop(NA,atop(\u0394\u03B4^{18}*'O','[\u2030]')))")

all_data_tseriesmil_grps <- c('TraCE @ NGP', 'TraCE @ EDC', 'TraCE @ NGP', 'TraCE @ EDC', 'NGP', 'EDC')

all_data_tseriesmil_u <- ACERtrcdatmilflt %>% 
  add_column(signal_name = all_data_tseriesmil_labs, 
             group = all_data_tseriesmil_grps) %>% 
  unnest(cols = c(tseries_flt))

all_data_tseriesmil_u$signal_name <- factor(all_data_tseriesmil_u$signal_name, 
                                            levels = c(all_data_tseriesmil_labs[5], all_data_tseriesmil_labs[3], all_data_tseriesmil_labs[1]))

all_data_tseriesmil_u$group <- factor(all_data_tseriesmil_u$group,
                                      levels = c(all_data_tseriesmil_grps[5], all_data_tseriesmil_grps[6], 
                                                 all_data_tseriesmil_grps[1], all_data_tseriesmil_grps[2]))

clrs <- c('#d7191c', '#e69f00', '#009e73', '#0072b3') # '#6e77a7','#cc79a7', 
alphas <- c(1,1,1,1,1,1)

# plot filtered time series and tweak the grobs
tseriesmil_plt <- ggplot(all_data_tseriesmil_u,
                         aes(color = group, x = 1*time/1000, y = signal, alpha = group)) + 
  geom_line() + 
  facet_wrap(. ~signal_name, scales = 'free_y', strip.position = 'left', ncol = 1, labeller = label_parsed, drop = F) + 
  scale_y_continuous(position = 'right', breaks = scales::pretty_breaks(n = 3)) +#sec.axis = sec_axis(~ .)) + 
  scale_x_continuous(labels = c(6,14,22), breaks = c(6,14,22), expand = c(0.01,0.01)) + 
  scale_alpha_manual(values = alphas, guide = FALSE) + 
  scale_color_manual(values = clrs, guide = guide_legend('Proxy/Simulation', override.aes = list(size = 4), title.position = 'top')) + 
  global_title_and_axis() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.y = element_blank(), axis.text = element_text(size = GLOBAL_FONT_SIZE + 3), axis.title.x = element_text(size = GLOBAL_FONT_SIZE + 3),
        legend.text = element_text(size = GLOBAL_FONT_SIZE + 3), legend.title = element_text(size = GLOBAL_FONT_SIZE + 3),
        strip.placement = 'outside', strip.background.y = element_blank(), strip.text = element_text(size = GLOBAL_FONT_SIZE + 11),
        panel.spacing.y = unit(0.02, 'npc'), legend.position = 'bottom', legend.direction = 'horizontal',
        panel.border = element_rect(size = 0.7), legend.box.background = element_rect(color = 'black', size = 0.7)) + #, panel.background = element_rect(color = 'black', inherit.blank = F, size = 1)) + 
  labs(x = 'time [ka BP]')

## for saving model series separately
#tseriesmil_grb <- ggplotGrob(tseriesmil_plt)
#tseriesmil_grb$grobs[[2]]$children[[5]] <- segmentsGrob(c(1,0),c(0,0),c(1,0),c(1,1))
#tseriesmil_grb$grobs[[3]]$children[[5]] <- segmentsGrob(c(1,0),c(0,0),c(1,0),c(1,1))
#tseriesmil_grb$grobs[[4]]$children[[5]] <- segmentsGrob(c(1,0,0),c(0,0,0),c(1,0,1),c(1,1,0))
#cmbplt_milseries <- gridExtra::arrangeGrob(tseriesmil_grb)
#ggsave(plot = cmbplt_milseries, file.path(DIR_FIGURES, 'PCN_ACERapTRACEtsfprts_mil2-8ka_6-22k_updnorm_EDC_NGP_onlytseriesmil.pdf'), 
#       width = 38, height = 22, units = 'cm', dpi = 'print', limitsize = FALSE, device = cairo_pdf)


# combine model series with networks and add more labels ----
tseriesmil_leg <- get_legend(tseriesmil_plt)
tseriesmil_grb <- ggplotGrob(tseriesmil_plt + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()))
tseriesmil_grb$grobs[[2]]$children[[5]] <- segmentsGrob(c(1,0),c(0,0),c(1,0),c(1,1))
tseriesmil_grb$grobs[[3]]$children[[5]] <- segmentsGrob(c(1,0),c(0,0),c(1,0),c(1,1))
tseriesmil_grb$grobs[[4]]$children[[5]] <- segmentsGrob(c(1,0,0),c(0,0,0),c(1,0,1),c(1,1,0))

anot <- annotate('label', label = c('NGP', 'EDC'), x = c(0.445, 0.69), y = c(0.88, 0.1), hjust = c(0,1),
                 fill = 'white', fontface = GLOBAL_FONT_FACE_TITLE, size = GLOBAL_FONT_SIZE - 10.5, label.size = 0.5)
pltst1 <- plot_grid(ggdraw(plts$ACER_ap + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank())) + anot + 
                      annotate('text', parse = T, label = "bold(atop(NA,atop('A | ACER AP', '   ice '*\u03B4^{18}*'O')))", 
                               x = 0.01, y = 0.89, vjust = 0, hjust = 0, size = GLOBAL_FONT_SIZE-5),
                    ggdraw(plts$TRACE_pr + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank())) + anot + 
                      annotate('text', parse = T, label = "bold(atop(NA,atop('B | TraCE PR', NA)))", 
                               x = 0.01, y = 0.89, vjust = 0, hjust = 0, size = GLOBAL_FONT_SIZE-5),
                    ggdraw(plts$TRACE_ts + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank())) + anot + 
                      annotate('text', parse = T, label = "bold(atop(NA,atop('C | TraCE TS', NA)))", 
                               x = 0.01, y = 0.89, vjust = 0, hjust = 0, size = GLOBAL_FONT_SIZE-5),
                    align = "hv", axis = "tblr", nrow = 2) #, 

plt <- ggdraw() + 
  draw_plot(pltst1, y = 0.1, height = 0.9) +
  draw_plot(plot = tseriesmil_grb,
            x = 0.538, y = 0.08, width = 0.4, height = 0.45) + 
  annotate('text', parse = T, label = "bold(atop(NA,atop('D', NA)))", 
           x = 0.51, y = 0.49, vjust = 0, hjust = 0, size = GLOBAL_FONT_SIZE-5) + 
  draw_grob(get_legend(plts$ACER_ap + theme(legend.justification = 'center',# legend.box = 'vertical', legend.box.just = 'left', legend.text = element_text(size = GLOBAL_FONT_SIZE + 4), legend.title = element_text(size = GLOBAL_FONT_SIZE + 8),
                                            legend.background = element_blank(), legend.box.background = element_rect(color = 'black', size = 0.7),
                                            legend.text = element_text(size = GLOBAL_FONT_SIZE + 3), legend.title = element_text(size = GLOBAL_FONT_SIZE + 3))),
            x = 0.0, y = -0.18, width = 0.5, height = 0.45) + 
  draw_grob(tseriesmil_leg,
            x = 0.5, y = -0.18, width = 0.5, height = 0.45)

# save the combined plot ----
ggsave(plot = plt, width = 63, height = 36, units = 'cm', dpi = 200,
       filename = file.path(DIR_FIGURES, 'PCN_ACERapTRACEtsfprts_mil2-8ka_6-22k_updnorm_EDC_NGP_wtseriesmil2-8ka_theta190.pdf'), limitsize = FALSE, device = cairo_pdf)
ggsave(plot = plt, width = 63, height = 35, units = 'cm', dpi = 150,
       filename = file.path(DIR_FIGURES, 'PCN_ACERapTRACEtsfprts_mil2-8ka_6-22k_updnorm_EDC_NGP_wtseriesmil2-8ka_theta190.png'), limitsize = FALSE)

