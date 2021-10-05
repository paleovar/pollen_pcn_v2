# Script to compute Figure S15 to S19 from the supplement
# - visualisation of paleoclimate networks on a map (individual links & aggregated by continent) of:
#   S17: ACER AP, LOVECLIM TS and PR (6.2-18 ka BP); 
#   S16: HadCM3 TSF, TS, and PR (6-22 ka BP); 
#   S15: TraCE TSF, TS, PR (6-22 ka BP, only individual links); 
#   emulated responses of TraCE TSF (S19) and HadCM3 TSF (S18) to forcing 
#   with TS, PR, CO2, and all three combined (all 6-22 ka BP).
#   Code is similar to figure 3 (see `create_fig3.R`)

#DONE & TESTED

# ACER AP and pseudo proxies from simulations ----
## initialize data (skips automatically if done previously in the session, see `main_pcns.R`) ----
source('main_pcns.R')
source('main_pcn_meas.R')

## Figure S15: TraCE TSF, TS, PR ----
ndg_list <- lapply(names(ACERtrc), function(nm) ACERtrc_ndg_raw %>% rename(ndg = !!sym(paste0('node_degree.',nm))) %>% select(window,site_id,ndg)) %>% 
  setNames(names(ACERtrc))

ndg_lims <- lapply(ndg_list, function(n) n %>% .$ndg) %>% unlist(use.names = FALSE)
ndg_lims <- c(min(ndg_lims),max(ndg_lims))

plts <- lapply(1:length(ACERtrc), function(i) {
  nm <- names(ACERtrc)[i]
  if (nm == 'ACER_ap') dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data' else dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data'
  pltloc <- run_network_plot(ACERtrc[[i]], 
                             orig_data = dat, link_width = 0.3, ndg_lims = ndg_lims,
                             zoom = c(-180, -65, 180, 65), #frac_abs_strongest = 1, 
                             bundled = T, save_plot = list(activate = F), color_strength = T, split_regions_by_sign = T,
                             filter_windows = c('6-22'), node_degree = T, return = T, bg = NULL)
  pltloc$corr_nw
}) %>% 
  setNames(names(ACERtrc))


## combine the 3 networks in the list into one plot
pltst1 <- plot_grid(plts$TRACE_apsb + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank(), plot.margin = margin()),
                    plts$TRACE_pr + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank(), plot.margin = margin()),
                    plts$TRACE_ts + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank(), plot.margin = margin()),
                    align = "hv", axis = "tblr", nrow = 2, 
                    labels = c('A | TraCE TSF', 'B | TraCE PR', 'C | TraCE TS'),
                    label_size = 24, label_fontface = 'bold', label_y = 1.0, hjust = 0, label_x = 0.09)

plt <- pltst1 + 
  draw_grob(get_legend(plts$ACER_ap + theme(legend.justification = 'center',
                                            legend.text = element_text(size = GLOBAL_FONT_SIZE + 4), legend.title = element_text(size = GLOBAL_FONT_SIZE + 8),
                                            legend.background = element_blank(), legend.box.background = element_rect(color = 'black', size = 0.8))),
            x = 0.52, y = 0.12, width = 0.46, height = 0.4)

## save combined plot to disk
## width = 61 for non cairo_pdf
ggsave(plot = plt, width = 62, height = 29.5, units = 'cm', dpi = 200,
       filename = file.path(DIR_FIGURES, 'PCN_TRACEtsfprts_raw_6-22k_updnorm.pdf'), device = cairo_pdf)
ggsave(plot = plt, width = 61, height = 29.5, units = 'cm', dpi = 150,
       filename = file.path(DIR_FIGURES, 'PCN_TRACEtsfprts_raw_6-22k_updnorm.png'))


## Figure S16: HadCM3 TSF, TS, PR ----
## prepare node degree data & limits
ndg_list <- lapply(names(ACERhcm), function(nm) ACERhcm_ndg_raw %>% rename(ndg = !!sym(paste0('node_degree.',nm))) %>% select(window,site_id,ndg)) %>% 
  setNames(names(ACERhcm))

ndg_lims <- lapply(ndg_list, function(n) n %>% .$ndg) %>% unlist(use.names = FALSE)
ndg_lims <- c(min(ndg_lims),max(ndg_lims))


## full network plot for all raw ACER AP links using the generic network plot function of the `PCNdata` S3 class
plts <- lapply(1:length(ACERhcm), function(i) { #1:length(ACERhcm)
  nm <- names(ACERhcm)[i]
  if (nm == 'ACER_ap') dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data' else dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data'
  pltloc <- run_network_plot(ACERhcm[[i]], 
                             orig_data = dat, link_width = 0.3, ndg_lims = ndg_lims, leg_opt = 'cmb',
                             zoom = c(-180, -65, 180, 65), #frac_abs_strongest = 1, 
                             bundled = T, save_plot = list(activate = F), color_strength = T, split_regions_by_sign = T,
                             filter_windows = c('6-22'), node_degree = T, return = T, return_nw = T, bg = NULL)
}) %>% 
  setNames(names(ACERhcm))


## aggregated network plot for ACER AP and HadCM3 networks 
grphs_agg <- lapply(names(ACERhcm), function(nm) {
  if (nm == 'ACER_ap') cnm <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data' else cnm <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data'
  testgraph <- create_bubble_network_from_corr_list(ACERhcm[[nm]][[cnm]])
}) %>% 
  setNames(names(ACERhcm))

width_lims <- lapply(grphs_agg, function(g) g %E>% as_tibble() %>% .$nlinks) %>% unlist(use.names = FALSE)
width_lims <- c(min(width_lims),max(width_lims))

site_data_list <- lapply(names(ACERhcm), function(nm) {
  ndg_list[[nm]] %>% full_join(ACERhcm[[nm]]$sites) %>% select(site_id,lat,long,ndg) %>% mutate(occurs = if_else(is.na(ndg) | ndg == 0,FALSE,TRUE),ndg = if_else(is.na(ndg),0,ndg))
}) %>% 
  setNames(names(ACERhcm))

plts_agg <- lapply(1:length(grphs_agg), function(i) {
  plot_bubble_network_spatial(grphs_agg[[i]],width_lims = width_lims,ndg_lims = ndg_lims,
                              leg_opt = 'cmb',
                              site_data = site_data_list[[i]])
}) %>% 
  setNames(names(ACERhcm))

## combine the 6 networks + 2 legends in one plot
leg1 <- get_legend(plts$BBC_ap$corr_nw + theme(legend.justification = 'center',
                                               legend.text = element_text(size = GLOBAL_FONT_SIZE + 4), legend.title = element_text(size = GLOBAL_FONT_SIZE + 8),
                                               legend.background = element_blank(), legend.spacing = unit(20, 'mm'))) #legend.box.background = element_rect(color = 'black', size = 0.8)))


leg2 <- get_legend(plts_agg$BBC_ap + theme(legend.justification = 'center', 
                                           legend.text = element_text(size = GLOBAL_FONT_SIZE + 4), legend.title = element_text(size = GLOBAL_FONT_SIZE + 8),
                                           legend.spacing = unit(11.1, 'mm')))

plt <- plot_grid(
  plot_grid(plts$BBC_ap$corr_nw + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            NULL,
            plts_agg$BBC_ap + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            NULL,NULL,NULL,
            plts$BBC_pr$corr_nw + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            NULL,
            plts_agg$BBC_pr + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            NULL,NULL,NULL,
            plts$BBC_ts$corr_nw + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            NULL,
            plts_agg$BBC_ts + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            align = "hv", axis = "tblr", ncol = 3, 
            labels = c('A | HadCM3 TSF LINKS', '', 'B | HadCM3 TSF', 
                       '','','',
                       'C | HadCM3 PR LINKS', '', 'D | HadCM3 PR', 
                       '','','',
                       'E | HadCM3 TS LINKS', '', 'F | HadCM3 TS'),
            label_size = 24, label_fontface = 'bold', label_y = 1.0, hjust = 0, label_x = 0.09,
            rel_heights = c(0.32,-0.022,0.32,-0.022,0.32),
            rel_widths = c(0.48,-0.02,0.48)),
  NULL, NULL,
  ncol = 1,
  rel_heights = c(0.77,-0.05,0.12)) + 
  draw_grob(leg1, x = -0.17, y = -0.445) + 
  draw_grob(leg2, x = 0.165, y = -0.442) + 
  draw_grob(grid::rectGrob(gp = grid::gpar(fill = NA, lineheight = 10, fontsize = 20)), 
            x = 0.168, y = 0.008, width = 0.66, height = 0.092)


## save combined plot to disk
## width = 61 for non cairo_pdf
ggsave(plot = plt, width = 61, height = 45.5, units = 'cm', dpi = 200,
       filename = file.path(DIR_FIGURES, 'PCNagg_HadCM3tsfprts_raw_6-22k_updnorm.pdf'), device = cairo_pdf)
ggsave(plot = plt, width = 61, height = 45.5, units = 'cm', dpi = 150,
       filename = file.path(DIR_FIGURES, 'PCNagg_HadCM3tsfprts_raw_6-22k_updnorm.png'))


## Figure S17: ACER AP (6.2-18 ka BP) and LOVECLIM TS, PR ----
## prepare node degree data & limits
ndg_list <- lapply(names(ACERlc), function(nm) ACERlc_ndg_raw %>% rename(ndg = !!sym(paste0('node_degree.',nm))) %>% select(window,site_id,ndg)) %>% 
  setNames(names(ACERlc))

ndg_lims <- lapply(ndg_list, function(n) n %>% .$ndg) %>% unlist(use.names = FALSE)
ndg_lims <- c(min(ndg_lims),max(ndg_lims))


## full network plot for all raw ACER AP links using the generic network plot function of the `PCNdata` S3 class
plts <- lapply(1:length(ACERlc), function(i) { #1:length(ACERlc)
  nm <- names(ACERlc)[i]
  if (nm == 'ACER_ap') dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data' else dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data'
  pltloc <- run_network_plot(ACERlc[[i]], 
                             orig_data = dat, link_width = 0.3, ndg_lims = ndg_lims, leg_opt = 'cmb',
                             zoom = c(-180, -65, 180, 65), #frac_abs_strongest = 1, 
                             bundled = T, save_plot = list(activate = F), color_strength = T, split_regions_by_sign = T,
                             filter_windows = c('6.2-18'), node_degree = T, return = T, return_nw = T, bg = NULL)
}) %>% 
  setNames(names(ACERlc))


## aggregated network plot for ACER AP and LOVECLIM networks 
grphs_agg <- lapply(names(ACERlc), function(nm) {
  if (nm == 'ACER_ap') cnm <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data' else cnm <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data'
  testgraph <- create_bubble_network_from_corr_list(ACERlc[[nm]][[cnm]])
}) %>% 
  setNames(names(ACERlc))

# use lims from above for same scaling
#width_lims <- lapply(grphs_agg, function(g) g %E>% as_tibble() %>% .$nlinks) %>% unlist(use.names = FALSE)
#width_lims <- c(min(width_lims),max(width_lims))

site_data_list <- lapply(names(ACERlc), function(nm) {
  ndg_list[[nm]] %>% full_join(ACERlc[[nm]]$sites) %>% select(site_id,lat,long,ndg) %>% mutate(occurs = if_else(is.na(ndg) | ndg == 0,FALSE,TRUE),ndg = if_else(is.na(ndg),0,ndg))
}) %>% 
  setNames(names(ACERlc))

plts_agg <- lapply(1:length(grphs_agg), function(i) {
  plot_bubble_network_spatial(grphs_agg[[i]],width_lims = width_lims,ndg_lims = ndg_lims,
                              leg_opt = 'cmb',
                              site_data = site_data_list[[i]])
}) %>% 
  setNames(names(ACERlc))

## combine the 6 networks + 2 legends in one plot
leg1 <- get_legend(plts$ACER_ap$corr_nw + theme(legend.justification = 'center',
                                               legend.text = element_text(size = GLOBAL_FONT_SIZE + 4), legend.title = element_text(size = GLOBAL_FONT_SIZE + 8),
                                               legend.background = element_blank(), legend.spacing = unit(20, 'mm'))) #legend.box.background = element_rect(color = 'black', size = 0.8)))


leg2 <- get_legend(plts_agg$ACER_ap + theme(legend.justification = 'center', 
                                           legend.text = element_text(size = GLOBAL_FONT_SIZE + 4), legend.title = element_text(size = GLOBAL_FONT_SIZE + 8),
                                           legend.spacing = unit(11.1, 'mm')))

plt <- plot_grid(
  plot_grid(plts$ACER_ap$corr_nw + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            NULL,
            plts_agg$ACER_ap + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            NULL,NULL,NULL,
            plts$LC_pr$corr_nw + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            NULL,
            plts_agg$LC_pr + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            NULL,NULL,NULL,
            plts$LC_ts$corr_nw + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            NULL,
            plts_agg$LC_ts + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            align = "hv", axis = "tblr", ncol = 3, 
            labels = c('A | ACER AP LINKS', '', 'B | ACER AP', 
                       '','','',
                       'C | LOVECLIM PR LINKS', '', 'D | LOVECLIM PR', 
                       '','','',
                       'E | LOVECLIM TS LINKS', '', 'F | LOVECLIM TS'),
            label_size = 24, label_fontface = 'bold', label_y = 1.0, hjust = 0, label_x = 0.09,
            rel_heights = c(0.32,-0.022,0.32,-0.022,0.32),
            rel_widths = c(0.48,-0.02,0.48)),
  NULL, NULL,
  ncol = 1,
  rel_heights = c(0.77,-0.05,0.12)) + 
  draw_grob(leg1, x = -0.17, y = -0.445) + 
  draw_grob(leg2, x = 0.165, y = -0.442) + 
  draw_grob(grid::rectGrob(gp = grid::gpar(fill = NA, lineheight = 10, fontsize = 20)), 
            x = 0.168, y = 0.008, width = 0.66, height = 0.092)


## save combined plot to disk
## width = 61 for non cairo_pdf
ggsave(plot = plt, width = 61, height = 45.5, units = 'cm', dpi = 200,
       filename = file.path(DIR_FIGURES, 'PCNagg_ACERapLCLIMprts_raw_6_2-18k_updnorm.pdf'), device = cairo_pdf)
ggsave(plot = plt, width = 61, height = 45.5, units = 'cm', dpi = 150,
       filename = file.path(DIR_FIGURES, 'PCNagg_ACERapLCLIMprts_raw_6_2-18k_updnorm.png'))


# Emulated response of TSF of HadCM3 and TraCE ----
## initialize data (skips automatically if done previously in the session, see `main_emulation_pcns.R`) ----
source('main_emulation_pcns.R')
source('main_emulation_pcn_meas.R')

## Figure S18: Emulated response of HadCM3 TSF to TS, PR, CO2, and combined forcing ----
## prepare node degree data & limits
ndg_list <- lapply(names(ACERhcm_em), function(nm) ACERhcm_em_ndg %>% rename(ndg = !!sym(paste0('node_degree.',nm))) %>% select(window,site_id,ndg)) %>% 
  setNames(names(ACERhcm_em))

ndg_lims <- lapply(ndg_list, function(n) n %>% .$ndg) %>% unlist(use.names = FALSE)
ndg_lims <- c(min(ndg_lims),max(ndg_lims))


## full network plot for all raw ACER AP links using the generic network plot function of the `PCNdata` S3 class
plts <- lapply(1:length(ACERhcm_em), function(i) { #1:length(ACERhcm_em)
  nm <- names(ACERhcm_em)[i]
  if (nm == 'ACER_ap') dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data' else dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data'
  pltloc <- run_network_plot(ACERhcm_em[[i]], 
                             orig_data = dat, link_width = 0.3, ndg_lims = ndg_lims, leg_opt = 'cmb',
                             zoom = c(-180, -65, 180, 65), #frac_abs_strongest = 1, 
                             bundled = T, save_plot = list(activate = F), color_strength = T, split_regions_by_sign = T,
                             filter_windows = c('6-22'), node_degree = T, return = T, return_nw = T, bg = NULL)
}) %>% 
  setNames(names(ACERhcm_em))


## aggregated network plot for ACER AP and HadCM3 networks 
grphs_agg <- lapply(names(ACERhcm_em), function(nm) {
  if (nm == 'ACER_ap') cnm <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data' else cnm <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data'
  testgraph <- create_bubble_network_from_corr_list(ACERhcm_em[[nm]][[cnm]])
}) %>% 
  setNames(names(ACERhcm_em))

width_lims <- lapply(grphs_agg, function(g) g %E>% as_tibble() %>% .$nlinks) %>% unlist(use.names = FALSE)
width_lims <- c(min(width_lims),max(width_lims))


site_data_list <- lapply(names(ACERhcm_em), function(nm) {
  ndg_list[[nm]] %>% full_join(ACERhcm_em[[nm]]$sites) %>% select(site_id,lat,long,ndg) %>% mutate(occurs = if_else(is.na(ndg) | ndg == 0,FALSE,TRUE),ndg = if_else(is.na(ndg),0,ndg))
}) %>% 
  setNames(names(ACERhcm_em))

plts_agg <- lapply(1:length(grphs_agg), function(i) {
  plot_bubble_network_spatial(grphs_agg[[i]],width_lims = width_lims,ndg_lims = ndg_lims,
                              leg_opt = 'cmb',
                              site_data = site_data_list[[i]])
}) %>% 
  setNames(names(ACERhcm_em))

## combine the 8 networks + 2 legends in one plot
leg1 <- get_legend(plts$forceall$corr_nw + theme(legend.justification = 'center',
                                               legend.text = element_text(size = GLOBAL_FONT_SIZE + 4), legend.title = element_text(size = GLOBAL_FONT_SIZE + 8),
                                               legend.background = element_blank(), legend.spacing = unit(20, 'mm'))) #legend.box.background = element_rect(color = 'black', size = 0.8)))


leg2 <- get_legend(plts_agg$forceall + theme(legend.justification = 'center', 
                                           legend.text = element_text(size = GLOBAL_FONT_SIZE + 4), legend.title = element_text(size = GLOBAL_FONT_SIZE + 8),
                                           legend.spacing = unit(11.1, 'mm')))

plt <- plot_grid(
  plot_grid(plts$forceall$corr_nw + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            NULL,
            plts_agg$forceall + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            NULL,NULL,NULL,
            plts$forcets$corr_nw + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            NULL,
            plts_agg$forcets + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            NULL,NULL,NULL,
            plts$forcepr$corr_nw + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            NULL,
            plts_agg$forcepr + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            NULL,NULL,NULL,
            plts$forceco2$corr_nw + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            NULL,
            plts_agg$forceco2 + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            align = "hv", axis = "tblr", ncol = 3, 
            labels = c('A | EM TSF ALL LINKS', '', 'B | EM TSF ALL', 
                       '','','',
                       'C | EM TSF TS LINKS', '', 'D | EM TSF TS', 
                       '','','',
                       'E | EM TSF PR LINKS', '', 'F | EM TSF PR',
                       '','','',
                       'G | EM TSF CO2 LINKS', '', 'H | EM TSF CO2'),
            label_size = 24, label_fontface = 'bold', label_y = 1.0, hjust = 0, label_x = 0.09,
            rel_heights = c(0.22,-0.018,0.22,-0.018,0.22,-0.018,0.22),
            rel_widths = c(0.48,-0.02,0.48)),
  NULL, NULL,
  ncol = 1,
  rel_heights = c(0.77,-0.05,0.11)) + 
  draw_grob(leg1, x = -0.17, y = -0.453) + 
  draw_grob(leg2, x = 0.165, y = -0.45) + 
  draw_grob(grid::rectGrob(gp = grid::gpar(fill = NA, lineheight = 10, fontsize = 20)), 
            x = 0.168, y = 0.008, width = 0.66, height = 0.076)


## save combined plot to disk
## width = 61 for non cairo_pdf
ggsave(plot = plt, width = 61, height = 59, units = 'cm', dpi = 200,
       filename = file.path(DIR_FIGURES, 'PCNagg_HadCM3tsf_em_raw_6-22k_updnorm.pdf'), device = cairo_pdf)
ggsave(plot = plt, width = 61, height = 59, units = 'cm', dpi = 150,
       filename = file.path(DIR_FIGURES, 'PCNagg_HadCM3tsf_em_raw_6-22k_updnorm.png'))


## Figure S19: Emulated response of TraCE TSF to TS, PR, CO2, and combined forcing ----
## prepare node degree data & limits
ndg_list <- lapply(names(ACERtrc_em), function(nm) ACERtrc_em_ndg %>% rename(ndg = !!sym(paste0('node_degree.',nm))) %>% select(window,site_id,ndg)) %>% 
  setNames(names(ACERtrc_em))

#ndg_lims <- lapply(ndg_list, function(n) n %>% .$ndg) %>% unlist(use.names = FALSE)
#ndg_lims <- c(min(ndg_lims),max(ndg_lims))


## full network plot for all raw ACER AP links using the generic network plot function of the `PCNdata` S3 class
plts <- lapply(1:length(ACERtrc_em), function(i) { #1:length(ACERtrc_em)
  nm <- names(ACERtrc_em)[i]
  if (nm == 'ACER_ap') dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data' else dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data'
  pltloc <- run_network_plot(ACERtrc_em[[i]], 
                             orig_data = dat, link_width = 0.3, ndg_lims = ndg_lims, leg_opt = 'cmb',
                             zoom = c(-180, -65, 180, 65), #frac_abs_strongest = 1, 
                             bundled = T, save_plot = list(activate = F), color_strength = T, split_regions_by_sign = T,
                             filter_windows = c('6-22'), node_degree = T, return = T, return_nw = T, bg = NULL)
}) %>% 
  setNames(names(ACERtrc_em))


## aggregated network plot for ACER AP and HadCM3 networks 
grphs_agg <- lapply(names(ACERtrc_em), function(nm) {
  if (nm == 'ACER_ap') cnm <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data' else cnm <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data'
  testgraph <- create_bubble_network_from_corr_list(ACERtrc_em[[nm]][[cnm]])
}) %>% 
  setNames(names(ACERtrc_em))

#width_lims <- lapply(grphs_agg, function(g) g %E>% as_tibble() %>% .$nlinks) %>% unlist(use.names = FALSE)
#width_lims <- c(min(width_lims),max(width_lims))

site_data_list <- lapply(names(ACERtrc_em), function(nm) {
  ndg_list[[nm]] %>% full_join(ACERtrc_em[[nm]]$sites) %>% select(site_id,lat,long,ndg) %>% mutate(occurs = if_else(is.na(ndg) | ndg == 0,FALSE,TRUE),ndg = if_else(is.na(ndg),0,ndg))
}) %>% 
  setNames(names(ACERtrc_em))

plts_agg <- lapply(1:length(grphs_agg), function(i) {
  plot_bubble_network_spatial(grphs_agg[[i]],width_lims = width_lims,ndg_lims = ndg_lims,
                              leg_opt = 'cmb',
                              site_data = site_data_list[[i]])
}) %>% 
  setNames(names(ACERtrc_em))

## combine the 8 networks + 2 legends in one plot
leg1 <- get_legend(plts$forceall$corr_nw + theme(legend.justification = 'center',
                                                 legend.text = element_text(size = GLOBAL_FONT_SIZE + 4), legend.title = element_text(size = GLOBAL_FONT_SIZE + 8),
                                                 legend.background = element_blank(), legend.spacing = unit(20, 'mm'))) #legend.box.background = element_rect(color = 'black', size = 0.8)))


leg2 <- get_legend(plts_agg$forceall + theme(legend.justification = 'center', 
                                             legend.text = element_text(size = GLOBAL_FONT_SIZE + 4), legend.title = element_text(size = GLOBAL_FONT_SIZE + 8),
                                             legend.spacing = unit(11.1, 'mm')))

plt <- plot_grid(
  plot_grid(plts$forceall$corr_nw + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            NULL,
            plts_agg$forceall + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            NULL,NULL,NULL,
            plts$forcets$corr_nw + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            NULL,
            plts_agg$forcets + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            NULL,NULL,NULL,
            plts$forcepr$corr_nw + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            NULL,
            plts_agg$forcepr + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            NULL,NULL,NULL,
            plts$forceco2$corr_nw + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            NULL,
            plts_agg$forceco2 + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
            align = "hv", axis = "tblr", ncol = 3, 
            labels = c('A | EM TSF ALL LINKS', '', 'B | EM TSF ALL', 
                       '','','',
                       'C | EM TSF TS LINKS', '', 'D | EM TSF TS', 
                       '','','',
                       'E | EM TSF PR LINKS', '', 'F | EM TSF PR',
                       '','','',
                       'G | EM TSF CO2 LINKS', '', 'H | EM TSF CO2'),
            label_size = 24, label_fontface = 'bold', label_y = 1.0, hjust = 0, label_x = 0.09,
            rel_heights = c(0.22,-0.018,0.22,-0.018,0.22,-0.018,0.22),
            rel_widths = c(0.48,-0.02,0.48)),
  NULL, NULL,
  ncol = 1,
  rel_heights = c(0.77,-0.05,0.11)) + 
  draw_grob(leg1, x = -0.17, y = -0.453) + 
  draw_grob(leg2, x = 0.165, y = -0.45) + 
  draw_grob(grid::rectGrob(gp = grid::gpar(fill = NA, lineheight = 10, fontsize = 20)), 
            x = 0.168, y = 0.008, width = 0.66, height = 0.076)


## save combined plot to disk
## width = 61 for non cairo_pdf
ggsave(plot = plt, width = 61, height = 59, units = 'cm', dpi = 200,
       filename = file.path(DIR_FIGURES, 'PCNagg_TraCEtsf_em_raw_6-22k_updnorm.pdf'), device = cairo_pdf)
ggsave(plot = plt, width = 61, height = 59, units = 'cm', dpi = 150,
       filename = file.path(DIR_FIGURES, 'PCNagg_TraCEtsf_em_raw_6-22k_updnorm.png'))


# remove objects to save some memory
rm(plts_agg,plt,leg1,leg2,site_data_list,grphs_agg,plts,ndg_list)
gc()
