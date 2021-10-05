# Script to compute Figure 3 from the manuscript
# - Visualization of paleoclimate network of ACER AP and TraCE TSF, PR, TS on a map.
#   ACER network with all links and continent-wise aggregation, TraCE with aggregation
#   Code is similar to supplementary Figures S15 to S19 (see `create_figS15-19.R`).

#DONE & TESTED

## initialize data (skips automatically if done previously in the session, see `main_pcns.R`) ----
source('main_pcns.R')
source('main_pcn_meas.R')


## prepare node degree data & limits ----
ndg_list <- lapply(names(ACERtrc), function(nm) ACERtrc_ndg_raw %>% rename(ndg = !!sym(paste0('node_degree.',nm))) %>% select(window,site_id,ndg)) %>% 
  setNames(names(ACERtrc))

ndg_lims <- lapply(ndg_list, function(n) n %>% .$ndg) %>% unlist(use.names = FALSE)
ndg_lims <- c(min(ndg_lims),max(ndg_lims))


## full network plot for all raw ACER AP links using the generic network plot function of the `PCNdata` S3 class ----
plts <- lapply(1:length(ACERtrc), function(i) { #1:length(ACERtrc)
  nm <- names(ACERtrc)[i]
  if (nm == 'ACER_ap') dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data' else dat <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data'
  pltloc <- run_network_plot(ACERtrc[[i]], 
                             orig_data = dat, link_width = 0.3, ndg_lims = ndg_lims, leg_opt = 'cmb',
                             zoom = c(-180, -65, 180, 65), #frac_abs_strongest = 1, 
                             bundled = T, save_plot = list(activate = F), color_strength = T, split_regions_by_sign = T,
                             filter_windows = c('6-22'), node_degree = T, return = T, return_nw = T, bg = NULL)
}) %>% 
  setNames(names(ACERtrc))


## aggregated network plot for ACER AP and TraCE networks ----
grphs_agg <- lapply(names(ACERtrc), function(nm) {
  if (nm == 'ACER_ap') cnm <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data' else cnm <- 'corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data'
  testgraph <- create_bubble_network_from_corr_list(ACERtrc[[nm]][[cnm]])
}) %>% 
  setNames(names(ACERtrc))

width_lims <- lapply(grphs_agg, function(g) g %E>% as_tibble() %>% .$nlinks) %>% unlist(use.names = FALSE)
width_lims <- c(min(width_lims),max(width_lims))


site_data_list <- lapply(names(ACERtrc), function(nm) {
  ndg_list[[nm]] %>% full_join(ACERtrc[[nm]]$sites) %>% select(site_id,lat,long,ndg) %>% mutate(occurs = if_else(is.na(ndg) | ndg == 0,FALSE,TRUE),ndg = if_else(is.na(ndg),0,ndg))
}) %>% 
  setNames(names(ACERtrc))

plts_agg <- lapply(1:length(grphs_agg), function(i) {
  plot_bubble_network_spatial(grphs_agg[[i]],width_lims = width_lims,ndg_lims = ndg_lims,
                              leg_opt = 'cmb',
                              site_data = site_data_list[[i]])
}) %>% 
  setNames(names(ACERtrc))

## combine the 5 networks + 2 legends in one plot ----
leg1 <- get_legend(plts$ACER_ap$corr_nw + theme(legend.justification = 'center',
                                                legend.text = element_text(size = GLOBAL_FONT_SIZE + 4), legend.title = element_text(size = GLOBAL_FONT_SIZE + 8),
                                                legend.background = element_blank(), legend.spacing = unit(20, 'mm'))) #legend.box.background = element_rect(color = 'black', size = 0.8)))
                                               

leg2 <- get_legend(plts_agg$ACER_ap + theme(legend.justification = 'center', 
                                            legend.text = element_text(size = GLOBAL_FONT_SIZE + 4), legend.title = element_text(size = GLOBAL_FONT_SIZE + 8),
                                            legend.spacing = unit(11.1, 'mm')))

plt <- plot_grid(plts$ACER_ap$corr_nw + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
                 NULL,
                 NULL,
                 NULL,NULL,NULL,
                 plts_agg$ACER_ap + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
                 NULL,
                 plts_agg$TRACE_apsb + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
                 NULL,NULL,NULL,
                 plts_agg$TRACE_pr + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
                 NULL,
                 plts_agg$TRACE_ts + theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank()),
                 align = "hv", axis = "tblr", ncol = 3, 
                 labels = c('A | ACER AP LINKS', '', '', 
                            '','','',
                            'B | ACER AP', '', 'C | TraCE TSF', 
                            '','','',
                            'D | TraCE PR', '', 'E | TraCE TS'),
                 label_size = 24, label_fontface = 'bold', label_y = 1.0, hjust = 0, label_x = 0.09,
                 rel_heights = c(0.32,-0.022,0.32,-0.022,0.32),
                 rel_widths = c(0.48,-0.02,0.48)) + 
  draw_grob(leg1, x = 0.246, y = 0.383) + 
  draw_grob(leg2, x = 0.2465, y = 0.29) + 
  draw_grob(grid::rectGrob(gp = grid::gpar(fill = NA, lineheight = 10, fontsize = 20)), 
            x = 0.58, y = 0.745, width = 0.33, height = 0.19)


## save combined plot to disk ----
## width = 61 for non cairo_pdf
ggsave(plot = plt, width = 62, height = 42, units = 'cm', dpi = 200,
       filename = file.path(DIR_FIGURES, 'PCNagg_ACERapTRACEtsfprts_raw_6-22k_updnorm.pdf'), device = cairo_pdf)
ggsave(plot = plt, width = 61, height = 42, units = 'cm', dpi = 150,
       filename = file.path(DIR_FIGURES, 'PCNagg_ACERapTRACEtsfprts_raw_6-22k_updnorm.png'))
