# Script to compute Figures S28 and S29 from the supplement
# - sensitivity of cross-link fractions and node degrees on site type (S28)
#   and on ACER AP transformation (S29, analogous to Figure 5A from the manuscript)
#   and in complement to Figure 7, panel B from the manuscript
# - sensitivity of link density (S29) on the ACER AP transformation (sqrt, identity; 
#   analogous to Figure 4 from the manuscript)

# load/compute data
source('main_pcn_sensitivity.R')

# sensitivity of ACER AP network on site type
## NDG and CLF split between terrestrial and marine sites ----
## NOTE: need to initialize functions and objects from `create_fig_S5_S9_S10.R`!
sens_suff <- list(
  'raw' = list(c('MARI','TERR'))
)

trc_matrs_cmb <- plot_nw_clf_and_ndg_matr1x1_nosig(ACERsens_clp_raw,
                                                   ACERsens_ndg_raw,
                                                   #plts_sep = TRUE,
                                                   suffix_pairs = sens_suff, legend = T, leg_option = 'not_shared_internal')[[1]]

trc_matrs_cmb_a <- trc_matrs_cmb + 
  draw_plot(cntpltx, x = 0.4975, y = 0.1, width = 0.48, height = 0.48, scale = 1) + 
  draw_plot(cntplty, x = 0.3, y = 0.24, width = 0.6, height = 0.49, scale = 1) + 
  draw_plot_label('MARI', x = 0.49, y = 0.75, vjust = 0.97, size = GLOBAL_FONT_SIZE+12) + 
  draw_plot_label('TERR', x = 0.88, y = 0.27, vjust = 0.97, size = GLOBAL_FONT_SIZE+12)
  
ggsave(plot = trc_matrs_cmb_a, filename = file.path(DIR_FIGURES, 'CLF_NDG_ACERTERRMARI_raw_updnorm_correctedNA_woanom_compare.pdf'), 
       width = 54, height = 45, units = 'cm', device = cairo_pdf) # height = 45; 54 45

# sensitivity of ACER AP network on AP transformation
## plot and save link density evaluation ----
ld_lims_trf <- round(max(abs(lds_trf_raw$link_dens)),3) + 0.001
lds_trf_raw_plt <- lds_trf_raw %>% 
  mutate(ld_sgn = if_else(ld_sgn == 'ldp', '+', '-')) %>% 
  mutate(link_dens = if_else(ld_sgn == '-', -1*link_dens, link_dens)) %>% 
  mutate_at(vars(signal), toupper)
lds_trf_raw_plt$ld_sgn <- factor(lds_trf_raw_plt$ld_sgn, levels = c('+','-'))
lds_trf_raw_plt$signal <- factor(lds_trf_raw_plt$signal, levels = c('PROBIT', 'SQRT', 'IDENTITY'))

plt_ld_trf <- ggplot(mapping = aes(x = signal)) + #color = ld_sgn, 
  geom_col(data = lds_trf_raw_plt, mapping = aes(y = link_dens, fill = ld_sgn), width = 0.8) + 
  scale_fill_manual(values = rev(GLOBAL_PN_COLORS), guide = guide_legend('correlation sign', title.position = 'top')) + 
  scale_y_continuous(limits = c(-1*ld_lims_trf,ld_lims_trf)) + 
  labs(y = 'LD') + 
  global_title_and_axis() + 
  theme(legend.direction = 'horizontal', legend.position = 'bottom', axis.title.x = element_blank(), panel.grid.minor = element_blank())

ggsave(plt_ld_trf, filename = file.path(DIR_FIGURES, 'ACER_AP_LD_vs_trafo.pdf'),
       width = 16, height = 14, units = 'cm')

## NDG and CLF split between terrestrial and marine sites ----
## NOTE: need to initialize functions and objects from `create_fig_S5_S9_S10.R`!
trf_suff <- list(
  'raw' = list(c('sqrt','identity'))
)

matrs_cmb <- plot_nw_clf_and_ndg_matr1x1_nosig(ACERtrf_clp_raw,
                                                   ACERtrf_ndg_raw,
                                                   #plts_sep = TRUE,
                                                   suffix_pairs = trf_suff, legend = T, leg_option = 'not_shared_internal')[[1]]

matrs_cmb_a <- matrs_cmb + 
  draw_plot(cntpltx, x = 0.4975, y = 0.1, width = 0.48, height = 0.48, scale = 1) + 
  draw_plot(cntplty, x = 0.3, y = 0.24, width = 0.6, height = 0.49, scale = 1) + 
  draw_plot_label('SQRT', x = 0.49, y = 0.75, vjust = 0.97, size = GLOBAL_FONT_SIZE+12) + 
  draw_plot_label('IDENTITY', x = 0.82, y = 0.27, vjust = 0.97, size = GLOBAL_FONT_SIZE+12)

ggsave(plot = matrs_cmb_a, filename = file.path(DIR_FIGURES, 'CLF_NDG_ACERTRAFO_raw_updnorm_correctedNA_woanom_compare.pdf'), 
       width = 54, height = 45, units = 'cm', device = cairo_pdf) # height = 45; 54 45

plt_cmb_tm <- matrs_cmb_a + 
  draw_plot(plt_ld_trf + theme(legend.position = 'none', axis.title = element_text(size = GLOBAL_FONT_SIZE + 12), axis.text = element_text(size = GLOBAL_FONT_SIZE + 10)),
            x = -0.15, y = 0.25, width = 0.47, height = 0.47) + 
  draw_plot_label('A', x = -0.15, y = 0.764, size = GLOBAL_FONT_SIZE + 12) + 
  draw_plot_label('B', x = 0.35, y = 0.764, size = GLOBAL_FONT_SIZE + 12)

ggsave(plt_cmb_tm,
       filename = file.path(DIR_FIGURES, 'LD_CLF_NDG_ACERTRAFO_raw.pdf'), 
       width = 54, height = 45, units = 'cm')
