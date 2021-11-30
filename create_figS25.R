# Script to compute Figure S25 from the Supplement
# - Visualisation of the kappa coefficient of the TraCE contingency matrix

#DONE & TESTED

## initialize data (skips automatically if done previously in the session, see `main_emulation_pcns.R`, `main_emulation_pcn_meas.R`) ----
source('main_pcns.R')
source('main_emulation_pcns.R')
source('main_emulation_pcn_meas.R')

# evaluate contingency table of realized links between networks & kappa statistics of the tables for TraCE TSF x TraCE & Emulated LPJ ----
## computations for full adjacency matrix ----
### 1 ~ positive
### 0 ~ not significant
### -1 ~ negative
cmb_cl_trc <- list('TSF' = ACERtrc$TRACE_apsb$corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data,
                   'PR'= ACERtrc$TRACE_pr$corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data,
                   'TS'= ACERtrc$TRACE_ts$corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data,
                   'EM TSF ALL'= ACERtrc_em$forceall$corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data,
                   'EM TSF TS'= ACERtrc_em$forcets$corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data,
                   'EM TSF PR'= ACERtrc_em$forcepr$corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data,
                   'EM TSF CO2'= ACERtrc_em$forceco2$corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data)
nms_trc <- names(cmb_cl_trc)
ct_trc_glb <- lapply(1, function(i) {
  lapply(1:length(cmb_cl_trc), function(j) {
    ct <- compute_nw_contingency_table(cmb_cl_trc[[i]],cmb_cl_trc[[j]]) %>% 
      mutate(a_name = nms_trc[i], b_name = nms_trc[j])
  }) %>% 
    bind_rows()
}) %>% 
  bind_rows() %>% 
  mutate(kappa = purrr::map(.$cont_table, function(ct) {
    kp <- cohens_kappa(ct)
    k <- tibble(kappa = kp[1],kappal = kp[2], kappah = kp[3])
    return(k)
  }),
  cont_table = purrr::map(.$cont_table, function(ct) {
    ct <- as_tibble(ct)
    colnames(ct) <- c(-1,0,1)
    ct <- ct %>% add_column(.,adj.a = c('-1','0','1')) %>% gather(-adj.a,key = 'adj.b', value = 'nlink')
    return(ct)
  })) 
kp_trc_glb <- ct_trc_glb %>% 
  select(-cont_table) %>% 
  unnest(kappa) %>% 
  filter(!(b_name == 'TSF' & a_name == 'TSF'))
ct_trc_glb <- ct_trc_glb %>% 
  select(-kappa) %>% 
  unnest(cont_table)
ct_trc_glb$adj.b <- factor(ct_trc_glb$adj.b,levels = c('1','0','-1'))

## computations for significant links in both networks only ----
### 1 ~ positive
### -1 ~ negative
adj_vals2 <- c(-1,1)
ct_trc_sgo <- lapply(1, function(i) {
  lapply(1:length(cmb_cl_trc), function(j) {
    ct <- compute_nw_contingency_table(cmb_cl_trc[[i]],cmb_cl_trc[[j]],adj_vals = adj_vals2) %>% 
      mutate(a_name = nms_trc[i], b_name = nms_trc[j])
  }) %>% 
    bind_rows()
}) %>% 
  bind_rows() %>% 
  mutate(kappa = purrr::map(.$cont_table, function(ct) {
    kp <- cohens_kappa(ct)
    k <- tibble(kappa = kp[1],kappal = kp[2], kappah = kp[3])
    return(k)
  }),
  cont_table = purrr::map(.$cont_table, function(ct) {
    ct <- as_tibble(ct)
    colnames(ct) <- c(-1,1)
    ct <- ct %>% add_column(.,adj.a = c('-1','1')) %>% gather(-adj.a,key = 'adj.b', value = 'nlink')
    return(ct)
  })) 
kp_trc_sgo <- ct_trc_sgo %>% 
  select(-cont_table) %>% 
  unnest(kappa) %>% 
  filter(!(b_name == 'TSF' & a_name == 'TSF'))
ct_trc_sgo <- ct_trc_sgo %>% 
  select(-kappa) %>% 
  unnest(cont_table)
ct_trc_sgo$adj.b <- factor(ct_trc_sgo$adj.b,levels = c('1','-1'))


## plot kappa + se ----
kappa_lims <- c(min(c(kp_trc_glb$kappal,kp_trc_sgo$kappal)),max(c(kp_trc_glb$kappah,kp_trc_sgo$kappah)))
plot_nw_kappa_trc_glb <- ggplot(kp_trc_glb %>% mutate(b_name = str_replace(b_name, '-DGVM', '')), aes(x=b_name,y=kappa)) + 
  geom_errorbar(aes(ymin=kappal,ymax=kappah),color = GLOBAL_GREY_DARK, width = 0.1) + 
  geom_point(color = 'black', size = 3) + 
  #scale_color_manual(values=c(GLOBAL_PN_COLORS[1],'black',GLOBAL_PN_COLORS[2])) +
  global_title_and_axis() + 
  scale_y_continuous(limits = kappa_lims) + 
  labs(y=expression(kappa),subtitle = 'A | All record pairs') + 
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank(),
        axis.title.y = element_text(family = GLOBAL_FONT_FAMILY,face = GLOBAL_FONT_FACE_TITLE,size = GLOBAL_FONT_SIZE+14))

plot_nw_kappa_trc_sgo <- ggplot(kp_trc_sgo %>% mutate(b_name = str_replace(b_name, '-DGVM', '')), aes(x=b_name,y=kappa)) + 
  geom_errorbar(aes(ymin=kappal,ymax=kappah),color = GLOBAL_GREY_DARK, width = 0.1) + 
  geom_point(color = 'black', size = 3) + 
  #scale_color_manual(values=c(GLOBAL_PN_COLORS[1],'black',GLOBAL_PN_COLORS[2])) +
  global_title_and_axis() + 
  scale_y_continuous(limits = kappa_lims) + 
  labs(y=expression(kappa),subtitle = 'B | Significant pairs') + 
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank(),
        axis.title.y = element_text(family = GLOBAL_FONT_FAMILY,face = GLOBAL_FONT_FACE_TITLE,size = GLOBAL_FONT_SIZE+14))

## combined ----
plt_cmb <- plot_grid(plot_nw_kappa_trc_glb,plot_nw_kappa_trc_sgo + theme(legend.position = 'none'),
                     #labels = c('A | All record pairs', 'B | Significant pairs'), 
                     #label_size = 24, label_fontface = 'bold', label_y = 1.0,
                     align = 'hv',axis = 'tblr',nrow = 2)
#plt_cmb
ggsave(file.path(DIR_FIGURES,paste0('kappa_LPJ_comb.pdf')),width = 32,height = 22,units = 'cm')

