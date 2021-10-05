# Script to compute Figure 4 from the manuscript
# - Link densities (LD) of paleoclimate networks of ACER AP and simulated TSF, TS, PR
#   for the raw signal and the signal filtered for variability on millennial time scales.
# - Link densities of the emulated LPJ-DGVM (TraCE) and TRIFFID (HadCM3) networks. 

# initialize data and compute/re-load network data
source('main_pcn_meas.R')
source('main_emulation_pcn_meas.R')

#DONE & TESTED

## processing of link densities for model, emulated, and ACER networks ----

# compile data frames for all models, ACER, emulations from data computed
# in `main_pcn_meas.R`, `main_emulation_pcns.R`
lds <- bind_rows(list(lds_trc_raw,lds_lc_raw,lds_hcm_raw,lds_trc_mil,lds_lc_mil)) %>% # lds_trc_milext
  mutate(ld = if_else(ld_sgn == 'ldp', link_dens, -1*link_dens)) %>% 
  filter(data == 'TRACE' | signal != 'ACER_ap') %>% 
  select(tscale,signal,data,ld_sgn,ld) %>% 
  mutate(data = if_else(signal == 'ACER_ap', 'ACER', data)) %>% 
  rowwise() %>% 
  mutate(signal = if_else(str_detect(signal,'_'), str_split(signal, '_')[[1]][2], signal)) %>% 
  mutate(signal = if_else(signal == 'apsb', 'ap', signal)) %>% 
  arrange(signal,factor(data, levels = c('ACER','TRACE','BBC','LCLIM'))) %>% 
  mutate(tscale = str_to_upper(tscale), signal = str_to_upper(signal))

lds_em <- bind_rows(list(lds_trc_em,lds_hcm_em)) %>%
  mutate(ld = if_else(ld_sgn == 'ldp', link_dens, -1*link_dens)) %>% 
  select(tscale,signal,data,ld_sgn,ld) %>% 
  rowwise() %>% 
  mutate(signal = if_else(str_detect(signal,'_'), str_split(signal, '_')[[1]][2], signal)) %>% 
  mutate(signal = if_else(signal == 'apsb', 'ap', signal)) %>% 
  arrange(signal,factor(data, levels = c('TRACE','BBC'))) %>% 
  mutate(tscale = str_to_upper(tscale), signal = str_to_upper(signal))

lds <- lds %>% 
  mutate(signal = if_else(signal == 'AP' & data != 'ACER', 'TSF', signal), 
         data = case_when(data == 'TRACE' ~ 'TraCE',
                          data == 'BBC' ~ 'HadCM3',
                          data == 'LCLIM' ~ 'LOVECLIM',
                          TRUE ~ data)) %>% 
  mutate(tscale = case_when(tscale == 'RAW' ~ 'A | RAW',
                            tscale == 'MIL' ~ 'C | MIL')) %>% 
  mutate(ld_sgn = case_when(ld_sgn == 'ldp' ~ '+',
                            ld_sgn == 'ldn' ~ '-')) %>% 
  # add dummy rows for plotting "window" in millennial to allow sharing axis
  bind_rows(tibble(tscale = rep('C | MIL',3),
                   signal = c('TSF', 'PR', 'TS'),
                   data = rep('HadCM3',3),
                   ld_sgn = rep('+',3),
                   ld = rep(0,3))) 
lds$data <- factor(lds$data, levels = c('ACER','TraCE','HadCM3','LOVECLIM'))

lds_em <- lds_em %>% 
  mutate(signal = if_else(signal == 'AP' & data != 'ACER', 'TSF', signal), 
         data = case_when(data == 'TRACE' ~ 'TraCE',
                          data == 'BBC' ~ 'HadCM3',
                          data == 'LCLIM' ~ 'LOVECLIM',
                          TRUE ~ data)) %>% 
  filter(signal %in% toupper(c('forceall','forcets','forcepr','forceco2'))) %>% 
  mutate(signal = str_split(signal, 'FORCE')[[1]][2]) %>% 
  mutate(tscale = case_when(tscale == 'RAW_EM' ~ 'B | RAW EM')) %>% 
  mutate(ld_sgn = case_when(ld_sgn == 'ldp' ~ '+',
                            ld_sgn == 'ldn' ~ '-'))

# define plot function for bar plot ----
# function has different options for depicting the LD bars and confidence level
# for the arguments used for the Figure in the manuscript see below
plot_nw_link_dens_bars <- function(dat, order = 'source', conf_bar = NULL, scales = 'free_x', ylims = NULL, leg = FALSE) {
  
  if (!is.null(conf_bar)) {
    dat <- dat %>% 
      mutate(cbar = conf_bar * abs(ld))
  }
  
  dat$tscale <- factor(dat$tscale, levels = unique(dat$tscale))
  dat$signal <- factor(dat$signal, levels = unique(dat$signal))
  if(!is.factor(dat$data)) dat$data <- factor(dat$data, levels = unique(dat$data))
  dat$ld_sgn <- factor(dat$ld_sgn, levels = unique(dat$ld_sgn))
  
  if (order == 'source') {
    mapp <- aes(fill = ld_sgn, y = ld, x = signal)
    erb <- geom_errorbar(aes(x = signal, ymin = ld - cbar, ymax = ld + cbar), size = 0.5, width = 0.1)
    fct <- facet_grid(cols = vars(data), rows = vars(tscale), space = 'free_x', scales = scales, switch = 'x')
  } else if (order == 'signal') {
    mapp <- aes(fill = ld_sgn, y = ld, x = data)
    erb <- geom_errorbar(aes(x = data, ymin = ld - cbar, ymax = ld + cbar), size = 0.5, width = 0.1)
    fct <- facet_grid(cols = vars(signal), rows = vars(tscale), space = 'free_x', scales = scales, switch = 'x')
  }
  
  plt <- ggplot(dat, mapp) + 
    geom_col(width = 0.75) 
  
  if (leg) {
    plt <- plt + 
      geom_hline(data = tibble(y = c(-0.025,0.025), lt = rep('LD = 0.025',2)),
                 mapping = aes(yintercept = y, linetype = lt),
                 color = 'grey20', alpha = 0.65, size = 0.8) + 
      scale_linetype_manual(values = 'dotted',
                            guide = guide_legend(title = '', direction = 'horizontal', title.position = 'left',
                                                 override.aes = list(size = 2,linetype="11"), order = 2)) +
      scale_fill_manual(values = c('#ED3537', '#1494E9'), 
                        guide = guide_legend(title = 'SIGN', direction = 'horizontal', title.position = 'left', order = 1,
                                             label.theme = element_text(size = GLOBAL_FONT_SIZE + 2, face = 'bold'))) # title = 'SIGN'
  } else {
    plt <- plt + 
      scale_fill_manual(values = c('#ED3537', '#1494E9'), guide = FALSE) + #c('red2','blue3'), c('#ff4040', '#499bd1'), c('#D9514EFF', '#2DA8D8FF')
      geom_hline(yintercept = c(-0.025,0.025), color = 'grey20', alpha = 0.65, linetype = 'dotted', size = 0.8)
  }
  
  if (!is.null(conf_bar)) {
    plt <- plt + 
      erb
  }
  
  plt <- plt +
    fct + 
    global_title_and_axis() + 
    labs(y = 'LD') + 
    #coord_cartesian(ylim = c(max(dat$ld),min(dat$ld)), expand = T, clip = 'off') + 
    theme(panel.grid.major.x = element_blank(),
          panel.spacing.x = unit(0, "lines"),
          #axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = GLOBAL_FONT_SIZE + 2),
          axis.text = element_text(size = GLOBAL_FONT_SIZE + 2),
          strip.background.x = element_blank(), strip.background.y = element_blank(),
          strip.text.x = element_text(hjust = 0.5, size = GLOBAL_FONT_SIZE + 2),
          strip.text.y = element_text(angle = 0, vjust = 1, size = GLOBAL_FONT_SIZE + 2),
          strip.placement = 'outside', 
          legend.box = 'horizontal',
          legend.text = element_text(size = GLOBAL_FONT_SIZE + 2),
          legend.title = element_text(size = GLOBAL_FONT_SIZE + 2),
          legend.background = element_blank(), legend.box.background = element_rect())
    
  if (!is.null(ylims)) {
    plt <- plt + 
      ylim(ylims)
  }
  
  if (order == 'signal') {
    plt <- plt + 
      theme(axis.text = element_text(size = GLOBAL_FONT_SIZE - 6))
  }
  return(plt)
}

# execute plots, combine, and save the figure to disk ----
ld_lims <- c(min(c(lds$ld,lds_em$ld)),max(c(lds$ld,lds_em$ld))) # limits for common y axis
plt_lds_raw <- plot_nw_link_dens_bars(lds %>% filter(tscale == 'A | RAW'),ylims = ld_lims) + 
  theme(strip.text.x = element_blank(),strip.text.y = element_blank())#,axis.text.x = element_blank(), axis.ticks.x = element_blank())
#plt_lds_raw
plt_lds_mil <- plot_nw_link_dens_bars(lds %>% filter(tscale == 'C | MIL'),ylims = ld_lims) + 
  theme(strip.text.y = element_blank())
#plt_lds_mil
plt_lds_em <- plot_nw_link_dens_bars(lds_em,ylims = ld_lims,leg = TRUE) + 
  theme(axis.title.y = element_blank(),strip.text.y = element_blank(),axis.ticks.y = element_blank(), axis.text.y = element_blank())
#plt_lds_em
plt_lds_leg <- get_legend(plt_lds_em)
plt_lds_em <- plt_lds_em + theme(legend.position = 'none') 

plt_lds_cmb <- plot_grid(NULL,NULL,NULL,
                         plt_lds_raw,NULL,plt_lds_em,
                         NULL,NULL,NULL,
                         plt_lds_mil,#plt_lds_leg,
                         align = 'hv', ncol = 3, axis = 'tb',
                         labels = c('','','',
                                    'A | RAW','','',#need to reposition 'C | RAW EM' below
                                    '','','',
                                    ' C | MIL'),
                         label_size = GLOBAL_FONT_SIZE + 2,
                         hjust = -0.83, vjust = -0.05,
                         rel_heights = c(0.03,0.4,-0.01,0.4),
                         rel_widths = c(0.5,0.0,0.4)) + 
  draw_grob(plt_lds_leg, x = 0.6775, y = 0.33, width = 0.2, height = 0.15) + 
  draw_text(x = 0.621, y = 0.976, size = GLOBAL_FONT_SIZE + 2, fontface = 'bold', text = 'B | RAW EM')


ggsave(plot = plt_lds_cmb,
       filename = file.path(DIR_FIGURES, 'LD_raw_mil2-8ka_emul_compare_alt2_updnorm_woVC_fixnm.pdf'),
       units = 'cm', width = 36, height = 20, device = cairo_pdf)#width = 20, height = 11)
