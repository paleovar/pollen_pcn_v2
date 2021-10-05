# Script to compute Figure S2 from the Supplement
# - uncategorized fraction through arboreal pollen classification based on harmonized taxa
# and to compute some related numbers that are explicitly reported in the manuscript/Supplement

# load/compute data
source('main_pcns.R')

# prepare data: median removed fraction taxa in 6-22 ka BP through unclassified arboreal pollen ----
ACERhere <- PCNdata() %>% 
  mix_sample_dating_and_err() %>% 
  harmonized_pollen_to_db() %>% 
  arboreal_pollen_to_db()

taxa_cat <- suppressMessages(read_xlsx(path = paste(DIR_ARB_POLLEN,
                                                    paste(paste('new_matchfile', 'harm_taxa', 2, 'arb_pollen_match', 0.01, 'min_occur_frac', sep = '_'), 'xlsx', sep = '.'), sep = '/'))) %>% 
  select(taxon_harmonized, is_arboreal_pollen)

sites_in_nws <- unique(c(ACERtrc$ACER_ap$corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data$from,
                         ACERtrc$ACER_ap$corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data$to))

ACERrem_dat <- ACERhere$harm_pollen_data %>% 
  inner_join(taxa_cat) %>% 
  select(-taxon_pcnt, -harm_pollen_data_id, -is_arboreal_pollen) %>% 
  group_by(site_id,sample_id) %>%
  summarise(count_cat_sum = sum(taxon_count, na.rm = TRUE), .groups = 'keep') %>% 
  full_join(ACERhere$harm_pollen_data %>% 
              select(-taxon_pcnt, -harm_pollen_data_id) %>% 
              group_by(site_id,sample_id) %>% 
              summarise(count_harm_sum = sum(taxon_count, na.rm = TRUE), .groups = 'keep')) %>% 
  mutate(frac_not_cat = (count_harm_sum - count_cat_sum)/count_harm_sum) %>% 
  inner_join(ACERhere$sample_dating %>% select(sample_id, mixed_age)) %>% 
  filter(mixed_age > WINDOWS_TRACE[[1]] & mixed_age < WINDOWS_TRACE[[2]]) %>% 
  filter(site_id %in% sites_in_nws) %>% 
  group_by(site_id) %>% 
  mutate(med_fnc = median(frac_not_cat, na.rm = TRUE))

site_63 <- filter_sites(ACERhere$sites, ACERhere$sample_dating, T, 'all','all')$site_id
ACERrem_dat_63 <- ACERhere$harm_pollen_data %>% 
  inner_join(taxa_cat) %>% 
  select(-taxon_pcnt, -harm_pollen_data_id, -is_arboreal_pollen) %>% 
  group_by(site_id,sample_id) %>%
  summarise(count_cat_sum = sum(taxon_count, na.rm = TRUE), .groups = 'keep') %>% 
  full_join(ACERhere$harm_pollen_data %>% 
              select(-taxon_pcnt, -harm_pollen_data_id) %>% 
              group_by(site_id,sample_id) %>% 
              summarise(count_harm_sum = sum(taxon_count, na.rm = TRUE), .groups = 'keep')) %>% 
  mutate(frac_not_cat = (count_harm_sum - count_cat_sum)/count_harm_sum) %>% 
  inner_join(ACERhere$sample_dating %>% select(sample_id, mixed_age)) %>% filter(mixed_age > WINDOWS_TRACE[[1]] & mixed_age < WINDOWS_TRACE[[2]]) %>% 
  filter(site_id %in% site_63) %>% 
  ungroup()

print(ACERrem_dat_63 %>% summarise(med_fnc = mean(frac_not_cat, na.rm = TRUE))) # mean over all samples (selected & time window); reported on manuscript
print(ACERrem_dat_63 %>% group_by(site_id) %>% mutate(med_fnc = median(frac_not_cat, na.rm = TRUE)) %>% distinct(site_id,med_fnc) %>% filter(med_fnc == 0)) # none lost
print(ACERrem_dat_63 %>% group_by(site_id) %>% mutate(med_fnc = median(frac_not_cat, na.rm = TRUE)) %>% summarise(mean_med_fnc = mean(med_fnc)) %>% arrange(desc(mean_med_fnc))) # highest uncat fractions (all 63)
print(ACERrem_dat %>% summarise(mean_med_fnc = mean(med_fnc))%>% arrange(desc(mean_med_fnc))) # highest uncat fractions (all 58)

### Figure S2: uncategorized fraction of all samples in 6-22 ka BP & median uncategorized fraction for all 58 records of ACER AP network ----
ACERrem_dat_plt <- ACERrem_dat %>% 
  ungroup() %>% 
  nest_by(site_id) %>% 
  arrange(site_id)
ACERrem_dat_plt$xid <- 1:length(ACERrem_dat_plt$site_id)
site_labs_x <- ACERrem_dat_plt$site_id %>% as.character() %>% str_replace('\\.', '\n\\/') %>% str_replace('\n\\/1','\n\\/100')
ACERrem_dat_plt$site_lab <- factor(site_labs_x, levels = site_labs_x)

ACERrem_dat_plt <- ACERrem_dat_plt %>%
  unnest(data) %>% 
  select(-count_cat_sum, -count_harm_sum) %>% 
  gather(key = 'key', value = 'val', frac_not_cat, med_fnc) %>% 
  distinct()

plt_cat <- ggplot(mapping = aes(x = xid, y = val, color = key, shape = key, group = site_id)) + 
  geom_point(data = ACERrem_dat_plt %>% filter(key == 'frac_not_cat') %>% mutate(key = 'sample') %>% distinct(), 
             alpha = 0.3, size = 2.5) + 
  geom_point(data = ACERrem_dat_plt %>% filter(key == 'med_fnc') %>% mutate(key = 'median') %>% distinct(site_id,key,val,xid,site_lab),
             alpha = 1, size = 3.5, stroke = 1) + 
  scale_color_manual(values = rev(c('black', 'red')), guide = guide_legend(override.aes = list(size = 5))) + 
  #scale_fill_manual(values = rev(c('black', 'red'))) + 
  scale_shape_manual(values = rev(c(21,23))) + 
  scale_x_continuous(breaks = ACERrem_dat_plt$xid, labels = ACERrem_dat_plt$site_lab, expand = c(0.01,0.01)) + 
  global_title_and_axis() + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), legend.title = element_blank(),
        legend.justification = c(1,1), legend.position = c(0.99,0.98)) + 
  labs(x = 'site id', y = 'uncategorized fraction')

ggsave(plt_cat, filename = file.path(DIR_FIGURES, 'uncat_tax_frac.pdf'), width = 40, height = 16, units = 'cm')


### for plotting the median uncategorized fraction on a map (not used in the manuscript/supplement) ----
#plot_ACER_sites_catfrac_map <- function(dat) {
#  projection <- 'robinson'
#  zoom <- c(-180, -60, 180, 75) 
#  graticules <- 30
#  
#  sites <- arrange(dat, site_id)
#  sites_transformed <- as_tibble(project(cbind(sites$long, sites$lat), proj = GLOBAL_CRS[[projection]])) %>% 
#    rename(long_new = V1, lat_new = V2) %>% 
#    bind_cols(select(sites, site_id))
#  
#  sites <- sites %>% 
#    full_join(sites_transformed,
#              by = 'site_id') %>% 
#    select(-lat, -long) %>% 
#    rename(long = long_new, lat = lat_new)
#  if (!is.null(zoom)) {
#    yaxislclip <- zoom[2]; xaxislclip <- zoom[1]; xaxisrclip <- zoom[3]
#  } else {
#    yaxislclip <- -90; xaxislclip <- -180; xaxisrclip <- 180
#  }
#  
#  data <- inner_join(sites, dat %>% select(site_id, catfrac))
#  
#  map_plot <- shapefile_map_plot_(ggobj = ggplot(), projection = projection, grat_labs = list(x = c(-180, -90, 0, 90, 180), y = c(-90, -60, -30, 0, 30, 60, 90)),
#                                  graticules = graticules, 
#                                  lgm_ice = F, cb_friendly = TRUE, zoom = zoom,
#                                  yaxislclip = yaxislclip, xaxisrclip = xaxisrclip, xaxislclip = xaxislclip)
#  map_plot <- map_plot + 
#    geom_point(data = data,
#               mapping = aes(x = long, y = lat,
#                             fill = catfrac),
#               size = 5,
#               shape = 21,
#               alpha = I(0.6),
#               stroke = 1) +
#    scale_fill_distiller(palette = 'Spectral', guide = guide_colorsteps(title = 'fraction uncategorized', order = 3, show.limits = T, ticks = T, direction = 'horizontal',
#                                                                        barheight = 0.75, barwidth = 25, title.position = 'top')) + 
#    theme(panel.background = element_blank(), 
#          axis.text = element_blank(), 
#          panel.grid = element_blank(), 
#          axis.ticks = element_blank(), 
#          panel.border = element_blank())
#  
#  return(map_plot)
#}
#
#ggsave(plot_ACER_sites_catfrac_map(ACERrem_dat_63 %>% filter(site_id %in% sites_in_nws) %>% group_by(site_id) %>% 
#                                     mutate(med_fnc = median(frac_not_cat, na.rm = TRUE)) %>% distinct(site_id,med_fnc) %>% 
#                                     inner_join(ACERhere$sites) %>% rename(catfrac = med_fnc)),
#       filename = file.path(DIR_FIGURES,'catfrac_map_58.pdf'),width = 34, height = 24, units = 'cm')
