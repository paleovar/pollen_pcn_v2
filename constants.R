library(viridisLite)


# file paths
STANDARD_WD <- '~/Documents/code/pollen_pcn_v2'

DIR_DATASETS <- './data_ACER/ACER_pollen_charcoal_2017-02-15_flat/all database tables' # previously './data/ACER_pollen_charcoal_2017-02-15_flat/all database tables'
DATAFILES_TYPE <- 'csv'

DIR_EVALUATIONS <- './data_out' # previously './eval' (rarely used, mostly superseded by `DIR_CACHE`)

DIR_DATASETS_SUPPLEMENTARY <- './data_in' # previously './data_supplementary'
DIR_PNV <- file.path(DIR_DATASETS_SUPPLEMENTARY, 'global potential vegetation distribution')
NAME_POT_VEG_DISTR <- 'PNV.MLR'

DIR_MODEL_DATA <- '/stacydata/data' #file.path(DIR_DATASETS_SUPPLEMENTARY)

NAT_EARTH_DATA_PATH <- file.path(DIR_DATASETS_SUPPLEMENTARY, 'naturalearthdata')
# for further projections see https://proj4.org/operations/projections/index.html
GLOBAL_CRS <- list(fixed = 'fixed', robinson = '+proj=robin', wintri = '+proj=wintri', azequidistant = '+proj=aeqd', area = '+proj=aea +lat_1=29.5 +lat_2=42.5')

# for time series from ACER
GLOBAL_SERIES_SIGNALS <- list(arboreal_pollen = list(signal_name = 'pcnt_arb_pollen', signal_id = 'arb_pollen_data_id'),
                              pollen = list(signal_name = 'taxon_pcnt', signal_id = 'pollen_data_id', signal_grouper = 'taxon'), 
                              harm_pollen = list(signal_name = 'taxon_pcnt', signal_id = 'harm_pollen_data_id', signal_grouper = 'taxon_harmonized'), 
                              pca_pollen = list(signal_name = 'pc_value', signal_id = 'pca_id', signal_grouper = 'pc'), 
                              hybrid_ice_core = list(signal_name = 'hic', signal_id = 'hic_id'))
GLOBAL_SERIES_WINDOWS <- list(0, 10000, 20000, 30000, 40000)
GLOBAL_SERIES_WINDOWS_LABELS <- list('0-10', '10-20', '20-30', '30-40')
GLOBAL_SERIES_DETREND_LABELS <- list(linear = 'linear detrend', gaussbp = 'gaussbandpass detrend')
GLOBAL_SERIES_WINDOWS_ <- list(0, 25000, 50000, 100000)
GLOBAL_SERIES_WINDOWS_LABELS_ <- list('0-25', '25-50', '50-100')
GLOBAL_SERIES_WINDOWS__ <- list(0, 25000, 65000, 110000)
GLOBAL_SERIES_WINDOWS_LABELS__ <- list('0-25', '25-65', '65-110')

DIR_FIGURES <- './figures'

DIR_CACHE = './cache'

# DIR_MEGABIOMES <- file.path(DIR_DATASETS_SUPPLEMENTARY, 'harmonization/megabiomes')
DIR_TAXA <- file.path(DIR_DATASETS_SUPPLEMENTARY, 'harmonization/taxa')
DIR_ARB_POLLEN <- file.path(DIR_DATASETS_SUPPLEMENTARY, 'harmonization/arboreal_pollen')
# null model data, was DIR_HIC_DATA <- file.path(DIR_DATASETS_SUPPLEMENTARY, 'hybrid_icecores')
DIR_NM_DATA <- file.path(DIR_DATASETS_SUPPLEMENTARY, 'nullmodel')

DIR_CORRECTED_DATA <- file.path(DIR_DATASETS_SUPPLEMENTARY, 'corrected_data')

# number of cores used for null model/PCN computation for parallelized code parts
# set to `NCORES = NULL` for automatic parallelization
NCORES_PCN <- 4

# coloured cat printing
cr <- crayon::yellow

# plotting options
GLOBAL_FONT_FACE_TITLE <- 'bold'
GLOBAL_FONT_FACE_TEXT <- 'plain'
GLOBAL_FONT_SIZE <- 17
GLOBAL_FONT_FAMILY <- 'sans'

GLOBAL_BLUE_LIGHT <- viridis(1, begin = 0.6)
GLOBAL_BLUE_LIGHT_ALPHA_LOW <- viridis(1, alpha = 0.45, begin = 0.6)
GLOBAL_BLUE_DARK <- viridis(1, begin = 0.4)
GLOBAL_GREEN_LIGHT <- viridis(1, begin = 0.9)
GLOBAL_GREEN_LIGHT_ALPHA_LOW <- viridis(1, alpha = 0.45, begin = 0.9)
GLOBAL_GREEN_DARK <- viridis(1, begin = 0.7)
GLOBAL_GREY_LIGHT <- grey.colors(1, start = 0.97, end = 0.99)
GLOBAL_GREY_LIGHT_ALPHA_LOW <- grey.colors(1, start = 0.97, end = 0.99, alpha = 0.45)
GLOBAL_GREY_DARK <- grey.colors(1, start = 0.6)
GLOBAL_GREY_MEDIUM <- grey.colors(1, start = 0.88)
GLOBAL_GREY_MEDIUM_LIGHT <- grey.colors(1, start = 0.95)
GLOBAL_RED_DARK <- viridis(1, begin = 0.6, option = 'B')
GLOBAL_RED_LIGHT <- viridis(1, begin = 0.7, option = 'B')
GLOBAL_BLUE_DDARK <- viridis(1, begin = 0.2, option = 'E')
GLOBAL_BLUE_DLIGHT <- viridis(1, begin = 0.4, option = 'E')
GLOBAL_ORANGE_DARK <- viridis(1, begin = 0.75, option = 'A')
GLOBAL_ORANGE_LIGHT <- viridis(1, begin = 0.85, option = 'A')
GLOBAL_YELLOW_DARK <- viridis(1, begin = 0.9, option = 'E')
GLOBAL_YELLOW_LIGHT <- viridis(1, begin = 0.98, option = 'E')
GLOBAL_VIOLET_DARK <- viridis(1, begin = 0.4, option = 'C')
GLOBAL_VIOLET_LIGHT <- viridis(1, begin = 0.5, option = 'C')
GLOBAL_PN_COLORS <- rev(c('#ED3537', '#1188d6'))

LABELED_SITES_FOR_LAT_ORDERED_PLOT <- c(10, 7, 19, 90, 29, 82, 17, 38, 28, 15)
SAMPLING_RATE_INTERVALS_SIMPLE <- list(list(start = 0, stop = 8000), list(start = 15000, stop = 73000))
SAMPLING_RATE_INTERVALS_DEFAULT <- list(list(start = 0, stop = 8000), list(start = 8001, stop = 18999),
                                        list(start = 19000, stop = 27000))
SAMPLING_RATE_INTERVALS_10k_WINDOW <- list(list(start = 0, stop = 10000), list(start = 10000, stop = 20000),
                                           list(start = 20000, stop = 30000), list(start = 30000, stop = 40000), 
                                           list(start = 40000, stop = 50000), list(start = 50000, stop = 60000), 
                                           list(start = 60000, stop = 70000))


# datasets constants
OLDEST_DATAPOINT_AGE <- 250955
BIOME_PERCENTAGES_ABBR <- data.frame(biome = c('Bo forest', 'Gr', 'Sav', 'SE Pine forest', ' Subtr forest', 'Te forest', 'Te mountain forest', 'Tr forest', 'WTe forest'), 
                                     biome_long = c('Boreal forest', 'Grass- and dry shrublands', 'Savannah', 'Southeastern pine forest', 'Subtropical forest',
                                                    'Temperate forest', 'Temperate mountain forest', 'Tropical forest', 'Warm temperate forest'))

REGIONS <- tibble(site_id = c(99,33,48,32,35,31,46,34,1,19,76,73,8,18,79,40,3,92,37,25,36,64,17,77,12,23,21,20,38,67,71,97,72,87,72.97,14,13,94,88,86,9,16,26,11,29,4,95,28,89,52,51,56,50,15,10,54,91,58,74,27,7,55,42,6,90,98,90.98,53,47,57,59,2,39,82,41,24,96,30,43,62,75,60,44,45,80,85,100,85.100,83,84,93,49,61,81,63,69,68,66,70), 
                  region  = c(rep('Europe', 22), rep('Africa', 13), rep('Asia', 11), rep('Australia', 8), rep('NAmerica', 23), rep('SAmerica', 22)))
REGION_LOCS <- tibble(region = c('Europe', 'Africa', 'Asia', 'Australia', 'NAmerica', 'SAmerica', 'earth'), 
                      lat = c(40, -50, 60, -55, 60, -40, 0), 
                      long = c(7, 28, 113, 143, -110, -48, 0))
REGION_LOCS_pn <- tibble(region = c('Europe.p', 'Africa.p', 'Asia.p', 'Australia.p', 'NAmerica.p', 'SAmerica.p', 'earth.p', 'Europe.n', 'Africa.n', 'Asia.n', 'Australia.n', 'NAmerica.n', 'SAmerica.n', 'earth.n'), 
                         lat = c(44, -50, 60, -55, 60, -40, 65, 32, -56, 59, -56, 59, -41, -60), 
                         long = c(7, 28, 113, 143, -110, -48, 0, -10, -25, 150, 150, -170, -160, 30))
REGION_LOCS_BUBBLE <- tibble(region = c('Europe', 'Africa', 'Asia', 'Australia', 'NAmerica', 'SAmerica', 'earth'), 
                             lat = c(50, -10, 40, -32, 40, -20, 0), 
                             long = c(9, 28, 125, 150, -120, -60, 0))

REGIONS <- bind_rows(REGIONS, tibble(site_id = c(101, 102), region = c('Antarctica', 'Greenland')))
REGION_LOCS <- bind_rows(REGION_LOCS, tibble(region = c('Antarctica', 'Greenland'), lat = c(-60, 60), long = c(90, -48)))
REGION_LOCS_pn <- bind_rows(REGION_LOCS_pn, tibble(region = c('Antarctica.p', 'Greenland.p', 'Antarctica.n', 'Greenland.n'), lat = c(-60, 60, -60, 60), long = c(90, -44, 90, -52)))

ZOOM_ATLANTIC <- c(c(-141,-56), c(30,56))
ZOOM_EUROPE <- c(c(-15,30), c(40,50))

NDG_RANGE <- c(3,9) # scale of node degrees

# time windows of available & used model data
WINDOWS_TRACE <- list(6000, 22000)
WINDOW_LABELS_TRACE <- list('6-22')
WINDOWS_LC <- list(6200, 18000)
WINDOW_LABELS_LC <- list('6.2-18')
WINDOWS_BBC <- list(6000, 22000)#list(0, 21000)
WINDOW_LABELS_BBC <- list('6-22')#list('0-21')
