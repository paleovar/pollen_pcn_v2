# helpers for taxa harmonization

# taxa harmonization --------------------
#' @rdname acer.harmonise.taxa
#' @param pollen_data tibble, ACER pollen data (sample wise)
#' @param harmonization_level integer, either 1 or 2 (1 clips all words depite the first one, 2 does the same with separation with slashes considered)
#'
#' @return tibble, conversion table
#'
#' @examples
harmonize_taxa_make_conversion_table <- function(pollen_data, harmonization_level = 2) {
  taxa <- select(pollen_data, taxon) %>% 
    distinct() %>% 
    arrange(taxon) %>% 
    rowwise() %>% 
    mutate(taxon_harmonized = if_else(str_split(taxon, ' ')[[1]][1] == 'cf',
                                      str_split(taxon, ' ')[[1]][2],
                                      str_split(taxon, ' ')[[1]][1])) %>% 
    {if (harmonization_level == 1) {.}else{transmute(., taxon = taxon, taxon_harmonized = str_split(taxon_harmonized, '/')[[1]][1])}} %>% 
    rownames_to_column(var = 'taxa_conversion_id')
  return(taxa)
}


#' @param harmonization_level integer, either 1 or 2 (1 clips all words despite the first one, 2 does the same with separation with slashes considered)
#' @param save_conversion_to_file logical
#'
#' @return tibble, ACER pollen data with harmonised taxonomy
#' @export
#'
#' @examples
harmonize_taxa_ax <- function(pollen_data,
                              harmonization_level = 2,
                              save_conversion_to_file = FALSE) {
  # level 2 harmonization stronger than level 1
  taxa <- harmonize_taxa_make_conversion_table(pollen_data = pollen_data, harmonization_level = harmonization_level)
  
  if (save_conversion_to_file) {
    library(writexl)
    write_xlsx(x = taxa,
               path = paste(DIR_TAXA, paste(paste('harmonized_taxa', harmonization_level, sep = '_'), 'xlsx', sep = '.'), sep = '/'),
               col_names = TRUE, format_headers = TRUE)
  }
  
  pollen_data <- inner_join(pollen_data, 
                            taxa,
                            by = 'taxon') %>% 
    select(-pollen_data_id, -taxon) %>% 
    group_by(site_id, sample_id, taxon_harmonized) %>% 
    summarise(taxon_count = sum(taxon_count, na.rm = TRUE), taxon_pcnt = sum(taxon_pcnt, na.rm = TRUE)) %>% 
    #group_by(site_id, sample_id) %>% 
    rownames_to_column(var = 'harm_pollen_data_id')
  
  return(pollen_data)
}


# arboreal pollen matching ---------------
# optimize_taxa_for_arb_pollen_by_taxa_minimal_occurence_fraction is to prepare and write the relevant individual taxa for arboreal pollen classification
# match_top_taxa_to_harmonized_arb_pollen_by_min_occ_frac is to read a "matchfile" i.e. the arboreal pollen classification
#' @rdname acer.harmonise.taxa
#' @param pollen_data_harm tibble, harmonised ACER pollen data
#' @param harmonisation_level integer, either 1 or 2 (1 clips all words depite the first one, 2 does the same with separation with slashes considered)
#' @param minimal_occurence_fraction numeric, between 0 and 1, fraction of minimal non-zero occurence of a taxon at a single site
#' @param save_to_file logical
#'
#' @return list, of shape list(taxa = tibble - taxa to consider under given occurrence fraction, taxa_stats = tibble - corresponding statistics)
#' @export
#'
#' @examples
optimize_taxa_for_arb_pollen_by_taxa_minimal_occurence_fraction <- function(pollen_data_harm, 
                                                                            harmonisation_level,
                                                                            minimal_occurence_fraction = 0.01, 
                                                                            save_to_file = FALSE) {
  taxa_considered <- filter(pollen_data_harm, taxon_pcnt > 0) %>% 
    group_by(site_id, taxon_harmonized) %>% 
    mutate(harm_pollen_occurence_at_site = n()) %>% 
    group_by(site_id) %>% 
    mutate(all_pollen_occurence_at_site = n()) %>% 
    mutate(fraction_pollen_occurence = harm_pollen_occurence_at_site / all_pollen_occurence_at_site) %>% 
    select(site_id, taxon_harmonized, fraction_pollen_occurence) %>% 
    filter(fraction_pollen_occurence >= minimal_occurence_fraction) %>% 
    ungroup() %>% 
    distinct(taxon_harmonized)
  
  stats <- calculate_taxa_to_arb_pollen_stats(pollen_data_harm = pollen_data_harm, 
                                              min_occurences = 0, #always!
                                              taxa_to_consider = taxa_considered)
  if (save_to_file) {
    library(writexl)
    write_xlsx(x = taxa_considered,
               path = paste(DIR_ARB_POLLEN,
                            paste(paste('harm_taxa', harmonisation_level, 'arb_pollen_match', minimal_occurence_fraction, 'min_occur_frac', sep = '_'), 'xlsx', sep = '.'), sep = '/'),
               col_names = TRUE, format_headers = TRUE)
    write_xlsx(x = stats,
               path = paste(DIR_ARB_POLLEN,
                            paste(paste('harm_taxa', harmonisation_level, 'arb_pollen_match_site_wise_stats', minimal_occurence_fraction, 'min_occur_frac', sep = '_'), 'xlsx', sep = '.'), sep = '/'),
               col_names = TRUE, format_headers = TRUE)
    print.xtable(xtable(taxa_considered, type = 'latex'),
                 file = paste(DIR_ARB_POLLEN,
                              paste(paste('harm_taxa', harmonisation_level, 'arb_pollen_match', minimal_occurence_fraction, 'min_occur_frac', sep = '_'), 'tex', sep = '.'), sep = '/'), 
                 include.rownames = FALSE)
    print.xtable(xtable(stats, type = 'latex'),
                 file = paste(DIR_ARB_POLLEN,
                              paste(paste('harm_taxa', harmonisation_level, 'arb_pollen_match_site_wise_stats', minimal_occurence_fraction, 'min_occur_frac', sep = '_'), 'tex', sep = '.'), sep = '/'), 
                 include.rownames = FALSE)
  }
  
  return(list(taxa = taxa_considered, taxa_stats = stats))
}



#' @rdname acer.harmonise.taxa
#' @param harm_pollen_data tibble, harmonised ACER pollen data
#' @param harmonization_level integer, either 1 or 2 (1 clips all words depite the first one, 2 does the same with separation with slashes considered)
#' @param min_occurence_fraction numeric, between 0 and 1, fraction of minimal non-zero occurence of a taxon at a single site
#'
#' @return tibble, ACER arboreal pollen data
#' @export
#' 
#' @note corresponding AP match file has to be provided
#'
#' @examples
match_top_taxa_to_harmonized_arb_pollen_by_min_occ_frac <- function(harm_pollen_data, 
                                                                    harmonization_level = 2,
                                                                    min_occurence_fraction = 0.01,
                                                                    incl_counts = FALSE){
  dir_matchfile <- paste(DIR_ARB_POLLEN,
                         paste(paste('new_matchfile', 'harm_taxa', harmonization_level, 'arb_pollen_match', min_occurence_fraction, 'min_occur_frac', sep = '_'), 'xlsx', sep = '.'), sep = '/')
  taxa <- suppressMessages(read_xlsx(path = dir_matchfile)) %>% 
    select(taxon_harmonized, is_arboreal_pollen) %>% 
    filter(is_arboreal_pollen == TRUE) %>% 
    select(-is_arboreal_pollen)
  
  pollen_data <- harm_pollen_data %>% 
    {if (!incl_counts) {select(., -taxon_count)} else {.}} %>% 
    inner_join(taxa, by = 'taxon_harmonized') %>% 
    group_by(site_id, sample_id) %>% 
    {if (!incl_counts) {summarise(., pcnt_arb_pollen = sum(taxon_pcnt))} else {summarise(., pcnt_arb_pollen = sum(taxon_pcnt), arb_pollen_count = sum(taxon_count))}} %>% 
    rownames_to_column(var = 'arb_pollen_data_id')
  
  cat(paste('matched taxa occuring with occurence fraction equal or larger than', min_occurence_fraction, '\n'))
  
  return(pollen_data)
}


calculate_taxa_to_arb_pollen_stats <- function(pollen_data_harm, 
                                               harmonization_level = 2, 
                                               min_occurences = 200, 
                                               add_summary_row = FALSE, 
                                               taxa_to_consider = NULL) {
  # addresses question: if global min_occurence limit is chosen, which percentage of samples of each site that are not zero is covered
  require(readxl)
  dir_default_matchfile <- paste(paste(DIR_ARB_POLLEN, paste('harmonized_taxa', harmonization_level, 'with_stats_200_min_occurences_with_arboreal_categorization', sep = '_'), sep = '/'), 'xlsx', sep = '.')
  taxa_considered <- {if(is.null(taxa_to_consider)){read_xlsx(path = dir_default_matchfile)}else {taxa_to_consider}} %>% 
    {if(is.null(taxa_to_consider)){filter(., taxon_occurence_harm >= min_occurences)}else {.}} %>% 
    {if(is.null(taxa_to_consider) | !('site_id' %in% names(taxa_to_consider))){select(., taxon_harmonized)}else {select(., site_id, taxon_harmonized)}}
  
  stats <- inner_join(filter(pollen_data_harm, taxon_pcnt > 0) %>% 
                        group_by(site_id) %>% 
                        summarise(nmb_harm_pollen_samples = n()), 
                      filter(inner_join(pollen_data_harm, taxa_considered, by = {if(is.null(taxa_to_consider) | !('site_id' %in% names(taxa_to_consider))){'taxon_harmonized'}else {c('taxon_harmonized', 'site_id')}}),
                             taxon_pcnt > 0) %>% 
                        group_by(site_id) %>% 
                        summarise(nmb_harm_arb_pollen_samples = n()), 
                      by = 'site_id') %>% 
    mutate(fraction_arb_harm_to_harm_pollen = nmb_harm_arb_pollen_samples / nmb_harm_pollen_samples) %>% 
    {if(add_summary_row){bind_rows(.,tibble(site_id = 0,
                                            nmb_harm_pollen_samples = (filter(pollen_data_harm, taxon_pcnt > 0) %>% count())$n, 
                                            nmb_harm_arb_pollen_samples = (filter(inner_join(pollen_data_harm, taxa_considered, by = 'taxon_harmonized'), taxon_pcnt > 0) %>% count())$n) %>% 
                                     mutate(fraction_arb_harm_to_harm_pollen = nmb_harm_arb_pollen_samples / nmb_harm_pollen_samples))}else {.}}
  
  return(stats)
}

