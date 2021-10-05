# Preliminary wrapper script to compute
# - principal component analysis of proxy and pseudo proxy time series.

USE_PCA_REPO_DATA <- TRUE

FILE_PCA_CACHE <- file.path(DIR_CACHE,'pre_processed_pca_data.RData')


# script flag preventing script to be executed if this has been done in the session before
# overwrite STATUS_PCA_RUN <- FALSE in global environment to force re-execution
if (!exists('STATUS_PCA_RUN')) STATUS_PCA_RUN <- FALSE

if (STATUS_PCA_RUN == FALSE) {
  message('main_pca.R has not been run in session, executing script with the defined flags')
} else {
  message('main_pca.R has been run in session, aborting executing the script. To force execution, set STATUS_PCA_RUN <- FALSE in global environment')
  invokeRestart('abort')
}


# if pre-processed simulation data should be used
if (USE_PCA_REPO_DATA) {
  if (file.exists(FILE_PCA_CACHE)) {
    message('Using pre-processed PCA data, therefore skipping computing them from scratch')
    load(FILE_PCA_CACHE)
    STATUS_PCA_RUN <- TRUE
    ## this loads main data objects into global environment:
    ## `pca_results`, `pca_results_filt`
    invokeRestart('abort')
  } else {
    message('USE_PCA_REPO_DATA == TRUE but no matching .RData file found, attempting to compute PCAs from scratch')
  }
}


# load/compute data
source('main_pcns.R') # contains ACER_ap in convenient format

# wrapper for PCA computation ----
## Note this uses the same functions as prior to computing similarity but outside all the wrappers used for that purpose
## (see the `PCNdata` class and it's generics)
compute_ap_pca <- function(data,record_window=c(6000,22000),interp_res=250,pca_limits=c(7000,21000),standardize=TRUE,trafo="none",filter=NULL,filter_sites=NULL) {
  data_as_zoos <- data$sample_dating %>% inner_join(data$arb_pollen_data, by = c('site_id', 'sample_id')) %>% dplyr::select(site_id, mixed_age,pcnt_arb_pollen) 
  if (!is.null(filter_sites)) data_as_zoos <- data_as_zoos %>% filter(site_id %in% filter_sites)
  data_as_zoos <- data_as_zoos %>% 
    filter(!is.na(mixed_age) & mixed_age >= 0) %>% group_by(site_id) %>% nest() %>% mutate(zoos = purrr::map(data, function(x) zoo(x$pcnt_arb_pollen,order.by = x$mixed_age)))
  # The windowing here can be done prior to the function
  data_as_zoos$zoos <- lapply(data_as_zoos$zoos,paleodata_windowing,record_window[1],record_window[2])
  data_as_zoos <- data_as_zoos[which(sapply(data_as_zoos$zoos,length)>=2),]
  if (trafo=="sqrt") {
    data_as_zoos$zoos <- lapply(data_as_zoos$zoos,sqrt)
  }
  if (trafo=="probit") {
    data_as_zoos$zoos <- lapply(data_as_zoos$zoos,function(x) zoo(qnorm(sapply(coredata(x),function(y) min(99.99,max(0.01,y)))/100),order.by=index(x)))
  }
  if (!is.null(filter)) {
    data_as_zoos$zoos <- lapply(data_as_zoos$zoos,function(x) gaussbandpass(x,per1=filter[1],per2=filter[2])$filt)
  }
  data_as_zoos <- data_as_zoos[which(sapply(data_as_zoos$zoos,length)>=2),]
  interp_seq <- seq(record_window[1],record_window[2],by=interp_res)
  interp_data <- t(sapply(1:length(data_as_zoos$zoos),function(i) rioja::interp.dataset(y=data.frame(y=coredata(data_as_zoos$zoos[[i]])),x=index(data_as_zoos$zoos[[i]]),xout=interp_seq,method="linear",rep.negt=FALSE)))
  ind_tmp <- which(sapply(1:dim(interp_data)[1],function(i) length(which(is.na(interp_data[i,which(pca_limits[1]<=interp_seq & pca_limits[2]>=interp_seq)]))))==0)
  interp_data_select <- interp_data[ind_tmp,which(pca_limits[1]<=interp_seq & pca_limits[2]>=interp_seq)]
  pca_result <- prcomp(t(interp_data_select),scale. = standardize)
  data_as_zoos <- data_as_zoos %>% inner_join(data$sites)
  expl_var <- apply(pca_result$x,2,var)
  expl_var <- expl_var/sum(expl_var)
  return(list(pca=pca_result,expl_var=expl_var,lon=data_as_zoos$long[ind_tmp],lat=data_as_zoos$lat[ind_tmp],time_axis=interp_seq[which(pca_limits[1]<=interp_seq & pca_limits[2]>=interp_seq)]))
}


#pca_result <- compute_ap_pca(ACERtrc$ACER_ap,trafo="probit")

#pca_result <- compute_ap_pca(ACERtrc$TRACE_pr)
# ...
#pca_result <- compute_ap_pca(ACERlc$LC_pr,pca_limits=c(7000,17000))
#pca_result <- compute_ap_pca(ACERlc$LC_ts,pca_limits=c(7000,17000))
# ...

# execute PCA ----
trafos <- c(TRACE_apsb = 'none', TRACE_ts = 'none', TRACE_pr = 'none', ACER_ap = 'probit')

# no site filtering
pca_results <- lapply(1:length(ACERtrc), function(i) compute_ap_pca(ACERtrc[[names(ACERtrc)[i]]],
                                                                    trafo = trafos[[names(ACERtrc)[i]]])) %>% 
  setNames(names(ACERtrc))

# with site filtering
sitesel <- filter_sites(sites = ACERtrc$ACER_ap$sites, dating = ACERtrc$ACER_ap$sample_dating,hres_only = TRUE, site_ids = 'all', site_lats = 'all')$site_id
pca_results_filt <- lapply(1:length(ACERtrc), function(i) compute_ap_pca(ACERtrc[[names(ACERtrc)[i]]],
                                                                         trafo = trafos[[names(ACERtrc)[i]]],
                                                                         filter_sites = sitesel)) %>% 
  setNames(names(ACERtrc))


# cache PCA if permitted ----
if (USE_PCA_REPO_DATA) {
  save(
    pca_results,pca_results_filt,
    file = FILE_PCA_CACHE
  )
}

# modifying run status of this script
STATUS_PCA_RUN <- TRUE

