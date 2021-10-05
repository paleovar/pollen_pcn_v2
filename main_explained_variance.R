# Preliminary wrapper script to compute
# - the explained variance of the ACER arboreal pollen signal
#   for the multivariate pollen record.
FILE_EXPLVAR_CACHE <- file.path(DIR_CACHE,'pre_processed_explvar_ACERAP.RData')

if (!exists('STATUS_EXPLVAR_RUN')) STATUS_EXPLVAR_RUN <- FALSE

if (STATUS_EXPLVAR_RUN == FALSE) {
  message('main_explained_variance.R has not been run in session, executing script with the defined flags')
} else {
  message('main_explained_variance.R has been run in session, aborting executing the script. To force execution, set STATUS_EXPLVAR_RUN <- FALSE in global environment')
  invokeRestart('abort')
}


# if pre-processed simulation data should be used
if (file.exists(FILE_EXPLVAR_CACHE)) {
  message('Using pre-processed explained variance data, therefore skipping aggregating them from scratch')
  load(FILE_EXPLVAR_CACHE)
  STATUS_EXPLVAR_RUN <- TRUE
  ## this loads explained variance objects into global environment:
  ## expl_var_raw, expl_var_orb, expl_var_mil, expl_var_smil
  invokeRestart('abort')
} else {
  message('no matching .RData file found, attempting to compute explained variance from scratch')
}


#### Prepare ACER data
ACERhere <- PCNdata() %>% 
  mix_sample_dating_and_err() %>% 
  harmonized_pollen_to_db() %>% 
  arboreal_pollen_to_db()
ACERdat_orig <- ACERhere$sites %>% 
  select(site_id) %>% 
  inner_join(ACERhere$pollen_data) %>% 
  inner_join(ACERhere$sample_dating) %>% 
  select(pollen_data_id,site_id,sample_id,mixed_age,taxon_pcnt,taxon) %>% 
  window_data(data = ., type_data = 'pollen', transform = NULL, windows = WINDOWS_TRACE, labels = WINDOW_LABELS_TRACE) %>% 
  ungroup() %>% 
  select(site_id, sample_id, mixed_age, taxon, taxon_pcnt) %>% 
  mutate(taxon_pcnt = if_else(sample_id == 20795 & taxon %in% c("Unknown", "Indeterminable"),0,taxon_pcnt)) # catch NA in sample_id = 20795 "Unknown" and "Indeterminable" which are internally handled for AP conversion

ACERdat_ap <- filter_sites(ACERhere$sites, ACERhere$sample_dating, hres_only = T, 'all', 'all') %>% 
  inner_join(ACERhere$arb_pollen_data) %>% 
  inner_join(ACERhere$sample_dating) %>% 
  select(arb_pollen_data_id,site_id,sample_id,mixed_age,pcnt_arb_pollen) %>% 
  qtransform_data(data = ., type_data = 'arboreal_pollen', transform = 'identity', sd_one = F) %>% 
  window_and_detrend(data = ., windows = WINDOWS_TRACE, labels = WINDOW_LABELS_TRACE, 
                     type_data = 'arboreal_pollen', transform = 'identity', detrend = list(activate = FALSE), sd_one = F) %>% 
  unnest()

rm(ACERhere)


#### ------------- EXPL. VAR -----------------------
# -------- Raw --------------
#load("/stacywork/fredweasley/Pollen/data_supplementary/model_data/6_22_TRACE/dat/6_22_TRACE_ACER_orig.RData")
acer_list <- ACERdat_orig %>% group_by(site_id) %>% nest() %>% mutate(new = purrr::map(data, function(x) group_by(x, mixed_age) %>% spread(key=taxon, value=taxon_pcnt) %>% ungroup())) %>% dplyr::select(-data) #%>% unnest(new)
acer_list$new <- acer_list$new #[[1]]
#load("/stacywork/fredweasley/Pollen/data_supplementary/model_data/6_22_TRACE/dat/6_22_TRACE_ACER_ap_proc.RData")
tmp_ind <- which(acer_list$site_id %in% ACERdat_ap$site_id)
acer_list <- acer_list[tmp_ind,]
#Remove NA rows
for (i in 1:length(acer_list$site_id)) {
  acer_list$new[[i]] <- acer_list$new[[i]][c(which(!is.na(rowSums(acer_list$new[[i]][,-1])))),]
}

# Compute RDA
rda_list <- list()
for (i in 1:length(acer_list$site_id)) {
  data <- as.matrix(acer_list$new[[i]][,3:dim(acer_list$new[[i]])[2]])
  data_mean <- array(rep(colMeans(data),each=dim(data)[1]),dim=dim(data))
  if (dim(data)[1] >= 4 & !is.null(ACERdat_ap$site_data[[which(ACERdat_ap$site_id == acer_list$site_id[i])]]$pcnt_arb_pollen)) {
    ap <- ACERdat_ap$site_data[[which(ACERdat_ap$site_id == acer_list$site_id[i])]]$pcnt_arb_pollen
    rda_list[[i]] <- rda(X=data,Y=ap)
  } else {
    rda_list[[i]] <- NA
  }
}

# Compute explained variance
expl_var <- c()
for (i in 1:length(acer_list$site_id)) {
  data <- as.matrix(acer_list$new[[i]][,3:dim(acer_list$new[[i]])[2]])
  ap <- ACERdat_ap$site_data[[which(ACERdat_ap$site_id == acer_list$site_id[i])]]$pcnt_arb_pollen
  if (dim(data)[1] >= 4 & !is.null(ACERdat_ap$site_data[[which(ACERdat_ap$site_id == acer_list$site_id[i])]]$pcnt_arb_pollen)) {
    expl_var[i] <- rda_list[[i]]$CCA$eig[1]/sum(rda_list[[i]]$CCA$eig,rda_list[[i]]$CA$eig)
  } else {
    expl_var[i] <- NA
  }
}

expl_var_raw <- tibble(site_id = acer_list$site_id,
                       expl_var = expl_var) %>% 
  mutate(expl_var = if_else(is.na(expl_var),0,expl_var))

# ----------- Orbital only --------------------
#load("/stacywork/fredweasley/Pollen/data_supplementary/model_data/6_22_TRACE/dat/6_22_TRACE_ACER_orig.RData")
acer_list <- ACERdat_orig %>% group_by(site_id) %>% nest() %>% mutate(new = purrr::map(data, function(x) group_by(x, mixed_age) %>% spread(key=taxon, value=taxon_pcnt) %>% ungroup())) %>% dplyr::select(-data) #%>% unnest(new)
acer_list$new <- acer_list$new #[[1]]
#load("/stacywork/fredweasley/Pollen/data_supplementary/model_data/6_22_TRACE/dat/6_22_TRACE_ACER_ap_proc.RData")
tmp_ind <- which(acer_list$site_id %in% ACERdat_ap$site_id)
acer_list <- acer_list[tmp_ind,]
#Remove NA rows and filter dataset
for (i in 1:length(acer_list$site_id)) {
  acer_list$new[[i]] <- acer_list$new[[i]][c(which(!is.na(rowSums(acer_list$new[[i]][,-1])))),]
  if (dim(as.matrix(acer_list$new[[i]][,3:dim(acer_list$new[[i]])[2]]))[1] >= 2) {
    for (j in 3:(dim(acer_list$new[[i]])[2])) {
      acer_list$new[[i]][,j] <- as.numeric(gaussdetr(zoo(acer_list$new[[i]][,j],order.by = acer_list$new[[i]]$mixed_age),tsc.in = 8000)$Xsmooth)
    }
  }
}

# Compute RDA and filter AP
rda_list <- list()
for (i in 1:length(acer_list$site_id)) {
  data <- as.matrix(acer_list$new[[i]][,3:dim(acer_list$new[[i]])[2]])
  data_mean <- array(rep(colMeans(data),each=dim(data)[1]),dim=dim(data))
  if (dim(data)[1] >= 4 & !is.null(ACERdat_ap$site_data[[which(ACERdat_ap$site_id == acer_list$site_id[i])]]$pcnt_arb_pollen)) {
    ap <- as.numeric(gaussdetr(zoo(ACERdat_ap$site_data[[which(ACERdat_ap$site_id == acer_list$site_id[i])]]$pcnt_arb_pollen,order.by = acer_list$new[[i]]$mixed_age),tsc.in = 8000)$Xsmooth)
    rda_list[[i]] <- rda(X=data,Y=ap)
  } else {
    rda_list[[i]] <- NA
  }
}

# Compute explained variance
expl_var <- c()
for (i in 1:length(acer_list$site_id)) {
  data <- as.matrix(acer_list$new[[i]][,3:dim(acer_list$new[[i]])[2]])
  ap <- ACERdat_ap$site_data[[which(ACERdat_ap$site_id == acer_list$site_id[i])]]$pcnt_arb_pollen
  if (dim(data)[1] >= 4 & !is.null(ACERdat_ap$site_data[[which(ACERdat_ap$site_id == acer_list$site_id[i])]]$pcnt_arb_pollen)) {
    expl_var[i] <- rda_list[[i]]$CCA$eig[1]/sum(rda_list[[i]]$CCA$eig,rda_list[[i]]$CA$eig)
  } else {
    expl_var[i] <- NA
  }
}

expl_var_orb <- tibble(site_id = acer_list$site_id,
                       expl_var = expl_var) %>% 
  mutate(expl_var = if_else(is.na(expl_var),0,expl_var))

# ----------- Millennial only --------------------
#load("/stacywork/fredweasley/Pollen/data_supplementary/model_data/6_22_TRACE/dat/6_22_TRACE_ACER_orig.RData")
acer_list <- ACERdat_orig %>% group_by(site_id) %>% nest() %>% mutate(new = purrr::map(data, function(x) group_by(x, mixed_age) %>% spread(key=taxon, value=taxon_pcnt) %>% ungroup())) %>% dplyr::select(-data) #%>% unnest(new)
acer_list$new <- acer_list$new #[[1]]
#load("/stacywork/fredweasley/Pollen/data_supplementary/model_data/6_22_TRACE/dat/6_22_TRACE_ACER_ap_proc.RData")
tmp_ind <- which(acer_list$site_id %in% ACERdat_ap$site_id)
acer_list <- acer_list[tmp_ind,]
#Remove NA rows and filter dataset
for (i in 1:length(acer_list$site_id)) {
  acer_list$new[[i]] <- acer_list$new[[i]][c(which(!is.na(rowSums(acer_list$new[[i]][,-1])))),]
  if (dim(as.matrix(acer_list$new[[i]][,3:dim(acer_list$new[[i]])[2]]))[1] >= 2) {
    for (j in 3:(dim(acer_list$new[[i]])[2])) {
      acer_list$new[[i]][,j] <- as.numeric(gaussbandpass(zoo(acer_list$new[[i]][,j],order.by = acer_list$new[[i]]$mixed_age),per1 = 2000, per2 = 8000)$filt)
    }
  }
}

# Compute RDA and filter AP
rda_list <- list()
for (i in 1:length(acer_list$site_id)) {
  data <- as.matrix(acer_list$new[[i]][,3:dim(acer_list$new[[i]])[2]])
  data_mean <- array(rep(colMeans(data),each=dim(data)[1]),dim=dim(data))
  if (dim(data)[1] >= 4 & !is.null(ACERdat_ap$site_data[[which(ACERdat_ap$site_id == acer_list$site_id[i])]]$pcnt_arb_pollen)) {
    ap <- as.numeric(gaussbandpass(zoo(ACERdat_ap$site_data[[which(ACERdat_ap$site_id == acer_list$site_id[i])]]$pcnt_arb_pollen,order.by = acer_list$new[[i]]$mixed_age),per1 = 2000, per2 = 8000)$filt)
    rda_list[[i]] <- rda(X=data,Y=ap)
  } else {
    rda_list[[i]] <- NA
  }
}

# Compute explained variance
expl_var <- c()
for (i in 1:length(acer_list$site_id)) {
  data <- as.matrix(acer_list$new[[i]][,3:dim(acer_list$new[[i]])[2]])
  ap <- ACERdat_ap$site_data[[which(ACERdat_ap$site_id == acer_list$site_id[i])]]$pcnt_arb_pollen
  if (dim(data)[1] >= 4 & !is.null(ACERdat_ap$site_data[[which(ACERdat_ap$site_id == acer_list$site_id[i])]]$pcnt_arb_pollen)) {
    expl_var[i] <- rda_list[[i]]$CCA$eig[1]/sum(rda_list[[i]]$CCA$eig,rda_list[[i]]$CA$eig)
  } else {
    expl_var[i] <- NA
  }
}

expl_var_mil <- tibble(site_id = acer_list$site_id,
                       expl_var = expl_var) %>% 
  mutate(expl_var = if_else(is.na(expl_var),0,expl_var))

# ----------- Submillennial only --------------------
#load("/stacywork/fredweasley/Pollen/data_supplementary/model_data/6_22_TRACE/dat/6_22_TRACE_ACER_orig.RData")
acer_list <- ACERdat_orig %>% group_by(site_id) %>% nest() %>% mutate(new = purrr::map(data, function(x) group_by(x, mixed_age) %>% spread(key=taxon, value=taxon_pcnt) %>% ungroup())) %>% dplyr::select(-data) #%>% unnest(new)
acer_list$new <- acer_list$new #[[1]]
#load("/stacywork/fredweasley/Pollen/data_supplementary/model_data/6_22_TRACE/dat/6_22_TRACE_ACER_ap_proc.RData")
tmp_ind <- which(acer_list$site_id %in% ACERdat_ap$site_id)
acer_list <- acer_list[tmp_ind,]
#Remove NA rows and filter dataset
for (i in 1:length(acer_list$site_id)) {
  acer_list$new[[i]] <- acer_list$new[[i]][c(which(!is.na(rowSums(acer_list$new[[i]][,-1])))),]
  if (dim(as.matrix(acer_list$new[[i]][,3:dim(acer_list$new[[i]])[2]]))[1] >= 2) {
    for (j in 3:(dim(acer_list$new[[i]])[2])) {
      acer_list$new[[i]][,j] <- as.numeric(gaussdetr(zoo(acer_list$new[[i]][,j],order.by = acer_list$new[[i]]$mixed_age),tsc.in = 2000)$detr)
    }
  }
}

# Compute RDA and filter AP
rda_list <- list()
for (i in 1:length(acer_list$site_id)) {
  data <- as.matrix(acer_list$new[[i]][,3:dim(acer_list$new[[i]])[2]])
  data_mean <- array(rep(colMeans(data),each=dim(data)[1]),dim=dim(data))
  if (dim(data)[1] >= 4) {
    ap <- as.numeric(gaussdetr(zoo(ACERdat_ap$site_data[[which(ACERdat_ap$site_id == acer_list$site_id[i])]]$pcnt_arb_pollen,order.by = acer_list$new[[i]]$mixed_age),tsc.in = 2000)$detr)
    rda_list[[i]] <- rda(X=data,Y=ap)
  } else {
    rda_list[[i]] <- NA
  }
}

# Compute explained variance
expl_var <- c()
for (i in 1:length(acer_list$site_id)) {
  data <- as.matrix(acer_list$new[[i]][,3:dim(acer_list$new[[i]])[2]])
  ap <- ACERdat_ap$site_data[[which(ACERdat_ap$site_id == acer_list$site_id[i])]]$pcnt_arb_pollen
  if (dim(data)[1] >= 4 & !is.null(ACERdat_ap$site_data[[which(ACERdat_ap$site_id == acer_list$site_id[i])]]$pcnt_arb_pollen)) {
    expl_var[i] <- rda_list[[i]]$CCA$eig[1]/sum(rda_list[[i]]$CCA$eig,rda_list[[i]]$CA$eig)
  } else {
    expl_var[i] <- NA
  }
}

expl_var_smil <- tibble(site_id = acer_list$site_id,
                        expl_var = expl_var) %>% 
  mutate(expl_var = if_else(is.na(expl_var),0,expl_var))

# cache explained variance ----
save(
  expl_var_raw,expl_var_orb,expl_var_mil,expl_var_smil,
  file = FILE_EXPLVAR_CACHE
)

# modifying run status of this script
STATUS_EXPLVAR_RUN <- TRUE
