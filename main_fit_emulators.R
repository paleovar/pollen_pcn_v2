# Wrapper script to fit emulators for the TSF response of TraCE and HadCM3
# to forcing by temperature, precipitation, and CO2, and all three combined.

USE_PREFITTED_EMUL <- TRUE
FIT_DATA_FRACTION <- 0.5 # use every second time slice for fitting (50%) and other half for validation

FILE_EMUL_FIT_CACHE <- file.path(DIR_CACHE,paste0('pre_fitted_emul_', FIT_DATA_FRACTION,'.RData'))


# script flag preventing script to be executed if this has been done in the session before
# overwrite STATUS_EMUL_FIT_RUN <- FALSE in global environment to force re-execution
if (!exists('STATUS_EMUL_FIT_RUN')) STATUS_EMUL_FIT_RUN <- FALSE

if (STATUS_EMUL_FIT_RUN == FALSE) {
  message('main_emulation_data.R has not been run in session, executing script with the defined flags')
} else {
  message('main_emulation_data.R has been run in session, aborting executing the script. To force execution, set STATUS_EMUL_FIT_RUN <- FALSE in global environment')
  invokeRestart('abort')
}


# if pre-fitted emulators should be used
if (USE_PREFITTED_EMUL) {
  if (file.exists(FILE_EMUL_FIT_CACHE)) {
    message('Using pre-fitted emulators, therefore skipping fitting them from scratch')
    load(FILE_EMUL_FIT_CACHE)
    STATUS_EMUL_FIT_RUN <- TRUE
    ## this loads main data objects into global environment:
    ## `trace_gam_all_te` - Emulator fit TraCE21ka simulation
    ## `hadcm3_gam_all_te` - Emulator fit HadCM3 BRIDGE simulation
    ## `trace_co2_ts`, `hadcm3_co2_ts` - interpolated CO2 forcing
    invokeRestart('abort')
  } else {
    message('USE_PREFITTED_EMUL == TRUE but no matching .RData file found, attempting to fit emulators from scratch. This will take a while.')
  }
}


# function for evaluating the emulator fit ----
evaluate_tsf_gam <- function(predict_data,data,scatterplot=TRUE) {
  cat("RMSE: ",sqrt(mean((predict_data-data$y)^2)),"\n",sep="")
  cat("MAE: ",mean(abs(predict_data-data$y)),"\n",sep="")
  cat("Expl.var: ",1-mean((predict_data-data$y)^2)/var(data$y),"\n",sep="")
  if (scatterplot==TRUE) {
    plot(data$y,predict_data,pch=20,cex=0.1)
  }
}


# Load datasets ----
## TraCE
forcings <- load_forcing_dataset()
nc <- nc_open(file.path(DIR_MODEL_DATA,'trace21ka/trace.01-36.22000BP.cam2.TREFHT.22000BP_decavg_400BCE.nc'))
trace_lon <- ncvar_get(nc,"lon")
trace_lat <- ncvar_get(nc,"lat")
trace_time <- ncvar_get(nc,"time")[1:1600]*1000
trace_ts <- ncvar_get(nc,"TREFHT")[,,1:1600]
nc <- nc_open(file.path(DIR_MODEL_DATA,'trace21ka/trace.01-36.22000BP.cam2.PRECT.22000BP_decavg_400BCE.nc'))
trace_pr <- ncvar_get(nc,"PRECT")[,,1:1600]*60*60*24*30*12*1000
nc <- nc_open(file.path(DIR_MODEL_DATA,'trace21ka/processed/trace.01-36.22000BP.PFTFRAC_decavg_400BCE.nc'))
trace_pfts <- ncvar_get(nc,"PFTFRAC")
trace_pfts[which(is.na(trace_pfts))] <- 0
trace_tsf <- apply(trace_pfts[,,,2:12],c(1,2,3),sum)/apply(trace_pfts[,,,2:15],c(1,2,3),sum)
trace_tsf <- aperm(trace_tsf,c(2,3,1))[,,1:1600]
rm(trace_pfts)
trace_co2_ts <- paleodata_interpolation(rev_time_axis(forcings$co2),interpolation_type="spline",interpolation_dates=trace_time)
trace_co2 <- array(rep(trace_co2_ts,each=length(trace_lon)*length(trace_lat)),dim=dim(trace_ts))
gc()


# HadCM3 ----
# 1) Broadleaf tree [Tropical], 2) Needleleaf tree, 3) C3 grass, 4) C4 grass, 5) Shrub, 6) Urban, 7) Inland water, 8) Bare soil, 9) Ice
nc <- nc_open('/stacywork/fredweasley/data/BRIDGE_work/21ka.temp_mm_1_5m.annual.nc')
hadcm3_lon <- ncvar_get(nc,"longitude")
hadcm3_lat <- ncvar_get(nc,"latitude")
hadcm3_time <- -22000:-6001
hadcm3_ts <- ncvar_get(nc,"temp_mm_1_5m")[,,21500:5501]
decs <- seq(1,dim(hadcm3_ts)[3],by=10)
hadcm3_time <- hadcm3_time[decs]
hadcm3_ts <- sapply(decs, function(t) apply(hadcm3_ts[,,t:t+9], c(1,2), mean)) %>% array(dim = c(dim(hadcm3_ts)[c(1,2)],length(decs)))
nc <- nc_open('/stacywork/fredweasley/data/BRIDGE_work/21ka.precip_mm_srf.annual.nc')
hadcm3_pr <- ncvar_get(nc,"precip_mm_srf")[,,21500:5501]*60*60*24*30*12
hadcm3_pr <- sapply(decs, function(t) apply(hadcm3_pr[,,t:t+9], c(1,2), mean)) %>% array(dim = c(dim(hadcm3_pr)[c(1,2)],length(decs)))
nc <- nc_open('/stacywork/fredweasley/data/BRIDGE_work/21ka.PFTsumnorm.annual.nc')
hadcm3_tsf <- ncvar_get(nc,"fracPFTs_mm_srf")[,,21500:5501]
hadcm3_tsf <- sapply(decs, function(t) apply(hadcm3_tsf[,,t:t+9], c(1,2), mean)) %>% array(dim = c(dim(hadcm3_tsf)[c(1,2)],length(decs)))
hadcm3_co2_ts <- paleodata_interpolation(rev_time_axis(forcings$co2),interpolation_type="spline",interpolation_dates=seq(-22000,-6000,by=1000))
hadcm3_co2 <- array(rep(rep(hadcm3_co2_ts,times=c(50,rep(100,times=15),50)),each=length(hadcm3_lon)*length(hadcm3_lat)),dim=dim(hadcm3_ts))
gc()

# Train and evaluate HadCM3 emulation ----
cal_timesteps <- seq(1,length(hadcm3_time),by=2)
val_timesteps <- (1:length(hadcm3_time))[!(1:length(hadcm3_time)) %in% cal_timesteps]
hadcm3_data <- data.frame(y=hadcm3_tsf[,,cal_timesteps][which(!is.na(hadcm3_tsf[,,cal_timesteps]))],x1=hadcm3_ts[,,cal_timesteps][which(!is.na(hadcm3_tsf[,,cal_timesteps]))],x2=hadcm3_pr[,,cal_timesteps][which(!is.na(hadcm3_tsf[,,cal_timesteps]))],x3=hadcm3_co2[,,cal_timesteps][which(!is.na(hadcm3_tsf[,,cal_timesteps]))])
hadcm3_data_val <- data.frame(y=hadcm3_tsf[,,val_timesteps][which(!is.na(hadcm3_tsf[,,val_timesteps]))],x1=hadcm3_ts[,,val_timesteps][which(!is.na(hadcm3_tsf[,,val_timesteps]))],x2=hadcm3_pr[,,val_timesteps][which(!is.na(hadcm3_tsf[,,val_timesteps]))],x3=hadcm3_co2[,,val_timesteps][which(!is.na(hadcm3_tsf[,,val_timesteps]))])
hadcm3_gam_all_te <- mgcv::gam(y ~ te(x1,x2,x3), family=binomial,data=hadcm3_data)
predict_hadcm3_gam_all_te <- predict(hadcm3_gam_all_te,type="response",newdata=hadcm3_data_val)
evaluate_tsf_gam(predict_hadcm3_gam_all_te,hadcm3_data_val,scatterplot = FALSE)
hadcm3_predict <- hadcm3_tsf[,,val_timesteps]
hadcm3_predict[which(!is.na(hadcm3_predict))] <- predict_hadcm3_gam_all_te
cor_field <- array(NA,dim=c(length(hadcm3_lon),length(hadcm3_lat)))
for (i in 1:length(hadcm3_lon)) {
  for (j in 1:length(hadcm3_lat)) {
    cor_field[i,j] <- cor(hadcm3_tsf[i,j,val_timesteps],hadcm3_predict[i,j,],use="pairwise.complete.obs")
  }
}
image.plot(c(hadcm3_lon[49:96]-360,hadcm3_lon[1:48]),rev(hadcm3_lat),cor_field[c(49:96,1:48),73:1],col=col_anomalies_ipcc_temperature(),zlim=c(-1,1),main="HadCM3 , correlation TSF vs. emulated TSF",xlab="",ylab="");maps::map(add=T)
gc()

# Train and evaluate TraCE emulation ----
cal_timesteps <- seq(1,length(trace_time),by=2)
val_timesteps <- (1:length(trace_time))[!(1:length(trace_time)) %in% cal_timesteps]
trace_data <- data.frame(y=trace_tsf[,,cal_timesteps][which(!is.na(trace_tsf[,,cal_timesteps]))],x1=trace_ts[,,cal_timesteps][which(!is.na(trace_tsf[,,cal_timesteps]))],x2=trace_pr[,,cal_timesteps][which(!is.na(trace_tsf[,,cal_timesteps]))],x3=trace_co2[,,cal_timesteps][which(!is.na(trace_tsf[,,cal_timesteps]))])
trace_data_val <- data.frame(y=trace_tsf[,,val_timesteps][which(!is.na(trace_tsf[,,val_timesteps]))],x1=trace_ts[,,val_timesteps][which(!is.na(trace_tsf[,,val_timesteps]))],x2=trace_pr[,,val_timesteps][which(!is.na(trace_tsf[,,val_timesteps]))],x3=trace_co2[,,val_timesteps][which(!is.na(trace_tsf[,,val_timesteps]))])
trace_gam_all_te <- mgcv::gam(y ~ te(x1,x2,x3), family=binomial,data=trace_data)
predict_trace_gam_all_te <- predict(trace_gam_all_te,type="response",newdata=trace_data_val)
evaluate_tsf_gam(predict_trace_gam_all_te,trace_data_val,scatterplot = FALSE)
trace_predict <- trace_tsf[,,val_timesteps]
trace_predict[which(!is.na(trace_predict))] <- predict_trace_gam_all_te
cor_field <- array(NA,dim=c(length(trace_lon),length(trace_lat)))
for (i in 1:length(trace_lon)) {
  for (j in 1:length(trace_lat)) {
    cor_field[i,j] <- cor(trace_tsf[i,j,val_timesteps],trace_predict[i,j,],use="pairwise.complete.obs")
  }
}
fields::image.plot(c(trace_lon[49:96]-360,trace_lon[1:48]),trace_lat,cor_field[c(49:96,1:48),],col=col_anomalies_ipcc_temperature(),zlim=c(-1,1),main="trace , correlation TSF vs. emulated TSF",xlab="",ylab="");maps::map(add=T)
gc()


# cache fitted emulators if permitted ----
if (USE_PREFITTED_EMUL) {
  save(
    trace_gam_all_te,
    hadcm3_gam_all_te,
    trace_co2_ts,
    hadcm3_co2_ts,
    file = FILE_EMUL_FIT_CACHE
  )
}

# modifying run status of this script
STATUS_EMUL_FIT_RUN <- TRUE
