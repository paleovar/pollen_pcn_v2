# Script to compute Figure S32 from the Supplement
# - correlation maps of mean annual surface air temperature and mean seasonal 
#   surface air temperatures in TraCE over the Last Deglaciation

nc <- nc_open(file.path(DIR_MODEL_DATA,"trace21ka/trace.01-36.22000BP.cam2.TREFHT.22000BP_decavg_400BCE.nc"))
trace_tann <- ncvar_get(nc,"TREFHT")
trace_lon <- ncvar_get(nc,"lon")
trace_lat <- ncvar_get(nc,"lat")
trace_time <- ncvar_get(nc,"time")
trace_time_ind <- which(trace_time>=-22 & trace_time <= -6)
nc <- nc_open(file.path(DIR_MODEL_DATA,"trace21ka/trace.01-36.22000BP.cam2.TREFHT.22000BP_decavgDJF_400BCE.nc"))
trace_djf <- ncvar_get(nc,"TREFHT")
nc <- nc_open(file.path(DIR_MODEL_DATA,"trace21ka/trace.01-36.22000BP.cam2.TREFHT.22000BP_decavgJJA_400BCE.nc"))
trace_jja <- ncvar_get(nc,"TREFHT")

pdf(file = file.path(DIR_FIGURES, 'TraCE_corr_TANN_TDJF_TJJA.pdf'), width = 7, height = 12)
par(mfrow = c(3,1),
    oma = c(1,0.5,1,-0.5) + 1,
    mar = c(1,0.5,1,-0.5) + 1)

# TANN vs. TDJF
fields::image.plot(c(trace_lon[49:96]-360,trace_lon[1:48]),
                   trace_lat,
                   t(sapply(1:96,
                            function(x) sapply(1:48, function(y) cor(trace_tann[x,y,trace_time_ind],trace_djf[x,y,trace_time_ind]))))[c(49:96,1:48),],
                   zlim=c(-1,1),
                   col=col_anomalies_ipcc_temperature(),
                   xlab="Latitude",ylab="Longitude",
                   main="Correlation T_ann vs. T_djf");maps::map(add=T,interior=FALSE)

# TANN vs. TJJA
fields::image.plot(c(trace_lon[49:96]-360,trace_lon[1:48]),
                   trace_lat,
                   t(sapply(1:96,
                            function(x) sapply(1:48, function(y) cor(trace_tann[x,y,trace_time_ind],trace_jja[x,y,trace_time_ind]))))[c(49:96,1:48),],
                   zlim=c(-1,1),
                   col=col_anomalies_ipcc_temperature(),
                   xlab="Latitude",ylab="Longitude",
                   main="Correlation T_ann vs. T_jja");maps::map(add=T,interior=FALSE)

# TJJA vs. TDJF
fields::image.plot(c(trace_lon[49:96]-360,trace_lon[1:48]),
                   trace_lat,
                   t(sapply(1:96,
                            function(x) sapply(1:48, function(y) cor(trace_djf[x,y,trace_time_ind],trace_jja[x,y,trace_time_ind]))))[c(49:96,1:48),],
                   zlim=c(-1,1),
                   col=col_anomalies_ipcc_temperature(),
                   xlab="Latitude",ylab="Longitude",
                   main="Correlation T_jja vs. T_djf");maps::map(add=T,interior=FALSE)
dev.off()
gc()
