# source("/home/ben/research/NOC/SRS_wave_analysis/CCI/L2/track_analysis/JMSE_coastal_extremes/plot_L2_buoy_series_MDPI.R")

# This R script reads in Hs observational data from data buoys (NDBC) and satellite (CCI L2P v1.1)
# and creates time series plots of the two data sets used in our JMSE publication (https://doi.org/10.3390/jmse8121039).
# Note that some annotation is not accurate (inherited from earlier revisions).

   library(ggplot2)
   library(grid)
   library(gridExtra)
   require(ncdf4)
   source("/home/ben/research/NOC/SRS_wave_analysis/analysis/functions/quantile_CI.R")
 
# Read track data.
# Dimensions:
# 1,2: lat,lon
#   3: buoy
#   4: month
   #attach("./output/buoy_array_2011_mpi.Robj")
   attach("./output/buoy_pairs_3x3/buoy_array_2011_mpi.Robj")
   array_list_buoy_data <- list_buoy_data[[2]]
# Meta data.
   i_year <- as.integer(list_buoy_data[[1]][[2]])
   buoy_list_output <- list_buoy_data[[1]][[3]]
# 44014 41001 41002 41010 41048 42003 46005 46059 46022 51004 51001
   buoy_idx <- 8
   detach()
   vec_datasets <- c("ers-1","topex","ers-2","gfo","jason-1","envisat","jason-2","cryosat-2","saral","jason-3")
   buoy_radius <- 50

#=======================================================================#
# Data loading, processing and QC for buoys.
#=======================================================================#
# Quality control threshold flag.
# 85%
   #flag_QC_thresh <- "H"
# 75%
   flag_QC_thresh <- "L"
   #flag_QC_thresh <- AA
# Title label and threshold.
   if (flag_QC_thresh == "H") {
      #lab_QC <- "<85% rejection"
      lab_QC <- "11/12 mths, 12 days"
      QC_mth_thresh <- 10
      QC_day_thresh <- 12
      QC_hour_thresh <- 12
   } else {
      #lab_QC <- "<75% rejection"
      lab_QC <- "10/12 mths, 8 days"
      QC_mth_thresh <- 9
      QC_day_thresh <- 8
      QC_hour_thresh <- 6
   }

# Load historical data for NDBC buoys.
   #buoy_list <- "41010"
   buoy_list <- buoy_list_output[buoy_idx]
   mat_buoy_obs <- matrix(0,nrow=500000,ncol=16)
   mat_data_dims <- matrix(0,ncol=2,nrow=length(buoy_list))
# Loop over buoys.
   for (b.idx in 1:length(buoy_list)) {
      #buoy_name <- c("46066")
      buoy_name <- buoy_list[b.idx]
# File path.
      buoy_data_dir <- paste("/home/ben/research/waves/buoy_data/NDBC",buoy_name,sep="")
      data_files <- list.files(path = buoy_data_dir, pattern = paste("^",buoy_name,"h",i_year,sep="") )
      no_files <- length(data_files)

      file_path <- 0
      data_dims <- matrix(0,nrow=no_files,ncol=2)
      mat_buoy <- 0

      for (i in 1:no_files) {
         file_path[i] <- paste(buoy_data_dir,"/",data_files[i],sep='')
         temp_table <- read.table(file_path[i],skip=1)
         data_dims[i,] <- dim(temp_table)
         mat_buoy <- rbind(mat_buoy,temp_table[,c(1:16)])
      }

      mat_data_dims[b.idx,] <- dim(mat_buoy)
      #array_buoy_obs[(1:dim(mat_buoy)[1]),,b.idx] <- as.matrix(mat_buoy)
      mat_buoy_obs[(1:dim(mat_buoy)[1]),] <- as.matrix(mat_buoy)
   }

# Remove missing data, but save hourly zeros.
   mat_buoy_obs[mat_buoy_obs==99.00] <- NA
   vec_hourly <- mat_buoy_obs[,4]
   mat_buoy_obs[mat_buoy_obs==0.00] <- NA
   mat_buoy_obs[,4] <- vec_hourly

## Adjust data for the change in data formats in 2005.
#   if ( buoy_name == "46005" ) {
#      hist_data2 <- seq(which(mat_buoy_obs[,1] == 2006)[1],dim(mat_buoy_obs)[1])
#   } else {
#      hist_data2 <- seq(which(mat_buoy_obs[,1] == 2005)[1],dim(mat_buoy_obs)[1])
#   }
#   hist_data1 <- seq(1,hist_data2[1]-1)
#   hist_buoy_hs <- c(mat_buoy_obs[hist_data1,8],mat_buoy_obs[hist_data2,9])
   if ( buoy_name == "46005" ) {
      if ( mat_buoy_obs[2,1] >= 2006 ) {
         hist_buoy_hs <- mat_buoy_obs[,9]
      } else {
         hist_buoy_hs <- mat_buoy_obs[,8]
      }
   } else {
      if ( mat_buoy_obs[2,1] >= 2005 ) {
         hist_buoy_hs <- mat_buoy_obs[,9]
      } else {
         hist_buoy_hs <- mat_buoy_obs[,8]
      }
   }

# Quantiles for 
   q_plot <- quantile(hist_buoy_hs,probs=c(0.5,0.9),na.rm=T)

# Bin data by year.
   if (mat_buoy_obs[3,1] < 100) {
      start_year <- mat_buoy_obs[3,1] + 1900
   } else {
      start_year <- mat_buoy_obs[3,1]
   }
   #seq_years <- start_year:2018
   seq_years <- start_year

# Create time series indices in NetCDF format for plotting.
   list_buoy_dates <- apply( X=mat_buoy_obs[1:mat_data_dims[1,1],],MAR=1,FUN=function(x) { strptime(paste(x[1],"-",x[2],"-",x[3]," ",x[4],":00:00",sep=""),format="%Y-%m-%d %H:%M:%S",tz="GMT") } )
   df_plot_buoy <- data.frame(date=rep(as.POSIXct(NA, origin = '1981-01-01', tz='GMT'),mat_data_dims[1,1]),hs=hist_buoy_hs[1:mat_data_dims[1,1]])
   for (t_idx in 1:mat_data_dims[1,1]) { df_plot_buoy[t_idx,1] <- as.POSIXct(list_buoy_dates[[t_idx]], origin = '1981-01-01', tz='GMT') }

#-----------------------------------------------------------------------#
# Data structures for monthly raw, daily means and daily 12 hourly means.
   mat_monthly_hs <- mat_monthly_hs1 <- mat_monthly_hs2 <- matrix(list(),length(seq_years),12)
# Loop over years to get monthly raw data.
   for (yy in 1:length(seq_years)) {
      if ( seq_years[yy] <= 1998 ) {
         for (mm in 1:12) {
            #mat_monthly_hs[[yy,mm]] <- hist_buoy_hs[which(mat_buoy_obs[,1] == (seq_years[yy]-1900) & mat_buoy_obs[,2] %in% mm)]
            month_hs_idx <- which(mat_buoy_obs[,1] == (seq_years[yy]-1900) & mat_buoy_obs[,2] %in% mm)
            mat_monthly_hs[[yy,mm]] <- cbind( hist_buoy_hs[month_hs_idx],mat_buoy_obs[month_hs_idx,3],mat_buoy_obs[month_hs_idx,4] )
         }
      } else {
         for (mm in 1:12) {
            month_hs_idx <- which(mat_buoy_obs[,1] == seq_years[yy] & mat_buoy_obs[,2] %in% mm)
            mat_monthly_hs[[yy,mm]] <- cbind( hist_buoy_hs[month_hs_idx],mat_buoy_obs[month_hs_idx,3],mat_buoy_obs[month_hs_idx,4] )
         }
      }
   }

# Data structure for annual.
   list_annual_hs <- list_annual_hs1 <- list_annual_hs2 <- list()

# Quality control monthly.
   mat_month_q_flag <- matrix(TRUE,length(seq_years),12)
# Quality control annual.
   year_q_flag <- logical(length(seq_years))
   year_q_flag[1:length(seq_years)] <- TRUE

# Loop over years for quality control and daily means.
   for (yy in 1:length(seq_years)) {
      for (mm in 1:12) {
         mat_day_hs <- mat_monthly_hs[[yy,mm]]
         #if ( (sum(!is.na(mat_day_hs[,1])) < 96) | (length(unique(mat_day_hs[,2])) < 13) ) {
# Require at least 13 days, each containing at least 13 measurements (typically one per hour).
# Assume equally weighted averages.
         if ( ! sum( sapply(X=unique(mat_day_hs[,2]),FUN=function(x) { sum(!is.na(mat_day_hs[mat_day_hs[,2]==x,1])) }) > QC_hour_thresh ) > QC_day_thresh ) {
            mat_month_q_flag[yy,mm] <- FALSE
         } else {
            mat_monthly_hs1[[yy,mm]] <- sapply(X=unique(mat_day_hs[,2]),FUN=function(x) { mean(mat_day_hs[mat_day_hs[,2]==x,1],na.rm=T) })
            mat_monthly_hs2[[yy,mm]] <- apply( X=cbind(rep(unique(mat_day_hs[,2]),2),rep(c(0,1),each=length(unique(mat_day_hs[,2])))) , MAR=1, FUN=function(x) { mean(mat_day_hs[mat_day_hs[,2]==x[1] & floor(mat_day_hs[,3] / 12)==x[2],1],na.rm=T) } )
         }
      }
      #print(paste("Year:",seq_years[yy],"SUM:",sum(unlist(sapply(X=c(1,2,3,10,11,12),FUN=function(x) { mat_month_q_flag[[yy,x]] })))))
      if ( sum(unlist(sapply(X=1:12,FUN=function(x) { mat_month_q_flag[[yy,x]] }))) < QC_mth_thresh ) {
         year_q_flag[yy] <- FALSE
         list_annual_hs[[yy]] <- NA
         list_annual_hs1[[yy]] <- NA
      } else {
# Annual raw data.
         hs_temp <- unlist(sapply(X=1:12,FUN=function(x) { mat_monthly_hs[[yy,x]][,1] }))
         list_annual_hs[[yy]] <- hs_temp[!is.na(hs_temp)]
# Daily mean data.
         hs_temp1 <- unlist(sapply(X=1:12,FUN=function(x) { mat_monthly_hs1[[yy,x]] }))
         list_annual_hs1[[yy]] <- hs_temp1[!is.na(hs_temp1)]
# Daily 12 hourly data.
         hs_temp2 <- unlist(sapply(X=1:12,FUN=function(x) { mat_monthly_hs2[[yy,x]] }))
         list_annual_hs2[[yy]] <- hs_temp2[!is.na(hs_temp2)]
      }
   }
# Hs 95%ile.
   #vec_buoy_24hr <- sort(list_annual_hs[[which(seq_years == i_year)]])
   #vec_buoy_12hr <- sort(list_annual_hs2[[which(seq_years == i_year)]])
   #vec_buoy_Q95_annual <- vec_buoy_12hr[quantile.CI(n=length(vec_buoy_12hr),q=0.95)$Interval]
   #vec_buoy_Q95_annual <- c(vec_buoy_Q95_annual[1],quantile(unlist(vec_buoy_12hr),probs=0.95),vec_buoy_Q95_annual[2])

#=======================================================================#
# Extract 1Hz data from 2x2 grid cell.
#=======================================================================#
#   swh_rejection_flags:flag_meanings = "not_water sea_ice swh_validity sigma0_validity waveform_validity ssh_validity swh_rms_outlier swh_outlier" ;
#   swh_rejection_flags:flag_masks = 1LL, 2LL, 4LL, 8LL, 16LL, 32LL, 64LL, 128LL ;
#-----------------------------------------------------------------------#

# Sample Hs within ~50 km radius of buoy.
      path_buoy_meta <- paste("/home/ben/research/waves/buoy_data/NDBC_metadata/",buoy_name,"_meta",sep="")
      com_b_lat <- paste("sed -ne '5p' ",path_buoy_meta," | cut -f1 -d' '",sep="")
      com_b_lon <- paste("sed -ne '5p' ",path_buoy_meta," | cut -f3 -d' '",sep="")
      df_buoy_loc <- data.frame(lat=as.numeric(system(com_b_lat,intern = TRUE)), lon=-as.numeric(system(com_b_lon,intern = TRUE)), name=buoy_name)
# Code to approximately convert degrees to km.
# Radius of the Earth in km.
      i_radius = 6371
# Function for radians.
      func_rads <- function(x) { x * pi / 180 }
# Function for distance.
      func_buoy_dist <- function(x) {
         fl_d_lat = func_rads(x[1]) - func_rads(df_buoy_loc$lat)
         fl_d_lon = func_rads(x[2]) - func_rads(df_buoy_loc$lon)
         fl_h = sin(fl_d_lat / 2) * sin(fl_d_lat / 2) + cos( func_rads(df_buoy_loc$lat) ) * cos( func_rads(x[1]) ) * sin(fl_d_lon / 2) * sin(fl_d_lon / 2)
         2 * i_radius * asin(sqrt(fl_h))
      }
# Function for distance (basic).
      func_buoy_hit <- function(x) { sqrt( (df_buoy_loc$lat - x[1])^2 + (df_buoy_loc$lon - x[2])^2 ) }

# Start loop over buoys (not used).
   #for (b_idx in 1:dim(array_list_buoy_data)[3]) {
   #for (b_idx in 2) {
   b_idx <- buoy_idx
# Grid cell boundaries.
      vec_lat_bound <- c( min(sapply(X=1:2,FUN=function(x) { floor(array_list_buoy_data[[x,1,b_idx,2]][[1]][3]) })),
                          max(sapply(X=1:2,FUN=function(x) { ceiling(array_list_buoy_data[[x,1,b_idx,2]][[1]][3]) })) )
      vec_lon_bound <- c( min(sapply(X=1:2,FUN=function(x) { floor(array_list_buoy_data[[1,x,b_idx,2]][[1]][4]) })),
                          max(sapply(X=1:2,FUN=function(x) { ceiling(array_list_buoy_data[[1,x,b_idx,2]][[1]][4]) })) )
# Identify cell corresponding to buoy.
      #array_list_buoy_data[[1,1,b_idx,2]][[1]]$lat_mid == floor(as.numeric(as.character(array_list_buoy_data[[1,1,b_idx,2]][[1]]$buoy_lat))) + 0.5
# Loop over months.
      #for (mm_idx in 1:dim(array_list_buoy_data)[4]) {
      plot_tracks_lat <- NULL
      plot_tracks_lon <- NULL
      list_lat <- list()
      list_lon <- list()
      list_swh <- list()
      list_time <- list()
      list_swh_med <- list()
      list_swh_med_time <- list()
      list_swh_qual <- list()
      list_swh_rej <- list()
      list_swh_mission <- list()
      list_buoy_swh <- list()
      list_buoy_swh_med <- list()
      for (mm_idx in 1:dim(array_list_buoy_data)[4]) {
         lat_temp <- NULL
         lon_temp <- NULL
         swh_temp <- NULL
         swh_time_temp <- NULL
         swh_qual_temp <- NULL
         swh_rej_temp <- NULL
         swh_mission_temp <- NULL
         swh_med_temp <- NULL
         swh_med_time_temp <- NULL
         for (lat_idx in 1:3) {
            for (lon_idx in 1:3) {
# Hs measurements.
               swh_temp <- c(swh_temp,array_list_buoy_data[[lat_idx,lon_idx,b_idx,mm_idx]][[5]])
               swh_time_temp <- c(swh_time_temp,array_list_buoy_data[[lat_idx,lon_idx,b_idx,mm_idx]][[4]])
# Lat / lon.
               lat_temp <- c(lat_temp,array_list_buoy_data[[lat_idx,lon_idx,b_idx,mm_idx]][[2]])
               lon_temp <- c(lon_temp,array_list_buoy_data[[lat_idx,lon_idx,b_idx,mm_idx]][[3]])
# Data quality, rejection and mission flags.
               swh_qual_temp <- c(swh_qual_temp,array_list_buoy_data[[lat_idx,lon_idx,b_idx,mm_idx]][[7]])
               swh_rej_temp <- c(swh_rej_temp,array_list_buoy_data[[lat_idx,lon_idx,b_idx,mm_idx]][[8]])
               swh_mission_temp <- c(swh_mission_temp,array_list_buoy_data[[lat_idx,lon_idx,b_idx,mm_idx]][[10]])
# Hs median measurements.
# Now medians calculated locally due to error flags.
               #swh_med_temp <- c(swh_med_temp,array_list_buoy_data[[lat_idx,lon_idx,b_idx,mm_idx]][[7]])
# Results for buoy cell.
               #if ( floor(array_list_buoy_data[[lat_idx,lon_idx,b_idx,mm_idx]][[1]][3]) == floor(array_list_buoy_data[[lat_idx,lon_idx,b_idx,mm_idx]][[1]][5]) &
               #     floor(array_list_buoy_data[[lat_idx,lon_idx,b_idx,mm_idx]][[1]][4]) == floor(array_list_buoy_data[[lat_idx,lon_idx,b_idx,mm_idx]][[1]][6]) ) {
               #   list_buoy_swh[[mm_idx]] <- swh_temp
               #   list_buoy_swh_med[[mm_idx]] <- swh_med_temp
               #}
# Hs median measurements.
## Satellite tracks and medians.
#               if (!is.null(array_list_buoy_data[[lat_idx,lon_idx,b_idx,mm_idx]][[6]])) {
#                  list_breaks <- array_list_buoy_data[[lat_idx,lon_idx,b_idx,mm_idx]][[6]]
#                  plot_tracks_lat <- c( plot_tracks_lat, unlist( lapply(X=list_breaks,FUN=function(x) { c(array_list_buoy_data[[lat_idx,lon_idx,b_idx,mm_idx]][[2]][x],NA) } ) ) )
#                  plot_tracks_lon <- c( plot_tracks_lon, unlist( lapply(X=list_breaks,FUN=function(x) { c(array_list_buoy_data[[lat_idx,lon_idx,b_idx,mm_idx]][[3]][x],NA) } ) ) )
## Medians based on Q=3 (Good) data.
#                  swh_med_temp <- c(swh_med_temp,unlist(
#                                          lapply(X=list_breaks,FUN=function(x) {
#                                             vec_dist <- apply(X=cbind(unlist(array_list_buoy_data[[lat_idx,lon_idx,b_idx,mm_idx]][[2]][x]),unlist(array_list_buoy_data[[lat_idx,lon_idx,b_idx,mm_idx]][[3]][x])),MAR=1,FUN=func_buoy_dist)
#                                             flag_50km_Q3 <- array_list_buoy_data[[lat_idx,lon_idx,b_idx,mm_idx]][[7]][x] == 3 & vec_dist < buoy_radius
#                                             median(array_list_buoy_data[[lat_idx,lon_idx,b_idx,mm_idx]][[5]][x][flag_50km_Q3]) } ) ) )
#                  swh_med_time_temp <- c(swh_med_time_temp,unlist( lapply(X=list_breaks,FUN=function(x) { mean(array_list_buoy_data[[lat_idx,lon_idx,b_idx,mm_idx]][[4]][x][array_list_buoy_data[[lat_idx,lon_idx,b_idx,mm_idx]][[7]][x]==3]) } ) ) )
               #} else {
               #   plot_tracks_lat <- c( plot_tracks_lat, c(array_list_buoy_data[[lat_idx,lon_idx,b_idx,mm_idx]][[2]],NA) )
               #   plot_tracks_lon <- c( plot_tracks_lon, c(array_list_buoy_data[[lat_idx,lon_idx,b_idx,mm_idx]][[3]],NA) )
               #}
            }
         }
         list_lat[[mm_idx]] <- lat_temp
         list_lon[[mm_idx]] <- lon_temp
         list_swh[[mm_idx]] <- swh_temp
         list_time[[mm_idx]] <- swh_time_temp
         list_swh_med[[mm_idx]] <- swh_med_temp
         list_swh_med_time[[mm_idx]] <- swh_med_time_temp
         list_swh_qual[[mm_idx]] <- swh_qual_temp
         list_swh_rej[[mm_idx]] <- swh_rej_temp
         list_swh_mission[[mm_idx]] <- swh_mission_temp
      }

# Sample Hs within ~50 km radius of buoy.
      path_buoy_meta <- paste("/home/ben/research/waves/buoy_data/NDBC_metadata/",buoy_name,"_meta",sep="")
      com_b_lat <- paste("sed -ne '5p' ",path_buoy_meta," | cut -f1 -d' '",sep="")
      com_b_lon <- paste("sed -ne '5p' ",path_buoy_meta," | cut -f3 -d' '",sep="")
      df_buoy_loc <- data.frame(lat=as.numeric(system(com_b_lat,intern = TRUE)), lon=-as.numeric(system(com_b_lon,intern = TRUE)), name=buoy_name)
# Code to approximately convert degrees to km.
# Radius of the Earth in km.
      i_radius = 6371
# Function for radians.
      func_rads <- function(x) { x * pi / 180 }
# Function for distance.
      func_buoy_dist <- function(x) {
         fl_d_lat = func_rads(x[1]) - func_rads(df_buoy_loc$lat)
         fl_d_lon = func_rads(x[2]) - func_rads(df_buoy_loc$lon)
         fl_h = sin(fl_d_lat / 2) * sin(fl_d_lat / 2) + cos( func_rads(df_buoy_loc$lat) ) * cos( func_rads(x[1]) ) * sin(fl_d_lon / 2) * sin(fl_d_lon / 2)
         2 * i_radius * asin(sqrt(fl_h))
      }
# Function for distance (basic).
      func_buoy_hit <- function(x) { sqrt( (df_buoy_loc$lat - x[1])^2 + (df_buoy_loc$lon - x[2])^2 ) }
      #vec_sat_buoy_dist <- apply(X=cbind(mat_sat_csv$lat,mat_sat_csv$lon),MAR=1,FUN=func_buoy_hit)
      vec_sat_buoy_dist <- apply(X=cbind(unlist(list_lat),unlist(list_lon)),MAR=1,FUN=func_buoy_dist)
      vec_flag_50km <- vec_sat_buoy_dist < buoy_radius
      vec_flag_50km_Q3 <- vec_sat_buoy_dist < buoy_radius & unlist(list_swh_qual) == 3

# Breaks for fixed radius.
      swh_breaks <- NULL
      vec_all_time_sort <- sort( unlist(list_time)[vec_flag_50km_Q3],index.return = TRUE )
      vec_all_time <- vec_all_time_sort[[1]]
      vec_all_time_idx <- vec_all_time_sort[[2]]
      for (i in 2:length(vec_all_time)) {
      #print(paste(" Difference:",( nc1_time_idx_cell[i] - nc1_time_idx_cell[i-1] )))
         if ( abs( vec_all_time[i] - vec_all_time[i-1] ) > 5 ) {
            swh_breaks <- c(swh_breaks,i)
            #print(paste(" Break before:",i))
         }
      }
      if ( is.null(swh_breaks) ) {
         #nc1_breaks <- c(1,length(nc1_time_idx_cell))
         #mat_nc1_breaks <- matrix(nc1_breaks,ncol=2)
         list_breaks <- 1:length(vec_all_time)
      } else {
# Matrix of 'block' indices for contiguous readings, and a list of indices.
         mat_breaks <- cbind(c(1,swh_breaks),c(swh_breaks-1,length(vec_all_time)))
         list_breaks <- list()
         for ( i in 1:dim(mat_breaks)[1] ) { list_breaks[[i]] <- mat_breaks[i,1]:mat_breaks[i,2] }
      }
# Hs track medians.
      vec_all_swh <- unlist(list_swh)[vec_flag_50km_Q3][vec_all_time_idx]
      vec_swh_med <- unlist( lapply(X=list_breaks,FUN=function(x) { median(vec_all_swh[x]) } ) )
# Time indices for medians.
      vec_swh_med_time <- unlist( lapply(X=list_breaks,FUN=function(x) { mean(vec_all_time[x]) } ) )
# Data frame for median sat.
      df_plot_sat_med <- data.frame(date=as.POSIXct(vec_swh_med_time, origin = '1981-01-01', tz='GMT'),hs=vec_swh_med)
# Data frame for raw sat.
      df_plot_sat <- data.frame(date=as.POSIXct(unlist(list_time), origin = '1981-01-01', tz='GMT'),hs=unlist(list_swh),rej=unlist(list_swh_rej))
# Replace zeros with NAs (for plotting).
      df_plot_sat[df_plot_sat[,3] == 0,3] <- NA
      #for (m_idx in unique(unlist(list_swh_mission))) {
      #   miss_idx <- unlist(list_swh_mission) == m_idx
      #   df_plot_mission <- data.frame(date=as.POSIXct(unlist(list_time)[miss_idx], origin = '1981-01-01', tz='GMT'),hs=unlist(list_swh)[miss_idx],rej=unlist(list_swh_rej)[miss_idx])
      #   print(paste("Dim:",dim(df_plot_mission)))
      #}

# Plotting.
      vec_cols <- c(NA,NA,NA,NA,"red","blue","green","blueviolet")
      vec_cols <- rainbow(8)
# Pastel colours.
      vec_cols <- rev(c("#1B9E77","#7570B3","#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666"))
      vec_cols <- rev(c("#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf","#e41a1c","#377eb8"))
      fig_file_name <- paste("./figures/time_series/series_",i_year,"_",buoy_name,"_",buoy_radius,"km_MDPI.png",sep="")

# Expand window.
      i_lim1 <- 1313798400
      i_lim2 <- i_lim1 + (88 * 24 * 3600)
# Open file.
      png(fig_file_name, width = 2600, height = 4000)
      par(mfrow=c(3,1),oma=c(4,5,4,4),mar=c(11,13,11,11),mgp=c(9,4,0))

# Panel A.
      plot(df_plot_sat[vec_flag_50km,c(1,2)],ylim=c(0,11),pch=19,cex=1.8,cex.lab=5,cex.axis=5,cex.main=5,main=paste("(a) Time series of Hs at Buoy: ",buoy_name," [",i_year,"]",sep=""),xlab="Date",ylab="Hs (m)")
      points(df_plot_buoy,pch=19,cex=0.8,col="red")
# Lines showing zoom.
      abline(v=i_lim1,lwd=6)
      abline(v=i_lim2,lwd=6)
      legend("topright",legend=c("Buoy","CCI L2 (1 Hz)","Rejection Flag #","[16] = waveform_validity","[64] = swh_rms_outlier","[128] = swh_outlier"),text.col=c("red","black","blue","blue","blue","blue"),pch=c(19,19,5,NA,NA,NA),col=c("red","black","blue",NA,NA,NA),cex=6)
      par(new=TRUE)
# Rejection flags.
      plot(df_plot_sat[vec_flag_50km,c(1,3)],axes=FALSE,xlab="",ylab="",ylim=c(0,150),pch=5,cex=6,col="blue")
      axis(4,at=c(16,64,80,128),ylim=c(0,150),col="blue",col.axis="blue",las=1,cex.axis=5,ylab="Rejection flag")
      abline(h=c(16,64,80,128),col="blue")

# Range for hurricane peaks at 41010 (2011), by mission.
      x_lim <- c(i_lim1,i_lim2)

# Panel B.
      plot(df_plot_sat_med[vec_flag_50km,c(1,2)],xlim=x_lim,ylim=c(0,11),cex.lab=5,cex.axis=5,cex.main=5,main="(b) Time series of Hs by satellite mission",xlab="Date",ylab="Hs (m)")
      for (m_idx in 1:8) {
         miss_idx <- unlist(list_swh_mission)[vec_flag_50km] == m_idx
         df_plot_mission <- data.frame(date=as.POSIXct(unlist(list_time)[vec_flag_50km][miss_idx], origin = '1981-01-01', tz='GMT'),hs=unlist(list_swh)[vec_flag_50km][miss_idx],rej=unlist(list_swh_rej)[vec_flag_50km][miss_idx])
         points(df_plot_mission[,c(1,2)],pch=19,cex=5.0,col=vec_cols[m_idx])
      }
      points(df_plot_buoy,pch=19,cex=0.8,col="red")
      legend("topright",legend=c("Buoy","Jason-1 (CCI 1 Hz)","Envisat (CCI 1 Hz)","Jason-2 (CCI 1 Hz)","Cryosat (CCI 1 Hz)"),pch=c(19,19,19,19,19),col=c("red",vec_cols[c(5,6,7,8)]),cex=6)

# Panel C.
# Range for hurricane peaks at 41010 (2011).
      plot(df_plot_sat_med[,c(1,2)],xlim=x_lim,ylim=c(0,11),pch='+',cex=7.0,cex.lab=5,cex.axis=5,cex.main=5,main="(c) Time series of satellite track median Hs",xlab="Date",ylab="Hs (m)")
      points(df_plot_buoy,pch=19,cex=0.8,col="red")
      legend("topright",legend=c("Buoy","Track section median"),pch=list(19,3),col=c("red","black"),cex=6)

      dev.off()
      system(paste("okular",fig_file_name,"&> /dev/null &"))

