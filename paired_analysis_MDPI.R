# source("/home/ben/research/NOC/SRS_wave_analysis/CCI/L2/track_analysis/JMSE_coastal_extremes/paired_analysis_MDPI.R")

# This R script reads in Hs observational data from data buoys (NDBC) and satellite (CCI L2P v1.1),
# creates consistent 1-hourly pairwise data series and scatterplots of the paired time series of thei
# two data sets used in our JMSE publication (https://doi.org/10.3390/jmse8121039).
# Note that some annotation is not accurate (inherited from earlier revisions).

# Set region idx.
   reg_idx <- 6
# Set sampling radius around buoy.
   buoy_radius <- 50
# Plotting flags.
   flag_qqplot <- FALSE
   flag_scatter <- TRUE
   flag_regression <- TRUE
   if ( flag_regression ) {
      require(lmodel2)
   }
   flag_extremes <- FALSE
   flag_extremes_thresh <- TRUE

# Define buoys.
   mat_buoys <- rbind( c("44005","44007"),
                       c("41002","41110"),
                       c("41010","41113"),
                       c("42002","42035"),
                       c("46059","46013"),
                       c("46005","46041") )
   list_years <- list( 1993:2018, 2008:2018, 2007:2018, 1993:2018, 1994:2018, 1994:2018 )
# Buoys and years.
   buoy_list <- mat_buoys[reg_idx,]
   vec_years <- list_years[[reg_idx]]
# Region labels.
   lab_reg <- paste("Region #",reg_idx,sep="")

# Filter by quality flag.
   flag_qual <- TRUE
# Take 1 deg track medians, not total data.
   flag_median <- TRUE
   lab_flag_med <- "MED"
# Set TRUE for single cell of sat data containing buoy, otherwise 2x2 cells are processed.
   flag_cell <- FALSE
# If false, list files to use.
   sat_list_idx <- 1:9

# Create a time series.
# 1 hourly series
   vec_pos_time <- as.POSIXlt((1*3600)*c(0:((( length(vec_years) * 365+6)*(24 / 1))-1)),origin = paste(vec_years[1],"-01-01",sep=""))
   lab_freq <- "1hr"
## 3 hourly series
#   vec_pos_time <- as.POSIXlt((3*3600)*c(0:((( length(vec_years) * 365+2)*(24 / 3))-1)),origin = paste(vec_years[1],"-01-01",sep=""))
#   lab_freq <- "3hr"
# 6 hourly series
#   vec_pos_time <- as.POSIXlt((6*3600)*c(0:((( length(vec_years) * 365+2)*(24 / 6))-1)),origin = paste(vec_years[1],"-01-01",sep=""))
# 12 hourly series
   #BB <- as.POSIXlt((12*3600)*c(1:(( length(vec_years) * 365+2)*(24 / 12))-1),origin = "1993-01-01")
# 24 hourly series
   #BB <- as.POSIXlt((24*3600)*c(1:(( length(vec_years) * 365+2)*(24 / 24))-1),origin = "1993-01-01")

# Satellite missions.
   vec_datasets <- c("ERS-1","TOPEX","ERS-2","GFO","JASON-1","ENVISAT","JASON-2","CRYOSAT-2","SARAL","JASON-3")

# Function for mode (for satellite mission determination).
   func_mode <- function(x) {
     if (!is.null(x)) {
        ux <- unique(x)
        ux[which.max(tabulate(match(x, ux)))]
      } else {
         return(NA)
      }
   }
#--------------------------------------------------------------#
# Load buoy data.
#--------------------------------------------------------------#
   list_buoy_hs_mean <- vector(mode = "list", length = length(buoy_list))
   list_buoy_hs_raw <- vector(mode = "list", length = length(buoy_list))
   for ( b_idx in 1:length(buoy_list) ) {
# Check for data already saved and load it.
      buoy_name_temp <- buoy_list[b_idx]
      master_file_name_buoy <- paste("/home/ben/research/NOC/SRS_wave_analysis/CCI/L2/track_analysis/buoy_pair_time_series/",buoy_name_temp,"_BUOY_",lab_freq,".Robj",sep="")
      if ( file.exists(master_file_name_buoy) ) {
         attach(master_file_name_buoy)
         list_buoy_hs_mean[[b_idx]] <- master_list[[3]]
         list_buoy_hs_mean[[b_idx]][list_buoy_hs_mean[[b_idx]] < 0.05] <- NA
         list_buoy_hs_raw[[b_idx]] <- master_list[[4]]
	 detach()
# Otherwise generate.
      } else {
         buoy_name_temp <- buoy_list[b_idx]
         buoy_data_file <- list.files(path = "/home/ben/research/waves/buoy_data/NDBC_complete_records/", pattern = paste("^",buoy_name_temp,sep="") )
         mat_buoy_csv <- read.csv(paste("/home/ben/research/waves/buoy_data/NDBC_complete_records/",buoy_data_file,sep=""))
         vec_buoy_time <-  strptime(as.character(mat_buoy_csv[,1]),format="%Y-%m-%d %H:%M:%S",tz="GMT")
         vec_buoy_hs <- mat_buoy_csv$hs
# Remove NAs.
         vec_buoy_time1 <- vec_buoy_time[!is.na(vec_buoy_time)]
         vec_buoy_hs1 <- vec_buoy_hs[!is.na(vec_buoy_time)]

         list_buoy_hs <- NULL
         for (y_idx in 1:length(vec_years)) {
# Find master record entries corresponding to year of analysis.
            vec_time_year <- vec_pos_time[format(vec_pos_time,"%Y") == as.character(vec_years[y_idx])]
            list_hs_temp <- vector(mode = "list", length = length(vec_time_year))
# Get buoy data.
            if ( sum(format(vec_buoy_time1,"%Y") == as.character(vec_years[y_idx])) > 1 ) {
# Extract relevant time stamps in the buoy data.
               vec_buoy_time_year <- vec_buoy_time1[format(vec_buoy_time1,"%Y") == as.character(vec_years[y_idx])]
# Extract relevant hs measurements in the buoy data.
               vec_buoy_hs_year <- vec_buoy_hs1[format(vec_buoy_time1,"%Y") == as.character(vec_years[y_idx])]
# Loop over relevant measurements and time slots to find corresponding time slot in the master record (vec_pos_time).
               for (t_idx in 1:length(vec_buoy_time_year)) {
                  mast_time_idx <- sum(vec_time_year < vec_buoy_time_year[t_idx]+1)
                  hs_temp1 <- c(unlist(list_hs_temp[mast_time_idx]),vec_buoy_hs_year[t_idx])
                  list_hs_temp[[mast_time_idx]] <- hs_temp1
               }
            }
            list_buoy_hs <- c(list_buoy_hs,list_hs_temp)
         }
# Store raw data (not averaged).
         list_buoy_hs_raw[[b_idx]] <- list_buoy_hs
# Create vector mean Hs.
         vec_buoy_hs_mean <- numeric(length(vec_pos_time))
         for (t_idx in 1:length(vec_buoy_hs_mean)) {
            vec_buoy_hs_mean[t_idx] <- suppressWarnings( mean(unlist(list_buoy_hs[t_idx]),na.rm=T) )
         }

         list_buoy_hs_mean[[b_idx]] <- vec_buoy_hs_mean
# Save output for specific buoy.
         master_meta <- list(buoy_name_temp,lab_freq,vec_years[1],vec_years[length(vec_years)],"time","buoy mean","buoy raw")
         master_list <- list(master_meta,vec_pos_time,list_buoy_hs_mean[[b_idx]],list_buoy_hs_raw[[b_idx]])
         save(master_list,file=master_file_name_buoy)
      }
   }

# Comparison between "raw" time series and x-hourly averaged time series.
   #X11(); qqplot(unlist(list_buoy_hs),list_buoy_hs_mean[[b_idx]])

#--------------------------------------------------------------#
# Load sat data.
#--------------------------------------------------------------#
# Load single cell corresponding to the buoy, or group of cells, and restrict by distance from the buoy.
   list_sat_hs_mean <- vector(mode = "list", length = length(buoy_list))
   list_sat_hs_raw <- vector(mode = "list", length = length(buoy_list))
   list_sat_hs_mission_mode <- vector(mode = "list", length = length(buoy_list))
   list_sat_hs_med <- vector(mode = "list", length = length(buoy_list))
   list_sat_hs_med_mean <- vector(mode = "list", length = length(buoy_list))
   list_sat_hs_time_med <- vector(mode = "list", length = length(buoy_list))
   list_sat_hs_mission_med_mode <- vector(mode = "list", length = length(buoy_list))
   list_sat_hs_mission_med_mode_mode <- vector(mode = "list", length = length(buoy_list))
   #list_sat_hs_qual <- vector(mode = "list", length = length(buoy_list))
   for ( b_idx in 1:length(buoy_list) ) {
      buoy_name <- buoy_list[b_idx]
# Check for data already saved and load it.
      master_file_name_sat <- paste("/home/ben/research/NOC/SRS_wave_analysis/CCI/L2/track_analysis/buoy_pair_time_series/",buoy_name,"_SAT_",buoy_radius,"km_",lab_freq,".Robj",sep="")
      if ( file.exists(master_file_name_sat) ) {
         attach(master_file_name_sat)
         list_sat_hs_mean[[b_idx]] <- master_list[[3]]
         list_sat_hs_raw[[b_idx]] <- master_list[[4]]
         list_sat_hs_mission_mode[[b_idx]] <- master_list[[5]] 
         list_sat_hs_med_mean[[b_idx]] <- master_list[[6]]
         list_sat_hs_mission_med_mode_mode[[b_idx]] <- master_list[[7]]
	 detach()
# Otherwise generate.
      } else {
         sat_file_path <- paste("/home/ben/research/NOC/SRS_wave_analysis/CCI/L2/track_analysis/buoy_records_3x3/",buoy_name,"/",sep="")
         if ( flag_cell ) {
            sat_data_file <- list.files(path = sat_file_path, pattern = paste("^.*BUOY",sep=""))
            mat_sat_csv <- read.csv(paste(sat_file_path,sat_data_file,sep=""))
         } else {
            sat_data_file <- list.files(path = sat_file_path, pattern = buoy_name)
            mat_sat_csv <- NULL
            for ( sat_file_idx in 1:length(sat_data_file) ) {
#print(paste("File:",sat_data_file[sat_file_idx]))
               mat_sat_csv_temp <- read.csv(paste(sat_file_path,sat_data_file[sat_file_idx],sep=""))
               mat_sat_csv <- rbind(mat_sat_csv,mat_sat_csv_temp)
            }
         }

# Sample Hs within ~50 km radius of buoy.
         path_buoy_meta <- paste("/home/ben/research/waves/buoy_data/NDBC_metadata/",buoy_name,"_meta",sep="")
	 if ( buoy_name == 41113 ) {
            com_b_lat <- paste("sed -ne '4p' ",path_buoy_meta," | cut -f1 -d' '",sep="")
            com_b_lon <- paste("sed -ne '4p' ",path_buoy_meta," | cut -f3 -d' '",sep="")
         } else {
            com_b_lat <- paste("sed -ne '5p' ",path_buoy_meta," | cut -f1 -d' '",sep="")
            com_b_lon <- paste("sed -ne '5p' ",path_buoy_meta," | cut -f3 -d' '",sep="")
         }
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
         #func_buoy_hit <- function(x) { sqrt( (df_buoy_loc$lat - x[1])^2 + (df_buoy_loc$lon - x[2])^2 ) }
         vec_sat_buoy_dist <- apply(X=cbind(mat_sat_csv$lat,mat_sat_csv$lon),MAR=1,FUN=func_buoy_dist)
         if (flag_qual) {
            mat_sat_csv1 <- mat_sat_csv[vec_sat_buoy_dist < buoy_radius & mat_sat_csv[,5]==3,]
         } else {
            mat_sat_csv1 <- mat_sat_csv[vec_sat_buoy_dist < buoy_radius,]
         }
# Extract variables of interest.
         vec_sat_time <-  strptime(as.character(mat_sat_csv1[,1]),format="%Y-%m-%d %H:%M:%S",tz="GMT")
         vec_sat_hs <- mat_sat_csv1$hs
         vec_mission <- mat_sat_csv1$mission
         vec_qual <- mat_sat_csv1$qual

# Code for computing track segment medians.
         #if (flag_median) {
            sat_hs_breaks <- 1
            vec_sat_hs_Nna_idx <- which(!is.na(vec_sat_hs))
            vec_sat_time_sort <- sort.int( as.numeric(vec_sat_time[vec_sat_hs_Nna_idx]),index.return = TRUE )

            vec_sat_time2 <- vec_sat_time[vec_sat_hs_Nna_idx][vec_sat_time_sort[[2]]]
            vec_sat_hs2 <- vec_sat_hs[vec_sat_hs_Nna_idx][vec_sat_time_sort[[2]]]
            vec_sat_mission2 <- vec_mission[vec_sat_hs_Nna_idx][vec_sat_time_sort[[2]]]
            vec_sat_qual2 <- vec_qual[vec_sat_hs_Nna_idx][vec_sat_time_sort[[2]]]

            for (i in 2:length(vec_sat_hs_Nna_idx)) {
            #print(paste(" Difference:",( nc1_time_idx_cell[i] - nc1_time_idx_cell[i-1] )))
               if ( abs( as.numeric(vec_sat_time2[i]) - as.numeric(vec_sat_time2[i-1]) ) > 5 ) {
                  sat_hs_breaks <- c(sat_hs_breaks,i)
                  #print(paste(" Break before:",i))
               }
            }

            #list_time_med <- vector(mode = "list", length = length(sat_hs_breaks)-1)
            vec_time_med <- strptime(as.character(mat_sat_csv1[,1]),format="%Y-%m-%d %H:%M:%S",tz="GMT")[1]
            vec_hs_med <- numeric(length(sat_hs_breaks)-1)
            vec_mission_med <- numeric(length(sat_hs_breaks)-1)
            #vec_qual_med <- numeric(length(sat_hs_breaks)-1)
# Loop over segments for find Hs medians and time means.
            for (i in 1:(length(sat_hs_breaks)-1)) {
               track_range <- sat_hs_breaks[i]:(sat_hs_breaks[i+1]-1)
# Hs.
               vec_hs_med[i] <- median(vec_sat_hs2[track_range],na.rm=T)
# Time.
               vec_time_med[i] <- mean(vec_sat_time2[track_range],na.rm=T)
               vec_mission_med[i] <- func_mode(vec_sat_mission2[track_range])
               #list_qual_med[[i]] <- func_mode(vec_sat_qual2[track_range])
            }
            list_sat_hs_med[[b_idx]] <- vec_hs_med
            list_sat_hs_time_med[[b_idx]] <- vec_time_med
            list_sat_hs_mission_med_mode[[b_idx]] <- vec_mission_med
         #}

# Process data into time series hourly averages.
# Remove NAs from all data.
         vec_sat_time1 <- vec_sat_time[!is.na(vec_sat_time)]
         vec_hs1 <- vec_sat_hs[!is.na(vec_sat_time)]
         vec_mission1 <- vec_mission[!is.na(vec_sat_time)]
         vec_qual1 <- vec_qual[!is.na(vec_sat_time)]

         list_sat_hs <- NULL
         list_sat_hs_mission <- NULL
         list_sat_med_hs <- NULL
         list_sat_hs_med_mission <- NULL
         for (y_idx in 1:length(vec_years)) {

# Find time indices matching the year.
            vec_time_year <- vec_pos_time[format(vec_pos_time,"%Y") == as.character(vec_years[y_idx])]

# All data: Loop over relevant measurements and time slots to find master record entries corresponding to year of analysis.
            list_hs_temp <- vector(mode = "list", length = length(vec_time_year))
            list_hs_mission_temp <- vector(mode = "list", length = length(vec_time_year))
            vec_sat_time_year <- vec_sat_time1[format(vec_sat_time1,"%Y") == as.character(vec_years[y_idx])]
            if ( length(vec_sat_time_year) > 1 ) {
# Match satellite data (Hs, mission, median) for year.
               vec_hs_year <- vec_hs1[format(vec_sat_time1,"%Y") == as.character(vec_years[y_idx])]
               vec_hs_mission_year <- vec_mission1[format(vec_sat_time1,"%Y") == as.character(vec_years[y_idx])]
# All data: Loop over relevant measurements and time slots to find corresponding time slot in the master record (vec_pos_time).
#TIME_IDX <- NULL
               for (t_idx in 1:length(vec_sat_time_year)) {
                  time_idx <- sum(vec_time_year < vec_sat_time_year[t_idx])
#TIME_IDX <- c(TIME_IDX,time_idx)
#print(paste("Time idx:",time_idx))
# Hs mean.
                  hs_temp1 <- c(unlist(list_hs_temp[time_idx]),vec_hs_year[t_idx])
                  list_hs_temp[[time_idx]] <- hs_temp1
# Mission mode.
                  hs_mission_temp1 <- c(unlist(list_hs_mission_temp[time_idx]),vec_hs_mission_year[t_idx])
                  list_hs_mission_temp[[time_idx]] <- hs_mission_temp1
               }
            }
# TEST
#for ( i in 1:length(unique(TIME_IDX)) ) { print(mean(list_hs_temp[[sort(unique(TIME_IDX))[i]]])) }
            list_sat_hs <- c(list_sat_hs,list_hs_temp)
            list_sat_hs_mission <- c(list_sat_hs_mission,list_hs_mission_temp)

# Median data: Loop over relevant measurements and time slots to find corresponding time slot in the master record (vec_pos_time).
            list_hs_temp <- vector(mode = "list", length = length(vec_time_year))
            vec_sat_time_year <- list_sat_hs_time_med[[b_idx]][format(list_sat_hs_time_med[[b_idx]],"%Y") == as.character(vec_years[y_idx])]
            list_hs_mission_temp <- vector(mode = "list", length = length(vec_time_year))
            if ( length(vec_sat_time_year) > 1 ) {
               vec_hs_med_year <- list_sat_hs_med[[b_idx]][format(list_sat_hs_time_med[[b_idx]],"%Y") == as.character(vec_years[y_idx])]
               for (t_idx in 1:length(vec_sat_time_year)) {
                  time_idx <- sum(vec_time_year < vec_sat_time_year[t_idx])
# Hs median.
                  hs_med_temp1 <- c(unlist(list_hs_temp[time_idx]),vec_hs_med_year[t_idx])
                  list_hs_temp[[time_idx]] <- hs_med_temp1
# Hs median mission mode.
                  hs_mission_temp1 <- c(unlist(list_hs_mission_temp[time_idx]),vec_hs_mission_year[t_idx])
                  list_hs_mission_temp[[time_idx]] <- hs_mission_temp1
               }
            }
            list_sat_med_hs <- c(list_sat_med_hs,list_hs_temp)
            list_sat_hs_med_mission <- c(list_sat_hs_med_mission,list_hs_mission_temp)
         }
# Store raw data (not averaged).
         list_sat_hs_raw[[b_idx]] <- list_sat_hs
# Create vectors for mean Hs, mean mission, mean of median Hs and mean of median mission.
         vec_sat_hs_mean <- numeric(length(list_sat_hs))
         vec_sat_hs_mission_mode <- numeric(length(list_sat_hs))
         vec_sat_hs_med_mean <- numeric(length(list_sat_hs))
         vec_sat_hs_mission_med_mode <- numeric(length(list_sat_hs))
         for (t_idx in 1:length(vec_sat_hs_mean)) {
            vec_sat_hs_mean[t_idx] <- suppressWarnings( mean(unlist(list_sat_hs[t_idx]),na.rm=T) )
            vec_sat_hs_mission_mode[t_idx] <- suppressWarnings( func_mode(unlist(list_sat_hs_mission[t_idx])) )
# Track medians.
            vec_sat_hs_med_mean[t_idx] <- suppressWarnings( mean(unlist(list_sat_med_hs[t_idx]),na.rm=T) )
            vec_sat_hs_mission_med_mode[t_idx] <- suppressWarnings( func_mode(unlist(list_sat_hs_med_mission[t_idx])) )
         }
         list_sat_hs_mean[[b_idx]] <- vec_sat_hs_mean
         list_sat_hs_mission_mode[[b_idx]] <- vec_sat_hs_mission_mode
         list_sat_hs_med_mean[[b_idx]] <- vec_sat_hs_med_mean
         list_sat_hs_mission_med_mode_mode[[b_idx]] <- vec_sat_hs_mission_med_mode
# Save files.
         master_meta <- list(buoy_name,buoy_radius,lab_freq,vec_years[1],vec_years[length(vec_years)],"time","satellite mean","satellite raw","satellite mission mode","satellite track median","satellite mission median mode")
         master_list <- list(master_meta,vec_pos_time,list_sat_hs_mean[[b_idx]],list_sat_hs_raw[[b_idx]],list_sat_hs_mission_mode[[b_idx]],list_sat_hs_med_mean[[b_idx]],vec_sat_hs_mission_med_mode)
         save(master_list,file=master_file_name_sat)
      }
   }

# Correlation.
#   print(paste("Correlation:",cor(cbind(vec_buoy_hs_mean[[1]],vec_sat_hs_mean[[1]]),use="pairwise.complete.obs")))

#==============================================================#
# Plotting.
#==============================================================#
   if (flag_scatter) {
   #vec_cols <- rep(rainbow(length(vec_years)),each=(4*365))
   vec_cols <- rainbow(10)
   #X11(); plot(list_buoy_hs_mean[[1]],list_buoy_hs_mean[[2]],xlim=c(0,12),ylim=c(0,12),col=vec_cols,xlab="Buoy 44005",ylab="Buoy 44007"); abline(a=0,b=1)
   #X11(); plot(list_buoy_hs_mean[[1]],vec_sat_hs_mean,xlim=c(0,12),ylim=c(0,12),col=vec_cols,xlab="Buoy 44005",ylab="Sat 44005"); abline(a=0,b=1)
# General case of two buoy locations in a loop.
   mat_pl_idx <- matrix(c(1,2,1,3,2,4,3,4,2,3,1,4),2)
   #list_all_data <- c(list_buoy_hs_mean[1],list_sat_hs_mean[1],list_buoy_hs_mean[2],list_sat_hs_mean[2])
   list_all_data <- list(list_buoy_hs_mean[[1]],list_sat_hs_med_mean[[1]],list_buoy_hs_mean[[2]],list_sat_hs_med_mean[[2]])
   vec_plot_labs <- c("Buoy","Sat","Buoy","Sat")
   x_lim <- y_lim <- c(0,ceiling(max(unlist(list_all_data),na.rm=T)))
   x_lim1 <- y_lim1 <- c( 0,ceiling(max(list_all_data[[1]][ !is.na(list_all_data[[1]]) & !is.na(list_all_data[[2]]) ])) )

   vec_plot_buoys <- c(paste(buoy_list[1],"(offshore)"),paste(buoy_list[1],"(offshore)"),paste(buoy_list[2],"(inshore)"),paste(buoy_list[2],"(inshore)"))
   if (flag_qual) {
      fig_file_name <- paste("./figures/buoys_3x3/",paste(buoy_list,collapse="-"),"_",paste(vec_years[1],"-",vec_years[length(vec_years)],sep=""),"_",buoy_radius,"km_",lab_flag_med,"_",lab_freq,"_Q3_MDPI.png",sep="")
   } else {
      fig_file_name <- paste("./figures/buoys_3x3/",paste(buoy_list,collapse="-"),"_",paste(vec_years[1],"-",vec_years[length(vec_years)],sep=""),"_",buoy_radius,"km_",lab_freq,".png",sep="")
   }
   #png(fig_file_name, width = 3600, height = 2400)
   png(fig_file_name, width = 3600, height = 2400)
   par(mfrow=c(2,3),oma=c(8,8,12,9),mar=c(12,14,9,9),mgp=c(9,4,0))
   #for (pl_A_idx in 1:3) {
   #   for (pl_B_idx in (pl_A_idx+1):4) {
   for (pl_idx in 1:6) {
      pl_A_idx <- mat_pl_idx[1,pl_idx]
      pl_B_idx <- mat_pl_idx[2,pl_idx]
      if ( pl_idx == 1 | pl_idx == 4 ) {
         N_buoy <- sum(!is.na(list_all_data[[pl_A_idx]]))
         N_sat <- sum(!is.na(list_all_data[[pl_B_idx]]))
         N_pairs <- sum( !is.na(list_all_data[[pl_A_idx]]) & !is.na(list_all_data[[pl_B_idx]]) )
         plot(list_all_data[[pl_A_idx]],list_all_data[[pl_B_idx]],
              xlim=x_lim1,ylim=y_lim1,pch=19,col=vec_cols[unlist(list_sat_hs_mission_mode[[as.integer(pl_B_idx/2)]])],cex=4,cex.lab=5,cex.axis=5,cex.main=6,
              #main=paste("N_buoy = ",sum(!is.na(list_all_data[[pl_A_idx]]))," ; N_sat = ",sum(!is.na(list_all_data[[pl_B_idx]])),sep=""),
              main=bquote(N[buoy] ~ "=" ~ .(N_buoy) * ";" ~ N[sat] ~ "=" ~ .(N_sat) * ";" ~ N[pairs] ~ "=" ~ .(N_pairs)),
              xlab=paste(vec_plot_labs[pl_A_idx],vec_plot_buoys[pl_A_idx]),ylab=paste(vec_plot_labs[pl_B_idx],vec_plot_buoys[pl_B_idx]))
         abline(a=0,b=1)
# Quadrants.
         abline(v=(x_lim1[2]/2))
         abline(h=(y_lim1[2]/2))
# Regression.
         if ( flag_regression ) {
            df_reg <- data.frame(buoy=list_all_data[[pl_A_idx]],sat=list_all_data[[pl_B_idx]])
            lm_1 <- lmodel2(sat ~ buoy, data=df_reg, "interval", "interval", 99)
# RMA intercept and slope.
            abline(a=lm_1$regression.results[4,2],b=lm_1$regression.results[4,3],col="red",lwd=6)
# RMA confidence intervals.
            #abline(a=lm_1$confidence.intervals[4,2],b=lm_1$confidence.intervals[4,4])
            #abline(a=lm_1$confidence.intervals[4,3],b=lm_1$confidence.intervals[4,5])
         }
      } else if ( pl_idx == 2 | pl_idx == 3 ) {
         #main_title <- c(NA,paste(vec_years[1],"-",vec_years[length(vec_years)]),NA)
         plot(list_all_data[[pl_A_idx]],list_all_data[[pl_B_idx]],
# No colours.
         #xlim=x_lim,ylim=y_lim,cex=4,cex.lab=4,cex.axis=4,cex.main=7,main=main_title[pl_idx],
         xlim=x_lim,ylim=y_lim,cex=4,cex.lab=4,cex.axis=4,cex.main=7,
         xlab=paste(vec_plot_labs[pl_A_idx],vec_plot_buoys[pl_A_idx]),ylab=paste(vec_plot_labs[pl_B_idx],vec_plot_buoys[pl_B_idx]))
         abline(a=0,b=1)
      } else if ( pl_idx == 5 | pl_idx == 6 ) {
         plot(list_all_data[[pl_A_idx]],list_all_data[[pl_B_idx]],
              xlim=x_lim,ylim=y_lim,pch=19,col=vec_cols[unlist(list_sat_hs_mission_mode[[as.integer(pl_B_idx-2)]])],cex=4,cex.lab=4,cex.axis=4,
              xlab=paste(vec_plot_labs[pl_A_idx],vec_plot_buoys[pl_A_idx]),ylab=paste(vec_plot_labs[pl_B_idx],vec_plot_buoys[pl_B_idx]))
         abline(a=0,b=1)
      }
      if ( pl_idx == 1 ) {
         legend(x=0,y=y_lim1[2],legend=vec_datasets[1:5],cex=5,col=vec_cols[1:5],pch=19)
         legend(x=(0.65*x_lim1[2]),y=(0.35*y_lim1[2]),legend=vec_datasets[6:10],cex=5,col=vec_cols[6:10],pch=19)
      }
   }
   mtext(text=paste(lab_reg,": ",vec_years[1],"-",vec_years[length(vec_years)]), side=3, line=3, adj=0.05, cex=5, outer=TRUE)
   mtext(text="(a)", side=3, line=-5, adj=0, cex=6, outer=TRUE); mtext(text="(b)", side=3, line=-5, adj=0.34, cex=6, outer=TRUE); mtext(text="(c)", side=3, line=-5, adj=0.68, cex=6, outer=TRUE)
   mtext(text="(d)", side=3, line=-123, adj=0, cex=6, outer=TRUE); mtext(text="(e)", side=3, line=-123, adj=0.34, cex=6, outer=TRUE); mtext(text="(f)", side=3, line=-123, adj=0.68, cex=6, outer=TRUE)
   dev.off()
   system(paste("okular",fig_file_name,"&> /dev/null &"))
   }

#--------------------------------------------------------------#
# Histograms with qqplots.
#--------------------------------------------------------------#
   if (flag_qqplot) {
# Labels.
      vec_plot_labs <- c("Buoy","Sat","Buoy","Sat")
      #vec_plot_buoys <- c(paste(buoy_list[1],"(offshore)"),paste(buoy_list[1],"(offshore)"),paste(buoy_list[2],"(inshore)"),paste(buoy_list[2],"(inshore)"))
      vec_plot_buoys <- c(paste(buoy_list[1],"(offshore)"),paste(buoy_list[2],"(inshore)"))

      source("/home/ben/research/NOC/SRS_wave_analysis/analysis/functions/quantile_CI.R")
      if (flag_qual) {
         fig_qq_file_name <- paste("./figures/buoys_3x3/QQ_",buoy_list[1],"_",paste(vec_years[1],"-",vec_years[length(vec_years)],sep=""),"_Q3_MDPI.png",sep="")
      } else {
         fig_qq_file_name <- paste("./figures/buoys_3x3/QQ_",buoy_list[1],"_",paste(vec_years[1],"-",vec_years[length(vec_years)],sep=""),".png",sep="")
      }

      png(fig_qq_file_name, width = 3200, height = 1600)
      par(mgp=c(9,4,0),mfrow=c(1,2),oma=c(3,3,3,3),mar=c(11,13,9,9))

      for (b_idx in 1:2) {
#   AA <- list_buoy_hs_mean[[1]]
         AA <- unlist(list_buoy_hs_raw[[b_idx]])
         AA <- AA[!is.na(AA)]
#   BB <- list_sat_hs_mean[[1]]
         BB <- unlist(list_sat_hs_raw[[b_idx]])
         BB <- BB[!is.na(BB)]

         x_lim <- y_lim <- c(0,ceiling(max(AA,BB)))

#   CC <- sample(AA,length(BB))
#   AA_plot <- AA
#
#   hist(AA_plot,breaks=50,xlim=c(2,2.5)); abline(v=quantile(AA_plot,probs=0.99),lwd=1.5); abline(v=sort(AA_plot)[quantile.CI(length(AA_plot),q=0.99)$Interval],lty=5,lwd=1.5)
#   DD <- 0; for (i in 1:5000) { DD[i] <- quantile(sample(AA,length(BB),replace = F),probs=0.99,type=8) }
#   #X11(); hist(DD,breaks=50,xlim=c(3,7)); abline(v=mean(DD),lwd=2,col="red"); abline(v=quantile(DD,probs=c(0.025,0.975)),lty=5,lwd=1.5)
         qqplot(AA,BB,xlab=paste(vec_plot_labs[1],vec_plot_buoys[b_idx]),xlim=x_lim,ylim=y_lim,ylab=paste(vec_plot_labs[2],vec_plot_buoys[b_idx]),cex=4,cex.lab=4,cex.axis=4)
         abline(a=0,b=1,lwd=5)
         mtext(text=paste(lab_reg,": ",vec_years[1],"-",vec_years[length(vec_years)]), side=3, line=-4, adj=0.05, cex=5, outer=TRUE)
#   abline(v=mean(DD),lwd=5,col="red")
#   abline(v=quantile(DD,probs=c(0.025,0.975)),lty=5,lwd=5)
#   abline(v=quantile(BB,probs=0.99,type=8),lwd=5,col="blue")
### Bootstrap over binomial estimation method.
##   mat_FF <- matrix(NA,nrow=10000,ncol=2); for (i in 1:10000) { DD <- sample(AA,10000,replace = T); mat_FF[i,] <- sort(DD)[quantile.CI(length(DD),q=0.99)$Interval] }
##   X11(); hist(DD,breaks=50,xlim=c(3,4)); abline(v=mean(DD),lwd=2,col="red"); abline(v=c(mean(mat_FF[,1]),mean(mat_FF[,2])),lty=5,lwd=1.5)
##
      }
         dev.off()
         system(paste("okular",fig_qq_file_name,"&> /dev/null &"))
   }

#--------------------------------------------------------------#
# 10-year return level.
#--------------------------------------------------------------#
   if (flag_extremes) {
# Libraries for extremes based on EV theory.
      require(extRemes)
      require(evd)

# Chiplot.
   #X11(); chiplot(cbind(unlist(list_buoy_hs_mean[1]),unlist(list_sat_hs_mean[1])),which = 1,ylim1=c(0,1)); abline(h=0.8)

# EV fitting (block maxima). Find annual maxima.
      mat_buoy_hs_an_max <- matrix(NA,nrow=length(vec_years),ncol=2)
      mat_sat_hs_an_max <- matrix(NA,nrow=length(vec_years),ncol=2)
      for (y_idx in 1:length(vec_years)) {
         for (b_idx in 1:2) {
            mat_buoy_hs_an_max[y_idx,b_idx] <- max( list_buoy_hs_mean[[b_idx]][format(vec_pos_time,"%Y") == as.character(vec_years[y_idx])],na.rm=T )
            mat_sat_hs_an_max[y_idx,b_idx] <- max( list_sat_hs_mean[[b_idx]][format(vec_pos_time,"%Y") == as.character(vec_years[y_idx])],na.rm=T )
            #print(paste("Sat length for year:",vec_years[y_idx],":",sum(!is.na(list_sat_hs_mean[[2]][format(vec_pos_time,"%Y") == as.character(vec_years[y_idx])]))))
         }
      }

# Compare annual max between buoy and satellite.
   #X11(); plot(mat_buoy_hs_an_max[,1],mat_sat_hs_an_max[,1],xlim=c(0,10),ylim=c(0,10),main="Annual Maxima"); abline(a=0,b=1)
# Fit EVD (block maxima) to buoy data.
#   return.level(fevd(mat_buoy_hs_an_max[,1]),return.period = c(2,5,10),do.ci=T)
# Resample buoy annual maxima given frequency of sat sampling *per year*.
      array_buoy_sub_hs_an_max <- array(NA,dim=c(1000,length(vec_years),2))
      array_buoy_sub_hs_an_max_thresh <- array(NA,dim=c(1000,length(vec_years),6,2))
      vec_samp_level <- c(0.01,0.05,0.1,0.25,0.5,0.75)
      for (b_idx in 1:2) {
         vec_buoy_an_len <- numeric(length(vec_years))
         vec_sat_an_len <- numeric(length(vec_years))
         for (y_idx in 1:length(vec_years)) {
            vec_buoy_an_len[y_idx] <- sum(!is.na(list_buoy_hs_mean[[b_idx]][format(vec_pos_time,"%Y") == as.character(vec_years[y_idx])]))
            vec_sat_an_len[y_idx] <- sum(!is.na(list_sat_hs_mean[[b_idx]][format(vec_pos_time,"%Y") == as.character(vec_years[y_idx])]))
         }
         for (y_idx in 1:length(vec_years)) {
            buoy_hs_an_temp <- list_buoy_hs_mean[[b_idx]][format(vec_pos_time,"%Y") == as.character(vec_years[y_idx])]
            buoy_hs_an_temp1 <- buoy_hs_an_temp[!is.na(buoy_hs_an_temp)]
            if (! length(buoy_hs_an_temp1) == 0) {
               if ( flag_extremes_thresh ) {
                  for (i_sl in 1:6) {
                     for (i_x in 1:1000) {
                        array_buoy_sub_hs_an_max_thresh[i_x,y_idx,i_sl,b_idx] <- max(sample(buoy_hs_an_temp1,floor(vec_samp_level[i_sl]*vec_buoy_an_len[y_idx])))
                     }
                  }
               }
               for (i_x in 1:1000) {
                  array_buoy_sub_hs_an_max[i_x,y_idx,b_idx] <- max(sample(buoy_hs_an_temp1,vec_sat_an_len[y_idx]))
               }
            }
         }
      }
      mat_buoy_sub_ret_level <- matrix(NA,nrow=1000,ncol=2)
      array_buoy_sub_ret_level <- array(NA,dim=c(1000,6,2))
      for (b_idx in 1:2) {
         vec_NA <- !is.na(array_buoy_sub_hs_an_max_thresh[1,,i_sl,b_idx])
         if ( flag_extremes_thresh ) {
            for (i_sl in 1:6) {
               for (i_x in 1:1000) {
                  array_buoy_sub_ret_level[i_x,i_sl,b_idx] <- return.level(fevd(array_buoy_sub_hs_an_max_thresh[i_x,vec_NA,i_sl,b_idx]),return.period = c(10),do.ci=F)[1]
               }
            }
         }
         for (i_x in 1:1000) {
            mat_buoy_sub_ret_level[i_x,b_idx] <- return.level(fevd(array_buoy_sub_hs_an_max[i_x,vec_NA,b_idx]),return.period = c(10),do.ci=F)[1]
         }
      }
# Plot histogram and return-level estimates.
      if (flag_qual) {
         fig_evd_file_name <- paste("./figures/buoys_3x3/EVD_",paste(buoy_list,collapse="-"),"_",paste(vec_years[1],"-",vec_years[length(vec_years)],sep=""),"_",lab_freq,"_Q3_MDPI.png",sep="")
      } else {
         fig_evd_file_name <- paste("./figures/buoys_3x3/EVD_",paste(buoy_list,collapse="-"),"_",paste(vec_years[1],"-",vec_years[length(vec_years)],sep=""),"_",lab_freq,".png",sep="")
      }
      png(fig_evd_file_name, width = 1400, height = 2000)
      #X11()
      par(mfrow=c(2,1),oma=c(5,5,5,5),mar=c(7,7,7,7),mgp=c(5,2,0))
      for (b_idx in 1:2) {
         hist(mat_buoy_sub_ret_level[,b_idx],breaks=25,xlim=c(floor(min(mat_buoy_sub_ret_level)),10),
              main=paste("Buoy ",buoy_list[b_idx],": Hs(",lab_freq,") 10-year return level (m) estimates",sep=""),
              xlab="Hs (m)",
              cex=3,cex.main=3,cex.axis=3,cex.lab=3)
# 10-year from bootstrap distribution.
         abline(v=mean(mat_buoy_sub_ret_level[,b_idx]),lwd=4,col="red")
         abline(v=quantile(mat_buoy_sub_ret_level[,b_idx],probs=c(0.025,0.927)),
                lwd=4,lty=5,col="red")
# 10-year from buoy.
         buoy_ret_level <- return.level(fevd(mat_buoy_hs_an_max[,b_idx][!is.infinite(mat_buoy_hs_an_max[,b_idx])]),return.period = c(10),do.ci=T)
         abline(v=buoy_ret_level[2],lwd=4,col="blue")
         abline(v=buoy_ret_level[c(1,3)],lwd=4,col="blue",lty=5)
# 10-year from satellite.
         sat_ret_level <- return.level(fevd(mat_sat_hs_an_max[,b_idx][!is.infinite(mat_sat_hs_an_max[,b_idx])]),return.period = c(10),do.ci=T)
         abline(v=sat_ret_level[2],lwd=4,col="darkorange")
         if (b_idx == 1) {
            legend(x=floor(min(mat_buoy_sub_ret_level)),y=150,legend=c(paste("Buoy (",lab_freq,")",sep=""),paste("Satellite (",lab_freq,")",sep=""),paste("Buoy resampled",sep="")),lty=c(1,1,1),lwd=c(3,3,3),col=c("Blue","Orange","Red"),cex=2)
         }
      }
      dev.off()
      system(paste("okular",fig_evd_file_name,"&> /dev/null &"))
   }

## Plot effect of resampling by threshold.
#   b_idx <- 2
#   fig_evd_thresh_file_name <- paste("./figures/buoys_3x3/EVD_thresh_test_",buoy_list[b_idx],".png",sep="")
#   mat_test <- apply(X=array_buoy_sub_ret_level[,,b_idx],MAR=2,FUN=function(x) { quantile(x[x < 20],probs=c(0.025,0.975)) } )
#   png(fig_evd_thresh_file_name, width = 2400, height = 2400)
#   #X11()
#   par(oma=c(5,5,5,5),mar=c(7,7,7,7),mgp=c(5,2,0))
#   plot(100*vec_samp_level,apply(X=array_buoy_sub_ret_level[,,b_idx],MAR=2,FUN=function(x) { mean(x[x < 20]) } ),
#        ylim=c(0,10),axes=F,cex=3,cex.main=3,cex.lab=3,
#        main=paste(buoy_list[b_idx],": Hs(1hr) 10-yr return level",sep=""),xlab="Buoy sub-sampling level (%)",ylab="Hs(m)")
#   lines(100*vec_samp_level,mat_test[1,],lwd=3)
#   lines(100*vec_samp_level,mat_test[2,],lwd=3)
#   #abline(h=sat_ret_level[2],lwd=3,col="darkorange")
#   #abline(h=sat_ret_level[1],lwd=3,col="darkorange",lty=5)
#   #abline(h=sat_ret_level[3],lwd=3,col="darkorange",lty=5)
#   axis(1,at=100*vec_samp_level,labels=TRUE,tick=TRUE,cex.axis=3,cex.lab=3)
#   axis(2,at=1:10,cex.axis=3)
#   dev.off()
#   system(paste("okular",fig_evd_thresh_file_name,"&> /dev/null &"))
#
### Synthetic data.
##   AA <- evd::rgev(10000, loc=0, scale=1, shape=0)
##   BB <- evd::rbvevd(10000,dep=0.3)

