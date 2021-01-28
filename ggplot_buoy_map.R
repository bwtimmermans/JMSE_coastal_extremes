# source("/home/ben/research/NOC/SRS_wave_analysis/CCI/L2/track_analysis/ggplot_buoy_map.R")
   library(ggplot2)
   library(grid)
   library(gridExtra)

## Read Irma data.
#   data_irma <- read.table("track.dat.irma",header=TRUE,skip=2)
#   path_irma <- cbind(data_irma[,c(2,3,5)],g=1,name="Irma")
## Read Maria data.
#   data_maria <- read.table("track.dat.maria",header=TRUE,skip=2)
#   path_maria <- cbind(data_maria[,c(2,3,5)],g=2,name="Maria")
## Read Ike data.
#   data_ike <- read.table("track.dat.ike",header=TRUE,skip=2)
#   path_ike <- cbind(data_ike[,c(2,3,5)],g=3,name="Ike")
## Read Katrina data.
#   data_katrina <- read.table("path_katrina.txt",header=TRUE,skip=0)
#   path_katrina <- cbind(data_katrina[,c(3,4,5)],g=4,name="Katrina")
#
#   path_ALL <- rbind(path_maria,path_irma)

# Atlantic (north).
   buoy_44005 <- data.frame(lat=43.201, lon=-69.128, name="44005", nudge=0)
   buoy_44007 <- data.frame(lat=43.525, lon=-70.141, name="44007", nudge=0)
   buoy_44008 <- data.frame(lat=40.504, lon=-69.248, name="44008", nudge=0)
   buoy_44013 <- data.frame(lat=42.346, lon=-70.651, name="44013", nudge=0)
   buoy_44011 <- data.frame(lat=41.070, lon=-66.588, name="44011", nudge=0)
   buoy_44017 <- data.frame(lat=40.693, lon=-72.049, name="44017", nudge=0)
   buoy_44018 <- data.frame(lat=42.206, lon=-70.143, name="44018", nudge=0)
   buoy_44020 <- data.frame(lat=41.493, lon=-70.279, name="44020", nudge=0)
   buoy_44008 <- data.frame(lat=40.504, lon=-69.248, name="44008", nudge=0)
   buoy_44065 <- data.frame(lat=40.369, lon=-73.703, name="44065", nudge=0)
   buoy_44066 <- data.frame(lat=39.618, lon=-72.644, name="44066", nudge=0)
   buoy_44009 <- data.frame(lat=38.457, lon=-74.702, name="44009", nudge=0)
   buoy_44014 <- data.frame(lat=36.601, lon=-74.834, name="44014", nudge=0)
   buoy_41001 <- data.frame(lat=34.724, lon=-72.317, name="41001", nudge=0)
   buoy_41025 <- data.frame(lat=35.025, lon=-75.363, name="41025", nudge=0)
   buoy_41004 <- data.frame(lat=32.501, lon=-79.099, name="41004", nudge=0)
   buoy_41008 <- data.frame(lat=31.400, lon=-80.866, name="41008", nudge=0)
   buoy_41009 <- data.frame(lat=28.508, lon=-80.185, name="41009", nudge=0)
   buoy_41010 <- data.frame(lat=28.878, lon=-78.485, name="41010", nudge=0)
   #buoy_ <- data.frame(lat=, lon=, name="", nudge=0)
# Atlantic (south).
   buoy_42060 <- data.frame(lat=16.406, lon=-63.188, name="42060", nudge=0)
   buoy_41002 <- data.frame(lat=31.760, lon=-74.84, name="41002", nudge=0)
   buoy_41040 <- data.frame(lat=14.559, lon=-53.073, name="41040", nudge=0)
   buoy_41043 <- data.frame(lat=21.132, lon=-64.856, name="41043", nudge=0)
   buoy_41044 <- data.frame(lat=21.575, lon=58.625, name="41044", nudge=0)
   buoy_41046 <- data.frame(lat=23.832, lon=-68.417, name="41046", nudge=0)
   buoy_41047 <- data.frame(lat=27.520, lon=-71.53, name="41047", nudge=0)
   buoy_41110 <- data.frame(lat=34.141, lon=-77.717, name="41110", nudge=0)
   buoy_41112 <- data.frame(lat=30.709, lon=-81.292, name="41112", nudge=0)
   buoy_41113 <- data.frame(lat=28.400, lon=-80.534, name="41113", nudge=0)
   buoy_41114 <- data.frame(lat=27.550, lon=-80.217, name="41114", nudge=0)
# Gulf.
   buoy_42001 <- data.frame(lat=25.897, lon=-89.668, name="42001", nudge=0)
   buoy_42002 <- data.frame(lat=26.091, lon=-93.758, name="42002", nudge=0)
   buoy_42003 <- data.frame(lat=26.007, lon=-85.648, name="42003", nudge=0)
   buoy_42019 <- data.frame(lat=27.907, lon=-95.352, name="42019", nudge=0)
   buoy_42036 <- data.frame(lat=28.501, lon=-84.516, name="42036", nudge=0)
   buoy_42035 <- data.frame(lat=29.232, lon=-94.413, name="42035", nudge=0)
   buoy_42039 <- data.frame(lat=28.788, lon=-86.008, name="42039", nudge=0.5)
   buoy_42040 <- data.frame(lat=29.208, lon=-88.226, name="42040", nudge=1)
   buoy_42055 <- data.frame(lat=22.120, lon=-93.96, name="42055", nudge=0)
   buoy_42059 <- data.frame(lat=15.252, lon=-67.483, name="42059", nudge=0)
   buoy_42085 <- data.frame(lat=17.869, lon=-66.532, name="42085", nudge=0)
# Pacific.
   buoy_46005 <- data.frame(lat=46.134,lon=-131.079, name="46005", nudge=0)
   buoy_46059 <- data.frame(lat=38.094,lon=-129.951, name="46059", nudge=0)
   buoy_46013 <- data.frame(lat=38.253,lon=-123.303, name="46013", nudge=0)
   buoy_46022 <- data.frame(lat=40.701,lon=-124.55, name="46022", nudge=0)
   buoy_46041 <- data.frame(lat=47.353,lon=-124.742, name="46041", nudge=0)
# Hawaii.
   buoy_51004 <- data.frame(lat=17.533,lon=-152.255, name="51004", nudge=0)
   buoy_51001 <- data.frame(lat=24.453,lon=-162.0, name="51001", nudge=0)

# Create dataframe.
   buoy_list_NA <- c("44007","44013","44011","44017","44008","44065","44066","44009","44014","41001","41025","41004","41008","41009","41010")
   buoy_list_SA <- c("42060","41002","41040","41043","41044","41046","41047")
   buoy_list_GM <- c("42001","42002","42003","42019","42036","42035","42039","42040","42055")
   #buoy_list_INSOFFS <- c("44005","44007","44008","44018","44020","41110","41112","41113",""41009,"41010","41002","42059","42085","42035","42019","42002")
   buoy_list_INSOFFS <- c("44005","44007","41002","41110","41010","41113","42002","42035","46059","46013","46005","46041")
   buoy_list_ALL <- c(buoy_list_NA,buoy_list_SA,buoy_list_GM)
   buoy_list <- buoy_list_ALL
   buoy_list <- buoy_list_INSOFFS
   dd <- NULL
   for (b_idx in 1:length(buoy_list)) {
      eval(parse(text=paste("dd <- rbind(buoy_",rev(buoy_list)[b_idx],",dd)",sep="")))
   }
   #dd <- rbind(buoy_42060,buoy_41040,buoy_41043,buoy_41046,buoy_41047,buoy_42001,buoy_42002,buoy_42003,buoy_42019,buoy_42035,buoy_42036,buoy_42039,buoy_42040,buoy_42055)
# Nudge data.
   x_nudge <- c( 1.5, 1.5,-1.5,-1.5, 1.5,-1.5,-1.50, 1.5, 2.0,1.5, 2.0, 2.0)
   y_nudge <- c(-0.75,0.5, 0.5, 0.5,-0.5, 0.5,-0.50,-0.5,-0.5,1.0,-0.5,-0.75)
   #x_nudge <- rep(2.0,length(buoy_list))
   #y_nudge <- rep(-1,length(buoy_list))
# Distance to coast.
   #vec_dist_coast <- c(67.8,6.3,332,9.7,178.4,5.2,350,31.8,538,20.9,495,31)
   vec_dist_coast <- c(68,6,332,10,178,5,350,32,538,21,495,31)
# Area rectangles and labels.
   df_rect <- data.frame(x1=c(-71,-81,-83.5,-97,-131,-132), x2=c(-65,-73,-75,-91.5,-120,-121), y1=c(40.5,31,27,24.5,36,44.5), y2=c(45,35.5,30,30,40,48), r=c("#1","#2","#3","#4","#5","#6"))

# Get time series durations.
   path_buoy_data <- "/home/ben/research/waves/buoy_data/NDBC_complete_records/"
   lab_range <- "Buoy   [Duration]      Dist to coast (km)\n"
   for (b_idx in 1:length(buoy_list)) {
   #for (b_idx in 1:5) {
      data_file <- list.files(path = path_buoy_data, pattern = paste("^",buoy_list[b_idx],"*",sep="") )
      range_start <- system(paste("sed -n 2p ",path_buoy_data,data_file," | cut -c2-5",sep=""),intern = TRUE)
      range_end <- system(paste("tail -1 ",path_buoy_data,data_file," | cut -c2-5",sep=""),intern = TRUE)
      if (b_idx == length(buoy_list)) {
         lab_range <- paste(lab_range,paste(buoy_list[b_idx]," [",range_start,"-",range_end,"]  ",vec_dist_coast[b_idx],sep=""),sep="")
      } else {
         lab_range <- paste(lab_range,paste(buoy_list[b_idx]," [",range_start,"-",range_end,"]  ",vec_dist_coast[b_idx],"\n",sep=""),sep="")
      }
   }

# Hurricane labels.
   #ka_st <- 29
   #ik_st <- 55
   #ir_st <- 51
   #ma_st <- 40
   #label_hur <- rbind(path_irma[ir_st,],path_maria[ma_st,],path_ike[ik_st,],path_katrina[ka_st,])

# Mapping.
   #myMap <- get_map(location =  c(-75,20), zoom = 4, maptype = "terrain")
   #png("./test.png",width = 4000, height = 2400)
   p1 <- ggplot(data = dd) + #scale_x_continuous(limits = c(-100, -50), expand = c(0,0)) + scale_y_continuous(limits = c(10, 50), expand = c(0,0)) +
	 ylab("Latitude\n") + xlab("\nLongitude") +
         coord_map(xlim = c(-98,-65), ylim = c(21,45), projection = "lambert", lat0 = 28, lat1 = 35) +
	 geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group), color = "#000000", fill = NA, size = 0.35) +
	 geom_polygon(data = map_data("state"), aes(x = long, y = lat, group = group), color = "#000000", fill = NA, size = 0.15) +
## Hurricane tracks.
#         geom_path(data=path_katrina[1:13,], aes(x=LON, y=LAT, group=g),colour="cornflowerblue",size=1.2,linetype=2) +
#         geom_path(data=path_katrina[13:ka_st,], aes(x=LON, y=LAT, group=g),colour="cornflowerblue",size=1.2) +
#         geom_path(data=path_ike[10:21,], aes(x=LON, y=LAT, group=g),colour="darkblue",size=1.2,linetype=2) +
#         geom_path(data=path_ike[21:ik_st,], aes(x=LON, y=LAT, group=g),colour="darkblue",size=1.2) +
#         geom_path(data=path_irma[5:22,], aes(x=LON, y=LAT, group=g),colour="darkred",size=1.2,linetype=2) +
#         geom_path(data=path_irma[22:ir_st,], aes(x=LON, y=LAT, group=g),colour="darkred",size=1.2) +
#         geom_path(data=path_maria[3:28,], aes(x=LON, y=LAT, group=g),colour="orange",size=1.2) +
#         geom_path(data=path_maria[28:ma_st,], aes(x=LON, y=LAT, group=g),colour="orange",size=1.2,linetype=2) +
#         geom_label(data=label_hur,aes(x=LON, y=LAT, label = name), label.padding = unit(0.40, "lines"), size = 7, nudge_x = 0.75, nudge_y = 0.75) +
# Buoy symbols.
         geom_point(aes(x=lon, y=lat), size=4, shape=23, fill="yellow") +
         geom_label(aes(x=lon, y=lat, label = name), label.padding = unit(0.30, "lines"), size = 4, nudge_x = x_nudge, nudge_y = y_nudge) +
         geom_label(aes(x=-100, y=37.5, label = lab_range), label.padding = unit(0.40, "lines"), size = 5, hjust=0) +
# Area boxes (in blue).
         geom_rect(data=df_rect,mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2), color="blue", alpha=0.0) +
         geom_text(data=df_rect,aes(x=x1+0.9*(x2-x1),y=y1+0.2*(y2-y1),label=r),size=6) +
# Theme stuff.
         #theme(axis.title.x = element_text(size = rel(1.8), angle = 00), axis.title.y = element_text(size = rel(1.8), angle = 90))
         theme(axis.title.x=element_blank(),
               axis.title.y=element_blank(),
               axis.text.x = element_text(size = 10),
               axis.text.y = element_text(size = 10),
	       #t,r,b,l
               plot.margin=unit(c(0,1,0,2),"cm") )

   p2 <- ggplot(data = dd) + #scale_x_continuous(limits = c(-100, -50), expand = c(0,0)) + scale_y_continuous(limits = c(10, 50), expand = c(0,0)) +
	 ylab("Latitude\n") + xlab("\nLongitude") +
         coord_map(xlim = c(-130,-116), ylim = c(20,50), projection = "lambert", lat0 = 28, lat1 = 35) +
	 geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group), color = "#000000", fill = NA, size = 0.35) +
	 geom_polygon(data = map_data("state"), aes(x = long, y = lat, group = group), color = "#000000", fill = NA, size = 0.15) +
# Buoy symbols.
         geom_point(aes(x=lon, y=lat), size=4, shape=23, fill="yellow") +
         geom_label(aes(x=lon, y=lat, label = name), label.padding = unit(0.30, "lines"), size = 4, nudge_x = x_nudge, nudge_y = y_nudge) +
# Area boxes (in blue).
         geom_rect(data=df_rect,mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2), color="blue", alpha=0.0) +
         geom_text(data=df_rect,aes(x=x1+0.9*(x2-x1),y=y1+0.2*(y2-y1),label=r),size=6) +
# Theme stuff.
         theme(axis.title.x=element_blank(),
               axis.title.y=element_blank(),
               axis.text.x = element_text(size = 10),
               axis.text.y = element_text(size = 10) )
 
   fig_file_name <- "./INSOFFS_MDPI.pdf"
   mat_lay <- cbind(c(1,1,1),matrix(rep(2,9),ncol=3))
   pdf(fig_file_name,width = 13, height = 7)
   grid.arrange(p2,p1,layout_matrix=mat_lay)
   dev.off()
   system(paste("okular",fig_file_name,"&> /dev/null &"))

