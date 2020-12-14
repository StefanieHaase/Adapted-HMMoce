##Geolocation model##
#based on HMMoce package by Camrin Braun 
#https://github.com/camrinbraun/HMMoce
#adapted for Eastern Baltic cod in the southern Baltic Sea by Stefanie Haase since 21.03.2019

###ARTIFICIALLY SIMULATED TRACK (GER)

rm(list=ls())
setwd("C:/Users/haase_s/Desktop/HMMoce")
#library(ggplot2)
#library(readr)
#library(dplyr)
#library(lubridate)
#library(StreamMetabolism)
#library(xts)
#library(reshape)
#library(scales)
#library(gridExtra)
#library(zoo)
#library(TTR)
#library(forecast)
#library(fUnitRoots)
#library(MetaCycle)
library(HMMoce)
#library(fields)
library(ncdf4) #newly added
#library(RNetCDF)
library(RColorBrewer)
library(SpaDES)
library(geosphere)
library(dplyr)
library(raster)
#------------
# LOAD THE TAG DATA
#------------

# SET INITIAL LOCATIONS (TAG AND POP-UP)
#release date, lat, long, recapture date, lat, long
iniloc <- data.frame(matrix(c(1, 2, 2016, 54.9730, 13.6750, 1, 2, 2017, 54.13151, 14.09495), nrow = 2, ncol = 5, byrow = TRUE)) #2222_noise



names(iniloc) <- list('day','month','year','lat','lon')
tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y', tz='UTC')
pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y', tz='UTC')

iniloc2 <- data.frame(date=c(tag,pop), lat=iniloc$lat, lon=iniloc$lon)
iniloc2$date <- as.POSIXct(iniloc2$date, format = '%d/%m/%Y', tz='UTC')
# VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
dateVec <- as.Date(seq(tag, pop, by = 'day'))
dates <- as.data.frame(dateVec); colnames(dates) <- c("date")
ptt<-2222
# READ IN DATA AS OUTPUT FROM WC PORTAL
# SST DATA

# DEPTH-TEMPERATURE DST DATA (format: Date, Temperature, Depth)

B2222_noise<-read.csv("B2222_noise2.csv", sep=",")


B2222_noise$date<-as.Date(B2222_noise$Date)
B2222_noise$date <- trunc(as.Date(B2222_noise$date, origin = "1900-01-01"))
B2222_noise$Date<-B2222_noise$date
#calculate maximum depth per day
max_depth_day<-aggregate(B2222_noise$Depth, by=list(date=B2222_noise$Date), FUN=max)
max_depth_day$date<-as.Date(max_depth_day$date)

 source('bin_TempTS.R')

 B2222_noise = B2222_noise[!is.nan(B2222_noise$Temperature), ]
 B2222_noise = B2222_noise[!is.nan(B2222_noise$Depth), ]


res=5  # res defines the binning size (res=5 means 5 meter bins)
btz = bin_TempTS(B2222_noise, res = res) 
saveRDS(btz, file = "btz_2222_noise_res5.rds")
btz<-readRDS(file = "btz_2222_noise_res5.rds")
udates = unique(as.POSIXct(as.character(btz$date), tz = 'UTC'))
btz$Date = btz$date #paste0(as.character(btz$date), ' 00:00:00')


pdt<-btz
pdt.udates <- udates

# SET SPATIAL LIMITS, IF DESIRED
# these are the lat/lon bounds of your study area (e.g. where you think the animal went)
sp.lim <- list(lonmin = 12, lonmax = 17.26852, latmin = 53.8, latmax = 56.17)

if (exists('sp.lim')){
  locs.grid <- setup.locs.grid(sp.lim)
} else{
  locs.grid <- setup.locs.grid(gpe2)
  sp.lim <- list(lonmin = min(locs.grid$lon[1,]), lonmax = max(locs.grid$lon[1,]),
                 latmin = min(locs.grid$lat[,1]), latmax = max(locs.grid$lat[,1]))
}


#------------
# GET ENVIRONMENTAL DATA
#------------ 

# YOU NEED SOME REPRESENTATION OF ENVIRONMENTAL DEPTH-TEMPERATURE and bathymetry
#Baltic ROM provided by Ulf Gräwe (IOW), ncdf4 format

bathy_name <- "D:/TaBaCod/z1m/WB600m.TaBaCod.20170101.z1m.nc"  #to plot bathymetry of the area
plot(raster(bathy_name, varname="bathymetry"))
bathy<-(raster(bathy_name, varname="bathymetry"))
values(bathy) <- values(bathy)*-1


#------------
# CALCULATE LIKELIHOODS
#------------

locs.grid <- setup.locs.grid(sp.lim)

# vector indicating which likelihoods to run (e.g. 1=light, 2=sst, 5=temperature-depth, 6=maximum depth)
# can be combined with if() statements around calc functions: if (any(likVec == 5) & !exists('L.5')){calc.hycom(...)}
likVec <- c(5) 

#maximum depth likelihood
source("C:/Users/haase_s/Desktop/HMMoce/calc_depth_par.R")
hycom.dir <- "D:/z1m.test"
 L.4 <- calc.depth.par(pdt, hycom.dir, focalDim = 9, dateVec = dateVec, max_depth_day = max_depth_day)#
 saveRDS(L.4, file = "L.4_2222_noise_res1res1.rds") #save to save time when rerunning code
L.4<-readRDS(file = "L.4_2222_noise_res1res1.rds")

#temperature-depth likelihood
source("C:/Users/haase_s/Desktop/HMMoce/calc_temp_par_ncdf4.R")
hycom.dir <- "D:/z3m.test"
 L.5 <- calc.temp.par.ncdf4(pdt, hycom.dir, focalDim = 9, dateVec = dateVec, use.se = T, max_depth_day=max_depth_day, res=res) #
 save.image() # good idea to save after these larger calculations in case the next one causes problems
 gc(); closeAllConnections()
 saveRDS(L.5, file = "L.5_2222_noise_res5res5.rds")
L.5<-readRDS(file = "L.5_2222_noise_res3res3.rds")
#------------
# PREPARE TO RUN THE MODEL
#------------
# create a list of the likelihood rasters just created
L.rasters <- mget(ls(pattern = 'L\\.')) # use with caution as all workspace items containing 'L.' will be listed. We only want the likelihood outputs calculated above
# resample them all to match the most coarse layer (typically light at 1/4 deg)
# this can be changed to use whatever resolution you choose
resamp.idx <- which.max(lapply(L.rasters, FUN=function(x) raster::res(x)[1]))  #sometimes doesnt work if you run the code for 2nd time and there are other L. things still active
L.res <- resample.grid(L.rasters, L.rasters[[resamp.idx]])

# Figure out appropriate L combinations
# use this if you have a vector (likVec) indicating which likelihoods you are calculating


likVec<-c(1,2)

if (length(likVec) > 2){
  L.idx <- c(utils::combn(likVec, 2, simplify=F), utils::combn(likVec, 3, simplify=F))
} else{
  L.idx <- utils::combn(likVec, 2, simplify=F)
}

# which of L.idx combinations do you want to run?
tt<-run.idx <- c(1)#

# vector of appropriate bounding in filter. see ?hmm.filter for more info
bnd<-bndVec <- c(NA)

# vector of appropriate migr kernel speed. see ?makePar for more info.
parVec <- c(seq(0.1,0.5,0.1)) #~70km/day = 0.81 m/s


#------------
# RUN THE MODEL
#------------

distance_df <- as.data.frame(parVec) #dataframe to save distances between k´true and modelled daily positions
distance_df$max<-distance_df$mean<-distance_df$median<-distance_df$min <- NA;


    for (i in 1:length(parVec)){
      par<-parVec[i]
      ptt="B_artificial_GER_res3"
      runName <- paste(ptt,'_idx',tt,'_bnd',bnd,'_par',par,sep='')
      
      # COMBINE LIKELIHOOD MATRICES
      # L.idx combination indicates likelihood surfaces to consider
      L <- make.L(L1 = L.res[[1]][L.idx[[tt]]],    #[L.idx[[tt]]
                  L.mle.res = L.res$L.mle.res, dateVec = dateVec,
                  locs.grid = locs.grid, iniloc = iniloc, bathy = bathy,
                  pdt = pdt)
      L.mle <- L$L.mle
      L <- L$L
      g <- L.res$g
      g.mle <- L.res$g.mle
      lon <- g$lon[1,]
      lat <- g$lat[,1]
      
      # GET MOVEMENT KERNELS AND SWITCH PROB FOR COARSE GRID
      par0 <- makePar(migr.spd=par, grid=g.mle, L.arr=L.mle, p.guess=c(.9,.9), calcP=T) #0.9
      P.final <- par0$P.final
      
      # GET MOVEMENT KERNELS AND SWITCH PROB FOR FINER GRID
      par0 <- makePar(migr.spd=par, grid=g, L.arr=L, p.guess=c(.9,.9), calcP=F)
      K1 <- par0$K1; K2 <- par0$K2
      
      # RUN THE FILTER STEP  # this step takes a while
      if(!is.na(bnd)){
        f <- hmm.filter(g, L, K1, K2, maskL=T, P.final, minBounds = bnd)
        maskL.logical <- TRUE
      } else{
        f <- hmm.filter(g, L, K1, K2, P.final, maskL=F)
        maskL.logical <- FALSE
      }
      nllf <- -sum(log(f$psi[f$psi>0])) # negative log-likelihood
      
      # RUN THE SMOOTHING STEP
      source("hmm.smoother.test.R")
      s <- hmm.smoother.test(f, K1, K2, L, P.final) #takes long
      
      # GET THE MOST PROBABLE TRACK
      tr <- calc.track(s, g, dateVec, iniloc, method="mean")

      
      gps_t<- read.csv("C:/Users/haase_s/Desktop/HMMoce/track_artificial_GER_test2.csv") #read known track
      #calculate distance between knwon and modelled daily positions:
      gps_t$date<- as.Date(gps_t$date, format='%Y-%m-%d')
      gps_t <- gps_t[,c(1,3,4)]
      colnames(gps_t)<- c("date", "lat_gps", "lon_gps")
      gps_t <- unique(gps_t)
      gps_t$dist<-NA
      for(j in 2:length(gps_t$date)){
        gps_t$dist[j]<-distm(c(gps_t$lon_gps[j], gps_t$lat_gps[j]), c(gps_t$lon_gps[j-1], gps_t$lat_gps[j-1]))
      }
      track_full <- full_join(gps_t,tr, by="date")
      track_full<-track_full[order(track_full$date),]
      track_full <- track_full[!is.na(track_full$lon),]
      track_full <- unique(track_full)
      library(geosphere)
      
      for(k in 1:length(track_full$date)){
        track_full$dist[k]<-distm (c(track_full$lon_gps[k], track_full$lat_gps[k]), c(track_full$lon[k], track_full$lat[k]))
      }
      dist <- track_full%>%
        group_by(date)%>%filter(is.na(dist)|dist==min(dist))
      dist<-distinct(dist)
      summary(dist$dist)
      distance_df$min[i]<-summary(dist$dist)[[1]]
      distance_df$median[i]<-summary(dist$dist)[[3]]
      distance_df$mean[i]<-summary(dist$dist)[[4]]
      distance_df$max[i]<-summary(dist$dist)[[6]]
      
      
      
      
      plotHMM(s, tr, dateVec, ptt=runName, save.plot = T, resid=T, behav.pts=T, known=iniloc2)
      cols <- brewer.pal(9, "Blues")
      pal <- colorRampPalette(rev(cols))
      jpeg(paste(runName, ".jpg", sep=""))
      par(mfrow=c(1,1))
      plot(bathy,col = pal(40))
      par(bg = "white")
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "light grey")
      plot(bathy,col = pal(40), add=T)
      points(dist$lon_gps, dist$lat_gps, col=3, pch=20, cex=0.6)
      lines(dist$lon_gps, dist$lat_gps, col=3, pch=20, cex=0.6)
      lines(tr$lon, tr$lat, col = "black", lwd=2)
      points(tr$lon[1], tr$lat[1], bg = "green", pch = 21)
      points(iniloc2$lon, iniloc2$lat, col="yellow", lwd=2, cex=1.5)
      TT <- length(tr$lon)
      points(tr$lon[TT], tr$lat[TT], bg = "red", pch = 21)
      dev.off()
      # WRITE OUT RESULTS
      outVec <- matrix(c(ptt=ptt, minBounds = bnd, migr.spd = par,
                         Lidx = paste(L.idx[[tt]],collapse=''), P1 = P.final[1,1], P2 = P.final[2,2],
                         spLims = sp.lim[1:4], resol = raster::res(L.rasters[[resamp.idx]]),
                         maskL = maskL.logical, NLL = nllf, name = runName), ncol=15)
      # write.table(outVec,paste(dataDir, 'outVec_results.csv', sep=''), sep=',', col.names=F, append=T)
      #names(outVec) <- c('ptt','bnd','migr.spd','Lidx','P1','P2','spLims','resol','maskL','nll','name')
      res <- list(outVec = outVec, s = s, g = g, tr = tr, dateVec = dateVec, iniloc = iniloc, grid = raster::res(L.res[[1]]$L.5)[1])
 

    } 
