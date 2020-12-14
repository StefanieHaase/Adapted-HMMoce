calc.temp.par.ncdf4 <- function (pdt, hycom.dir, focalDim = 9, dateVec, use.se = TRUE, ncores = NULL, max_depth_day, res) #
{
 library(ncdf4)
  library(doParallel)
   options(warn = -1)
  t0 <- Sys.time()
  print(paste("Starting Hycom profile likelihood calculation..."))
  if (is.null(ncores)) 
    ncores <- ceiling(parallel::detectCores() * 0.9)
  if (is.na(ncores) | ncores < 0) 
    ncores <- ceiling(as.numeric(system("nproc", intern = T)) * 
                        0.9)
  pdt$MidTemp <- (pdt$MaxTemp + pdt$MinTemp)/2
  dateVec = lubridate::parse_date_time(dateVec, "%Y-%m-%d")
  udates <- unique(lubridate::parse_date_time(pdt$Date, orders = "%Y-%m-%d"))
  T <- length (udates)    #length(dir(hycom.dir, full.names = T))
  print(paste0("Generating profile likelihood for ", udates[1], 
               " through ", udates[length(udates)]))
 
  nc1 = ncdf4::nc_open(dir(hycom.dir, full.names = T)[1], write=FALSE, readunlim=TRUE, verbose=FALSE,
          auto_GMT=TRUE, suppress_dimvals=FALSE )  #takes the first file from direction 

  ncnames = NULL
  nmax <- nc1$ndims 
  for (ii in 1:nmax) ncnames[ii] <- attributes(nc1$dim)$names[ii]  
 
  depth <- ncdf4::ncvar_get( nc1, attributes(nc1$dim)$names[3]) 
  lon  <- ncdf4::ncvar_get( nc1, attributes(nc1$dim)$names[1]) 
  if (length(dim(lon)) == 2) 
    lon <- lon[, 1]
  if (!any(lon < 180)) 
    lon <- lon - 360
  lat <- ncdf4::ncvar_get(nc1, attributes(nc1$dim)$names[2]) 
  if (length(dim(lat)) == 2) 
    lat <- lat[1, ]
  L.hycom <- array(0, dim = c(length(lon), length(lat), length(dateVec)))  
  print("Processing in parallel... ")
  cl = parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl, cores = ncores)
  #ans = foreach::foreach(i = 1:T) %dopar% {
  ans1<-list()
  for(i in 1:T){
      print(i)
    time <- as.Date(udates[i])
    pdt.i <- pdt[which(pdt$Date == time), ]
    nc <- ncdf4::nc_open(paste(hycom.dir, "/depthtemp_", as.Date(time), ".nc", sep = ""), write=FALSE, readunlim=TRUE, verbose=FALSE,
                         auto_GMT=TRUE, suppress_dimvals=FALSE)
    

  
    dat <- ncdf4::ncvar_get(nc, (nc$var)$temp)

    y <- pdt.i$Depth[!is.na(pdt.i$Depth)]
    y[y < 0] <- 0
    
    x <- pdt.i$MidTemp[!is.na(pdt.i$Depth)]
    depIdx = unique(apply(as.data.frame(pdt.i$Depth), 1, 
                          FUN = function(x) which.min((x - depth)^2)))
    hycomDep <- depth[depIdx]
    suppressWarnings(fit.low <- locfit::locfit(pdt.i$MinTemp ~ 
                                                 pdt.i$Depth))
    suppressWarnings(fit.high <- locfit::locfit(pdt.i$MaxTemp ~ 
                                                  pdt.i$Depth))
    n = length(hycomDep)
    pred.low = stats::predict(fit.low, newdata = hycomDep, 
                              se = T, get.data = T)
    pred.high = stats::predict(fit.high, newdata = hycomDep, 
                               se = T, get.data = T)
    if (use.se) {
      df = data.frame(low = pred.low$fit - pred.low$se.fit * 
                        sqrt(n), high = pred.high$fit + pred.high$se.fit * 
                        sqrt(n), depth = hycomDep)
    }
    else {
      df = data.frame(low = pred.low$fit, high = pred.high$fit, depth = hycomDep)
    }

    sd.i = array(NA, dim = c(dim(dat)[1:2], length(depIdx)))  #geändert, passt dsa so?
    for (ii in 1:length(depIdx)) {
      r = raster::raster(t(dat[, , depIdx[ii]]))
      f1 = raster::focal(r, w = matrix(1, nrow = focalDim,  ncol = focalDim), fun = function(x) stats::sd(x, na.rm = T))
      f1 = t(raster::as.matrix(f1))
      sd.i[, , ii] = f1
    }
    didx = base::match(udates, dateVec)
    lik.pdt = array(NA, dim = c(dim(dat)[1:2], length(depIdx)))   #auch geändert
    for (b in 1:length(depIdx)) {
      lik.try <- try(likint3(dat[, , depIdx[b]], sd.i[, , b], df[b, 1], df[b, 2]), TRUE)  #sidn hier die Dimensionen schon falsch?
      class.try <- class(lik.try)
      if (!any(which(lik.try > 0))) 
        class.try <- "try-error"
      if (class.try == "try-error" & use.se == FALSE) {
        df[b, 1] <- pred.low$fit[b] - pred.low$se.fit[b] * 
          sqrt(n)
        df[b, 2] <- pred.high$fit[b] - pred.high$se.fit[b] * 
          sqrt(n)
        lik.try <- try(likint3(dat[, , depIdx[b]], sd.i[, , b], df[b, 1], df[b, 2]), TRUE)
        class.try <- class(lik.try)
        if (!any(which(lik.try > 0))) 
          class.try <- "try-error"
        if (class.try == "try-error") {
          lik.try <- dat[, , depIdx[b]] * 0
          warning(paste("Warning: likint3 failed after trying with and without SE prediction of depth-temp profiles. This is most likely a divergent integral for ", 
                        time, "...", sep = ""))
        }
      }
      else if (class.try == "try-error" & use.se == TRUE) {
        lik.try <- dat[, , depIdx[b]] * 0
        warning(paste("Warning: likint3 failed after trying with and without SE prediction of depth-temp profiles. This is most likely a divergent integral for ", 
                      time, "...", sep = ""))
      }
      lik.pdt[, , b] <- lik.try
    }
    lik.pdt0 <- lik.pdt
    lik.pdt0[is.na(lik.pdt0)] <- 0
    use.idx <- unique(which(lik.pdt0 != 0, arr.ind = T)[, 3])  # find depth
    lik.pdt <- apply(lik.pdt[, , use.idx], 1:2, FUN = function(x) prod(x, na.rm = F))  
  

    ans1[[i]]<-lik.pdt
    }  
 # parallel::stopCluster(cl)
  didx <- base::match(udates, dateVec)
  lik.pdt <- lapply(ans1, function(x) x/max(x, na.rm = T))   
  ii = 1
  for (i in didx) {
 
  # L.hycom[, , i] = lik.pdt[[ii]]   #safe likelihood in L-hycom vector for each day
    
    
    time <- as.Date(dateVec[i])
    max_depth <- max_depth_day[which(max_depth_day$date == time), ]
    bathy_name <- "D:/TaBaCod/z1m/WB600m.TaBaCod.20170101.z1m.nc"
    bathy<-(raster(bathy_name, varname="bathymetry"))
    values(bathy) <- as.numeric(values(bathy) >=max_depth$x )
    L.hycom[, , i] = lik.pdt[[ii]]   #safe likelihood in L-hycom vector for each day
    
    ii = ii + 1
  }
  
  
  didx_diff<-setdiff(1:length(dateVec), didx)
   #
  
  #for days were there are no temperature data (e.g. when fish did not cover 3*res)
   for (i in 1:length(didx_diff)) {
     time <- as.Date(dateVec[didx_diff[i]])
     max_depth <- max_depth_day[which(max_depth_day$date == time), ]
     bathy_name <- "D:/TaBaCod/z1m/WB600m.TaBaCod.20170101.z1m.nc"
     bathy<-(raster(bathy_name, varname="bathymetry"))
     values(bathy) <- as.numeric(values(bathy) >=max_depth$x & values(bathy) <=(max_depth$x+3*res))
     values(bathy)[values(bathy) == 0] = NA
     L.hycom[, , didx_diff[i]] = t(values(bathy))   #safe likelihood in L-hycom vector for each day
   }
  
  #

  # #only filter for depth with reasonable distance to bottom (maxdepth)
  #  for (i in 1:length(didx)) {
  #    time <- as.Date(dateVec[didx[i]])
  #    max_depth <- max_depth_day[which(max_depth_day$date == time), ]
  #    bathy_name <- "E:/TaBaCod/z1m/WB600m.TaBaCod.20170101.z1m.nc"
  #    bathy<-(raster(bathy_name, varname="bathymetry"))
  #    values(bathy) <- as.numeric(values(bathy) <=max_depth$x & values(bathy) >=max_depth$x-20)
  #   L.hycom[, , didx[i]] = L.hycom[, , didx[i]]*values(bathy)   #safe likelihood in L-hycom vector for each day
  #  }

  print(paste("Making final likelihood raster..."))
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  L.hycom <- raster::brick(L.hycom, xmn = min(lon), xmx = max(lon),   #stack data
                           ymn = min(lat), ymx = max(lat), transpose = T, crs)

  names(L.hycom) = as.character(dateVec)
  t1 <- Sys.time()
  print(paste("Hycom profile calculations took ", round(as.numeric(difftime(t1, 
                                                                            t0, units = "mins")), 2), "minutes..."))
  options(warn = 2)
  return(L.hycom)
}