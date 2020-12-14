calc.depth.par <- function (pdt, hycom.dir, focalDim = 9, dateVec, ncores = NULL, max_depth_day) #
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
  udates <- unique(lubridate::parse_date_time(pdt$Date, orders = "%Y-%m-%d %H%:%M:%S"))
  T <- length (udates)    
  print(paste0("Generating profile likelihood for ", udates[1], 
               " through ", udates[length(udates)]))
  
  nc1 = ncdf4::nc_open(dir(hycom.dir, full.names = T)[1], write=FALSE, readunlim=TRUE, verbose=FALSE,
                       auto_GMT=TRUE, suppress_dimvals=FALSE )  #here we load the first file! 

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
  

  ii = 1
  for (i in 1:length(dateVec)) {
    

    time <- as.Date(dateVec[i])
    print(time)
    max_depth <- max_depth_day[which(max_depth_day$date == time), ]
    bathy_name <- "D:/TaBaCod/z1m/WB600m.TaBaCod.20170101.z1m.nc"
    bathy<-(raster(bathy_name, varname="bathymetry"))
    #filter for bathymetry which is within 5 to the maximum depth recorded in the DST of that day
    values(bathy) <- (values(bathy) >=max_depth$x & values(bathy) <=max_depth$x+5) 
    values(bathy)[values(bathy) == 0] = NA
    L.hycom[, , i] = t(values(bathy))    #safe likelihood in L-hycom vector for each day
  }
  

  
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