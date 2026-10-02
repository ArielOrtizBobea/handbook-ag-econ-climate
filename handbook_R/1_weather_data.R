#===============================================================================
# Description: This script files performs a series of tasks related to weather
# data. Here is a description of each section of the code:
# 1- Runs preliminary loading packages, directories and functions
# 2- Import gridded and station weather, shapefiles, land cover
# 3- Aggregate land cover to the PRISM grid
# 4- Create transformation "P" matrix for gridded data aggregation
# 5- Interpolate weather station data to counties
# 6- Interpolate weather station data to PRISM grid
# 7- Compute temperature exposure bins
# 8- Aggregate all gridded data to the county level
# 9- Create figures 1-7 in the chapter
#
# By default (full <- FALSE) section 7 only computes the bins for August 2020
# (needed for Figure 7) and section 8 is skipped: the county-level data for
# 1981-2020 that these sections produce were downloaded by 0_download_data.R.
# Set full <- TRUE to rebuild them from the daily PRISM data, after downloading
# these data with full <- TRUE in 0_download_data.R.
#===============================================================================

# ==== BLOCK: setup ====

#===============================================================================
# 1). Preliminary ------
#===============================================================================

# Clean up workspace
  rm(list=ls())
  gc()

# Clean and load packages
  wants <- c("terra","sf","RColorBrewer","Matrix","data.table","R.utils","readr","spam","FNN")
  needs <- wants[!(wants %in% installed.packages()[,"Package"])]
  if(length(needs)) install.packages(needs, repos="https://cloud.r-project.org")
  lapply(wants, function(i) require(i, character.only=TRUE))
  rm(needs,wants)
  sf_use_s2(FALSE) # planar geometry on longitude/latitude, as in the 2021 code

# Rebuild the county-level data for 1981-2020? (see description above)
  full <- FALSE

# Directories
  dir <- list()
  dir$root <- dirname(getwd())
  dir$map <- paste0(dir$root,"/data/gz_2010_us_050_00_20m")
  dir$world <- paste0(dir$root,"/data/naturalearth")
  dir$prism <- paste0(dir$root,"/data/PRISM")
  dir$prism2 <- paste0(dir$root,"/data2/PRISM")
  dir$stations <- paste0(dir$root,"/data/GHCN")
  dir$data <- paste0(dir$root,"/data2/PRISM_co")
  dir$bins <- paste0(dir$root,"/data2/PRISM/bins")
  dir$landcover <- paste0(dir$root,"/data/NLCD")
  dir$figures <- paste0(dir$root,"/figures")
  lapply(dir, function(i) dir.create(i, recursive = T, showWarnings = F))

# Functions
  source("functions.R")

#===============================================================================
# 2). Import data ------
#===============================================================================

# 1. PRISM land and ID masks ------

  # Create land mask for PRISM raster (any daily file works; all share the same grid)
  s <- list.files(paste0(dir$prism,"/daily/tmax"), pattern="20200816[.]bil$", full.names=T, recursive=T)
  s <- rast(s)
  mask <- rast(s)
  values(mask) <- as.numeric(!is.na(values(s, mat=FALSE)))
  plot(mask)

  # Create a "grid cell ID" raster only for cells over land (NA otherwise)
  id <- mask
  values(id) <- ifelse(values(mask, mat=FALSE)==1, 1:ncell(mask), NA)
  plot(id)
  names(id) <- names(mask) <- "cellid"


# 2. Polygons of US counties -----

  map <- st_read(paste0(dir$map,"/gz_2010_us_050_00_20m.shp"), quiet=TRUE)
  map <- map[!(map$STATE %in% c("15","02","72")), ] # drop Hawaii, Alaska and Puerto Rico
  map <- st_transform(map, crs(s)) # make sure map has same projection as raster
  map$fips <- as.numeric(map$STATE)*10^3 + as.numeric(map$COUNTY) # county FIPS code
  world <- list.files(dir$world, pattern="^ne_10m_admin_0_countries[.](gpkg|shp)$", full.names=T)
  world <- st_read(world[1], quiet=TRUE) # map of the world
  world <- st_transform(world, crs(map))
  us <- st_union(map)

# 3. GHCND weather station data for 2020 -----

  # Stations and their coordinates

  # From readme.txt
  #IV. FORMAT OF "ghcnd-stations.txt"
  #- - - - - - - - - - - - - - - -
  #  Variable   Columns   Type
  #- - - - - - - - - - - - - - - - -
  #ID            1-11   Character
  #LATITUDE     13-20   Real
  #LONGITUDE    22-30   Real
  #ELEVATION    32-37   Real
  #STATE        39-40   Character
  #NAME         42-71   Character
  #GSN FLAG     73-75   Character
  #HCN/CRN FLAG 77-79   Character
  #WMO ID       81-85   Character
  pos <- fwf_positions(start=c(1,13,22,32,39,42,73,77,81),
                       end=c(11,20,30,37,40,71,75,79,85),
                       col_names = c("id","latitude","longitude","elevation","state","name","gsn flag","hcn/crn flag","wmo id"))
  Stations <- read_fwf(file=paste0(dir$stations,"/ghcnd-stations.txt"), skip=0, col_positions = pos, show_col_types = FALSE)
  Stations <- as.data.frame(Stations)
  Stations <- st_as_sf(Stations, coords=c("longitude","latitude"), crs=st_crs(map), remove=FALSE)
  temp <- lengths(st_intersects(Stations, map)) > 0
  stations <- Stations[temp,] # Subset to stations falling over CONUS
  rm(temp)

  # Weather station data for August 2020 (TMAX and PRCP records, see 0_download_data.R)
  Weather <- fread(paste0(dir$stations,"/ghcnd_2020-08_tmax_prcp.csv.gz"), header=FALSE)
  Weather <- as.data.frame(Weather)
  Weather <- Weather[Weather$V1 %in% stations$id, ] # Subset to CONUS
  weather <- Weather[Weather$V2=="20200816",] # subset to one of hottest days of the year
  weather <- weather[weather$V3=="TMAX",] # subset to tmax
  weather$V4 <- weather$V4/10 # decimal point is missing

# ==== END BLOCK ====


# ==== BLOCK: landcover ====

#===============================================================================.
# 3). Aggregate fine-scale 30m NLCD land cover data to PRISM grid ---------
#===============================================================================.


# Aggregation of land cover to grid and save to disk ----

  # File name for weights
  rname <- paste(dir$prism2,"/nlcd_prism_weights.tif",sep="")

  # Compute or read weights object from disk
  if(file.exists(rname)) {

    rstack <- rast(rname)

  } else { # Takes about 30 minutes

    # NLCD land cover raster (30 m), e.g. NLCD_2016_Land_Cover_L48_20190424.img,
    # placed in data/NLCD. See README for where to find it.
    nlcd <- list.files(dir$landcover, pattern="[.](img|tif)$", recursive=T, full.names=T)
    nlcd <- rast(nlcd[1])

    # Split the ID raster into 12 x 12 tiles. The PRISM grid cells of each tile
    # are converted to polygons in the projection of the land cover data,
    # which are then used to compute "zonal" statistics of land cover for each
    # PRISM grid cell
    rows <- split(1:nrow(id), cut(1:nrow(id), 12, labels=F))
    cols <- split(1:ncol(id), cut(1:ncol(id), 12, labels=F))
    tiles <- expand.grid(r=1:12, c=1:12)

    # Loop over tiles
    system.time({
      lout <- lapply(1:nrow(tiles), function(t) {

        print(paste("tile ",t,"/",nrow(tiles),sep=""))

        # 1. Polygons of the PRISM grid cells over land in the tile
        r <- id[min(rows[[tiles$r[t]]]):max(rows[[tiles$r[t]]]), min(cols[[tiles$c[t]]]):max(cols[[tiles$c[t]]]), drop=FALSE]
        if (all(is.na(values(r)))) return(NULL) # tile without land
        p <- as.polygons(r, aggregate=FALSE, na.rm=TRUE)
        p <- vect(st_transform(st_as_sf(p), crs(nlcd))) # sf rather than terra::project, whose datum shift moves the polygons by ~0.2 m

        # 2. Crop land cover to tile
        d <- crop(nlcd, p, snap="out")

        # 3. Rasterize tile of IDs on the land cover grid
        idraster <- rasterize(p, d, field="cellid")

        # 4. Count land cover pixels of each class in each PRISM grid cell
        dat <- data.table(id=values(idraster, mat=FALSE), value=values(d, mat=FALSE))
        dat[!is.na(id) & !is.na(value), .N, by=.(id, value)]

      }) # tile loop
    })

    # Merge and compute shares (%)
    out <- rbindlist(lout)
    out[, share := round(N/sum(N)*100,3), by=id]

    # Double check correct, sum to 100
    range(out[, sum(share), by=id]$V1)

    # Store in a raster, one layer per land cover class
    classes <- setdiff(sort(unique(out$value)), 0) # 0 is "unclassified"
    rstack <- rast(lapply(classes, function(k) {
      v <- rep(NA_real_, ncell(mask))
      v[unique(out$id)] <- 0
      v[out$id[out$value==k]] <- out$share[out$value==k]
      setValues(rast(mask), v)
    }))
    names(rstack) <- paste0("nlcd_",classes)

    # Save to disk
    writeRaster(rstack, rname, datatype="FLT8S", gdal=c("COMPRESS=DEFLATE","PREDICTOR=3"))

  }

# Cultivated land ----

  # Get cultivated land, pasture/hay and grasslands
  # See classes here: https://www.mrlc.gov/data/legends/national-land-cover-database-2016-nlcd2016-legend
  grass <- rstack[["nlcd_71"]]
  past  <- rstack[["nlcd_81"]]
  crop  <- rstack[["nlcd_82"]]
  all <- (past+crop+grass)
  names(all) <- "layer"

# ==== END BLOCK ====


# ==== BLOCK: pmatrix DEPS: landcover ====

#===============================================================================
# 4). Create transformation "P" matrix to aggregate raster to counties ------
#===============================================================================

# This section will be useful to generate data shown in Figure 7 of the chapter.

# 1. Extract cell IDs and cultivation weights for all counties ----

  # Prepare extraction of cell IDs and cropland+pasture+grassland weights
  # Takes ~ 90s
  system.time({
  temp <- c(id,all)
  info <- fun$extract_weights(temp, map)
  names(info) <- as.character(map$fips)
  })

  # Scale weights to 1 within each county
  info <- lapply(names(info), function(i) {
    #print(i)
    df <- info[[paste(i)]]
    df <- as.data.frame(df)
    df$w <- df$layer/sum(df$layer,na.rm=T) # make them add to 1
    df <- df[!is.na(df$w),, drop=F] # drop=F, to make sure object stays as a matrix, even with 1 row
    df$fips.order <- match(i, names(info))
    df
  })
  names(info) <- as.character(map$fips)

  # Double check
  unique(sapply(info, function(x) sum(x$w))) # weights add to 1?
  range(sapply(info, function(x) nrow(x))) # range of grid cell numbers inside counties

  # Convert to dataframe that will go into our transformation "P" matrix
  info <- do.call("rbind",info)
  info <- info[,c("cellid","fips.order","w")]

# 2. Create projection "P" matrix to extract weather information ------

  # We can use the info object to create weights to extract weather information
  # for each country using matrix algebra. We need a transformation/projection
  # matrix (named P) so that we can do the following:

  # Say you have a
  # A <- P %*% G

  # where:
  # G is a matrix with all the climate data, obtained as G <- values(g), where g is a raster stack
  # M has a dimension n x t, where n is the number of grid cells in the climate data, and t the number of time periods or raster layers
  # A is a matrix with N rows (number of counties) and t columns (number of layers in g)
  # P is our weight/projection matrix we need to construct
  # P dimensions? n rows (number of grid cells) by N columns (number of countries)

  # Let's start
  g <- c(s,s,s) # test with 3 raster layers
  G <- values(g)
  dim(G) # should be 872505 (number of grid cell in PRISM) x 3 (number of layers)

  # Create T matrix
  P <- sparseMatrix(i=info$cellid,
                    j=info$fips.order,
                    x=info$w,
                    dims=c(ncell(g), nrow=length(unique(info$fips.order))))
  colnames(P) <- as.character(map$fips)

  # Double check that all columns sum to 1
  unique(colSums(P))

  # Test
  A <- t(P) %*% G
  dim(A)
  head(A)

  # Save projection matrix to disk
  fname <- paste0(dir$prism2,"/p.RDS")
  saveRDS(P, fname)

# ==== END BLOCK ====


# ==== BLOCK: stations_counties ====

#===============================================================================
# 5). Interpolate weather stations to counties ------
#===============================================================================

# This section generates data that is shown in Figure 3 of the paper.

# The aggregation of weather station data follows a very similar pattern to
# the aggregation of gridded data. The whole point is to represent the source
# data as a matrix and to construct a sparse matrix that projects that matrix
# into the target dataset which will be interpolated. Each row of the projection
# matrix performs the local interpolation.


# Prepare data
  stations2 <- stations[stations$id %in% weather[,1],] # stations with weather data
  centroids <- fun$labpt(map) # county centroids
  weather <- weather[match(stations2$id,weather[,1]),]
  identical(weather[,1], stations2$id)

# Aggregation weights

# Stations within 0.5 degrees
  delta1 <- 3 # distance around centroid of each county in degrees
  nclose <- 1
  dist1 <- nearest.dist(st_coordinates(stations2), centroids, method = "greatcircle", delta=delta1)
  dist1 <- 1/dist1
  dist1 <- apply(dist1, 2, FUN=function(x) ifelse(rank(-x) %in% 1:nclose,x,0)) # get 5 closest distances
  dist1 <- t(apply(dist1, 2, FUN=function(x) x/sum(x,na.rm=T))) # scale to 1
  colnames(dist1) <- stations2$id
  rownames(dist1) <- map$fips
  P1 <- Matrix(dist1)

# 5 closest stations
  delta2 <- 3 # distance around centroid of each county in degrees
  nclose <- 5
  dist2 <- nearest.dist(st_coordinates(stations2), centroids, method = "greatcircle", delta=delta2)
  dist2 <- 1/dist2
  dist2 <- apply(dist2, 2, FUN=function(x) ifelse(rank(-x) %in% 1:nclose,x,0)) # get 5 closest distances
  dist2 <- t(apply(dist2, 2, FUN=function(x) x/sum(x,na.rm=T))) # scale to 1
  colnames(dist2) <- stations2$id
  rownames(dist2) <- map$fips
  P2 <- Matrix(dist2)

# Stations within 1 degree
  delta3 <- 1 # distance around centroid of each county in degrees
  dist3 <- nearest.dist(st_coordinates(stations2), centroids, method = "greatcircle", delta=delta3)
  dist3 <- 1/dist3
  dist3 <- t(apply(dist3, 2, FUN=function(x) x/sum(x,na.rm=T))) # scale to 1
  colnames(dist3) <- stations2$id
  rownames(dist3) <- map$fips
  P3 <- Matrix(dist3)


# Double check distance
  if (F) {
  i <- 55
  plot(st_geometry(map), lwd=.1)
  plot(st_geometry(map[i,]), add=T, col="black")
  plot(st_geometry(stations2[which(dist2[i,]>0),]), add=T, pch=16, cex=.3, col="cyan")
  plot(st_geometry(stations2[which(dist1[i,]>0),]), add=T, pch=16, cex=.3, col="blue")
  plot(st_geometry(stations2[which(dist3[i,]>0),]), add=T, pch=16, cex=.3, col="red")
  }

# Aggregate weather station to counties
  # target    <- projection x source
  weather.co1 <- as.vector(P1 %*% cbind(weather$V4))
  weather.co2 <- as.vector(P2 %*% cbind(weather$V4))
  weather.co3 <- as.vector(P3 %*% cbind(weather$V4))

# ==== END BLOCK ====


# ==== BLOCK: stations_grid ====

#===============================================================================
# 6). Interpolate weather stations to PRISM grid ------
#===============================================================================

# This section generates data that is shown in Figure 4 of the paper.

# Loop over weather variables
# Takes ~ 1 min
system.time({
slist <- lapply(c("tmax","ppt"), function(var) {

  print("----------------")
  print(var)
  var2 <- ifelse(var=="tmax","TMAX","PRCP")

  # Import monthly PRISM data
    flist <- list.files(paste0(dir$prism,"/monthly/",var), pattern="[.]bil$", full.names = T, recursive = T)
    flist <- flist[grepl("202008",flist)]
    m <- rast(flist)

  # Import weather stations data
    weather <- Weather[grepl("202008",Weather$V2),] # subset to month
    weather <- weather[weather$V3==var2,] # subset to variable
    weather$V4 <- weather$V4/10 # decimal point is missing

  # Form a balanced panel of weather stations
    table(weather$V1) # not all stations have 31 observations
    hist(table(weather$V1), main=var) # distribution
    keep <- table(weather$V1)==max(table(weather$V1))
    keep <- names(keep)[keep]
    length(keep)  # number of stations in the panel
    stations.balanced <- stations[stations$id %in% keep, ]

  # Get PRISM raster grid centroids
    centroids <- xyFromCell(m, 1:ncell(m))
    centroids <- centroids[values(mask, mat=FALSE)==1,]
    coords.stations <- st_coordinates(stations.balanced)

  # Get 5 nearest stations to each PRISM gridcell
    k <- 5
    nn <- get.knnx(coords.stations, centroids, k=k)
    # nn[[1]][1,] # index of 5 nearest stations
    # plot(coords.stations, pch=16, cex=.5)
    # points(coords.stations[nn[[1]][1,],], col="red")

  # Compute the inverse distance weight for each one
    w <- 1/nn[[2]] # inverse distance
    w <- w/rowSums(w) # scaled to 1
    rowSums(w)

  # Create inverse distance weighting matrix
    dist <- sparseMatrix(i= rep(1:nrow(centroids),each=k), # row index - gridcells
                         j= c(t(nn[[1]])), # column index - stations
                         x= c(t(w)),
                         dims=c(nrow(centroids),nrow(coords.stations)))
    colnames(dist) <- stations.balanced$id

  # Loop over days of the month
    s <- lapply(unique(weather$V2), function(day) { # day <- 20200801

      print(day)

      # Subset weather to day
      d <- weather[weather$V2==day,]
      d <- d[match(colnames(dist),d$V1),]

      # Map to PRISM raster grids
      g <- dist %*% cbind(d$V4)
      g <- c(as.matrix(g))

      # Store in a raster
      v <- rep(NA, ncell(mask))
      v[values(mask, mat=FALSE)==1] <- g
      r <- setValues(rast(mask), v)
      plot(r, main=paste(var,"-",day))
      r
    })
    s <- rast(s)

  # Remove interpolated mean and add monthly PRISM mean
  if (var=="ppt") {
    stot <- sum(s)
    scale <- m/stot
    s <- s*scale
    test <- sum(s)
    fill <- which(is.na(values(test, mat=FALSE)))
    fill <- fill[fill %in% which(values(mask, mat=FALSE)==1)]
    v <- values(s)
    v[fill,] <- 0
    values(s) <- v
    #plot(sum(s))
    #plot(m)
  } else {
    savg <- mean(s)
    s <- s - savg # Remove interpolated mean
    s <- s + m    # Add monthly PRISM mean
  }

  # Return
  names(s) <- unique(weather$V2)
  s
})
})
names(slist) <- c("tmax","ppt")

# ==== END BLOCK ====


# ==== BLOCK: bins ====

#===============================================================================
# 7). Create monthly temperature exposure bins ------
#===============================================================================

# This section generates data that is shown in Figure 7 of the paper.

# With full <- TRUE this takes about a day and 10GB of disk for 1981-2020.
# Otherwise it only processes August 2020 (~1 min).

# Create monthly list
  flist <- list.files(paste0(dir$prism,"/daily"), pattern="[.]bil$", full.names = T, recursive = T)
  months <- if (full) rep(1981:2020, each=12)*100 + 1:12 else 202008
  missing0 <- values(mask, mat=FALSE)==0 # # Create a logical vector to subset land
  trez <- 15/60  # the time resolution of the interpolation
  bins <- -10:50 # bins to that will be saved
  binsbase <- seq(-100.5,100.5,1) # base bins over which computation are done

# Loop over months
system.time({
lapply(months, function(mo) {

  print("--------------------------")
  timer <- system.time({
  # 1. Prepare data

    # Get files
    minfiles <- flist[grepl(mo,flist) & grepl("tmin",flist)]
    maxfiles <- flist[grepl(mo,flist) & grepl("tmax",flist)]
    ndays <- length(minfiles)
    print(paste(mo,"-",ndays,"days"))
    if (length(minfiles)!=length(maxfiles)) stop("missing daily files")
    # Add first day of the following month (or repeat the last day at the end of the record)
    mo2 <- ifelse(mo %% 100 == 12, (mo %/% 100 + 1)*100 + 1, mo + 1)
    nextmin <- flist[grepl(mo2*100+1,flist) & grepl("tmin",flist)]
    nextmax <- flist[grepl(mo2*100+1,flist) & grepl("tmax",flist)]
    if (length(nextmin)==1 & length(nextmax)==1) {
      minfiles <- c(minfiles, nextmin)
      maxfiles <- c(maxfiles, nextmax)
    } else {
      minfiles <- c(minfiles,minfiles[length(minfiles)])
      maxfiles <- c(maxfiles,maxfiles[length(maxfiles)])
    }

    # Read daily data to memory as matrices (land grid cells only)
    system.time({
    tmin <- values(rast(minfiles))[!missing0,]
    tmax <- values(rast(maxfiles))[!missing0,]
    })


  # 2. Interpolate and compute bins

    # Process 20,000 grid cells at a time, which keeps memory use low
    # (the interpolated matrix for all cells is ~11GB)
    index <- ceiling((1:nrow(tmin))/20000)
    system.time({
    m2 <- lapply(unique(index), function(j) {

      # Interpolate temperature every 15 minutes
      m <- fun$sine.approx(tmin=tmin[index==j,,drop=F], tmax=tmax[index==j,,drop=F], trez=trez)

      # Remove last 12 hours to get exactly ndays days
      remove <- (ncol(m)-((1/trez)*12)+1):ncol(m) # remove half day
      m <- m[,-remove,drop=F]
      ncol(m)/4 # 744 hours for a month with 31 days with each column representing 15 mins

      # Compute bins
      fun$exposuretobins(m=m, bins=bins, binsbase=binsbase, trez=trez)
    })
    m2 <- do.call("rbind",m2)
    })
    rm(tmin, tmax)
    gc()

  # 3. Store bins in a raster file

    # Fill matrix
    mout  <- matrix(NA, nrow=length(missing0), ncol=length(bins))
    mout[!missing0,] <- m2

    # Create raster with one layer per bin
    r <- rast(mask, nlyrs=length(bins))
    values(r) <- mout
    names(r) <- paste0("bin",bins)

    # Check it out
    # plot(r[[41]])

    # Save to disk
    fname <- paste0(dir$bins,"/bins_",mo,".tif")
    writeRaster(r, fname, overwrite=TRUE, gdal=c("COMPRESS=DEFLATE"))

  })
  print(paste(round(timer[[3]]/60,1),"min"))

})
})

# ==== END BLOCK ====


# ==== BLOCK: aggregate DEPS: pmatrix ====

#===============================================================================
# 8). Aggregate raster variables to the county level ------
#===============================================================================

# This section generates data that is used in the regression analysis in the chapter

# Loop over variables: ~ 35 mins (only with full <- TRUE)
if (full) {
system.time({
lapply(c("tmin","tmax","ppt","bins"), function(var) {

  print("---------------")
  print(var)

  # Read file list
  if (var=="bins") {
    flist <- list.files(dir$bins, pattern="[.]tif$", full.names = T)
  } else {
    flist <- list.files(paste0(dir$prism,"/monthly/",var), pattern="[.]bil$", full.names = T, recursive = T)
  }

  # Loop over file list
  flist <- flist[order(basename(flist))]
  d <- lapply(flist, function(f) { # f <- flist[1]

    #print(f)
    # Read data to a matrix
    r <- rast(f)
    G <- values(r)

    # Aggregate to county level
    A <- t(P) %*% G
    fips <- as.numeric(rownames(A))
    A <- as.matrix(A)
    if (var!="bins") colnames(A) <- var

    # Arrange for export
    ym <- regmatches(basename(f), regexpr("[0-9]{6}", basename(f)))
    year <- as.numeric(substr(ym,1,4))
    month <- as.numeric(substr(ym,5,6))
    cat(c(year*100+month), sep="...")
    o <- data.frame(fips=fips, year=year, month=month, A)

  })

  # Combine aggregated data from all files
  d <- do.call("rbind", d)
  d <- d[order(d$year,d$month,d$fips),]


  # Write to disk
  fname <- paste0(dir$data,"/prism_county_",var,"_",paste(range(d$year), collapse="-"),".csv")
  write.csv(d, fname, row.names = F)

})
})
}

# ==== END BLOCK ====

#===============================================================================
# 9). Visualizations ------
#===============================================================================

# ==== BLOCK: fig1 ====

# Figure 1 -----
# Distribution of weather stations

  fname <- paste0(dir$figures,"/fig1_map_stations.png")
  png(fname, 1800, 2200, pointsize = 45)
  par(mar=c(0,0,2,0), mfrow=c(2,1))
  plot(st_geometry(world), lwd=1, border="grey70", col="grey95")
  plot(st_geometry(Stations), pch=16, cex=.15, col="red3", add=T)
  mtext("A", side=3, adj=.02, font=1, cex=2)
  mtext("Global", side=3, font=3, cex=1.5)
  plot(st_geometry(map), lwd=1, border="grey70", col="grey95")
  plot(st_geometry(stations), pch=16, cex=.2, col="red3", add=T)
  mtext("B", side=3, adj=.02, font=1, cex=2)
  mtext("Contiguous US", side=3, font=3, cex=1.5)
  dev.off()

# ==== END BLOCK ====

# ==== BLOCK: fig2 DEPS: landcover ====

# Figure 2 -----
# Plot raster and gridcell

  # Prepare data and plotting objects
  flist <- list.files(paste0(dir$prism,"/daily/tmax"), pattern="[.]bil$", full.names = T, recursive = T)
  s <- flist[grepl(20200816,flist)]
  s   <- rast(s)
  st <- map[map$STATE=="06",] # CA
  s2 <- crop(s,st) # crop weather raster over CA
  l <- crop(crop,st) # crop cropland raster over CA
  inside <- rasterize(vect(st), s2, field=1) # get raster cells inside the state
  s2[is.na(inside)] <- NA
  l[is.na(inside)] <- NA

  # Settings for plot
  breaks <- seq(5,55,5)
  colors <- colorRampPalette(rev(brewer.pal(9,"Spectral")))(length(breaks)-1)
  breaks2 <- seq(0,1,.1)
  colors2 <- c(colorRampPalette(brewer.pal(9,"YlGn"))(length(breaks2)-1))
  fname <- paste0(dir$figures,"/fig2_map_grid.png")

  # Visualize
  png(fname, width = 1000*2, height = 1000, pointsize = 35)
  par(mar=c(2,2,2,2), oma=c(0,0,0,3), mfrow=c(1,2))
  fun$plot_raster(s2, col=colors, axes=F, box=F, legend.width=2, legend.args=list(text='Maximum temperature (°C)', side=4, font=1, line=2.5, cex=1))
  plot(st_geometry(st), lwd=1, add=T)
  mtext("A", side=3, adj=-.1, font=1, cex=2)
  fun$plot_raster(l/100, breaks=breaks2, col=colors2, axes=F, box=F, legend.width=2, legend.args=list(text='Fraction of cropland in PRISM gridcell', side=4, font=1, line=2.5, cex=1))
  plot(st_geometry(st), lwd=1, add=T)
  mtext("B", side=3, adj=-.1, font=1, cex=2)
  dev.off()

# ==== END BLOCK ====

# ==== BLOCK: fig3 DEPS: stations_counties ====

# Figure 3 -----
# Weather station and interpolated values of Tmax for August 16 2020


# Plot settings
  bgcol <- "grey30"
  breaks <- seq(5,55,5)
  bucket <- findInterval(weather$V4, breaks)
  colors <- colorRampPalette(rev(brewer.pal(9,"Spectral")))(length(breaks)-1)
  o <- match(stations2$id,weather[,1])
  colvec <- colors[bucket][o]
  leg <- paste(breaks[-length(breaks)],breaks[-1], sep="-")
  colvec1 <- colors[findInterval(weather.co1, breaks)]
  colvec1[is.na( colvec1)] <- bgcol
  colvec2 <- colors[findInterval(weather.co2, breaks)]
  colvec2[is.na( colvec2)] <- bgcol
  colvec3 <- colors[findInterval(weather.co3, breaks)]
  colvec3[is.na( colvec3)] <- bgcol

# Plot
  fname <- paste0(dir$figures,"/fig3_map_us_stations.png")
  png(fname, width = 1800, height = 900*4, pointsize =  70)
  par(mar=c(0,0,1,3), oma=c(0,0,1.5,0), xpd=T, mfrow=c(4,1))
  # Stations
  plot(st_geometry(map), lwd=.5, col=bgcol)
  plot(st_geometry(stations2), pch=21, cex=.5, bg=colvec, add=T)
  mtext("A", side=3, adj=.02, font=1, cex=2)
  mtext("Weather stations", side=3, font=3, cex=1)
  # Closest station
  plot(st_geometry(map), lwd=.5, col=colvec1)
  mtext("B", side=3, adj=.02, font=1, cex=2)
  mtext("Nearest station to centroid", side=3, font=3, cex=1)
  # 5 closest stations
  plot(st_geometry(map), lwd=.5, col=colvec2)
  mtext("C", side=3, adj=.02, font=1, cex=2)
  mtext("Interpolation of 5 nearest stations", side=3, font=3, cex=1)
  # Stations within 1 degree
  plot(st_geometry(map), lwd=.5, col=colvec3)
  mtext("D", side=3, adj=.02, font=1, cex=2)
  mtext("Interpolation of stations within 1°", side=3, font=3, cex=1)
  legend("bottomright", rev(c("n/a",leg)), pt.bg=rev(c(bgcol,colors)), pch=21, inset=c(-.1,0), title="Tmax (°C)", bty="n", cex=1, pt.cex=2)
  dev.off()

# ==== END BLOCK ====

# ==== BLOCK: fig4 DEPS: stations_grid ====

# Figure 4 -----

  # Import daily PRISM
  flist <- list.files(paste0(dir$prism,"/daily/tmax"), pattern="[.]bil$", full.names = T, recursive = T)
  stmax <- flist[grepl(202008,flist) & grepl("tmax",flist)]
  length(stmax)==31 # make sure this is 31, otherwise some daily files are missing
  flist <- list.files(paste0(dir$prism,"/daily/ppt"), pattern="[.]bil$", full.names = T, recursive = T)
  sppt  <- flist[grepl(202008,flist) & grepl("ppt",flist)]
  length(sppt)==31 # make sure this is 31, otherwise some daily files are missing
  stmax <- rast(stmax)
  sppt <- rast(sppt)
  diff.tmax <- slist$tmax[["20200816"]]-stmax[[16]]
  diff.ppt <- slist$ppt[["20200816"]]-sppt[[16]]

  # Colors
  breaks <- seq(5,55,5)
  colors <- colorRampPalette(rev(brewer.pal(9,"Spectral")))(length(breaks)-1)
  breaks2 <- c(seq(-14,-2,4),seq(2,14,4))
  colors2 <- colorRampPalette((brewer.pal(7,"RdBu")))(length(breaks2)-1)
  breaks3 <- seq(0,275,25)
  breaks3 <- c(0,1,breaks3[-1])
  colors3 <- c("white",colorRampPalette(c("aliceblue","blue","magenta"))(length(breaks3)-2))
  breaks4 <- seq(-275,275,50)
  colors4 <- colorRampPalette((brewer.pal(9,"RdBu")))(length(breaks4)-1)

  # Plot
  figname <- paste0(dir$figures,"/fig4_interpol_prism.png")
  png(figname,2400,900, pointsize=30)
  par(mfrow=c(2,3), oma=c(0,7,2,2), mar=c(0,1,1,1))

  # TMAX
  fun$plot_raster(slist$tmax[["20200816"]], breaks=breaks, col=colors, axes=F, box=F, legend.width=1.5)
  plot(us, add=T, lwd=1)
  mtext("Daily interpolation", side=3, las=1, font=3, line=0)
  mtext("TMAX (°C)", side=2, las=2, line=1, at=42, font=3)
  mtext("A", side=3, adj=.02, font=1, cex=2)

  fun$plot_raster(stmax[[16]], breaks=breaks, col=colors, axes=F, box=F, legend.width=1.5)
  plot(us, add=T, lwd=1)
  mtext("Daily PRISM", side=3, font=3, line=0)
  mtext("B", side=3, adj=.02, font=1, cex=2)

  fun$plot_raster(diff.tmax,  breaks=breaks2, col=colors2, axes=F, box=F, legend.width=1.5)
  plot(us, add=T, lwd=1)
  mtext("Difference", side=3, font=3, line=0)
  mtext("C", side=3, adj=.02, font=1, cex=2)

  # PPT
  fun$plot_raster(slist$ppt[["20200816"]], breaks=breaks3, col=colors3, axes=F, box=F, legend.width=1.5)
  plot(us, add=T, lwd=1)
  mtext("PPT (mm)", side=2, las=2, line=1, at=42, font=3)
  mtext("D", side=3, adj=.02, font=1, cex=2)

  fun$plot_raster(sppt[[16]], breaks=breaks3, col=colors3, axes=F, box=F, legend.width=1.5)
  plot(us, add=T, lwd=1)
  mtext("E", side=3, adj=.02, font=1, cex=2)

  fun$plot_raster(diff.ppt, breaks=breaks4, col=colors4, axes=F, box=F, legend.width=1.5)
  plot(us, add=T, lwd=1)
  mtext("F", side=3, adj=.02, font=1, cex=2)

  dev.off()

# ==== END BLOCK ====

# ==== BLOCK: fig5 DEPS: landcover ====

# Figure 5 -----
# Land cover fraction of cropland+pasture+grassland over each PRISM gridcell

  # Settings
  cuts <- seq(0,1,.1)
  pal <- c(colorRampPalette(brewer.pal(9,"YlGn"))(length(cuts)-1))

  # Plot
  fname <- paste0(dir$figures,"/fig5_map_landcover_prism.png")
  png(fname, width = 2000, height = 1000, pointsize = 30)
  par(mar=c(0,0,0,0), xpd=T)
  fun$plot_raster(all/max(values(all), na.rm=T), axes=F, box=F, breaks=cuts, col=pal, legend.width=1.5, legend.args = list(text = 'Fraction of cropland, pasture and grasslands in PRISM gridcell', side = 4, font = 1, line = 2.5, cex = 1))
  plot(st_geometry(map), add=T, lwd=.75)
  dev.off()

# ==== END BLOCK ====


# ==== BLOCK: fig6 ====

# Figure 6 -----
# Show a sequence of daily Tmax and Tmin over a month as well as the interpolated
# temperatures based on a double sine curve.

  # Select time period and get files
  mo <- 202008
  flist <- list.files(paste0(dir$prism,"/daily"), pattern="[.]bil$", full.names = T, recursive = T)
  minfiles <- flist[grepl(mo,flist) & grepl("tmin",flist)]
  maxfiles <- flist[grepl(mo,flist) & grepl("tmax",flist)]

  # Read daily data to memory as matrices
  tmin <- values(rast(minfiles))
  tmax <- values(rast(maxfiles))
  missing0 <- values(mask, mat=FALSE)==0 # # Create a logical vector to subset land
  trez <- 15/60  # the time resolution of the interpolation
  bins <- -10:50 # bins to that will be saved
  binsbase <- seq(-100.5,100.5,1) # base bins over which computation are done

  # Pick a random grid cell over land
  set.seed(8675309)
  i <- sample(1:sum(!missing0),1)
  tmin <- tmin[!missing0,][i,,drop=F]
  tmax <- tmax[!missing0,][i,,drop=F]

  # Interpolate (for this grid cell only)
  m <- fun$sine.approx(tmin=tmin, tmax=tmax, trez=trez)

  # Compute bins
  m2 <- fun$exposuretobins(m=m, bins=bins, binsbase=binsbase, trez=trez)

  # Settings
  cols <- c("red1","royalblue","grey50","green3")
  ylim=c(0,35)
  binrange <- c(floor(min(tmin[1,c(1:31)])): floor(max(tmax[1,c(1:31)])))

  # Plot
  if (T) {
  fname <- paste0(dir$figures,"/fig6_temperature_bins.png")
  png(fname,1400,800, pointsize=30)
  par(mar=c(3.5,3.5,1,9))
  plot(seq(1,31.5,length.out=ncol(m)),m[1,], type="p", ylab="", xlab="", axes=F, pch=16, cex=.3, ylim=ylim, col=cols[3])
  abline(h=bins, lty=2, lwd=.5)
  points(1:31,tmin[1,c(1:31)], bg=cols[2], pch=21, cex=.8)
  points(1:31+.5,tmax[1,c(1:31)], bg=cols[1], pch=21, cex=.8)
  # Outer
  axis(1, at=1:32, labels=NA, tck=-0.015)
  axis(1, at=c(1,5,10,15,20,25,31))
  axis(2, las=2)
  #axis(4, labels=F, tck=0.025)
  mtext("Temperature time distribution", col=cols[4], side=4, line=6.5, at=mean(range(ylim+c(10,0))))
  mtext("(31 days = 744 hours total)", col=cols[4], side=4, line=7.5, at=mean(range(ylim+c(10,0))), font=3)
  mtext("Day of the month", side=1, line=2)
  mtext("Temperature (°C)", side=2, line=2.5)
  box()
  legend("bottom", c("Daily tmax","Daily tmin","Sine interpolation (15 mins)","1°C bin intervals"),
         col=c(1,1,cols[3],1), pt.bg=c(cols[c(1:2)],NA,NA),
         lwd=c(NA,NA,NA,1), lty=c(NA, NA, NA, 2), pch=c(21,21,16, NA), pt.cex=c(1,1,.8,NA),
         inset=.02, ncol=2, cex=1, x.intersp=0)
  # Histogram
  par(xpd=T)
  lapply(binrange, function(x) rect(xleft = 33.25, xright = 33.25+m2[1,paste(x)]/7, ybottom = x, ytop = x+1, col=cols[4]))
  par(xpd=F)
  dev.off()
  }

# ==== END BLOCK ====


# ==== BLOCK: fig7 ====

# Figure 7 -----
# Exposure above 30C for the month of August 2020

  # Prepare data
  #r <- rast(paste0(dir$bins,"/bins_198107.tif"))
  r <- rast(paste0(dir$bins,"/bins_202008.tif")) # created in section 7
  bins <- -10:50
  above30C <- sum(r[[which(bins>=30)]])
  breakpoints <- seq(0,744,24) # 24-hour intervals up to 31 days
  colors <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(length(breakpoints)-1)

  # Plot map
  fname <- paste0(dir$figures,"/fig7_map_above30C.png")
  png(fname,2000,1000, pointsize=30)
  par(mar=c(0,0,0,0))
  fun$plot_raster(above30C,breaks=breakpoints,col=colors, box=F, axes=F, legend.width=1.5, legend.shrink=.9, legend.args=list(text='Exposure >30°C (hours)', side=4, font=2, line=2.5, cex=1))
  plot(st_geometry(map), add=T, lwd=.5)
  dev.off()

# ==== END BLOCK ====


# The end
