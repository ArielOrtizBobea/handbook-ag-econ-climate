#===============================================================================
# Description: Build the data files attached to the GitHub release "data-v1"
# from the 2021 replication package archived at CISER
# (https://doi.org/10.6077/fb1a-c376).
#
# Users do not need to run this script: 0_download_data.R fetches its output.
# It is kept here to document exactly how each archived file was produced.
#
# Files written to ../release:
# - prism_county_{tmin,tmax,ppt,bins}_1981-2020.csv.gz
#     County-level PRISM aggregates (output of 1_weather_data.R, sections 7-8),
#     gzipped copies of the 2021 CSV files.
# - corn_usa_1981-1999.csv.gz, corn_usa_2000-2020.csv.gz
#     USDA NASS Quick Stats exports of county corn yields (May 2021).
# - ghcnd-stations.txt.gz, ghcnd_2020-08_tmax_prcp.csv.gz
#     GHCN-Daily station list and the August 2020 TMAX/PRCP records from the
#     2020 file downloaded from NOAA in May 2021.
# - nlcd_prism_weights.tif
#     Share (%) of each NLCD 2016 land cover class within each PRISM grid cell
#     (output of 1_weather_data.R, section 3), 16 layers named nlcd_<class>.
#===============================================================================

# Packages (raster is only needed to read the 2021 RasterStack object)
  wants <- c("data.table","R.utils","terra","raster")
  needs <- wants[!(wants %in% installed.packages()[,"Package"])]
  if(length(needs)) install.packages(needs)
  invisible(lapply(wants, function(i) require(i, character.only=TRUE)))

# Directories
  pkg <- path.expand("~/Downloads/ru-2021-ortizbobea-1") # unzipped 2021 package
  out <- paste0(dirname(getwd()),"/release")
  dir.create(out, showWarnings = F)

# 1. County-level PRISM aggregates -----
  for (v in c("tmin","tmax","ppt","bins")) {
    f <- paste0(pkg,"/data2/PRISM_co/prism_county_",v,"_1981-2020.csv")
    gzip(f, destname=paste0(out,"/",basename(f),".gz"), remove=FALSE, overwrite=TRUE)
  }

# 2. NASS corn yields -----
  for (f in list.files(paste0(pkg,"/data/yields"), full.names=T)) {
    gzip(f, destname=paste0(out,"/",basename(f),".gz"), remove=FALSE, overwrite=TRUE)
  }

# 3. GHCN-Daily -----
  gzip(paste0(pkg,"/data/GHCN/ghcnd-stations.txt"), destname=paste0(out,"/ghcnd-stations.txt.gz"),
       remove=FALSE, overwrite=TRUE)
  w <- fread(paste0(pkg,"/data/GHCN/2020.csv"), header=FALSE, colClasses="character")
  w <- w[substr(V2,1,6)=="202008" & V3 %in% c("TMAX","PRCP")]
  fwrite(w, paste0(out,"/ghcnd_2020-08_tmax_prcp.csv.gz"), col.names=FALSE)

# 4. Land cover shares on the PRISM grid -----
  r <- readRDS(paste0(pkg,"/data2/PRISM/nlcd_prism_weights.RDS")) # RasterStack, layer i = NLCD class i
  classes <- c(11,12,21,22,23,24,31,41,42,43,52,71,81,82,90,95)
  r <- rast(r[[classes]])
  names(r) <- paste0("nlcd_",classes)
  writeRaster(r, paste0(out,"/nlcd_prism_weights.tif"), datatype="FLT8S", overwrite=TRUE,
              gdal=c("COMPRESS=DEFLATE","PREDICTOR=3"))

# Sizes
  print(file.info(list.files(out, full.names=T))[,"size",drop=F])

# The end
