#===============================================================================
# Description: Download the data used in the other script files
#===============================================================================

# This script file downloads into ../data and ../data2:
#   - US county boundaries (US Census Bureau, 2010 cartographic boundary file)
#   - Country boundaries (Natural Earth, 1:10m)
#   - GHCN-Daily station list and August 2020 station records (NOAA)
#   - County corn yields, 1981-2020 (USDA NASS Quick Stats)
#   - PRISM daily tmax, tmin and ppt for August 2020 and monthly tmax and ppt
#     for August 2020 (PRISM Climate Group)
#   - County-level PRISM data for 1981-2020 and land cover shares over the PRISM
#     grid, which are outputs of 1_weather_data.R
#
# Default download: about 600 MB, a few minutes.
#
# Files that are not downloaded from their original provider are attached to
# the GitHub release "data-v1" of this repository. They are the versions used
# for the chapter in 2021 (see README and _make_release_assets.R):
#   - The county-level PRISM data took about 120 GB of daily PRISM files and a
#     day of computing to produce (sections 7-8 of 1_weather_data.R).
#   - The land cover shares were computed from the NLCD 2016 release of April
#     2019, which MRLC no longer distributes.
#   - The yields are a Quick Stats export (the Quick Stats API requires a key).
#   - NOAA keeps revising GHCN-Daily, so today's files differ slightly from
#     the ones used in 2021 (set ghcn_current <- TRUE to use them anyway).
#
# Set full <- TRUE to also download the PRISM files needed to rebuild the
# county-level data with 1_weather_data.R (daily tmin and tmax and monthly
# tmin, tmax and ppt for 1981-2020: about 30,000 files and 120 GB of disk).
#
# The PRISM web service allows each file to be downloaded twice per day from
# the same IP address. Files already on disk are skipped, so the script can be
# re-run after an interruption.

#===============================================================================
# 1). Preliminary ------
#===============================================================================

# Clean up workspace
  rm(list=ls())
  gc()

# Clean and load packages
  wants <- c("prism","rnaturalearth","data.table","R.utils")
  needs <- wants[!(wants %in% installed.packages()[,"Package"])]
  if(length(needs)) install.packages(needs)
  lapply(wants, function(i) require(i, character.only=TRUE))
  rm(needs,wants)

# Settings
  full <- FALSE          # TRUE to also download the PRISM files for 1981-2020
  ghcn_current <- FALSE  # TRUE to use today's GHCN-Daily files instead of the 2021 ones
  release <- "https://github.com/ArielOrtizBobea/handbook-ag-econ-climate/releases/download/data-v1"
  options(timeout = 3600) # some files are large

# Directories
  dir <- list()
  dir$root <- dirname(getwd())
  dir$map <- paste0(dir$root,"/data/gz_2010_us_050_00_20m")
  dir$world <- paste0(dir$root,"/data/naturalearth")
  dir$stations <- paste0(dir$root,"/data/GHCN")
  dir$yields <- paste0(dir$root,"/data/yields")
  dir$prism <- paste0(dir$root,"/data/PRISM")
  dir$prism2 <- paste0(dir$root,"/data2/PRISM")
  dir$data <- paste0(dir$root,"/data2/PRISM_co")
  lapply(dir, function(i) dir.create(i, recursive = T, showWarnings = F))

# Download a file unless it is already on disk, and decompress .gz files
  download <- function(url, folder, gunzip=TRUE) {
    dest <- paste0(folder,"/",basename(url))
    out <- if (gunzip) sub("[.]gz$","",dest) else dest
    if (!file.exists(out)) {
      download.file(url, dest, mode="wb")
      if (gunzip & grepl("[.]gz$",dest)) gunzip(dest, overwrite=TRUE)
    }
    out
  }


#===============================================================================
# 2). Boundaries ------
#===============================================================================

# US counties (2010 cartographic boundary file, 1:20 million)
  f <- download("https://www2.census.gov/geo/tiger/GENZ2010/gz_2010_us_050_00_20m.zip", dir$map)
  unzip(f, exdir=dir$map)

# Countries of the world (Natural Earth 1:10m, for Figure 1)
  if (length(list.files(dir$world, pattern="^ne_10m_admin_0_countries"))==0) {
    ne_download(scale=10, type="countries", category="cultural", destdir=dir$world, load=FALSE)
  }


#===============================================================================
# 3). Weather stations (GHCN-Daily) ------
#===============================================================================

  if (ghcn_current) {
    # Today's station list and 2020 records from NOAA (~170 MB), subset to the
    # August 2020 TMAX and PRCP records used in 1_weather_data.R
    download("https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt", dir$stations)
    f <- download("https://www.ncei.noaa.gov/pub/data/ghcn/daily/by_year/2020.csv.gz", dir$stations, gunzip=FALSE)
    w <- fread(f, header=FALSE, colClasses="character")
    w <- w[substr(V2,1,6)=="202008" & V3 %in% c("TMAX","PRCP")]
    fwrite(w, paste0(dir$stations,"/ghcnd_2020-08_tmax_prcp.csv.gz"), col.names=FALSE)
    rm(w)
  } else {
    # Versions downloaded from NOAA in May 2021
    download(paste0(release,"/ghcnd-stations.txt.gz"), dir$stations)
    download(paste0(release,"/ghcnd_2020-08_tmax_prcp.csv.gz"), dir$stations, gunzip=FALSE)
  }


#===============================================================================
# 4). Corn yields (USDA NASS) ------
#===============================================================================

# Quick Stats export of county yields for "CORN, GRAIN - YIELD, MEASURED IN BU / ACRE",
# 1981-2020, in two files because of the export size limit (see README)
  download(paste0(release,"/corn_usa_1981-1999.csv.gz"), dir$yields)
  download(paste0(release,"/corn_usa_2000-2020.csv.gz"), dir$yields)


#===============================================================================
# 5). PRISM grids for August 2020 ------
#===============================================================================

# Daily tmax and tmin for August 2020 and September 1 (the following day is
# needed to interpolate temperature over the last night of the month)
  for (var in c("tmax","tmin")) {
    prism_set_dl_dir(paste0(dir$prism,"/daily/",var))
    get_prism_dailys(type=var, minDate="2020-08-01", maxDate="2020-09-01", keepZip=FALSE)
  }

# Daily ppt for August 2020
  prism_set_dl_dir(paste0(dir$prism,"/daily/ppt"))
  get_prism_dailys(type="ppt", minDate="2020-08-01", maxDate="2020-08-31", keepZip=FALSE)

# Monthly tmax and ppt for August 2020
  for (var in c("tmax","ppt")) {
    prism_set_dl_dir(paste0(dir$prism,"/monthly/",var))
    get_prism_monthlys(type=var, years=2020, mon=8, keepZip=FALSE)
  }


#===============================================================================
# 6). Outputs of 1_weather_data.R archived in 2021 ------
#===============================================================================

# Share (%) of each NLCD 2016 land cover class in each PRISM grid cell (section 3)
  download(paste0(release,"/nlcd_prism_weights.tif"), dir$prism2)

# County-level PRISM data, 1981-2020 (sections 7-8)
  if (!full) {
    for (var in c("tmin","tmax","ppt","bins")) {
      download(paste0(release,"/prism_county_",var,"_1981-2020.csv.gz"), dir$data)
    }
  }


#===============================================================================
# 7). PRISM grids for 1981-2020 (only if full <- TRUE) ------
#===============================================================================

  if (full) {

    # Monthly tmin, tmax and ppt (1,440 files)
    for (var in c("tmin","tmax","ppt")) {
      prism_set_dl_dir(paste0(dir$prism,"/monthly/",var))
      get_prism_monthlys(type=var, years=1981:2020, mon=1:12, keepZip=FALSE)
    }

    # Daily tmin and tmax (29,220 files, ~120 GB, many hours)
    for (var in c("tmin","tmax")) {
      prism_set_dl_dir(paste0(dir$prism,"/daily/",var))
      get_prism_dailys(type=var, minDate="1981-01-01", maxDate="2020-12-31", keepZip=FALSE)
    }

  }

# The end
