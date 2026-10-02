#===============================================================================
# Description: Robustness checks. The script file creates figure 14 in the
# chapter.
#===============================================================================

# ==== BLOCK: setup ====

#===============================================================================
# 1). Preliminary ------
#===============================================================================

# Clean up workspace
  rm(list=ls())

# Clean and load packages
  wants <- c("sf","RColorBrewer","fixest","data.table")
  needs <- wants[!(wants %in% installed.packages()[,"Package"])]
  if(length(needs)) install.packages(needs)
  lapply(wants, function(i) require(i, character.only=TRUE))
  rm(needs,wants)

# Standard errors: use the small-sample corrections of lfe::felm, which
# produced the figures in the chapter (singleton counties are kept, which
# also matters for the adjusted R2 used to sort the models)
  setFixest_ssc(ssc(G.df = "conventional"), c("cluster","twoway"))
  setFixest_estimation(fixef.rm = "none")

# Directories
  dir <- list()
  dir$root <- dirname(getwd())
  dir$yields <- paste0(dir$root,"/data/yields")
  dir$prism <- paste0(dir$root,"/data2/PRISM_co")
  dir$figures <- paste0(dir$root,"/figures")
  dir$map <- paste0(dir$root,"/data/gz_2010_us_050_00_20m")
  lapply(dir, function(i) dir.create(i, recursive = T, showWarnings = F))

# Functions
  source("functions.R")

#===============================================================================
# 2). Import and manipulate data ------
#===============================================================================

# Map of US counties
  map <- st_read(paste0(dir$map,"/gz_2010_us_050_00_20m.shp"), quiet=TRUE)
  map <- map[!(map$STATE %in% c("15","02","72")), ] # drop Hawaii, Alaska and Puerto Rico
  map <- st_transform(map, "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs") # nice projection
  map$fips <- as.numeric(map$STATE)*10^3 + as.numeric(map$COUNTY)
  # Note: the map is projected, so the centroid x-coordinate is in meters and
  # "> -100" keeps the counties east of about 96°W (see README)
  eastfips <- map$fips[ fun$labpt(map)[,1] > -100 ]
  #plot(st_geometry(map), col=ifelse(map$fips %in% eastfips, "green3","grey90"), lwd=.5)

# Corn yield data
  flist <- list.files(dir$yields, full.names=T)
  yields <- lapply(flist, function(f) {
    d <- read.csv(f)
    d$statefips <- d$State.ANSI
    d$fips <-  d$statefips*10^3 +  d$County.ANSI
    d$dist <-  d$statefips*10^3 +  d$Ag.District.Code
    d$yield <- d$Value
    d$year <- d$Year
    d <- d[,c("statefips","dist","fips","year","yield")]
    d
  })
  yields <- do.call("rbind",yields)
  yields <- yields[order(yields$fips, yields$year),]

# County-level PRISM data
  flist <- list.files(dir$prism, full.names = T, recursive = T)
  weather <- lapply(flist, function(f) as.data.frame(fread(f)))
  names(weather) <- sapply(strsplit(basename(flist),"_"), function(x) rev(x)[2])
  identical(weather[[1]][,1:3],weather[[2]][,1:3])
  identical(weather[[1]][,1:3],weather[[3]][,1:3])
  identical(weather[[1]][,1:3],weather[[4]][,1:3])
  weather <- data.frame(weather$ppt, tmax=weather$tmax[,-c(1:3)], tmin=weather$tmin[,-c(1:3)], weather$bins[,-c(1:3)])
  names(weather) <- gsub("[.]","-",names(weather))
  weather$tmean <- (weather$tmax + weather$tmin)/2
  weather <- weather[,c("fips", "year", "month","ppt","tmax","tmin","tmean",names(weather)[grepl("bin",names(weather))])]

# ==== END BLOCK ====

# ==== BLOCK: specifications ====

#===============================================================================
# 3). Estimation of non-linear effects ------
#===============================================================================

# 1. List models -----

  models <- expand.grid(tvar    = c("tmax","tmean","tmin"),
                        precip  = c(TRUE,FALSE),
                        form    = c("quad","cubic"),
                        season  = c("maraug","aprsep","annual"),
                        trend   = c("pooled","state"), stringsAsFactors = F)

# 2. Loop over models (~30s)
  system.time({
  out <- lapply(1:nrow(models), function(i) {

    print("------------------------------")
    print(paste(i,"/",nrow(models)))
    print(models[i,])

    # A. Prepare regression formula -----

      # Dep variable
      lhs <- "log(yield)"
      # Regressors
      tvars <- tvar <- models$tvar[i]
      if (models$form[i]=="quad")  tvars <- unlist(lapply(tvars, function(x) c(x,paste0(x,"_sq"))))
      if (models$form[i]=="cubic") tvars <- unlist(lapply(tvars, function(x) c(x,paste0(x,"_sq"),paste0(x,"_cu"))))
      if (models$prec[i] & models$form[i]=="quad") rhs <- c(tvars,"ppt","ppt_sq")
      if (models$prec[i] & models$form[i]=="cubic") rhs <- c(tvars,"ppt","ppt_sq","ppt_cu")
      if (!models$prec[i]) rhs <- tvars
      rhs <- paste(rhs, collapse=" + ")
      # Trend
      if (models$trend[i]=="pooled") trend <- "I(year) + I(year^2)"
      if (models$trend[i]=="state") trend <- "I(year):as.factor(statefips) + I(year^2):as.factor(statefips)"
      # Fixed effects (standard errors clustered by state below)
      fe <- "| fips"
      # Put formula together
      f <- as.formula(paste(lhs,"~",rhs,"+",trend,fe))

    # B. Aggregate data to the right season ------

      # Season months. Note: as in the 2021 code behind Figure 14, "maraug"
      # covers April-August and "aprsep" covers March-September, while the
      # labels of the figure read March-August and April-September (see README)
      if (models$season[i]=="aprsep") season <- 3:9 # April-September
      if (models$season[i]=="maraug") season <- 4:8 # March-August
      if (models$season[i]=="annual") season <- 1:12 # Full year

      # Seasonal data aggregation
      w <- weather[weather$month %in% season, ]
      w <- as.data.table(w)
      w1 <- w[, lapply(.SD, sum), by=.(fips, year), .SDcols=c("ppt")]
      w2 <- w[, lapply(.SD, mean), by=.(fips, year), .SDcols=c("tmax","tmin","tmean")]
      w <- merge(w2,w1)
      rm(w1,w2)
      w <- as.data.frame(w)

      # Data with 2C warming
      w$tmax2  <- w$tmax + 2
      w$tmin2  <- w$tmin + 2
      w$tmean2 <- w$tmean + 2

      # Compute cubics and quadratics for all temperature and preciptiation variables
      for (j in c("tmax","tmean","tmin","ppt","tmax2","tmean2","tmin2")) {
        w[,paste0(j,"_sq")] <- w[,j]^2
        w[,paste0(j,"_cu")] <- w[,j]^3
      }
      w <- as.data.table(w)

      # Compute climatology
      wclim <- w[, lapply(.SD, mean), by=.(fips), .SDcols=names(w)[-c(1,2)]]
      w <- as.data.frame(w)
      wclim <- as.data.frame(wclim)

      # Regression data
      regdata <- merge(yields,w, by=c("fips","year"))
      regdata <- regdata[regdata$fips %in% eastfips, ] # eastern counties only
      sum(regdata$yield==0)       # Drop yield equal to zero (suspicious)
      regdata <-   regdata[regdata$yield>0,]

    # C. Run regression ------

      # Regression
      reg <- feols(f, regdata, cluster = ~statefips)
      summary(reg)

      # Temperature coefficients
      sel <- grepl(tvar,names(coef(reg)))
      beta <- coef(reg)[sel]
      beta <- t(t(beta))
      vcov <- vcov(reg)[sel,sel]

    # D. Post-estimation -----

      # Climates
      clim0 <- wclim[,tvars]
      clim2 <- wclim[,gsub(tvar,paste0(tvar,"2"),tvars)]
      delta <- as.matrix(clim2 - clim0)

      # Impact at county level
      impact <- delta %*% beta
      meanimp <- colMeans(impact)
      se <- delta %*% vcov %*% t(delta)
      se <- sqrt(diag(se))
      up <- impact + 1.96*se
      dw <- impact - 1.96*se
      mean(se)

      # Convert to percentage
      impact<- (exp(mean(impact))-1)*100
      up    <- (exp(mean(up))-1)*100
      dw    <- (exp(mean(dw))-1)*100

    # E. Export -----

      c(impact=impact, up=up, dw=dw, r2adj=r2(reg, "ar2"))

  })
  })

  out <- data.frame(do.call("rbind",out))
  names(out) <- c("impact","up","dw","r2adj")

# ==== END BLOCK ====

# ==== BLOCK: fig14 DEPS: specifications ====

#===============================================================================
# 4). Plot results ------
#===============================================================================

# Arrange data in propoer format for the function

  o1 <- model.matrix(~ factor(models$tvar) - 1)
  colnames(o1) <- sapply(strsplit(colnames(o1),")"), function(x) rev(x)[1])

  o2 <- model.matrix(~ factor(models$precip) - 1)
  colnames(o2) <- c("noprecip","precip")
  o2 <- o2[,c(2,1)]

  o3 <- model.matrix(~ factor(models$form) - 1)
  colnames(o3) <- sapply(strsplit(colnames(o3),")"), function(x) rev(x)[1])
  o3 <- o3[,c(2,1)]

  o4 <- model.matrix(~ factor(models$season) - 1)
  colnames(o4) <- sapply(strsplit(colnames(o4),")"), function(x) rev(x)[1])
  o4 <- o4[,3:1]

  o5 <- model.matrix(~ factor(models$trend) - 1)
  colnames(o5) <- sapply(strsplit(colnames(o5),")"), function(x) rev(x)[1])


  o <- cbind(o1,o2,o3,o4,o5)
  o <- o==1
  rm(o1,o2,o3,o4,o5)

  data <- cbind(out[,c("impact","up","dw")],o)

# Create list of labels

  labels <- list("Temperature:" = c("Tmax","Tmean","Tmin"),
                 "Precipitation:" = c("Yes","No"),
                 "Functional form:" = c("Quadratic","Cubic"),
                 "Growing season:" = c("March-August","April-September","Full year"),
                 "Time trend:" = c("Pooled","By state"))

# First try
  par(oma=c(1,0,1,1))
  fun$schart(data, index.ci=c(2,3))

# Try with labels
  fun$schart(data, index.ci=c(2,3), labels=labels)

# Try with labels
  fun$schart(data, index.ci=c(2,3), labels=labels, ylab="Impact of a 2°C warming (%)")

# Highlight baseline model
  fun$schart(data, index.ci=c(2,3), highlight=1, labels=labels, ylab="Impact of a 2°C warming (%)")


# Sort data by R2
  data2 <- data[order(out$r2adj),] # sorted dataframe
  ref <- which(rownames(data[order(out$r2adj),])=="2") # position of reference model in original dataframe

# Write to disk
  fname <- paste0(dir$figures,"/fig14_specchart.png")
  png(fname, width = 1600, height = 1000, pointsize = 25)
  par(oma=c(1,0,1,3)) #, lwd=3)
  fun$schart(data2, highlight=ref, index.ci=c(2,3), n=9, labels=labels, ylab="Impact of a 2°C warming (%)",
             lwd.est=5,lwd.symbol=1, pch.est=16, col.est=c("steelblue","red2"),
             order="asis", adj=0, offset=c(13,12), leftmargin=7)
  axis(4, las=2)
  dev.off()

# ==== END BLOCK ====


# The end
