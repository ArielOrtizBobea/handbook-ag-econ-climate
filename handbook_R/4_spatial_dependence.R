#===============================================================================
# Description: Accounting for spatial dependence when estimating non-linear
# effects of temperature exposure on US corn yields. The script file creates
# figure 13 in the chapter.
#===============================================================================

# ==== BLOCK: setup ====

#===============================================================================
# 1). Preliminary ------
#===============================================================================

# Clean up workspace
  rm(list=ls())

# Clean and load packages
  wants <- c("sf","RColorBrewer","splines","fixest","Matrix","data.table","spdep","splm")
  needs <- wants[!(wants %in% installed.packages()[,"Package"])]
  if(length(needs)) install.packages(needs, repos="https://cloud.r-project.org")
  lapply(wants, function(i) require(i, character.only=TRUE))
  rm(needs,wants)

# Standard errors: use the small-sample corrections of lfe::felm, which
# produced the figures in the chapter (each clustering dimension gets its own
# G/(G-1) correction, and singleton counties are kept)
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
  map$fips <- as.numeric(map$STATE)*10^3 + as.numeric(map$COUNTY)
  map$lat <- fun$labpt(map)[,2]
  map$lon <- fun$labpt(map)[,1]
  eastfips <- map$fips[ map$lon > -100 ]
  #plot(st_geometry(map), col=ifelse(map$fips %in% eastfips, "green3","grey90"), lwd=.5, main="CONUS counties east of 100th meridian west")

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

#===============================================================================
# 3). Estimation ------
#===============================================================================

# A. Prepare regression data -----

  # Settings
  season <- 4:9 # April-September
  tails  <- 0.01 # minimum density of extreme bins (top/bottom coding)

  # Data aggregation
  w <- weather[weather$month %in% season, ]
  w <- as.data.table(w)
  w1 <- w[, lapply(.SD, sum), by=.(fips, year), .SDcols=c("ppt",names(w)[grepl("bin",names(w))])]
  w2 <- w[, lapply(.SD, mean), by=.(fips, year), .SDcols=c("tmax","tmin","tmean")]
  w <- merge(w2,w1)
  rm(w1,w2)
  w <- as.data.frame(w)

  # Regression data
  regdata <- merge(yields,st_drop_geometry(map)[,c("fips","lon","lat")], by=c("fips"))
  regdata <- merge(regdata,w, by=c("fips","year"))
  regdata <- regdata[regdata$fips %in% eastfips, ] # eastern counties only

  # Drop yield equal to zero (suspicious)
  sum(regdata$yield==0) # 212 cases
  regdata <- regdata[regdata$yield>0,]

  # Prepare data
  bins <- 0:38 # bins
  df <- 7 # degrees of freedom for spline
  B <- ns(bins,df, intercept = T)
  bindata <- regdata[,grepl("bin",names(regdata))]/24 # in days
  center <- bindata[,paste0("bin",bins[-c(1,length(bins))])]
  left   <- bindata[,1:(match(paste0("bin",bins[1]), names(bindata))-1)]
  right  <- bindata[,(match(paste0("bin",bins[length(bins)]), names(bindata))+1):ncol(bindata)]
  bindata <- cbind(rowSums(left),center,rowSums(right))
  names(bindata) <- paste0("bin",bins)
  bdata <- as.matrix(bindata) %*% B
  colnames(bdata) <- paste0("ns",1:df)
  rdata <- regdata[,c("statefips","dist","lon","lat","fips","year","yield","ppt","tmin","tmax","tmean")]
  rdata <- cbind(rdata,bdata)
  rdata$year <- as.numeric(rdata$year) # solves some issues (year being taken as a factor instead of a number by some pacakges)
  dens <- colMeans(bindata)
  names(dens) <- bins
  dens <- dens/sum(dens)

# ==== END BLOCK ====

# ==== BLOCK: spatial ====

# B. Run models and store results -----

# 1. Naive SEs -----

  # Formula
  f <- paste("log(yield) ~",paste(colnames(bdata), collapse=" + "), "+ ppt + I(ppt^2) + year + I(year^2)")
  f <- paste(f,"| fips") # only FIPS fixed effects

  # Run model
  reg <- feols(as.formula(f), rdata, vcov="iid")

  # Plot error term in space
  e <- residuals(reg)
  year <- 1993
  e <- e[rdata$year==year]
  p <- map[map$fips %in% rdata$fips[rdata$year==year], ]
  p <- p[order(p$fips),]
  brks <- quantile(e, seq(0,1,.1)) # deciles
  bucket <- findInterval(e, brks)
  colors <- colorRampPalette(brewer.pal(9,"RdBu"))(length(unique(bucket))-1)
  plot(st_geometry(p), col=colors[bucket])

  # Test for spatial dependence in the error term
  coords <- fun$labpt(p)
  nb  <- knn2nb(knearneigh(coords, k = 5), row.names = p$fips) # 5-nearest neighbors
  dlist  <- nbdists(nb, coords, longlat=T)
  idlist <- lapply(dlist, function(x) 1/x) # inverse distance function
  nbw    <- nb2listw(nb, glist=idlist, style="W")
  moran.test(e, listw=nbw, alternative="two.sided") # yes, indeed, small p-val


  # Get marginal effects and SEs
  me <- B %*% cbind(coef(reg)[1:df])
  me <- me - sum(me*(dens/sum(dens))) # Exposure-center marginal effects
  se1 <- B %*% vcov(reg)[1:df,1:df] %*% t(B) # Get confidence bands
  se1 <- sqrt(diag(se1))

  # Quick view
  if (T) {
    plot(bins, me, type="l", ylim=c(-.1,.05), main= "naive")
    lines(bins, me+1.96*se1)
    lines(bins, me-1.96*se1)
  }

# 2. Heteroscedasticity Robust SEs -----

  # Get SEs (HC0, without small-sample correction)
  vcov <- vcov(reg, vcov="hetero", ssc=ssc(K.adj=FALSE))
  se2 <- B %*% vcov[1:df,1:df] %*% t(B)
  se2 <- sqrt(diag(se2))

  # Quick view
  if (T) {
    plot(bins, me, type="l", ylim=c(-.1,.05), main= "Hteroscedasticity robust")
    lines(bins, me+1.96*se2)
    lines(bins, me-1.96*se2)
  }

# 3. Clustered by state -----

  # Get SEs (FIPS fixed effects and cluster by state)
  se3 <- B %*% vcov(reg, cluster=~statefips)[1:df,1:df] %*% t(B) # Get confidence bands
  se3 <- sqrt(diag(se3))

  # Quick view
  if (T) {
    plot(bins, me, type="l", ylim=c(-.1,.05), main="clustered by state")
    lines(bins, me+1.96*se3)
    lines(bins, me-1.96*se3)
  }

# 4. Clustered by state and by year -----

  # Get SEs (FIPS fixed effects and cluster by state and year)
  se4 <- B %*% vcov(reg, cluster=~statefips+year)[1:df,1:df] %*% t(B) # Get confidence bands
  se4 <- sqrt(diag(se4))

  # Quick view
  if (T) {
    plot(bins, me, type="l",  ylim=c(-.1,.05), main="cluster by state and year")
    lines(bins, me+1.96*se4)
    lines(bins, me-1.96*se4)
  }

# 5. Conley 500 miles -----

  # Run model, keeping the demeaned regressors needed for the spatial HAC
  rdata <- rdata[order(rdata$year,rdata$fips),]
  reg <- feols(as.formula(f), rdata, vcov="iid", demeaned=TRUE)

  # Compute new covaraince matrix for various distance cutoffs (in miles): 0, 500 and 1000
  vcov.conley <- fun$conley(reg, idvec=rdata$fips, timevec=rdata$year, latvec=rdata$lat, lonvec=rdata$lon,
                        kernel = "bartlett", dist_cutoff = c(0,500,1000))

  # This gives you a list of covariance matrices
  class(vcov.conley)
  names(vcov.conley)

  # The "-1" element is simply the naive SEs
  setemp <- B %*% vcov.conley[["-1"]][1:df,1:df] %*% t(B) # Get confidence bands
  setemp <- sqrt(diag(setemp))
  head(cbind(setemp,se1,setemp-se1)) # virtually the same

  # The "0" mile cutoff is the same as heteroscedasticity robust SE
  setemp <- B %*% vcov.conley[["0"]][1:df,1:df] %*% t(B) # Get confidence bands
  setemp <- sqrt(diag(setemp))
  head(cbind(setemp,se2,setemp-se2)) # virtually the same

  # Quick view Conley errors at difference distances
  if (T) {
    lapply(names(vcov.conley)[-c(1)], function(d) {
      setemp <- B %*% vcov.conley[[paste(d)]][1:df,1:df] %*% t(B) # Get confidence bands
      setemp <- sqrt(diag(setemp))
      plot(bins, me, type="l",  ylim=c(-.1,.05), main=paste("Conley:",d,"miles"))
      lines(bins, me+1.96*setemp)
      lines(bins, me-1.96*setemp)
      setemp
    })
  }


# 6. Spatial error model (SEM) -----

  # This does not work for unbalanced panels: keep counties observed in all years
  rdata2 <- rdata[rdata$year %in% 1981:2020,]
  rdata2 <- rdata2[rdata2$fips %in% names(which(table(rdata2$fips)==length(unique(rdata2$year)))),]
  table(table(rdata2$fips)) # there's only 599 counties with 40 years
  nrow(rdata2)

  # Create spatial neighboring object
  samplemap <- map[map$fips %in% unique(rdata2$fips),]
  # Important to SORT counties in the shapefile so that they match the regression data
  # Otherwise the function will match counties with wrong neighbors
  samplemap <- samplemap[order(samplemap$fips),]
  coords <- fun$labpt(samplemap)
  nb  <- knn2nb(knearneigh(coords, k = 5), row.names = samplemap$fips) # 5-nearest neighbors
  dlist  <- nbdists(nb, coords, longlat=T)
  idlist <- lapply(dlist, function(x) 1/x) # inverse distance function
  nbw    <- nb2listw(nb, glist=idlist, style="W")

  # View neighboring relationships
  plot(st_geometry(samplemap))
  plot(nb, coords, col="red", add=T, points=F)

  # Important to have "fips" and "year" as first columns so the function works
  rdata2 <- rdata2[order(rdata2$fips, rdata2$year), c("fips","year","yield","ppt",colnames(bdata))]

  # Formula
  f <- paste("log(yield) ~",paste(colnames(bdata), collapse=" + "), "+ ppt + I(ppt^2) + I(as.numeric(year)) + I(as.numeric(year)^2)")

  # Run model
  reg <- spml(formula = as.formula(f), data = rdata2,
              listw = nbw, lag = FALSE, model="within", effect="individual", spatial.error = "b")

  # Check results
  summary(reg) # note how high is "rho" showing that the is substantial spatial dependence

  # Test for spatial dependence in the error term
  e <- residuals(reg)
  e <- e[rdata2$year==2000] # residuals for year 2000
  moran.test(e, listw=nbw, alternative="two.sided") # yes, indeed, small p-val

  # Get marginal effects and SEs
  mesp <- B %*% cbind(coef(reg)[-1][1:df])
  mesp <- mesp - sum(mesp*(dens/sum(dens))) # Exposure-center marginal effects
  vcovsp <- vcov(reg)[-1,-1] # remove covariance element for spatial autocorrelation parameter
  se6 <- B %*% vcovsp[1:df,1:df] %*% t(B) # Get confidence bands
  se6 <- sqrt(diag(se6))

  # Quick view
  if (T) {
    plot(bins, mesp, type="l",  ylim=c(-.1,.05), main="SEM")
    lines(bins, mesp+1.96*se6)
    lines(bins, mesp-1.96*se6)
  }


# 7. Organize previous results in a list -----

  # Organize ME and SEs
  rlist <- list(naive     = list(me=me, se=se1),
                robust    = list(me=me, se=se2),
                clu.st    = list(me=me, se=se3),
                clu.styr  = list(me=me, se=se4),
                conley500 = list(me=me, se=sqrt(diag(B %*% vcov.conley[["500" ]][1:df,1:df] %*% t(B)))),
                conley1000= list(me=me, se=sqrt(diag(B %*% vcov.conley[["1000"]][1:df,1:df] %*% t(B)))),
                sem       = list(me=mesp, se=se6))



  # Labels
  labels <- list(naive     = "Naive i.i.d.",
                 robust    = "Heteroscedasticity Robust",
                 clu.st    = "Clustered (State)",
                 clu.styr  = "Clustered (State + Year)",
                 conley500 = "Spatial HAC (500 mi)",
                 conley1000= "Spatial HAC (1000 mi)",
                 sem       = "SEM")

# ==== END BLOCK ====

# ==== BLOCK: fig13 DEPS: spatial ====

#===============================================================================
# 4). Visualize  ------
#===============================================================================

# Settings
  colors <- c("darkblue", adjustcolor("steelblue", alpha.f = .4))
  ylim <- c(-.09,.03)
  lwd <- 3

# Figure 13

  # Header
  fname <- paste0(dir$figures,"/fig13_spdep.png")
  png(fname, height = 1600, width = 2600, pointsize = 40)
  par(mfrow=c(2,4), mar=c(2,2,3,2), oma=c(2,2,0,0), cex.axis=1.2, lwd=2)

  # Loop over models
  lapply(names(rlist), function(type) {

    print(type)

    # Get confidence bands an marginal effects
    m <- c(rlist[[type]]$me)
    s <- c(rlist[[type]]$se)
    bot  <- c(m) - 1.96 * s
    top  <- c(m) + 1.96 * s
    bot2 <- c(m) - 2.58 * s
    top2 <- c(m) + 2.58 * s

    # Marginal effects + confidence
    plot(bins, m, type="l", axes=F, xlab="", ylab="", col=colors[1], ylim=ylim, lwd=lwd)
    polygon(x=c(bins,rev(bins)), y=c(top2,rev(bot2)), col=colors[2], border=NA)
    polygon(x=c(bins,rev(bins)), y=c(top ,rev(bot )), col=colors[2], border=NA)
    abline(h=0, lty=2)
    axis(1)
    axis(2, las=2)
    box()

    # Add distributions
    scale <- max(dens) * 5
    for (i in 1:length(dens)) rect(xleft=bins[i]-.5, xright = bins[i]+.5, ybottom = min(ylim), ytop = min(ylim) + dens[i]*scale, col="green3", border="green4")
    mtext("Growing-season distribution", side=1, font=1, cex=.8, line=-5, col="green4")
    if (type %in% c("conley500","conley1000","sem"))  mtext("Temperature bin (°C)", side=1, font=1, cex=1, line=2.5)
    mtext(labels[[type]], side=3, font=2, cex=1, line=.5)
    mtext(toupper(letters)[match(type,names(rlist))], side=3, adj=-.2, font=1, cex=2.5)

  })
  dev.off()

# ==== END BLOCK ====

# The end
