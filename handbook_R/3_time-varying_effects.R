#===============================================================================
# Description: Estimate time-varying non-linear effects of temperature exposure
# on US corn yields. The script file creates figures 11-12 in the chapter.
#===============================================================================

# ==== BLOCK: setup ====

#===============================================================================
# 1). Preliminary ------
#===============================================================================

# Clean up workspace
  rm(list=ls())

# Clean and load packages
  wants <- c("terra","sf","RColorBrewer","splines","fixest","Matrix","data.table")
  needs <- wants[!(wants %in% installed.packages()[,"Package"])]
  if(length(needs)) install.packages(needs)
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

# ==== BLOCK: timevarying ====

#===============================================================================
# 3). Estimation of time-varying non-linear effects ------
#===============================================================================

# 1. Preliminary steps -----

  # Settings
  bins1 <- 0:35 # temperature bins
  season <- 4:10 ##4:10 #3:9
  bins2 <- 1:length(season) # temporal bins are months
  df1 <- 6 # degrees of freedom in the temperature dimension
  df2 <- 3 # degrees of freedom in the temporal dimension

  # Bindata
  bindata <- weather[,grepl("bin",names(weather))]/24 # in days
  center <- bindata[,paste0("bin",bins1[-c(1,length(bins1))])]
  left   <- bindata[,1:(match(paste0("bin",bins1[1]), names(bindata))-1)]
  right  <- bindata[,(match(paste0("bin",bins1[length(bins1)]), names(bindata))+1):ncol(bindata)]
  bindata <- cbind(rowSums(left),center,rowSums(right))
  names(bindata) <- paste0("bin",bins1)

  # Now arrange so that monthly data is formatted in the "wide" format
  bindata$fips <- weather$fips
  bindata$year <- weather$year
  bindata$month <- weather$month
  bindata2 <- reshape(bindata, v.names=paste0("bin",bins1), timevar="month", idvar=c("fips","year"), direction="wide",sep="_m")

  # Double check
  bindata[bindata$fips==17001 & bindata$year==2020,c("month","bin20")]
  t(bindata2[bindata2$fips==17001 & bindata2$year==2020,paste0("bin20_m",1:12)])
  bindata <- bindata2
  rm(bindata2)

  # Create new tensor basis matrix
  basis1 <- as.matrix(ns(bins1, df1, intercept = T))
  basis2 <- as.matrix(ns(bins2, df2, intercept = T))
  rownames(basis1) <- paste0("bin1.",bins1)
  rownames(basis2) <- paste0("bin2.",bins2)
  basis <- basis2 %x% basis1
  dim(basis1)
  dim(basis2)
  dim(basis)

  # Drop irrelevant month "bins"
  keep <- c("fips","year", unlist(lapply(paste0("_m",season), function(x) names(bindata)[grepl(x, names(bindata))] )))
  bindata <- bindata[, names(bindata) %in% keep]
  # Modify data for 1 county and year
  mat <- bindata[bindata$fips==17001 & bindata$year==2020,-c(1,2)]
  mat <- as.matrix(mat)
  dim(mat)
  # Arrange in a rectangular format, with months as rows and temperature bins as columns
  mat2 <- matrix(mat, ncol=length(bins2), nrow=length(bins1))
  mat2 <- t(mat2)
  # Double check first month of the season
  x1 <- unlist(mat2[1,])
  x2 <- unlist(mat[,paste0("bin",bins1,"_m",season[1])])
  cbind(x1,x2) # the same, so we're ok

  # Transforming the bindata to lower dimensional data can be done as:
  bdata <- t(basis2) %*% mat2 %*% basis1
  colnames(bdata) <- paste0("h.",colnames(bdata)) # heat (h) dimension reduced from 38 to df1
  rownames(bdata) <- paste0("t.",rownames(bdata)) # temporal (t) dimension reduced from 12 to df2
  dim(bdata)
  bdata

  # Note this can be done directly with the tensor basis matrix obtained with the kronecker product
  bdata2 <- mat %*% basis
  dim(bdata2)
  bdata2

  # Compare, double check
  cbind(c(bdata2), c(t(bdata)))

  # This means that one can name the columns of bdata2
  colnames(bdata2) <- paste0("h",1:df1,".t", rep(1:df2, each=df1))
  bdata2

  # Confirm naming is correct
  t(bdata2)
  bdata

# 2. Prepare regressions data -----

  # Do it for all variables
  bindata2 <- as.matrix(bindata[,-c(1,2)]) %*% basis
  dim(bindata2) # we now have df1-by-df2 regressors
  colnames(bindata2) <- paste0("h",1:df1,".t", rep(1:df2, each=df1)) # assign name
  bindata2 <- cbind(data.frame(fips=bindata$fips, year=bindata$year), bindata2)

  # Compute density
  dens <- colMeans(bindata[,-c(1,2)])

  # Regression data
  regdata <- merge(yields,bindata2, by=c("fips","year"))
  regdata <- regdata[regdata$fips %in% eastfips, ] # eastern counties only
  regdata <- regdata[regdata$statefips %in% c(17,18,19,39), ] # Illinois, Indiana, Iowa and Ohio only

  # Drop yield equal to zero
  sum(regdata$yield==0) # 212 cases
  regdata <-   regdata[regdata$yield>0,]

# 3. Regression -----

  # Formula
  f <- paste("log(yield) ~",paste(names(bindata2)[-c(1,2)], collapse=" + "), "+ year + I(year^2)")
  f <- paste(f,"| fips")
  f <- as.formula(f)

  # Run model (county fixed effects, standard errors clustered by state and year)
  reg <- feols(f, regdata, cluster = ~statefips + year)
  summary(reg)


# 4. Post-estimation exploration -----

  # Get marginal effects and covariance
  beta <- coef(reg)[names(bindata2)[-c(1,2)]]
  me <- basis %*% beta
  var <- vcov(reg)[names(bindata2)[-c(1,2)],names(bindata2)[-c(1,2)]]
  var <- basis %*% var %*% t(basis)
  se <- sqrt(diag(var))

  # But note the marginal effects are organized in a vector
  # We need to rearrange them in the 2-D space for plotting
  me2d <- matrix(me, nrow=length(bins2), ncol=length(bins1), byrow = T)
  se2d <- matrix(se, nrow=length(bins2), ncol=length(bins1), byrow = T)
  dens2d <- matrix(dens, nrow=length(bins2), ncol=length(bins1), byrow = T)

  # Add names to dimensions
  dim(me2d)
  rownames(me2d) <- paste0("month",season)
  colnames(me2d) <- paste0("bin",bins1)
  rownames(se2d) <- paste0("month",season)
  colnames(se2d) <- paste0("bin",bins1)
  rownames(dens2d) <- paste0("month",season)
  colnames(dens2d) <- paste0("bin",bins1)

  # Quick view
  image(me2d)
  image(se2d)
  image(dens2d)

  # Re-center response by exposure
  for (i in 1:nrow(me2d)) {
    me2d[i,] <-  me2d[i,] - sum(me2d[i,]*(dens2d[i,]/sum(dens2d[i,])))
  }

  # Recheck marginal effects in 2D
  image(me2d)

  # Quick visualization - marginal effects + SEs
  lapply(season, function(sea) {
    s <- paste0("month",sea)
    plot(bins1, me2d[s,], type="l", lwd=2, ylim=c(-.2,.2), main=month.abb[sea], col=2)
    lines(bins1, me2d[s,]+1.96*se2d[s,], type="l", lwd=1, col=2)
    lines(bins1, me2d[s,]-1.96*se2d[s,], type="l", lwd=1, col=2)
    abline(h=0, lty=2)
    lines((dens2d[s,]/25)-.2, col="green3") # density
  })

  # Quick visualization in 2D
  breaks <- c(seq(-.1,.1-.025,.025),.4)
  colors <- colorRampPalette(brewer.pal(9,"RdBu"))(length(breaks)-1)
  image(me2d, axes=F, breaks=breaks, col=colors)
  axis(1, (bins2-1)/max(bins2-1), month.abb[season])
  axis(2, (bins1-1)/max(bins1-1), bins1, las=2)
  box()

  # Add a black block over effects that are not statistically significant
  signif <- (me2d+se2d*1.96 < 0) | (me2d-se2d*1.96 > 0)
  #signif[signif] <- NA
  #image(!signif*1, add=T, col=adjustcolor("black", alpha.f = .7))

  # Density
  image(dens2d, axes=F)
  axis(1, (bins2-1)/max(bins2-1), month.abb[season])
  axis(2, (bins1-1)/max(bins1-1), bins1, las=2)
  box()

# ==== END BLOCK ====

# ==== BLOCK: fig11 DEPS: timevarying ====

# 5. Figure 11 -----

  # Settings
  a <- t(me2d[,ncol(me2d):1])
  d <- t(dens2d[,ncol(dens2d):1])
  breaks <- seq(-.1,.1,.025)
  colors <- colorRampPalette(brewer.pal(9,"RdBu"))(length(breaks)-1)
  sel <- bins1 %in% seq(0,35,5)
  breaks2 <- seq(0,2.5,.25)
  colors2 <- c("grey95",colorRampPalette(brewer.pal(9,"YlGn"))(length(breaks2)-2))
  xvec <- bins2/max(bins2)
  yvec <- (bins1)/(max(bins1)+1)

  # Header
  fname <- paste0(dir$figures,"/fig11_timevarying.png")
  png(fname, height = 1600, width = 1100, pointsize = 28)
  par(mar=c(3,5,1,2), oma=c(0,.5,2,2.5), mfcol=c(2,1), lwd=2, cex.axis=1.2)
  if (T) {

  # Top panel: marginal effects
  fun$plot_raster(rast(a, extent=ext(0,1,0,1)), axes=F, breaks=breaks, col=colors, legend.width=1.5, legend.shrink=.9, legend.args=list(text='Marginal effect', side=4, font=1, line=4.5, cex=1.5))
  axis(1, (bins2-.5)/max(bins2), month.abb[season])
  axis(2, (bins1[sel]+.5)/(max(bins1)+1), bins1[sel], las=2)
  mtext("Temperature bin (°C)", side=2, line=3, cex=1.5)
  mtext("A", side=3, adj=-.3, font=1, cex=3)
  abline(v=xvec, lwd=.5)
  abline(h=yvec, lwd=.5)
  box()

  # Add significant effects
  for (i in 1:nrow(signif)) {
    pos <- which((signif[i,])==T)
    x <- xvec[i]-diff(xvec)[1]/2
    x <- rep(x, length(pos))
    if (length(x)>0) points(x,yvec[pos]+diff(yvec)[1]/2, pch="*", cex=1)
  }

  # Bottom panel: density
  fun$plot_raster(rast(d, extent=ext(0,1,0,1)), axes=F, breaks=breaks2, col=colors2, legend.width=1.5, legend.shrink=.9, legend.args=list(text='Density (days)', side=4, font=1, line=4.5, cex=1.5))
  axis(1, (bins2-.5)/max(bins2), month.abb[season])
  axis(2, (bins1[sel]+.5)/(max(bins1)+1), bins1[sel], las=2)
  mtext("Temperature bin (°C)", side=2, line=3, cex=1.5)
  mtext("B", side=3, adj=-.3, font=1, cex=3)
  abline(v=xvec, lwd=.5)
  abline(h=yvec, lwd=.5)
  box()
  }
  dev.off()

# ==== END BLOCK ====

# ==== BLOCK: fig12 DEPS: timevarying ====

# 6. Figure 12 -----

  # Settings
  colors <- c("darkblue", adjustcolor("steelblue", alpha.f = .4))
  ylim <- c(-.15,.15)

  # Header
  fname <- paste0(dir$figures,"/fig12_timevarying.png")
  png(fname, height = 1400, width = 2600, pointsize = 40)
  par(mfrow=c(2,4), mar=c(1,2,2,2), oma=c(3,2,0,0), cex.axis=1.2, lwd=2)

  # Loop over months of the season
  lapply(season, function(sea) {

    s <- paste0("month",sea)
    top <- me2d[s,]+1.96*se2d[s,]
    bot <- me2d[s,]-1.96*se2d[s,]
    top2 <- me2d[s,]+2.58*se2d[s,]
    bot2 <- me2d[s,]-2.58*se2d[s,]

    # Marginal effect + confidence
    plot(bins1, me2d[s,], type="l", lwd=4, axes=F, xlab="", ylab="", col=colors[1], ylim=ylim)
    polygon(x=c(bins1,rev(bins1)), y=c(top2,rev(bot2)), col=colors[2], border=NA)
    polygon(x=c(bins1,rev(bins1)), y=c(top ,rev(bot )), col=colors[2], border=NA)
    abline(h=0, lty=2)
    axis(1)
    axis(2, las=2)
    box()
    if (sea %in% season[(ceiling(length(season)/2)+1):length(season)]) mtext("Temperature bin (°C)", side=1, font=1, cex=1, line=2.5)

    # Density
    scale <- max(dens2d) * .01
    for (i in 1:length(dens2d[s,])) rect(xleft=bins1[i]-.5, xright = bins1[i]+.5, ybottom = min(ylim), ytop = min(ylim) + dens2d[s,i]*scale, col="green3", border="green4")
    mtext("Growing-season distribution", side=1, font=1, cex=.8, line=-4, col="green4")
    mtext(month.name[sea], side=3, line=-2, font=3)
  })
  plot(1, axes=F, type="n", xlab="", ylab="")
  dev.off()

# ==== END BLOCK ====

# The end
