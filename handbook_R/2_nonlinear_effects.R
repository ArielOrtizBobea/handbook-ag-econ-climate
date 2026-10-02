#===============================================================================
# Description: Estimate non-linear effects of temperature exposure on US corn
# yields based on step functions, natural cubic splines and Chebyshev polynomials.
# The script file creates figures 8-10 in the chapter.
#===============================================================================

# ==== BLOCK: setup ====

#===============================================================================
# 1). Preliminary ------
#===============================================================================

# Clean up workspace
  rm(list=ls())

# Clean and load packages
  wants <- c("sf","RColorBrewer","splines","fixest","Matrix","data.table")
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
  # "> -100" keeps the 2,127 counties east of about 96°W. This is the sample
  # used for the figures in the chapter (see README).
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

#===============================================================================
# 3). Estimation of non-linear effects ------
#===============================================================================

# 0. Prepare data -----

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
  regdata <- merge(yields,w, by=c("fips","year"))
  regdata <- regdata[regdata$fips %in% eastfips, ] # easter counties only

  # Drop yield equal to zero (suspicious)
  sum(regdata$yield==0) # 212 cases
  regdata[which(regdata$yield==0),1:8]
  regdata <-   regdata[regdata$yield>0,]

# ==== END BLOCK ====

# ==== BLOCK: step ====

# 1. Step functions (+ Fig 8) ------

# Parameters
  bins <- 0:38 # bins
  dfs <- c(1,3,7) # width of the step functions

  # Prepare data, run regressions, and export marginal effects
  out1 <- lapply(dfs, function(df) {

    print(df)

    # Prepare data
    steps <- seq(min(bins),max(bins),df)
    idx <- rep(1:length(steps),each=df)
    idx <- idx[1:length(bins)]
    # If only one small step function at the end, aggregate with the previous step
    if (table(idx)[max(idx)]==1) idx[length(idx)] <- idx[length(idx)-1]
    B <- sparseMatrix(i=1:length(bins), j=idx, x=1)
    bindata <- regdata[,grepl("bin",names(regdata))]/24 # in days
    center <- bindata[,paste0("bin",bins[-c(1,length(bins))])]
    left   <- bindata[,1:(match(paste0("bin",bins[1]), names(bindata))-1)]
    right  <- bindata[,(match(paste0("bin",bins[length(bins)]), names(bindata))+1):ncol(bindata)]
    bindata <- cbind(rowSums(left),center,rowSums(right))
    names(bindata) <- paste0("bin",bins)
    bdata <- as.matrix(bindata) %*% B
    bdata <- as.matrix(bdata)
    colnames(bdata) <- paste0("step",unique(idx))
    rdata <- regdata[,c("statefips","dist","fips","year","yield","ppt","tmin","tmax","tmean")]
    rdata <- cbind(rdata,bdata)
    dens <- colMeans(bindata)
    names(dens) <- bins
    dens <- dens/sum(dens)

    # Formula
    f <- paste("log(yield) ~",paste(colnames(bdata), collapse=" + "), "+ ppt + I(ppt^2) + year + I(year^2)")
    f <- paste(f,"| fips")
    f <- as.formula(f)

    # Run model (county fixed effects, standard errors clustered by state and year)
    reg <- feols(f, rdata, cluster = ~statefips + year)
    summary(reg)

    # Get marginal effects
    me <- B %*% cbind(coef(reg)[1:length(unique(idx))])
    me <- as.matrix(me)

    # Exposure-center marginal effects
    me <- me - sum(me*(dens/sum(dens)))

    # Get confidence bands
    se <- B %*% vcov(reg)[1:length(unique(idx)),1:length(unique(idx))] %*% t(B)
    se <- sqrt(diag(se))

    # Quick view
    if (F) {
    plot(bins, me, type="s")
    lines(bins, me+1.96*se, type="s")
    lines(bins, me-1.96*se, type="s")
    }

    # Export
    list(df=df, me=me, se=se, dens=dens, steps=steps, idx=idx)

  })
  names(out1) <- dfs

# ==== END BLOCK ====

# ==== BLOCK: fig8 DEPS: step ====

# Figure 8
if (T) {
  # Header
  fname <- paste0(dir$figures,"/fig8_step.png")
  png(fname, width = 2200, height = 800, pointsize = 40)
  par(mar=c(1,2,2,2), oma=c(3,3,2,0), mfcol=c(1,3), lwd=2, cex.axis=1.2)

  # Loop over degrees of freedom
  lapply(dfs, function(df) {

    print(df)

    # Parameters
    colors2 <- c("darkblue", adjustcolor("steelblue", alpha.f = .4))
    lwd <- 3
    bot <- out1[[paste(df)]]$me - 1.96 * out1[[paste(df)]]$se
    top <- out1[[paste(df)]]$me + 1.96 * out1[[paste(df)]]$se
    bot2 <- out1[[paste(df)]]$me - 2.58 * out1[[paste(df)]]$se
    top2 <- out1[[paste(df)]]$me + 2.58 * out1[[paste(df)]]$se
    ylim2 <- c(-.15,.05)
    dens <- out1[[paste(df)]]$dens
    # Remove repeated values
    top <- top[!duplicated(top)]
    top2 <- top2[!duplicated(top2)]
    bot <- bot[!duplicated(bot)]
    bot2 <- bot2[!duplicated(bot2)]
    bins2 <- c(out1[[paste(df)]]$steps, max(bins))

    # Panel
    plot(bins, out1[[paste(df)]]$me, type="s", axes=F, xlab="", ylab="", col=colors2[1], ylim=ylim2, lwd=lwd+1)
    for(i in 1:(length(bins2)-1)) {
      polygon(x=c(bins2[c(i,i+1)],rev(bins2[c(i,i+1)])), y=c(top2[c(i,i)],rev(bot2[c(i,i)])), col=colors2[2], border=NA)
      polygon(x=c(bins2[c(i,i+1)],rev(bins2[c(i,i+1)])), y=c(top [c(i,i)],rev(bot[c(i,i)])), col=colors2[2], border=NA)
    }
    abline(h=0, lty=2)
    axis(1)
    axis(2, las=2)
    box()
    if (df==dfs[1]) mtext("Response function", side=2, cex=1, font=3, line=3.5)
    mtext("Temperature bin (°C)", side=1, font=1, cex=1, line=2.5)
    mtext(paste0(df,"°C step"), side=3, font=3)
    # Add distributions
    scale <- max(dens) * 10
    for (i in 1:length(dens)) rect(xleft=bins[i]-.5, xright = bins[i]+.5, ybottom = min(ylim2), ytop = min(ylim2) + dens[i]*scale, col="green3", border="green4")
    mtext("Growing-season distribution", side=1, font=1, cex=.8, line=-5, col="green4")

  })
  mtext("Step function", side=3, outer=T, cex=1.2, font=2)
  dev.off()
}

# ==== END BLOCK ====

# ==== BLOCK: spline ====

# 2. Natural cubic spline (+ Fig 9)------

# Parameters
  bins <- 0:38 # bins
  dfs <- c(3,7,12) # degrees of freedom

# Prepare data, run regressions, and export marginal effects
  out2 <- lapply(dfs, function(df) {

    print(df)

    # Prepare data
    B <- ns(bins,df, intercept = T)
    bindata <- regdata[,grepl("bin",names(regdata))]/24 # in days
    center <- bindata[,paste0("bin",bins[-c(1,length(bins))])]
    left   <- bindata[,1:(match(paste0("bin",bins[1]), names(bindata))-1)]
    right  <- bindata[,(match(paste0("bin",bins[length(bins)]), names(bindata))+1):ncol(bindata)]
    bindata <- cbind(rowSums(left),center,rowSums(right))
    names(bindata) <- paste0("bin",bins)
    bdata <- as.matrix(bindata) %*% B
    colnames(bdata) <- paste0("ns",1:df)
    rdata <- regdata[,c("statefips","dist","fips","year","yield","ppt","tmin","tmax","tmean")]
    rdata <- cbind(rdata,bdata)
    dens <- colMeans(bindata)
    names(dens) <- bins
    dens <- dens/sum(dens)

    # Formula
    f <- paste("log(yield) ~",paste(colnames(bdata), collapse=" + "), "+ ppt + I(ppt^2) + year + I(year^2)")
    f <- paste(f,"| fips")
    f <- as.formula(f)

    # Run model (county fixed effects, standard errors clustered by state and year)
    reg <- feols(f, rdata, cluster = ~statefips + year)
    summary(reg)

    # Get marginal effects
    me <- B %*% cbind(coef(reg)[1:df])

    # Exposure-center marginal effects
    me <- me - sum(me*(dens/sum(dens)))

    # Get confidence bands
    se <- B %*% vcov(reg)[1:df,1:df] %*% t(B)
    se <- sqrt(diag(se))

    # Quick view
    #plot(bins, me, type="l")
    #lines(bins, me+1.96*se)
    #lines(bins, me-1.96*se)

    # Export
    list(df=df, me=me, se=se, dens=dens)

  })
  names(out2) <- dfs

# ==== END BLOCK ====

# ==== BLOCK: fig9 DEPS: spline ====

# Figure 9
if (T) {
  # Header
  fname <- paste0(dir$figures,"/fig9_spline.png")
  png(fname, width = 2200, height = 1400, pointsize = 40)
  par(mar=c(1,2,2,2), oma=c(3,3,2,0), mfcol=c(2,3), lwd=2, cex.axis=1.2)

  # Loop over degrees of freedom
  lapply(dfs, function(df) {

    print(df)

    # Parameters
    colors  <- colorRampPalette(brewer.pal(name="Blues", n=9))(df)
    colors2 <- c("darkblue", adjustcolor("steelblue", alpha.f = .4))
    B <- ns(bins,df, intercept = T)
    lwd <- 3
    ncol <- ifelse(df>=dfs[2],ceiling(length(colors)/2),length(colors))
    bot <- out2[[paste(df)]]$me - 1.96 * out2[[paste(df)]]$se
    top <- out2[[paste(df)]]$me + 1.96 * out2[[paste(df)]]$se
    bot2 <- out2[[paste(df)]]$me - 2.58 * out2[[paste(df)]]$se
    top2 <- out2[[paste(df)]]$me + 2.58 * out2[[paste(df)]]$se
    ylim1 <- c(-.8,.8)
    ylim2 <- c(-.15,.05)
    dens <- out2[[paste(df)]]$dens

    # Top panel
    plot(bins, B[,1], type="l", axes=F, xlab="", ylab="", col=colors[1], ylim=ylim1, lwd=lwd)
    for (i in 2:length(colors)) lines(bins, B[,i], col=colors[i], lwd=lwd)
    abline(h=0, lty=2)
    axis(1)
    axis(2, las=2)
    box()
    legend("bottom", paste(1:length(colors)), col=colors, lwd=lwd+1, ncol=ncol, inset=.02, title="Basis matrix column:", cex=.8, x.intersp=.25)
    mtext(paste0("df = ",df), side=3, font=3)
    if (df==dfs[1]) mtext("Basis matrix", side=2, cex=1, font=3, line=3.5)

    # Bottom panel
    plot(bins, out2[[paste(df)]]$me, type="l", axes=F, xlab="", ylab="", col=colors2[1], ylim=ylim2, lwd=lwd+1)
    polygon(x=c(bins,rev(bins)), y=c(top2,rev(bot2)), col=colors2[2], border=NA)
    polygon(x=c(bins,rev(bins)), y=c(top ,rev(bot )), col=colors2[2], border=NA)
    abline(h=0, lty=2)
    axis(1)
    axis(2, las=2)
    box()
    if (df==dfs[1]) mtext("Response function", side=2, cex=1, font=3, line=3.5)
    mtext("Temperature bin (°C)", side=1, font=1, cex=1, line=2.5)
    # Add distributions
    scale <- max(dens) * 10
    for (i in 1:length(dens)) rect(xleft=bins[i]-.5, xright = bins[i]+.5, ybottom = min(ylim2), ytop = min(ylim2) + dens[i]*scale, col="green3", border="green4")
    mtext("Growing-season distribution", side=1, font=1, cex=.8, line=-5, col="green4")

  })
  mtext("Spline degrees of freedom", side=3, outer=T, cex=1.2, font=2)
  dev.off()
}

# ==== END BLOCK ====

# ==== BLOCK: poly ====

# 3. Chebyshev polynomial (+ Fig 10) ------

# Parameters
  bins <- 0:38 # bins
  dfs <- c(3,7,12) # polynomial degrees

# Prepare data, run regressions, and export marginal effects
  out3 <- lapply(dfs, function(df) {

    print(df)

    # Prepare data
    B <- poly(bins,df)
    bindata <- regdata[,grepl("bin",names(regdata))]/24 # in days
    center <- bindata[,paste0("bin",bins[-c(1,length(bins))])]
    left   <- bindata[,1:(match(paste0("bin",bins[1]), names(bindata))-1)]
    right  <- bindata[,(match(paste0("bin",bins[length(bins)]), names(bindata))+1):ncol(bindata)]
    bindata <- cbind(rowSums(left),center,rowSums(right))
    names(bindata) <- paste0("bin",bins)
    bdata <- as.matrix(bindata) %*% B
    colnames(bdata) <- paste0("poly",1:df)
    rdata <- regdata[,c("statefips","dist","fips","year","yield","ppt","tmin","tmax","tmean")]
    rdata <- cbind(rdata,bdata)
    dens <- colMeans(bindata)
    names(dens) <- bins
    dens <- dens/sum(dens)

    # Formula
    f <- paste("log(yield) ~",paste(colnames(bdata), collapse=" + "), "+ ppt + I(ppt^2) + year + I(year^2)")
    f <- paste(f,"| fips")
    f <- as.formula(f)

    # Run model (county fixed effects, standard errors clustered by state and year)
    reg <- feols(f, rdata, cluster = ~statefips + year)
    summary(reg)

    # Get marginal effects
    me <- B %*% cbind(coef(reg)[1:df])

    # Exposure-center marginal effects
    me <- me - sum(me*(dens/sum(dens)))

    # Get confidence bands
    se <- B %*% vcov(reg)[1:df,1:df] %*% t(B)
    se <- sqrt(diag(se))

    # Quick view
    #plot(bins, me, type="l")
    #lines(bins, me+1.96*se)
    #lines(bins, me-1.96*se)

    # Export
    list(df=df, me=me, se=se, dens=dens)

  })
  names(out3) <- dfs

# ==== END BLOCK ====

# ==== BLOCK: fig10 DEPS: poly ====

# Figure 10
if (T) {
  # Header
  fname <- paste0(dir$figures,"/fig10_poly.png")
  png(fname, width = 2200, height = 1400, pointsize = 40)
  par(mar=c(1,2,2,2), oma=c(3,3,2,0), mfcol=c(2,3), lwd=2, cex.axis=1.2)

  # Loop over degrees of freedom
  lapply(dfs, function(df) {

    print(df)

    # Parameters
    colors  <- colorRampPalette(brewer.pal(name="Blues", n=9))(df)
    colors2 <- c("darkblue", adjustcolor("steelblue", alpha.f = .4))
    B <- poly(bins,df)
    lwd <- 3
    ncol <- ifelse(df>=dfs[2],ceiling(length(colors)/2),length(colors))
    bot <- out3[[paste(df)]]$me - 1.96 * out3[[paste(df)]]$se
    top <- out3[[paste(df)]]$me + 1.96 * out3[[paste(df)]]$se
    bot2 <- out3[[paste(df)]]$me - 2.58 * out3[[paste(df)]]$se
    top2 <- out3[[paste(df)]]$me + 2.58 * out3[[paste(df)]]$se
    ylim1 <- c(-.8,.8)
    ylim2 <- c(-.15,.05)
    dens <- out3[[paste(df)]]$dens

    # Top panel
    plot(bins, B[,1], type="l", axes=F, xlab="", ylab="", col=colors[1], ylim=ylim1, lwd=lwd)
    for (i in 2:length(colors)) lines(bins, B[,i], col=colors[i], lwd=lwd)
    abline(h=0, lty=2)
    axis(1)
    axis(2, las=2)
    box()
    legend("bottom", paste(1:length(colors)), col=colors, lwd=lwd+1, ncol=ncol, inset=.02, title="Basis matrix column:", cex=.8, x.intersp=.25)
    mtext(paste0("degree = ",df), side=3, font=3)
    if (df==dfs[1]) mtext("Basis matrix", side=2, cex=1, font=3, line=3.5)

    # Bottom panel
    plot(bins, out3[[paste(df)]]$me, type="l", axes=F, xlab="", ylab="", col=colors2[1], ylim=ylim2, lwd=lwd+1)
    polygon(x=c(bins,rev(bins)), y=c(top2,rev(bot2)), col=colors2[2], border=NA)
    polygon(x=c(bins,rev(bins)), y=c(top ,rev(bot )), col=colors2[2], border=NA)
    abline(h=0, lty=2)
    axis(1)
    axis(2, las=2)
    box()
    if (df==dfs[1]) mtext("Response function", side=2, cex=1, font=3, line=3.5)
    mtext("Temperature bin (°C)", side=1, font=1, cex=1, line=2.5)
    # Add distributions
    scale <- max(dens) * 10
    for (i in 1:length(dens)) rect(xleft=bins[i]-.5, xright = bins[i]+.5, ybottom = min(ylim2), ytop = min(ylim2) + dens[i]*scale, col="green3", border="green4")
    mtext("Growing-season distribution", side=1, font=1, cex=.8, line=-5, col="green4")

  })
  mtext("Polynomial degree", side=3, outer=T, cex=1.2, font=2)
  dev.off()
}

# ==== END BLOCK ====


# The end
