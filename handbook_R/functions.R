#===============================================================================
# Functions used in the other script files. They are stored in a list named
# "fun" and called as fun$name().
#===============================================================================

fun <- list()

# Fit sine curve over min and max values
fun$sine.approx <- function( tmin, tmax, trez ) {
  B <- (2*pi)/24  # period = 24 hours
  C <- pi/2 # horizontal shift
  tmin <- t( tmin )
  tmax <- t( tmax )
  idx <- seq.int(  24 * 4 * nrow( tmax ) - 12 * 4 )
  xout <- ( idx - 1 ) * trez
  cycles <- sin( B * xout - C )
  mag <- matrix( NA, ncol=ncol( tmax ), nrow=length( idx ) )
  idxup <- 2 * seq.int( nrow( tmax ) ) - 1
  idxdn <- idxup[ -1 ] - 1
  idxups <- rep( seq.int( nrow( tmax ) ), each=48 )
  idxdns <- rep( seq.int( nrow( tmax ) - 1 ), each=48 )
  idxupm <- c( outer( 1:48, 48*( idxup - 1 ), FUN="+" ) )
  idxdnm <- c( outer( 1:48, 48*( idxdn - 1 ), FUN="+" ) )
  mag[ idxupm, ] <- ( tmax - tmin )[ idxups, , drop=FALSE ] / 2
  mag[ idxdnm, ] <- ( tmax[ -nrow( tmax ), , drop=FALSE ]
                      - tmin[ -1, , drop=FALSE ] )[ idxdns, , drop=FALSE ] / 2
  mag <- mag * rep( cycles, times=ncol( mag ) )
  mag[ idxupm, ] <- mag[ idxupm, ] + ( tmax + tmin )[ idxups, , drop=FALSE ] / 2
  mag[ idxdnm, ] <- mag[ idxdnm, ] + ( tmax[ -nrow( tmax ), , drop=FALSE ]
                                       + tmin[ -1, , drop=FALSE ] )[ idxdns, , drop=FALSE ] / 2
  t( mag )
}

# Compute exposure to bins for a matrix
fun$exposuretobins <- function(m, bins, binsbase, trez){
  # set up bins
  nbins <- length(binsbase) - 1L
  binsmid <- (binsbase[-1] + binsbase[-length(binsbase)])/2 # bin midpoints
  # compute bin counts
  M <- array(findInterval(m, binsbase,rightmost.closed=TRUE), dim = dim(m))
  M <- t(apply(M, 1, tabulate, nbins = nbins))
  M <- M * trez # scale to hours
  # aggregate extreme bins
  left <- rowSums(M[,1: match(bins[1], binsmid), drop=FALSE])
  right<- rowSums(M[,match(bins[length(bins)], binsmid): ncol(M), drop=FALSE])
  mid  <- M[,match(bins[2], binsmid):match(bins[length(bins)-1], binsmid), drop=FALSE]
  # export
  mout <- cbind(left,mid,right)
  colnames(mout) <- bins
  mout
}

# Spatial HAC (Conley)
fun$conley <- function(reg, idvec, timevec, latvec, lonvec, kernel = "bartlett",dist_cutoff = 500, ncores=1) {
  # Syntax: conley(reg, data$id, data$year, data$lat, data$lon, kernel = "bartlett",dist_cutoff = seq(0,500,100), ncores=1)
  # reg: a fixest object estimated with feols(..., demeaned = TRUE)
  # idvec: a vector of location IDs (n x t), in the rows of the estimation data
  # timevec: a vector of time IDs (n x t)
  # latvec: a vector of latitute (n x t)
  # lonvec: a vector of longitude (n x t)
  # kernel: one of "uniform" and "bartlett"
  # dist_cutoff: a vector of distance cutoffs (miles)

  # Internal functions (note: distance is in miles)
  weightfun <- list(uniform=function(dist, cut) {dist<=cut}, bartlett=function(dist, cut) {w <- dist<=cut ; (1-dist/(cut+1)) * w } )
  iterateObs   <-function(e1,X1,fordist,cutoff=500) {
    # Distances and kernel weight
    distances2 <- as.matrix(dist(fordist)) * 69 # in miles
    weights2  <- apply(distances2, 1, function(x) weightfun[[kernel]](dist=x, cut=cutoff))
    E1 <- t(t(e1)) %*% t(e1)
    XeeXhs <- t(X1) %*% (E1 * weights2) %*% X1
  }

  # Demeaned regressors and residuals of the estimation sample
  if (!inherits(reg, "fixest") || is.null(reg$X_demeaned)) {
    stop("reg must be estimated with feols(..., demeaned = TRUE)")
  }
  sel <- obs(reg) # rows of the data used in the estimation
  X <- reg$X_demeaned
  e <- residuals(reg)
  dat <- data.frame(id=idvec[sel], time=timevec[sel], lat=latvec[sel], lon=lonvec[sel])
  olsvcov <- vcov(reg, vcov="iid")
  n <- nrow(X)

  # Correct for spatial correlation:
  timeUnique <- unique(dat$time)

  # Loop over cutoffs
  out <- lapply(dist_cutoff, function(d) {

    print(paste(round(d,1),"miles"))

    # Compute
    XeeXhs <- parallel::mclapply(timeUnique, mc.cores=ncores, FUN=function(t) {
      i <- dat$time==t
      iterateObs(e[i], X[i,,drop=F], as.matrix(dat[i,c("lon","lat")]), cutoff=d)
    })
    XeeX <- Reduce("+",  XeeXhs)

    # Generate VCE for only cross-sectional spatial correlation:
    invXX <- solve(t(X) %*% X) * n
    V <- invXX %*% (XeeX / n) %*% invXX / n
  })

  # Export
  return_list <- c( list(olsvcov), out)
  names(return_list) <- c(-1, dist_cutoff)
  return(return_list)
}

# County centroids
fun$labpt <- function(x) {
  # Returns a 2-column matrix with the centroid of the largest polygon of each
  # feature of an sf object, in the units of its coordinates. This is the
  # point that sp::coordinates() returned for polygons in the 2021 code.
  g <- sf::st_geometry(x)
  out <- t(vapply(seq_along(g), function(i) {
    rings <- lapply(unclass(sf::st_cast(g[i], "POLYGON")), function(p) p[[1]]) # outer rings
    area <- vapply(rings, function(r) {
      n <- nrow(r)
      abs(sum(r[-n,1]*r[-1,2] - r[-1,1]*r[-n,2]))/2
    }, 0)
    r <- rings[[which.max(area)]]
    n <- nrow(r)
    cr <- r[-n,1]*r[-1,2] - r[-1,1]*r[-n,2]
    a <- sum(cr)/2
    c(sum((r[-n,1]+r[-1,1])*cr)/(6*a), sum((r[-n,2]+r[-1,2])*cr)/(6*a))
  }, numeric(2)))
  colnames(out) <- c("x","y")
  out
}

# Raster cells covered by polygons
fun$extract_weights <- function(x, p) {
  # For each polygon of the sf object p, returns a data.frame with the values of
  # the SpatRaster x in the cells covered by the polygon and the approximate
  # fraction of each cell covered ("weight", normalized to sum to 1). This
  # follows raster::extract(x, p, weights=TRUE), used in the 2021 code: each
  # cell is split into 10 x 10 sub-cells (20 x 20 or 100 x 100 when the
  # polygon spans fewer than 17 or 5 cells) and a cell is kept if the polygon
  # contains the center of at least one sub-cell. Polygons too small to contain
  # any sub-cell center get the cells under their vertices.
  v <- terra::vect(p)
  template <- terra::rast(x[[1]])
  half <- max(terra::res(x))/2
  lapply(seq_len(nrow(v)), function(i) {
    e <- as.vector(terra::ext(v[i])) + c(-half, half, -half, half)
    rc <- terra::crop(template, terra::ext(e), snap="near")
    nc <- terra::ncell(rc)
    f <- if (nc < 5) 100 else if (nc < 17) 20 else 10
    cover <- terra::rasterize(v[i], terra::disagg(rc, f), background=0)
    cover <- terra::values(terra::aggregate(cover, f, mean), mat=FALSE)
    cells <- terra::cellFromXY(template, terra::xyFromCell(rc, which(cover > 0)))
    weight <- cover[cover > 0]
    if (length(cells) == 0) { # small polygon
      xy <- sf::st_coordinates(sf::st_cast(sf::st_geometry(p)[i], "POLYGON"))
      xy <- xy[xy[,"L1"]==1, c("X","Y"), drop=FALSE] # outer rings
      cells <- unique(terra::cellFromXY(template, xy))
      weight <- rep(1, length(cells))
    }
    data.frame(terra::extract(x, cells), weight=weight/sum(weight))
  })
}

# Plot a raster with a vertical color legend
fun$plot_raster <- function(x, col, breaks=NULL, legend.width=0.6, legend.shrink=0.5,
                            legend.mar=5.1, legend.args=NULL, axis.args=NULL,
                            maxpixels=5e5, box=TRUE, ...) {
  # Draws a single-layer SpatRaster (terra) the way raster::plot() did in the
  # 2021 code, so the maps keep their original look: rasters with more than
  # maxpixels cells are sampled on a regular grid, and the map region is
  # narrowed to make room for a legend bar on its right.
  # x: single-layer SpatRaster
  # col: vector of colors
  # breaks: break points for the colors (if NULL, colors span the data range)
  # legend.width, legend.shrink, legend.mar: legend size and position (as in raster::plot)
  # legend.args: list of arguments passed to mtext() to label the legend
  # axis.args: list of arguments passed to axis() for the legend
  # ...: other arguments passed to plot() for the map (e.g. axes=F)

  # 1. Values as a matrix, sampled on a regular grid of at most maxpixels cells
  m <- matrix(terra::values(x, mat=FALSE), nrow(x), ncol(x), byrow=TRUE)
  if (length(m) > maxpixels) {
    z <- sqrt(length(m)/maxpixels)
    nr <- max(1, floor(nrow(m)/z))
    nc <- max(1, floor(ncol(m)/z))
    rows <- unique(round(nrow(m)/nr * 1:nr - 0.5*nrow(m)/nr))
    cols <- unique(round(ncol(m)/nc * 1:nc - 0.5*ncol(m)/nc))
    m <- m[rows[rows>0], cols[cols>0], drop=FALSE]
  }
  e <- as.vector(terra::ext(x)) # xmin, xmax, ymin, ymax

  # 2. Colors
  zrange <- range(m, breaks, na.rm=TRUE)
  tocol <- function(z) {
    if (!is.null(breaks)) {
      k <- as.numeric(cut(z, breaks, include.lowest=TRUE))
    } else {
      k <- round((z - zrange[1])/(zrange[2] - zrange[1]) * (length(col) - 1) + 1)
    }
    col[k]
  }
  img <- matrix(tocol(c(m)), nrow(m), ncol(m))

  # 3. Map and legend regions (fractions of the figure region)
  old.par <- par(no.readonly=TRUE)
  char.size <- par("cin")[1]/par("din")[1]
  offset <- char.size * par("mar")[4]
  smallplot <- old.par$plt
  smallplot[2] <- 1 - legend.mar * char.size
  smallplot[1] <- smallplot[2] - legend.width * char.size
  pr <- (smallplot[4] - smallplot[3]) * ((1 - legend.shrink)/2)
  smallplot[3:4] <- smallplot[3:4] + c(pr, -pr)
  bigplot <- old.par$plt
  bigplot[2] <- min(bigplot[2], smallplot[1] - offset)
  smallplot[1:2] <- min(bigplot[2] + offset, smallplot[1]) + c(0, diff(smallplot[1:2]))

  # 4. Map
  par(plt=bigplot)
  lonlat <- isTRUE(terra::is.lonlat(x, perhaps=TRUE, warn=FALSE))
  asp <- ifelse(lonlat, 1/cos(mean(e[3:4]) * pi/180), 1)
  plot(NA, NA, xlim=e[1:2], ylim=e[3:4], type="n", xaxs="i", yaxs="i", asp=asp, xlab="", ylab="", ...)
  rasterImage(as.raster(img), e[1], e[3], e[2], e[4], interpolate=FALSE)
  big.par <- par(no.readonly=TRUE)

  # 5. Legend
  par(new=TRUE, pty="m", plt=smallplot, err=-1)
  plot(NA, NA, xlim=c(0,1), ylim=zrange, type="n", xlab="", ylab="", xaxs="i", yaxs="i", axes=FALSE)
  axis.args <- c(list(side=4, mgp=c(3,1,0), las=2), axis.args)
  if (!is.null(breaks)) {
    bar <- rev(tocol(seq(zrange[1], zrange[2], by=diff(zrange)/100)))
    if (is.null(axis.args$at)) axis.args$at <- breaks
  } else {
    mult <- round(max(1, 100/length(col)))
    z <- ((mult*length(col)):1)/mult
    bar <- col[round((z - min(z))/(max(z) - min(z)) * (length(col) - 1) + 1)]
  }
  rasterImage(as.raster(matrix(bar, ncol=1)), 0, zrange[1], 1, zrange[2], interpolate=FALSE)
  do.call(axis, axis.args)
  box()
  if (!is.null(legend.args)) do.call(mtext, legend.args)

  # 6. Return to the map region so that later calls (e.g. plot(..., add=TRUE)) draw on the map
  mfg <- par("mfg")
  par(big.par)
  par(plt=big.par$plt, xpd=FALSE)
  par(mfg=mfg, new=FALSE)
  if (box) box()
  invisible()
}

# Specification chart function
fun$schart <- function(data, labels=NA, highlight=NA, n=1, index.est=1, index.se=2, index.ci=NA,
                   order="asis", ci=.95, ylim=NA, axes=T, heights=c(1,1), leftmargin=11, offset=c(0,0), ylab="Coefficient", lwd.border=1, horizontal=T,
                   lwd.est=4, pch.est=21, lwd.symbol=2, ref=0, lwd.ref=1, lty.ref=2, col.ref="black", band.ref=NA, col.band.ref=NA,length=0,
                   col.est=c("grey60", "red3"), col.est2=c("grey80","lightcoral"), bg.est=c("white", "white"),
                   col.dot=c("grey60","grey95","grey95","red3"),
                   bg.dot=c("grey60","grey95","grey95","white"),
                   pch.dot=c(22,22,22,22), fonts=c(2,1), adj=c(1,1),cex=c(1,1)) {

  # Authors: Ariel Ortiz-Bobea (ao332@cornell.edu).
  # Version: March 22, 2021
  # If you like this function and use it, please send me a note. It might motivate
  # me improve it or write new ones and share them.

  # Description of arguments

  # Data:
  # data: data.frame with data, ideally with columns 1-2 with coef and SE, then logical variables.
  # labels: list of labels by group. Can also be a character vector if no groups. Default is rownames of data.
  # index.est: numeric indicating position of the coefficient column.
  # index.se: numeric indicating position of the SE column.
  # index.ci: numeric vector indicating position of low-high bars for SE. Can take up to 2 CI, so vector can be up to length 4

  # Arrangement and basic setup:
  # highlight: numeric indicating position(s) of models (row) to highlight in original dataframe.
  # n: size of model grouping. n=1 removes groupings. A vector yields arbitrary groupings.
  # order: whether models should be sorted or not. Options: "asis", "increasing", "decreasing"
  # ci: numeric indicating level(s) of confidence. 2 values can be indicated.
  # ylim: if one wants to set an arbitrary range for Y-axis

  # Figure layout:
  # heights: Ratio of top/bottom panel. Default is c(1,1) for 1 50/50 split
  # (Note: ratio for left/right panel if horizontal=F)
  # leftmargin: amount of space on the left margin
  # offset: vector of numeric with offset for the group and specific labels
  # ylab: Label on the y-axis of top panel. Default is "Coefficient"
  # lwd.border: width of border and other lines
  # horizontal:  should the plot be horizontal? (default is TRUE)

  # Line and symbol styles and colors:
  # lwd.est: numeric indicating the width of lines in the top panel
  # ref: numeric vector indicating horizontal reference lines(s) Default is 0.
  # lty.ref: Style of reference lines. Default is dash line (lty=2).
  # lwd.ref. Width of reference lines. Default is 1.
  # col.ref: vector of colors of reference lines. Default is black.
  # band.ref: vector of 2 numerics indicating upper abdn lower height for a band
  # col.band.ref: color of this band
  # col.est: vector of 2 colors indicating for "other" and "highlighted" models
  # col.est2: same for outer confidence interval if more than 1 confidence interval
  # col.dot: vector of 4 colors indicating colors for borders of symbol in bottom panel for "yes", "no", "NA", and "yes for highlighted model"
  # bg.dot : vector of 4 colors indicating colors for background of symbol in bottom panel for "yes", "no", "NA", and "yes for highlighted model"
  # pch.dot: style of symbols in bottom panel for "yes", "no", "NA", and "yes for highlighted model"
  # length: length of the upper notch on th vertical lines. default is 0.

  # Letter styles
  # fonts: numeric vector indicating font type for group (first) and other labels (second) (e.g. 1:normal, 2:bold, 3:italic)
  # adj: numeric vector indicating alignment adjustment for text label: 0 is left, .5 is center, 1 is right.
  # cex: numeric vector for size of fonts for top panel (first) and bottom panel (Second)

  # 1. Set up
  if (T) {
    # Arrange data
    d <- data
    rownames(d) <- 1:nrow(d)

    # Create ordering vector
    if (order=="asis")       o <- 1:length(d[,index.est])
    if (order=="increasing") o <- order(d[,index.est])
    if (order=="decreasing") o <- order(-d[,index.est])
    if (!is.numeric(d[,index.est])) {warning("index.est does not point to a numeric vector.") ; break}
    d <- d[o,]
    est <- d[,index.est] # Estimate
    if (length(index.ci)>1) {
      l1 <- d[,index.ci[1]]
      h1 <- d[,index.ci[2]]
      if (length(index.ci)>2) {
        l2 <- d[,index.ci[3]]
        h2 <- d[,index.ci[4]]
      }
    } else {
      if (!is.numeric(d[,index.se]))  {warning("index.se does not point to a numeric vector.") ; break}
      se  <- d[,index.se] # Std error
      ci <- sort(ci)
      a <- qnorm(1-(1-ci)/2)
      l1 <- est - a[1]*se
      h1 <- est + a[1]*se
      if (length(ci)>1) {
        l2 <- est - a[2]*se
        h2 <- est + a[2]*se
      }
    }

    # Table
    if (length(index.ci)>1) remove.index <- c(index.est,index.ci) else remove.index <- c(index.est,index.se)
    remove.index <- remove.index[!is.na(remove.index)]
    tab <- t(d[,-remove.index]) # get only the relevant info for bottom panel
    if (!is.list(labels) & !is.character(labels)) labels <- rownames(tab)

    # Double check we have enough labels
    if ( nrow(tab) != length(unlist(labels))) {
      print("Warning: number of labels don't match number of models.")
      labels <- rownames(tab)
    }

    # Plotting objects
    xs <- 1:nrow(d) # the Xs for bars and dots
    if (n[1]>1 & length(n)==1) xs <- xs + ceiling(seq_along(xs)/n) - 1 # group models by n
    if (length(n)>1) {
      if (sum(n) != nrow(d) ) {
        warning("Group sizes don't add up.")
      } else {
        idx <- unlist(lapply(1:length(n), function(i) rep(i,n[i])))
        xs <- xs + idx - 1
      }
    }
    h <- nrow(tab) + ifelse(is.list(labels),length(labels),0) # number of rows in table
    # Location of data and labels
    if (is.list(labels)) {
      index <- unlist(lapply(1:length(labels), function(i) rep(i, length(labels[[i]])) ))
      locs <- split(1:length(index),index)
      locs <- lapply(unique(index), function(i) {
        x <- locs[[i]]+i-1
        x <- c(x,max(x)+1)
      })
      yloc  <- unlist(lapply(locs, function(i) i[-1])) # rows where data points are located
      yloc2 <- sapply(locs, function(i) i[1]) # rows where group lables are located
    } else {
      yloc <- 1:length(labels)
    }

    # Range
    if (is.na(ylim[1]) | length(ylim)!=2) {
      if (length(index.ci)>2 | length(ci)>1) {
        ylim <- range(c(l2,h2,ref)) # range that includes reference lines
      } else {
        ylim <- range(c(l1,h1,ref))
      }
      ylim <- ylim + diff(ylim)/10*c(-1,1) # and a bit more
    }
    xlim <- range(xs) #+ c(1,-1)
  }

  # 2. Plot
  if (T) {
    #par(mfrow=c(2,1), mar=c(0,leftmargin,0,0), oma=oma, xpd=F, family=family)
    if (horizontal) {
      layout(t(t(2:1)), height=heights, widths=1)
      par(mar=c(0,leftmargin,0,0), xpd=F)
    } else  {
      layout(t(1:2), height=1, widths=heights)
      par(mar=c(leftmargin,0,0,0), xpd=F)
    }

    # Bottom panel (plotted first)
    if (horizontal)  plot(1:nrow(tab), xlab="", ylab="", axes=F, type="n", ylim=c(h,1), xlim=xlim)
    if (!horizontal) plot(1:nrow(tab), xlab="", ylab="", axes=F, type="n", ylim=xlim, xlim=c(1,h))
    lapply(1:nrow(tab), function(i) {
      # Get colors and point type
      type <- ifelse(is.na(tab[i,]),3,ifelse(tab[i,]==TRUE,1, ifelse(tab[i,]==FALSE,2,NA)))
      type <- ifelse(names(type) %in% paste(highlight) & type==1,4,type) # replace colors for baseline model
      col <- col.dot[type]
      bg  <- bg.dot[type]
      pch <- as.numeric(pch.dot[type])
      sel <- is.na(pch)
      # Plot points
      if (horizontal) {
        points(xs, rep(yloc[i],length(xs)), col=col, bg=bg, pch=pch, lwd=lwd.symbol)
        points(xs[sel], rep(yloc[i],length(xs))[sel], col=col[sel], bg=bg[sel], pch=pch.dot[3]) # symbol for missing value
      } else {
        points(rep((yloc)[i],length(xs)),xs, col=col, bg=bg, pch=pch, lwd=lwd.symbol)
        points(rep((yloc)[i],length(xs))[sel],xs[sel], col=col[sel], bg=bg[sel], pch=pch.dot[3]) # symbol for missing value
      }

    })
    par(xpd=T)
    if (is.list(labels)) {
      if (horizontal) {
        text(-offset[1], yloc2, labels=names(labels), adj=adj[1], font=fonts[1], cex=cex[2])
      } else {
        text((yloc2), -offset[1], labels=names(labels), adj=adj[1], font=fonts[1], cex=cex[2], srt=90)
      }
    }
    # Does not accomodate subscripts
    if (horizontal) {
      text(-rev(offset)[1], yloc , labels=unlist(labels), adj=rev(adj)[1], font=fonts[2], cex=cex[2])
    } else {
      text(yloc,-rev(offset)[1], labels=unlist(labels), adj=rev(adj)[1], font=fonts[2], cex=cex[2], srt=90)
    }
    # Accomodates subscripts at the end of each string
    if (F) {
      labels1 <- unlist(labels)
      lapply(1:length(labels1), function(i) {
        a  <- labels1[i]
        a1 <- strsplit(a,"\\[|\\]")[[1]][1]
        a2 <- rev(strsplit(a,"\\[|\\]")[[1]])[1]
        if (identical(a1,a2))  a2 <- NULL
        text(-rev(offset)[1], yloc[i], labels=bquote(.(a1)[.(a2)]), adj=adj[2], font=fonts[2], cex=cex[2])
      })
    }
    par(xpd=F)

    # Top panel (plotted second)
    colvec  <- ifelse(colnames(tab) %in% paste(highlight), col.est[2], col.est[1])
    bg.colvec  <- ifelse(colnames(tab) %in% paste(highlight), bg.est[2], bg.est[1])
    colvec2 <- ifelse(colnames(tab) %in% paste(highlight),col.est2[2], col.est2[1])
    if (horizontal)   plot(est, xlab="", ylab="", axes=F, type="n", ylim=ylim, xlim=xlim)
    if (!horizontal)  plot(est, xlab="", ylab="", axes=F, type="n", ylim=xlim, xlim=ylim)
    # Band if present
    if (!is.na(band.ref[1])) {
      if (horizontal) {
        rect(xleft=min(xlim)-diff(xlim)/10, ybottom=band.ref[1], xright=max(xlim)+diff(xlim)/10, ytop=band.ref[2],
             col=col.band.ref, border=NA)
      } else {
        rect(ybottom=min(xlim)-diff(xlim)/10, xleft=band.ref[1], ytop=max(xlim)+diff(xlim)/10, xright=band.ref[2],
             col=col.band.ref, border=NA)
      }
    }
    # Reference lines
    if (horizontal) {
      abline(h=ref, lty=lty.ref, lwd=lwd.ref, col=col.ref)
    } else {
      abline(v=ref, lty=lty.ref, lwd=lwd.ref, col=col.ref)
    }
    # Vertical bars
    if (horizontal) {
      if (length(ci)>1 | length(index.ci)>2) arrows(x0=xs, y0=l2, x1=xs, y1=h2, length=length, code=3, lwd=rev(lwd.est)[1], col=colvec2, angle=90)
      arrows(x0=xs, y0=l1, x1=xs, y1=h1, length=length, code=3, lwd=lwd.est[1]     , col=colvec, angle=90)
      points(xs, est, pch=pch.est, lwd=lwd.symbol, col=colvec, bg=bg.colvec)
    } else {
      if (length(ci)>1 | length(index.ci)>2) arrows(y0=xs, x0=l2, y1=xs, x1=h2, length=length, code=3, lwd=rev(lwd.est)[1], col=colvec2, angle=90)
      arrows(y0=xs, x0=l1, y1=xs, x1=h1, length=length, code=3, lwd=lwd.est[1]     , col=colvec, angle=90)
      points(est,xs, pch=pch.est, lwd=lwd.symbol, col=colvec, bg=bg.colvec)
    }
    # Axes
    if (axes) {
      axis(ifelse(horizontal,2,1), las=2, cex.axis=cex[1], lwd=lwd.border)
      axis(ifelse(horizontal,4,3), labels=NA, lwd=lwd.border)
    }
    mtext(ylab, side=ifelse(horizontal,2,1), line=3.5, cex=cex[1])
    box(lwd=lwd.border)

  }

}
