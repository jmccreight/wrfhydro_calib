# Namelist read function

ReadNamelist <- function(nlist) {
   source(nlist)
   load(obsFile)
}

# Convert to daily flow

Convert2Daily <- function(str) {
   str$Date <- rwrfhydro::CalcDateTrunc(str$POSIXct)
   setkey(str, Date)
   str.d <- str[, list(q_cms=mean(q_cms, na.rm=TRUE)), by = "Date"]
   str.d$POSIXct <- as.POSIXct(paste0(str.d$Date, " 00:00"), tz="UTC")
   str.d
}

# Read streamflow from netcdf file

ReadChFile <- function(file, idList){
    nc <- ncdf4::nc_open(file)
    output <- data.frame(q_cms = ncdf4::ncvar_get(nc, varid = "streamflow", start = idList , count =1),
                         POSIXct = as.POSIXct(strsplit(basename(file),"[.]")[[1]][1], format = "%Y%m%d%H%M", tz = "UTC"))
    ncdf4::nc_close(nc)
    return(output)
}


# DDS parameter selection function

DDS.sel <- function(i, m, r, xnames, x_min, x_max, x_best) {

   # Set parameter set
   P_i <- 1-log(i)/log(m)
   sel <- c()
   for (d in 1:length(xnames)) {
      sel[d] <- sample(c(1,0), 1, prob=c(P_i, 1-P_i))
   }
   N <- xnames[as.logical(sel)]
   if (length(N) < 1) N <- sample(xnames, 1)

   # Set new values for selected parameters
   x_new <- x_best
   for (j in N) {
      xj_min <- x_min[[j]]
      xj_max <- x_max[[j]]
      xj_best <- x_best[[j]]
      sigj <- r * (xj_max - xj_min)
      x_new[[j]] <- xj_best + sigj*rnorm(1)
      if (x_new[[j]] < xj_min) {
         x_new[[j]] <- xj_min + (xj_min - x_new[[j]])
         if (x_new[[j]] > xj_max) {
            x_new[[j]] <- xj_min
         }
      }
      if (x_new[[j]] > xj_max) {
         x_new[[j]] <- xj_max - (x_new[[j]] - xj_max)
         if (x_new[[j]] < xj_min) {
            x_new[[j]] <- xj_max
         }
      }
   }
   x_new

}

# Kge
Kge <- function (m, o, na.rm=TRUE, s.r=1, s.alpha=1, s.beta=1) {
  use <- if(na.rm) 'pairwise.complete.obs' else 'everything'
  r     <- cor(m, o, use=use)
  alpha <- sd(m, na.rm=na.rm) / sd(o, na.rm=na.rm)
  beta  <- mean(m, na.rm=na.rm) / mean(o, na.rm=na.rm)
  kge = sqrt( (s.r*(1-r))^2 + (s.alpha*(1-alpha))^2 + (s.beta*(1-beta))^2 )
  kge
}

# Flow duration curve calcs
CalcFdc <- function(strDf, strCol="q_cms") {
    if (data.table::is.data.table(strDf)) {
        tmp <- rank(-strDf[, strCol, with=FALSE], na.last="keep")
        strDf[, paste0(strCol,".fdc") := NA]
        strDf[, paste0(strCol,".fdc") := tmp/(sum(!is.na(strDf[,strCol,with=FALSE]))+1)]
    } else {
        tmp <- rank(-strDf[,strCol],na.last="keep")
        strDf[,paste0(strCol,".fdc")] <- NA
        strDf[,paste0(strCol,".fdc")] <- tmp/(sum(!is.na(strDf[,strCol]))+1)
    }
    strDf
}

# ggplot color palette
gg_color_hue <- function(n) {
   hues = seq(15, 375, length = n + 1)
   hcl(h = hues, l = 65, c = 100)[1:n]
 }

# Multiplot
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
