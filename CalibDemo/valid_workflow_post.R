.libPaths("/glade/u/home/adugger/system/R/Libraries/R3.2.2")
library(rwrfhydro)
library(data.table)
library(ggplot2)
library(plyr)
library(gridExtra)

#########################################################
# SETUP
#########################################################

source("calib_utils.R")

# Multi-core
parallelFlag <- TRUE
ncores <- 8
if (parallelFlag && ncores>1) {
        library(doParallel)
        cl <- makeForkCluster(ncores)
        registerDoParallel(cl)
}

# Calib/valid date split
validDate <- as.POSIXct("2013-10-01", format="%Y-%m-%d", tz="UTC")

# Initialize results
metrics <- c("obj", "nse", "nselog", "cor", "rmse", "bias", "kge")
results <- as.data.frame(matrix(, nrow=4, ncol=length(metrics)+3))
names(results) <- c("site_no", "run", "period", metrics)


#########################################################
# MAIN CODE
#########################################################

ReadNamelist("namelist.calib")

# Load obs
load(obsFile)
obsDT$q_cms <- NULL
obsDT <- obsDT[!is.na(obs),]

# Find the index of the gage
rtLink <- ReadRouteLink(rtlinkFile)
rtLink <- data.table(rtLink)
linkId <- which(trimws(rtLink$gages) %in% siteId)

   ### VALIDATION

   # Read model out and calculate performance metric
   outPath <- paste0(runDir, "/RUN.VALID/OUTPUT/")
   print(outPath)

   # Read files
   message("Reading validation model out files.")
   system.time({
   filesList <- list.files(path = outPath,
                          pattern = glob2rx("*.CHRTOUT_DOMAIN*"),
                          full.names = TRUE)
   filesListDate <- as.POSIXct(unlist(plyr::llply(strsplit(basename(filesList),"[.]"), '[',1)), format = "%Y%m%d%H%M", tz = "UTC")
   whFiles <- which(filesListDate >= startDate)
   filesList <- filesList[whFiles]
   if (length(filesList) == 0) stop("No matching files in specified directory.")
   chrt.valid <- as.data.table(plyr::ldply(filesList, ReadChFile, linkId, .parallel = parallelFlag))
   })

   # Convert to daily
   chrt.valid.d <- Convert2Daily(chrt.valid)
   chrt.valid.d[, site_no := siteId]

   save.image("proj_data_VALID.Rdata")

   ### CONTROL

   # Sead model out and calculate performance metric
   outPath <- paste0(runDir, "/RUN.CONTROL/OUTPUT/")
   print(outPath)

   # Read files
   message("Reading control model out files.")
   system.time({
   filesList <- list.files(path = outPath,
                          pattern = glob2rx("*.CHRTOUT_DOMAIN*"),
                          full.names = TRUE)
   filesListDate <- as.POSIXct(unlist(plyr::llply(strsplit(basename(filesList),"[.]"), '[',1)), format = "%Y%m%d%H%M", tz = "UTC")
   whFiles <- which(filesListDate >= startDate)
   filesList <- filesList[whFiles]
   if (length(filesList) == 0) stop("No matching files in specified directory.")
   chrt.cont <- as.data.table(plyr::ldply(filesList, ReadChFile, linkId, .parallel = parallelFlag))
   })

   # Convert to daily
   chrt.cont.d <- Convert2Daily(chrt.cont)
   chrt.cont.d[, site_no := siteId]

   if (parallelFlag) stopCluster(cl)

   save.image("proj_data_VALID.Rdata")

   ### MERGE

   # Merge
   setkey(chrt.valid.d, "site_no", "POSIXct")
   setkey(chrt.cont.d, "site_no", "POSIXct")
   setkey(obsDT, "site_no", "POSIXct")
   chrt.valid.d <- merge(chrt.valid.d, obsDT, by=c("site_no", "POSIXct"), all.x=FALSE, all.y=FALSE)
   chrt.valid.d[, run := "validation"]
   chrt.cont.d <- merge(chrt.cont.d, obsDT, by=c("site_no", "POSIXct"), all.x=FALSE, all.y=FALSE)
   chrt.cont.d[, run := "control"]
   tmpobs <- chrt.cont.d[,c("POSIXct", "obs", "run"), with=FALSE]
   names(tmpobs)[2] <- "q_cms"
   tmpobs[, run := "observed"]
   chrt.all.d <- rbindlist(list(chrt.cont.d[,c("POSIXct", "q_cms", "obs", "run"), with=FALSE], 
                              chrt.valid.d[,c("POSIXct", "q_cms", "obs", "run"), with=FALSE],
                              tmpobs), fill=TRUE)

   ### METRICS

   CalcStats <- function(chrt.d) {
      F_new <- objFn(chrt.d$q_cms, chrt.d$obs)
      statNse <- rwrfhydro::Nse(chrt.d$q_cms, chrt.d$obs)
      statNseLog <- rwrfhydro::NseLog(chrt.d$q_cms, chrt.d$obs)
      statCor <- cor(chrt.d$q_cms, chrt.d$obs)
      statRmse <- rwrfhydro::Rmse(chrt.d$q_cms, chrt.d$obs)
      statBias <- sum(chrt.d$q_cms - chrt.d$obs, na.rm=TRUE)/sum(chrt.d$obs, na.rm=TRUE) * 100
      statKge <- Kge(chrt.d$q_cms, chrt.d$obs)
      c(F_new, statNse, statNseLog, statCor, statRmse, statBias, statKge)
   }

   # Control
   results[1,] <- c(siteId, "control", "full", CalcStats(chrt.cont.d))
   results[2,] <- c(siteId, "control", "calib", CalcStats(chrt.cont.d[POSIXct < validDate,]))
   results[3,] <- c(siteId, "control", "valid", CalcStats(chrt.cont.d[POSIXct >= validDate,]))

   # Validation
   results[4,] <- c(siteId, "valid", "full", CalcStats(chrt.valid.d))
   results[5,] <- c(siteId, "valid", "calib", CalcStats(chrt.valid.d[POSIXct < validDate,]))
   results[6,] <- c(siteId, "valid", "valid", CalcStats(chrt.valid.d[POSIXct >= validDate,]))

   ### PLOTS

   writePlotDir <- paste0(runDir, "/plots")
   dir.create(writePlotDir)

   # Hydrographs
   gg <- ggplot(data=chrt.all.d[run != "observed",], aes(x=POSIXct, y=q_cms, color=run)) +
              geom_line(lwd=0.6) +
              geom_line(data=chrt.all.d[run == "observed",], aes(x=POSIXct, y=q_cms, color='observed'), lwd=0.4) +
              geom_vline(xintercept=as.numeric(validDate), lwd=1.8, col=alpha('grey70', 0.7), lty=2) +
              ggtitle(paste0("Model Validation Hydrograph: ", siteId)) +
              scale_color_manual(name="", values=c('dodgerblue', 'orange', 'black'), 
                                limits=c('control','validation','observed'),
                                label=c('control', 'validation', 'observed')) +
              labs(x="", y="Streamflow (m3/s)") +
              theme_bw()

   ggsave(filename=paste0(writePlotDir, "/", siteId, "_valid_hydrogr.png"),
              plot=gg, units="in", width=16, height=8, dpi=300)
   ggsave(filename=paste0(writePlotDir, "/", siteId, "_valid_hydrogr_log.png"),
              plot=gg+scale_y_log10(), units="in", width=16, height=8, dpi=300)

   # Scatterplots
   maxval <- max(chrt.all.d$q_cms)
   gg1 <- ggplot(data=chrt.all.d[run != "observed",], aes(x=obs, y=q_cms, color=run)) +
              geom_point(shape=1, size=3) +
              scale_shape_discrete(solid=FALSE) +
              geom_abline(intercept=0, slope=1, col='black', lty=1) +
              ggtitle(paste0("Full Period (WY2009-2016): ", siteId)) +
              scale_color_manual(name="", values=c('dodgerblue', 'orange'),
                                limits=c('control','validation'),
                                label=c('control', 'validation')) +
              labs(x="", y="modeled streamflow (m3/s)") +
              theme_bw() + theme(legend.position="none") +
              xlim(0,maxval) + ylim(0,maxval)
   gg2 <- ggplot(data=chrt.all.d[run != "observed" & POSIXct < validDate,], aes(x=obs, y=q_cms, color=run)) +
              geom_point(shape=1, size=3) +
              scale_shape_discrete(solid=FALSE) +
              geom_abline(intercept=0, slope=1, col='black', lty=1) +
              ggtitle(paste0("Calibration Period (WY2009-2013): ", siteId)) +
              scale_color_manual(name="", values=c('dodgerblue', 'orange'),
                                limits=c('control','validation'),
                                label=c('control', 'validation')) +
              labs(x="observed streamflow (m3/s)", y="") +
              theme_bw() + theme(legend.position="none") +
              xlim(0,maxval) + ylim(0,maxval)
   gg3 <- ggplot(data=chrt.all.d[run != "observed" & POSIXct >= validDate,], aes(x=obs, y=q_cms, color=run)) +
              geom_point(shape=1, size=3) +
              scale_shape_discrete(solid=FALSE) +
              geom_abline(intercept=0, slope=1, col='black', lty=1) +
              ggtitle(paste0("Validation Period (WY2014-2016): ", siteId)) +
              scale_color_manual(name="", values=c('dodgerblue', 'orange'),
                                limits=c('control','validation'),
                                label=c('control', 'validation')) +
              labs(x="", y="") +
              theme_bw() +
              xlim(0,maxval) + ylim(0,maxval)
   gg.all <- grid.arrange(gg1, gg2, gg3, ncol=3)

   ggsave(filename=paste0(writePlotDir, "/", siteId, "_valid_scatter.png"),
              plot=gg.all, units="in", width=16, height=8, dpi=300)


   # Stats Barplots
   results.plot <- melt(results[,c("run", "period", metrics)], id=c("period", "run"))
   results.plot$period <- factor(results.plot$period, levels=c("calib", "valid", "full"))
   results.plot <- results.plot[order(results.plot$variable, results.plot$period, results.plot$run),]
   results.plot$value <- as.numeric(results.plot$value)
   gg <- ggplot(data=results.plot, aes(x=factor(period), y=value, fill=factor(run))) + 
         geom_bar(stat="identity", position="dodge") + 
         facet_wrap(~variable, scales="free_y") +
         scale_fill_manual(name="", values=c('dodgerblue', 'orange'), 
             limits=c('control','valid'),
             label=c('control', 'validation')) +
         ggtitle(paste0("Model Validation Performance Metrics: ", siteId)) +
         labs(x="run period", y="value") +
         theme_bw()
   ggsave(filename=paste0(writePlotDir, "/", siteId, "_valid_metrics.png"),
        plot=gg, units="in", width=16, height=8, dpi=300)

   # Save and exit
   save.image("proj_data_VALID.Rdata")
   quit("no")



