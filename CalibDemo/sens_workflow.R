.libPaths("/glade/u/home/adugger/system/R/Libraries/R3.2.2")
library(rwrfhydro)
library(data.table)
library(ggplot2)
library(plyr)

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

# Metrics
metrics <- c("obj", "nse", "nselog", "cor", "rmse", "bias", "kge", "fdc")


#########################################################
# MAIN CODE
#########################################################

# First loop check
if (file.exists("proj_data_SENS.Rdata")) { 
   load("proj_data_SENS.Rdata")
} else {
   # First run so need to initialize
   ReadNamelist("namelist.calib")
   cyclecount <- 0

   # Load obs so we have them for next iteration
   load(obsFile)
   obsDT$q_cms <- NULL

   # Find the index of the gage
   rtLink <- ReadRouteLink(rtlinkFile)
   rtLink <- data.table(rtLink)
   linkId <- which(trimws(rtLink$gages) %in% siteId)

   # Setup value lists from paramBnds
   message("Setup lists")
   xnames <- paramBnds$param
   x0 <- paramBnds$ini
   names(x0) <- xnames
   x_min <- paramBnds$min
   names(x_min) <- xnames
   x_max <- paramBnds$max
   names(x_max) <- xnames

   # Setup parameter sets
   message("Setup parameter sets")
   x_all <- data.frame(t(x0))
   x_all$tag <- "control"
   x_all$bound <- "control"
   cnt <- 2
   for (ivar in names(x0)) {
      message(ivar)
      x_all[cnt,] <- c(x0, "", "")
      x_all[cnt,ivar] <- paramBnds$min[paramBnds$param==ivar]
      x_all[cnt,"tag"] <- ivar
      x_all[cnt,"bound"] <- "Low"
      x_all[cnt+1,] <- c(x0, "", "")
      x_all[cnt+1,ivar] <- paramBnds$max[paramBnds$param==ivar]
      x_all[cnt+1,"tag"] <- ivar
      x_all[cnt+1,"bound"] <- "High"
      cnt <- cnt+2
   }

   # Initialize chrtout
   chrt.d.all <- data.table()

   # Initialize parameter archive DF
   message("Initialize parameter archive")
   x_archive <- as.data.frame(matrix(, nrow=1, ncol=length(xnames)+11))
   names(x_archive) <- c("id", xnames, "tag", "bound", metrics)

   # Output parameter set
   message("Output parameter set")
   x_new <- x0
   cyclecount <- 1

   x_new_out <- c(cyclecount, x_new)
   names(x_new_out)[1] <- "id"
   write.table(data.frame(t(x_new_out)), file="params_new.txt", row.names=FALSE, sep=" ")

   # Save and exit
   save.image("proj_data_SENS.Rdata")
   quit("no")
}

if (cyclecount > 0) {

   file.rename("params_new.txt", paste0("SENS_RESULTS/OUTPUT", cyclecount, "/params_new.txt"))

   # Read model out and calculate performance metric
   outPath <- paste0(runDir, "/SENS_RESULTS/OUTPUT", cyclecount)
   print(outPath)

   # Read files
   message("Reading model out files.")
   system.time({
   filesList <- list.files(path = outPath,
                          pattern = glob2rx("*.CHRTOUT_DOMAIN*"),
                          full.names = TRUE)
   filesListDate <- as.POSIXct(unlist(plyr::llply(strsplit(basename(filesList),"[.]"), '[',1)), format = "%Y%m%d%H%M", tz = "UTC")
   whFiles <- which(filesListDate >= startDate)
   filesList <- filesList[whFiles]
   if (length(filesList) == 0) stop("No matching files in specified directory.")
   chrt <- as.data.table(plyr::ldply(filesList, ReadChFile, linkId, .parallel = parallelFlag))
   })

   # Convert to daily
   chrt.d <- Convert2Daily(chrt)
   chrt.d[, site_no := siteId]
   # Merge
   setkey(chrt.d, "site_no", "POSIXct")
   setkey(obsDT, "site_no", "POSIXct")
   chrt.d <- merge(chrt.d, obsDT, by=c("site_no", "POSIXct"), all.x=FALSE, all.y=FALSE)
   chrt.d$id <- cyclecount
   chrt.d$tag <- x_all$tag[cyclecount]
   chrt.d.all <- rbindlist(list(chrt.d.all, chrt.d))

   # Assess performance
   F_new <- objFn(chrt.d$q_cms, chrt.d$obs)
   statNse <- rwrfhydro::Nse(chrt.d$q_cms, chrt.d$obs)
   statNseLog <- rwrfhydro::NseLog(chrt.d$q_cms, chrt.d$obs)
   statCor <- cor(chrt.d$q_cms, chrt.d$obs)
   statRmse <- rwrfhydro::Rmse(chrt.d$q_cms, chrt.d$obs)
   statBias <- sum(chrt.d$q_cms - chrt.d$obs, na.rm=TRUE)/sum(chrt.d$obs, na.rm=TRUE) * 100
   statKge <- Kge(chrt.d$q_cms, chrt.d$obs)
   chrt.d <- CalcFdc(chrt.d, "q_cms")
   chrt.d <- CalcFdc(chrt.d, "obs")
   statFdc <- integrate(splinefun(as.data.frame(chrt.d)[,"q_cms.fdc"], as.data.frame(chrt.d)[,"q_cms"], method='natural'), 0.05, 0.95, subdivisions=1000)$value

   # Archive results
   x_archive[cyclecount,] <- c(cyclecount, x_new, x_all$tag[cyclecount], x_all$bound[cyclecount], F_new, statNse, statNseLog, statCor, statRmse, statBias, statKge, statFdc)
   print(x_archive[cyclecount,])

   if (cyclecount < nrow(x_all)) {
      # Select next parameter set
      cyclecount <- cyclecount+1
      x_new <- as.numeric(x_all[cyclecount,1:length(xnames)])
      names(x_new) <- xnames

      # Output next parameter set
      x_new_out <- c(cyclecount, x_new)
      names(x_new_out)[1] <- "id"
      write.table(data.frame(t(x_new_out)), file="params_new.txt", row.names=FALSE, sep=" ")

   } else {

      # Summary stats
      for (i in names(x_archive)[18:25]) { x_archive[,i] <- as.numeric(x_archive[,i]) }
      x_archive$obj.diff <- x_archive$obj - x_archive$obj[1]
      x_archive$nse.diff <- x_archive$nse - x_archive$nse[1]
      x_archive$nselog.diff <- x_archive$nselog - x_archive$nselog[1]
      x_archive$cor.diff <- x_archive$cor - x_archive$cor[1]
      x_archive$rmse.diff <- x_archive$rmse - x_archive$rmse[1]
      x_archive$bias.diff <- x_archive$bias - x_archive$bias[1]
      x_archive$kge.diff <- x_archive$kge - x_archive$kge[1]
      x_archive$fdc.diff <- x_archive$fdc - x_archive$fdc[1]

      # Summary plots

      # Plot setup
      ggPalette <- gg_color_hue(14)

      plotGroups <- list(soil=c('bexp', 'dksat', 'smcmax', 'refkdt', 'slope', 'RETDEPRTFAC', 'LKSATFAC'),
                   other=c('Zmax', 'Expon', 'CWPVT', 'VCMX25', 'MP', 'HVT', 'MFSNO'))

      writePlotDir <- paste0(runDir, "/plots")
      dir.create(writePlotDir)

      # Hydrographs
      gg <- ggplot(data=chrt.d.all[tag %in% plotGroups[["soil"]],], aes(x=POSIXct, y=q_cms, color=tag, group=id)) +
              geom_line(lwd=0.6) +
              scale_y_log10() +
              facet_wrap(~ tag) +
              annotate(geom='line', x=chrt.d.all[tag=="control",]$POSIXct, y=chrt.d.all[tag=="control",]$q_cms, color='black', lwd=0.3) +
              ggtitle(paste0("Model Sensitivity: ", chrt.d.all$site_no[1]))
      ggsave(filename=paste0(writePlotDir, "/", chrt.d.all$site_no[1], "_soilp_sens.png"),
              plot=gg, units="in", width=16, height=8, dpi=300)


      gg <- ggplot(data=chrt.d.all[tag %in% plotGroups[["other"]],], aes(x=POSIXct, y=q_cms, color=tag, group=id)) +
              geom_line(lwd=0.6) +
              scale_y_log10() +
              facet_wrap(~ tag) +
              annotate(geom='line', x=chrt.d.all[tag=="control",]$POSIXct, y=chrt.d.all[tag=="control",]$q_cms, color='black', lwd=0.3) +
              ggtitle(paste0("Model Sensitivity: ", chrt.d.all$site_no[1]))
      ggsave(filename=paste0(writePlotDir, "/", chrt.d.all$site_no[1], "_otherp_sens.png"),
              plot=gg, units="in", width=16, height=8, dpi=300)

     # Metric plots
     tmp<-x_archive[,c("tag", metrics)]
     tmp<-melt(tmp, id='tag')
     tmpmax <- ddply(tmp, .(tag, variable), function(x) {max(x$value)})
     tmpmin <- ddply(tmp, .(tag, variable), function(x) {min(x$value)})
     tmpall <- merge(tmpmin, tmpmax, by=c("tag", "variable"))
     names(tmpall)[3:4]<-c("min", "max")
     tmpall <- plyr::join(tmpall, data.frame(tag=unique(tmpall$tag), id=seq(1:length(unique(tmpall$tag)))), by="tag")

     gg <- ggplot(data=tmpall, aes(x=tag, fill=tag)) + 
        geom_rect(aes(x=tag, xmin=id-0.45, xmax=id+0.45, ymin=min, ymax=max, color=tag), alpha=0.8) + 
        facet_wrap(~variable, scales="free_y") +
        ggtitle(paste0("Metric Sensitivity: ", chrt.d.all$site_no[1])) +
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
     ggsave(filename=paste0(writePlotDir, "/", chrt.d.all$site_no[1], "_metric_sens.png"),
              plot=gg, units="in", width=16, height=8, dpi=300)

   }

   # Save and exit
   if (parallelFlag) stopCluster(cl)
   save.image("proj_data_SENS.Rdata")
   quit("no")

}

