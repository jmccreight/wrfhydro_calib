.libPaths("/glade/u/home/adugger/system/R/Libraries/R3.2.2")
library(rwrfhydro)
library(data.table)
library(ggplot2)
library(plyr)

#########################################################
# SETUP
#########################################################

source("calib_utils.R")

# Metrics
metrics <- c("obj", "nse", "nselog", "cor", "rmse", "bias", "kge", "fdc")


#########################################################
# MAIN CODE
#########################################################

# First run so need to initialize
ReadNamelist("namelist.calib")
cyclecount <- 0

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
      x_all[cnt,] <- c(as.numeric(x0), "", "")
      x_all[cnt,ivar] <- as.numeric(paramBnds$min[paramBnds$param==ivar])
      x_all[cnt,"tag"] <- ivar
      x_all[cnt,"bound"] <- "Low"
      x_all[cnt+1,] <- c(x0, "", "")
      x_all[cnt+1,ivar] <- as.numeric(paramBnds$max[paramBnds$param==ivar])
      x_all[cnt+1,"tag"] <- ivar
      x_all[cnt+1,"bound"] <- "High"
      cnt <- cnt+2
}
for (i in xnames) { x_all[,i] <- as.numeric(x_all[,i]) }

# Initialize parameter archive DF
message("Initialize parameter archive")
x_archive <- as.data.frame(matrix(, nrow=1, ncol=length(xnames)+11))
names(x_archive) <- c("id", xnames, "tag", "bound", metrics)

# Output parameter set
message("Output parameter set")
x_all_out <- cbind(data.frame(id=seq(1:nrow(x_all))), x_all[,1:14])
write.table(x_all_out, file="params_new.txt", row.names=FALSE, sep=" ")

# Save and exit
save.image("proj_data_SENS.Rdata")
quit("no")

