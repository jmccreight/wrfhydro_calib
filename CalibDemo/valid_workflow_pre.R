.libPaths("/glade/u/home/adugger/system/R/Libraries/R3.2.2")
library(rwrfhydro)
library(data.table)
library(ggplot2)
library(plyr)

#########################################################
# SETUP
#########################################################


#########################################################
# MAIN CODE
#########################################################

if (file.exists("proj_data.Rdata")) {
      load("proj_data.Rdata")
      sel <- x_archive
      for (j in names(x_best)) {
         sel <- subset(sel, round(sel[,j], 4) == round(x_best[j], 4))
      }
      x_new_out <- c(sel$id, x_best)
      names(x_new_out)[1] <- "id"
      write.table(data.frame(t(x_new_out)), file="params_best.txt", row.names=FALSE, sep=" ")
} else {
      message("No proj_data.Rdata file found.")
}

quit("no")

