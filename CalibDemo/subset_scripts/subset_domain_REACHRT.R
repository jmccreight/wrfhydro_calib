################ INPUTS #####################

### CHANGE THESE VALUES ###

# Specify your new domain file directory
myPath <- "/glade/scratch/adugger/CONUS_CALIB/CalibDemo/RUN.TEMPLATE/DOMAIN/"

# Specify the clip bounding coordinates
        y_south <- 611800
        y_north <- 696570
        x_west <- -1098460
        x_east <- -1045560

# Projection for boundaing coordinates. This needs to be a PROJ4 string 
# (e.g., "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs").
        coordProj <- "+proj=lcc +lat_1=30 +lat_2=60 +lat_0=40.0000076293945 +lon_0=-97 +x_0=0 +y_0=0 +a=6370000 +b=6370000 +units=m +no_defs"

# Multiplier between routing grid and LSM grid
# (e.g., 1-km LSM and 250-m routing means a value of 4)
        dxy <- 4
# Number of cells to buffer
        cellBuff <- 2


# Specify the ORIGINAL (full extent) domain files:

	# Routing domain file
	fullHydFile <- "/glade/p/ral/RHAP/adugger/WH_REPO/DOMAINS/CONUS_IOC/Fulldom_hires_netcdf_file_250m_CALIB3_lksatfac.nc"

	# Geogrid domain file
	fullGeoFile <- "/glade/p/ral/RHAP/adugger/WH_REPO/DOMAINS/CONUS_IOC/geo_em.d01.nc.conus_1km_nlcd11_glacfix"

	# Wrfinput domain file
	fullWrfFile <- "/glade/p/ral/RHAP/adugger/WH_REPO/DOMAINS/CONUS_IOC/wrfinput_d01_1km_nlcd11_glacfix_soilwatfix"

	# Route link file
	fullRtlinkFile <- "/glade/p/ral/RHAP/adugger/WH_REPO/DOMAINS/CONUS_IOC/RouteLink_2016_04_07.nudgingOperational2016-04-08_chanparm3_mann_BtmWdth_XYCOORD.nc"

	# Spatial weights file
	fullSpwtFile <- "/glade/p/ral/RHAP/adugger/WH_REPO/DOMAINS/CONUS_IOC/spatialweights_IOC_all_basins_250m_2015_12_30.nc"

	# GW bucket parameter file
	fullGwbuckFile <- "/glade/p/ral/RHAP/adugger/WH_REPO/DOMAINS/CONUS_IOC/GWBUCKPARM_OCONUS_2016_01_03_CALIB3.nc"

	# Soil parameter file
	fullSoilparmFile <- "/glade/p/ral/RHAP/adugger/WH_REPO/DOMAINS/CONUS_IOC/soil_properties_CONUS_v4_ter_rtg_CALIB3_nlcd11_SandDwsatFix.nc"

	# Lake parameter file
	fullLakeparmFile <- "/glade/p/ral/RHAP/adugger/WH_REPO/DOMAINS/CONUS_IOC/LAKEPARM_ks_scaled_max_el_Feb_12_2016_1d_WeirC_OrificeC_min_alt_1.0_WeirL_35.0_RELABEL_goodlakes1260.nc"


# Specify the NEW (subset extent) domain files:

	# Routing domain file
	subHydFile <- paste0(myPath, "/Fulldom.nc")

	# Geogrid domain file
	subGeoFile <- paste0(myPath, "/geo_em.nc")

	# Wrfinput domain file
	subWrfFile <- paste0(myPath, "/wrfinput.nc")

	# Route link file
	subRtlinkFile <- paste0(myPath, "/RouteLink.nc")

	# Spatial weights file
	subSpwtFile <- paste0(myPath, "/spatialweights.nc")

	# GW bucket parameter file
	subGwbuckFile <- paste0(myPath, "/GWBUCKPARM.nc")

	# Soil parameter file
	subSoilparmFile <- paste0(myPath, "/soil_properties.nc")

	# Lake parameter file
	subLakeparmFile <- paste0(myPath, "/LAKEPARM.nc")

	# Coordinate parameter text file
	subCoordParamFile <- paste0(myPath, "/params.txt")

	# Forcing clip script file
	subScriptFile <- paste0(myPath, "/script_forcing_subset.txt")



################ PROCESSING #####################

################ CALCULATE INDICES

library(rwrfhydro)
library(ncdf4)
source("Utils_ReachFiles.R")

# Setup coordinates df
coords <- data.frame(id=c(1,2,3,4), lat=c(y_south, y_north, y_north, y_south),
        lon=c(x_west, x_west, x_east, x_east))

# Create temp geogrid tif
tmpfile <- tempfile(fileext=".tif")
ExportGeogrid(fullGeoFile, "HGT_M", tmpfile)
geohgt <- raster::raster(tmpfile)
file.remove(tmpfile)

# Generate spatial coords
sp <- sp::SpatialPoints(data.frame(x=coords[,"lon"], y=coords[,"lat"]))
raster::crs(sp) <- coordProj
sp_proj <- sp::spTransform(sp, crs(geohgt))
geoindex <- as.data.frame(raster::rowColFromCell(geohgt, raster::cellFromXY(geohgt, sp_proj)))
geoindex$we <- geoindex$col

# Change row count from N->S to S->N
geoindex$sn <- dim(geohgt)[1] - geoindex$row + 1
geoindex$id <- coords[,"id"]

# Get subsetting dimensions
geo_w <- min(geoindex[,"we"])
geo_e <- max(geoindex[,"we"])
geo_s <- min(geoindex[,"sn"])
geo_n <- max(geoindex[,"sn"])
hyd_w <- (geo_w-1)*dxy+1
hyd_e <- geo_e*dxy
hyd_s <- (geo_s-1)*dxy+1
hyd_n <- geo_n*dxy
hyd_min <- (min(geoindex$row)-1)*dxy+1
hyd_max <- max(geoindex$row)*dxy
geo_min <- min(geoindex$row)
geo_max <- max(geoindex$row)

# Get relevant real coords for new bounds
geo_min_col <- min(geoindex[,"col"])
geo_max_col <- max(geoindex[,"col"])
geo_min_row <- min(geoindex[,"row"])
geo_max_row <- max(geoindex[,"row"])
rowcol_new <- data.frame(id=c(1,2,3,4), row=c(geo_max_row, geo_min_row, geo_min_row, geo_max_row),
        col=c(geo_min_col, geo_min_col, geo_max_col, geo_max_col))
rowcol_new_buff <- data.frame(id=c(1,2,3,4), row=c(geo_max_row+cellBuff, geo_min_row-cellBuff, geo_min_row-cellBuff, geo_max_row+cellBuff),
        col=c(geo_min_col-cellBuff, geo_min_col-cellBuff, geo_max_col+cellBuff, geo_max_col+cellBuff))
sp_new_proj <- xyFromCell(geohgt, raster::cellFromRowCol(geohgt, rowcol_new$row, rowcol_new$col), spatial=TRUE)
sp_new_proj_buff <- xyFromCell(geohgt, raster::cellFromRowCol(geohgt, rowcol_new_buff$row, rowcol_new_buff$col), spatial=TRUE)
sp_new_wrf <- sp::coordinates(sp::spTransform(sp_new_proj, "+proj=longlat +a=6370000 +b=6370000 +no_defs"))
sp_new_nad83 <- sp::coordinates(sp::spTransform(sp_new_proj, "+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"))
sp_new_buff_wrf <- sp::coordinates(sp::spTransform(sp_new_proj_buff, "+proj=longlat +a=6370000 +b=6370000 +no_defs"))
sp_new_buff_nad83 <- sp::coordinates(sp::spTransform(sp_new_proj_buff, "+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"))


################# SUBSET DOMAINS

# NCO starts with 0 index for dimensions, so we have to subtract 1

# ROUTING GRID

cmd <- paste0("ncks -d x,", hyd_w-1, ",", hyd_e-1, " -d y,", hyd_min-1, ",", hyd_max-1, " ", fullHydFile, " ", subHydFile)
print(cmd)
system(cmd)

# GEO GRID

# Dimension subsetting
cmd <- paste0("ncks -d west_east,", geo_w-1, ",", geo_e-1, " -d south_north,", geo_s-1, ",", geo_n-1, " ", fullGeoFile, " ", subGeoFile)
print(cmd)
system(cmd)
# Remove stagger dimensions since not used
cmd <- paste0("ncwa -O -a west_east_stag ", subGeoFile, " ", subGeoFile)
system(cmd)
cmd <- paste0("ncwa -O -a south_north_stag ", subGeoFile, " ", subGeoFile)
system(cmd)
# Attribute updates
cmd <- paste0("ncatted -h -a WEST-EAST_GRID_DIMENSION,global,o,l,", geo_e-geo_w+2, " ", subGeoFile)
system(cmd)
cmd <- paste0("ncatted -h -a SOUTH-NORTH_GRID_DIMENSION,global,o,l,", geo_n-geo_s+2, " ", subGeoFile)
system(cmd)
cmd <- paste0("ncatted -h -a WEST-EAST_PATCH_END_UNSTAG,global,o,l,", geo_e-geo_w+1, " ", subGeoFile)
system(cmd)
cmd <- paste0("ncatted -h -a SOUTH-NORTH_PATCH_END_UNSTAG,global,o,l,", geo_n-geo_s+1, " ", subGeoFile)
system(cmd)
cmd <- paste0("ncatted -h -a WEST-EAST_PATCH_START_STAG,global,d,,, ", subGeoFile)
system(cmd)
cmd <- paste0("ncatted -h -a SOUTH-NORTH_PATCH_START_STAG,global,d,,, ", subGeoFile)
system(cmd)
cmd <- paste0("ncatted -h -a WEST-EAST_PATCH_END_STAG,global,d,,, ", subGeoFile)
system(cmd)
cmd <- paste0("ncatted -h -a SOUTH-NORTH_PATCH_END_STAG,global,d,,, ", subGeoFile)
system(cmd)
cmd <- paste0("ncatted -h -a i_parent_end,global,o,l,", geo_e-geo_w+2, " ", subGeoFile)
system(cmd)
cmd <- paste0("ncatted -h -a j_parent_end,global,o,l,", geo_n-geo_s+2, " ", subGeoFile)
system(cmd)
cmd <- paste0("ncatted -h -a corner_lons,global,o,f,", sp_new_wrf[1,1], " ", subGeoFile)
system(cmd)
cmd <- paste0("ncatted -h -a corner_lats,global,o,f,", sp_new_wrf[1,2], " ", subGeoFile)
system(cmd)

# WRFINPUT GRID

cmd <- paste0("ncks -d west_east,", geo_w-1, ",", geo_e-1, " -d south_north,", geo_s-1, ",", geo_n-1, " ", fullWrfFile, " ", subWrfFile)
print(cmd)
system(cmd)
# Attribute updates
cmd <- paste0("ncatted -h -a WEST-EAST_GRID_DIMENSION,global,o,l,", geo_e-geo_w+2, " ", subWrfFile)
system(cmd)
cmd <- paste0("ncatted -h -a SOUTH-NORTH_GRID_DIMENSION,global,o,l,", geo_n-geo_s+2, " ", subWrfFile)
system(cmd)


################# SUBSET PARAMS

# Identify catchments to keep

fullWts <- ReadWtFile(fullSpwtFile)
keepIdsPoly <- subset(fullWts[[1]], fullWts[[1]]$i_index >= hyd_w & fullWts[[1]]$i_index <= hyd_e &
			fullWts[[1]]$j_index >= hyd_s & fullWts[[1]]$j_index <= hyd_n)
keepIdsPoly <- unique(keepIdsPoly$IDmask)

fullRtlink <- ReadLinkFile(fullRtlinkFile)
#keepIdsLink <- subset(fullRtlink, fullRtlink$lon >= min(sp_new_buff_nad83[,1]) & fullRtlink$lon <= max(sp_new_buff_nad83[,1]) &
#			fullRtlink$lat >= min(sp_new_buff_nad83[,2]) & fullRtlink$lat <= max(sp_new_buff_nad83[,2]))
keepIdsLink <- subset(fullRtlink, fullRtlink$x >= x_west & fullRtlink$x <= x_east &
                      fullRtlink$y >= y_south & fullRtlink$y <= y_north)
keepIdsLink <- unique(keepIdsLink$link)

keepIds <- unique(c(keepIdsPoly, keepIdsLink))

# SPATIAL WEIGHT

subWts <- SubsetWts(fullWts, keepIdsPoly, hyd_w, hyd_e, hyd_s, hyd_n)
file.copy(fullSpwtFile, subSpwtFile)
UpdateWtFile(subSpwtFile, subWts[[1]], subWts[[2]], subDim=TRUE)

# GWBUCK PARAMETER

fullGwbuck <- GetNcdfFile(fullGwbuckFile, quiet=TRUE)
subGwbuck <- subset(fullGwbuck, fullGwbuck$ComID %in% keepIdsPoly)
subGwbuck$Basin <- seq(1, nrow(subGwbuck), 1)
file.copy(fullGwbuckFile, subGwbuckFile)
UpdateGwbuckFile(subGwbuckFile, subGwbuck, subDim=TRUE)

# ROUTE LINK

subRtlink <- subset(fullRtlink, fullRtlink$link %in% keepIds)
subRtlink$to <- ifelse(subRtlink$to %in% unique(subRtlink$link), subRtlink$to, 0)
file.copy(fullRtlinkFile, subRtlinkFile)
UpdateLinkFile(subRtlinkFile, subRtlink, subDim=TRUE)

# SOIL PARAMETER

cmd <- paste0("ncks -d west_east,", geo_w-1, ",", geo_e-1, " -d south_north,", geo_s-1, ",", geo_n-1, " ", fullSoilparmFile, " ", subSoilparmFile)
system(cmd)

# LAKE PARAMETER

# nothing for now...
file.copy(fullLakeparmFile, subLakeparmFile)


################# CREATE SCRIPT FILES

# Save the coordinate parameter file
coordsExport <- data.frame(grid=c("hyd_sn", "hyd_ns", "geo_sn", "geo_ns"),
                            imin=c(hyd_w, hyd_w, geo_w, geo_w),
                            imax=c(hyd_e, hyd_e, geo_e, geo_e),
                            jmin=c(hyd_s, hyd_min, geo_s, geo_min),
                            jmax=c(hyd_n, hyd_max, geo_n, geo_max),
                            index_start=c(1,1,1,1))
write.table(coordsExport, file=subCoordParamFile, row.names=FALSE, sep="\t")

# Save the forcing subset script file
#ncksCmd <- paste0("ncks -d west_east,", geo_w-1, ",", geo_e-1, " -d south_north,", geo_s-1, ",", geo_n-1, " ${OLDFORCPATH}/${i} ${NEWFORCPATH}/${i}")
ncksCmd <- paste0("ncks -d west_east,", geo_w-1, ",", geo_e-1, " -d south_north,", geo_s-1, ",", geo_n-1, " ${i} ${NEWFORCPATH}/${i##*/}")
fileConn <- file(subScriptFile)
writeLines(c("#!/bin/bash",
		"OLDFORCPATH='PATH_TO_OLD_FORCING_DATA_FOLDER'",
		"NEWFORCPATH='PATH_TO_NEW_FORCING_DATA_FOLDER'",
		"for i in `ls $OLDFORCPATH`; do",
		"echo ${i##*/}",
		ncksCmd,
		"done"),
		fileConn)
close(fileConn)

ncksCmd <- paste0("ncks -d ncl0,", geo_s-1, ",", geo_n-1, " -d ncl1,", geo_w-1, ",", geo_e-1, " -d ncl2,", geo_s-1, ",", geo_n-1, " -d ncl3,", geo_w-1, ",", geo_e-1, " -d ncl4,", geo_s-1, ",", geo_n-1, " -d ncl5,", geo_w-1, ",", geo_e-1, " -d ncl6,", geo_s-1, ",", geo_n-1, " -d ncl7,", geo_w-1, ",", geo_e-1, " -d ncl8,", geo_s-1, ",", geo_n-1, " -d ncl9,", geo_w-1, ",", geo_e-1, " -d ncl10,", geo_s-1, ",", geo_n-1, " -d ncl11,", geo_w-1, ",", geo_e-1, " -d ncl12,", geo_s-1, ",", geo_n-1, " -d ncl13,", geo_w-1, ",", geo_e-1, " -d ncl14,", geo_s-1, ",", geo_n-1, " -d ncl15,", geo_w-1, ",", geo_e-1, " ${i} ${NEWFORCPATH}/${i##*/}")
fileConn <- file(paste0(subScriptFile, "_NWM_REALTIME"))
writeLines(c("#!/bin/bash",
                "OLDFORCPATH='PATH_TO_OLD_FORCING_DATA_FOLDER'",
                "NEWFORCPATH='PATH_TO_NEW_FORCING_DATA_FOLDER'",
                "for i in `ls $OLDFORCPATH`; do",
                "echo ${i##*/}",
                ncksCmd,
                "done"),
                fileConn)
close(fileConn)


#quit("no")

