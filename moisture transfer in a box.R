# This program calculates moisture transfer over vertical or horisontal line
# Transfer is based on wind and humidity on 19 levels in the atmosphere stack

# we consider rectangle outside of the border which is researched
# in case of this particular grid, we have only odd degrees, therfore
# actual border is the even degree between outer definded below and inner 
# (2 degrees smaller in each direction)
# LONGITUDE 1 & 2 - X degree to the east
# LATITUDE 1 & 2 - Y degree to the north 

LONGITUDE1 = 34 # east
LONGITUDE2 = 60
LATITUDE1 = 46 # north
LATITUDE2 = 60

# the only library required is NETCDF opener
library("ncdf4")
# source data folder
setwd("E:/ncep_ncar20C")

Levels <- 19 # count of levels in a column

# Take first NC file to create header data
ncfile <- nc_open("special_humidity/shum.1901.nc")
nc.lat <- ncvar_get(ncfile, "lat")
Lat1 <- which(abs(nc.lat-LATITUDE1)<0.5)-1 # shift 1 step top (inside)
Lat2 <- which(abs(nc.lat-LATITUDE2)<0.5) # Lat2 < Lat1 because data in NETCDF starts from 90 and descends

nc.lon <- ncvar_get(ncfile, "lon")
Lon1 <- which(abs(nc.lon-LONGITUDE1)<0.5)
Lon2 <- which(abs(nc.lon-LONGITUDE2)<0.5)-1 # shift 1 step to left (inside)

nc.time <- ncvar_get(ncfile, "time")
nc.levels <- ncvar_get(ncfile, "level") # 19 levels from 1000 with step of 50
nc_close(ncfile)

# month lengths (each day has 4 measures)
# this function below takes in account visokosny year, however, they are not in the file
month.intervals <- function(year, tm.len){
	month.sroks <- array(dim=tm.len)
	mn2 <- 4*31
	month.sroks[1:mn2] <- 1
	mn1 <- mn2 + 1
	mn2 <- mn1 - 1 + ifelse((year - floor(year/4)*4) == 0, 29*4, 28*4) 
	month.sroks[mn1:mn2] <- 2
	mn1 <- mn2 + 1
	mn2 <- mn1 - 1 + 31*4
	month.sroks[mn1:mn2] <- 3
	mn1 <- mn2 + 1
	mn2 <- mn1 - 1 + 30*4
	month.sroks[mn1:mn2] <- 4 # April
	mn1 <- mn2 + 1
	mn2 <- mn1 - 1 + 31*4
	month.sroks[mn1:mn2] <- 5
	mn1 <- mn2 + 1
	mn2 <- mn1 - 1 + 30*4
	month.sroks[mn1:mn2] <- 6
	mn1 <- mn2 + 1
	mn2 <- mn1 - 1 + 31*4
	month.sroks[mn1:mn2] <- 7
	mn1 <- mn2 + 1
	mn2 <- mn1 - 1 + 31*4
	month.sroks[mn1:mn2] <- 8
	mn1 <- mn2 + 1
	mn2 <- mn1 - 1 + 30*4
	month.sroks[mn1:mn2] <- 9
	mn1 <- mn2 + 1
	mn2 <- mn1 - 1 + 31*4
	month.sroks[mn1:mn2] <- 10
	mn1 <- mn2 + 1
	mn2 <- mn1 - 1 + 30*4
	month.sroks[mn1:mn2] <- 11
	mn1 <- mn2 + 1
	mn2 <- mn1 - 1 + 31*4
	month.sroks[mn1:mn2] <- 12
	
	return(month.sroks)
}



# Transfer through horisontal border for particular srok
srok.h.transfer <- function(h, w){
	trans <- 0
	steps <- dim(h)[1]
	for(l in 1:steps){ #for each l on horizontal border
		m.left <- h[l,1,]
		m.right <- h[l,2,]
		wind.left <- w[l,1,]
		# replace wind with 0 if it blows out of border (leave only those wich really transfers moisture)
		wind.left[wind.left>0] <- 0
		wind.right <- w[l,2,]
		wind.right[wind.right<0] <- 0
		#humidity volume at 1000h level
		level <- 1
		trans <- trans + wind.left[level]  * 25 * ((m.left[level]+m.left[level+1])/2 + m.left[level])
		trans <- trans + wind.right[level] * 25 * ((m.right[level]+m.right[level+1])/2+m.right[level])
		
		#volume of humidity on all other levels
		for(level in 2:18){
			trans <- trans + wind.left[level]  * 25 * (3*m.left[level]+(m.left[level+1] + m.left[level-1])/2)
			trans <- trans + wind.right[level] * 25 * (3*m.right[level]+(m.right[level+1] + m.right[level-1])/2)
		}
		
		#volume at 100 level
		level <- 19
		trans <- trans + wind.left[level]  * 50 * ((m.left[level]+m.left[level-1])/2+m.left[level])/2
		trans <- trans + wind.right[level] * 50 * ((m.right[level]+m.right[level-1])/2+m.right[level])/2
	}
	trans
}	

# Transfer through vertical border	
srok.v.transfer <- function(h, w){
	trans <- 0
#	cat(dim(h))
#	cat(dim(w))
	steps <- dim(h)[2]
	for(l in 1:steps){ #for each l on vertical border
		m.left <- h[1,l,]  #left 
		m.right <- h[2,l,]
		wind.left <- w[1,l,]
		# replace wind with 0 if it blows out of border (leave only those wich really transfers moisture)
		wind.left[wind.left<0] <- 0
		wind.right <- w[2,l,]
		wind.right[wind.right>0] <- 0
		#humidity volume at 1000h level
		level <- 1
		trans <- trans + wind.left[level]  * 25 * ((m.left[level]+m.left[level+1])/2 + m.left[level])
		trans <- trans + wind.right[level] * 25 * ((m.right[level]+m.right[level+1])/2+m.right[level])
		
		#volume of humidity on all other levels
		for(level in 2:18){
			trans <- trans + wind.left[level]  * 25 * (3*m.left[level]+(m.left[level+1] + m.left[level-1])/2)
			trans <- trans + wind.right[level] * 25 * (3*m.right[level]+(m.right[level+1] + m.right[level-1])/2)
		}
		
		#volume at 100 level
		level <- 19
		trans <- trans + wind.left[level]  * 50 * ((m.left[level]+m.left[level-1])/2+m.left[level])/2
		trans <- trans + wind.right[level] * 50 * ((m.right[level]+m.right[level-1])/2+m.right[level])/2
	}
	trans
}	



DimLat <- Lat1-Lat2
DimLon <- Lon2-Lon1

# Output file 
fileconnector <- file(description = "C:/Users/alexe/Documents/Climat/Wind-and-humidity/monthly box transfer.txt", open="wt")
write(c("YEAR", 1:12), file=fileconnector, append=TRUE, ncolumns=13, sep = "\t")

year <- 1901
for(year in 1901:2012) {
	cat("year", year)
	tm.len <- 4 * ifelse((year - floor(year/4)*4) == 0, 366, 365)
	month.sroks <- month.intervals(year, tm.len)
	# Humidity data
	ncfileh <- nc_open(paste0("special_humidity/shum.", year, ".nc"))
	ncfileu <- nc_open(paste0("uwnd/uwnd.", year, ".nc"))
	ncfilev <- nc_open(paste0("vwnd/vwnd.", year, ".nc"))
	#print(ncfile)
	# read full year, 4 blocks

	# block L - left border of a box + 1 vertical to the right (east)
	# humidity
	hL <- ncvar_get(ncfileh, "shum", start=c(Lon1, Lat2+1, 1, 1), count=c(2, DimLat, 19, tm.len) )
	# u-wind component, m/s 
	wind.uL <- ncvar_get(ncfileu, "uwnd", start=c(Lon1, Lat2+1, 1, 1), count=c(2, DimLat, 19, tm.len ) )

	# block R - right border of a box + 1 vertical to the right
	# humidity
	hR <- ncvar_get(ncfileh, "shum", start=c(Lon2, Lat2+1, 1, 1), count=c(2, DimLat, 19, tm.len ) )
	# u-wind component, m/s 
	wind.uR <- ncvar_get(ncfileu, "uwnd", start=c(Lon2, Lat2+1, 1, 1), count=c(2, DimLat, 19, tm.len) )
	
	# block T - top border of a box + 1 horisontal to the bottom (south)
	# humidity
	hT <- ncvar_get(ncfileh, "shum", start=c(Lon1+1, Lat2, 1, 1), count=c(DimLon, 2, 19, tm.len) )
	# v-wind component, m/s 
	wind.vT <- ncvar_get(ncfilev, "vwnd", start=c(Lon1+1, Lat2, 1, 1), count=c(DimLon, 2, 19, tm.len) )
	
	# block B - bottom border of a box + 1 horisontal to the bottom (south)
	# humidity
	hB <- ncvar_get(ncfileh, "shum", start=c(Lon1+1, Lat1, 1, 1), count=c(DimLon, 2, 19, tm.len) )
	# v-wind component, m/s 
	wind.vB <- ncvar_get(ncfilev, "vwnd", start=c(Lon1+1, Lat1, 1, 1), count=c(DimLon, 2, 19, tm.len) )

	nc_close(ncfileh)
	nc_close(ncfileu)
	nc_close(ncfilev)

	transfer <- array(dim=12)
	cat(" Data taken\nprocessing...")
	for(month in 1:12) {
		cat(month.abb[month], " ")
		transfer[month] <- 0
		for(time in which(month.sroks==month, arr.ind = TRUE)){ #for each moment of time
			transfer[month] <- transfer[month] + srok.v.transfer(hL[,,,time],wind.uL[,,,time])
			transfer[month] <- transfer[month] + srok.v.transfer(hR[,,,time],wind.uR[,,,time])
			transfer[month] <- transfer[month] + srok.h.transfer(hT[,,,time],wind.vT[,,,time])
			transfer[month] <- transfer[month] + srok.h.transfer(hB[,,,time],wind.vB[,,,time])
		} # end of loop by all sroks within month
		
	} # end of loop by month
	
	write(c(year, transfer), file=fileconnector, append=TRUE, ncolumns=13, sep = "\t")
	cat("year done\n")
} # end of loop by year
write(paste("Outer rectangle: LONGITUDE = [", LONGITUDE1, ",", LONGITUDE2, "]\tLATITUDE = [", LATITUDE1, ",", LATITUDE2, "]"), file=fileconnector, append=TRUE)

close(fileconnector)

