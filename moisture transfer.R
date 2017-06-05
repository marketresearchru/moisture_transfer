library("ncdf4")
setwd("E:/ncep_ncar20C")
Dim1 <- 3	# take Dim1 cells by vertical dimension LAT

Lon <- 18 # index, lon to the right is the base for calculations of transfer
Levels <- 19 # count of levels in a column

# Take first NC file to create header data
ncfile <- nc_open("special_humidity/shum.1901.nc")
nc.lat <- ncvar_get(ncfile, "lat", start = 21, count = Dim1)
nc.lon <- ncvar_get(ncfile, "lon")
nc.time <- ncvar_get(ncfile, "time")
nc.levels <- ncvar_get(ncfile, "level") # 19 levels from 1000 with step of 50
nc_close(ncfile)

# Output file
fileconnector <- file(description = paste0("c:/temp/monthly transfer via ", nc.lon[Lon]+1, ".txt"), open="wt")
write(c("YEAR", "Vertical", 1:12), file=fileconnector, append=TRUE, ncolumns=14, sep = "\t")

for(year in 1901:2012) {
	cat("year ", year, "\n")
	tm.len <- 4 * ifelse((year - floor(year/4)*4) == 0, 366, 365)
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

	# Humidity data
	ncfile <- nc_open(paste0("special_humidity/shum.", year, ".nc"))
	#print(ncfile)
	# read full year - 2 vertical coords (LONs) from 10th
	h1 <- ncvar_get(ncfile, "shum", start=c(Lon, 4, 1, 1), count=c(2, Dim1, 19, tm.len ) )
	nc_close(ncfile)
	
	# u-wind component, m/s 
	ncfile <- nc_open(paste0("uwnd/uwnd.", year, ".nc"))
	wind.u <- ncvar_get(ncfile, "uwnd", start=c(Lon, 4, 1, 1), count=c(2, Dim1, 19, tm.len ) )
	nc_close(ncfile)
	
	transfer <- array(dim=12)
	cat("months: ")
	for(month in 1:12) {
		cat(month)
		transfer[month] <- 0
		for(time in which(month.sroks==month, arr.ind = TRUE)){ #for each moment of time
			for(lat in 1:Dim1){ #for each vertical step
				m.left <- h1[1,lat,1:Levels, time]
				m.right <- h1[2,lat,1:Levels, time]
				wind.left <- wind.u[1,lat,1:Levels, time]
				wind.right <- wind.u[2,lat,1:Levels, time]
				#humidity volume at 1000h level
				level <- 1
				if(wind.left[level] > 0 ){
					# wind to East, account amount in left column
					transfer[month] <- transfer[month] + wind.left[level] * 25 * ((m.left[level]+m.left[level+1])/2 + m.left[level])
				} 
				if(wind.right[level] < 0 ){
					# wind to West, account amount in right column as negative
					transfer[month] <- transfer[month] + wind.right[level] * 25 * ((m.right[level]+m.right[level+1])/2+m.right[level])
				}
				
				#volume of humidity on all other levels
				for(level in 2:18){
					if(wind.left[level] > 0 ){
						# wind to East, account amount in left column
						transfer[month] <- transfer[month] + wind.left[level] * 25 * (3*m.left[level]+(m.left[level+1] + m.left[level-1])/2)
					} 
					if(wind.right[level] < 0 ){
						# wind to West, account amount in right column as negative
						transfer[month] <- transfer[month] + wind.right[level] * 25 * (3*m.right[level]+(m.right[level+1] + m.right[level-1])/2)
					}
				}

				#volume at 100 level
				level <- 19
				if(wind.left[level] > 0 ){
					# wind to East, account amount in left column
					transfer[month] <- transfer[month] + wind.left[level] * 50 * ((m.left[level]+m.left[level-1])/2+m.left[level])/2
				} 
				if(wind.right[level] < 0 ){
					# wind to West, account amount in right column as negative
					transfer[month] <- transfer[month] + wind.right[level] * 50 * ((m.right[level]+m.right[level-1])/2+m.right[level])/2
				}
			}
		} # end of loop by all sroks within month
		
	} # end of loop by month
	write(c(year, nc.lon[Lon], transfer), file=fileconnector, append=TRUE, ncolumns=14, sep = "\t")
	cat(" year complete\n")
} # end of loop by year

close(fileconnector)

