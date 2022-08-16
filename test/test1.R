# me trying stuff

install.packages("ecmwfr")
install.packages("rNOMADS")
install.packages("ncdf4")
install.packages("ncdf")
install.packages("tidync")
install.packages("rgdal")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(rNOMADS)
x <- ReadGrib("adaptor.mars.internal-1650164883.1168177-27172-11-b22d2b89-43d3-46de-a7cd-35cb5080697d.grib",
              file.type = "grib")





# ERA-interim
# https://www.r-bloggers.com/2018/09/access-to-climate-reanalysis-data-from-r/

if(!require("reticulate")) install.packages("reticulate")
if(!require("ncdf4")) install.packages("ncdf4") #to manage netCDF format

#load packages
library(reticulate)
library(ncdf4)

#install the python ECMWF API
py_install("ecmwf-api-client")
## 
## Installation complete.





# ERA5
# https://www.r-bloggers.com/2018/09/access-to-climate-reanalysis-data-from-r/

#load libraries 
library(sf)
library(ncdf4)
library(tidyverse)
library(lubridate)
library(reticulate)

#install the CDS API
conda_install("r-reticulate","cdsapi",pip=TRUE)

#import python CDS-API
cdsapi <- import('cdsapi')

#for this step there must exist the file .cdsapirc
server = cdsapi$Client() #start the connection

#we create the query
query <- r_to_py(list(
  variable= "2m_temperature",
  product_type= "reanalysis",
  year= "2018",
  month= "07", #formato: "01","01", etc.
  day= str_pad(1:31,2,"left","0"),   
  time= str_c(0:23,"00",sep=":")%>%str_pad(5,"left","0"),
  format= "netcdf",
  area = "45/-20/35/5" # North, West, South, East
))

#query to get the ncdf
server$retrieve("reanalysis-era5-single-levels",
                query,
                "era5_ta_2018.nc")

#open the connection with the file
nc <- nc_open("era5_ta_2018.nc")

#extract lon, lat
lat <- ncvar_get(nc,'latitude')
lon <- ncvar_get(nc,'longitude')
dim(lat);dim(lon)
## [1] 41
## [1] 101
#extract time
t <- ncvar_get(nc, "time")

#time unit: hours from 1900-01-01
ncatt_get(nc,'time')
## $units
## [1] "hours since 1900-01-01 00:00:00.0"
## 
## $long_name
## [1] "time"
## 
## $calendar
## [1] "gregorian"
#we convert the hours into date+time 
#as_datetime from lubridate needs seconds
timestamp <- as_datetime(c(t*60*60),origin="1900-01-01")

#temperatures in K from july 2018
head(timestamp)
## [1] "2018-07-01 00:00:00 UTC" "2018-07-01 01:00:00 UTC"
## [3] "2018-07-01 02:00:00 UTC" "2018-07-01 03:00:00 UTC"
## [5] "2018-07-01 04:00:00 UTC" "2018-07-01 05:00:00 UTC"
#import temperature data
data <- ncvar_get(nc,"t2m")

#plot 2018-07-01
filled.contour(data[,,1])

#time serie plot for a pixel
plot(data.frame(date=timestamp,
                ta=data[1,5,]),
     type="l")

#close the conection with the ncdf file
nc_close(nc)







# ERA5 multiple levels

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#load libraries 
library(sf)
library(ncdf4)
library(tidyverse)
library(lubridate)
library(reticulate)

#install the CDS API
#conda_install("r-reticulate","cdsapi",pip=TRUE)

#import python CDS-API
cdsapi <- import('cdsapi')

#for this step there must exist the file .cdsapirc
server = cdsapi$Client() #start the connection

#we create the query
query <- r_to_py(list(
  variable= c("u_component_of_wind", "v_component_of_wind"),
  pressure_level= c("1000","850"),
  product_type= "reanalysis",
  year= "2018",
  month= "07", #formato: "01","01", etc.
  day= str_pad(1:31,2,"left","0"),   
  time= str_c(0:23,"00",sep=":")%>%str_pad(5,"left","0"),
  format= "netcdf",
  area = "17/-91/7/-81" # North, West, South, East
))

#query to get the ncdf
server$retrieve("reanalysis-era5-pressure-levels",
                query,
                "era5_test1.nc")

#open the connection with the file
nc <- nc_open("era5_test1.nc")

#extract lon, lat
lat <- ncvar_get(nc,'latitude')
lon <- ncvar_get(nc,'longitude')
dim(lat);dim(lon)
## [1] 41
## [1] 101
#extract time
t <- ncvar_get(nc, "time")

#time unit: hours from 1900-01-01
ncatt_get(nc,'time')
## $units
## [1] "hours since 1900-01-01 00:00:00.0"
## 
## $long_name
## [1] "time"
## 
## $calendar
## [1] "gregorian"
#we convert the hours into date+time 
#as_datetime from lubridate needs seconds
timestamp <- as_datetime(c(t*60*60),origin="1900-01-01")

lev <- ncvar_get(nc,"level")

#temperatures in K from july 2018
head(timestamp)
## [1] "2018-07-01 00:00:00 UTC" "2018-07-01 01:00:00 UTC"
## [3] "2018-07-01 02:00:00 UTC" "2018-07-01 03:00:00 UTC"
## [5] "2018-07-01 04:00:00 UTC" "2018-07-01 05:00:00 UTC"
#import temperature data
data <- ncvar_get(nc,"u")

#plot 2018-07-01
filled.contour(data[,,1,1])

#time serie plot for a pixel
plot(data.frame(date=timestamp,
                ua=data[1,5,1,]),
     type="l")

#close the conection with the ncdf file
nc_close(nc)
