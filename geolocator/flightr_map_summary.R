# playing around with FlightR for plotting map and determining stopover locations for Purple Martin geolocator data
# followed vignette for main FLightR functions from Rakhimberdiev et al. 2017 (https://cran.r-project.org/web/packages/FLightR/vignettes/FLightR_with_black-tailed_godwit_vignette_from_MEE_2017.html)
# modified code and set up for my data needs - Nov. 2020

library(FLightR)
library(devtools)
library(tidyverse)
library(ggmap)
library(grid)
library(ggplot2)

#-----------input info----------

#set working directory
setwd("C:/Users/eveli/Dropbox/PUMA/Geolocator/FlightR/newbatch") 

#put in band number and year info
id="158135707_12-13" 

#input breeding colony coordinates
longitude=-112.863565 
latitude=53.011129
  
#calibration period (last two weeks of breeding season, this will be different for each colony location)
calibration.start<- as.POSIXct("2012-07-20 00:00:00", "%Y-%m-%d %H:%M:%S", tz = "GMT")
calibration.stop<- as.POSIXct("2012-08-03 00:00:00", "%Y-%m-%d %H:%M:%S", tz = "GMT")

#define start and end date for migration period (can look at beginning and end of frametime twilights timeline)
start_date <- as.POSIXct("2012-07-17 05:00:00", "%Y-%m-%d %H:%M:%S", tz = "GMT")
end_date <- as.POSIXct("2013-06-10 00:00:00", "%Y-%m-%d %H:%M:%S", tz = "GMT")

#enter Google API key if you want to make map later
register_google(key="##############################")
has_google_key() #double check that is it "TRUE"

#setting up function for loading lig file
readLig <- function(file,skip=1) {
  # Read csv file and add column names
  d <- read.csv(file,header=F,skip=skip,
                col.names=c("Valid","Date","Julian","Light"),
                colClasses=c("character","character","numeric","integer"))
  # Parse date
  d$Date <- as.POSIXct(strptime(d$Date,"%d/%m/%y %H:%M:%S",tz="GMT"))
  d
}

#load lig file
raw_lightdata <- readLig(file.choose())


#-----------time to run the things!----------

#subset only "Date" and "Light" column
raw_lightdata <- subset(raw_lightdata, select=c("Date", "Light"))

#load twilight time data from bastag/geolight analysis
twilights <- read.csv(paste(id,"_frametime.csv",sep=""))

#putting dates and times into POSIX format
twilights$tFirst <- as.POSIXct(twilights$tFirst, "%Y-%m-%d %H:%M", tz="GMT")
twilights$tSecond <- as.POSIXct(twilights$tSecond, "%Y-%m-%d %H:%M", tz="GMT")

#converting data into a FLIGHTR data object
TAGS_twilights <- GeoLight2TAGS(raw_lightdata, twilights, threshold=32)
write.csv(TAGS_twilights, paste(id,".TAGS_twilights.csv",sep=""))
flightr_data <- get.tags.data(paste(id,".TAGS_twilights.csv",sep=""), start_date, end_date, measurement.period = 120)

#setting up calibration
calibration.periods <- data.frame(calibration.start, calibration.stop, lon=longitude,lat=latitude)
print(calibration.periods) #just to double check
calibration <- make.calibration(flightr_data, calibration.periods)

#assign spatial extent, can adjust as desired
grid <-make.grid(left = -115, bottom = -35, right = -35, top = 70,
                 distance.from.land.allowed.to.use = c(-Inf, Inf),
                 distance.from.land.allowed.to.stay = c(-Inf, Inf))

#combining the input data together to prep data for more flightr things. This can take ~45 min to run.
combine_all <- make.prerun.object(flightr_data, grid, start=c(longitude,latitude), Calibration=calibration)

#run the particle filter step, can take a while maybe ~90-120min.
geo_result <- run.particle.filter(combine_all,
                              threads = -1,
                              nParticles = 1e6,
                              known.last = TRUE,
                              precision.sd = 25,
                              check.outliers = FALSE)

save(geo_result, file=paste(id,".flightr.particle.filter.RData",sep=""))
dev.copy(png,paste(id,".particlefilter.map.png",sep=""))
dev.off()


#plot graph of longitude and latitude over track timeline
png(file=paste(id,".lonlatplot.png",sep=""))
plot_lon_lat(geo_result)
dev.off()

#basic plot map
#map.FLightR.ggmap(geo_result)

#plot map with some more specifications, adjust as desired)
map <- map.FLightR.ggmap(geo_result, zoom=3, seasonal.colors=TRUE, seasonal.donut.location="bottomleft", seasonal.donut.proportion=0.4, save=FALSE)
ggsave(paste(id,".ggmap.png",sep=""))


#derive time when bird arrived/departed area
index <- which(geo_result$Spatial$Grid[,2]>(latitude))
breeding_arrival <- find.times.distribution(geo_result,index)
print(breeding_arrival)
write.csv(breeding_arrival, paste(id,".find.times.distribution.csv",sep=""))


#summarize residency, stopovers, and migration periods (can take ~15min). 
#creates list of 3 (Stationary.periods, Potential_stat_periods, and Potential_movement_periods)
summary <- stationary.migration.summary(geo_result, prob.cutoff = 0.1, min.stay = 3)

capture.output(summary$Stationary.periods, file=paste(id,".Stationary.periods.txt",sep=""))
capture.output(summary$Potential_stat_periods, file=paste(id,".Potential.stat.periods.txt",sep=""))
capture.output(summary$Potential_movement_periods, file=paste(id,".Potential.movement.periods.txt",sep=""))
