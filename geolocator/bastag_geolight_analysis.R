#Script to analyze geolocator data from raw light level files (.lig) to obtain estimates of daily locations

library(BAStag) 
library(tidyverse) 
library(useful)


#Importing lig file:
readLig <- function(file,skip=1) {
  # Read csv file and add column names
  d <- read.csv(file,header=F,skip=skip,
                col.names=c("Valid","Date","Julian","Light"),
                colClasses=c("character","character","numeric","integer"))
  # Parse date
  d$Date <- as.POSIXct(strptime(d$Date,"%d/%m/%y %H:%M:%S",tz="GMT"))
  d
}

PumaRaw <- readLig(file.choose()) #upload the lig file for specific geolocator
head(PumaRaw) #double check that it loaded correctly


#Determining light levels - work through interactive stages in preprocessLight to edit twilights.
#light threshold for definition of twilight = 32 for purple martin
PUMA_transitions <- preprocessLight(tagdata=PumaRaw, threshold=32, offset=19, zlim=c(0,12), dark.min=240)


#Formatting twilight data into readable format for GeoLight
#shift the Twilight3 column to make it into the tSecond column later
PUMA_transitions <- shift.column(data=PUMA_transitions, columns="Twilight3", len=1, up=TRUE) 

#making "Rise" into  "type" column
PUMA_transitions <- PUMA_transitions %>%
  dplyr::mutate(type = ifelse(Rise == TRUE, "1", "2"))

#removing rows of deleted twilights
PUMA_transitions <- subset(PUMA_transitions, Deleted!="TRUE")

#Selecting the 3 columns we want
PUMA_transitions <- subset(PUMA_transitions, select=c("Twilight", "Twilight3.Shifted", "type"))

#renaming columns
PUMA_transitions <- select(PUMA_transitions, tFirst=Twilight, tSecond=Twilight3.Shifted, type=type)

#making tFirst and tSecond confirmed POSIXct format
PUMA_transitions$tFirst <- as.POSIXct(PUMA_transitions$tFirst, "%Y-%m-%d %H:%M", tz="GMT")
PUMA_transitions$tSecond <- as.POSIXct(PUMA_transitions$tSecond, "%Y-%m-%d %H:%M", tz="GMT")

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Next using GeoLight. Code based from Amanda Shave, and on tutorial: http://scbi-migbirds.github.io/Geolocator_GeoLight.html by Michael T. Hallworth.

library(maps) 
library(GeoLight)

#Calculating sun elevation angle. 
#Need to select appropriate interval and correct coordinates for calibrating the sun elevation angle. This example is using last week of July till first week of August in Manitoba. 

SunElev<-getElevation(tFirst=PUMA_transitions[23:51,1],
                      tSecond=PUMA_transitions[23:51,2],
                      type=PUMA_transitions[23:51,3],
                      known.coord=c(-112.863565,53.011129),
                      plot=TRUE)

SunElev #double check that it worked


#Estimating location from consecutive twilights. The tol argument in this coord function defines how many positions will be discarded around the equinox period
PUMA_Locations<-coord(tFirst=PUMA_transitions[,1],
                       tSecond=PUMA_transitions[,2],
                       type=PUMA_transitions[,3], 
                       degElevation=SunElev, 
                       tol=0.13)  

head(PUMA_Locations)

#Can make a quick preliminary plot to make sure it looks reasonable
plot(PUMA_Locations, pch="*", col="red", xlab="Longitude", ylab="Latitude")
map("world",add=TRUE)


#Using changeLight fuction in GeoLight to distinguish between residency and movement periods. Adjust number of days based on how you want to define how long minimum would be residency period

stop<-changeLight(tFirst=PUMA_transitions[,1], tSecond=PUMA_transitions[,2], type=2, twl=PUMA_transitions, quantile=0.9, rise.prob=NA, set.prob=NA, days=1, plot=TRUE, summary=TRUE)

#This is an example of the output to expect, should probably copy and paste & save somewhere
#Probability threshold(s):

#  Sunrise:  0.03653	Sunset:  0.04464

#Migration schedule table:
  
#  Site             Arrival           Departure  Days     P.start       P.end Days.1   P.start.1
#1     a 2016-07-06 16:49:46 2016-07-08 17:17:19   2.0 0.023958333 0.000000000    2.0 0.023958333
#2     b 2016-07-12 06:32:53 2016-07-29 06:34:12  17.0 0.000000000 0.021896660   17.0 0.000000000
#3     c 2016-07-30 06:30:05 2016-08-02 18:33:48   3.5 0.000000000 0.000000000    3.5 0.000000000
#4     d 2016-08-03 17:42:36 2016-08-11 06:42:16   7.5 0.004333333 0.000000000    7.5 0.004333333
#5     e 2016-08-12 06:48:16 2016-08-25 06:39:46  13.0 0.000000000 0.000000000   13.0 0.000000000
#6     f 2016-08-26 06:44:01 2016-08-28 18:43:24   2.5 0.004166667 0.004666667    2.5 0.004166667
#7     g 2016-08-31 18:22:08 2016-09-03 18:22:04   3.0 0.000000000 0.015740741    3.0 0.000000000
#8     h 2016-09-04 18:06:54 2016-09-11 05:34:47   6.5 0.000000000 0.000000000    6.5 0.000000000
#9     i 2016-09-15 04:47:35 2016-09-19 04:26:30   4.0 0.000000000 0.000000000    4.0 0.000000000
#10    j 2016-09-20 04:26:51 2017-04-11 04:02:01 203.0 0.003703704 0.000000000  203.0 0.003703704
#11    k 2017-04-17 05:39:36 2017-04-20 17:58:10   3.5 0.000000000 0.004166667    3.5 0.000000000
#12    l 2017-04-24 18:31:29 2017-04-28 06:18:05   3.5 0.000000000 0.000000000    3.5 0.000000000
#13    m 2017-05-02 06:21:25 2017-05-06 18:01:58   4.5 0.000000000 0.025925926    4.5 0.000000000
#14    n 2017-05-11 04:08:25                <NA>    NA 0.000000000 0.000000000     NA 0.000000000



#Additional mapping options:

#Plotting positions and trip on map
tripMap(PUMA_Locations, equinox = TRUE, map.range = c("America"), legend = TRUE)

#Plotting sites of residency
siteMap(PUMA_Locations, stop$site, type = "points", quantiles = c(0.25, 0.75), hull = T, map.range = c("America"))


#Convert location and timing data to dataframes. When merging these together, can look at daily location estimates.
framelocation<-data.frame(PUMA_Locations) #location file
write.csv(framelocation, "250147420_16-17_framelocation.csv")

frametime<-data.frame(PUMA_transitions) #timing file
write.csv(frametime, "250147420_16-17_frametime.csv")  
