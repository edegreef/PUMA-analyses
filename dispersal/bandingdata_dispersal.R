#Script for calculating and mapping dispersal distances using bird banding encounter records. It's not a very elegant script, but contains steps for when I examined juvenile dispersal. 
#Oct 2019- Evelien de Greef

setwd("C:/Users/Evelien de Greef/Dropbox/PUMA/Encounter data")

#------------------PREPPING AND FILTERING DATA--------------------

#Load data file
PUMAraw <- read.csv("./PUMAbands.csv") #Data obtained from NABBP, here I had already filtered the file to only include juveniles, locals, and HY birds

library(tidyverse)

#Renaming some of the column names
PUMA <- select(PUMAraw,
               Band = Band.Number,
               B_AOU = B.AOU, #PUMA code
               B_permit = B.Permit,
               S,
               B_date = B.Date,
               B_10minBLk = B.10.min.Blk,
               B_lat = B.Lat,
               B_long = B.Lon,
               B_CP = B.CP, #coordinate precision
               B_D = B.D, #direction code
               B_reg = B.Reg, #region
               B_status = B.Bird.Status,
               B_age = B.Age,
               B_sex = B.Sex,
               E_date = E.Date,
               E_10minBlk = E.10.min.Blk,
               E_lat = E.Lat,
               E_long = E.Lon,
               E_CP = E.CP, #coordinate precision
               E_D = E.D, #direction code
               E_reg = E.Reg, #region
               E_how = E.How, #how recaptured
               E_PC = E.PC, #present condition
               E_who = E.Who,
               E_repmthd = E.Rep.Mthd,
               HSS, #hunting season survived/ approx age of bird at encounter
               Other_bands = Other.Bands)

#Converting the B (banding) and E (encounter) dates into date format, this will turn all the inexact dates into NAs
#The "inexact" dates are in Bandit code format, but I only want to extract month and year in the end
PUMA$B_date <- as.Date(PUMA$B_date, format="%Y-%m-%d")
PUMA$E_date <- as.Date(PUMA$E_date, format="%Y-%m-%d")

#Making a separate dataframe of the E (encounter) dates that have NA - will use this to fill in the month and/or year later
PUMA_inexact_Edate <- subset(PUMA, is.na(PUMA$E_date))

#Removing NAs from main dataframe so it will only contain exact dates and removed the 4 NAs with no lat/long for E location
PUMA <- na.omit(PUMA)

#Creating new columns with year and month separated for B_dates and E_dates
PUMA$B_month <- as.numeric(format(PUMA$B_date, "%m"))
PUMA$B_year <- as.numeric(format(PUMA$B_date, "%Y"))
PUMA$E_month <- as.numeric(format(PUMA$E_date, "%m"))
PUMA$E_year <- as.numeric(format(PUMA$E_date, "%Y"))

#Sswitching back to factor for E_date column to able to merge with the subsetted data later
PUMA$E_date <- as.factor(format(PUMA$E_date))

#Going back to the dataframe with inexact encounter dates and creating list of just band numbers for individuals with inexact dates.
PUMA_inexact_Edate_bands <- as.data.frame(PUMA_inexact_Edate[,1])
PUMA_inexact_Edate_bands$`PUMA_inexact_Edate[, 1]` <- as.character(PUMA_inexact_Edate_bands$`PUMA_inexact_Edate[, 1]`)

#If band matches raw data, then enter in the E_date (from PUMAraw to PUMA_inexact_Edate)
PUMAraw$Band.Number <- as.character(PUMAraw$Band.Number)
PUMA_inexact_Edate_extract <- left_join(PUMA_inexact_Edate_bands, PUMAraw, by = c("PUMA_inexact_Edate[, 1]" = "Band.Number"))

#Extracting the PUMA_inexact_Edate_extract as a .csv to for manual adjustment (for ones we can salvage month and year)
PUMA_dates_subset <-data.frame(PUMA_inexact_Edate_extract)
write.csv(PUMA_dates_subset, "PUMA_dates_subset.csv")  #save for backup

#Removing rows that have NA in the E_month column
PUMA_dates_subset <- PUMA_dates_subset[!is.na(PUMA_dates_subset$E_month),]

PUMA_dates_subset$B_date <- as.Date(PUMA_dates_subset$B_date, format="%Y-%m-%d")
PUMA_dates_subset$B_month <- as.numeric(format(PUMA_dates_subset$B_date, "%m"))
PUMA_dates_subset$B_year <- as.numeric(format(PUMA_dates_subset$B_date, "%Y"))

#Combining PUMA and PUMA_dates_subset together
PUMA <- rbind(PUMA, PUMA_dates_subset)


#Extracting the individuals with multiple/duplicate band encounter records
PUMA_dup <- PUMA[duplicated(PUMA$Band)|duplicated(PUMA$Band, fromLast=TRUE),]
write.csv(PUMA_dup, "PUMA_dup.csv") #saving as backup

#I manually edited the PUMA_dup.csv file to include a column with "Y" or "N" for duplicate or not. It was not always the first record I wanted to keep- it depended on the enounter. I wanted to examine dispersal from natal site to following or first-known breeding year, while excluding the extra encounter records within the same individual.

PUMA_dup_noted <- read.csv("./PUMA_dup.csv") 

#Editing main dataset to remove all duplicates, and then adding in the individuals back in but with only 1 record each
PUMA <- PUMA[! PUMA$Band %in% unique(PUMA[duplicated(PUMA$Band), "Band"]), ] #Taking out all band# that are in more than once
PUMA$Dup <- NA #creating empty column
PUMA_dups_noted <- filter(PUMA_dups_noted, Dup == "N") #Keeping only the band records I want from the duplicated records
PUMA <- rbind(PUMA, PUMA_dups_noted)
length(unique(PUMA$Band)) #double checking number of records to make sure it adds up


#Filtering out encounters in same year as banding, and filtering months timespan to avoid using banding or encounters on migration
PUMA <- subset(PUMA, PUMA$B_year!=PUMA$E_year)
PUMA_filtered <- filter(PUMA, B_month >= 5 & B_month <= 7, E_month >= 5 & E_month <=7) #May, June, July

#Removing one bird encounted in Utah, as this one will be considered an outlier (was found dead substantially out of range):
PUMA_filtered_noUTAH <- PUMA_filtered[-c(463),]
write.csv(PUMA_filtered_noUTAH, "PUMA_bands_2060.csv") #saving as backup


#------------------MEASURING DISPERSAL DISTANCES--------------------
library(geosphere)

PUMA_bands_ready <- PUMA_filtered_noUTAH
PUMA_bands_ready <- select(PUMA_bands_ready, -X) #only necessary if loading data back in and there is an extra "X" column with row number


#Test run before running loop. Distm is in longitude,latitude format. Here, column 8 = longitude of banding, 7 = latitude of banding, 18 = longitude of encounter, 17 = latitude of encounter.
dist_test <- distm(c(PUMA_bands_ready[1,8],PUMA_bands_ready[1,7]), c(PUMA_bands_ready[1,18],PUMA_bands_ready[1,17])) #dist_test should be 176852m for this dataset

#Creating empty matrix with 1 column, which will contain the dispersal distances
dist_dispersal <- matrix(nrow = nrow(PUMA_bands_ready),ncol = 1)

#For loop to run distm for each row in PUMA. Note that distm outputs results in meters
for (row in 1:nrow(PUMA_bands_ready)){
  dist_dispersal[row,1] <- distm(c(PUMA_bands_ready[row,8],PUMA_bands_ready[row,7]), c(PUMA_bands_ready[row,18],PUMA_bands_ready[row,17]))[1,1]
}

dist_dispersal <- round(dist_dispersal, digits = 0) #rounding distances to whole number

mean(dist_dispersal) #calculating average dispersal distance (m)
range(dist_dispersal) #calculating range
hist(dist_dispersal, breaks=200) #histogram

#Playing around with quantifying dispersal across all the individuals. Here I am counting number of individuals with a specific dispersal distance range.
sum(dist_dispersal < 1000) #Cunting number of individuals who dispersed under 1000 meters)
sum(dist_dispersal >= 0 & dist_dispersal <= 19000)
sum(dist_dispersal >= 19000 & dist_dispersal <= 38000)
sum(dist_dispersal >= 38000 & dist_dispersal <= 57000)
sum(dist_dispersal >= 57000 & dist_dispersal <= 76000)
sum(dist_dispersal >= 76000 & dist_dispersal <= 95000)
sum(dist_dispersal >= 95000 & dist_dispersal <= 114000)
sum(dist_dispersal >= 114000 & dist_dispersal <= 133000)
sum(dist_dispersal >= 133000 & dist_dispersal <= 152000)
sum(dist_dispersal >= 133000 & dist_dispersal <= 152000)
sum(dist_dispersal >= 152000 & dist_dispersal <= 171000)
sum(dist_dispersal >= 171000 & dist_dispersal <= 190000)
sum(dist_dispersal >= 190000 & dist_dispersal <= 209000)
sum(dist_dispersal >= 209000 & dist_dispersal <= 228000)
sum(dist_dispersal >= 228000 & dist_dispersal <= 247000)
sum(dist_dispersal >= 247000 & dist_dispersal <= 266000)
sum(dist_dispersal >= 266000 & dist_dispersal <= 285000)
sum(dist_dispersal >= 266000 & dist_dispersal <= 285000)
sum(dist_dispersal >= 285000 & dist_dispersal <= 304000)
sum(dist_dispersal >= 304000 & dist_dispersal <= 323000)
sum(dist_dispersal >= 323000 & dist_dispersal <= 342000)
sum(dist_dispersal >= 342000 & dist_dispersal <= 361000)
sum(dist_dispersal >= 361000 & dist_dispersal <= 380000)
sum(dist_dispersal >= 380000 & dist_dispersal <= 399000)
sum(dist_dispersal >= 399000 & dist_dispersal <= 418000)
sum(dist_dispersal >= 418000 & dist_dispersal <= 437000)
sum(dist_dispersal >= 437000 & dist_dispersal <= 456000)
sum(dist_dispersal >= 456000 & dist_dispersal <= 475000)
sum(dist_dispersal >= 475000 & dist_dispersal <= 494000)
sum(dist_dispersal >= 494000 & dist_dispersal <= 513000)
sum(dist_dispersal >= 513000 & dist_dispersal <= 532000)
sum(dist_dispersal >= 532000 & dist_dispersal <= 1000000)
sum(dist_dispersal >= 1000000 & dist_dispersal <= 1500000)
sum(dist_dispersal >= 1500000 & dist_dispersal <= 2000000)
sum(dist_dispersal >= 2000000 & dist_dispersal <= 2500000)
sum(dist_dispersal >= 2500000 & dist_dispersal <= 3000000)

#Dividing dist_dispersals into 19km blocks
dist_dispersal_10min <- dist_dispersal/19000
dist_dispersal_far <- subset(dist_dispersal >= 19000)

#Using larger chunks/intervals:
sum(dist_dispersal >= 0 & dist_dispersal <= 20000)
sum(dist_dispersal >= 20000 & dist_dispersal <= 50000)
sum(dist_dispersal >= 50000 & dist_dispersal <= 100000)
sum(dist_dispersal >= 100000 & dist_dispersal <= 200000)
sum(dist_dispersal >= 200000 & dist_dispersal <= 500000)
sum(dist_dispersal >= 500000 & dist_dispersal <= 1000000)
sum(dist_dispersal > 1000000)

#Sorting data by order of dispersal (to ID the large ones easily, at the bottom of the list)
PUMA_with_dist <- cbind(PUMA_bands_ready, dist_dispersal)
PUMA_with_dist <- PUMA_with_dist[order(dist_dispersal),]

top_dispersers <- tail(PUMA_with_dist, 34)
top_7 <- tail(PUMA_with_dist, 7)

#------------------MAPPING DISPERSAL--------------------

#https://www.r-spatial.org/r/2018/10/25/ggplot2-sf.html
library(rnaturalearth)
library(ggspatial)

canada <- ne_states(country="canada", returnclass="sf")
usa <- ne_states(country="united states of america", returnclass="sf")
mexico <- ne_states(country="mexico", returnclass="sf")
#countries <-ne_countries(returnclass = "sf")

library(ggplot2)
countryborders <- map_data("world", c("usa", "Canada", "Mexico"))

#Mapping initial banding (red) and encounter (blue) locations for each record, and drawing line between them.
png("PUMA_bands2060_1000x1200_p6.png", height=1000, width=1200)
ggplot() +
  geom_sf(data=canada, color="gray90", size=1, fill="white") +
  geom_sf(data=usa, color="gray90", size=1, fill="white") +
  geom_sf(data=mexico, color="gray90",size=1, fill="white") +
  geom_polygon(data=countryborders, aes(x=long, y=lat, group=group), color="black", size=1, fill=NA)+
  annotation_scale(location="bl") +
  geom_point(data=PUMA_bands_ready, aes(x=E_long, y=E_lat), size=6, shape=20, col="blue", alpha=0.5) +
  geom_point(data=PUMA_bands_ready, aes(x=B_long, y=B_lat), size=6, shape=20, col="red", alpha=0.5) +
  coord_sf(xlim = c(-130, -60), ylim = c(20, 60), expand = FALSE) +
  geom_segment(data=PUMA_bands_ready, mapping=aes(x=E_long, y=E_lat, xend=B_long, yend=B_lat), size=1, color="orange", alpha=0.5)+
  theme(panel.background=element_rect(fill="#E1F1FF")) +
  labs(title="PUMA juvenile dispersal", y="Latitude", x="Longitude", family="arial") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


#Mapping the top 34 dispersers
ggplot() +
  geom_sf(data=canada, color="gray90", size=1, fill="white") +
  geom_sf(data=usa, color="gray90", size=1, fill="white") +
  geom_sf(data=mexico, color="gray90",size=1, fill="white") +
  geom_polygon(data=countryborders, aes(x=long, y=lat, group=group), color="black", size=1, fill=NA)+
  annotation_scale(location="bl") +
  geom_point(data=top_dispersers, aes(x=E_long, y=E_lat), size=5, shape=20, col="blue", alpha=0.5) +
  geom_point(data=top_dispersers, aes(x=B_long, y=B_lat), size=5, shape=20, col="red", alpha=0.5) +
  coord_sf(xlim = c(-130, -60), ylim = c(20, 60), expand = FALSE) +
  geom_segment(data=top_dispersers, mapping=aes(x=E_long, y=E_lat, xend=B_long, yend=B_lat), size=1, color="orange", alpha=0.5)+
  theme(panel.grid.major=element_line(color=gray(.5), linetype="dashed", size=0.5), panel.background=element_rect(fill="lightskyblue1")) +
  labs(title="PUMA top dispersers", y="Latitude", x="Longitude", family="arial") +
  theme(plot.title = element_text(hjust = 0.5))

#Mapping the top 7 dispersers:
ggplot() +
  geom_sf(data=canada, color="gray90", size=1, fill="white") +
  geom_sf(data=usa, color="gray90", size=1, fill="white") +
  geom_sf(data=mexico, color="gray90",size=1, fill="white") +
  geom_polygon(data=countryborders, aes(x=long, y=lat, group=group), color="black", size=1, fill=NA)+
  annotation_scale(location="bl") +
  geom_point(data=top_7, aes(x=E_long, y=E_lat), size=5, shape=20, col="blue", alpha=0.5) +
  geom_point(data=top_7, aes(x=B_long, y=B_lat), size=5, shape=20, col="red", alpha=0.5) +
  coord_sf(xlim = c(-130, -60), ylim = c(20, 60), expand = FALSE) +
  geom_segment(data=top_7, mapping=aes(x=E_long, y=E_lat, xend=B_long, yend=B_lat), size=1, color="orange", alpha=0.5)+
  theme(panel.grid.major=element_line(color=gray(.5), linetype="dashed", size=0.5), panel.background=element_rect(fill="lightskyblue1")) +
  labs(title="PUMA top 7 dispersers", y="Latitude", x="Longitude", family="arial") +
  theme(plot.title = element_text(hjust = 0.5))


#---------------EXAMINING LATITUDE VS LONGITUDE PATTERN-----------------
library(geosphere)

PUMA_bands_2060 <- PUMA_bands_ready

#Creating empty matrix with 1 column to fill in later with dispersal distance on the north-south direction
dist_dispersal_south <- matrix(nrow = nrow(PUMA_bands_2060),ncol = 1)

#Caluculating dispersal distance on north-south axis. Measuring from Long and Lat of initial banding location, and Long of initial banding with Lat of encounter location.
for (row in 1:nrow(PUMA_bands_2060)){
  dist_dispersal_south[row,1] <- distm(c(PUMA_bands_2060[row,8],PUMA_bands_2060[row,7]), c(PUMA_bands_2060[row,8],PUMA_bands_2060[row,17]))[1,1]
}
dist_dispersal_south <- round(dist_dispersal_south, digits = 0) #rounding distances to whole number

#Same for movements west-east
dist_dispersal_east <- matrix(nrow = nrow(PUMA_bands_2060),ncol = 1)

#Calculating dispersal distance on west-east axis. Measuring from Long of initial banding with Lat of encounter, and Long and Lat of encounter location.
for (row in 1:nrow(PUMA_bands_2060)){
  dist_dispersal_east[row,1] <- distm(c(PUMA_bands_2060[row,8],PUMA_bands_2060[row,17]), c(PUMA_bands_2060[row,18],PUMA_bands_2060[row,17]))[1,1]
}
dist_dispersal_east <- round(dist_dispersal_east, digits = 0) #rounding distances to whole number

#Merging the info together with the main dataframe
updated <- cbind(PUMA_bands_2060, dist_dispersal, dist_dispersal_south, dist_dispersal_east)

updated$ratio <- (updated$dist_dispersal_east / updated$dist_dispersal_south)
#ones with ratio > 1 means LONGITUTDINAL
#ones with ratio < 1 means LATITUDINAL

#Sorting by ratio and spot checking lat/long for band and encounter, making sure it worked properly
sorted <- updated[order(updated$ratio),]

#Filtering out ones that did not disperse
library(tidyverse)
dispersers_only <- subset(updated, updated$dist_dispersal != 0)

ratio_more <- subset(dispersers_only, dispersers_only$ratio > 1)
ratio_less <- subset(dispersers_only, dispersers_only$ratio < 1)

write.csv(ratio_more, "longitudinal_dispersers.csv")
write.csv(ratio_less, "latitudinal_dispersers.csv")

#if need to read back in
long <- read.csv("longitudinal_dispersers.csv", header=T)
lat <- read.csv("latitudinal_dispersers.csv", header=T)

long <- ratio_more
lat <- ratio_less

#Plotting map of latitudinal and longitudinal dispersal

#Same mapping set-up as before
library(rnaturalearth)
library(ggspatial)

canada <- ne_states(country="canada", returnclass="sf")
usa <- ne_states(country="united states of america", returnclass="sf")
mexico <- ne_states(country="mexico", returnclass="sf")
#countries <-ne_countries(returnclass = "sf")

library(ggplot2)
countryborders <- map_data("world", c("usa", "Canada", "Mexico"))

#latitudinal dispersers only, can do the same for longitudinal too using other dataframe
png("latitudinal dispersers.png", height=700, width=1400)
ggplot() +
  geom_sf(data=canada, color="gray90", size=1, fill="white") +
  geom_sf(data=usa, color="gray90", size=1, fill="white") +
  geom_sf(data=mexico, color="gray90",size=1, fill="white") +
  geom_polygon(data=countryborders, aes(x=long, y=lat, group=group), color="black", size=1, fill=NA)+
  annotation_scale(location="bl") +
  geom_point(data=lat, aes(x=E_long, y=E_lat), size=5, shape=20, col="blue", alpha=0.5) +
  geom_point(data=lat, aes(x=B_long, y=B_lat), size=5, shape=20, col="red", alpha=0.5) +
  coord_sf(xlim = c(-130, -60), ylim = c(20, 60), expand = FALSE) +
  geom_segment(data=lat, mapping=aes(x=E_long, y=E_lat, xend=B_long, yend=B_lat), size=1, color="orange", alpha=0.5)+
  theme(panel.grid.major=element_line(color=gray(.5), linetype="dashed", size=0.5), panel.background=element_rect(fill="lightskyblue1")) +
  labs(title="latitudinal dispersers", y="Latitude", x="Longitude", family="arial") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#sorting by overal dispersal distance
sorted2 <- dispersers_only[order(dispersers_only$dist_dispersal),]


#Trying something out for just the longer distance ones
long$dist_dispersal <- (long$dist_dispersal / 1000)
lat$dist_dispersal <- (lat$dist_dispersal / 1000)

hist(long$dist_dispersal, breaks=seq(0,3000, 10))
hist(lat$dist_dispersal, breaks=seq(0,3000, 10))

long$type <- "Long"
lat$type <- "Lat"

filter_lat <- subset(lat, dist_dispersal >=100)
filter_long <- subset(long, dist_dispersal >=100)
hist(filter_long$dist_dispersal, breaks=seq(0,3000, 10))
hist(filter_lat$dist_dispersal, breaks=seq(0,3000, 10))

longandlat <- rbind(long, lat)
filter_500km <- subset(longandlat, dist_dispersal >= 500)
filter_lat500 <- subset(lat, dist_dispersal >=500)
filter_long500 <- subset(long, dist_dispersal >=500)

#Mapping all dispersal records, highlight latitudinal and longitudinal distances in separate colors.
png("both dispersers.png", height=700, width=1400)
ggplot() +
  geom_sf(data=canada, color="gray90", size=1, fill="white") +
  geom_sf(data=usa, color="gray90", size=1, fill="white") +
  geom_sf(data=mexico, color="gray90",size=1, fill="white") +
  geom_polygon(data=countryborders, aes(x=long, y=lat, group=group), color="black", size=1, fill=NA)+
  annotation_scale(location="bl") +
  geom_point(data=longandlat, aes(x=E_long, y=E_lat), size=5, shape=20, col="blue", alpha=0.5) +
  geom_point(data=longandlat, aes(x=B_long, y=B_lat), size=5, shape=20, col="red", alpha=0.5) +
  coord_sf(xlim = c(-130, -60), ylim = c(20, 60), expand = FALSE) +
  geom_segment(data=long, mapping=aes(x=E_long, y=E_lat, xend=B_long, yend=B_lat), size=1, color="orange", alpha=0.5)+
  geom_segment(data=lat, mapping=aes(x=E_long, y=E_lat, xend=B_long, yend=B_lat), size=1, color="green", alpha=0.5)+
  theme(panel.grid.major=element_line(color=gray(.5), linetype="dashed", size=0.5), panel.background=element_rect(fill="lightskyblue1")) +
  labs(title="both lat and long dispersers", y="Latitude", x="Longitude", family="arial") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


#Mapping only dispersals over 500km (top dispersers) with latitudinal and longitudinal colors.
png("PUMA_bands2060_500KM_lat_and_long_1000x1200_p6.png", height=1000, width=1200)
ggplot() +
  geom_sf(data=canada, color="gray90", size=1, fill="white") +
  geom_sf(data=usa, color="gray90", size=1, fill="white") +
  geom_sf(data=mexico, color="gray90",size=1, fill="white") +
  geom_polygon(data=countryborders, aes(x=long, y=lat, group=group), color="black", size=1, fill=NA)+
  annotation_scale(location="bl") +
  geom_point(data=filter_500km, aes(x=E_long, y=E_lat), size=6, shape=20, col="blue", alpha=0.7) +
  geom_point(data=filter_500km, aes(x=B_long, y=B_lat), size=6, shape=20, col="red", alpha=0.7) +
  coord_sf(xlim = c(-130, -60), ylim = c(20, 60), expand = FALSE) +
  geom_segment(data=filter_long500, mapping=aes(x=E_long, y=E_lat, xend=B_long, yend=B_lat), size=1.3, color="orange", alpha=0.7)+
  geom_segment(data=filter_lat500, mapping=aes(x=E_long, y=E_lat, xend=B_long, yend=B_lat), size=1.3, color="green", alpha=0.7)+
  theme(panel.background=element_rect(fill="#E1F1FF")) +
  labs(title="lat and long dispersers top 15%", y="Latitude", x="Longitude", family="arial") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

