#### GLS paper Southern Rockhopper penguins - Bahia Franklin - Dodino et al.2024 #####
###  Geolocators analysis - Original code from Lisovski et al. 2019 DOI: 10.1111/1365-2656.13036 (https://geolocationmanual.vogelwarte.ch/probGLS.html)

ls()
rm(list=ls())
ls()

#library(devtools)
#install_github("SLisovski/GeoLocTools")
library(tidyverse)
library(TwGeos)
#remotes::install_github("benjamin-merkel/probGLS", force=TRUE)
library(probGLS)
library(ncdf4)
library(GeoLight)
library(SGAT)
#install.packages ("paletteer")
library (paletteer)

library(devtools)
#install_github("SLisovski/GeoLocTools")
library(GeoLocTools)
setupGeolocation()

######## probGLS ########
wd <- "C:/Users/User/Desktop/GLS projectos"

ID  <- "BV341_04Jan21_180945driftadj"
tag <- "BV341_04Jan21_180945"
Species <- "Rockhopper_franklin"

#load("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper/prob_track/BV317_04Jan21_173451driftadj_probGLS_premolt.RData")


### Deployment
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

map <- ggplot(data=world) +
  geom_sf(fill = gray(0.3), colour = NA)+ #colour = NA removes all borders
  coord_sf(xlim = c(-170,-120), ylim = c(35, 65), expand = FALSE) +
  geom_point(aes(x = -146.326007, y = 59.438033), size = 2, 
             shape = 21, fill = "darkred")

lon.calib <- -64.2500
lat.calib <- -54.7833

## Raw light data / twilight annotation
raw <- readMTlux(paste0(wd,"/RawData/", Species, "/", ID, ".lux")) %>% mutate(Light = log(Light))
head(raw)
# plot(raw[1:3000,])

offset <- 17
lightImage(raw, offset = offset, zlim = c(0,15))
tsimageDeploymentLines(raw$Date, lon.calib, lat.calib, offset = offset, lwd = 1.5, col = "orange")
##############  twilight annotation ###########
#######two windows will appear:
#Move them so they are not on top of each other and you can see both. 
#They should look like a big black blob. This identifies the “nightime” period over time. 
#The top of the blob shows all the sunrises and the bottom of blob shows all the sunsets. 
#You can note for instance that the days get longer (and thus the nights shorter) at the end of the time series, 
#because the blob gets thinner. 
#You may even note changes in the light image that relate to changes in activity patterns or breeding behavior.

###Step 1: Click on “Select subset”. 
#With the left mouse button choose where you want the start of the dataset to be, and right mouse button to choose the end
#Once you are happy with the start and end of the timeseries press “a” on the keyboard to accept and move to next step.

###Step 2: Click on "Find twilights"
#All you need to do here is click in the dark part (in the zoomed in image i.e. the one not entitled “Find twilights”)
#of the image and this will identify all the sunrises (orange) and sunsets (blue) based on the threshold defined in the previous section. 
#Press “a” on the keyboard to accept and move to next step.

#Final step
#Then close the windows with “q”.

twl <- preprocessLight(raw, threshold = 1, lmax = 15, offset = offset)
write.csv(twl, file = glue::glue("C:/Users/User/Desktop/GLS projectos/Results/{Species}/{ID}_twl_2.csv"))

twl <- read.csv(glue::glue("C:/Users/User/Desktop/GLS projectos/Results/{Species}/{ID}_twl.csv")) %>%
  mutate(Twilight = as.POSIXct(Twilight, tz = "GMT")) %>% twilightAdjust(., 5*60)

      
twl <- twilightEdit(twilights = twl, #revisar valores de filtrado
                    offset = offset,
                    window = 4,           # two days before and two days after
                    outlier.mins = 45,    # difference in mins
                    stationary.mins = 25, # are the other surrounding twilights within 25 mins of one another
                    plot = TRUE)

head(twl)
write.csv(twl, file = glue::glue("C:/Users/User/Desktop/GLS projectos/Results/{Species}/{ID}_twl_prueba.csv"))


## Calibration
tw <- twilight_error_estimation(shape = 2.49, scale = 0.94, delay = 0) ## Shoudl be adjusted



## Immersion temperature and tag data
aux <- read.table(glue::glue("C:/Users/User/Desktop/GLS projectos/RawData/{Species}/{tag}.deg"), skip = 20) %>%
  setNames(c("Date", "Time", "Tmin", "Tmax", "Tmean", "wetdry", "wetdry2","wetdry3","wetdry4","Cond")) %>%
  unite(col = "Datetime", c(Date, Time), sep = "", remove = TRUE) %>%
  mutate(Datetime = as.POSIXct(Datetime, format = "%d/%m/%Y %H:%M:%S", tz = "GMT"))%>%
  filter(wetdry >= 2)


## determine daily SST value recorded by the logger
sst_daily <- with(aux %>% dplyr::select("Datetime", "Tmin", "Tmax", "Tmean") %>% pivot_longer(cols = c("Tmin", "Tmax", "Tmean")),
                  sst_deduction(datetime = Datetime, temp = value, temp.range = c(-2,19)))


## act data
act <- aux %>% dplyr::select("Datetime", "wetdry") %>% setNames(c("dtime", "wetdry"))

plot(unique(as.Date(act$dtime)), tapply(act$wetdry,as.Date(act$dtime),sum)/288,
     type="l",col=grey(0.5),ylab="daily proportion of saltwater immersion")
points(unique(as.Date(act$dtime)),tapply(act$wetdry,as.Date(act$dtime),sum)/288,
       pch=19,cex=0.8)

speed_wet  = c(0.27, 0.2, 1.9)# for pre-molt
                             

## ProbGLS
particle.num  = 2000
iteration.num = 200

trn <- export2GeoLight(twl)

## initial track
library(SGAT)
x <- thresholdPath(twl$Twilight, twl$Rise)
plot(x$x)


xlim <- range(x$x[,1]+c(-5,5))
ylim <- range(x$x[,2]+c(-5,5))


#for pre-molt
speed_wet  = c(0.27, 0.2, 1.9)
prob_track  <- prob_algorithm(trn                         = trn,
                              sensor                      = sst_daily[sst_daily$SST.remove==F,],
                              act                         = act,
                              tagging.date                = "2020-01-26",
                              retrieval.date              = "2020-02-22",
                              loess.quartile              = NULL,
                              tagging.location            = c(lon.calib, lat.calib),
                              particle.number             = particle.num,
                              iteration.number            = iteration.num,
                              sunrise.sd                  = tw,
                              sunset.sd                   = tw,
                              range.solar                 = c(-7, -1),
                              boundary.box                = c(-80, -40, -65, -45), # min lon, max lon, min lat and max lat
                              #speed.dry                   = speed_dry,
                              speed.wet                   = speed_wet,
                              sst.sd                      = 0.5,
                              max.sst.diff                = 3,
                              east.west.comp              = T,
                              land.mask                   = T,
                              ice.conc.cutoff             = 0.9,
                              wetdry.resolution           = median(as.numeric(diff(act$dtime))*60),
                              NOAA.OI.location = "C:/Users/User/Desktop/GLS projectos/SST") #downloand from: https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html



## summary

summary(prob_track)

sm <- prob_track[[2]] %>% st_drop_geometry()
plot(sm$dtime, sm$tag.sst)
# abline(v = max_phenotrunc)
points(sm$dtime, sm$mean.sat.sst, type = "o", pch = 16, cex = 0.5, col = "blue")

## Select best track
best_track <- prob_track$`most probable track`

## Transform best track in sf object 
best_track_sf <- st_as_sf(best_track)
best_track_sf <- 
  best_track_sf %>% 
  mutate(lon = st_coordinates(best_track_sf)[, 1],
         lat = st_coordinates(best_track_sf)[, 2])

##Save best track as .RData 
save(best_track_sf, file= paste0(wd, "/Results/", Species, "/Rdata_loop_pre_molt/", ID, "_probGLS_premolt_prueba.RData"))

#for premolt
save(prob_track,file= paste0(wd, "/Results/", Species, "/prob_track/", ID, "_probGLS_premolt_.RData"))

write.csv(best_track_sf,file = glue::glue("C:/Users/User/Desktop/GLS projectos/Results/{Species}/best_track/{ID}_premolt_OP2.csv"))


### Save ALL tracks
load("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/Rdata_loop_pre_molt/GLS_BV316_franklin_premolt.RData")

all_track <- prob_track$`all tracks`

## Transform all track in sf object 
all_track_sf <- st_as_sf(all_track)
all_track_sf <- 
  all_track_sf %>% 
  mutate(lon = st_coordinates(all_track_sf)[, 1],
         lat = st_coordinates(all_track_sf)[, 2])

write.csv(all_track_sf,file = glue::glue("C:/Users/User/Desktop/GLS projectos/Results/{Species}/all_track/BV316_premolt_all_track_2.csv"), row.names = FALSE)


##### ALL TRACKS in one file ######

library(dplyr)


##### merge .csv from best_track in one file #####

BV316 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/best_track/BV316_04Jan21_182652driftadj_premolt_FINAL.csv")
BV318 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/best_track/BV318_04Jan21_185213driftadj_premolt2.csv")
BV319 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/best_track/BV319_04Jan21_182141driftadj_premolt2.csv")
BV328 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/best_track/BV328_04Jan21_171830driftadj_premolt2.csv")
BV330 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/best_track/BV330_04Jan21_175626driftadj_premolt2.csv")
BV332 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/best_track/BV332_04Jan21_180343driftadj_premolt_OP2.csv")
BV334 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/best_track/BV334_04Jan21_183207driftadj_premolt2.csv")
BV335 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/best_track/BV335_04Jan21_173821driftadj_premolt_OP2.csv")
BV337 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/best_track/BV337_04Jan21_174458driftadj_premolt_OP2.csv")
BV338 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/best_track/BV338_04Jan21_175036driftadj_premolt2.csv")
BV339 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/best_track/BV339_04Jan21_183819driftadj_premolt_OP2.csv")
BV341 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/best_track/BV341_04Jan21_180945driftadj_premolt2.csv")
BV342 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/best_track/BV342_04Jan21_184313driftadj_premolt2.csv")
BV344 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/best_track/BV344_04Jan21_185847driftadj_premolt_OP2.csv")


BV316 <- BV316 %>% mutate(id = "BV316")
BV318 <- BV318 %>% mutate(id = "BV318")
BV319 <- BV319 %>% mutate(id = "BV319")
BV328 <- BV328 %>% mutate(id = "BV328")
BV330 <- BV330 %>% mutate(id = "BV330")
BV332 <- BV332 %>% mutate(id = "BV332")
BV334 <- BV334 %>% mutate(id = "BV334")
BV335 <- BV335 %>% mutate(id = "BV335")
BV337 <- BV337 %>% mutate(id = "BV337")
BV338 <- BV338 %>% mutate(id = "BV338")
BV339 <- BV339 %>% mutate(id = "BV339")
BV341 <- BV341 %>% mutate(id = "BV341")
BV342 <- BV337 %>% mutate(id = "BV342")
BV344 <- BV344 %>% mutate(id = "BV344")


BV316 <- BV316[, c("dtime", "lon", "lat", "id")]
BV318 <- BV318[, c("dtime", "lon", "lat", "id")]
BV319 <- BV319[, c("dtime", "lon", "lat", "id")]
BV328 <- BV328[, c("dtime", "lon", "lat", "id")]
BV330 <- BV330[, c("dtime", "lon", "lat", "id")]
BV332 <- BV332[, c("dtime", "lon", "lat", "id")]
BV334 <- BV334[, c("dtime", "lon", "lat", "id")]
BV335 <- BV335[, c("dtime", "lon", "lat", "id")]
BV337 <- BV337[, c("dtime", "lon", "lat", "id")]
BV338 <- BV338[, c("dtime", "lon", "lat", "id")]
BV339 <- BV339[, c("dtime", "lon", "lat", "id")]
BV341 <- BV341[, c("dtime", "lon", "lat", "id")]
BV342 <- BV342[, c("dtime", "lon", "lat", "id")]
BV344 <- BV344[, c("dtime", "lon", "lat", "id")]


merge_tracks <- rbind(BV316, BV318, BV319,
                       BV328, BV330,BV332, 
                       BV334, BV335,BV337,
                       BV338, BV339,BV341,
                       BV342, BV344)


sex_data <- data.frame(
  id = c("BV316", "BV318", "BV319", "BV328", "BV330", "BV332", "BV334", "BV335", "BV337", "BV338", "BV339", "BV341", "BV342", "BV344"),
  sex = c("M", "M", "F", 
          "M", "M","M",
          "M", "F", "M", 
          "F", "F", "M", 
          "F", "M")
)


merge_tracks <- left_join(merge_tracks, sex_data, by = "id")



write.csv(merge_tracks, "merge_tracks_pre-molt_best_track_withSex.csv", row.names = FALSE)

###### merge .csv of all tracks in one file #####

BV316 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/all_track/BV316_premolt_all_track.csv")
BV318 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/all_track/BV318_premolt_all_track.csv")
BV319 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/all_track/BV319_premolt_all_track.csv")
BV328 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/all_track/BV328_premolt_all_track.csv")
BV330 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/all_track/BV330_premolt_all_track.csv")
BV332 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/all_track/BV332_premolt_all_track.csv")
BV334 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/all_track/BV334_premolt_all_track.csv")
BV335 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/all_track/BV335_premolt_all_track.csv")
BV337 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/all_track/BV337_premolt_all_track.csv")
BV338 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/all_track/BV338_premolt_all_track.csv")
BV339 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/all_track/BV339_premolt_all_track.csv")
BV341 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/all_track/BV341_premolt_all_track.csv")
BV342 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/all_track/BV342_premolt_all_track.csv")
BV344 <- read.csv("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/all_track/BV344_premolt_all_track.csv")


library(dplyr)

BV316 <- BV316 %>% mutate(id = "BV316")
BV318 <- BV318 %>% mutate(id = "BV318")
BV319 <- BV319 %>% mutate(id = "BV319")
BV328 <- BV328 %>% mutate(id = "BV328")
BV330 <- BV330 %>% mutate(id = "BV330")
BV332 <- BV332 %>% mutate(id = "BV332")
BV334 <- BV334 %>% mutate(id = "BV334")
BV335 <- BV335 %>% mutate(id = "BV335")
BV337 <- BV337 %>% mutate(id = "BV337")
BV338 <- BV338 %>% mutate(id = "BV338")
BV339 <- BV339 %>% mutate(id = "BV339")
BV341 <- BV341 %>% mutate(id = "BV341")
BV342 <- BV337 %>% mutate(id = "BV342")
BV344 <- BV344 %>% mutate(id = "BV344")


BV316 <- BV316[, c("type", "lon", "lat", "id")] 
BV318 <- BV318[, c("type", "lon", "lat", "id")] 
BV319 <- BV319[, c("type", "lon", "lat", "id")]
BV328 <- BV328[, c("type", "lon", "lat", "id")]
BV330 <- BV330[, c("type", "lon", "lat", "id")]
BV332 <- BV332[, c("type", "lon", "lat", "id")]
BV334 <- BV334[, c("type", "lon", "lat", "id")]
BV335 <- BV335[, c("type", "lon", "lat", "id")]
BV337 <- BV337[, c("type", "lon", "lat", "id")]
BV338 <- BV338[, c("type", "lon", "lat", "id")]
BV339 <- BV339[, c("type", "lon", "lat", "id")]
BV341 <- BV341[, c("type", "lon", "lat", "id")]
BV342 <- BV342[, c("type", "lon", "lat", "id")]
BV344 <- BV344[, c("type", "lon", "lat", "id")]


merge_tracks <- rbind(BV316, BV318, BV319,
                       BV328, BV330,BV332, 
                       BV334, BV335,BV337,
                       BV338, BV339,BV341,
                       BV342, BV344)


sex_data <- data.frame(
  id = c("BV316", "BV318", "BV319", "BV328", "BV330", "BV332", "BV334", "BV335", "BV337", "BV338", "BV339", "BV341", "BV342", "BV344"),
  sex = c("M", "M", "F", 
          "M", "M","M",
          "M", "F", "M", 
          "F", "F", "M", 
          "F", "M")
)


merge_tracks <- left_join(merge_tracks, sex_data, by = "id")

write.csv(merge_tracks, "merge_tracks_pre-molt_all_track.csv", row.names = FALSE)


