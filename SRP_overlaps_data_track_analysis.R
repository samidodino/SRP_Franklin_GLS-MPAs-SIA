##### Data track analysis and MPAs/Polar Front overlaps ####
#### GLS paper Southern Rockhopper penguins - Bahia Franklin - Dodino et al.2024 #####

library(dplyr)
library(nlme)
library(lme4)
library(MuMIn)

####latitudinal range --> using best track ###

rm(list = ls())
ls()

tabla<-read.csv("C:/Users/User/Desktop/GLS projectos/PPA_Franklin/Results/datos_agrupados_pre-molt_withSex.csv",header=T)
names(tabla)
str(tabla)

latitudinal_range <- tabla %>%
  group_by(id) %>%
  summarize(northernmost_lat = max(lat),
            southernmost_lat = min(lat)) %>%
  mutate(latitudinal_range = northernmost_lat - southernmost_lat)


tabla <- tabla %>%
  left_join(latitudinal_range, by = "id") %>%
  mutate(latitudinal_range_km = latitudinal_range * (40075 / 360))

resumen_latitudinal <- tabla %>%
  group_by(id) %>%
  summarize(mean(latitudinal_range_km))

View(resumen_latitudinal)


#### Differences between sexes in duration (days) / distance (km) in dispersion --->best track ####

tabla2<-read.csv("C:/Users/User/Desktop/GLS projectos/PPA_Franklin/Results/ppa_gls_data_track_analysis.csv",header=T)

M.lm.1<-lm(dias~sex,data=tabla2, na.action=na.fail)
plot(M.lm.1) #qqplot normality
M.lm.2<-lm(km~sex,data=tabla2, na.action=na.fail)
plot(M.lm.2) #qqplot normality


vf2 <- varIdent(form= ~ 1 | sex)

M.gls1<-gls(dias~sex,weights = vf2, data=tabla2, method="ML", na.action=na.fail)
summary(M.gls1)
anova(M.gls1, test = "F")

M.gls2<-gls(km~sex,weights = vf2, data=tabla2, method="ML", na.action=na.fail)
summary(M.gls2)

anova(M.gls2, test = "F")

###### Estimation of overlap with frontal zone --> all tracks ######

library(sf)
library(dplyr)
library(tmap)
library(tmaptools)
library(data.table)

fronts <- st_read("C:/Users/User/Desktop/GLS projectos/Fronts_of_the_Antarctic_Circumpolar_Current_-_GIS_data/shapefile/antarctic_circumpolar_current_fronts.shp")
str(frentes2$NAME)

frontal_zone<-dplyr::filter(fronts, NAME %in% c("Polar Front (PF)"))

# cast to polygons 
polygon_FZ <- frente_polar %>%  
  st_cast("POLYGON") 

FZ_sf <- st_as_sf(polygon_FZ)

long1<- -76.35
lat1<-  -61.72

long2<- -72.8
lat2<- -65.1

long3<- -50
lat3<-  -61.72

long4<- -50
lat4<-  -65.1

latlong<- data.table(long= c(long1, long2, long3, long4), lat= c(lat1, lat2, lat3, lat4))

latlong_sf <- st_as_sf(latlong, coords = c("long", "lat"), crs = 4326)

polygon_complement<- st_convex_hull(st_union(st_geometry(latlong_sf)))


#### points as a spatial object --> all tracks

tabla2<-read.csv("C:/Users/User/Desktop/GLS projectos/PPA_Franklin/Results/overlaps_with_all_track/datos_agrupados_pre-molt_all_track.csv",header=T)

# add a column called d2 with a continuous sequence of numbers, so that each point has a unique ID and facilitates subsequent analysis.
tabla2 <- tabla2 %>%
  mutate(id2 = row_number())

write.csv(tabla2, "C:/Users/User/Desktop/GLS projectos/PPA_Franklin/Results/overlaps_with_all_track/datos_agrupados_pre-molt_id2.csv", row.names = FALSE)

tabla2<-read.csv("C:/Users/User/Desktop/GLS projectos/PPA_Franklin/Results/overlaps_with_all_track/datos_agrupados_pre-molt_id2.csv",header=T)

#Dataframe to spatial object

ppa_locations<- st_as_sf(tabla2, coords = c("lon", "lat"), crs=4326)

# Test polygon dimension

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

map <- ggplot(data=world) +
  geom_sf(fill = gray(0.3), colour = NA)+ #colour = NA removes all borders
  coord_sf(xlim = c(-170,-120), ylim = c(35, 65), expand = FALSE) +
  geom_point(aes(x = -146.326007, y = 59.438033), size = 2, 
             shape = 21, fill = "darkred")

lon.calib <- -64.2500
lat.calib <- -54.7833


ggplot(data = world) +
  geom_sf(fill = gray(0.3), colour = NA) +
  geom_sf(data = polygon_complement,col=NA, fill = "red", alpha = 0.5) +
  geom_sf(data = FP_sf,col=NA, fill = "red", alpha = 0.5) +
  geom_sf(data = ppa_locations)+
  geom_point(aes(x = lon.calib, y = lat.calib), size = 3, 
             shape = 21, fill = "darkred") +
  geom_sf(data = frentes2 %>% filter(NAME %in% c("Polar Front (PF)")), aes(color = "blue"), size = 3) +
  coord_sf(xlim = c(-90,-40), ylim = c(-70, -45), expand = FALSE) +
  theme_light() +
  scale_x_continuous(name="") +
  scale_y_continuous(name="")  


# Calculate the percentage of points within the polygon.

intersection_polygon<-st_intersection(ppa_locations,polygon_complement)

percentage_inside <- nrow(intersection_polygon) / nrow(ppa_locations) * 100

intersection_FZ<-st_intersection(ppa_locations, FZ_sf)

percentage_inside2 <- nrow(intersection_FZ) / nrow(ppa_locations) * 100

sum(percentage_inside,percentage_inside2) # 61.3% of points around the frontal zone with best_tracks
                                          # 61.8% of points around the frontal zone with all_tracks


## assign zeros and ones to the days that fall beyond the frontal zone and to those that do not ###


## intersections to data.frame

intersection_FP_df <- as.data.frame(intersection_FP)


intersection_polygon_df <- as.data.frame(intersection_polygon)


# Create a logical vector that indicates whether id2 is present in the smaller DataFrames
present_in_small1 <- tabla2$id2 %in% intersection_FP_df$id2
present_in_small2 <- tabla2$id2 %in% intersection_polygon_df$id2

# Use ifelse to create the new column in original df
tabla2$intersection_FP <- ifelse(present_in_small1 | present_in_small2, 1, 0)

# Save
write.csv(tabla2, file = "data_merge_pre-molt_plusInterac_FZ_all_tracks.csv", row.names = FALSE)


### Number of days at the frontal zone for each individual #######

data<-read.csv("C:/Users/User/Desktop/GLS projectos/PPA_Franklin/Results/overlaps_with_all_track/datos_agrupados_pre-molt_plusInterac_PF_all_tracks.csv",header=T)

data$dtime <- as.Date(data$dtime, format = "%d/%m/%Y")

# Calculate date range for each id where intersection equals 1
results <- data %>%
  filter(intersection_FP == 1) %>%
  group_by(id) %>%
  summarize(from = min(dtime), to = max(dtime))%>%  
  mutate(days_at_frotal_zone = as.numeric(to - from, units = "days"))

# Calculate the total days for each id
total_days <- data %>%
  group_by(id) %>%
  summarise(total_days_trip = n_distinct(dtime))

# Unir los resultados anteriores y calcular el porcentaje de dias en el frente polar
results <- results %>%
  left_join(total_days, by = "id") %>%
  mutate(percent_of_thetrip = (days_at_frontal_zone / total_days_trip) * 100)

# Mostrar los resultados finales
print(results)


write.csv(results, file = "days_at_frontal_zone_all_tracks.csv", row.names = FALSE)


###### Estimation of overlap between Marine Protected Areas and tracks ######

library(sf)

##Shapefiles with marine protected areas

yaganesI<-st_read("C:/Users/User/Desktop/GLS projectos/MarineAreas/Yaganes I.shp")
yaganesII<-st_read("C:/Users/User/Desktop/GLS projectos/MarineAreas/Yaganes II.shp")

burdwood<-st_read("C:/Users/User/Desktop/GLS projectos/MarineAreas/BB_I-II/BB_I_II_limite_zonificacion.shp")

diego_ramirez<-st_read("C:/Users/User/Desktop/GLS projectos/MarineAreas/Diego_Ramirez/WDPA_WDOECM_Sep2023_Public_555643508_shp-polygons.shp")


# points as a spatial object
# using best track

# tabla1<-read.csv("C:/Users/User/Desktop/GLS projectos/PPA_Franklin/Results/datos_agrupados_pre-molt_withSex.csv",header=T)

# using all tracks
tabla2<-read.csv("C:/Users/User/Desktop/GLS projectos/PPA_Franklin/Results/overlaps_with_all_track/datos_agrupados_pre-molt_id2.csv",header=T)

ppa_locations<- st_as_sf(tabla2, coords = c("lon", "lat"), crs=4326)

# Percent of points inside each marine area.

intersection_yaganesI<-st_intersection(ppa_locations,yaganesI)

percentage_yaganesI <- nrow(intersection_yaganesI) / nrow(ppa_locations) * 100

intersection_yaganesII<-st_intersection(ppa_locations, yaganesII)

percentage_yaganesII <- nrow(intersection_yaganesII) / nrow(ppa_locations) * 100

sum(percentage_yaganesI,percentage_yaganesII) # 10.3% using best tracks
                                              # 7.3 % using all tracks 

intersection_burdwood<-st_intersection(ppa_locations,burdwood)

percentage_burdwood <- nrow(intersection_burdwood) / nrow(ppa_locations) * 100 # 3.4% using best tracks
                                                                               # 2.9% using all tracks
intersection_diego_ramirez<-st_intersection(ppa_locations, diego_ramirez)

percentage_diego_ramirez <- nrow(intersection_diego_ramirez) / nrow(ppa_locations) * 100 # 2.4% using best tracks
                                                                                         # 4.0% using all tracks
#percent total i all the marine areas:

sum(percentage_yaganesI,percentage_yaganesII,percentage_burdwood,percentage_diego_ramirez) # 16% using best tracks
                                                                                           # 14.2% using all tracks


intersection_yaganesI_df <- as.data.frame(intersection_yaganesI)
intersection_yaganesII_df <- as.data.frame(intersection_yaganesII)
intersection_burdwood_df <- as.data.frame(intersection_burdwood)
intersection_diego_ramirez_df <- as.data.frame(intersection_diego_ramirez)

present_in_small3 <- tabla2$id2 %in% intersection_yaganesI_df$id2
present_in_small4 <- tabla2$id2 %in% intersection_yaganesII_df$id2
present_in_small5 <- tabla2$id2 %in% intersection_burdwood_df$id2
present_in_small6 <- tabla2$id2 %in% intersection_diego_ramirez_df$id2


tabla2$intersection_Yaganes <- ifelse(present_in_small3 | present_in_small4, 1, 0)
tabla2$intersection_Burdwood <- ifelse(present_in_small5, 1, 0)
tabla2$intersection_DR <- ifelse(present_in_small6, 1, 0)


# Save
write.csv(tabla2, file = "datos_agrupados_pre-molt_Interac_MPAs_PF_all_tracks.csv", row.names = FALSE)



### Number of days in the MPAs for each individual #######

data2 <- read.csv("C:/Users/User/Desktop/GLS projectos/PPA_Franklin/Results/overlaps_with_all_track/datos_agrupados_pre-molt_Interac_MPAs_PF_all_tracks.csv", header = TRUE)


data2$dtime <- as.Date(data2$dtime, format = "%d/%m/%Y")

cols_interes <- c("intersection_FP", "intersection_Yaganes", "intersection_Burdwood", "intersection_DR")

results2 <- data2 %>%
  group_by(id) %>%
  summarise(across(all_of(cols_interes), ~mean(. == 1, na.rm = TRUE) * 100))

print(results2)

write.csv(results2, file = "days_at_MPAs_FZ_all_tracks.csv", row.names = FALSE)


# Add a column that sums the percentages of Yaganes, Burdwood and DR by id to the summary table and general table
results2<-read.csv("C:/Users/User/Desktop/GLS projectos/PPA_Franklin/Results/overlaps_with_all_track/days_at_MPAs_PF_all_tracks.csv", header = TRUE)

results2$Sum_MPAs <- rowSums(results2[ , c(3:5)], na.rm=TRUE)
 
write.csv(results2, file = "days_at_MPAs_FZ_all_tracks.csv", row.names = FALSE)

data2$Sum_MPAs <- rowSums(data2[ , c(8:10)], na.rm=TRUE)

write.csv(data2, file = "days_at_MPAs_FZ_all_tracks_raw.csv", row.names = FALSE)

data2$dtime <- as.Date(data2$dtime, format = "%d/%m/%Y")

results <- data2 %>%
  filter(Sum_MPAs == 1) %>%
  group_by(id) %>%
  summarize(from = min(dtime), to = max(dtime))%>%  
  mutate(days_at_MPAs = as.numeric(to - from, units = "days"))

# total days for each id
total_days <- data %>%
  group_by(id) %>%
  summarise(total_days_trip = n_distinct(dtime))

# Unir los resultados anteriores y calcular el porcentaje de dias en el frente polar
resultados <- resultados %>%
  left_join(total_days, by = "id") %>%
  mutate(percent_of_thetrip = (days_at_frontal_zone / total_days_trip) * 100)

# Mostrar los resultados finales
print(results)


write.csv(results, file = "days_at_frontal_zone.csv", row.names = FALSE)


###### Differences between sexes in the MPAs and Polar Front overlap ####
### using all tracks! ###

library(MuMIn)
library(dplyr)
library(lme4)
library(nlme)
library(sjPlot)

data4 <- read.csv("C:/Users/User/Desktop/GLS projectos/PPA_Franklin/Results/overlaps_with_all_track/days_at_MPAs_PF_all_tracks_raw.csv", header = TRUE)
str(data4)

# sex vs polar front
GLMM.bino.PF<-glmer(intersection_FP~sex+(1|id),data=data4,family="binomial",na.action=na.fail)
summary(GLMM.bino.PF)
#plot(GLMM.bino.PF)

ModelSel.BSR<-dredge(GLMM.bino.PF,rank="AICc")
ModelSel.BSR
MMI.NDI.Pss<-model.avg(ModelSel.BSR,revised.var=TRUE)
summary(MMI.NDI.Pss)
coef(MMI.NDI.Pss)
confint(MMI.NDI.Pss)

tab_model(GLMM.bino.PF)

#sex vs MPAs

GLMM.bino.MPAs<-glmer(Sum_MPAs~sex+(1|id),data=data4,family="binomial",na.action=na.fail)
summary(GLMM.bino.MPAs)
#plot(GLMM.bino.MPAs)

ModelSel.BSR<-dredge(GLMM.bino.MPAs,rank="AICc")
ModelSel.BSR
MMI.NDI.Pss<-model.avg(ModelSel.BSR,revised.var=TRUE)
summary(MMI.NDI.Pss)
coef(MMI.NDI.Pss)
confint(MMI.NDI.Pss)

tab_model(GLMM.bino.MPAs)


