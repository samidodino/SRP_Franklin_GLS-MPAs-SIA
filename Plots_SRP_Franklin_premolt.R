#######  Plots codes#####
### Dodino et al. 2024 - GLS Southern Rockhoppers - Bahia Franklin ####

####### Plot by sex --> best tracks #####

library(ggplot2)
library(dplyr)
library(ggspatial)
library(patchwork)
library(sf)
library(raster)
library(RColorBrewer)
library(marmap)

### best tracks
tabla1<-read.csv("C:/Users/User/Desktop/GLS projectos/PPA_Franklin/Results/overlaps_with_best_track/datos_agrupados_pre-molt_withSex.csv",header=T)
str(tabla1)

# MPAs + Polar Front

yaganesI<-st_read("C:/Users/User/Desktop/GLS projectos/MarineAreas/Yaganes I.shp")
yaganesII<-st_read("C:/Users/User/Desktop/GLS projectos/MarineAreas/Yaganes II.shp")
burdwood<-st_read("C:/Users/User/Desktop/GLS projectos/MarineAreas/BB_I-II/BB_I_II_limite_zonificacion.shp")
diego_ramirez<-st_read("C:/Users/User/Desktop/GLS projectos/MarineAreas/Diego_Ramirez/WDPA_WDOECM_Sep2023_Public_555643508_shp-polygons.shp")
limite_chile<-st_read("C:/Users/User/Desktop/GLS projectos/Limite_con_Chile/linea_de_limite_FA004.shp")
fronts <- st_read("C:/Users/User/Desktop/GLS projectos/Fronts_of_the_Antarctic_Circumpolar_Current_-_GIS_data/shapefile/antarctic_circumpolar_current_fronts.shp")
str(fronts$NAME)
polar_front<-dplyr::filter(frentes, NAME %in% c("Polar Front (PF)"))


####### Plot all tracks + Base map ######

# Transform to gls points to sf object

premolt_sf <- st_as_sf(tabla1, coords = c("lon", "lat"), crs = 4326)

### get base map from NOAA

pt.lim = data.frame(ylim=c(-65, -50), xlim=c(-78, -48))
pt.bbox <- st_bbox(c(xmin=pt.lim$xlim[1],
                     xmax=pt.lim$xlim[2],
                     ymin=pt.lim$ylim[1],
                     ymax=pt.lim$ylim[2]))
pt.baty = getNOAA.bathy(
  lon1 = pt.lim$xlim[1], 
  lon2 = pt.lim$xlim[2], 
  lat1 = pt.lim$ylim[1], 
  lat2 = pt.lim$ylim[2], resolution = 1)



# Base map
baty_plot <- ggplot() +
  geom_raster(data = pt.baty, aes(x = x, y = y, fill = z)) +
  scale_fill_etopo() +
  geom_contour(data = pt.baty, aes(x = x, y = y, z = z),
               breaks = c(0, -10, -20, -50, -100, -200, -1000), colour = "grey", linewidth = 0.2) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_sf(xlim=c(-78, -48),ylim=c(-65, -50), expand = FALSE)

# Colony location

lon.calib <- -64.2500
lat.calib <- -54.7833

# MPAs + Polar Front

plot_premolt_puntos<-baty_plot +
  geom_sf(data = premolt_sf, aes(color = id)) + #pone los colores por id, random eleccion
  geom_point(data = tabla1, aes(x = lon.calib, y = lat.calib), size = 3, shape = 21, fill = "darkred") +
  geom_sf(data = yaganesI, fill = NA, color = "black", size = 12, alpha = 0.5) +
  geom_sf(data = yaganesII, fill = NA, color = "black", size = 12, alpha = 0.5) +
  geom_sf(data = burdwood, fill = NA, color = "black", size = 12, alpha = 0.5) +
  geom_sf(data = diego_ramirez, fill = NA, color = "black", size = 12, alpha = 0.5)+
  geom_sf(data = limite_chile, fill = NA, color = "black", size = 12, alpha = 0.5)+
  geom_sf(data = polar_front, fill=NA, color="black",size=12, alpha=0.5) +
  annotation_scale(location = "bl") +
  coord_sf(xlim=c(-78, -48),ylim=c(-65, -50), expand = FALSE) +
  theme(axis.text = element_text(size = 10))


print(plot_premolt_puntos)

plot_premolt_puntos<-baty_plot +
  geom_sf(data = premolt_sf, aes(color = sex)) +
  scale_color_manual(values = c("M" = "violetred", "F" = "springgreen"))+
  geom_point(data = tabla1, aes(x = lon.calib, y = lat.calib), size = 3, shape = 23, fill = "black") +
  geom_sf(data = yaganesI, fill = NA, color = "black", size = 13, alpha = 0.5) +
  geom_sf(data = yaganesII, fill = NA, color = "black", size = 13, alpha = 0.5) +
  geom_sf(data = burdwood, fill = NA, color = "black", size = 13, alpha = 0.5) +
  geom_sf(data = diego_ramirez, fill = NA, color = "black", size = 13, alpha = 0.5)+
  geom_sf(data = limite_chile, fill = NA, color = "black", size = 12, alpha = 0.5)+
  geom_sf(data = polar_front, fill=NA, color="black",size=12, alpha=0.5) +
  annotation_scale(location = "bl") +
  coord_sf(xlim=c(-78, -48),ylim=c(-65, -50), expand = FALSE) +
  theme(axis.text = element_text(size = 10))


print(plot_premolt_puntos)



##### Density map by sex ######

install.packages("rnaturalearthdata")
library(rnaturalearth)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")


# Males, normalized (fill=after_stat(ncount))
data_M_BT <- tabla1 %>% filter(sex == "M")

plot_males_BT<-ggplot(data = world) +
  geom_bin2d(data = data_M_BT, aes(x = lon, y = lat, fill=after_stat(ncount)), binwidth = 1) +
  geom_sf(fill = gray(0.3), colour = NA) +
  geom_point(aes(x = lon.calib, y = lat.calib), size = 1.5,
             shape = 21, fill = "darkred") +
  geom_sf(data = polar_front,col="black", size=5,fill = NA, alpha = 0.5) +
  geom_sf(data = yaganesI,col="black", size=5,fill = NA, alpha = 0.5) +
  geom_sf(data = yaganesII,col="black", size=5,fill = NA, alpha = 0.5) +
  geom_sf(data = burdwood,col="black", size=5,fill = NA, alpha = 0.5) +
  geom_sf(data = diego_ramirez,col="black", size=5,fill = NA, alpha = 0.5)+ 
  coord_sf(xlim=c(-78, -48),ylim=c(-65, -50), expand = FALSE)+
  scale_fill_continuous(type = "viridis") +
  labs(title = "Density map - Males - Best track", x = NULL, y = NULL)+
  theme_light() +
  theme(axis.text = element_text(size = 10))+
  annotation_scale(location = "bl")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Females, normalized (fill=after_stat(ncount))
data_F_BT <- tabla1 %>% filter(sex == "F")

plot_females_BT<-ggplot(data = world) +
  geom_bin2d(data = data_F_BT, aes(x = lon, y = lat, fill=after_stat(ncount)), binwidth = 1) +
  geom_sf(fill = gray(0.3), colour = NA) +
  geom_point(aes(x = lon.calib, y = lat.calib), size = 1.5,
             shape = 21, fill = "darkred") +
  geom_sf(data = polar_front,col="black", size=5,fill = NA, alpha = 0.5) +
  geom_sf(data = yaganesI,col="black", size=5,fill = NA, alpha = 0.5) +
  geom_sf(data = yaganesII,col="black", size=5,fill = NA, alpha = 0.5) +
  geom_sf(data = burdwood,col="black", size=5,fill = NA, alpha = 0.5) +
  geom_sf(data = diego_ramirez,col="black", size=5,fill = NA, alpha = 0.5) +
  coord_sf(xlim=c(-78, -48),ylim=c(-65, -50), expand = FALSE)+
  scale_fill_continuous(type = "viridis") +
  labs(title = "Density map - Females - Best track", x = NULL, y = NULL)+
  theme_light() +
  theme(axis.text = element_text(size = 10))+
  annotation_scale(location = "bl") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


# Combined plot
combined_plot_BT <- plot_females_BT + plot_males_BT
combined_plot_BT


### Density both sex ####

ggplot(data = world) +
  geom_bin2d(data = tabla1, aes(x = lon, y = lat, fill=after_stat(ncount)), binwidth = 1.5) +
  geom_sf(fill = gray(0.3), colour = NA) +
  geom_point(aes(x = lon.calib, y = lat.calib), size = 2,
             shape = 21, fill = "darkred") +
  geom_sf(data = frentes2 %>% filter(NAME %in% c("Polar Front (PF)")), size = 4) +
  geom_sf(data = yaganesI,col="black", size=4,fill = NA, alpha = 0.5) +
  geom_sf(data = yaganesII,col="black", size=4,fill = NA, alpha = 0.5) +
  geom_sf(data = burdwood,col="black", size=4,fill = NA, alpha = 0.5) +
  geom_sf(data = diego_ramirez,col="black", size=4,fill = NA, alpha = 0.5)+ 
  coord_sf(xlim=c(-78, -48),ylim=c(-65, -50), expand = FALSE)+
  scale_fill_continuous(type = "viridis") +
  labs(title = "Densisty map - ALL - Best track", x = NULL, y = NULL)+
  theme_light() +
  theme(axis.text = element_text(size = 10))+
  annotation_scale(location = "bl")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))







