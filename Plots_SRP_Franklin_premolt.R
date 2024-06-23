#######  Plots codes#####
### Dodino et al. 2024 - GLS Southern Rockhoppers - Bahia Franklin ####

####### Plot by sex --> best tracks #####

library(ggplot2)
library(dplyr)
library(ggspatial)
library(sf)
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

plot_premolt_points<-baty_plot +
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


print(plot_premolt_points)

plot_premolt_points<-baty_plot +
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


print(plot_premolt_points)

#####  GLS days from start for each individual + polar front+ marine area ########
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
load("C:/Users/User/Desktop/GLS projectos/Results/Rockhopper_franklin/prob_track/BV342_04Jan21_184313driftadj_probGLS_premolt.RData")

## confident intervals shapes

prj <- paste0("+proj=aeqd +lat_0=", median(st_coordinates(prob_track[[1]])[, 2]), " +lon_0=", 
              median(st_coordinates(prob_track[[1]])[, 1]), " +units=km")
polys <- lapply(c(95, 50, 25), function(x) {
  probGLS:::calc_mcp(prob_track[[1]][, "step"], percent = x) %>% st_transform(prj) %>% st_union() %>%
    st_transform(4326)
})
doy_colours <- colorRampPalette(
  paletteer:::paletteer_d("ggthemes::Classic_Cyclic", n = 12)[c(1:12, 1)])(366)

#plot
ggplot(data = world) +
  geom_sf(fill = gray(0.3), colour = NA) +
  geom_sf(data = polys[[1]], aes(geometry = geometry),
          fill = grey(0.9, alpha = 0.6), col = "transparent") +
  geom_sf(data = polys[[2]], aes(geometry = geometry),
          fill = grey(0.8, alpha = 0.6), col = "transparent") +
  geom_sf(data = polys[[3]], aes(geometry = geometry),
          fill = grey(0.7, alpha = 0.6), col = "transparent") +
  
  geom_sf(data = yaganesI,col=NA, fill = "orange", alpha = 0.5) +
  geom_sf(data = yaganesII,col=NA, fill = "orange", alpha = 0.5) +
  geom_sf(data = burdwood,col=NA, fill = "orange", alpha = 0.5) +
  geom_sf(data = diego_ramirez,col=NA, fill = "orange", alpha = 0.5) +
  geom_sf(data = polar_front, fill=NA, color="black",size=12, alpha=0.5)
  geom_sf(data = prob_track[[2]] %>% summarise(do_union = FALSE) %>%
            st_cast("LINESTRING"), mapping = aes(geometry = geometry), linewidth = 0.4) +
  geom_sf(data = prob_track[[2]] %>% rowid_to_column(var = "nr"), mapping = aes(geometry = geometry, size = mean.rel_weight, 
                                                                                fill = as.numeric(nr)/2), shape = 21) +
  scale_size(range = c(2, 4), name = "Mean.rel_weight") +
  scale_fill_continuous(type = "viridis", name = "day from start") +
  
  geom_point(aes(x = lon.calib, y = lat.calib), size = 3, 
             shape = 21, fill = "darkred") +
  coord_sf(xlim = c(-90,-40), ylim = c(-70, -45), expand = FALSE) +
  theme_light() +
  scale_x_continuous(name="") +
  scale_y_continuous(name="") +
  theme(axis.text = element_text(size = 10)) + 
  annotation_north_arrow(location = "br", which_north = "true", pad_x = unit(0.5, "in"), pad_y =unit(0.5,"in"), style = north_arrow_fancy_orienteering)+
  annotation_scale(location = "bl")






