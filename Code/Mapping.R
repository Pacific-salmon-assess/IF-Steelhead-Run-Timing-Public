# #Mapping script
# install.packages("sf")
# install.packages("rnaturalearth")
# install.packages("rnaturalearthdata")
# install.packages("remotes")
# remotes::install_github("bcgov/bcmaps")

library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(bcmaps)
library(tidyverse)
library(cowplot)


#Download world map data
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
#Subset to us and canada
us_can <- world[world$admin == "United States of America" | world$admin == "Canada",]
ggplot(data = us_can)+
  geom_sf()

#View layers available in bcmaps package
layers <- available_layers()
#Get bc boundary 
bc <- bc_bound()
#Get projection of bc data - BC Albers
st_crs(bc)$proj4string

#Get drainages from bcmaps 
wsc <- wsc_drainages()
#subset to the Fraser river drainages 
fr <- wsc[wsc$SUB_DRAINAGE_AREA_NAME %in% c("Upper Fraser","Lower Fraser","Thompson", "Nechako"),]
#dissolve to make one watershed
fr <- fr %>% group_by(OCEAN_DRAINAGE_AREA_NAME) %>% summarize() 

#Create a list of point locations of interest: 
#Look up these coordinates manually 
#Barnston island = 1240005, 470241
#Albion = 1250562, 468840
#Harrison = 1295302, 476754
#Hope = 1330249, 495701
#Nahatlatch = 1321736, 561483
#Lytton = 1314665, 589166
#Stein = 1311231, 595511
#Bridge = 1286752, 645427
#Chilcotin = 1248204, 753809

#Create points data frame 
points <- data.frame(name = c("Barnston Island", "Albion", "Harrison", "Hope", "Nahatlatch", "Lytton", "Stein", "Bridge", "Chilcotin"),
                     long = c(1240005, 1250562, 1295302,1330249, 1321736, 1314665,1311231,1286752,1248204),
                     lat = c(470241, 468840, 476754, 495701, 561483, 589166, 595511, 645427, 753809))
#Remove the river points for this map
points_sub <- data.frame(name = c("Barnston Island", "Albion", "Hope", "Lytton"),
                         long = c(1240005, 1250562,1330249,  1314665),
                         lat = c(470241, 468840,495701,589166))
#Convert to sf object
points_geo <- st_as_sf(points_sub, coords = c("long", "lat"))
#Set proj to BC albers EPSG 3153
st_crs(points_geo) <- 3153

#Look at Fraser river geodatabase
st_layers("C:\\Users\\POTAPOVAA\\Documents\\fwa_network_fraser_historic.gdb\\fwa_network_fraser_historic.gdb")
#Read in the layer inside the geodatabase
fc <- sf::st_read("C:\\Users\\POTAPOVAA\\Documents\\fwa_network_fraser_historic.gdb\\fwa_network_fraser_historic.gdb", layer = "fwa_network_fraser_rollup_pnwnamet_1980_2010_adj_m3s")
#Subset to stream order > 5 for mapping 
fc <- fc[fc$STREAM_ORDER > 5,]

#Set the bounding box for main map and inset here
bbox <-c(x1 = 1152862, x2 = 1421635, y1 = 467000, y2 = 780000)

#Create map 
map <- ggplot(data = bc)+
  #geom_sf()+
  geom_sf(data = us_can)+
  geom_sf(data = fr)+
  geom_sf(data = fc)+
  geom_sf(data = points_geo, color = "black")+
  #geom_sf_text(data = points_geo, aes(label = name))+
  theme_bw()+
  xlim(c(bbox["x1"],bbox["x2"]))+
  ylim(c(bbox["y1"],bbox["y2"]))+
  xlab("Longitude")+
  ylab("Latitude")
map
ggsave(file = "Outputs/3_Plots/Map_noinset.png", height = 8, width = 6, units = "in")


#Vancouver island
van <- ggplot(data = bc)+
  #geom_sf()+
  geom_sf(data = us_can)+
  geom_sf(data = fr)+
  geom_sf(data = fc)+
  #geom_sf_text(data = points_geo, aes(label = name))+
  theme_bw()+
  xlim(c(835000,1210000))+
  ylim(c(354543,690000))+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))
van
ggsave(file = "Outputs/3_Plots/Map_van.png", height = 3, width = 3, units = "in")

#Create inset of province with box
#Need a box for vancouver island as well - can't quite get it right because of axis limits so just draw it manually
inset <- ggplot(data = bc)+
  geom_sf(data = us_can)+
  geom_sf(data = bc)+
  geom_sf(data = fr)+
  geom_rect(aes(xmin = bbox["x1"], xmax = bbox["x2"], ymin = bbox["y1"], ymax = bbox["y2"]), color = "black", fill = NA) +
  theme_bw()+
  ylim(391459,1677475)+
  xlim(449589,1836761)+
  theme(axis.text = element_blank(), axis.ticks.length = unit(0, "pt"), axis.ticks = element_blank(),plot.margin = margin(0, 0, 0, 0, "cm"))
inset
ggsave(file = "Outputs/3_Plots/Map_inset.png", height = 3, width = 3, units = "in")

#Manually add labels for the rivers in inkscape

