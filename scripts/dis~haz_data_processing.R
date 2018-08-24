### This R script creates a big data frame with difffrent attributes for all watersheds. It uses Watersheds 
### wich have experienced a natural hazard event during the last 30 years (watersheds1) and those which have not (watersheds0)
### and merges them to a complete dataset for the model analysis.
### Author J. Sebald 
### Data from C.Scheidl, M. Heiser and D. Pflugmacher

{

library(dplyr)# version 0.7.5
library(raster)# version 2.6-7
library(sp)# version 1.3-1
library(rgdal)# version 1.2-20
library(igraph)# version 1.2.1
library(ggplot2)# version 3.0.0
library(tidyverse)
library(data.table)
#rm(list=ls())

#server <- "/home/jsebald/upload_to_server"
local <- "D:/JULIUS/PhD/Projects/disturbances_and_natural_hazards/"

setwd(local)
#setwd(server)

rm(list=ls())
}


# Load raw data --------------------------------------------------------------

# shapefile of watersheds with event

shp_ws1 <- raster::shapefile("materials/raw_data/data_natural_hazards/GroupEvents_final.shp")

# shapefile of watersheds without event

shp_ws0 <- raster::shapefile("materials/raw_data/data_natural_hazards/GroupNOEvents_final.shp")

# landcover data

landcover <- raster("materials/raw_data/data_landcover/europe_landcover2015_lc3_3year_austria.tif") 

# disturbance data

disturbance <- raster("materials/raw_data/data_disturbances/austria_ensemble_aggregate32_19862016_firstdisturbance_laea_mmu6.tif")
disturbance <- projectRaster(disturbance, landcover, method = "ngb") 
writeRaster(disturbance, "methods/r/data_processed/rasters/disturbance.tif")

disturbance <- raster("methods/r/data_processed/rasters/disturbance.tif")

# DEM austria

dem_austria <- raster("methods/r/data_processed/rasters/dem_austria.tif")

# natural hazard events

events <- read.csv("materials/raw_data/data_natural_hazards/GroupEvents.csv") 

# ecological units

eco <- raster::shapefile("materials/raw_data/shapefiles_ecological_units/WLamPoly.shp")

# Clean shp_ws1 and shp_ws0 to build rasterfile of complete studyarea ----------------------------

# transform column WLK_ID to numeric, otherwise rasterize function will get confused

shp_ws1$WLK_ID <- as.integer(as.character(shp_ws1$WLK_ID))
shp_ws0$WLK_ID <- as.integer(as.character(shp_ws0$WLK_ID))

# create raster of watersheds with event (needed later)

ras <- raster(ncol = 180, nrow = 180)
extent(ras) <- extent(shp_ws1)
res(ras) <- res(landcover)
crs(ras) <- crs(landcover)
rast_ws1 <- rasterize(shp_ws1, ras, field = "WLK_ID")
raster_ws1 <- projectRaster(rast_ws1, landcover, method = "ngb")

# i don`t no why but obviously the rasterize function eats 8 out of the 2204 watersheds and gives them the IDs of the surrounding watershed. 
# after some back and forth i decided to simply filter them out after the area is included in the surrounding watershed anyways
# the result is a shapefile and a new dataframe without these 8 watersheds

# have a look at the diffrent WLK_IDs which exist in the shapefile and the rasterlayer

wlks_raster <- unique(values(raster_ws1))
wlks_shp <- shp_ws1$WLK_ID

# check if there are WLK_IDs which are in the shapefile but aren´t in the raster

diff <- setdiff(wlks_shp, wlks_raster)

# subset the watersheds which are eaten by the rasterize function

eaten <- subset(shp_ws1, WLK_ID %in% diff)

# built raster from the subsetted polygons

rast_eaters <- crop(rast_ws1, extent(eaten)) 
raster_eaters <-  mask(rast_eaters, eaten)

# get new WLK_IDs of the watersheds which where eaten by the rasterize function 

wlks_eaters <- sort(unique(values(raster_eaters)))

# create polygon layer which shows the watersheds which have eaten the missing WLK-IDs 

eaters <- subset(shp_ws1, WLK_ID %in% wlks_eaters)

# show everything with plots

plot(eaters)
plot(eaten, col = "red", add = T)

# create new watersheds-shapefile wich dosen`t contain the eaten watersheds anymore

without_eaten <- subset(shp_ws1, ! WLK_ID %in% diff)

# update watersheds1 so that it contains only events between 1986 and 2018 
# and that the events of the eaten watersheds get the WLK_ID from those which have eaten them (the surrounding ones)


# create vector of WLK-IDs without the eaten watersheds

without_eaten_wlks <- without_eaten$WLK_ID


# get events data.frame, filter out all events between 1980 and 1985 
# and filter out the 8 events which happend in the eaten watersheds

events_clean <- events %>%
  filter(ereignis_jahr %in% c(1986:2018)) %>%
  filter(WLK_ID %in% without_eaten_wlks)

# create subset of events dataframe which contains only the events of the eaten watersheds 

events_in_eaten <- events %>%
  filter(! WLK_ID %in% without_eaten_wlks)

# create vector of the eater-watersheds. These are the ones which are surrounding the eaten watersheds
# i created this vector manually because two small watersheds are eaten by the same big watersheds and to make sure the order
# of the WLK-IDs is correct. I did this in QGIS by eyeballing.

wlk_eaters <- c(104593,105457,103284,104258,104270,104544,104097,104097)

# assign wlk-id of the eaters watersheds to the events which happend in the eaten watersheds

events_in_eaten$WLK_ID <- wlk_eaters

# bring back events of eaten watersheds to the events dataframe (now with wlk_ids of surrounding "eaters" watersheds)

events_clean_new <- bind_rows(events_clean, events_in_eaten)


# Change process names

events_ws <- events_clean_new %>%   #CHECK WHAT IS ACTUALLY NEEDED FOR THE MODEL!!!!
  rename(year = ereignis_jahr) %>%
  mutate(Event = case_when(
    prozessart == "Fluviatiler Feststofftransport" ~ "FST",
    prozessart ==  "Murartiger Feststofftransport" ~ "DFLOOD",
    prozessart == "Murgang" ~ "DFLOW")) %>%
  dplyr::select(WLK_ID, year, Event) %>%
  group_by(Event, WLK_ID) %>%
  summarize(n_eve = n()) %>%
  spread(key = Event, value = n_eve) %>%
  mutate_at(.vars = vars(DFLOOD, DFLOW, FST), function(x) ifelse(is.na(x), 0, x))

summary(events_ws)

write_csv(events_ws, "methods/r/data_processed/dataframes/events_ws.csv")

events_ws <- read.csv("methods/r/data_processed/dataframes/temp/events_ws.csv") 

# create vector of WLK_IDs of shp_ws1_new

wlks_ws1_new <- events_ws$WLK_ID


# subset shp_ws1 so that it doesn`t contain watersheds which have expierienced an event only between 1980 and 1985 
# and watersheds which were eaten by the rasterize function anymore


shp_ws1_new <- subset(shp_ws1, WLK_ID %in% wlks_ws1_new)


# join together shp_ws1_new and shp_ws0 to get final shapefile of study-area


temp_ws1 <- shp_ws1_new[,1]
temp_ws0 <- shp_ws0[,1]

shp_ws <- rbind(temp_ws1, temp_ws0)

writeOGR(shp_ws, "methods/r/data_processed/shapefiles", layer = "shp_ws", driver = "ESRI Shapefile")

#

# create vector of WLK_IDs of final studyarea

wlks_ws <- as.integer(shp_ws$WLK_ID)

# built raster from shapefile of final studyarea

ras_ws <- raster(ncol = 180, nrow = 180)
extent(ras_ws) <- extent(shp_ws)
res(ras_ws) <- res(landcover)
crs(ras_ws) <- crs(landcover)
rast_ws <- rasterize(shp_ws, ras_ws, field = "WLK_ID")
raster_ws <- projectRaster(rast_ws, landcover, method = "ngb")

writeRaster(raster_ws, "methods/r/data_processed/rasters/raster_ws.tif")


# Data spatial ---------------------------------------------------

# load rasterfile of all watersheds

raster_ws <- raster("methods/r/data_processed/rasters/raster_ws.tif")


# Geomorphology -----------------------------------------------------------

df_ws1 <- as_data_frame(shp_ws1) %>%
  dplyr::select(WLK_ID, area, Melton___, Elevation, Circularit, Relief_rat, Elongation, Form_facto) %>% 
  rename (Melton = Melton___)

df_ws0 <- as_data_frame(shp_ws0) %>%
  dplyr::select(- ART_IDX) %>%
  rename(Melton = MELTON, Elevation = Err)

geomorphology_ws <- bind_rows(df_ws1, df_ws0) %>%
  filter(WLK_ID %in% wlks_ws) 

summary(geomorphology_ws)

write_csv(geomorphology_ws, "methods/r/data_processed/dataframes/old/geomorpholgy_ws.csv")

geomorphology_ws <- read.csv("methods/r/data_processed/dataframes/temp/geomorpholgy_ws.csv") 
  
# Landcover and topography ---------------------------------------------------------------

# create seperate rasters for every landuse-category; assaign 1 for pixels with
# a specific category an 0 for others, they are needed later to execute zonal statistics for the diffrent categories

artif <- mask(landcover == 1, landcover, maskvalue = 2 , updatevalue= 0)
forest <- mask(landcover %in% c(4:6), landcover, maskvalue = 1 , updatevalue= 0) 
clumps <- clump(forest, directions = 8, gaps = F)
slope <- terrain(dem_austria, opt = "slope", neighbours = 8, unit = "degrees") 

# store rasters, otherwise environemnt will get very large very soon

writeRaster(artif, "methods/r/data_processed/rasters/artif.tiff")
writeRaster(forest, "methods/r/data_processed/rasters/forest.tiff")
writeRaster(clumps, "methods/r/data_processed/rasters/clumps.tiff")
writeRaster(slope, "methods/r/data_processed/rasters/slope.tiff")

# after creating and storeing the landcover rasters they can be loaded in from the harddrive

artif <- raster("methods/r/data_processed/rasters/artif.tif")
forest <- raster("methods/r/data_processed/rasters/forest.tif")
clumps <- raster("methods/r/data_processed/rasters/clumps.tif")
dem_austria <- raster("methods/r/data_processed/rasters/dem_austria.tif")

# calculate forest-cover, share of artfical landcover (buidlings, infrastructure) and forest-patchdensity for every watershed

values_landcover <- data.table(WLK_ID = values(raster_ws),
                               forest = values(forest),
                               artifical = values(artif),
                               clumps = values(clumps)) %>%
  filter(!is.na(WLK_ID)) 
    

landcover_ws <- values_landcover %>%
  group_by(WLK_ID) %>%
  summarise( forest = sum(forest) / length(forest),
             artifical = sum(artifical) / length(artifical),
             clumps = length(unique(clumps))
             ) 

summary(landcover_ws)

write_csv(landcover_ws,"methods/r/data_processed/dataframes/old/landcover_ws.csv")

landcover_ws <- read.csv("methods/r/data_processed/dataframes/temp/landcover_ws.csv")


# calculate mean slope and mean elevation for every watershed

values_topography <- data.table(WLK_ID = values(raster_ws), 
                                elevation = values(dem_austria)) %>%
  filter(!is.na(WLK_ID))  


topography_ws <- values_topography %>%
  group_by(WLK_ID) %>%
  summarise(h_mean = mean(elevation)) 


summary(topography_ws)

write_csv(topography_ws, "methods/r/data_processed/dataframes/temp/topography_ws.csv")
topography_ws <- read.csv("methods/r/data_processed/dataframes/temp/topography_ws.csv")



# Ecological units --------------------------------------------------------

eco$Wuchsge1 <- as.integer(as.character(eco$Wuchsge1))

ras <- raster(ncol = 180, nrow = 180)
extent(ras) <- extent(eco)
res(ras) <- res(landcover)
crs(ras) <- crs(landcover)
rast_eco <- rasterize(eco, ras, field = "Wuchsge1")
eco_units <- projectRaster(rast_eco, landcover, method = "ngb")

writeRaster(eco_units, "methods/r/data_processed/rasters/eco_units.tif")

eco_units <- raster("methods/r/data_processed/rasters/eco_units.tif")

values_eco <- data.table(WLK_ID = values(raster_ws),
                         eco_unit = values(eco_units)) %>%
  filter(!is.na(WLK_ID)) %>%
  filter(!is.na(eco_unit))

summary(values_eco)

eco_ws <- values_eco %>%
  group_by(WLK_ID) %>%
  summarise(eco_unit = raster::modal(eco_unit))


summary(eco_ws)

write_csv(eco_ws, "methods/r/data_processed/dataframes/temp/eco_ws.csv")
eco_ws <- read.csv("methods/r/data_processed/dataframes/temp/eco_ws.csv")


# Disturbances ------------------------------------------------------------


values_disturbance <- data.table(WLK_ID = values(raster_ws),
                      forest = values(forest),
                      disturbance = values(disturbance)) %>%
  filter(!is.na(WLK_ID)) 

disturbance_ws <- values_disturbance %>%
  group_by(WLK_ID) %>%
  summarise(severity = sum(disturbance > 0) / sum(forest),
            frequency = DescTools::Gini(table(factor(disturbance[disturbance > 0], levels = 1986:2016)))) %>%
  mutate_all(function(x) ifelse(is.na(x) | is.nan(x), 0, x))

summary(disturbance_ws)

write_csv(disturbance_ws, "methods/r/data_processed/dataframes/temp/disturbance_ws.csv")
disturbance_ws <- read.csv("methods/r/data_processed/dataframes/temp/disturbance_ws.csv")



# Join together final dataframe -------------------------------------------


data_for_model <- geomorphology_ws %>%
  left_join(topography_ws, by = "WLK_ID") %>%
  left_join(landcover_ws, by = "WLK_ID") %>%
  left_join(disturbance_ws, by = "WLK_ID") %>%
  left_join(events_ws, by = "WLK_ID") %>%
  left_join(eco_ws, by = "WLK_ID") %>%
  mutate(patchdensity = clumps / (area*100)) %>%
  dplyr::select(-clumps) %>%
  mutate_at(.vars = vars(DFLOOD, DFLOW, FST), function(x) ifelse(is.na(x), 0, x))


summary(data_for_model)

write_csv(data_for_model, "methods/r/data_processed/dataframes/data_for_model_08222018.csv")

