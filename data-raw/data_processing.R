library(tidyverse)
library(sf)
library(raster)
library(magrittr)
library(INLA)
library(fasterize)

# Set coordinate reference system with units specified as km instead of m
# sdmTMB doesn't like coordinates in m
ncrs <- "+proj=utm +zone=18 +datum=WGS84 +units=km +no_defs"

# load VIMS shapefile for Hog Island
load(here::here("data/HogIslandBaySeagrass.rdata"))
hi_sg %<>% 
  filter(YEAR == 2018) %>% 
  st_transform(ncrs) %>% 
  st_union()

# load adjusted shapefile that extends the meadow to low density or zero sites
hi_sg_adjusted <- st_read(here::here("data/hi_adjusted.kml")) %>% 
  st_transform(ncrs)

# load shoot density data from VCR LTER and process
hi_2018_init <- readxl::read_xlsx(here::here("data/hog_island_2018.xlsx")) %>% 
  mutate(site = paste("HI",site))

# read in synoptic seagrass sampling data and point locations
syn_locs <- st_read(here::here("data/Syn_Site_Points_2018")) %>% 
  dplyr::rename(site = SiteName) %>% 
  mutate(site = str_replace(site, "-", " "))

all(hi_2018_init$site %in% syn_locs$site)

# bind the locations and site names, transform to new CRS
hi_2018 <- hi_2018_init %>% 
  left_join(.,syn_locs) %>% 
  st_as_sf(.) %>% 
  st_transform(st_crs(hi_sg)) %>% 
  mutate(longitude = Easting__X/1000,
         latitude = Northing__/1000) %>% 
  dplyr::select(-Easting__X, -Northing__) 

hi_2018_df <- 
  hi_2018 %>% 
  st_set_geometry(NULL)

ggplot() +
  geom_sf(data = hi_sg_adjusted, fill = "transparent") +
  geom_sf(data = hi_sg, color = "red", fill = "transparent") +
  geom_sf(data = hi_2018) +
  theme_bw()

# Define a raster based on the adjusted meadow polygon - for predicting over
# once we fit the model
rast_res <- 0.025
r <- raster(extent(hi_sg_adjusted),
            res = rast_res,
            crs = crs(ncrs))
pred_df <- fasterize::fasterize(hi_sg_adjusted, raster = r) %>% 
  as("SpatialPixelsDataFrame") %>% 
  as.data.frame() %>% 
  dplyr::select(longitude = x, latitude = y)

# Create triangular mesh using the adjusted meadow polygon
max.edge = 0.1
bound.outer = 0.05
mesh <- inla.mesh.2d(boundary = hi_sg_adjusted %>% 
                       st_zm() %>% 
                       as_Spatial(),
                     max.edge = 0.1)

# View the mesh
plot(mesh)

save(mesh, r, pred_df, rast_res, hi_2018, hi_sg_adjusted, hi_2018_df,
     file = here::here("data/sg_mod_prereqs.rdata"))

