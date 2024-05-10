library(INLA)
library(tidyverse)
library(sdmTMB)
library(glmmTMB)
library(broom)
library(ggbreak)
library(sf)
library(gstat)
library(rnaturalearth)
library(raster)


# sample_sc <- cgi::cpue_sc %>%
#   filter(year %in% 2010:2015,
#          mat_sex == "legal male")
# write.csv(sample_sc %>% 
#             dplyr::select(-yw), 
#           file = here::here("data/sample_snow_crab.csv"))

sample_sc <- read_csv(here::here("data/sample_snow_crab.csv"))

# Set coordinate reference system
ncrs <- "+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"

# polyon for cropping objects
ebs_large <- st_read(here::here("data/ebs_large.kml")) %>% 
  st_transform(., crs = ncrs)

# shoreline
shoreline <- 
  ne_states(country = "United States of America") %>% 
  filter(name == "Alaska") %>% 
  st_union() %>% 
  st_transform(., crs = ncrs) %>% 
  st_intersection(., ebs_large) %>% 
  st_zm()

plot(shoreline)

# survey area
survey_buffed_st <-
  sample_sc %>% 
  st_as_sf(., coords = c('lon','lat'), crs = ncrs) %>% 
  concaveman::concaveman(.) %>% 
  st_union() %>%
  st_buffer(dist = 15) %>%
  st_difference(.,shoreline) %>% # excise the islands
  {. ->> survey_buffed_sf_st} %>%
  as_Spatial()

ggplot() + 
  geom_sf(data = shoreline) +
  geom_sf(data= survey_buffed_sf_st) + 
  geom_point(data = sample_sc,
             aes(y = lat, x = lon))
