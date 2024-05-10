library(INLA)
library(tidyverse)
library(sdmTMB)
library(glmmTMB)
library(broom)
library(sf)
library(gstat)
library(rnaturalearth)
library(raster)
library(concaveman)

# Set coordinate reference system
ncrs <- "+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"

# polyon for cropping objects
ebs <- st_read(here::here("data/ebs_friendly.kml")) %>% 
  st_transform(., crs = ncrs)

# shoreline
shoreline <- 
  ne_states(country = "United States of America") %>% 
  filter(name == "Alaska") %>% 
  st_union() %>% 
  st_transform(., crs = ncrs) %>% 
  st_intersection(., ebs) 

# sample data to work with
# write.csv(cgi::cpue_bbrkc %>% 
#             filter(mat_sex == "legal male",
#                    year %in% 2010:2011) %>% 
#             dplyr::select(-yw), file = here::here("data/sample_bbrkc.csv"))


sample <- read_csv(file = here::here("data/sample_bbrkc.csv")) %>% 
  filter(year == 2010) %>% 
  mutate(lon_group = factor(
           cut_number(lon, n = 10))
         )

# Question 1: What is the mean density of RKC in the EBS?

# Example 1: Tweedie GLM. We estimate crab densities using an intercept-only
# GLM

# Tweedie GLM----
m1 <- sdmTMB(cpue ~ 1, data = sample, family = tweedie(), spatial = "off")

# Extract predictions and SEs. In the intercept only model, the intercept is the 
# estimate of the mean density
glm_pred_df <- predict(m1, se_fit = T) %>% 
  mutate(glm_est = exp(est),
         glm_upr = exp(est + 1.96 * est_se),
         glm_lwr = exp(est - 1.96 * est_se)) %>% 
  dplyr::select(-1)


ggplot(glm_pred_df) + 
  geom_histogram(aes(cpue)) + 
  geom_vline(aes(xintercept = glm_est), color = "purple", linewidth = 1) + 
  geom_vline(aes(xintercept = glm_lwr), color = "purple", linetype = 2, linewidth = 1) +
  geom_vline(aes(xintercept = glm_upr), color = "purple", linetype = 2, linewidth = 1) + 
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_minimal() + 
  labs(y = "Frequency",
       x = "RKC per sq. km") + 
  theme(axis.title = element_text(size = 16))

# Look at residual patterning
sample$glm_resids <- residuals(m1, type = "mle-mvn")

# visualize residuals across the EBS
ggplot(sample) + 
  geom_sf(data = shoreline) +
  geom_point(aes(y = lat, x = lon, color = glm_resids), size = 5) +
  scale_color_gradient2(midpoint = 0) + 
  labs(y = "Latitude",
       x = "Longitude") +
  theme_gray() + 
  theme(axis.title = element_text(size = 16)) + 
  coord_sf(ylim = c(550, 1000),
           xlim = c(-850, -250))

# visualize using a variogram
v1_out <- tibble()
for (i in 1:500){
  sample$glm_resids <- residuals(m1, type = "mle-mvn")
  ex1_sp <- sample
  coordinates(ex1_sp) <- c("lon", "lat")
  
  v1 <- variogram(object = glm_resids ~ 1,
                  data = ex1_sp,
                  cressie = T,
                  cutoff = 500) %>% 
    mutate(draw = i)
  assign("v1_out", rbind(v1, v1_out))
}

ggplot(v1_out) + 
  geom_boxplot(aes(y = gamma,
                   x = dist,
                   group = dist)) +
  labs(y = "Sample Variogram",
       x = "Distance (km)") + 
  theme_gray() + 
  theme(axis.title = element_text(size = 15))

# Example 2:
# Tweedie GLMM -----
ggplot(sample) + 
  geom_sf(data = shoreline) +
  geom_point(aes(y = lat, x = lon, color = lon_group), size = 5) +
  labs(y = "Latitude",
       x = "Longitude") +
  scale_color_viridis_d() +
  theme_gray() + 
  theme(axis.title = element_text(size = 16)) + 
  coord_sf(ylim = c(550, 1000),
           xlim = c(-850, -250))

# add random intercepts that cluster observations nearby in space
m2 <- sdmTMB(cpue ~ 1 + (1|lon_group), 
             data = sample,
             family = tweedie(), 
             spatial = "off")

# GLMM estimate of mean crab density after accounting for 
# group-level variability ("population-level estimate")
glmm_pop_pred <- predict(m2, re_form_iid = NA, se_fit = T) %>% 
  mutate(glmm_est = exp(est),
         glmm_upr = exp(est + 1.96 * est_se),
         glmm_lwr = exp(est - 1.96 * est_se))

# Group-level predictions
glmm_ind_pred <- predict(m2) %>% 
  mutate(glmm_ind_est = exp(est)) %>% 
  dplyr::select(.,glmm_ind_est)

# combine
glmm_pred_df <- bind_cols(glm_pred_df, 
                          glmm_pop_pred %>% 
                              dplyr::select(
                                glmm_est,
                                glmm_upr,
                                glmm_lwr
                              ),
                          glmm_ind_pred)

ggplot(glmm_pred_df) + 
  geom_histogram(aes(cpue),
                 binwidth = 20) + 
  # GLM
  geom_vline(aes(xintercept = glm_est), color = "purple", linewidth = 2) + 
  geom_vline(aes(xintercept = glm_lwr), color = "purple", linetype = 2, linewidth  = 2) +
  geom_vline(aes(xintercept = glm_upr), color = "purple", linetype = 2, linewidth = 2) + 
  
  # GLMM
  geom_vline(aes(xintercept = glmm_est), color = "darkorange", linewidth = 2) + 
  geom_vline(aes(xintercept = glmm_lwr), color = "darkorange", linetype = 2, linewidth = 2) +
  geom_vline(aes(xintercept = glmm_upr), color = "darkorange", linetype = 2, linewidth = 2) + 
  
  labs(y = "Frequency", x = "CPUE")+
  
  scale_y_continuous(expand = c(0.01, 0.01)) +
  scale_x_continuous(limits = c(-10, 1000)) +
  theme(axis.title = element_text(size = 15))

# Look at residual patterning
sample$glmm_resids <- residuals(m2, type = "mle-mvn")

# visualize residuals across the EBS
ggplot(sample) + 
  geom_sf(data = shoreline) +
  geom_point(aes(y = lat, x = lon, color = glmm_resids), size = 5) +
  scale_color_gradient2(midpoint = 0) + 
  labs(y = "Latitude",
       x = "Longitude") +
  theme_gray() + 
  theme(axis.title = element_text(size = 16)) + 
  coord_sf(ylim = c(550, 1000),
           xlim = c(-850, -250))


# visualize using a variogram
v2_out <- tibble()
for (i in 1:500){
  sample$glmm_resids <- residuals(m2, type = "mle-mvn")
  ex2_sp <- sample
  coordinates(ex2_sp) <- c("lon", "lat")
  
  v2 <- variogram(object = glmm_resids ~ 1,
                  data = ex2_sp,
                  cressie = T,
                  cutoff = 500) %>% 
    mutate(draw = i)
  assign("v2_out", rbind(v2, v2_out))
}

ggplot(v2_out) + 
  geom_boxplot(aes(y = gamma,
                   x = dist,
                   group = dist)) +
  labs(y = "Sample Variogram",
       x = "Distance (km)") + 
  theme_gray() + 
  theme(axis.title = element_text(size = 15))

# The Tweedie GLMM does a much better job capturing spatial dependencies in the data,
# and is therefore more useful for our original purposes. However, by assigning random
# intercepts to various groups, we discretize a continuous spatial process, and we can 
# actually model this underlying spatial process explicitly using spatially 
# correlated random intercepts. These intercepts allow for capturing spatial variability,
# similar to the random intercept from the GLMM, but differ in a few key ways: 1) the spatial 
# random effect is indexed at each observation rather than at groups of observations. 2) We no
# longer assume independence in how these intercepts relate to one another. Instead, we assume
# that they are correlated, and also that their correlations decay as a function of distance. 

# Ultimately, this allows us to generate predictions with uncertainty at unmeasured locations. 
# This is also called "kriging".

ggplot(glmm_pred_df) + 
  geom_sf(data = shoreline) +
  geom_point(aes(y = lat, 
                 x = lon,
                 color = glmm_ind_est,
                 size = glmm_ind_est)) +
  labs(y = "Latitude",
       x = "Longitude") +
  theme_gray() + 
  theme(axis.title = element_text(size = 16)) + 
  coord_sf(ylim = c(550, 1000),
           xlim = c(-850, -250))


# the shortcomings become more clear when we predict in 1D
ggplot(glmm_pred_df) +
  geom_point(aes(y = cpue, 
                 x = lon)) +
  geom_point(aes(y = glmm_ind_est,
                x = lon,
                group = lon_group), color = "purple",size = 1.5) +
  geom_hline(aes(yintercept = glmm_est),
             color = "darkorange") +
  labs(y = "CPUE",
       x = "Longitude") +
  theme(axis.title = element_text(size = 15))
  
# Tweedie GLMM w/ spatial RE-----

# create buffer region around survey locs to use in mesh
survey_buffed <-
  sample %>% 
  st_as_sf(., coords = c('lon','lat'), crs = ncrs) %>% 
  concaveman::concaveman(.) %>%
  st_union() %>%
  st_buffer(dist = 15) %>%
  {. ->> survey_buffed_sf} %>% 
  as_Spatial()

plot(survey_buffed)

# rasterize the polygon to create the prediction grid
r <- raster(extent(survey_buffed),
            res = 10,
            crs = crs(ncrs))

# prediction grid 
pred_grid <- 
  fasterize::fasterize(survey_buffed_sf, r) %>% 
  as("SpatialPixelsDataFrame") %>%
  as_tibble() %>% 
  dplyr::select(lon = x, 
                lat = y)
  
## INLA mesh----
mesh <- INLA::inla.mesh.2d(boundary = survey_buffed,
                           max.edge = c(20, 50), # max allowed triangle edge length inside and outside boundary
                           offset = c(10, 40), # extension of inner and outer boundary
                           cutoff = 20) #minimum distance between points
plot(mesh)

## convert for sdmTMB----
ebs_spde <- sdmTMB::make_mesh(data = sample,
                              xy_cols = c("lon", "lat"),
                              mesh = mesh)
plot(ebs_spde)

# fit a spatial GLMM
m3 <- sdmTMB(cpue ~ 1, 
             spatial = "on",
             mesh = ebs_spde,
             family = tweedie(),
             data = sample)

# evaluate residuals
# Look at residual patterning
sample$s_glmm_resids <- residuals(m3, type = "mle-mvn")

# visualize residuals across the EBS
ggplot(sample) + 
  geom_sf(data = shoreline) +
  geom_point(aes(y = lat, x = lon, color = s_glmm_resids), size = 5) +
  scale_color_gradient2(midpoint = 0) + 
  labs(y = "Latitude",
       x = "Longitude") +
  theme_gray() + 
  theme(axis.title = element_text(size = 16)) + 
  coord_sf(ylim = c(550, 1000),
           xlim = c(-850, -250))


# visualize using a variogram
v3_out <- tibble()
for (i in 1:500){
  sample$s_glmm_resids <- residuals(m3, type = "mle-mvn")
  ex3_sp <- sample
  coordinates(ex3_sp) <- c("lon", "lat")
  
  v3 <- variogram(object = s_glmm_resids ~ 1,
                  data = ex3_sp,
                  cressie = T,
                  cutoff = 500) %>% 
    mutate(draw = i)
  assign("v3_out", rbind(v3, v3_out))
}

ggplot(v3_out) + 
  geom_boxplot(aes(y = gamma,
                   x = dist,
                   group = dist)) +
  labs(y = "Sample Variogram",
       x = "Distance (km)") + 
  theme_gray() + 
  theme(axis.title = element_text(size = 15))

# predict from the model to extract spatial random field (LOG SCALE)
pred_df <- predict(m3, newdata = pred_grid)

# model predictions (global intercept + spatial random field)
ggplot(pred_df) +
  geom_raster(aes(y = lat, x = lon, fill = est)) + 
  geom_sf(data = shoreline) +
  scale_fill_gradient2() + 
  labs(y = "Latitude",
       x = "Longitude",
       fill = "Log RKC per sq. km") +
  theme_gray() + 
  theme(axis.title = element_text(size = 16),
        legend.position = "bottom") + 
  coord_sf(ylim = c(550, 1000),
           xlim = c(-850, -250))

# Spatial random field - spatial interpolation of u_s
# not different from above because we only have intercept + u_s in model
ggplot(pred_df) +
  geom_raster(aes(y = lat, x = lon, fill = est_rf)) + 
  geom_sf(data = shoreline) +
  scale_fill_gradient2() + 
  labs(y = "Latitude",
       x = "Longitude",
       fill = "Log RKC per sq. km") +
  theme_gray() + 
  theme(axis.title = element_text(size = 16),
        legend.position = "bottom") + 
  coord_sf(ylim = c(550, 1000),
           xlim = c(-850, -250))


# Tweedie GLMM w/ spatial and spatiotemporal RE----
sample_st <- read_csv(file = here::here("data/sample_bbrkc.csv"))

survey_buffed_st <-
  sample_st %>% 
  st_as_sf(., coords = c('lon','lat'), crs = ncrs) %>% 
  concaveman::concaveman(.) %>%
  st_union() %>%
  st_buffer(dist = 15) %>%
  {. ->> survey_buffed_sf_st} %>% 
  as_Spatial()

plot(survey_buffed_st)

# rasterize the polygon to create the prediction grid
r <- raster(extent(survey_buffed_st),
            res = 10,
            crs = crs(ncrs))

# prediction grid 
pred_grid_st <- 
  fasterize::fasterize(survey_buffed_sf_st, r) %>% 
  as("SpatialPixelsDataFrame") %>%
  as_tibble() %>% 
  dplyr::select(lon = x, 
                lat = y) %>% 
  expand_grid(., # expand prediction grid across years
              year = 2010:2011)

## INLA mesh----
mesh_st <- INLA::inla.mesh.2d(boundary = survey_buffed_st,
                           max.edge = c(20, 50), # max allowed triangle edge length inside and outside boundary
                           offset = c(10, 40), # extension of inner and outer boundary
                           cutoff = 20) #minimum distance between points
plot(mesh_st)

## convert for sdmTMB----
ebs_spde_st <- sdmTMB::make_mesh(data = sample_st,
                              xy_cols = c("lon", "lat"),
                              mesh = mesh_st)
plot(ebs_spde_st)

# fit a GLMM with spatial and IID spatiotemporal random effects
m4 <- sdmTMB(cpue ~ 1, 
             spatial = "on",
             time = "year",
             spatiotemporal = "IID",
             mesh = ebs_spde_st,
             family = tweedie(),
             data = sample_st)

# predict from the model to extract spatial random field (LOG SCALE)
pred_df_st <- predict(m4, newdata = pred_grid_st)

# model predictions (global intercept + 
# spatial random field + spatiotemporal random field)
ggplot(pred_df_st) +
  geom_raster(aes(y = lat, x = lon, fill = est)) + 
  geom_sf(data = shoreline) +
  facet_wrap(~year) +
  scale_fill_gradient2() + 
  labs(y = "Latitude",
       x = "Longitude",
       fill = "Log RKC per sq. km") +
  theme_gray() + 
  theme(axis.title = element_text(size = 16),
        legend.position = "bottom") + 
  coord_sf(ylim = c(550, 1000),
           xlim = c(-850, -250))

