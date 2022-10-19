library(sdmTMB)
library(tidyverse)
library(sf)
library(raster)
library(magrittr)
library(INLA)


load(here::here("data/sg_mod_prereqs.rdata"))

# LM demonstration----
# fit a lm using sdmTMB by turning off the spatial random fields component
linear_mod <- 
  sdmTMB(shoots_per_m ~ 1, 
       spatial = "off",
       data = hi_2018)
summary(linear_mod)

# Spatial model fitting 1-----
#1. Convert custom INLA mesh to mesh sdmTMB can handle
# Besides mesh, input data need to be data.frame format, not sf
mesh2 <- make_mesh(hi_2018_df, c("longitude", "latitude"), mesh = mesh)
plot(mesh2)

# Fit a Gaussian GLMM with spatial random field
fit1 <- sdmTMB(
  shoots_per_m ~ 1,
  data = hi_2018_df,
  mesh = mesh2,
  spatial = "on"
)

# we can use the sanity() function to check whether we have a well-fitting model
sanity(fit1)

# another method to check results is via residual simulation
simulate(fit1, nsim = 500) %>% 
  dharma_residuals(fit1)

# Nope!

# Shoot densities have many zeros and can't be <0. Need to generalize response family


# Spatial model fitting 2----
# Fit a Tweedie GLMM with spatial random field
fit2 <- sdmTMB(
  shoots_per_m ~ 1,
  data = hi_2018_df,
  mesh = mesh2,
  family = tweedie(),
  spatial = "on"
)

# Check again:
sanity(fit2)
simulate(fit2, nsim = 500) %>% 
  dharma_residuals(fit2)

# Much better!!!

# Model visualization----

# Residuals over space----
hi_2018_df$res <- residuals(fit2)

# Clusters of residuals of similar magnitude means spatial correlation wasn't
# dealt with effectively by the spatial random effect. This looks pretty good
ggplot(hi_2018_df) +
  geom_sf(data = hi_sg_adjusted) +
  geom_point(aes(y = latitude, x = longitude, 
                 size = res, color = res))

# Generate a shoot density surface----
pred_out <- 
  predict(fit2, newdata = pred_df) %>% 
  mutate(est = exp(est),
         est_rf = exp(est_rf))

ggplot(pred_out) +
  geom_sf(data = hi_sg_adjusted) +
  geom_raster(data = pred_out, aes(y = latitude, x = longitude, fill = est)) +
  geom_sf(data = hi_sg, fill = "transparent", color = "red") +
  scale_fill_viridis_c() +
  labs(fill = "Shoot density (shoot m^-2)") +
  theme_bw()

ggplot(pred_out) +
  geom_sf(data = hi_sg_adjusted) +
  geom_raster(data = pred_out, aes(y = latitude, x = longitude, fill = est)) +
  scale_fill_viridis_b()


