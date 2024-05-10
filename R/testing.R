library(tidyverse)
library(sf)
library(sdmTMB)
library(sp)
library(gstat)
library(rnaturalearth)
library(raster)

ebs <- st_read(here::here("data/ebs_friendly.kml")) %>% 
  st_transform(., crs = crabber::ak_crs)

shoreline <- 
  ne_states(country = "United States of America") %>% 
  filter(name == "Alaska") %>% 
  st_union() %>% 
  st_transform(., crs = crabber::ak_crs) %>% 
  st_intersection(., ebs) 



ex1_df <- crabber::cpue_bbrkc %>% 
  filter(mat_sex == "legal male",
         year %in% 2010) %>% 
  mutate(station = factor(station))

# visualize across longitude
ggplot(ex1_df) +
  geom_point(aes(y = lat, x = lon,
                 color = station),
             show.legend = F) +
  geom_sf(data = shoreline) + 
  coord_sf(ylim = c(550, 1000),
           xlim = c(-900, -150))

# Example 1 (GLM) ------

# fit intercept-only Tweedie GLM
m1 <- sdmTMB(cpue ~ 1,
             spatial = "off",
             data = ex1_df,
             family = tweedie())

# extract predicted mean density
density_prediction <- 
  tidy(m1) %>% 
    mutate(estimate_resp = exp(estimate),
           upr = exp(estimate + 2 * std.error ),
           lwr = exp(estimate - 2 * std.error))

# same as predict(m1) %>% mutate(est = exp(est)) %>% ...

# visualize prediction
ggplot(density_prediction) + 
  geom_histogram(data  = ex1_df,
                 aes(cpue)) +
  geom_vline(aes(xintercept = estimate_resp),
             color = "purple",
             linewidth = 2) +
  geom_vline(aes(xintercept = lwr),
             color = "purple",
             linewidth = 1,
             linetype = 2) +
  geom_vline(aes(xintercept = upr),
             color = "purple",
             linewidth = 1,
             linetype = 2) +
  labs(y = "Frequency",
       x = "RKC per sq. km") +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(axis.title = element_text(size = 14)) +
  theme_minimal()

# extract standardized residuals
ex1_df$residuals <- residuals(m1, type = "mle-eb")

# visualize spatial patterning of residuals
ggplot(ex1_df) + 
  geom_point(aes(y = lat, x = lon,
                 color = residuals),
             size = 2.5,
             alpha = 0.8) +
  scale_color_gradient2(midpoint = 0) +
  dream::theme_fade()

# Evaluate residuals spatial autocorrelation using a variogram
ex1_df_sp <- ex1_df
coordinates(ex1_df_sp) <- c("lon", "lat")

v1 <- variogram(object = residuals ~ 1,
          data = ex1_df_sp,
          alpha = 290,
          cressie = T,
          cutoff = 1000)

ggplot(v1) + 
  geom_point(aes(y = gamma,
                 x = dist)) +
  labs(y = "Sample Variogram",
        x = "Distance (km)")

# Example 2 (Individual level GLMM)

m2 <- sdmTMB(cpue ~ 1 + (1|station),
             spatial = "off",
             data = ex1_df,
             family = tweedie())

# extract standardized residuals
ex1_df$residuals2 <- residuals(m2, type = "mle-eb")

# visualize spatial patterning of residuals
ggplot(ex1_df) + 
  geom_point(aes(y = lat, x = lon,
                 color = residuals2),
             size = 2.5,
             alpha = 0.8) +
  scale_color_gradient2(midpoint = 0)

ggplot(ex1_df) + 
  geom_point(aes(y = lat, x = lon,
                 color = residuals),
             size = 2.5,
             alpha = 0.8) +
  scale_color_gradient2(midpoint = 0)

density_prediction2 <- 
  tidy(m2) %>% 
  mutate(estimate_resp = exp(estimate),
         upr = exp(estimate + 2 * std.error ),
         lwr = exp(estimate - 2 * std.error))

# same as predict(m1) %>% mutate(est = exp(est)) %>% ...

# visualize prediction
ggplot(density_prediction) + 
  geom_histogram(data  = ex1_df,
                 aes(cpue)) +
  geom_vline(aes(xintercept = estimate_resp),
             color = "purple",
             linewidth = 2) +
  geom_vline(aes(xintercept = lwr),
             color = "purple",
             linewidth = 1,
             linetype = 2) +
  geom_vline(aes(xintercept = upr),
             color = "purple",
             linewidth = 1,
             linetype = 2) +
  
  geom_vline(data = density_prediction2,
             aes(xintercept = estimate_resp),
             color = "red",
             linewidth = 2) +
  geom_vline(data = density_prediction2,
             aes(xintercept = lwr),
             color = "red",
             linewidth = 1,
             linetype = 2) +
  geom_vline(data = density_prediction2,
             aes(xintercept = upr),
             color = "red",
             linewidth = 1,
             linetype = 2) +
  labs(y = "Frequency",
       x = "RKC per sq. km") +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme(axis.title = element_text(size = 14)) +
  theme_minimal()



predict(m2, type = "response") %>% 
  ggplot() +
    geom_point(aes(y = cpue, x = station)) +
    geom_point(aes(y = est, x = station), color = "red") + 
    theme(axis.text = element_text(angle = 90))


# Residual autocorrelation is resolved via the incorporation
# of station-level random intercept terms

# Example 3 (Spatial GLMM) -----
df3 <- crabber::cpue_bbrkc %>% 
  filter(mat_sex == "legal male",
         year %in% 2010) %>% 
  mutate(station = factor(station))

mesh <- make_mesh(df3, c("lon","lat"), cutoff = 20)
plot(mesh)

# fit intercept-only Tweedie GLMM
m3 <- sdmTMB(cpue ~ 1,
             spatial = "on",
             mesh = mesh,
             data = df3,
             anisotropy = F,
             family = tweedie())

m3

# extract standardized residuals
df3$residuals <- residuals(m3, type = "mle-eb")

# visualize spatial patterning of residuals
ggplot(df3) + 
  geom_point(aes(y = lat, x = lon,
                 color = residuals),
             size = 2.5,
             alpha = 0.8) +
  scale_color_gradient2(midpoint = 0) +
  dream::theme_fade()

# spatial anisotropy - when the range of the correlation differs directionally


# Evaluate residuals spatial autocorrelation using a variogram
coordinates(df3) <- c("lon", "lat")

v3 <- variogram(object = residuals ~ 1,
                data = df3,
                cressie = T,
                cutoff = 1000)
# 
# ggplot(v2) + 
#   geom_point(aes(y = gamma,
#                  x = dist)) +
#   geom_point(data = v1,
#              aes(y = gamma,
#                  x = dist), color = "red") +
ggplot()+
  geom_point(data = v3,
             aes(y = gamma,
                 x = dist), color = "blue") +
  labs(y = "Sample Variogram",
       x = "Distance (km)")
