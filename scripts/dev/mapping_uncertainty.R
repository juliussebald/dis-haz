
#install.packages("colorspace", repos = "http://R-Forge.R-project.org")
#devtools::install_github("clauswilke/multiscales")

library(ggplot2)
library(multiscales)
library(tidyverse)
library(sf)

load("/Users/corneliussenf/Downloads/for_cornelius.RData")

shape <- sf::read_sf("/Users/corneliussenf/Dropbox/Research/European Forests/Hazards/shp_ws.shp")

pred <- data_frame(WLK_ID = as.character(r$data$WLK_ID),
                   pred = apply(r$predictions_all, 2, mean),
                   sd = apply(r$predictions_all, 2, sd)) %>%
  mutate(., varcof = sd / pred)

shape <- shape %>%
  left_join(pred, by = "WLK_ID")

# shp <- r$shp_fortiy %>%
#   mutate(WLK_ID = as.integer(id)) %>%
#   left_join(pred, by = "WLK_ID")

colors <- scales::colour_ramp(
  colors = c(blue = "#2265A3", purple = "#740280", red = "#AC202F")
)((0:7)/7)

pred_cutoff <- quantile(shape$pred, 0.95, na.rm = TRUE)
varcof_cutoff <- quantile(shape$varcof, 0.95, na.rm = TRUE)

shape %>%
  mutate(pred = ifelse(pred >= pred_cutoff, pred_cutoff, pred),
         varcof = ifelse(varcof >= varcof_cutoff, varcof_cutoff, varcof)) %>%
  ggplot(.) + 
  geom_sf(aes(fill = zip(pred, varcof)), color = NA, size = 0.2) + 
  coord_sf(datum = NA) +
  bivariate_scale("fill",
                  pal_vsup(values = colors, max_desat = 0.6, pow_desat = 1, max_light = 0.5, pow_light = 1),
                  name = c("P(Event)", "uncertainty"),
                  limits = list(quantile(shape$pred, c(0, 0.95), na.rm = TRUE), 
                                quantile(shape$varcof, c(0, 0.95), na.rm = TRUE)),
                  breaks = list(quantile(shape$pred, c(0, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE), 
                                quantile(shape$varcof, c(0, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)),
                  guide = "colourfan"
  ) +
  theme_void() +
  theme(
    legend.key.size = grid::unit(0.8, "cm"),
    legend.title.align = 0.5,
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  )
