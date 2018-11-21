### This is for plotting the results of dis~haz_analysis 
### written by J.Sebald & C. Senf

# load packages, set working directory to either local (pc/macbook) or server (valenor)

{
  library(raster)# version 2.6-7
  library(sp)# version 1.3-1
  library(rgdal)# version 1.2-20
  library(igraph)# version 1.2.1
  library(tidyverse)# version 1.2.1
  library(projpred)
  library(ggthemes)
  library(rstanarm)
  library(sf)
  library(ggridges)
  library(gridExtra)
  
  rm(list=ls())
  
}



# Load data ---------------------------------------------------------------

load(file = "../../../../../results/results.RData")

data <- read.csv("../data/data_for_model.csv")

shp_ws <- raster::shapefile("../data/shp_ws.shp")
shp_ws_fortify <- broom::tidy(shp_ws, region = "WLK_ID")

# Set names ---------------------------------------------------------------

processes <- c("DFLOOD", "DFLOW", "FST")

vars_ws <- data.frame(varname = c("h_mean", "Melton", "Elevation", "Circularit", "Elongation", "artifical", "forest", "area", "patchdensity", 
                                  "severity", 
                                  "frequency",
                                  "severity:frequency"), 
                      name = c("Elevation", "Melton ratio", "Elevation ratio", "Circularity", "Elongtion", "Artificial", "Forest", "Area", "Patch density", 
                               "Severity", 
                               "Frequency",
                               "Severity x Frequency"),
                      stringsAsFactors = FALSE)

# Plotting ----------------------------------------------------------------

subset_size <- results %>% 
  map(~ .[[2]])

p_varsel <- results %>% 
  map(~ .[[4]]$data) %>%
  set_names(processes) %>%
  bind_rows(.id = "process") %>%
  ggplot(., aes(x = size, y = value)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = value - se, ymax = value + se), width = 0) +
  facet_wrap(~ process, scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank()) +
  labs(x = "Number of predictors", y = "Expected log predictive density") +
  geom_vline(data = data.frame(process = processes, n_select = unlist(subset_size)),
             aes(xintercept = n_select), col = scales::muted("red"), linetype = "dashed")

ggsave("variable_selection.pdf", p_varsel, path = "results/", width = 7.5, height = 2.5)
ggsave("variable_selection.png", p_varsel, path = "results/", width = 7.5, height = 2.5)


predictors_selected <- results %>% 
  map(~ .[[3]]$name) %>%
  map(~ paste(., collapse = ", "))

elpd <- results %>% 
  map(~ .[[6]])

auc <- results %>% 
  map(~ .[[7]])

watershed_model_results_table <- data.frame(process = processes,
                                            subset_size = unlist(subset_size),
                                            predictors_selected = unlist(predictors_selected),
                                            elpd = unlist(elpd),
                                            auc = unlist(auc))

write_csv(watershed_model_results_table, "results/watershed_model_results_table.csv")

p_roc <- results %>% 
  map(~ .[[8]]$data) %>%
  set_names(processes) %>%
  bind_rows(.id = "process") %>%
  ggplot(., aes(x = fpr, y = tpr)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", col = scales::muted("red")) +
  facet_wrap(~ process, scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank()) +
  labs(x = "False positive rate", y = "True positive rate")

ggsave("roc_curves.pdf", p_roc, path = "results/", width = 7.5, height = 2.5)
ggsave("roc_curves.png", p_roc, path = "results/", width = 7.5, height = 2.5)


p_estimate <- results %>% 
  map(~ .[[9]]) %>%
  bind_rows() %>%
  ggplot(., aes(x = fct_rev(name), y = value)) +
  geom_violin(fill = "grey") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank()) +
  coord_flip() +
  theme(strip.background = element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed", col = scales::muted("red")) +
  labs(y = "Posterior probability distribution of parameter estimates", x = NULL, fill = "Process") +
  scale_fill_brewer(palette = "Greys", direction = -1) +
  facet_wrap(~process)

ggsave("estimates.pdf", p_estimate, path = "results/", width = 7.5, height = 2.5)
ggsave("estimates.png", p_estimate, path = "results/", width = 7.5, height = 2.5)


p_maps <- results %>% 
  map(~ .[[10]]$data) %>%
  set_names(processes) %>%
  bind_rows(.id = "process") %>%
  ggplot(.) +
  aes(long, lat, group = group, fill = prob_hazard_mean) +
  geom_polygon() +
  coord_equal() +
  scale_fill_gradient(low = "#fff5f0", high = "#67000d") +
  labs(x = NULL, y = NULL, fill = paste0("Probability of event")) +
  ggthemes::theme_map() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        strip.background = element_blank()) +
  guides(fill = guide_colorbar(barwidth = 10, barheight = 0.5,
                               direction = "horizontal", title.position = "top")) +
  facet_wrap(~process)

ggsave("probability_maps.pdf", p_maps, path = "results/", width = 7.5, height = 2.5)
ggsave("probability_maps.png", p_maps, path = "results/", width = 7.5, height = 2.5)


p_maps_difference <- results %>% 
  map(~ .[[11]]$data) %>%
  set_names(processes) %>%
  bind_rows(.id = "process") %>%
  ggplot(.) +
  aes(long, lat, group = group, fill = prob_hazard_mean) +
  geom_polygon() +
  coord_equal() +
  scale_fill_gradient2() +
  labs(x = NULL, y = NULL, fill = paste0("Change in event probability after including disturbances")) +
  ggthemes::theme_map() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        strip.background = element_blank()) +
  guides(fill = guide_colorbar(barwidth = 10, barheight = 0.5,
                               direction = "horizontal", title.position = "top")) +
  facet_wrap(~process)

ggsave("probability_maps_difference.pdf", p_maps_difference, path = "results/", width = 7.5, height = 2.5)
ggsave("probability_maps_difference.png", p_maps_difference, path = "results/", width = 7.5, height = 2.5)


p_response_classic <- results %>% 
  map(~ .[[12]]$data) %>%
  set_names(processes) %>%
  bind_rows(.id = "process") %>%
  #mutate(frequency = factor(frequency, labels = c("High (+1SD)", "Average", "Low (-1SD)"))) %>%
  ggplot(., aes(x = severity, y = prob_mean)) +
  geom_ribbon(aes(ymin = prob_mean - prob_sd * 2, ymax = prob_mean + prob_sd * 2, fill = frequency), alpha = 0.3) +
  geom_line(aes(col = frequency)) +
  geom_point(data = sample_n(data %>% mutate(severity = as.double(scale(severity))) %>% filter(severity < quantile(severity, 0.99)), 1000), 
             aes(x = severity, y = -0.05), shape = 124, alpha = 0.3) +
  theme_bw() +
  theme(#legend.position = c(0, 1),
    #legend.justification = c(0, 1),
    legend.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank()) +
  labs(x = "Disturbance severity", y = paste0("Probability of event"), 
       col = "Disturbance frequency", fill = "Disturbance frequency") +
  scale_color_manual(values = c(scales::muted("blue"), "grey", scales::muted("red"))) +
  scale_fill_manual(values = c(scales::muted("blue"), "grey", scales::muted("red"))) +
  facet_wrap(~process)

ggsave("responsecurves_classic.pdf", p_response_classic, path = "results/", width = 7.5, height = 2.25)
ggsave("responsecurves_classic.png", p_response_classic, path = "results/", width = 7.5, height = 2.25)



p_response_heatmap <- results %>% 
  map(~ .[[13]]$data) %>%
  set_names(processes) %>%
  bind_rows(.id = "process") %>%
  ggplot(., aes(x = severity, y = frequency * (-1))) +
  geom_tile(aes(fill = prob_mean)) +
  # geom_point(data = data %>% mutate(severity = as.double(scale(severity)),
  #                                   frequency = as.double(scale(frequency))) %>%
  #              sample_n(1000), 
  #                 aes(x = severity, y = frequency * (-1)), alpha = 0.3) +
  scale_fill_continuous(low = "#fff5f0", high = "#67000d", limit = c(0, 1)) +
  theme_bw() +
  theme(legend.position = c(0, 0),
        legend.justification = c(-0.05, -0.1),
        legend.box.background = element_rect(colour = "black"),
        legend.margin = margin(2, 6, 2, 6),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        panel.grid = element_blank(),
        strip.background = element_blank()) +
  labs(x = "Disturbance severity", fill = paste0("Probability of event"), 
       y = "Disturbance frequency") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill = guide_colorbar(barwidth = 6, barheight = 0.5,
                               direction = "horizontal", title.position = "top")) +
  facet_wrap(~process)


### Correlations

prob_count <- results %>%
  map(~ posterior_linpred(.[[1]], transform = TRUE), re.form = NA) %>%
  map2(.y = list(data$DFLOOD, data$DFLOW, data$FST),
       ~ data.frame(pred = apply(.x, 2, mean), count = .y)) %>%
  set_names(processes) %>%
  bind_rows(.id = "process") %>%
  ggplot(., aes(x = factor(count), y = pred)) +
  geom_violin(fill = "grey") +
  facet_wrap(~process)


ggsave("correlations.pdf", p_response_classic, path = "results/", width = 7.5, height = 2.25)
ggsave("correlations.png", p_response_classic, path = "results/", width = 7.5, height = 2.25)



# temporal trend in natural hazard events -------------------------------------------------------

events <- read.csv("../../../../../materials/raw_data/data_natural_hazards/GroupEvents.csv")


summary(events)

data <- events %>%
  dplyr::select(prozessart, ereignis_jahr) %>%
  rename(year = ereignis_jahr, process = prozessart) %>%
  filter(year %in% 1984:2018) %>%
  mutate(process = case_when(
    process == "Fluviatiler Feststofftransport" ~ "FST",
    process ==  "Murartiger Feststofftransport" ~ "DFLOOD",
    process == "Murgang" ~ "DFLOW")) %>%
  group_by(year, process) %>%
  summarize(events = n()) %>%
  mutate(process = factor(process, levels = c("DFLOW","DFLOOD","FST")))


p_trend <- ggplot(data, aes(x = year, y = events, col = process)) +
  geom_point() +
  labs(title = "", col = "") +
  geom_smooth(method = "gam", formula = y ~ s(x)) +
  facet_wrap(~process, scales = "free") +
  theme_few() +
  theme(legend.position = "none",
        axis.title = element_blank()) 

ggsave("p_trend.pdf", p_trend, path = "../../../../../results/figures/", width = 6, height = 3)
ggsave("p_trend.png", p_trend, path = "../../../../../results/figures/", width = 6, height = 3)




# Maps of study area --------------------------------------------------------------------

library(maptools)
library(maps)

# study area

austria <- read_sf("../../../../../materials/raw_data/shapfiles_austria/austria_EPSG_3035_ETRS89.shp") %>%
  mutate(country = "austria") %>%
  group_by(country) %>%
  summarize()

studyarea <- read_sf("../data/shapes/shp_ws.shp") %>%
  mutate(country = "austria") %>%
  group_by(country) %>%
  summarize()

p_study_area <- ggplot(austria) +
  geom_sf(aes(fill = " Austria")) +
  geom_sf(data = studyarea, aes(fill = " Study Area")) +
  scale_fill_manual(values = c("gray82", "gray42")) +
  theme_bw() +
  theme(legend.title = element_blank())


ggsave("../../../../../results/figures/p_study_area.png", p_study_area, width = 8, height = 5)


worldmap <- rnaturalearth::ne_download(scale = 50,
                                       type = "countries",
                                       category = "cultural",
                                       destdir = tempdir(),
                                       load = TRUE,
                                       returnclass = "sf") %>%
  mutate(austria = ifelse(NAME == "Austria", 1, 0))
  


ggplot(worldmap) +
  geom_sf(aes(fill = factor(austria))) +
  ylim(30, 70) +
  xlim(-20, 40) +
  theme_bw() +
  scale_fill_manual(values = c("gray82", "gray42")) +
  theme(legend.position = "none")



overview_map <- ggplot(worldmap) +
  geom_sf(aes(fill = factor(austria))) +
  geom_sf(data = studyarea, aes(fill = "Study area")) +
  ylim(37, 70) +
  xlim(-11, 35) +
  theme_bw() +
  scale_fill_manual(labels = c("Europe", "Austria", "Study area"), values = c("gray94","gray65", "gray42")) +
  guides(fill = guide_legend(title = "")) 

ggsave("../../../../../results/figures/overview_map.png", overview_map, width = 7, height = 5)
ggsave("../../../../../results/figures/overview_map.pdf", overview_map, width = 8, height = 6, dpi = 300)

 



# boxplots of disturbance metrics --------



data <- read.csv("../data/tables/data_for_model.csv") %>%
  mutate(process = case_when(
    FST >= 1 ~ "FST",
    DFLOOD >= 1  ~ "DFLOOD",
    DFLOW >= 1 ~ "DFLOW", 
    TRUE ~ "NO EVENT"),
    event = ifelse(process != "NO EVENT" , "yes", "no")) %>%
  mutate(process = factor(process, levels = c("NO EVENT","DFLOW","DFLOOD","FST"))) %>%
  mutate(press = 1 - pulse) %>%
  mutate(mean_extent = (extent/32)*100) %>%
  mutate(sum = DFLOW+DFLOOD+FST) %>%
  mutate(intens = case_when(sum > 5 ~ "> 5",
                           sum == 0 ~ "0",#
                           TRUE ~ "1 - 5"))


ggplot(data, aes(x = factor(process), y = press)) +
  geom_boxplot(fill = c("#d9f0d3","#a6dba0", "#5aae61", "#1b7837" )) +
  theme_bw() +
  labs(y = "Disturbance press") +
  theme(axis.title.x = element_blank())

ggplot(data, aes(x = factor(process), y = mean_extent)) +
  geom_boxplot(fill = c("#d9f0d3","#a6dba0", "#5aae61", "#1b7837" )) +
  theme_bw() +
  labs(y = "% forest area disturbed / year") +
  theme(axis.title.x = element_blank())




# Point plotsm boxplots and ridges of forest related predictors  -------------------------------


data <- read.csv("../data/tables/data_for_model.csv") %>%
  mutate(process = case_when(
    FST >= 1 ~ "FST",
    DFLOOD >= 1  ~ "DFLOOD",
    DFLOW >= 1 ~ "DFLOW", 
    TRUE ~ "NO EVENT"),
    event = ifelse(process != "NO EVENT" , "yes", "no")) %>%
  mutate(process = factor(process, levels = c("NO EVENT","DFLOW","DFLOOD","FST"))) %>%
  mutate(press = 1 - pulse) %>%
  mutate(mean_extent = (extent/32)*100) %>%
  mutate(sum = DFLOW+DFLOOD+FST) %>%
  mutate(intens = case_when(sum > 5 ~ "> 5",
                            sum == 0 ~ "0",#
                            TRUE ~ "1 - 5"))

data_long <- data %>%
  select(WLK_ID, forest, patchdensity, extent, press, intens) %>%
  mutate_at(.vars = vars(forest:press), function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)) %>%
  gather(key = predictor, value = value, forest:press) %>%
  mutate(predictor = factor(predictor, levels = c("forest", "patchdensity", "extent", "press"))) %>%
  mutate(intens = factor(intens, levels = c("0", "1 - 5", "> 5")))

data_group <- data_long %>%
  group_by(predictor, intens) %>%
  summarize(mea = mean(value),
            se = sd(value)) %>%
  ungroup() %>%
  mutate(predictor = factor(predictor, levels = c("forest", "patchdensity", "extent", "press"))) %>%
  mutate(intens = factor(intens, levels = c("> 5", "1 - 5", "0")))


forest_and_disturbances_point <- ggplot(data_group, aes(x = mea, y = predictor, colour = intens)) +
         geom_point(size = 6, shape = 17) +
  labs(y = "", x = "") +
  theme_bw() +
  guides(col = guide_legend(title = "Hazards")) +
  scale_colour_manual(values = c("#c51b7d", "#f1b6da","#7fbc41")) +
  scale_y_discrete(labels = c("Disturbance Press", "Canopy disturbed [%]", "Forest patches [n/km²]", "Forest [%]"))
 

  
ggsave("../../../../../results/figures/forest_and_disturbances_point.png", forest_and_disturbances_point, width = 6, height = 3)
ggsave("../../../../../results/figures/forest_and_disturbances_point.pdf", forest_and_disturbances_point, width = 6, height = 3)

forest_and_disturbances_boxplot <- ggplot(data_long, aes(x = predictor, y = value, fill = intens)) +
  geom_boxplot() +
  labs(y = "", x = "") +
  ylim(-2.5, 2.5) +
  theme_bw() +
  guides(fill = guide_legend(title = "Hazards")) +
  scale_fill_manual(values = c("#7fbc41", "#f1b6da", "#c51b7d")) +
  scale_x_discrete(labels = c("Forest share [%]","Forest patches [n/km²]", "Canopy disturbed [%]", "Disturbance Press")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave("../../../../../results/figures/forest_and_disturbances_boxplot.png", forest_and_disturbances_boxplot, width = 6, height = 4)
ggsave("../../../../../results/figures/forest_and_disturbances_boxplot.pdf", forest_and_disturbances_boxplot, width = 6, height = 3)


forest_and_disturbances_ridges <- ggplot(data_long, aes(x = value, y = predictor,  fill = intens)) +
  geom_density_ridges(scale = 0.95, alpha = 0.8) +
  labs(y = "", x = "") +
  theme_bw() +
  xlim(-3, 4)+
  guides(fill = guide_legend(title = "Hazards")) +
  scale_fill_manual(values = c("#c51b7d", "#f1b6da","#7fbc41")) +
  scale_y_discrete(labels = c("Disturbance Press","Canopy disturbed [%]",  "Forest patches [n/km²]", "Forest [%]"))

ggsave("../../../../../results/figures/forest_and_disturbances_ridges.png", forest_and_disturbances_ridges, width = 6, height = 3)
ggsave("../../../../../results/figures/forest_and_disturbances_ridges.pdf", forest_and_disturbances_ridges, width = 9, height = 6)



# Point plot of final predictors ------------------------------



p_predictors_points <- read.csv("../data/tables/data_for_model.csv") %>%
  select(WLK_ID,area, h_mean, Elevation, Melton, Elongation, Circularit, forest, patchdensity, extent, pulse) %>%
  mutate(press = 1 - pulse,
         forest = forest*100,
         extent = extent*100) %>%
  select(- pulse) %>%
  rename("Area [km²]" = area, 
         "Elevation [m]" = h_mean, 
         "Elevation ratio" = Elevation,
         "Meltion ratio" = Melton,
         Circularity = Circularit,
         "Forest cover [%]" = forest,
         "Forest patches [n/km²]" = patchdensity,
         "Forest area disturbed [%]" = extent,
         "Disturbance press" = press) %>%
  gather(key = predictor, value = value, "Area [km²]":"Disturbance press") %>%
  group_by(predictor) %>%
  summarize(mean = mean(value),
            min = min(value),
            max = max(value)) %>%
  gather(key = term, value = value, mean:max) %>%
  ggplot(., aes(x = value, y = "", colour = term)) +
  geom_point(shape = 17, size = 4) +
  facet_wrap( ~predictor, scales = "free") +
  labs(y = "", x = "") +
  theme_few() +
  guides(col = guide_legend(title = "")) +
  scale_colour_manual(values = c("#c51b7d","#7fbc41", "#f1b6da"))


ggsave("../../../../../results/figures/p_predictors_points.png", p_predictors_points, width = 8, height = 6)



# Violoin plot of final predictors -----------------------------------------


p_predictors_violin_gen <- read.csv("../data/tables/data_for_model.csv") %>%
  dplyr::select(WLK_ID,area, h_mean) %>%
  rename("Area [km²]" = area, 
         "Mean elevation [m]" = h_mean) %>%
  gather(key = predictor, value = value, "Area [km²]":"Mean elevation [m]") %>%
  mutate(predictor = factor(predictor, levels = c("Area [km²]", "Mean elevation [m]"))) %>%
  ggplot(., aes(x = "", y = value)) +
  geom_violin(fill = "#ffffbf") +
  facet_wrap(~ predictor, scales = "free") +
  theme_few() +
  labs(y = "", x = "", title = "General")  +
  guides(fill = FALSE)


p_predictors_violin_geo <- read.csv("../data/tables/data_for_model.csv") %>%
  dplyr::select(WLK_ID, Elevation, Melton, Elongation, Circularit) %>%
  rename("Elevation ratio" = Elevation,
         "Meltion ratio" = Melton,
         "Circularity" = Circularit) %>%
  gather(key = predictor, value = value, "Elevation ratio":"Circularity") %>%
  mutate(predictor = factor(predictor, levels = c("Elevation ratio", "Meltion ratio", "Elongation", "Circularity"))) %>%
  ggplot(., aes(x = "", y = value)) +
  geom_violin(fill = "#4393c3") +
  facet_wrap(~ predictor, scales = "free") +
  theme_few() +
  labs(y = "", x = "", title = "Geomorphological")  +
  guides(fill = FALSE)

p_predictors_violin_forest <- read.csv("../data/tables/data_for_model.csv") %>%
  dplyr::select(WLK_ID, forest, patchdensity, extent, pulse) %>%
  mutate(press = 1 - pulse,
         forest = forest*100,
         extent = extent*100) %>%
  dplyr::select(- pulse) %>%
  rename("Forest cover [%]" = forest,
         "Forest patches [n/km²]" = patchdensity,
         "Canopy disturbed [%]" = extent,
         "Disturbance press" = press) %>%
  gather(key = predictor, value = value, "Forest cover [%]":"Disturbance press") %>%
  mutate(predictor = factor(predictor, levels = c("Forest cover [%]", "Forest patches [n/km²]", "Canopy disturbed [%]", "Disturbance press" ))) %>%
  ggplot(., aes(x = "", y = value)) +
  geom_violin(fill = "#276419") +
  facet_wrap(~ predictor, scales = "free") +
  theme_few() +
  labs(y = "", x = "", title = "Forest related")  +
  guides(fill = FALSE)


p_predictors_arranged <- grid.arrange( p_predictors_violin_gen, p_predictors_violin_geo, p_predictors_violin_forest, ncol = 3)


ggsave("../../../../../results/figures/p_predictors_arranged.pdf", p_predictors_arranged, width = 11.5, height = 5, dpi = 300)
ggsave("../../../../../results/figures/p_predictors_arranged.png", p_predictors_arranged, width = 11, height = 5)




# Numbers of final dataset ------------------------------------------------


data <- read.csv("../data/tables/data_for_model.csv") 

study_area_km2 <- sum(data$area)
austria_km2 <- 83879

study_area_km2/austria_km2 *100

forest_km2 <- study_area_km2*(mean(data$forest))

austria_forest <- 40000

forest_km2/austria_forest

forest_km2*mean(data$extent)

events <- sum(data$DFLOW) + sum(data$DFLOOD) + sum(data$FST)

watersheds <- nrow(data)

 