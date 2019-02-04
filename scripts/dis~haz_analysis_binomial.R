
# Load packages and set options -------------------------------------------

{
  Sys.setenv(USE_CXX14 = 1)
  library(raster) # version 2.6-7
  library(sp) # version 1.3-1
  library(rgdal) # version 1.3-4
  library(igraph) # version 1.2.2
  library(tidyverse) # version 1.2.1
  library(rstanarm) # version 2.17.4
  library(projpred )# version 0.8.0
  library(multiscales) # version 0.1.0 # devtools::install_github("clauswilke/multiscales")
  library(patchwork) # version 0.0.1
  library(rstan) # version 2.17.4
  library(sf) # version 0.6-3
  library(gridExtra) # version 2.3
  
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  
}

# Load data ---------------------------------------------------------------

data <- read.csv("../data/data_for_model.csv") 

processes <- c("MFL", "FST")
processes_names <- c("Debris-flow", "Sediment-transport")

# Selecting a model -----------------------------------------------

# Create dataframe with variable short names and variable long names

vars_ws <- data.frame(varname = c("h_mean", "Melton", "Elevation", "Circularit", "Elongation", 
                                  "artifical", "forest", "area", "patchdensity", 
                                  "eco_region", 
                                  "extent", 
                                  "type",
                                  "extent:type"), 
                      name = c("Elevation", "Melton ratio", "Elevation ratio", "Circularity", "Elongation", 
                               "Artificial", "Forest", "Area", "Patch density", 
                               "Ecological region", 
                               "Extent", 
                               "Type",
                               "Extent x Type"),
                      stringsAsFactors = FALSE)

# Loop through processes and calibrate varying models

models <- vector("list", length = length(processes))

k <- 0

for (process in processes) {
  
  k <- k + 1
  
  # Bring data into form
  
  vars_nointeraction <- vars_ws %>% 
    filter(!varname %in% c("extent:type", "eco_region"))
  
  data_model <- data
  
  # make sure that have not excpierienced any disturbances have NA as "disturbance type"
  
  data_model[data_model$extent == 0, "type"] <- NA
  
  # z-transform predictors
  
  data_model <- data_model %>%
    mutate_at(.vars = vars(c(vars_nointeraction$varname)), function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  
  # devide dataset in watershed which have had a torrential hazard event (1) and watersheds which had no torrential hazard event (0)
  
  data_model[, "response"] <- ifelse(data_model[, process] > 0, 1, 0)
  
  # make sure that watersheds which have not expierienced any disturbances are set to disturbance type "avarage"
  
  data_model[data_model$extent == min(data_model$extent), "type"] <- 0
  
  # transform eco region predictor to factor
  
  data_model$eco_region <- as.factor(data_model$eco_region)
  
  # Fit null model
  
  fit_null <- stan_glm(response ~ 1,
                       data = data_model,
                       family = binomial,
                       prior = normal(0, 0.5, autoscale = TRUE),
                       prior_intercept = normal(0, 0.5, autoscale = TRUE),
                       prior_aux = exponential(rate = 1, autoscale = TRUE))
  
  # Calculate loo for null model
  
  loo_fit_null <- loo(fit_null)
  
  # Fit general model
  
  fit_general <- update(fit_null, . ~ . + area + h_mean + artifical + eco_region, QR = TRUE)
  
  # Calculate loo for general model
  
  loo_fit_general <- loo(fit_general)
  
  # Fit general+geomorph model
  
  fit_general_geomorph <- update(fit_general, . ~ . + Elevation + Melton + Circularit + Elongation, QR = TRUE)
  
  # Calculate loo for general+geomorph model
  
  loo_fit_general_geomorph <- loo(fit_general_geomorph)
  
  # Fit general+geomorph+forest model
  
  fit_general_geomorph_forest <- update(fit_general_geomorph, . ~ . + forest + patchdensity, QR = TRUE)
  
  # Calculate loo for general+geomorph+forest model
  
  loo_fit_general_geomorph_forest <- loo(fit_general_geomorph_forest)
  
  # Fit general+geomorph+forest-disturbance model
  
  fit_general_geomorph_forest_disturbances <- update(fit_general_geomorph_forest, . ~ . + extent * type)
  
  # Calculate loo for full model
  
  loo_fit_general_geomorph_forest_disturbances <- loo(fit_general_geomorph_forest_disturbances)
  
  # compare watershed model and full model
  
  model_comparison <- loo::compare(loo_fit_null,
                                   loo_fit_general,
                                   loo_fit_general_geomorph,
                                   loo_fit_general_geomorph_forest,
                                   loo_fit_general_geomorph_forest_disturbances)
  
  model_comparison_direct <- list(loo::compare(loo_fit_null,
                                               loo_fit_general),
                                  loo::compare(loo_fit_general,
                                               loo_fit_general_geomorph),
                                  loo::compare(loo_fit_general_geomorph,
                                               loo_fit_general_geomorph_forest),
                                  loo::compare(loo_fit_general_geomorph_forest,
                                               loo_fit_general_geomorph_forest_disturbances))
  
  model_comparison_direct <- model_comparison_direct %>%
    map(~ as.vector(.)) %>%
    set_names(c("General",
                "Geomorph",
                "Forest",
                "Disturbances")) %>%
    bind_rows() %>%
    t() %>%
    as.data.frame() %>%
    rename_("ELPD difference" = "V1", "SE" = "V2")
    
    
  # Store everything in a list
  
  models[[k]] <- list(list(fit_null,
                           fit_general,
                           fit_general_geomorph,
                           fit_general_geomorph_forest,
                           fit_general_geomorph_forest_disturbances),
                      list(loo_fit_null,
                           loo_fit_general,
                           loo_fit_general_geomorph,
                           loo_fit_general_geomorph_forest,
                           loo_fit_general_geomorph_forest_disturbances),
                      model_comparison,
                      model_comparison_direct)
                       
}

save(models, file = "../results/binomial/models_binomial.RData")

load(file =  "../results/binomial/models_binomial.RData")

# Model evaluation -------------------------------------------------------------

### Model comparison

elpd <- models %>%
  map(., ~ .[[3]]) %>%
  map(., ~ data.frame(ELPD = .[,2])) %>%
  map(., ~ rownames_to_column(., var = "predictor")) %>%
  set_names(processes_names) %>%
  bind_rows(.id = "process") %>%
  mutate(predictor = gsub("loo_fit_", "", predictor)) %>%
  mutate(ELPD = round(ELPD, 2))

elpd_differeences <- models %>% 
  map(., ~ .[[4]]) %>%
  map(., ~ rownames_to_column(., var = "predictor")) %>%
  set_names(processes_names) %>%
  bind_rows(.id = "process") %>%
  dplyr::rename("ELPD" = "ELPD difference") %>%
  mutate(Difference = paste(round(ELPD, 2), round(SE, 2), sep = "±")) %>%
  dplyr::select(-SE, -ELPD) %>%
  mutate(predictor = case_when(predictor == "General" ~ "general",
                               predictor == "Geomorph" ~ "general_geomorph",
                               predictor == "Forest" ~ "general_geomorph_forest",
                               predictor == "Disturbances" ~ "general_geomorph_forest_disturbances"))

model_comparsion <- left_join(elpd, elpd_differeences, by = c("process", "predictor")) %>%
  mutate(value = paste0(ELPD, " (", Difference, ")")) %>%
  dplyr::select( -ELPD, -Difference) %>%
  spread(key = predictor, value = value) %>%
  dplyr::select(process, null, general, general_geomorph, general_geomorph_forest, general_geomorph_forest_disturbances)

write_csv(model_comparsion, "../results/binomial/model_comparsion_binomial.csv")

### Extract final model

final_models <- models %>% map(., ~ .[[1]][[5]])

### Posterior predictive checks

pred_posterior_full <- final_models %>%
  map(~ posterior_predict(., draws = 100))

ppc_mean <- pred_posterior_full %>%
  map2(.y = list(data$MFL, data$FST),
       ~ bayesplot::ppc_stat(y = ifelse(.y > 0, 1, 0),
                             yrep = .x,
                             stat = "mean")) %>%
  map2(.y = processes_names, ~ . +
         theme_bw() +
         theme(panel.grid = element_blank(),
               strip.background = element_blank(),
               legend.position = "none") +
         scale_fill_manual(values = "grey") +
         scale_color_manual(values = "black") +
         labs(title = .y, x = "Mean probability", y = "Count")) %>%
  patchwork::wrap_plots(.)

ggsave("ppc_binomial.pdf", ppc_mean, path = "../results/binomial/", width = 5.5, height = 2.5)

# Extract and plot estimates ----------------------------------------------

estimates <- final_models %>%
  map(~ as.data.frame(.) %>%
        dplyr::select(-matches("Intercept")) %>% # Select everything that is not an intercept
        gather(key = varname, value = value) %>%
        left_join(vars_ws, by = "varname")) %>%
  set_names(processes_names) %>%
  bind_rows(.id = "process") %>%
  filter(!varname %in% c(paste0("eco_region", 1:9))) %>% 
  mutate(name = factor(name, levels = c("Area", "Artificial", "Elevation", "Elevation ratio", "Circularity", 
                                        "Melton ratio", "Elongation", "Forest", "Patch density", "Extent", 
                                        "Type", "Extent x Type"))) %>%
  mutate(type = case_when(name %in% c("Area", "Elevation", "Artificial") ~ "General",
                          name %in% c ("Elevation ratio", "Circularity", "Melton ratio", "Elongation") ~ "Geomorphology",
                          name %in% c("Forest", "Patch density", "Extent", "Type", "Extent x Type") ~ "Forest")) %>%
  mutate(type = factor(type, levels = c("General", "Geomorphology", "Forest")))

# Effect of forest cover on probability in % ("devide by four rule")

debris <- filter(estimates, process == "Debris-flow" & name == "Forest")

sediment <- filter(estimates, process == "Sediment-transport" & name == "Forest")

(median(debris$value)/4)*100

(median(sediment$value)/4)*100

# Plot effect sizes and directions

p_estimates <- ggplot(estimates, aes(x = fct_rev(name), y = value)) +
  geom_violin(aes(fill = paste0("  ", type))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank()) +
  coord_flip() +
  theme(strip.background = element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed", col = scales::muted("red")) +
  labs(y = "Effect size", x = NULL) +
  #scale_fill_brewer(palette = "Greys", direction = -1) +
  facet_wrap(~process) +
  scale_fill_manual(values = c("#276419","#ffffbf", "#4393c3"), breaks = c("  General", "  Geomorphology", "  Forest" )) +
  theme(legend.title = element_blank())

estimates$model <- "binomial"
write_csv(estimates, "../results/binomial/estimates_binomial.csv")

ggsave("estimates_binomial.pdf", p_estimates, path = "../results/binomial/", width = 5.5, height = 2.5, dpi = 300)

# Extract and plot eco region effects -----------------------------------------

# Estimates

ecoregion_effects <- final_models %>%
  map(~ as.data.frame(.) %>%
        dplyr::select(matches("Intercept"), eco_region2:eco_region9) %>% # Select everything that is an intercept
        mutate(draw = 1:4000) %>%
        gather(key = varname, value = value, -draw) %>%
        mutate(varname = ifelse(varname == "(Intercept)", "intercept", varname))) %>%
  set_names(processes_names) %>%
  bind_rows(.id = "process") %>%
  spread(key = varname, value = value) %>%
  gather(key = eco_region, value = value, -process, -draw, -intercept) %>%
  mutate(value = value + intercept)

p_eff <- ggplot(ecoregion_effects, aes(x = "", y = exp(value))) +
  geom_violin(aes(fill = fct_rev(eco_region))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank()) +
  coord_flip() +
  theme(strip.background = element_blank()) +
  geom_hline(data = ecoregion_effects %>% group_by(process) %>% summarize(m = mean(exp(intercept))),
             aes(yintercept = m), linetype = "dashed", col = scales::muted("red")) +
  labs(y = "Effect size", x = "", fill = "Ecological region") +
  scale_fill_brewer(palette = "Set1", labels = c("2","3","4","5","6","7","8","9"), breaks = c("eco_region2", "eco_region3", "eco_region4", "eco_region5", "eco_region6", "eco_region7", "eco_region8", "eco_region9")) +
  facet_wrap(~process) +
  theme(legend.position = "none") 

ecoregion_effects$model <- "binomial"
write_csv(ecoregion_effects, "../results/binomial/ecoregion_effects_binomial.csv")

ggsave("ecoregion_effects_binomial.pdf", p_eff, path = "../results/binomial", width = 7, height = 2.5)

# Map 

# load original shape file

raw_eco <- read_sf("../../../../../materials/raw_data/shapefiles_ecological_units/WLamPoly.shp")

p_map <- raw_eco %>%
  mutate(main = Wuchsge1) %>%
  mutate(main = substring(main, 1,1)) %>%
  st_set_precision(100) %>%
  group_by(main) %>%
  summarize() %>%
  ggplot(., aes( fill = factor(main))) +
  geom_sf() +
  theme_bw() +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  guides(fill = guide_legend(title = "Ecological region")) 

ggsave("../../../../../results/supplement/ecoregion_map_binomial.pdf", p_map, width = 7.5, height = 5)
ggsave("../../../../../results/supplement/ecoregion_map_binomial.png", p_map, width = 7.5, height = 5)

ggsave("../results/binomial/ecoregion_map_binomial.pdf", p_map, width = 7.5, height = 5)


# Create response curve plots ---------------------------------------------

# DFL

response_disturbance <- expand.grid(eco_region = factor(1),
                                    h_mean = 0,    
                                    Circularit = 0,
                                    Elongation = 0,
                                    artifical = 0,
                                    area = 0,
                                    patchdensity = 0,
                                    forest = 0,
                                    Elevation = 0,
                                    Melton = 0,
                                    extent = seq(quantile(scale(data$extent), 0), 
                                                 quantile(scale(data$extent), 0.99), 
                                                 length.out = 100),
                                    type = c(quantile(data_model$type, 0.05), 0, quantile(data_model$type, 0.95)))  # viele EZGs haben Pulisty kleiner als -1

predictions <- final_models[[1]] %>%
  posterior_linpred(., newdata = response_disturbance, transform = TRUE)

response_disturbance[, "mean"] <- apply(predictions, 2, mean)
response_disturbance[, "sd"] <- apply(predictions, 2, sd)  

p_response_dfl <- response_disturbance %>%
  mutate(type = factor(type, labels =  c("Press", "Average", "Pulse"))) %>%
  ggplot(., aes(x = extent, y = mean)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = type), alpha = 0.3) +
  geom_line(aes(col = type)) +
  geom_point(data = sample_n(data %>% mutate(extent = as.double(scale(extent))) %>% filter(extent < quantile(extent, 0.99)), 1000), 
             aes(x = extent, y = -0.01), shape = 124, alpha = 0.3) +
  theme_bw() +
  theme(legend.background = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size = 10),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.title = element_text(size = 9)) +
  scale_color_brewer(palette = "Set1", breaks = c("Press", "Average", "Pulse")) +
  scale_fill_brewer(palette = "Set1", breaks = c("Press", "Average", "Pulse")) +
  labs(x = "Disturbance extent", y = "Probability of occurrence", 
       fill = "Disturbance type", col = "Disturbance type",
       title = "Debris-flow") +
  guides(fill = guide_legend(ncol = 1, 
                             keywidth = 0.1,
                             keyheight = 0.1,
                             default.unit = "inch"))

ggsave("response_curve_binomial.pdf", p_response_dfl, path = "../results/binomial/", width = 3.5, height = 3.5)


