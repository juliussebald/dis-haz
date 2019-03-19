
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
 
  rm(list=ls()) 
}

# Load data ---------------------------------------------------------------

data <- read.csv("../data/data_for_model.csv") 

processes <- c("MFL", "FST")
processes_names <- c("Debris flow", "Flood")

# Selecting a model -----------------------------------------------

# Create dataframe with variable short names and variable long names

vars_ws <- data.frame(varname = c("h_mean", "Melton", "Elevation", "Circularit", "Elongation", 
                                  "artifical", "forest", "area", "patchdensity", 
                                  "eco_region", 
                                  "extent", 
                                  "type",
                                  "extent:type"), 
                      name = c("Elevation", "Melton ratio", "Elevation ratio", "Circularity", "Elongation", 
                               "Infrastructure", "Forest cover", "Area", "Patch density", 
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
  
  data_model[, "response"] <- data_model[, process]
  
  # make sure that watersheds which have not expierienced any disturbances are set to disturbance type "avarage"
  
  data_model[data_model$extent == min(data_model$extent), "type"] <- 0
  
  # transform eco region predictor to factor
  
  data_model$eco_region <- as.factor(data_model$eco_region)
  
  # Fit null model
  
  fit_null <- stan_glm(response ~ 1,
                       data = data_model,
                       family = neg_binomial_2,
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

save(models, file = "../results/count/models_count.RData")
load(file =  "../results/count/models_count.RData")

# Model evaluation -------------------------------------------------------------

### Model comparison

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

write_csv(model_comparsion, "../results/count/model_performances_count.csv")

### Extract final model

final_models <- models %>% map(., ~ .[[1]][[5]])

### Posterior predictive checks

pred_posterior_full <- final_models %>%
  map(~ posterior_predict(., draws = 100))

ppc_mean <- pred_posterior_full %>%
  map2(.y = list(data$FST), # TODO: Add data$DCOMBINE
       ~ bayesplot::ppc_bars(y = .y,
                             yrep = .x)) %>%
  map2(.y = processes_names, ~ . +
         theme_bw() +
         theme(panel.grid = element_blank(),
               strip.background = element_blank(),
               legend.position = "none") +
         labs(title = .y) +
         scale_size(limits = c(0, 0.5)) +
         scale_y_log10() +
         labs(x = "Number of events") +
         xlim(-1, 15)) %>%
  patchwork::wrap_plots(.)

ggsave("ppc_count.pdf", ppc_mean, path = "../results/count/", width = 5.5, height = 3.5)
ggsave("ppc_count.png", ppc_mean, path = "../../../../../results/supplement/", width = 5.5, height = 3.5)
ggsave("ppc_count.pdf", ppc_mean, path = "../../../../../results/supplement/", width = 5.5, height = 3.5)


# Extract and plot estimates ----------------------------------------------

estimates <- final_models %>%
  map(~ as.data.frame(.) %>%
        dplyr::select(-matches("Intercept")) %>% # Select everything that is not an intercept
        gather(key = varname, value = value) %>%
        left_join(vars_ws, by = "varname")) %>%
  set_names(processes_names) %>%
  bind_rows(.id = "process") %>%
  filter(!varname %in% c(paste0("eco_region", 1:9), "reciprocal_dispersion")) %>% 
  mutate(name = factor(name, levels = c("Area", "Infrastructure", "Elevation", "Elevation ratio", "Circularity", 
                                         "Melton ratio", "Elongation", "Forest cover", "Patch density", "Extent", 
                                         "Type", "Extent x Type"))) %>%
  mutate(type = case_when(name %in% c("Area", "Elevation", "Infrastructure") ~ "General",
                          name %in% c ("Elevation ratio", "Circularity", "Melton ratio", "Elongation") ~ "Geomorphology",
                          name %in% c("Forest cover", "Patch density", "Extent", "Type", "Extent x Type") ~ "Forest")) %>%
  mutate(type = factor(type, levels = c("General", "Geomorphology", "Forest")))

p_estimates <- ggplot(estimates, aes(x = fct_rev(name), y = value)) +
  geom_violin(aes(fill = paste0("  ", type))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank()) +
  coord_flip() +
  theme(strip.background = element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed", col = scales::muted("red")) +
  labs(y = "Posterior probability distribution of parameter estimates", x = NULL) +
  facet_wrap(~process) +
  scale_fill_manual(values = c("#276419","#ffffbf", "#4393c3"), breaks = c("  General", "  Geomorphology", "  Forest" )) +
  theme(legend.title = element_blank())

estimates$model <- "count"
write_csv(estimates, "../results/count/estimates_count.csv")

ggsave("estimates_count.pdf", p_estimates, path = "../results/count/", width = 5.5, height = 2.5)
ggsave("estimates_count.png", p_estimates, path = "../../../../../results/supplement/", width = 5.5, height = 2.5)

# Extract and plot eco_region effects -----------------------------------------

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

p_ecoregion_effects <- ggplot(ecoregion_effects, aes(x = fct_rev(eco_region), y = exp(value))) +
  geom_violin(fill = "grey") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank()) +
  coord_flip() +
  theme(strip.background = element_blank()) +
  geom_hline(data = ecoregion_effects %>% group_by(process) %>% summarize(m = mean(exp(intercept))),
             aes(yintercept = m), linetype = "dashed", col = scales::muted("red")) +
  labs(y = "Effect size", x = "Ecological region", fill = "Process") +
  scale_fill_brewer(palette = "Greys", direction = -1) +
  facet_wrap(~process) +
  scale_x_discrete(labels = c("9","8","7","6","5","4","3","2","1"))

ecoregion_effects$model <- "count"

write_csv(ecoregion_effects, "../results/count/ecoregion_effects_count.csv")

ggsave("ecoregion_effects_count.pdf", p_ecoregion_effects, path = "../results/count/", width = 5, height = 2.5)

# Expected counts plot ---------------------------------------------

# Re-create data_model data frame (same as for modeling above)

vars_nointeraction <- vars_ws %>% 
  filter(!varname %in% c("extent:type", "eco_region"))
data_model <- data
data_model[data_model$extent == 0, "type"] <- NA
data_model <- data_model %>%
  mutate_at(.vars = vars(c(vars_nointeraction$varname)), function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
data_model[data_model$extent == min(data_model$extent), "type"] <- 0

max_event <- 7

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
                                    extent = c((0.5 - mean(data$extent)) / sd(data$extent), 
                                               (0.1 - mean(data$extent)) / sd(data$extent)),

                                    type = c(quantile(data_model$type, 0.05), 0, quantile(data_model$type, 0.95)))  # viele EZGs haben Pulisty kleiner als -1

predictions_debris <- final_models[[1]] %>%
  posterior_predict(., newdata = response_disturbance, transform = TRUE)

plotdata_debris <- predictions_debris %>%
  as.data.frame(.) %>%
  mutate(draw = as.integer(1:4000)) %>%
  gather(., key = key, value = count, -draw) %>%
  group_by(key, count) %>%
  summarize(n = length(count)) %>%
  group_by(key) %>%
  mutate(p = (n / sum(n))) %>%
  ungroup() %>%
  dplyr::select(-n) %>%
  spread(key = count, value = p) %>%
  bind_cols(response_disturbance) %>%
  gather(key = count, value = p, -key, -(eco_region:type)) %>%
  filter(count %in% 1:max_event) %>%
  mutate(type = factor(type, labels =  c("Press", "Average", "Pulse"))) %>%
  mutate(extent = factor(extent, labels =  c("Small extent (10%)", "Large extent (50%)")))

p_debris <- ggplot(plotdata_debris, aes(x = factor(count, levels = 1:max_event), y = (p / length(1986:2018)), fill = factor(type))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~extent) +

  theme_bw() +
  theme(legend.background = element_blank(),
        legend.position = "right",
        legend.justification = c(1, 1),
        #legend.text = element_text(size = 9),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size = 11)) +

  labs(x = "Number of events", y = bquote("Probability (%"*yr^-1*")"), 
       fill = "Disturbance type", title = "a) Debris flow") +
  #scale_fill_manual(values = RColorBrewer::brewer.pal(6, name = "Reds")[c(2, 4, 6)]) +

  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(labels = function(x) round(x * 100, 3))

# Some numbers


plotdata_debris %>% 
  filter(count == 1) %>%
  mutate(prop_report = p / length(1986:2018) * 100) %>%
  dplyr::select(type, extent, prop_report)
  

# FST

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
                                    extent = 0,
                                    type = c(quantile(data_model$type, 0.05), 0, quantile(data_model$type, 0.95)))  # viele EZGs haben Pulisty kleiner als -1

predictions_flood <- final_models[[2]] %>%
  posterior_predict(., newdata = response_disturbance, transform = TRUE)

plotdata_flood <- predictions %>%
  as.data.frame(.) %>%
  mutate(draw = as.integer(1:4000)) %>%
  gather(., key = key, value = count, -draw) %>%
  group_by(key, count) %>%
  summarize(n = length(count)) %>%
  group_by(key) %>%
  mutate(p = (n / sum(n))) %>%
  ungroup() %>%
  dplyr::select(-n) %>%
  spread(key = count, value = p) %>%
  bind_cols(response_disturbance) %>%
  gather(key = count, value = p, -key, -(eco_region:type)) %>%
  filter(count %in% 1:max_event) %>%
  mutate(type = factor(type, labels =  c("Press", "Average", "Pulse")))

p_flood <- ggplot(plotdata, aes(x = factor(count, levels = 1:max_event), y = (p / length(1986:2018)), fill = factor(type))) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(legend.background = element_blank(),
        legend.position = "right",
        legend.justification = c(1, 1),
        legend.text = element_text(size = 9),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size = 11)) +
  labs(x = "Number of events", y = bquote("Probability (%"*yr^-1*")"), 
       fill = "Disturbance type", title = "b) Flood") +
  #scale_fill_manual(values = RColorBrewer::brewer.pal(6, name = "Reds")[c(2, 4, 6)]) +
  scale_fill_brewer(palette = "Set1") +

  scale_y_continuous(labels = function(x) round(x * 100, 3)) +
  guides(fill = guide_legend(keyheight = 0.4, keywidth = 0.6))


# Some numbers

plotdata_flood %>% 
  filter(count == 1) %>%
  mutate(prop_report = p / length(1986:2018) * 100) %>%
  dplyr::select(type, prop_report)


# Combine plots

p_response <- p_debris + theme(legend.position = "none") +
  p_flood + 
  plot_layout(ncol = 2, width = c(2, 1))


ggsave("expected_counts.pdf", p_response, path = "../results/count/", width = 7.5, height = 2.5)
ggsave("expected_counts.png", p_response, path = "../../../../../results/figures/", width = 7.5, height = 2.5)

