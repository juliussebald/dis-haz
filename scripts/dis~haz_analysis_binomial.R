
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
  
  # Fit watershed-only model
  
  fit_ws_only <- stan_glm(as.formula(paste0("response ~ ", paste0(paste(vars_ws$varname[-which(vars_ws$varname %in% c("extent", "type", "extent:type"))], collapse = "+")))),
                          data = data_model,
                          family = binomial,
                          prior = normal(0, 0.5, autoscale = TRUE),
                          prior_intercept = normal(0, 0.5, autoscale = TRUE),
                          prior_aux = exponential(rate = 1, autoscale = TRUE),
                          QR = TRUE)
  
  # Calculate loo for watershed-pnly model
  
  loo_fit_ws_only <- loo(fit_ws_only)
  
  # Include disturbance predictors
  
  fit_full_exp <- update(fit_ws_only, . ~ . + extent * type)
  
  # Calculate loo for full model
  
  loo_fit_full_exp <- loo(fit_full_exp)
  
  # compare watershed model and full model
  
  model_comparison_ws_full <- loo::compare(loo_fit_ws_only,
                                           loo_fit_full_exp)
  
  # Store everything in a list
  
  models[[k]] <- list(fit_ws_only, #1
                      fit_full_exp, #2
                      list(loo_fit_ws_only,
                           loo_fit_full_exp), #3
                      model_comparison_ws_full #4
  )
                       
}

save(models, file = "../results/binomial/models_binomial.RData")

load(file =  "../results/binomial/models_binomial.RData")

# Model evaluation -------------------------------------------------------------

### Model comparison

elpds_ws_only_model <- models %>% 
  map(., ~ as.data.frame(.[[3]][[1]][[1]]) %>%
        rownames_to_column(., var = "model")) %>%
  set_names(processes_names) %>%
  bind_rows(.id = "process") %>%
  filter(model == "elpd_loo") %>%
  dplyr::select(process, elpd_ws_only_model = Estimate)

elpds_full_model <- models %>% 
  map(., ~ as.data.frame(.[[3]][[2]][[1]]) %>%
        rownames_to_column(., var = "model")) %>%
  set_names(processes_names) %>%
  bind_rows(.id = "process") %>%
  filter(model == "elpd_loo") %>%
  dplyr::select(process, elpd_full_model = Estimate)

elpds_difference <- models %>% 
  map(., ~ data.frame(elpd_dif_ws_full = .[[4]][[1]], 
                      elpd_dif_se_ws_full = .[[4]][[2]])) %>%
  set_names(processes_names) %>%
  bind_rows(.id = "process")

# get_auc <- function(model) { # Function for calculating AUC of model via loo
#   loo <- loo(model, save_psis = TRUE)
#   preds <- posterior_linpred(model, transform = TRUE, re.form = NA)
#   ploo <- loo::E_loo(preds, psis_object = loo$psis_object, type = "mean", log_ratios = -log_lik(model))$value
#   auc <- AUC::auc(AUC::roc(ploo, factor(model$y)))
#   return(auc)
# }

# auc_ws_only <- models %>%
#   map(., ~ .[[1]]) %>%
#   map(~ get_auc(.)) %>%
#   map(., ~ data.frame(auc_ws_only = .)) %>%
#   set_names(processes_names) %>%
#   bind_rows(.id = "process")
 
# auc_full <- models %>%
#   map(., ~ .[[2]]) %>%
#   map(~ get_auc(.)) %>%
#   map(., ~ data.frame(auc_full = .)) %>%
#   set_names(processes_names) %>%
#   bind_rows(.id = "process")

model_performances <- elpds_ws_only_model %>%
  left_join(elpds_full_model, by = "process") %>%
  left_join(elpds_difference, by = "process")

write_csv(model_performances, "../results/binomial/model_performances_binomial.csv")

write_csv(model_performances, "../../../../../results/tables/model_performances_binomial.csv")


### Extract final model

final_models <- models %>% map2(.y = c(2, 2), ~ .[[.y]])

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

ggsave("ppc_binomial.pdf", ppc_mean, path = "../../../../../results/supplement/", width = 5.5, height = 2.5)
ggsave("ppc_binomial.png", ppc_mean, path = "../../../../../results/supplement/", width = 5.5, height = 2.5)


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

ggsave("p_eff.pdf", p_eff, path = "../../../../../results/supplement/", width = 5.5, height = 2.5)
ggsave("p_eff.png", p_eff, path = "../../../../../results/supplement/", width = 5.5, height = 2.5)

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

ggsave("response_curve_binomial.png", p_response_dfl, path = "../../../../../results/figures", width = 3.5, height = 3.5)
ggsave("response_curve_binomial.pdf", p_response_dfl, path = "../../../../../results/figures", width = 3.5, height = 3.5)


