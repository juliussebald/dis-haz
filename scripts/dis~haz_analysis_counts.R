
# Load packages and set options -------------------------------------------

{
  Sys.setenv(USE_CXX14 = 1)
  library(raster)# version 2.6-7
  library(sp)# version 1.3-1
  library(rgdal)# version 1.2-20
  library(igraph)# version 1.2.1
  library(tidyverse)# version 1.2.1
  library(rstanarm)# version 2.1.7.4
  library(projpred)# version 0.8.0
  library(multiscales)#devtools::install_github("clauswilke/multiscales")
  library(patchwork)
  library(rstan)
  
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  
}

# Load data ---------------------------------------------------------------

data <- read.csv("../data/tables/data_for_model.csv") %>% 
  mutate(press = 1 - pulse) %>%
  mutate(DCOMBINE = DFLOOD + DFLOW)

data$eco_unit <- as.factor(data$eco_unit) 

processes <- c("DCOMBINE", "FST")
processes_names <- c("Mud-flow", "Sediment-flow")


# Selecting a model -----------------------------------------------

# Create dataframe with variable short names and variable long names

vars_ws <- data.frame(varname = c("h_mean", "Melton", "Elevation", "Circularit", "Elongation", 
                                  "artifical", "forest", "area", "patchdensity", 
                                  "eco_unit", 
                                  "extent", 
                                  "press",
                                  "extent:press"), 
                      name = c("Elevation", "Melton ratio", "Elevation ratio", "Circularity", "Elongation", 
                               "Artificial", "Forest", "Area", "Patch density", 
                               "Forest region", 
                               "Extent", 
                               "Press",
                               "Extent x Press"),
                      stringsAsFactors = FALSE)

# Loop through processes and calibrate varying models

models <- vector("list", length = length(processes))

k <- 0

for (process in processes) {
  
  k <- k + 1
  
  # Bring data into form
  
  vars_nointeraction <- vars_ws %>% 
    filter(!varname %in% c("extent:press", "eco_unit"))
  
  data_model <- data
  data_model[data_model$extent == 0, "press"] <- NA
  data_model <- data_model %>%
    mutate_at(.vars = vars(c(vars_nointeraction$varname)), function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  data_model[, "response"] <- data_model[, process]
  data_model[data_model$extent == min(data_model$extent), "press"] <- 0
  data_model$eco_unit <- as.factor(data_model$eco_unit)
  
  # Fit watershed-only model
  
  fit_ws_only <- stan_glm(as.formula(paste0("response ~ ", paste0(paste(vars_ws$varname[-which(vars_ws$varname %in% c("extent", "press", "extent:press"))], collapse = "+")))),
                          data = data_model,
                          family = neg_binomial_2,
                          prior = normal(0, 0.5, autoscale = TRUE),
                          prior_intercept = normal(0, 0.5, autoscale = TRUE),
                          prior_aux = exponential(rate = 1, autoscale = TRUE),
                          QR = TRUE)
  
  loo_fit_ws_only <- loo(fit_ws_only)
  
  # Include disturbance predictors
  
  fit_full_exp <- update(fit_ws_only, . ~ . + extent * press)
  
  # calculate loo for full model
  
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

save(models, file = "../results/two_processes/count/models_count.RData")
load(file =  "../results/two_processes/count/models_count.RData")

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

model_performances <- elpds_ws_only_model %>%
  left_join(elpds_full_model, by = "process") %>%
  left_join(elpds_difference, by = "process")

write_csv(model_performances, "../results/two_processes/count/model_performances_count.csv")

### Extract final model

final_models <- models %>% map2(.y = c(2, 2), ~ .[[.y]])

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

ggsave("ppc_count.pdf", ppc_mean, path = "../results/two_processes/count/", width = 5.5, height = 2.5)

# Extract and plot estimates ----------------------------------------------

estimates <- final_models %>%
  map(~ as.data.frame(.) %>%
        dplyr::select(-matches("Intercept")) %>% # Select everything that is not an intercept
        gather(key = varname, value = value) %>%
        left_join(vars_ws, by = "varname")) %>%
  set_names(processes_names) %>%
  bind_rows(.id = "process") %>%
  filter(!varname %in% c(paste0("eco_unit", 1:9), "reciprocal_dispersion")) %>% 
  mutate(name = factor(name, levels = c("Area", "Artificial", "Elevation", "Elevation ratio", "Circularity", 
                                         "Melton ratio", "Elongation", "Forest", "Patch density", "Extent", 
                                         "Press", "Extent x Press"))) %>%
  mutate(type = case_when(name %in% c("Area", "Elevation", "Artificial") ~ "General",
                          name %in% c ("Elevation ratio", "Circularity", "Melton ratio", "Elongation") ~ "Geomorphological",
                          name %in% c("Forest", "Patch density", "Extent", "Press", "Extent x Press") ~ "Forest related")) %>%
  mutate(type = factor(type, levels = c("General", "Geomorphological", "Forest related")))

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
  scale_fill_manual(values = c("#276419","#ffffbf", "#4393c3"), breaks = c("  General", "  Geomorphological", "  Forest related" )) +
  theme(legend.title = element_blank())

estimates$model <- "count"
write_csv(estimates, "../results/two_processes/count/estimates_count.csv")

ggsave("estimates_count.png", p_estimates, path = "../results/two_processes/count/", width = 5.5, height = 2.5)
ggsave("estimates_count.pdf", p_estimates, path = "../results/two_processes/count/", width = 5.5, height = 2.5)

# Extract and plot random effects -----------------------------------------

ecounit_effects <- final_models %>%
  map(~ as.data.frame(.) %>%
        dplyr::select(matches("Intercept"), eco_unit2:eco_unit9) %>% # Select everything that is an intercept
        mutate(draw = 1:4000) %>%
        gather(key = varname, value = value, -draw) %>%
        mutate(varname = ifelse(varname == "(Intercept)", "intercept", varname))) %>%
  set_names(processes_names) %>%
  bind_rows(.id = "process") %>%
  spread(key = varname, value = value) %>%
  gather(key = eco_unit, value = value, -process, -draw, -intercept) %>%
  mutate(value = value + intercept)

p_ecounit_effects <- ggplot(ecounit_effects, aes(x = fct_rev(eco_unit), y = exp(value))) +
  geom_violin(fill = "grey") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank()) +
  coord_flip() +
  theme(strip.background = element_blank()) +
  geom_hline(data = ecounit_effects %>% group_by(process) %>% summarize(m = mean(exp(intercept))),
             aes(yintercept = m), linetype = "dashed", col = scales::muted("red")) +
  labs(y = "Posterior probability distribution of parameter estimates", x = "Ecological unit", fill = "Process") +
  scale_fill_brewer(palette = "Greys", direction = -1) +
  facet_wrap(~process) +
  scale_x_discrete(labels = c("9","8","7","6","5","4","3","2","1"))

ecounit_effects$model <- "count"
write_csv(ecounit_effects, "../results/two_processes/count/ecounit_effects_count.csv")

ggsave("ecounit_effects_count.pdf", p_ecounit_effects, path = "../results/two_processes/count/", width = 5, height = 2.5)
ggsave("ecounit_effects_count.png", p_ecounit_effects, path = "../results/two_processes/count/", width = 5, height = 2.5)

# Expected counts plot ---------------------------------------------

# DCOMBINE

response_disturbance <- expand.grid(eco_unit = factor(1),
                                    h_mean = 0,    
                                    Circularit = 0,
                                    Elongation = 0,
                                    artifical = 0,
                                    area = 0,
                                    patchdensity = 0,
                                    forest = 0,
                                    Elevation = 0,
                                    Melton = 0,
                                    extent = c(quantile(data_model$extent, 0.95), 0, quantile(data_model$extent, 0.05)),
                                    press = c(quantile(data_model$press, 0.95), 0, quantile(data_model$press, 0.05)))  # viele EZGs haben Pulisty kleiner als -1

predictions <- final_models[[1]] %>%
  posterior_linpred(., newdata = response_disturbance, transform = TRUE)

response_disturbance[, "mean"] <- apply(predictions, 2, mean)
response_disturbance[, "sd"] <- apply(predictions, 2, sd)  

p_response_dcombine <- response_disturbance %>%
  mutate(press = factor(press, labels =  c("Pulse", "Average", "Press"))) %>%
  mutate(extent = factor(extent, labels =  c("Low extent", "Average extent", "Large extent"))) %>%
  split(list(.$press, .$extent)) %>%
  map(~ data.frame(count = 1:3, 
                   prop = dpois(x = 1:3, lambda = .$mean),
                   prop_lwr = dpois(x = 1:3, lambda = .$mean - .$sd),
                   prop_upr = dpois(x = 1:3, lambda = .$mean + .$sd))) %>%
  bind_rows(.id = "id") %>%
  separate("id", c("press", "extent"), "\\.") %>%
  mutate(press = factor(press, levels =  c("Press", "Average", "Pulse"))) %>%
  mutate(extent = factor(extent, levels =  c("Low extent", "Average extent", "High extent"))) %>%
  ggplot(., aes(x = as.factor(count), y = prop, fill = press)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = prop_lwr, ymax = prop_upr), position = position_dodge(width = 0.9), width = 0.25) +
  theme_bw() +
  theme(legend.background = element_blank(),
        legend.position = "none",
        legend.justification = c(1, 1),
        panel.grid = element_blank(),
        strip.background = element_blank()) +
  labs(x = "Number of events", y = "Probability", 
       fill = "Disturbance type", title = "a) Mud-flow") +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~extent) +
  ylim(0, 0.15)

# FST

response_disturbance <- expand.grid(eco_unit = factor(1),
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
                                    press = c(quantile(data_model$press, 0.95), 0, quantile(data_model$press, 0.05)))  # viele EZGs haben Pulisty kleiner als -1

predictions <- final_models[[2]] %>%
  posterior_linpred(., newdata = response_disturbance, transform = TRUE)

response_disturbance[, "mean"] <- apply(predictions, 2, mean)
response_disturbance[, "sd"] <- apply(predictions, 2, sd)  

p_response_fst <- response_disturbance %>%
  mutate(press = factor(press, labels =  c("Pulse", "Average", "Press"))) %>%
  split(.$press) %>%
  map(~ data.frame(count = 1:3, 
                   prop = dpois(x = 1:3, lambda = .$mean),
                   prop_lwr = dpois(x = 1:3, lambda = .$mean - .$sd),
                   prop_upr = dpois(x = 1:3, lambda = .$mean + .$sd))) %>%
  bind_rows(.id = "press") %>%
  mutate(press = factor(press, levels =  c("Press", "Average", "Pulse"))) %>%
  ggplot(., aes(x = as.factor(count), y = prop, fill = press)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = prop_lwr, ymax = prop_upr), position = position_dodge(width = 0.9), width = 0.25) +
  theme_bw() +
  theme(legend.background = element_blank(),
        legend.position = "right",
        legend.justification = c(1, 1),
        panel.grid = element_blank(),
        strip.background = element_blank()) +
  labs(x = "Number of events", y = "Probability", 
       fill = "Disturbance type", title = "b) Sediment-flow") +
  scale_fill_brewer(palette = "Set1") +
  guides(fill = guide_legend(ncol = 1, 
                             keywidth = 0.1,
                             keyheight = 0.1,
                             default.unit = "inch"))

# Combine plots

p_response <- p_response_dcombine + p_response_fst + plot_layout(ncol = 2, widths = c(3, 1.1))

ggsave("expected_counts.pdf", p_response, path = "../results/two_processes/count/", width = 7.5, height = 2.25)
ggsave("expected_counts.png", p_response, path = "../results/two_processes/count/", width = 7.5, height = 2.25)

