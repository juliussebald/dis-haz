
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


# Selecting a model -----------------------------------------------

# Create dataframe with variable short names and variable long names

vars_ws <- data.frame(varname = c("h_mean", "Melton", "Elevation", "Circularit", "Elongation", "artifical", "forest", "area", "patchdensity", 
                                  "extent", 
                                  "press",
                                  "extent:press"), 
                      name = c("Elevation", "Melton ratio", "Elevation ratio", "Circularity", "Elongation", "Artificial", "Forest", "Area", "Patch density", 
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
    filter(varname != c("extent:press"))
  
  data_model <- data
  data_model[data_model$extent == 0, "press"] <- NA
  data_model <- data_model %>%
    mutate_at(.vars = vars(c(vars_nointeraction$varname)), function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  data_model[, "response"] <- ifelse(data_model[, process] > 0, 1, 0)
  data_model[data_model$extent == min(data_model$extent), "press"] <- 0
  
  # # Null-model
  # 
  # fit_null_rf <- stan_glmer(response ~ 1 + (1|eco_unit),
  #                           data = data_model,
  #                           family = binomial(link = "logit"),
  #                           prior_intercept = normal(0|1), prior = normal(0|1))
  # 
  # loo_fit_null_rf <- loo(fit_null_rf)
  
  # Fit watershed-only model
  
  fit_ws_only <- stan_glm(as.formula(paste0("response ~ ", paste0(paste(vars_ws$varname[-which(vars_ws$varname %in% c("extent", "press", "extent:press"))], collapse = "+")))),
                          data = data_model, 
                          family = binomial(link = "logit"), 
                          prior_intercept = normal(0|1), prior = normal(0|1))
  
  # Select important watershed predictors
  
  ws_varsel <- varsel(fit_ws_only)
  
  subset_size <- suggest_size(ws_varsel, alpha = 0.05)
  ws_predictors_selected <- names(ws_varsel$varsel$vind[1:subset_size])
  
  # Include ecological units as random effects 
  
  fit_ws_only_reduced_rf <- stan_glmer(as.formula(paste0("response ~ ", paste0(paste(ws_predictors_selected, collapse = "+"), "+ (1|eco_unit)"))),
                                       data = data_model,
                                       family = binomial(link = "logit"),
                                       prior_intercept = normal(0|1), prior = normal(0|1))
  
  # calculate loo of reduced watershed model with ecological units as random effects 
  
  loo_fit_ws_only_reduced_rf <- loo(fit_ws_only_reduced_rf)
  
  # Include disturbance predictors
  
  fit_full_exp <- update(fit_ws_only_reduced_rf, . ~ . + extent * press)
  
  # calculate loo for full model
  
  loo_fit_full_exp <- loo(fit_full_exp)
  
  # compare watershed model and full model
  
  # model_comparison_null_ws <- loo::compare(loo_fit_null_rf,
  #                                          loo_fit_ws_only_reduced_rf)
  
  model_comparison_ws_full <- loo::compare(loo_fit_ws_only_reduced_rf,
                                           loo_fit_full_exp)

  # Store everything in a list
  
  models[[k]] <- list(fit_ws_only_reduced_rf, #1
                      fit_full_exp, #2
                      list(loo_fit_ws_only_reduced_rf,
                           loo_fit_full_exp), #3
                      model_comparison_ws_full #4
                      )
                       
}




save(models, file = "../results/two_processes/random_effects_models_two_processes.RData")
load(file =  "../results/two_processes/random_effects_models_two_processes.RData")


# Select final model -------------------------------------------------------------

### Compare models


# elpds_null_model <- models %>% 
#   map(., ~ as.data.frame(.[[4]][[1]][[1]]) %>%
#         rownames_to_column(., var = "model")) %>%
#   set_names(processes) %>%
#   bind_rows(.id = "process") %>%
#   filter(model == "elpd_loo") %>%
#   dplyr::select(process, elpd_null_model = Estimate)


elpds_ws_only_model <- models %>% 
  map(., ~ as.data.frame(.[[3]][[1]][[1]]) %>%
        rownames_to_column(., var = "model")) %>%
  set_names(processes) %>%
  bind_rows(.id = "process") %>%
  filter(model == "elpd_loo") %>%
  dplyr::select(process, elpd_ws_only_model = Estimate)


elpds_full_model <- models %>% 
  map(., ~ as.data.frame(.[[3]][[2]][[1]]) %>%
        rownames_to_column(., var = "model")) %>%
  set_names(processes) %>%
  bind_rows(.id = "process") %>%
  filter(model == "elpd_loo") %>%
  dplyr::select(process, elpd_full_model = Estimate)


elpds_difference <- models %>% 
  map(., ~ data.frame(elpd_dif_ws_full = .[[4]][[1]], 
                      elpd_dif_se_ws_full = .[[4]][[2]])) %>%
  set_names(processes) %>%
  bind_rows(.id = "process")


get_auc <- function(model) { # Function for calculating AUC of model via loo
  loo <- loo(model, save_psis = TRUE)
  preds <- posterior_linpred(model, transform = TRUE, re.form = NA)
  ploo <- loo::E_loo(preds, psis_object = loo$psis_object, type = "mean", log_ratios = -log_lik(model))$value
  auc <- AUC::auc(AUC::roc(ploo, factor(model$y)))
  return(auc)
}



auc_ws_only <- models %>%
  map(., ~ .[[1]]) %>%
  map(~ get_auc(.)) %>%
  map(., ~ data.frame(auc_ws_only = .)) %>%
  set_names(processes) %>%
  bind_rows(.id = "process")

auc_full <- models %>%
  map(., ~ .[[2]]) %>%
  map(~ get_auc(.)) %>%
  map(., ~ data.frame(auc_full = .)) %>%
  set_names(processes) %>%
  bind_rows(.id = "process")

model_performances <- elpds_ws_only_model %>%
  left_join(elpds_full_model, by = "process") %>%
  left_join(elpds_difference, by = "process") %>%
  left_join(auc_ws_only, by = "process") %>%
  left_join(auc_full, by = "process")

write_csv(model_performances, "../results/two_processes/model_performances.csv")

### Decide for final model

final_models <- models %>% map2(.y = c(2, 2), ~ .[[.y]])

### Posterior predictive checks

pred_posterior_full <- final_models %>%
  map(~ posterior_predict(., draws = 100))

ppc_mean <- pred_posterior_full %>%
  map2(.y = list(data$DCOMBINE, data$FST),
       ~ bayesplot::ppc_stat(y = ifelse(.y > 0, 1, 0),
                             yrep = .x,
                             stat = "mean")) %>%
  map2(.y = processes, ~ . +
         theme_bw() +
         theme(legend.position = "none") +
         scale_fill_manual(values = "grey") +
         scale_color_manual(values = "black") +
         labs(title = .y)) %>%
  patchwork::wrap_plots(.)

ggsave("ppc_mean.pdf", ppc_mean, path = "../results/two_processes/", width = 5.5, height = 2.5)

# Create plots and tables ---------------------------------------------------------



# Extract and plot estimates ----------------------------------------------

estimates <- final_models %>%
  map(~ as.data.frame(.) %>%
        dplyr::select(-matches("Intercept")) %>% # Select everything that is not an intercept
        gather(key = varname, value = value) %>%
        left_join(vars_ws, by = "varname")) %>%
  set_names(processes) %>%
  bind_rows(.id = "process") %>%
  mutate(process = factor(process, levels = c("DCOMBINE","FST"), labels = c("mudflow", "sediment transport"))) %>%
  mutate(name = factor(name, levels = c("Area", "Elevation", "Elevation ratio", "Melton ratio", "Elongation", "Forest", "Patch density", "Extent", "Press", "Extent x Press"))) %>%
  mutate(type = case_when(name %in% c("Area", "Elevation") ~ "General",
                          name %in% c ("Elevation ratio", "Melton ratio", "Elongation") ~ "Geomorphological",
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
  #scale_fill_brewer(palette = "Greys", direction = -1) +
  facet_wrap(~process) +
  scale_fill_manual(values = c("#276419","#ffffbf", "#4393c3"), breaks = c("  General", "  Geomorphological", "  Forest related" )) +
  theme(legend.title = element_blank())


ggsave("estimates.png", p_estimates, path = "../results/two_processes/", width = 5.5, height = 2.5)
ggsave("estimates.pdf", p_estimates, path = "../results/two_processes/", width = 5.5, height = 2.5, dpi = 300)

# Extract and plot random effects -----------------------------------------

randomeffects <- final_models %>%
  map(~ as.data.frame(.) %>%
        dplyr::select(matches("Intercept")) %>% # Select everything that is an intercept
        gather(key = varname, value = value)) %>%
  set_names(processes) %>%
  bind_rows(.id = "process") %>%
  mutate(process = factor(process, levels = c("DCOMBINE","FST")))

randomeffects <- randomeffects %>%
  filter(!varname %in% c("(Intercept)", "Sigma[eco_unit:(Intercept),(Intercept)]")) %>%
  mutate(name = substr(varname, 15, 24)) 

p_randomeffects <- ggplot(randomeffects, aes(x = fct_rev(name), y = value)) +
  geom_violin(fill = "grey") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank()) +
  coord_flip() +
  theme(strip.background = element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed", col = scales::muted("red")) +
  labs(y = "Posterior probability distribution of parameter estimates", x = "Ecological unit", fill = "Process") +
  scale_fill_brewer(palette = "Greys", direction = -1) +
  facet_wrap(~process) +
  scale_x_discrete(labels = c("9","8","7","6","5","4","3","2","1"))

ggsave("randomeffects.pdf", p_randomeffects, path = "../results/two_processes/", width = 5.5, height = 2.5)
ggsave("randomeffects.png", p_randomeffects, path = "../results/two_processes/", width = 5.5, height = 2.5)

# Create response curve plots ---------------------------------------------

response_disturbance <- expand.grid(h_mean = 0,    
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
                                    press = c(quantile(data_model$press, 0.95), 0, quantile(data_model$press, 0.05)))  # viele EZGs haben Pulisty kleiner als -1

predictions <- final_models %>%
  map(~ posterior_linpred(., newdata = response_disturbance, transform = TRUE, re.form = NA))

for (i in 1:length(predictions)) {
  response_disturbance[, paste0(processes[i], "_mean")] <- apply(predictions[[i]], 2, mean)
  response_disturbance[, paste0(processes[i], "_sd")] <- apply(predictions[[i]], 2, sd)  
}


p_response <- response_disturbance %>%
  mutate(press = factor(press, labels =  c("Pulse", "Average", "Press"))) %>% # Warum muss ich das hier verkehrt herum reinschreiben????
  gather(key = process, value = prop, -(h_mean:press)) %>%
  separate(process, c("process", "metric"), "_") %>%
  spread(key = metric, value = prop) %>%
  mutate(process = factor(process, levels = c("DCOMBINE","FST"), labels = c("mudflow", "sediment transport"))) %>%
  ggplot(., aes(x = extent, y = mean)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = press), alpha = 0.3) +
  geom_line(aes(col = press)) +
  geom_point(data = sample_n(data %>% mutate(extent = as.double(scale(extent))) %>% filter(extent < quantile(extent, 0.99)), 1000), 
             aes(x = extent, y = -0.01), shape = 124, alpha = 0.3) +
  theme_bw() +
  theme(#legend.position = c(0, 1),
    #legend.justification = c(0, 1),
    legend.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank()) +
  scale_color_manual(values = c(scales::muted("blue"), "grey", scales::muted("red")), breaks = c("Press", "Average", "Pulse")) +
  scale_fill_manual(values = c(scales::muted("blue"), "grey", scales::muted("red")), breaks = c("Press", "Average", "Pulse")) +
  facet_wrap(~process) +
  labs(x = "Disturbance extent", y = "Probability of event", 
       fill = "Disturbance type", col = "Disturbance type")

ggsave("response_curves.png", p_response, path = "../results/two_processes/", width = 5.5, height = 2.5)
ggsave("response_curves.pdf", p_response, path = "../results/two_processes/", width = 5.5, height = 2.5, dpi = 300)

# Difference in risk w/o disturbances -------------------------------------

# predict_posterior_ws_only <- models %>%
#   map(., ~ .[[1]]) %>%
#   map(., ~ posterior_linpred(., draws = 200, transform = TRUE))
# 
# predict_posterior_full_model <- models %>%
#   map(., ~ .[[2]]) %>%
#   map(., ~ posterior_linpred(., draws = 200, transform = TRUE))
# 
# predict_posterior_difference <- predict_posterior_ws_only %>%
#   map2(.y = predict_posterior_full_model, ~ .y - .x)
# 
# predict_posterior_difference <- predict_posterior_difference %>%
#   map(., ~ apply(., 2, mean)) %>%
#   map2(.y = models, ~ data.frame(WLK = .y[[1]]$data$WLK_ID,
#                       difference = .)) %>%
#   set_names(processes) %>%
#   bind_rows(.id = "process")
# 
# predict_posterior_difference <- predict_posterior_difference %>%
#   left_join(data, by = c("WLK" = "WLK_ID")) %>%
#   mutate(press_cut = cut(press, quantile(press, c(0, 0.5, 1)), include.lowest = TRUE, labels = c("lp", "hp")),
#          extent_cut = cut(extent, quantile(extent, c(0, 0.5, 1)), include.lowest = TRUE, labels = c("le", "he"))) %>%
#   mutate(disturbance_type = paste(press_cut, extent_cut, sep = "-"))
# 
# ggplot(predict_posterior_difference, aes(x = as.integer(reorder(WLK, difference)), y = difference, col = disturbance_type)) +
#   geom_point(alpha = 0.3) +
#   geom_vline(xintercept = 0, col = "red") +
#   facet_wrap(~process, scales = "free_x") +
#   coord_flip()
# 
# predict_posterior_difference %>%
#   group_by(disturbance_type, process) %>%
#   mutate(diff_mean = mean(difference),
#          diff_sd = sd(difference)) %>%
#   ungroup() %>%
#   ggplot(., aes(x = process, y = diff_mean, col = disturbance_type)) +
#   geom_point(position = position_dodge(width = 0.25)) +
#   geom_errorbar(aes(ymin = diff_mean - diff_sd, ymax = diff_mean + diff_sd), width = 0,
#                 position = position_dodge(width = 0.25)) +
#   coord_flip() +
#   geom_hline(yintercept = 0)
# 
# predict_posterior_ws_only_counts <- models %>%
#   map(., ~ .[[1]]) %>%
#   map(., ~ posterior_predict(., draws = 500)) %>%
#   map(., ~ apply(., 1, sum))
# 
# predict_posterior_full_model_counts <- models %>%
#   map(., ~ .[[2]]) %>%
#   map(., ~ posterior_predict(., draws = 500)) %>%
#   map(., ~ apply(., 1, sum))
# 
# predict_posterior_difference <- predict_posterior_ws_only_counts %>%
#   map2(.y = predict_posterior_full_model_counts, 
#        ~ .y - .x) %>%
#   map2(.y = models, ~ data.frame(draws = 1:500,
#                                  difference = .)) %>%
#   set_names(processes) %>%
#   bind_rows(.id = "process")
# 
# predict_posterior_difference %>%
#   group_by(process) %>%
#   summarize(diff = median(difference))
# 
# ggplot(predict_posterior_difference, aes(x = difference)) +
#   geom_histogram() +
#   facet_wrap(~process, scales = "free") +
#   geom_vline(xintercept = 0, col = "red")

