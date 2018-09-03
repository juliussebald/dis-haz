
# Load packages and set options -------------------------------------------

{
  library(raster)# version 2.6-7
  library(sp)# version 1.3-1
  library(rgdal)# version 1.2-20
  library(igraph)# version 1.2.1
  library(tidyverse)# version 1.2.1
  library(rstanarm)# version 2.1.7.4
  library(projpred)# version 0.8.0
  library(multiscales)#devtools::install_github("clauswilke/multiscales")
  
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = 4) 
  
<<<<<<< HEAD
=======
  rm(list=ls())  
>>>>>>> juliussebald/master
}

# Load data ---------------------------------------------------------------

data <- read.csv("../data/tables/data_for_model.csv")

data$eco_unit <- as.factor(data$eco_unit) 

processes <- c("DFLOOD", "DFLOW", "FST")

shp_ws <- sf::read_sf("../data/shapes/shp_ws.shp")

# Selecting a model -----------------------------------------------

# Create dataframe with variable short names and variable long names

vars_ws <- data.frame(varname = c("h_mean", "Melton", "Elevation", "Circularit", "Elongation", "artifical", "forest", "area", "patchdensity", 
                                  "extent", 
                                  "pulse",
<<<<<<< HEAD
                                  "extent:pulse"), 
                      name = c("Elevation", "Melton ratio", "Elevation ratio", "Circularity", "Elongtion", "Artificial", "Forest", "Area", "Patch density", 
                               "Extent", 
                               "Pulse",
                               "Extent x Pulse"),
                      stringsAsFactors = FALSE)

# Loop through processes and calibrate varying models

=======
                                  "rel_years",
                                  "extent:pulse",
                                  "extent:rel_years"), 
                      name = c("Elevation", "Melton ratio", "Elevation ratio", "Circularity", "Elongtion", "Artificial", "Forest", "Area", "Patch density", 
                               "Extent", 
                               "Pulse",
                               "rel_years",
                               "Extent x Pulse",
                               "Extent x Disturbance years"),
                      stringsAsFactors = FALSE)

>>>>>>> juliussebald/master
models <- vector("list", length = 3)

k <- 0

for (process in processes) {
  
  k <- k + 1
  
  # Bring data into form
  
  vars_nointeraction <- vars_ws %>% 
<<<<<<< HEAD
    filter(varname != c("extent:pulse"))
=======
    filter(varname != c("extent:pulse", "extent:rel_years"))
>>>>>>> juliussebald/master
  
  data_model <- data
  data_model[data_model$extent == 0, "pulse"] <- NA
  data_model <- data_model %>%
    mutate_at(.vars = vars(c(vars_nointeraction$varname)), function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  data_model[, "response"] <- ifelse(data_model[, process] > 0, 1, 0)
  data_model[data_model$extent == min(data_model$extent), "pulse"] <- 0
  
  # Fit watershed-only model
  
<<<<<<< HEAD
  fit_ws_only <- stan_glm(as.formula(paste0("response ~ ", paste0(paste(vars_ws$varname[-which(vars_ws$varname %in% c("extent", "pulse", "extent:pulse"))], collapse = "+")))),
=======
  # Watershed-only model
  
  fit_ws_only <- stan_glm(as.formula(paste0("response ~ ", paste0(paste(vars_ws$varname[-which(vars_ws$varname %in% c( "rel_years","extent", "pulse", "extent:pulse", "extent:rel_years"))], collapse = "+")))),
>>>>>>> juliussebald/master
                          data = data_model, 
                          family = binomial(link = "logit"), 
                          prior_intercept = normal(0|1), prior = normal(0|1))
  
  # Select important watershed predictors
  
  ws_varsel <- varsel(fit_ws_only)
  
  subset_size <- suggest_size(ws_varsel, alpha = 0.05)
  ws_predictors_selected <- names(ws_varsel$varsel$vind[1:subset_size])
  
<<<<<<< HEAD
  # Fit reduced watershed-only model
  
  fit_ws_only_reduced <- stan_glm(as.formula(paste0("response ~ ", paste0(paste(ws_predictors_selected, collapse = "+")))),
                                  data = data_model,
                                  family = binomial(link = "logit"),
                                  prior_intercept = normal(0|1), prior = normal(0|1))
  
  # Test if including Wuchsgebiete as random effect improves model
  
  fit_ws_only_reduced_rf <- stan_glmer(as.formula(paste0("response ~ ", paste0(paste(ws_predictors_selected, collapse = "+"), "+ (1|eco_unit)"))),
                                       data = data_model,
                                       family = binomial(link = "logit"),
                                       prior_intercept = normal(0|1), prior = normal(0|1))
=======
  fit_ws_only_reduced <- stan_glmer(as.formula(paste0("response ~ ", paste0(paste(ws_predictors_selected, collapse = "+"), "+ (1|eco_unit)"))),
                                    data = data_model,
                                    family = binomial(link = "logit"),
                                    prior_intercept = normal(0|1), prior = normal(0|1))
  
  # Full models
  
  fit_full_exp <- update(fit_ws_only_reduced, . ~ . + extent * pulse)
  fit_full_e_p <- update(fit_ws_only_reduced, . ~ . + extent + pulse)
  fit_full_exr <- update(fit_ws_only_reduced, . ~ . + extent * rel_years)
  fit_full_e_r <- update(fit_ws_only_reduced, . ~ . + extent + rel_years)
  
  
  
  # Null model
  
  fit_null <- update(fit_full_exp, . ~ 1 + ( 1 | eco_unit))
  
  ### Compare models
  
  # Using LOO-ELPD
>>>>>>> juliussebald/master
  
  loo_fit_ws_only_reduced <- loo(fit_ws_only_reduced)
<<<<<<< HEAD
  loo_fit_ws_only_reduced_rf <- loo(fit_ws_only_reduced_rf)
  
  comp <- loo::compare(loo_fit_ws_only_reduced, loo_fit_ws_only_reduced_rf)
  
  if ((comp[1] + comp[2]) > 0) {
    loo_fit_ws_only_reduced <- loo_fit_ws_only_reduced_rf
    fit_ws_only_reduced <- fit_ws_only_reduced_rf
  }
  
  # Include disturbance predictors and compare to watershed-only model
  
  fit_full_exp <- update(fit_ws_only_reduced, . ~ . + extent * pulse)
  fit_full_e_p <- update(fit_ws_only_reduced, . ~ . + extent + pulse)
  
  loo_fit_full_exp <- loo(fit_full_exp)
  loo_fit_full_e_p <- loo(fit_full_e_p)
  
  model_comparison <- loo::compare(loo_fit_ws_only_reduced,
                                   loo_fit_full_exp,
                                   loo_fit_full_e_p)
  
  # Store everything in a list
  
  models[[k]] <- list(fit_ws_only_reduced, #1
                      fit_full_exp, #2
                      fit_full_e_p, #3
                      list(loo_fit_ws_only_reduced,
                           loo_fit_full_exp,
                           loo_fit_full_e_p), #4
                      model_comparison) #5
                       
}

save(models, file = "../results/models.RData")
#load(file = "../results/models.RData")
=======
  loo_fit_full_exp <- loo(fit_full_exp)
  loo_fit_full_e_p <- loo(fit_full_e_p)
  loo_fit_full_exr <- loo(fit_full_exr)
  loo_fit_full_e_r <- loo(fit_full_e_r)
  loo_fit_null <- loo(fit_null)
  
  
elpds <-  rbind(null = loo_fit_null$estimates[1,], 
                ws_only = loo_fit_ws_only_reduced$estimates[1,],
                exp = loo_fit_full_exp$estimates[1,],
                e_p = loo_fit_full_e_p$estimates[1,],
                exr = loo_fit_full_exr$estimates[1,],
                e_r = loo_fit_full_e_r$estimates[1,])

  
  
  models[[k]] <- list(fit_full_exp,#1
                      fit_full_e_p,#2
                      fit_full_exr,#3
                      fit_full_e_r,#4
                      elpds)#5
                       
}

save(models, file = "../results/models.RData")



# Final model -------------------------------------------------------------

null_vs_ws_only <- loo::compare(loo_fit_null, loo_fit_ws_only_reduced)
ws_only_vs_full <- loo::compare(loo_fit_ws_only_reduced, loo_fit_full)

model_comparison <- as.data.frame(rbind(null_vs_ws_only, ws_only_vs_full))
model_comparison$comparison <- rownames(model_comparison)
rownames(model_comparison) <- NULL

# Calculate LOO AUC for final model

get_auc <- function(model) { # Function for calculating AUC of model via loo
  loo <- loo(model, save_psis = TRUE)
  preds <- posterior_linpred(model, transform = TRUE, re.form = NA)
  ploo <- loo::E_loo(preds, psis_object = loo$psis_object, type = "mean", log_ratios = -log_lik(model))$value
  auc <- AUC::auc(AUC::roc(ploo, factor(data_model$response)))
  return(auc)
}

auc_final <- get_auc(models[[1]][[2]])

# ### Plot roc curves -> I grayed this out as we do not really need those plots
#
# roc <- AUC::roc(ploo, factor(data_model$response))
# roc <- data.frame(cutoffs = roc$cutoffs, fpr = roc$fpr, tpr = roc$tpr)
# 
# p_roc <- ggplot(roc, aes(x = fpr, y = tpr)) +
#   geom_line() +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed", col = scales::muted("red")) +
#   theme_bw() +
#   labs(x = "False positive rate", y = "True positive rate", title = process)

### Extract estimates


# store estimates of model in format which is suitable to plot with ggplot
# the matrix contains 4000 draws from the posterior distribution of every variable which is included in
# the model

estimates <- as.matrix(models[[1]][[1]]) %>%
  as.data.frame() %>%
  dplyr::select(-matches("Intercept")) %>% # Select everything that is not an intercept
  gather(key = varname, value = value) %>%
  left_join(vars_ws, by = "varname") %>%
  mutate(process = process)

# Do the same stuff for the random effects + the scale parameter of the random effect

randomeffects <- as.matrix(fit_full) %>%
  as.data.frame() %>%
  dplyr::select(matches("Intercept")) %>% # Select everything that is an intercept
  gather(key = varname, value = value) %>%
  mutate(process = process)

### Create prediction for mapping

pred_posterior_full <- posterior_predict(fit_full)
pred_posterior_ws_only <- posterior_predict(fit_ws_only_reduced)

pred_posterior <- data_frame(WLK_ID = as.character(data_model$WLK_ID),
                             pred = apply(pred_posterior_full, 2, mean),
                             sd = apply(pred_posterior_full, 2, sd),
                             pred_diff = apply(pred_posterior_full - pred_posterior_ws_only, 2, mean),
                             sd_diff = apply(pred_posterior_full - pred_posterior_ws_only, 2, sd)) %>%
  mutate(., 
         varcof = sd / pred,
         varcof_diff = sd_diff / pred_diff)

shape <- shp_ws %>%
  left_join(pred_posterior, by = "WLK_ID")

#### Create plots of response curves

# bring all variables except disturbances to 0 to show the effect of disturbances on 
# natural hazard probability

response_disturbance <- expand.grid(h_mean = 0,    
                                    Circularit = 0,
                                    Elongation = 0,
                                    artifical = 0,
                                    area = 0,
                                    patchdensity = 0,
                                    forest = 0,
                                    Elevation = 0,
                                    Melton = 0,
                                    extent = seq(min(data_model$extent), quantile(data_model$extent, 0.99), length.out = 100),
                                    pulse = c(-1, 0, 1))

# create predictons from the linear predictor with new data but on the basis of the fitted model

predictions <- posterior_linpred(fit_full, newdata = response_disturbance, transform = TRUE, re.form = NA)

# Calculate mean and sd of posterior predictions and add to new data

response_disturbance$prob_mean <- apply(predictions, 2, mean)
response_disturbance$prob_sd <- apply(predictions, 2, sd)

>>>>>>> juliussebald/master


# Select final model -------------------------------------------------------------

### Compare models

model_comparison <- models %>% 
  map(~ as.data.frame(.[[5]]) %>%
        rownames_to_column(., var = "model")) %>%
  set_names(processes) %>%
  bind_rows(.id = "process")

write_csv(model_comparison, "../results/model_comparison.csv")

models %>%
  map(~ loo::compare(.[[4]][[3]], .[[4]][[2]]))

models %>%
  map(~ loo::compare(.[[4]][[1]], .[[4]][[2]]))

### Decide for final model

final_models <- models %>% map2(.y = c(2, 3, 3), ~ .[[.y]])

final_models <- models %>% map2(.y = c(2, 3, 3), ~ .[[2]])

### Calculate LOO AUC for final model

get_auc <- function(model) { # Function for calculating AUC of model via loo
  loo <- loo(model, save_psis = TRUE)
  preds <- posterior_linpred(model, transform = TRUE, re.form = NA)
  ploo <- loo::E_loo(preds, psis_object = loo$psis_object, type = "mean", log_ratios = -log_lik(model))$value
  auc <- AUC::auc(AUC::roc(ploo, factor(model$y)))
  return(auc)
}

auc <- final_models %>%
  map(~ get_auc(.))

auc %>%
  map(., ~ data.frame(auc = .)) %>%
  set_names(processes) %>%
  bind_rows(.id = "process") %>%
  write_csv(., "../results/model_auc.csv")

### Posterior predictive checks

pred_posterior_full <- final_models %>%
  map(~ posterior_predict(., draws = 100))

ppc_mean <- pred_posterior_full %>%
  map2(.y = list(data$DFLOOD, data$DFLOW, data$FST),
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

ggsave("ppc_mean.pdf", ppc_mean, path = "../results/", width = 7.5, height = 2.5)

# Create plots and tables ---------------------------------------------------------

### Extract and plot estimates

estimates <- final_models %>%
  map(~ as.data.frame(.) %>%
        dplyr::select(-matches("Intercept")) %>% # Select everything that is not an intercept
        gather(key = varname, value = value) %>%
        left_join(vars_ws, by = "varname")) %>%
  set_names(processes) %>%
  bind_rows(.id = "process")

p_estimates <- ggplot(estimates, aes(x = fct_rev(name), y = value)) +
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

ggsave("estimates.pdf", p_estimates, path = "../results/", width = 7.5, height = 2.5)


### Extract and plot random effects

randomeffects <- final_models %>%
  map(~ as.data.frame(.) %>%
        dplyr::select(matches("Intercept")) %>% # Select everything that is an intercept
        gather(key = varname, value = value)) %>%
  set_names(processes) %>%
<<<<<<< HEAD
  bind_rows(.id = "process")

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
  labs(y = "Posterior probability distribution of parameter estimates", x = NULL, fill = "Process") +
  scale_fill_brewer(palette = "Greys", direction = -1) +
  facet_wrap(~process)

ggsave("randomeffects.pdf", p_randomeffects, path = "../results/", width = 7.5, height = 2.5)


### Create response curve plots

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
                                    pulse = c(-1, 0, 1))

predictions <- final_models %>%
  map(~ posterior_linpred(., newdata = response_disturbance, transform = TRUE, re.form = NA))

for (i in 1:length(predictions)) {
  response_disturbance[, paste0(processes[i], "_mean")] <- apply(predictions[[i]], 2, mean)
  response_disturbance[, paste0(processes[i], "_sd")] <- apply(predictions[[i]], 2, sd)  
}

p_response <- response_disturbance %>%
  mutate(pulse = factor(pulse, labels = c("Low (-1SD)", "Average", "High (+1SD)"))) %>%
  gather(key = process, value = prop, -(h_mean:pulse)) %>%
  separate(process, c("process", "metric"), "_") %>%
  spread(key = metric, value = prop) %>%
  ggplot(., aes(x = extent, y = mean)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = pulse), alpha = 0.3) +
  geom_line(aes(col = pulse)) +
  geom_point(data = sample_n(data %>% mutate(extent = as.double(scale(extent))) %>% filter(extent < quantile(extent, 0.99)), 1000), 
             aes(x = extent, y = -0.01), shape = 124, alpha = 0.3) +
  #geom_point(data = data %>% 
  #             mutate(extent = as.double(scale(extent))) %>% 
  #             filter(extent < quantile(extent, 0.99)), 
  #           aes(x = extent, y = ifelse(DFLOOD > 0, 1, 0)),
  #           alpha = 0.3) +
  theme_bw() +
  theme(#legend.position = c(0, 1),
        #legend.justification = c(0, 1),
        legend.background = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank()) +
=======
  bind_rows(.id = "process") %>%
  mutate(process = factor(process, levels = c("FST", "DFLOOD", "DFLOW"))) %>%
  mutate(pulse = factor(pulse, labels = c("High (+1SD)", "Average", "Low (-1SD)"))) %>%
  ggplot(., aes(x = extent, y = prob_mean)) +
  geom_ribbon(aes(ymin = prob_mean - prob_sd, ymax = prob_mean + prob_sd, fill = pulse), alpha = 0.3) +
  geom_line(aes(col = pulse)) +
  geom_point(data = sample_n(data %>% mutate(extent = as.double(scale(extent))) %>% filter(extent < quantile(extent, 0.99)), 1000), 
             aes(x = extent, y = -0.05), shape = 124, alpha = 0.3) +
  theme_bw() +
  theme(#legend.position = c(0, 1),
    #legend.justification = c(0, 1),
    legend.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank()) +
  labs(x = "Disturbance extent", y = paste0("Probability of event"), 
       col = "Disturbance pulse", fill = "Disturbance pulse") +
>>>>>>> juliussebald/master
  scale_color_manual(values = c(scales::muted("blue"), "grey", scales::muted("red"))) +
  scale_fill_manual(values = c(scales::muted("blue"), "grey", scales::muted("red"))) +
  facet_wrap(~process) +
  labs(x = "Disturbance extent", y = "Probability of event", 
       fill = "Pulsity", col = "Pulsity")

ggsave("response_curves_test.pdf", p_response, path = "../results/", width = 7.5, height = 2.5)


# Posterior predictive checks with modulated predictor

response_disturbance <- expand.grid(h_mean = 0,    
                                    Circularit = 0,
                                    Elongation = 0,
                                    artifical = 0,
                                    area = 0,
                                    patchdensity = 0,
                                    forest = 0,
                                    Elevation = 0,
                                    Melton = 0,
                                    extent = seq(min(data$extent), quantile(data$extent, 0.99), length.out = 100),
                                    pulse = c(-1, 0, 1))

predictions <- final_models %>%
  map(~ posterior_linpred(., newdata = response_disturbance, transform = TRUE, re.form = NA))


# Create prediction maps

pred_posterior_full <- posterior_predict(fit_full)
pred_posterior_ws_only <- posterior_predict(fit_ws_only_reduced)

pred_posterior <- data_frame(WLK_ID = as.character(data_model$WLK_ID),
                             pred = apply(pred_posterior_full, 2, mean),
                             sd = apply(pred_posterior_full, 2, sd),
                             pred_diff = apply(pred_posterior_full - pred_posterior_ws_only, 2, mean),
                             sd_diff = apply(pred_posterior_full - pred_posterior_ws_only, 2, sd)) %>%
  mutate(., 
         varcof = sd / pred,
         varcof_diff = sd_diff / pred_diff)

shape <- shp_ws %>%
  left_join(pred_posterior, by = "WLK_ID")

for (i in 1:3) {
  
  shape <- results[[i]][[7]]
  
  pred_cutoff <- quantile(shape$pred, 0.99, na.rm = TRUE)
  pred_diff_cutoff <- quantile(shape$pred_diff, 0.99, na.rm = TRUE)
  varcof_cutoff <- quantile(shape$varcof, 0.99, na.rm = TRUE)
  varcof_diff_cutoff <- quantile(shape$varcof_diff, 0.99, na.rm = TRUE)
  
  normalize <- function(x, ...) (x - min(x, ...)) / (max(x, ...) - min(x, ...))
  
  map_prop <- shape %>%
    #sample_n(150) %>%
    mutate(pred = ifelse(pred >= pred_cutoff, pred_cutoff, pred),
           varcof = ifelse(varcof >= varcof_cutoff, varcof_cutoff, varcof)) %>%
    ggplot(.) + 
    geom_sf(aes(fill = pred, alpha = normalize((varcof^2 * -1), na.rm = TRUE)), color = NA, size = 0.2) + 
    coord_sf(datum = NA) +
    scale_fill_gradient(low = "#fff5f0", high = "#67000d") +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    labs(x = NULL, y = NULL, 
         fill = paste0("Probability of ", processes[i]), 
         #alpha = "Certainty") +
         alpha = NULL) +
    ggthemes::theme_map() +
    scale_alpha(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("Low", "", "", "", "High")) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.box = "horizontal",
          panel.grid = element_blank(),
          strip.background = element_blank()) +
    guides(fill = guide_colorbar(barwidth = 7.5, barheight = 0.25,
                                 direction = "horizontal", title.position = "top"),
           alpha = FALSE)
  
  map_diff <- shape %>%
    #sample_n(150) %>%
    ggplot(.) + 
    geom_sf(aes(fill = pred_diff), color = NA, size = 0.2) + 
    coord_sf(datum = NA) +
    scale_fill_gradient2(limits = c(-0.2, 0.2), low = scales::muted("blue"), high = scales::muted("red")) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    labs(x = NULL, y = NULL, 
         fill = paste0("Difference in probability when\nincluding disturbances")) +
    ggthemes::theme_map() +
    scale_alpha(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("Low", "", "", "", "High")) +
    theme(legend.position = "top",
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.box = "vertical",
          panel.grid = element_blank(),
          strip.background = element_blank()) +
    guides(fill = guide_colorbar(barwidth = 7.5, barheight = 0.25,
                                 direction = "horizontal", title.position = "top"))
  
  assign(paste0("map_prop_", i), map_prop)
  assign(paste0("map_diff_", i), map_diff)
  
}

p_map <- map_prop_3 + 
  map_prop_1 + 
  map_prop_2 + 
  map_diff_3 + 
  (map_diff_1 + theme(legend.position = "none")) + 
  (map_diff_2 + theme(legend.position = "none")) +
  plot_layout(ncol = 3)

ggsave("probability_maps.pdf", p_map, path = "../results", width = 7.5, height = 5.5)

<<<<<<< HEAD
=======
### Random effects

p_estimate <- results %>% 
  map(~ .[[5]]) %>%
  bind_rows() %>%
  mutate(process = factor(process, levels = c("FST", "DFLOOD", "DFLOW"))) %>%
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

ggsave("estimates.pdf", p_estimate, path = "../results/", width = 7.5, height = 2.5)
ggsave("estimates.png", p_estimate, path = "../results/", width = 7.5, height = 2.5)
>>>>>>> juliussebald/master
