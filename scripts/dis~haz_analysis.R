### This r-script analysis the pre-processed data from dis~hat_data_processing_20180816.R
### written by J.Sebald 

# load packages, set working directory to either local (pc/macbook) or server (valenor)

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

  rm(list=ls())  
}

# Load data ---------------------------------------------------------------

data <- read.csv("../data/data_for_model.csv")

data$eco_unit <- as.factor(data$eco_unit) 

processes <- c("DFLOOD", "DFLOW", "FST")

shp_ws <- sf::read_sf("../data/shp_ws.shp")

# Watershed risk model with disturbances -----------------------------------------------

# create dataframe with variable short names and variable long names

vars_ws <- data.frame(varname = c("h_mean", "Melton", "Elevation", "Circularit", "Elongation", "artifical", "forest", "area", "patchdensity", 
                                  "severity", 
                                  "frequency",
                                  "severity:frequency"), 
                      name = c("Elevation", "Melton ratio", "Elevation ratio", "Circularity", "Elongtion", "Artificial", "Forest", "Area", "Patch density", 
                               "Severity", 
                               "Frequency",
                               "Severity x Frequency"),
                      stringsAsFactors = FALSE)

results <- vector("list", length = 3)

k <- 0

for (process in processes) {
  
  k <- k + 1
  
  # Bring data into form
  
  vars_nointeraction <- vars_ws %>% 
    filter(varname != "severity:frequency")
  
  data_model <- data %>%
    mutate_at(.vars = vars(c(vars_nointeraction$varname)), function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  data_model[, "response"] <- ifelse(data_model[, process] > 0, 1, 0)
  
  ### Calibrate models
  
  # Watershed-only model
  
  fit_ws_only <- stan_glm(as.formula(paste0("response ~ ", paste0(paste(vars_ws$varname[-which(vars_ws$varname %in% c( "severity", "frequency", "severity:frequency"))], collapse = "+")))),
                          data = data_model, 
                          family = binomial(link = "logit"), 
                          prior_intercept = normal(0|1), prior = normal(0|1))
  
  ws_varsel <- varsel(fit_ws_only)
  
  subset_size <- suggest_size(ws_varsel, alpha = 0.05)
  ws_predictors_selected <- names(ws_varsel$varsel$vind[1:subset_size])
  
  fit_ws_only_reduced <- stan_glmer(as.formula(paste0("response ~ ", paste0(paste(ws_predictors_selected, collapse = "+"), "+ (1|eco_unit)"))),
                                    data = data_model,
                                    family = binomial(link = "logit"),
                                    prior_intercept = normal(0|1), prior = normal(0|1))
  
  # Full model
  
  fit_full <- update(fit_ws_only_reduced, . ~ . + severity * frequency)
  
  # Watershed-only model
  
  fit_null <- update(fit_full, . ~ 1 + ( 1 | eco_unit))

  ### Compare models
  
  # Using LOO-ELPD
  
  loo_fit_ws_only <- loo(fit_ws_only)
  loo_fit_ws_only_reduced <- loo(fit_ws_only_reduced)
  loo_fit_full <- loo(fit_full)
  loo_fit_null <- loo(fit_null)
  
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
  
  auc_final <- get_auc(fit_full)
  
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
  
  estimates <- as.matrix(fit_full) %>%
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
                                      severity = seq(min(data_model$severity), quantile(data_model$severity, 0.99), length.out = 100),
                                      frequency = c(-1, 0, 1))

  # create predictons from the linear predictor with new data but on the basis of the fitted model
  
  predictions <- posterior_linpred(fit_full, newdata = response_disturbance, transform = TRUE, re.form = NA)
  
  # Calculate mean and sd of posterior predictions and add to new data
  
  response_disturbance$prob_mean <- apply(predictions, 2, mean)
  response_disturbance$prob_sd <- apply(predictions, 2, sd)
  
  ### Gather results and add to list
  
  results[[k]] <- list(fit_full, #1
                       model_comparison, #2
                       auc_final, #3
                       estimates, #4
                       randomeffects, #5
                       response_disturbance, #6
                       shape) #7
  
}

save(results, file = "../results/results.RData")

# Create plots and tables ---------------------------------------------------------

### Model comparison table

model_comparison <- results %>%
  map(~ .[[2]]) %>%
  set_names(processes) %>%
  bind_rows(.id = "process")

### Model estimates

p_estimate <- results %>% 
  map(~ .[[4]]) %>%
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

### Response curves

p_response <- results %>% 
  map(~ .[[6]]) %>%
  set_names(processes) %>%
  bind_rows(.id = "process") %>%
  mutate(process = factor(process, levels = c("FST", "DFLOOD", "DFLOW"))) %>%
  mutate(frequency = factor(frequency, labels = c("High (+1SD)", "Average", "Low (-1SD)"))) %>%
  ggplot(., aes(x = severity, y = prob_mean)) +
  geom_ribbon(aes(ymin = prob_mean - prob_sd, ymax = prob_mean + prob_sd, fill = frequency), alpha = 0.3) +
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

ggsave("responsecurves.pdf", p_response, path = "../results/", width = 7.5, height = 2.5)
ggsave("responsecurves.png", p_response, path = "../results/", width = 7.5, height = 2.5)

### Maps

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

