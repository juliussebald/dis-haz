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

  
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 4) 

rm(list=ls())  
}

# Load data ---------------------------------------------------------------

data <- read.csv("../data/data_for_model.csv")

data$eco_unit <- as.factor(data$eco_unit) 

processes <- c("DFLOOD", "DFLOW", "FST")

shp_ws <- raster::shapefile("../data/shp_ws.shp")
shp_ws_fortify <- broom::tidy(shp_ws, region = "WLK_ID")


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
  
  # Calibrate model
  
  fit_ws <- stan_glmer(as.formula(paste0("response ~ ", paste0(paste(vars_ws$varname, collapse = "+"), "+ (1|eco_unit)"))),
                     data = data_model, family = binomial(link = "logit"), prior_intercept = normal(0|1), prior = normal(0|1))
  

  # Reduce predictors by projection predictive variable selection
  
  fit_ws_varsel <- varsel(fit_ws, method = "forward")
  
  subset_size <- suggest_size(fit_ws_varsel, alpha = 0.05)
  
  p_varsel <- varsel_plot(fit_ws_varsel) +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          legend.position = "none") +
    labs(y = "Expected log predictive density", title = process)
  
  # Make sure that severity and frequency are included also as single variables when interaction is included
  
  predictors_selected <- names(fit_ws_varsel$varsel$vind[1:subset_size])
  if ("severity:frequency" %in% predictors_selected) {
    if (!"severity" %in% predictors_selected) {
      predictors_selected <- c(predictors_selected, "severity")
      subset_size <- subset_size + 1
    }
    if (!"frequency" %in% predictors_selected) {
      predictors_selected <- c(predictors_selected, "frequency")
      subset_size <- subset_size + 1
    }
  }
  
  # create table of seleceted predictors and join them with their long names which occur in the plots later
  # and add their relativ importance for the predictive power of the model
  
  predictors_selected_names <- data.frame(varname = predictors_selected, 
                                          importance = 1:length(predictors_selected)) %>%
    left_join(vars_ws, by = "varname")
  
  # calibrate model with reduced number of variables
  
  fit_ws_reduced <- stan_glm(as.formula(paste0("response ~ ", paste(predictors_selected, collapse = "+"))),
                             data = data_model, family = binomial(link = "logit"), 
                             prior_intercept = normal(0|1), prior = normal(0|1))
  

  # store estimates of model in format which is suitable to plot with ggplot
  # the matrix contains 4000 draws from the posterior distribution of every variable which is included in
  # the model
  
  estimates <- as.matrix(fit_ws_reduced) %>%
    as.data.frame() %>%
    dplyr::select(-1) %>%
    dplyr::select(.dots = vars$varnames) %>%
    gather(key = varname, value = value) %>%
    left_join(vars_ws, by = "varname") %>%
    mutate(process = process)
  
  randomeffects <- as.matrix(fit_ws_reduced) %>%
    as.data.frame() %>%
    dplyr::select(-1) %>%
    dplyr::select(.dots = -(vars$varnames)) %>%
    gather(key = varname, value = value) %>%
    left_join(vars_ws, by = "varname") %>%
    mutate(process = process)
  
  
  # Assessing the model fit: loo = leave one out cross-validation, elpd = expected log predictive density
  
  loo <- loo(fit_ws_reduced, save_psis = TRUE)
  elpd <- loo$estimates[1,1]
  
  # compute matrix of predictions 
  
  preds <- posterior_linpred(fit_ws_reduced, transform = TRUE, re.form = NA)
  
  # expected mean of probablity distribution for all watersheds
  
  ploo <- loo::E_loo(preds, psis_object = loo$psis_object, type = "mean", log_ratios = -log_lik(fit_ws_reduced))$value
  
  # calculate auc for model performance, auc of 0,5 would mean the model prediction is totaly random
  # a value of 1 would mean a perfect model performance
  
  auc <- AUC::auc(AUC::roc(ploo, factor(data_model$response)))
  
  # compute data_frame to plot roc curves; the dataframe contains cutoffs, fpr = false positive rate
  # and tpr = true positive rate
  
  roc <- AUC::roc(ploo, factor(data_model$response))
  roc <- data.frame(cutoffs = roc$cutoffs, fpr = roc$fpr, tpr = roc$tpr)
  
  # plot roc curve
  
  p_roc <- ggplot(roc, aes(x = fpr, y = tpr)) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", col = scales::muted("red")) +
    theme_bw() +
    labs(x = "False positive rate", y = "True positive rate", title = process)
  
  # Create prediction for mapping
  
  # fit model without disturbances to get general hazard exposure of watershed, based on 
  # geomorphological and topographic variables
  
  fit_ws_reduced_nodist <- update(fit_ws_reduced, . ~ . - severity - frequency - severity:frequency)
  
  # compute matrix of predictions when disturbances are not included in the model
  
  pred_nodist <- posterior_linpred(fit_ws_reduced_nodist, transform = TRUE, re.form = NA)
  
  # create dataframe of predictions with disturbances included and without disturbances
  
  predictions <- data.frame(WLK_ID = data_model$WLK_ID,
                            prob_hazard_mean = apply(preds, 2, mean),
                            prob_hazard_sd = apply(preds, 2, sd),
                            prob_hazard_nodist_mean = apply(pred_nodist, 2, mean),
                            prob_hazard_nodist_sd= apply(pred_nodist, 2, sd))
  
  # add predictions to shapefile of watersheds

  shp_ws_fortify_predictions <- shp_ws_fortify %>%
    mutate(WLK_ID = as.integer(id)) %>%
    left_join(predictions, by = "WLK_ID")
  
  # create map of natural hazard event probability 
  
  p_probabilitymap_all <- shp_ws_fortify_predictions %>%
    ggplot(.) +
    aes(long, lat, group = group, fill = prob_hazard_mean) +
    geom_polygon() +
    coord_equal() +
    scale_fill_gradient(low = "#fff5f0", high = "#67000d", limit = c(0, quantile(shp_ws_fortify_predictions$prob_hazard_mean, 0.99, na.rm = TRUE))) +
    labs(x = NULL, y = NULL, fill = paste0("P(", process, ")")) +
    ggthemes::theme_map() +
    theme(legend.justification = c(0, 1),
          legend.position = c(0, 1),
          legend.background = element_blank()) +
    guides(fill = guide_colorbar(barwidth = 6, barheight = 0.5,
                                 direction = "horizontal", title.position = "top"))
  
  # create map of difference in natural hazard probability when disturbances are included in the model
  
  p_probabilitymap_difference <- shp_ws_fortify_predictions %>%
    ggplot(.) +
    aes(long, lat, group = group, fill = prob_hazard_mean - prob_hazard_nodist_mean) +
    geom_polygon() +
    coord_equal() +
    scale_fill_gradient2() +
    labs(x = NULL, y = NULL, fill = paste0("P(", process, ")")) +
    ggthemes::theme_map() +
    theme(legend.justification = c(0, 1),
          legend.position = c(0, 1),
          legend.background = element_blank()) +
    guides(fill = guide_colorbar(barwidth = 6, barheight = 0.5,
                                 direction = "horizontal", title.position = "top"))
    
  # Create plots of response curves classical
  
  # bring all variables except disturbances to 0 to show the effect of disturbances on 
  # natural hazard probability
  
  newdata <- expand.grid(h_mean = 0,    
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

  # create predictons with new data but on the basis of the fitted model
  
  predictions <- posterior_linpred(fit_ws, newdata = newdata, transform = TRUE, re.form = NA)

  newdata$prob_mean <- apply(predictions, 2, mean)
  newdata$prob_sd <- apply(predictions, 2, sd)
  
  # plot response curve

  p_responsecurve_classical <- newdata %>%
    mutate(frequency = factor(frequency, labels = c("High (+1SD)", "Average", "Low (-1SD)"))) %>%
    ggplot(., aes(x = severity, y = prob_mean)) +
    geom_ribbon(aes(ymin = prob_mean - prob_sd * 2, ymax = prob_mean + prob_sd * 2, fill = frequency), alpha = 0.3) +
    geom_line(aes(col = frequency)) +
    geom_point(data = sample_n(data_model %>% filter(severity < quantile(severity, 0.99)), 1000), aes(x = severity, y = -0.01), shape = 124, alpha = 0.3) +
    theme_bw() +
    theme(legend.position = c(0, 1),
          legend.justification = c(0, 1),
          legend.background = element_blank(),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8)) +
    labs(x = "Disturbance severity", y = paste0("P(", process, ")"), 
         col = "Disturbance frequency", fill = "Disturbance frequency") +
    scale_color_brewer(palette = "Set1")
  
  # create plots of response curve as heatmap
  
  newdata <- expand.grid(h_mean = 0,    
                         Circularit = 0,
                         Elongation = 0,
                         artifical = 0,
                         area = 0,
                         patchdensity = 0,
                         forest = 0,
                         Elevation = 0,
                         Melton = 0,
                         severity = seq(min(data_model$severity), quantile(data_model$severity, 0.99), length.out = 100),
                         frequency = seq(min(data_model$frequency), quantile(data_model$frequency, 0.99), length.out = 100))

  predictions <- posterior_linpred(fit_ws, newdata = newdata, transform = TRUE, re.form = NA)

  newdata$prob_mean <- apply(predictions, 2, mean)
  newdata$prob_sd <- apply(predictions, 2, sd)
  
  p_responsecurve_heatmap <- ggplot(newdata, aes(x = severity, y = frequency, fill = prob_mean)) +
    geom_tile() +
    scale_fill_continuous(low = "#fff5f0", high = "#67000d", limit = c(0, 1)) +
    theme_bw() +
    theme(legend.position = c(0, 1),
          legend.justification = c(-0.05, 1.1),
          legend.box.background = element_rect(colour = "black"),
          legend.margin = margin(2, 6, 2, 6),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8)) +
    labs(x = "Disturbance severity", fill = paste0("P(", process, ")"), 
         y = "Disturbance frequency") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    guides(fill = guide_colorbar(barwidth = 6, barheight = 0.5,
                                 direction = "horizontal", title.position = "top"))
  
  # Gather results and add to list
  
  results[[k]] <- list(fit_ws, #1
                       subset_size, #2 
                       predictors_selected_names, #3
                       p_varsel, #4
                       fit_ws_reduced, #5 
                       elpd, #6
                       auc, #7
                       p_roc, #8
                       estimates, #9
                       randomeffects, #10
                       p_probabilitymap_all, #11
                       p_probabilitymap_difference, #12
                       p_responsecurve_classical, #13
                       p_responsecurve_heatmap) #14
  
}



save(results, file = "../../../../../results/results.RData")

