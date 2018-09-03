
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
                                  "rel_years",
                                  "extent:rel_years"), 
                      name = c("Elevation", "Melton ratio", "Elevation ratio", "Circularity", "Elongtion", "Artificial", "Forest", "Area", "Patch density", 
                               "Extent", 
                               "Frequency",
                               "Extent x Frequency"),
                      stringsAsFactors = FALSE)

# Loop through processes and calibrate varying models

models <- vector("list", length = 3)

k <- 0

for (process in processes) {
  
  k <- k + 1
  
  # Bring data into form
  
  vars_nointeraction <- vars_ws %>% 
    filter(varname != c("extent:rel_years"))
  
  data_model <- data
  data_model[data$extent == 0, "rel_years"] <- mean(data_model$rel_years)
  data_model <- data_model %>%
    mutate_at(.vars = vars(c(vars_nointeraction$varname)), function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  data_model[, "response"] <- ifelse(data_model[, process] > 0, 1, 0)
  
  # Fit watershed-only model
  
  fit_ws_only <- stan_glm(as.formula(paste0("response ~ ", paste0(paste(vars_ws$varname[-which(vars_ws$varname %in% c("extent", "rel_years", "extent:rel_years"))], collapse = "+")))),
                          data = data_model, 
                          family = binomial(link = "logit"), 
                          prior_intercept = normal(0|1), prior = normal(0|1))
  
  # Select important watershed predictors
  
  ws_varsel <- varsel(fit_ws_only)
  
  subset_size <- suggest_size(ws_varsel, alpha = 0.05)
  ws_predictors_selected <- names(ws_varsel$varsel$vind[1:subset_size])
  
  # Fit reduced watershed-only model
  
  fit_ws_only_reduced <- stan_glm(as.formula(paste0("response ~ ", paste0(paste(ws_predictors_selected, collapse = "+")))),
                                  data = data_model,
                                  family = binomial(link = "logit"),
                                  prior_intercept = normal(0|1), prior = normal(0|1))
  
  # Test if including forestry areas as random effect improves model
  
  fit_ws_only_reduced_rf <- stan_glmer(as.formula(paste0("response ~ ", paste0(paste(ws_predictors_selected, collapse = "+"), "+ (1|eco_unit)"))),
                                       data = data_model,
                                       family = binomial(link = "logit"),
                                       prior_intercept = normal(0|1), prior = normal(0|1))
  
  loo_fit_ws_only_reduced <- loo(fit_ws_only_reduced)
  loo_fit_ws_only_reduced_rf <- loo(fit_ws_only_reduced_rf)
  
  comp <- loo::compare(loo_fit_ws_only_reduced, loo_fit_ws_only_reduced_rf)
  
  if ((comp[1] + comp[2]) > 0) {
    loo_fit_ws_only_reduced <- loo_fit_ws_only_reduced_rf
    fit_ws_only_reduced <- fit_ws_only_reduced_rf
  }
  
  # Include disturbance predictors and compare to watershed-only model
  
  fit_full_exp <- update(fit_ws_only_reduced, . ~ . + extent * rel_years)
  fit_full_e_p <- update(fit_ws_only_reduced, . ~ . + extent + rel_years)
  
  loo_fit_full_exp <- loo(fit_full_exp)
  loo_fit_full_e_p <- loo(fit_full_e_p)
  
  model_comparison <- loo::compare(loo_fit_ws_only_reduced,
                                   loo_fit_full_exp,
                                   loo_fit_full_e_p)
  
  # Store everything in a list
  
  models[[k]] <- list(fit_full_exp, #1
                      fit_full_e_p, #2
                      list(loo_fit_ws_only_reduced,
                           loo_fit_full_exp,
                           loo_fit_full_e_p), #3
                      model_comparison) #4
                       
}

save(models, file = "../results/rel_years_model/models_rel_years.RData")


# Select final model -------------------------------------------------------------

# Compare models

model_comparison <- models %>% 
  map(~ as.data.frame(.[[4]]) %>%
        rownames_to_column(., var = "model")) %>%
  set_names(processes) %>%
  bind_rows(.id = "process")

write_csv(model_comparison, "../results/rel_years_model/model_comparison_rel_years.csv")

models %>%
  map(~ loo::compare(.[[3]][[3]], .[[3]][[2]]))

# Decide for final model

final_models <- models %>% map2(.y = c(1, 2, 2), ~ .[[.y]])

# Calculate LOO AUC for final model

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
  write_csv(., "../results/model_auc_rel_years.csv")


# Create plots and tables ---------------------------------------------------------


# Extract and plot estimates

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

ggsave("estimates_rel_years.pdf", p_estimates, path = "../results/rel_years_model/", width = 7.5, height = 2.5)


# Extract and plot random effects

randomeffects <- final_models %>%
  map(~ as.data.frame(.) %>%
        dplyr::select(matches("Intercept")) %>% # Select everything that is an intercept
        gather(key = varname, value = value)) %>%
  set_names(processes) %>%
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

ggsave("randomeffects.pdf", p_randomeffects, path = "../results/rel_years_model/", width = 7.5, height = 2.5)


# Create response curve plots

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
                                    rel_years = c(-1, 0, 1))

predictions <- final_models %>%
  map(~ posterior_linpred(., newdata = response_disturbance, transform = TRUE, re.form = NA))

for (i in 1:length(predictions)) {
  response_disturbance[, paste0(processes[i], "_mean")] <- apply(predictions[[i]], 2, mean)
  response_disturbance[, paste0(processes[i], "_sd")] <- apply(predictions[[i]], 2, sd)  
}

p_response <- response_disturbance %>%
  mutate(rel_years = factor(rel_years, labels = c("Low (-1SD)", "Average", "High (+1SD)"))) %>%
  gather(key = process, value = prop, -(h_mean:rel_years)) %>%
  separate(process, c("process", "metric"), "_") %>%
  spread(key = metric, value = prop) %>%
  ggplot(., aes(x = extent, y = mean)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = rel_years), alpha = 0.3) +
  geom_line(aes(col = rel_years)) +
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
  scale_color_manual(values = c(scales::muted("blue"), "grey", scales::muted("red"))) +
  scale_fill_manual(values = c(scales::muted("blue"), "grey", scales::muted("red"))) +
  facet_wrap(~process) +
  labs(x = "Disturbance extent", y = "Probability of event", 
       fill = "Frequency", col = "Frequency")

ggsave("response_curves_rel_years.pdf", p_response, path = "../results/rel_years_model/", width = 7.5, height = 2.5)

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

