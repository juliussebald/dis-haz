
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
  library(patchwork)
  
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = 4) 
  
}

# Load data ---------------------------------------------------------------

data <- read.csv("../data/tables/data_for_model.csv")

processes <- c("DFLOOD", "DFLOW", "FST")


# Selecting a model -----------------------------------------------

# Create dataframe with variable short names and variable long names

vars_ws <- data.frame(varname = c("h_mean", "Melton", "Elevation", "Circularit", "Elongation", "artifical", "forest", "area", "patchdensity", 
                                  "extent", 
                                  "pulse",
                                  "extent:pulse"), 
                      name = c("Elevation", "Melton ratio", "Elevation ratio", "Circularity", "Elongtion", "Artificial", "Forest", "Area", "Patch density", 
                               "Extent", 
                               "Pulse",
                               "Extent x Pulse"),
                      stringsAsFactors = FALSE)

# Loop through processes and calibrate varying models

models <- vector("list", length = 3)

k <- 0

for (process in processes) {
  
  k <- k + 1
  
  # Bring data into form
  
  vars_nointeraction <- vars_ws %>% 
    filter(varname != c("extent:pulse"))
  
  data_model <- data
  data_model[data_model$extent == 0, "pulse"] <- NA
  data_model <- data_model %>%
    mutate_at(.vars = vars(c(vars_nointeraction$varname)), function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  data_model[, "response"] <- ifelse(data_model[, process] > 0, 1, 0)
  data_model[data_model$extent == min(data_model$extent), "pulse"] <- 0
  
  # Fit watershed-only model
  
  fit_ws_only <- stan_glm(as.formula(paste0("response ~ ", paste0(paste(vars_ws$varname[-which(vars_ws$varname %in% c("extent", "pulse", "extent:pulse"))], collapse = "+")))),
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
  
  
  # Include disturbance predictors and compare to watershed-only model
  
  fit_full_exp <- update(fit_ws_only_reduced, . ~ . + extent * pulse)

  loo_fit_full_exp <- loo(fit_full_exp)
  loo_fit_ws_only_reduced  <- loo(fit_ws_only_reduced)

  model_comparison <- loo::compare(loo_fit_ws_only_reduced,
                                   loo_fit_full_exp
                                   )
  
  # Store everything in a list
  
  models[[k]] <- list(fit_ws_only_reduced, #1
                      fit_full_exp, #2
                      list(loo_fit_ws_only_reduced,
                           loo_fit_full_exp), #3
                      model_comparison) #4
                       
}


save(models, file = "../results/models.RData")
#load(file = "../results/models.RData")


# Select final model -------------------------------------------------------------

### Compare models

model_comparison <- models %>% 
  map(~ as.data.frame(.[[4]]) %>%
        rownames_to_column(., var = "model")) %>%
  set_names(processes) %>%
  bind_rows(.id = "process")

write_csv(model_comparison, "../results/model_comparison.csv")

### Decide for final model

final_models <- models %>% map(., ~.[[2]])


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
                                    pulse = c(quantile(data_model$pulse, 0.05), 0, quantile(data_model$pulse, 0.95)))  # viele EZGs haben Pulisty kleiner als -1

predictions <- final_models %>%
  map(~ posterior_linpred(., newdata = response_disturbance, transform = TRUE, re.form = NA))

for (i in 1:length(predictions)) {
  response_disturbance[, paste0(processes[i], "_mean")] <- apply(predictions[[i]], 2, mean)
  response_disturbance[, paste0(processes[i], "_sd")] <- apply(predictions[[i]], 2, sd)  
}

p_response <- response_disturbance %>%
  mutate(pulse = factor(pulse, labels = c("High Press", "Average", "High Pulse"))) %>%
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
  scale_color_manual(values = c(scales::muted("red"), "grey", scales::muted("blue"))) +
  scale_fill_manual(values = c(scales::muted("red"), "grey", scales::muted("blue"))) +
  facet_wrap(~process) +
  labs(x = "Disturbance extent", y = "Probability of event", 
       fill = "Disturbance type", col = "Disturbance type")

ggsave("response_curves_20182008.png", p_response, path = "../results/", width = 7.5, height = 2.5)







