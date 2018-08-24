### This is for plotting the results of dis~haz_analysis 
### written by J.Sebald & C. Senf

# load packages, set working directory to either local (pc/macbook) or server (valenor)

{
  library(raster)# version 2.6-7
  library(sp)# version 1.3-1
  library(rgdal)# version 1.2-20
  library(igraph)# version 1.2.1
  library(tidyverse)# version 1.2.1
  library(projpred)
  library(ggthemes)
  library(rstanarm)
  
  rm(list=ls())
  
}



# Load data ---------------------------------------------------------------

load(file = "../../../../../results/results.RData")

data <- read.csv("../data/data_for_model.csv")

shp_ws <- raster::shapefile("../data/shp_ws.shp")
shp_ws_fortify <- broom::tidy(shp_ws, region = "WLK_ID")

# Set names ---------------------------------------------------------------

processes <- c("DFLOOD", "DFLOW", "FST")

vars_ws <- data.frame(varname = c("h_mean", "Melton", "Elevation", "Circularit", "Elongation", "artifical", "forest", "area", "patchdensity", 
                                  "severity", 
                                  "frequency",
                                  "severity:frequency"), 
                      name = c("Elevation", "Melton ratio", "Elevation ratio", "Circularity", "Elongtion", "Artificial", "Forest", "Area", "Patch density", 
                               "Severity", 
                               "Frequency",
                               "Severity x Frequency"),
                      stringsAsFactors = FALSE)

# Plotting ----------------------------------------------------------------

subset_size <- results %>% 
  map(~ .[[2]])

p_varsel <- results %>% 
  map(~ .[[4]]$data) %>%
  set_names(processes) %>%
  bind_rows(.id = "process") %>%
  ggplot(., aes(x = size, y = value)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = value - se, ymax = value + se), width = 0) +
  facet_wrap(~ process, scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank()) +
  labs(x = "Number of predictors", y = "Expected log predictive density") +
  geom_vline(data = data.frame(process = processes, n_select = unlist(subset_size)),
             aes(xintercept = n_select), col = scales::muted("red"), linetype = "dashed")

ggsave("variable_selection.pdf", p_varsel, path = "results/", width = 7.5, height = 2.5)
ggsave("variable_selection.png", p_varsel, path = "results/", width = 7.5, height = 2.5)


predictors_selected <- results %>% 
  map(~ .[[3]]$name) %>%
  map(~ paste(., collapse = ", "))

elpd <- results %>% 
  map(~ .[[6]])

auc <- results %>% 
  map(~ .[[7]])

watershed_model_results_table <- data.frame(process = processes,
                                            subset_size = unlist(subset_size),
                                            predictors_selected = unlist(predictors_selected),
                                            elpd = unlist(elpd),
                                            auc = unlist(auc))

write_csv(watershed_model_results_table, "results/watershed_model_results_table.csv")

p_roc <- results %>% 
  map(~ .[[8]]$data) %>%
  set_names(processes) %>%
  bind_rows(.id = "process") %>%
  ggplot(., aes(x = fpr, y = tpr)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", col = scales::muted("red")) +
  facet_wrap(~ process, scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank()) +
  labs(x = "False positive rate", y = "True positive rate")

ggsave("roc_curves.pdf", p_roc, path = "results/", width = 7.5, height = 2.5)
ggsave("roc_curves.png", p_roc, path = "results/", width = 7.5, height = 2.5)


p_estimate <- results %>% 
  map(~ .[[9]]) %>%
  bind_rows() %>%
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

ggsave("estimates.pdf", p_estimate, path = "results/", width = 7.5, height = 2.5)
ggsave("estimates.png", p_estimate, path = "results/", width = 7.5, height = 2.5)


p_maps <- results %>% 
  map(~ .[[10]]$data) %>%
  set_names(processes) %>%
  bind_rows(.id = "process") %>%
  ggplot(.) +
  aes(long, lat, group = group, fill = prob_hazard_mean) +
  geom_polygon() +
  coord_equal() +
  scale_fill_gradient(low = "#fff5f0", high = "#67000d") +
  labs(x = NULL, y = NULL, fill = paste0("Probability of event")) +
  ggthemes::theme_map() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        strip.background = element_blank()) +
  guides(fill = guide_colorbar(barwidth = 10, barheight = 0.5,
                               direction = "horizontal", title.position = "top")) +
  facet_wrap(~process)

ggsave("probability_maps.pdf", p_maps, path = "results/", width = 7.5, height = 2.5)
ggsave("probability_maps.png", p_maps, path = "results/", width = 7.5, height = 2.5)


p_maps_difference <- results %>% 
  map(~ .[[11]]$data) %>%
  set_names(processes) %>%
  bind_rows(.id = "process") %>%
  ggplot(.) +
  aes(long, lat, group = group, fill = prob_hazard_mean) +
  geom_polygon() +
  coord_equal() +
  scale_fill_gradient2() +
  labs(x = NULL, y = NULL, fill = paste0("Change in event probability after including disturbances")) +
  ggthemes::theme_map() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        strip.background = element_blank()) +
  guides(fill = guide_colorbar(barwidth = 10, barheight = 0.5,
                               direction = "horizontal", title.position = "top")) +
  facet_wrap(~process)

ggsave("probability_maps_difference.pdf", p_maps_difference, path = "results/", width = 7.5, height = 2.5)
ggsave("probability_maps_difference.png", p_maps_difference, path = "results/", width = 7.5, height = 2.5)


p_response_classic <- results %>% 
  map(~ .[[12]]$data) %>%
  set_names(processes) %>%
  bind_rows(.id = "process") %>%
  #mutate(frequency = factor(frequency, labels = c("High (+1SD)", "Average", "Low (-1SD)"))) %>%
  ggplot(., aes(x = severity, y = prob_mean)) +
  geom_ribbon(aes(ymin = prob_mean - prob_sd * 2, ymax = prob_mean + prob_sd * 2, fill = frequency), alpha = 0.3) +
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

ggsave("responsecurves_classic.pdf", p_response_classic, path = "results/", width = 7.5, height = 2.25)
ggsave("responsecurves_classic.png", p_response_classic, path = "results/", width = 7.5, height = 2.25)



p_response_heatmap <- results %>% 
  map(~ .[[13]]$data) %>%
  set_names(processes) %>%
  bind_rows(.id = "process") %>%
  ggplot(., aes(x = severity, y = frequency * (-1))) +
  geom_tile(aes(fill = prob_mean)) +
  # geom_point(data = data %>% mutate(severity = as.double(scale(severity)),
  #                                   frequency = as.double(scale(frequency))) %>%
  #              sample_n(1000), 
  #                 aes(x = severity, y = frequency * (-1)), alpha = 0.3) +
  scale_fill_continuous(low = "#fff5f0", high = "#67000d", limit = c(0, 1)) +
  theme_bw() +
  theme(legend.position = c(0, 0),
        legend.justification = c(-0.05, -0.1),
        legend.box.background = element_rect(colour = "black"),
        legend.margin = margin(2, 6, 2, 6),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        panel.grid = element_blank(),
        strip.background = element_blank()) +
  labs(x = "Disturbance severity", fill = paste0("Probability of event"), 
       y = "Disturbance frequency") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill = guide_colorbar(barwidth = 6, barheight = 0.5,
                               direction = "horizontal", title.position = "top")) +
  facet_wrap(~process)


### Correlations

prob_count <- results %>%
  map(~ posterior_linpred(.[[1]], transform = TRUE), re.form = NA) %>%
  map2(.y = list(data$DFLOOD, data$DFLOW, data$FST),
       ~ data.frame(pred = apply(.x, 2, mean), count = .y)) %>%
  set_names(processes) %>%
  bind_rows(.id = "process") %>%
  ggplot(., aes(x = factor(count), y = pred)) +
  geom_violin(fill = "grey") +
  facet_wrap(~process)


ggsave("correlations.pdf", p_response_classic, path = "results/", width = 7.5, height = 2.25)
ggsave("correlations.png", p_response_classic, path = "results/", width = 7.5, height = 2.25)


