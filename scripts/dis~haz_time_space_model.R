
library(tidyverse)
library(raster)
library(data.table)
library(rstanarm)

options(mc.cores = parallel::detectCores())


# Prepare data ------------------------------------------------------------

# Load disturbance, watershep-id and forest raster

disturbance <- raster("../data/rasters/disturbance.tif")
watersheds <- raster("../data/rasters/raster_ws.tif")
forest <- raster("../data/rasters/forest.tif")

# Load watershed characteristics and remove aggregated event information (i.e., event counts per ws)

ws_characteristics <- read.csv("../data/tables/data_for_model.csv") %>%
  dplyr::select(-DFLOOD, -DFLOW, -FST)

# Load raw event file and aggregate to number of events per year per watershed

events <- read.csv("../rawdata/GroupEvents.csv") %>%
  rename(., year = ereignis_jahr) %>%
  mutate(., event = case_when(
    prozessart == "Fluviatiler Feststofftransport" ~ "FST",
    prozessart ==  "Murartiger Feststofftransport" ~ "DFLOOD",
    prozessart == "Murgang" ~ "DFLOW")) %>%
  dplyr::select(WLK_ID, year, event) %>%
  group_by(WLK_ID, year, event) %>%
  summarize(n = n()) %>%
  spread(key = event, value = n)

# Create dataset containing disturbances, forestcover and watershed-ids

dat <- data.table(disturbance = values(disturbance),
                  WLK_ID = values(watersheds),
                  forest = values(forest)) %>%
  na.omit()

# Aggregate disturbances

disturbance <- dat %>%
  filter(disturbance > 0) %>%
  group_by(WLK_ID, year = disturbance) %>%
  summarize(disturbance = n())

# Aggregate forest dover

forest <- dat %>%
  group_by(WLK_ID) %>%
  summarize(forest_count = sum(forest))

# Combine everthing into one data frame, join watershed characteristics, calculate disturbance rate, join events

dat_full <- expand.grid(year = 1986:2016,
                        WLK_ID = unique(dat$WLK_ID)) %>%
  left_join(disturbance, by = c("WLK_ID", "year")) %>%
  left_join(forest, by = "WLK_ID") %>%
  left_join(ws_characteristics, by = "WLK_ID") %>%
  mutate(disturbance = ifelse(is.na(disturbance), 0, disturbance),
         disturbance_rate = disturbance / forest_count) %>%
  left_join(events, by = c("year", "WLK_ID")) %>%
  mutate_at(.vars = vars(DFLOOD, DFLOW, FST), ~ ifelse(is.na(.), 0, .))

# Calculate lagged disturbance values

dat_full <- dat_full %>%
  split(.$WLK_ID) %>%
  map(~ mutate(., disturbance_rate_lag1 = lag(disturbance_rate, 1),
               disturbance_rate_lag2 = lag(disturbance_rate, 2),
               disturbance_rate_lag3 = lag(disturbance_rate, 3),
               disturbance_rate_lag4 = lag(disturbance_rate, 4),
               disturbance_rate_lag5 = lag(disturbance_rate, 5),
               disturbance_rate_lag6 = lag(disturbance_rate, 6),
               disturbance_rate_lag7 = lag(disturbance_rate, 7),
               disturbance_rate_lag8 = lag(disturbance_rate, 8),
               disturbance_rate_lag9 = lag(disturbance_rate, 9),
               disturbance_rate_lag10 = lag(disturbance_rate, 10),
               disturbance_rate_lag11 = lag(disturbance_rate, 11),
               disturbance_rate_lag12 = lag(disturbance_rate, 12),
               disturbance_rate_lag13 = lag(disturbance_rate, 13),
               disturbance_rate_lag14 = lag(disturbance_rate, 14),
               disturbance_rate_lag15 = lag(disturbance_rate, 15))) %>%
  bind_rows()

# Calculate summed disturbance rate over variable lag distances

for (i in 1:15) {
  dat_full[, paste0("disturbance_rate_lag", 0, "_to_lag", i)] <- rowSums(as.matrix(dat_full[, paste0("disturbance_rate_lag", 1:i)])) + dat_full$disturbance_rate
}  

# Write to disc

write_csv(dat_full, "../data/tables/data_for_model_temporal.csv")


# Model -------------------------------------------------------------------

# Load data

dat_full <- read.csv("../data/tables/data_for_model_temporal.csv")

# Standardize predictors and create binary response

dat_full <- dat_full %>%
  mutate_at(.vars = vars(area, Melton, Circularit, Elongation, Elevation, h_mean, forest, patchdensity, disturbance_rate:disturbance_rate_lag0_to_lag15),
            .funs = function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)) %>%
  mutate_at(.vars = vars(DFLOOD, DFLOW, FST), .funs = function(x) ifelse(x > 0, 1, 0))

# Loop through processes and calculate models for variable lag distances

#process <- "DFLOOD"
processes <- c("DFLOOD", "DFLOW", "FST")

results <- vector("list", 3)

k <- 0

for (process in processes) {
  
  k <- k + 1
  
  # Sub-sample response to maximum number of 1 and double the amount for 0
  
  dat_full_sample <- dat_full %>%
    split(.[,process]) %>%
    map2(.y = (table(dat_full[,process])[2] * c(2, 1)), ~ sample_n(., .y)) %>%
    bind_rows()
  
  # Load important predictors from spatial analysis
  
  # -> TBD
  
  predictors <- c("area + Melton + Circularit + Elongation + Elevation + h_mean + forest + patchdensity")
  
  # List for storing the posterior of disturbance and LOO-ELDP with variable lags
  
  posterior_disturbance <- vector("list", 15)
  elpd_disturbance <- vector("list", 15)
  
  # Loop through lags
  
  #start.time <- Sys.time()
  
  for (i in 1:15) {
    
    print(paste0("Lag: ", i))
    
    # Fit models for variable lags
    
    formular <- as.formula(paste0(process, " ~ ", predictors, " + disturbance_rate_lag0_to_lag", i," + (1 | WLK_ID)"))
    fit <- stan_glmer(formular, data = dat_full_sample, family = binomial(link = "logit"))
    
    # Extract posterior and store in list
    
    posterior <- as.matrix(fit)
    posterior_disturbance[[i]] <- posterior[, paste0("disturbance_rate_lag0_to_lag", i)]
    
    # Calculate LOO-ELDP and store in list
    #elpd <- kfold(fit, K = 10)
    #elpd_disturbance[[i]] <- elpd
    
  }
  
  #end.time <- Sys.time()
  #print(paste0("Total runinning time for 15 lags: ", end.time - start.time))
  
  # Create data frame storing posteriors
  
  posterior_disturbance <- posterior_disturbance %>%
    map(~ data.frame(value = .)) %>%
    set_names(1:15) %>%
    bind_rows(.id = "lag")
  
  # Store results in list
  
  results[[k]] <- list(posterior_disturbance
                       #elpd
                       )
  
}

save(results, file = "../results/results_temporal.RData")

posterior_disturbance_lags <- results %>%
  map(~.[[1]]) %>%
  set_names(processes) %>%
  bind_rows(.id = "process")

posterior_disturbance_lags <- posterior_disturbance_lags %>%
  group_by(lag, process) %>%
  summarize(above_zero = mean(value > 0)) %>%
  right_join(posterior_disturbance_lags, by = c("lag", "process"))

p <- posterior_disturbance_lags %>%
  mutate(process = factor(process, levels = c("FST", "DFLOOD", "DFLOW"))) %>%
  ggplot(., 
       aes(x = factor(lag, levels = 1:15), y = value)) +
  geom_violin(fill = "grey") +
  geom_hline(yintercept = 0) +
  facet_wrap(~process, scales = "free") +
  theme_bw() +
  labs(x = "Number of years before event year",
       y = "Posterior probability of\ndisturbance extent")
  
ggsave("temporal_analysis_disturbances.pdf", p, path = "../results/", width = 7.5, height = 2.5)
  
  

