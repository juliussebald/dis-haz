
library(tidyverse)
library(rstanarm)

options(mc.cores = parallel::detectCores())

setwd("Dropbox/Research/European Forests/Hazards/")

data <- read.csv("data_for_model_08162018.csv")

data <- data %>%
  mutate_at(.vars = vars(area, Elevation, h_mean, forest, severity, frequency),
            .funs = function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))

fit_poisson <- stan_glm(FST ~ area + Elevation + h_mean + forest + severity * frequency, 
                data = data, family = poisson)

fit_negbinom <- stan_glm(FST ~ area + Elevation + h_mean + forest + severity * frequency, 
                data = data, family = neg_binomial_2)

prop_zero <- function(y) mean(y == 0)

prop_zero_poisson <- pp_check(fit_poisson, plotfun = "stat", stat = "prop_zero")
prop_zero_negbinom <- pp_check(fit_negbinom, plotfun = "stat", stat = "prop_zero")

gridExtra::grid.arrange(prop_zero_poisson + ggtitle("Poisson"), 
                        prop_zero_negbinom + ggtitle("Negative Binomial"), 
                        ncol = 2)

loo::compare(loo(fit_poisson), loo(fit_negbinom))

bayesplot::mcmc_areas(as.matrix(fit_negbinom))

summary(fit_negbinom)

varsel <- projpred::varsel(fit_poisson)

projpred::varsel_plot(varsel)
varsel$varsel$vind

projpred::varsel_stats(varsel)

projpred::suggest_size(varsel)


