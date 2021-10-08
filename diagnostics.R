# Libraries ---------------------------------------------------------------

library(mvnfast)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

# Posterior samples -------------------------------------------------------

posterior_samples <- read_rds("data/posterior_samples.rds")

Lout <- posterior_samples$Lout
Pout <- posterior_samples$Pout
Oout <- posterior_samples$Oout

# Traceplots --------------------------------------------------------------

as_tibble(Lout) %>% 
  mutate(iter = row_number()) %>% 
  pivot_longer(-iter) %>% 
  ggplot(aes(iter, value)) +
  geom_line() +
  facet_wrap(~name) +
  labs(title = "Factor Loadings") +
  theme_bw()

ggsave("figures/factor_loadings_traceplot.png")

as_tibble(Pout) %>% 
  mutate(iter = row_number()) %>% 
  pivot_longer(-iter) %>% 
  ggplot(aes(iter, value)) +
  geom_line() +
  facet_wrap(~name) +
  labs(title = "Error Precision") +
  theme_bw()

ggsave("figures/error_precision_traceplot.png")

# Posterior summaries -----------------------------------------------------

as_tibble(Lout) %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(name, value)) +
  geom_boxplot() +
  geom_point(
    data = tibble(truth = c(Lt), name = str_c("V", 1:length(truth))),
    aes(name, truth),
    color = "red"
  ) +
  labs(title = "Factor Loadings") +
  theme_bw()

ggsave("figures/factor_loadings_posterior.png")

as_tibble(1 /  Pout) %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(name, value)) +
  geom_boxplot() +
  geom_point(
    data = tibble(truth = diag(St), name = str_c("V", 1:length(truth))),
    aes(name, truth),
    color = "red"
  ) +
  labs(title = "Error Precision") +
  theme_bw()

ggsave("figures/error_precision_posterior.png")

as_tibble(Oout) %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(name, value)) +
  geom_boxplot() +
  geom_point(
    data = tibble(truth = c(Ot), name = str_c("V", 1:length(truth))),
    aes(name, truth),
    color = "red"
  ) +
  labs(title = "Omega") +
  theme_bw()

ggsave("figures/Omega_posterior.png")

