library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(broom)

### 1. Load inferred and observed bipartite food-web data
# inferred_foodweb_properties: per-site network metrics produced by the food-web metrics script.
# observed_foodweb_properties: per-site network metrics produced by the food-web metrics script.
# The inferred (metaweb-based) and independently observed food-web properties are each tested against the same environmental covariates.
inferred_foodweb_properties <- read.csv("data/02_processed/fw_properties.csv") %>% filter(type == "Bipartite")
observed_foodweb_properties <- read.csv("data/02_processed/observed_fw_properties.csv")
covariates <- read.csv("data/01_raw/covariates.csv")

inferred_foodweb_properties <- merge(inferred_foodweb_properties, covariates)
observed_foodweb_properties <- merge(observed_foodweb_properties, covariates)

combined_foodweb_properties <- bind_rows(mutate(inferred_foodweb_properties, data_source = "Metaweb"),
                                         mutate(observed_foodweb_properties, data_source = "Observed"))

### 2. Define structural metrics and environmental predictors
foodweb_metrics <- c("connectance",
                     "nestedness",
                     "modularity",
                     "r50",
                     "median_degree_z",
                     "skew_degree_z",
                     "cv_degree_z")

environmental_covariates <- c("impervious_1000",
                              "no.plants",
                              "prop.nonnative",
                              "prop.wood",
                              "maint_disturbance",
                              "maint_time")

### 3. Fit metric-covariate models for inferred and observed food webs
# For each food-web data source, metric, and covariate, fit a linear model and  extract the covariate coefficient and 95% confidence interval. 
coefficient_data <- combined_foodweb_properties %>%
  pivot_longer(cols = all_of(foodweb_metrics),
               names_to = "metric",
               values_to = "metric_value") %>%
  pivot_longer(cols = all_of(environmental_covariates),
               names_to = "covariate",
               values_to = "covariate_value") %>%
  mutate(covariate = factor(covariate, levels = rev(environmental_covariates)), 
         metric = factor(metric, levels = foodweb_metrics)) %>%
  group_by(covariate, data_source) %>%
  mutate(covariate_scaled = as.numeric(scale(covariate_value))) %>%
  group_by(data_source, metric, covariate) %>%
  do(tidy(lm(metric_value ~ covariate_scaled, data = .),
          conf.int = TRUE,
          conf.level = 0.95)) %>%
  filter(term == "covariate_scaled") %>%
  ungroup() %>%
  mutate(significance = case_when(p.value < 0.001 ~ "***",
                                  p.value < 0.01 ~ "**",
                                  p.value < 0.05 ~ "*",
                                  p.value < 0.1 ~ "·",
                                  TRUE ~ ""),
         significance_category = if_else(p.value < 0.1, "Significant", "NS"))

### 4. Compare effect directions and statistical significance
# Each metric-covariate relationship is classified according to whether inferred and observed networks show the same effect direction and significance pattern.
direction_significance_data <- coefficient_data %>%
  mutate(significant_raw = p.value < 0.05,
         direction = case_when(estimate > 0 ~ "Positive",
                               estimate < 0 ~ "Negative",
                               estimate == 0 ~ "Zero")) %>%
  select(metric, covariate, data_source, estimate, p.value, significant_raw, direction) %>%
  pivot_wider(names_from = data_source,
              values_from = c(estimate, p.value, significant_raw, direction),
              names_sep = "_") %>%
  mutate(direction_agreement = case_when(is.na(direction_Metaweb) | is.na(direction_Observed) ~ NA_character_,
                                         direction_Metaweb == direction_Observed ~ "Same trend",
                                         TRUE ~ "Different trend"),
         significance_pattern = case_when(significant_raw_Metaweb & significant_raw_Observed ~ "both significant",
                                          significant_raw_Observed & !significant_raw_Metaweb ~ "observed-only significant",
                                          significant_raw_Metaweb & !significant_raw_Observed ~ "inferred-only significant",
                                          !significant_raw_Metaweb & !significant_raw_Observed ~ "Neither significant"),
         agreement_class = case_when(significance_pattern == "Neither significant" ~ "Neither significant",
                                     TRUE ~ paste(direction_agreement, significance_pattern, sep = ", ")))

direction_significance_summary <- direction_significance_data %>%
  filter(!is.na(agreement_class)) %>%
  count(agreement_class) %>%
  mutate(agreement_class = factor(agreement_class,
                                  levels = c("Neither significant",
                                             "Same trend, both significant",
                                             "Same trend, observed-only significant",
                                             "Same trend, inferred-only significant",
                                             "Different trend, both significant",
                                             "Different trend, observed-only significant",
                                             "Different trend, inferred-only significant")))

### 5. Figure S7: agreement between inferred and observed effects
ggplot(direction_significance_summary, aes(x = agreement_class, y = n)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.4, size = 5) +
  theme_pubr(base_size = 14) +
  labs(x = NULL,
       y = "Number of metric-covariate tests") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   vjust = 1))
