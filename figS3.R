library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(broom)

### 1. Load food-web properties and environmental covariates
# foodweb_properties: per-site network metrics produced by the food-web metrics script.
# covariates: per-site covariates.
foodweb_properties <- read.csv("data/02_processed/fw_properties.csv")
covariates <- read.csv("data/01_raw/covariates.csv")

foodweb_properties <- merge(foodweb_properties, covariates)

# Shared colour scale for the three focal network types.
type_colors <- c("Combined" = "black",
                 "Bipartite" = "#2166ac",
                 "Multitrophic" = "#b2182b")

### 2. Fit standardized degree-metric responses to impervious cover
# For every network type, degree metric, and spatial scale of impervious cover, fit a linear model using standardized predictor and response values.
# Extract the standardized coefficient, its 95% confidence interval, and significance.
coefficient_data <- foodweb_properties %>%
  select(ugs, type, impervious_1000, impervious_500, impervious_150, cv_degree_z, skew_degree_z, median_degree_z) %>%
  pivot_longer(cols = ends_with("_z"),
               names_to = "metric",
               values_to = "metric_value") %>%
  pivot_longer(cols = c(impervious_1000, impervious_500, impervious_150),
               names_to = "covariate",
               values_to = "covariate_value") %>%
  mutate(covariate = factor(covariate, levels = c("impervious_150", "impervious_500", "impervious_1000")),
         metric = factor(metric, levels = c("median_degree_z", "skew_degree_z", "cv_degree_z")),
         covariate = recode(covariate, impervious_1000 = "% impervious\n(1000 m)", impervious_500 = "% impervious\n(500 m)", impervious_150 = "% impervious\n(150 m)"),
         metric = recode(metric, cv_degree_z = "Coefficient of variation (z)", median_degree_z = "Median (z)", skew_degree_z = "Skewness (z)")) %>%
  group_by(covariate, type) %>%
  mutate(covariate_scaled = scale(covariate_value)) %>%
  group_by(type, metric, covariate) %>%
  do(tidy(lm(scale(metric_value) ~ covariate_scaled, data = .),
          conf.int = TRUE,
          conf.level = 0.95)) %>%
  filter(term == "covariate_scaled") %>%
  mutate(significance = case_when(p.value < 0.001 ~ "***",
                                  p.value < 0.01 ~ "**",
                                  p.value < 0.05 ~ "*",
                                  p.value < 0.1 ~ "·",
                                  TRUE ~ ""),
    significance_category = if_else(p.value < 0.1, "Significant", "NS"),
    star_x = if_else(estimate > 0, conf.high, conf.low),
    star_hjust = if_else(estimate > 0, -0.5, 1.5),
    type = factor(type, levels = c("Combined", "Bipartite", "Multitrophic")))

# 3. Figure S3: effects of impervious cover on degree distributions
ggplot(coefficient_data, aes(x = estimate, y = covariate, color = type, fill = type, alpha = significance_category)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
                width = 0,
                position = position_dodge(width = 0.6)) +
  geom_text(aes(label = significance, x = star_x, hjust = star_hjust),
            position = position_dodge(width = 0.6),
            vjust = 0.7,
            size = 5,
            show.legend = FALSE) +
  scale_alpha_manual(values = c("Significant" = 1, "NS" = 0.2), name = NULL) +
  scale_color_manual(values = type_colors, name = NULL) +
  scale_fill_manual(values = type_colors, name = NULL) +
  facet_wrap(~ metric, scales = "free_x") +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  theme_pubr(base_size = 16) +
  labs(y = "Standardized coefficient and 95% CI",
       x = NULL) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
