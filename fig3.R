library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(broom)

### 1. Load and merge data
# foodweb_properties: per-site network metrics produced by the food-web metrics script.
# covariates: per-site covariates.
foodweb_properties <- read.csv("data/02_processed/fw_properties.csv")
covariates <- read.csv("data/01_raw/covariates.csv")
foodweb_properties <- merge(foodweb_properties, covariates)

# 2. Plot aesthetics
type_colors <- c("Combined" = "black",
                 "Bipartite" = "#2166ac",
                 "Multitrophic" = "#b2182b")

covariate_labels <- c(impervious_1000 = "% impervious\n(1000 m)",
                      no.plants = "Plant\nrichness",
                      prop.nonnative = "% non-native\nplants",
                      prop.wood = "% woody\nplants",
                      maint_disturbance = "Disturbance\nintensity",
                      maint_time = "Time\ninvestment")

### 3. Build the standardized-coefficient table
# For every combination of network type (Combined/Bipartite/Multitrophic), degree-distribution metric (median/skewness/CV, as z-scores), 
# and site covariate, fit metric ~ covariate (both standardized) and extract the slope, its 95% CI, and its significance. 
# This produces the data behind a forest-plot-style figure of standardized effect sizes.
coef_df <- foodweb_properties %>%
  select(ugs, type, names(covariate_labels), cv_degree_z, skew_degree_z, median_degree_z) %>%
  pivot_longer(cols = ends_with("_z"), names_to = "metric", values_to = "metric_value") %>%
  pivot_longer(cols = all_of(names(covariate_labels)),
               names_to = "covariate",
               values_to = "cov_value") %>%
  mutate(covariate = factor(covariate, levels = rev(names(covariate_labels))),
         metric = factor(metric, levels = c("median_degree_z", "skew_degree_z", "cv_degree_z"))) %>%
  mutate(covariate = recode(covariate, !!!covariate_labels),
         metric = recode(metric,
                         median_degree_z = "Median (z)",
                         skew_degree_z = "Skewness (z)",
                         cv_degree_z = "Coefficient of variation (z)"),
         type = factor(type, levels = c("Combined", "Bipartite", "Multitrophic"))) %>%
  group_by(covariate, type) %>% # Standardize the covariate within each covariate x network-type group
  mutate(cov_scaled = scale(cov_value)) %>%
  group_by(type, metric, covariate) %>% # Fit one standardized linear model per network type x metric x covariate.
  do(tidy(lm(scale(metric_value) ~ cov_scaled, data = .),
          conf.int = TRUE,
          conf.level = 0.95)) %>%
  filter(term == "cov_scaled") %>% # Keep only the covariate's slope term, drop the intercept row.
  mutate(sig = case_when(p.value < 0.001 ~ "***",   # Significance stars and a coarse significant/NS flag (used for point/CI opacity in the plot).
                         p.value < 0.01 ~ "**",
                         p.value < 0.05 ~ "*",
                         p.value < 0.1 ~ "·",
                         TRUE ~ ""),
         significant = ifelse(p.value < 0.1, "Significant", "NS")) %>%
  mutate(star_x = ifelse(estimate > 0, conf.high, conf.low),   # Position each significance star just past the CI end that's furthest from zero, on whichever side the estimate points.
         hjust_star = ifelse(estimate > 0, -0.5, 1.5))

### 4. Figure 3: standardized coefficients by covariate, metric, and network type
ggplot(coef_df, aes(x = estimate, y = covariate, color = type, fill = type, alpha = significant)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3, position = position_dodge(0.6)) +
  geom_errorbar(aes(xmin = conf.low,
                    xmax = conf.high),
                width = 0, position = position_dodge(0.6)) +
  geom_text(aes(label = sig, x = star_x, hjust = hjust_star),
            position = position_dodge(0.6), vjust = 0.7,
            size = 5,
            show.legend = FALSE) +
  scale_alpha_manual(values = c("Significant" = 1, "NS" = 0.2), name = "") +
  scale_color_manual(values = type_colors, name = "") +
  scale_fill_manual(values = type_colors, name = "") +
  facet_wrap(~ metric, scales = "free_x") +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  theme_pubr(base_size = 16) +
  labs(x = "Standardized coefficient and 95% CI", y = NULL) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 16), 
        legend.position = "right")
