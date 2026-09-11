library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

### 1. Load food-web properties
# foodweb_properties: per-site network metrics produced by the food-web metrics script.
# Each metric has an observed value and a corresponding `.rand` null-model value. Their distributions are compared for each focal network type.
foodweb_properties <- read.csv("data/02_processed/fw_properties.csv")

### 2. Reshape observed and null-model metrics for plotting
metric_labels <- c(modularity = "Modularity",
                   nestedness = "Nestedness",
                   r50 = "Robustness",
                   median_degree = "Median (degree)",
                   skew_degree = "Skewness (degree)",
                   cv_degree = "CV (degree)")

foodweb_properties_long <- foodweb_properties %>%
  pivot_longer(cols = all_of(c(names(metric_labels), paste0(names(metric_labels), "_null_mean"))),
               names_to = "metric_name",
               values_to = "value") %>%
  mutate(model_type = if_else(grepl("\\_null_mean$", metric_name), "Null", "Observed"),
    metric = sub("\\_null_mean$", "", metric_name),
    model_type = factor(model_type, levels = c("Observed", "Null")),
    type = factor(type, levels = c("Combined", "Bipartite", "Multitrophic")))

### 3. Figure S9: observed versus null-model metric distributions
ggplot(foodweb_properties_long, aes(x = value, fill = model_type)) +
  geom_density(alpha = 0.5) +
  facet_wrap(metric ~ type,
             scales = "free",
             nrow = 6,
             ncol = 3,
             axis.labels = "margins",
             labeller = labeller(.multi_line = FALSE,
                                 metric = metric_labels)) +
  labs(x = NULL, y = "Density", fill = NULL) +
  theme_pubr(base_size = 16) +
  theme(strip.background = element_blank(),
        strip.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        panel.spacing.x = unit(1, "lines"),
        legend.position = "top")
