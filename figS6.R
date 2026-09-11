library(dplyr)
library(ggplot2)
library(ggpubr)

### 1. Load inferred and observed food-web properties
# inferred_foodweb_properties: per-site network metrics produced by the food-web metrics script.
# observed_foodweb_properties: per-site network metrics produced by the food-web metrics script.
# The inferred bipartite food webs are compared with independently observed food-web properties to assess whether structural metrics are consistently estimated across sites.
inferred_foodweb_properties <- read.csv("data/02_processed/fw_properties.csv") %>% filter(type == "Bipartite")
observed_foodweb_properties <- read.csv("data/02_processed/observed_fw_properties.csv")

### 2. Match sites and select food-web metrics for validation
# Joining by `ugs` ensures inferred and observed values from the same study site are paired before their correlations are calculated.
metric_labels <- c(connectance = "Connectance",
                   nestedness = "Nestedness",
                   modularity = "Modularity",
                   r50 = "Robustness",
                   median_degree_z = "Median degree",
                   skew_degree_z = "Skewness degree",
                   cv_degree_z = "CV degree")

validation_data <- inferred_foodweb_properties %>%
  select(ugs, all_of(names(metric_labels))) %>%
  inner_join(observed_foodweb_properties %>% select(ugs, all_of(names(metric_labels))),
             by = "ugs", suffix = c("_inferred", "_observed"))

### 3. Calculate correlations between inferred and observed metrics
correlation_summary <- lapply(names(metric_labels), function(metric) {
  
  inferred_values <- validation_data[[paste0(metric, "_inferred")]]
  observed_values <- validation_data[[paste0(metric, "_observed")]]

  correlation_test <- cor.test(inferred_values, observed_values, use = "complete.obs")
  
  data.frame(metric = metric_labels[[metric]],
             correlation = unname(correlation_test$estimate),
             p_value = correlation_test$p.value)
  }) %>%
  bind_rows() %>%
  mutate(metric = factor(metric, levels = c("Median degree", "Skewness degree", "CV degree", "Connectance", "Modularity", "Nestedness", "Robustness")))

### 4. Figure S6: correlations between inferred and observed metrics
ggplot(correlation_summary, aes(x = metric, y = correlation)) +
  geom_point(size = 12, shape = 1) +
  geom_text(aes(label = round(correlation, 2)), size = 4) +
  labs(x = NULL, y = "Correlation") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_pubclean(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16))
