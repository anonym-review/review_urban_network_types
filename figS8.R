library(dplyr)
library(ggplot2)

### 1. Load food-web properties and environmental covariates
# foodweb_properties: per-site network metrics produced by the food-web metrics script.
# covariates: per-site covariates.
foodweb_properties <- read.csv("data/02_processed/fw_properties.csv")
covariates <- read.csv("data/01_raw/covariates.csv")

foodweb_properties <- merge(foodweb_properties, covariates)

metric_labels <- c(median_degree_z = "Median (z)",
                   skew_degree_z = "Skewness (z)",
                   cv_degree_z = "Coefficient of variation (z)")

### 2. Calculate correlations among degree-distribution metrics
# For each network type, calculate pairwise Pearson correlations between median, skewness, and coefficient-of-variation degree metrics. 
calculate_metric_correlations <- function(data, metrics) {
  
  metric_pairs <- combn(metrics, 2, simplify = FALSE)
  
  bind_rows(lapply(metric_pairs, function(metric_pair) {
    
    metric_1 <- metric_pair[1]
    metric_2 <- metric_pair[2]
    
    complete_rows <- complete.cases(data[[metric_1]], data[[metric_2]])
    
    correlation_test <- cor.test(data[[metric_1]][complete_rows],
                                 data[[metric_2]][complete_rows],
                                 method = "pearson")
    
    data.frame(metric_1 = metric_1,
               metric_2 = metric_2,
               correlation = unname(correlation_test$estimate),
               p_value = correlation_test$p.value)
  }))
}

correlation_data <- foodweb_properties %>%
  group_by(type) %>%
  group_modify(~ calculate_metric_correlations(.x, names(metric_labels))) %>%
  ungroup() %>%
  mutate(significance = case_when(p_value < 0.001 ~ "***",
                                  p_value < 0.01 ~ "**",
                                  p_value < 0.05 ~ "*",
                                  TRUE ~ ""),
    metric_1 = recode(metric_1, !!!metric_labels),
    metric_2 = recode(metric_2, !!!metric_labels),
    metric_1 = factor(metric_1, levels = unname(metric_labels)),
    metric_2 = factor(metric_2, levels = unname(metric_labels)),
    type = factor(type, levels = c("Bipartite", "Multitrophic", "Combined")),
    correlation_label = if_else(p_value < 0.05, paste0(round(correlation, 2), significance), ""))

# 3. Figure S8: correlations among degree-distribution metrics
ggplot(correlation_data, aes(x = metric_1, y = metric_2, fill = correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = correlation_label), size = 5) +
  scale_fill_gradient2(limits = c(-1, 1),
                       low = "#2166ac",
                       mid = "white",
                       high = "#b2182b") +
  facet_wrap(~ type) +
  labs(fill = "Correlation") +
  coord_fixed() +
  theme_minimal(base_size = 16) +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(size = 16))
