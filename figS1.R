library(dplyr)
library(tidyr)
library(ggplot2)

### 1. Load and standardize environmental covariates
# covariates: per-site covariates.
# Covariates are standardized before correlation analysis so their scales are comparable. 
# Pearson correlations are calculated using complete observations for each individual covariate pair.
covariates <- read.csv("data/01_raw/covariates.csv")

covariate_labels <- c(impervious_1000 = "% impervious (1000 m)",
                      no.plants = "Plant richness",
                      prop.nonnative = "% non-native plants",
                      prop.wood = "% woody species",
                      maint_disturbance = "Disturbance intensity",
                      maint_time = "Time investment")

standardized_covariates <- covariates %>%
  select(all_of(names(covariate_labels))) %>%
  scale() %>%
  as.data.frame()

### 2. Calculate pairwise Pearson correlations and p-values
calculate_correlation_p_values <- function(data) {
  
  variable_names <- colnames(data)
  p_value_matrix <- matrix(NA_real_, nrow = length(variable_names), ncol = length(variable_names), dimnames = list(variable_names, variable_names))
  
  for (row_index in seq_along(variable_names)) {
    for (column_index in seq_along(variable_names)) {
      
      variable_1 <- data[[row_index]]
      variable_2 <- data[[column_index]]
      complete_rows <- complete.cases(variable_1, variable_2)
      
      p_value_matrix[row_index, column_index] <- cor.test(variable_1[complete_rows],
                                                          variable_2[complete_rows],
                                                          method = "pearson")$p.value
    }
  }
  
  p_value_matrix
}

correlation_matrix <- cor(standardized_covariates,
                          use = "pairwise.complete.obs",
                          method = "pearson")

p_value_matrix <- calculate_correlation_p_values(standardized_covariates)

### 3. Prepare the lower-triangle correlation matrix for plotting
axis_levels <- recode(colnames(standardized_covariates), !!!covariate_labels)

correlation_data <- as.data.frame(correlation_matrix) %>%
  tibble::rownames_to_column("covariate_1") %>%
  pivot_longer(-covariate_1,
               names_to = "covariate_2",
               values_to = "correlation") %>%
  left_join(as.data.frame(p_value_matrix) %>%
              tibble::rownames_to_column("covariate_1") %>%
              pivot_longer(-covariate_1,
                           names_to = "covariate_2",
                           values_to = "p_value"),
            by = c("covariate_1", "covariate_2")) %>%
  mutate(axis_x = factor(recode(covariate_1, !!!covariate_labels), levels = axis_levels),
         axis_y = factor(recode(covariate_2, !!!covariate_labels), levels = axis_levels)) %>%
  filter(as.integer(axis_y) >= as.integer(axis_x))

### 4. Figure S1: correlations among environmental covariates
ggplot(correlation_data, aes(x = axis_x, y = axis_y)) +
  geom_vline(aes(xintercept = axis_x),
             color = "grey",
             linewidth = 0.1) +
  geom_hline(aes(yintercept = axis_y),
             color = "grey",
             linewidth = 0.1) +
  geom_point(data = filter(correlation_data, p_value <= 0.05),
             aes(color = correlation),
             size = 9,
             shape = 16) +
  geom_point(data = filter(correlation_data, p_value > 0.05),
             size = 7,
             shape = 16,
             color = NA) +
  geom_text(data = filter(correlation_data, p_value <= 0.05),
            aes(label = sprintf("%.2f", correlation)),
            fontface = "bold",
            size = 3.3,
            color = "black") +
  scale_color_gradient2(low = "#2166ac",
                        mid = "white",
                        high = "#b2182b",
                        limits = c(-1, 1),
                        name = "Correlation") +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 16) +
  theme(axis.text.y = element_text(hjust = 1),
        axis.text.x = element_text(angle = 45, hjust = 1))
