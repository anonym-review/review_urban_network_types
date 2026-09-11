library(piecewiseSEM)
library(dplyr)
library(qgraph)

# 1. Load data and split by network type
foodweb_properties <- read.csv("data/02_processed/fw_properties.csv")
covariates <- read.csv("data/01_raw/covariates.csv")

foodweb_properties <- merge(foodweb_properties, covariates)

combined_df <- foodweb_properties %>% filter(type == "Combined")
bipartite_df <- foodweb_properties %>% filter(type == "Bipartite")
multitrophic_df <- foodweb_properties %>% filter(type == "Multitrophic")

# Shared color scale for the three network types, used in the final figure.
type_colors <- c("Combined" = "black",
                 "Bipartite" = "#2166ac",
                 "Multitrophic" = "#b2182b")

### 2. Fit one structural equation model (SEM) per predictor x network type
# Every predictor is tested against the same 5-equation causal chain:
# predictor -> species richness -> connectance -> modularity, nestedness -> robustness (r50), with the predictor also allowed a direct path into every downstream node.
# Only the predictor variable changes between models, so it's factored out into one helper instead of retyping the 5-model block 18 times.
build_foodweb_sem <- function(data, predictor) {
  psem(do.call(MASS::glm.nb, list(formula = reformulate(predictor, response = "species_richness"), data = data)),
       do.call(lm, list(formula = reformulate(c("species_richness", predictor), response = "connectance"), data = data)),
       do.call(lm, list(formula = reformulate(c("species_richness", "connectance", predictor), response = "modularity"), data = data)),
       do.call(lm, list(formula = reformulate(c("species_richness", "connectance", predictor), response = "nestedness"), data = data)),
       do.call(lm, list(formula = reformulate(c("modularity", "nestedness", "species_richness", "connectance", predictor), response = "r50"), data = data)))
}

# One row per predictor: the column name used in the model formulas, and the label semEff's effect table reports it under
# semEff renders some of these with a "." instead of "_", so the two need to be tracked separately.
predictor_specs <- list(list(id = "impervious", formula_var = "impervious_1000", effect_var = "impervious.1000", label = "% impervious\n(1000 m)"),
                        list(id = "plant_richness", formula_var = "no.plants", effect_var = "no.plants", label = "Plant richness"),
                        list(id = "nonative", formula_var = "prop.nonnative", effect_var = "prop.nonnative", label = "% non-native plants"),
                        list(id = "wood", formula_var = "prop.wood", effect_var = "prop.wood", label = "% woody plants"),
                        list(id = "disturbance", formula_var = "maint_disturbance", effect_var = "maint.disturbance", label = "Disturbance frequency"),
                        list(id = "time", formula_var = "maint_time", effect_var = "maint.time", label = "Time investment"))

network_data <- list(Combined = combined_df,
                     Bipartite = bipartite_df,
                     Multitrophic = multitrophic_df)

# Fit all 18 predictor x network-type SEMs, keyed as e.g. "Bipartite_wood".
sem_models <- list()
for (spec in predictor_specs) {
  for (type_label in names(network_data)) {
    model_id <- paste(type_label, spec$id, sep = "_")
    sem_models[[model_id]] <- build_foodweb_sem(network_data[[type_label]], spec$formula_var)
  }
}

# Inspect each model's fit (Fisher's C, separation tests, path coefficients, R^2, etc.).
for (model_id in names(sem_models)) {
  cat("\n====", model_id, "====\n")
  print(summary(sem_models[[model_id]]))
}


### 3. Create path diagrams for fitted SEMs
# Edge width represents the absolute standardized coefficient. Positive effects are blue, negative effects are red, and non-significant effects are grey.
# Solid edges indicate p < 0.05; dashed edges indicate 0.05 <= p < 0.10.
create_sem_plot <- function(sem_model,
                            predictor_label,
                            response_labels,
                            node_size = c(43, 23),
                            edge_label_size = 2.2) {
  
  coefficient_data <- summary(sem_model, .progressBar = FALSE)$coefficients
  
  names(coefficient_data)[names(coefficient_data) == ""] <- "significance"
  
  coefficient_data <- coefficient_data %>%
    mutate(edge_style = case_when(P.Value < 0.05 ~ "solid",
                                  P.Value < 0.1 ~ "dashed",
                                  TRUE ~ "solid"),
           edge_color = case_when(P.Value >= 0.1 ~ "grey",
                                  Std.Estimate < 0 ~ "#d7301f",
                                  TRUE ~ "#0570b0"),
           standardized_estimate = if_else(P.Value >= 0.1, NA_real_, round(Std.Estimate, 3)),
           edge_background = case_when(edge_color == "#d7301f" ~ "#fddbc7",
                                       edge_color == "#0570b0" ~ "#d1e5f0",
                                       TRUE ~ NA_character_))
  
  sem_graph <- igraph::graph_from_adjacency_matrix(getDAG(sem_model))
  edge_list <- igraph::as_edgelist(sem_graph)
  
  # Match coefficient rows to the directed edges returned by the SEM DAG.
  coefficient_data <- coefficient_data[match(paste(edge_list[, 1], edge_list[, 2]),
                                             paste(coefficient_data$Predictor, coefficient_data$Response)), ]
  
  node_x <- c(1.5, 0, 3, 0, 3, 1.5)
  node_y <- c(2.2, 1.6, 1.6, 0.6, 0.6, 0.1)
  
  node_layout <- matrix(c(node_x, node_y), ncol = 2)
  
  sem_plot <- qgraph(edge_list,
                     edge.color = coefficient_data$edge_color,
                     edge.labels = coefficient_data$standardized_estimate,
                     edge.label.position = 0.25,
                     edge.label.cex = edge_label_size,
                     edge.label.margin = 0.007,
                     edge.label.bg = coefficient_data$edge_background,
                     edge.label.color = coefficient_data$edge_color,
                     lty = coefficient_data$edge_style,
                     layout = node_layout,
                     shape = "ellipse",
                     label.cex = 1.2,
                     label.scale = FALSE,
                     borders = FALSE,
                     vsize = node_size[1],
                     vsize2 = node_size[2],
                     edge.width = if_else(is.na(coefficient_data$standardized_estimate), 1,abs(coefficient_data$standardized_estimate) * 10),
                     mar = c(2, 3, 2, 4.5),
                     DoNotPlot = TRUE)
  
  # The predictor occupies the first node; all remaining nodes display their
  # corresponding R-squared values from the piecewise SEM summary.
  node_labels <- bind_rows(data.frame(label = predictor_label, r_squared = ""),
                           data.frame(label = response_labels,
                                      r_squared = paste0("R2: ", summary(sem_model)$R2$R.squared)))
  
  sem_plot$graphAttributes$Nodes$labels <- paste(node_labels$label,
                                                 node_labels$r_squared,
                                                 sep = "\n")
  
  sem_plot
}

response_labels <- c(
  "Consumer\nrichness",
  "Connectance",
  "Modularity",
  "Nestedness",
  "Robustness"
)

### 4. Plot the SEM diagrams
# One diagram is produced for every predictor x network-type combination.

network_plot_order <- c("Bipartite", "Multitrophic", "Combined")

par(mfrow = c(6, 3))

for (predictor_spec in predictor_specs) {
  for (network_type in network_plot_order) {
    
    model_id <- paste(network_type, predictor_spec$id, sep = "_")
    
    plot(create_sem_plot(sem_model = sem_models[[model_id]],
                         predictor_label = predictor_spec$label,
                         response_labels = response_labels))
  }
}
