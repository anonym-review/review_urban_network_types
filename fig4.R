library(semEff)
library(piecewiseSEM)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggpattern)

### 1. Load data and split by network type
# foodweb_properties: per-site network metrics produced by the food-web metrics script.
# covariates: per-site covariates.
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
predictor_specs <- list(list(id = "impervious", formula_var = "impervious_1000",   effect_var = "impervious.1000"),
                        list(id = "herb",       formula_var = "no.plants",         effect_var = "no.plants"),
                        list(id = "nonnative",  formula_var = "prop.nonnative",    effect_var = "prop.nonnative"),
                        list(id = "wood",       formula_var = "prop.wood",         effect_var = "prop.wood"),
                        list(id = "disturb",    formula_var = "maint_disturbance", effect_var = "maint.disturbance"),
                        list(id = "time",       formula_var = "maint_time",        effect_var = "maint.time"))

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

### 3. Bootstrap standardized direct/indirect/total effects
R <- 1000
seed <- 123

boot_models <- list()
for (model_id in names(sem_models)) {
  boot_models[[model_id]] <- bootEff(sem_models[[model_id]],
                                     seed = seed, R = R,
                                     parallel = "multicore", type = "nonparametric")
}

# Pull the effect table for one predictor at the 95% CI level, then attach the 90/99/99.9% CI bounds as extra columns (used for the significance stars below).
get_sem_table <- function(boot_obj, predictor, type_label) {
  
  eff90 <- semEff(boot_obj, predictors = predictor, ci.type = "norm", ci.conf = 0.90)$Effects$Table
  eff95 <- semEff(boot_obj, predictors = predictor, ci.type = "norm", ci.conf = 0.95)$Effects$Table
  eff99 <- semEff(boot_obj, predictors = predictor, ci.type = "norm", ci.conf = 0.99)$Effects$Table
  eff999 <- semEff(boot_obj, predictors = predictor, ci.type = "norm", ci.conf = 0.999)$Effects$Table
  
  eff95$lower90 <- eff90$lower_ci
  eff95$upper90 <- eff90$upper_ci
  eff95$lower99 <- eff99$lower_ci
  eff95$upper99 <- eff99$upper_ci
  eff95$lower999 <- eff999$lower_ci
  eff95$upper999 <- eff999$upper_ci
  eff95$type <- type_label
  
  return(eff95)
}

effect_tables <- list()
for (spec in predictor_specs) {
  for (type_label in names(network_data)) {
    model_id <- paste(type_label, spec$id, sep = "_")
    effect_tables[[model_id]] <- get_sem_table(boot_models[[model_id]], spec$effect_var, type_label)
  }
}

# 4. Combine all effect tables and prepare labels for plotting
sem_eff <- bind_rows(effect_tables) %>%
  filter(!effect_type == "mediators") %>% 
  mutate(response = recode(response,
                           species.richness = "Species richness",
                           connectance = "Connectance",
                           modularity = "Modularity",
                           nestedness = "Nestedness",
                           r50 = "Robustness")) %>%
  mutate(response = factor(response, levels = c("Robustness", "Modularity", "Nestedness", "Connectance", "Species richness"))) %>%
  mutate(predictor = recode(predictor,
                            impervious.1000 = "% impervious\n(1000 m)",
                            no.plants = "Plant\nrichness",
                            prop.nonnative = "% non-native\nplants",
                            prop.wood = "% woody\nplants",
                            maint.disturbance = "Disturbance\nintensity",
                            maint.time = "Time\ninvestment")) %>%
  mutate(predictor = factor(predictor, levels = c("% impervious\n(1000 m)",
                                                  "Plant\nrichness",
                                                  "% non-native\nplants",
                                                  "% woody\nplants",
                                                  "Disturbance\nintensity",
                                                  "Time\ninvestment"))) %>%
  mutate(stars = case_when(lower999 > 0 | upper999 < 0 ~ "***", # Significance stars from the tightest CI level that still excludes zero, and a coarse significant/non-significant flag for point/CI opacity.
                           lower99 > 0 | upper99 < 0 ~ "**",
                           lower_ci > 0 | upper_ci < 0 ~ "*",
                           lower90 > 0 | upper90 < 0 ~ "·",
                           TRUE ~ ""),
         signif_cat = case_when(stars != "" ~ "significant",
                                TRUE ~ "ns")) %>%
  mutate(ci_width = upper_ci - lower_ci,
         star_x = ifelse(effect > 0, upper_ci + 0.1, lower_ci - 0.1))


# 5. Figure 4: direct, indirect, and total standardized effects
ggplot(sem_eff, aes(y = response, x = effect, fill = type, color = type)) +
  geom_col_pattern(data = sem_eff %>% filter(effect_type %in% c("direct", "indirect")),
                   aes(pattern = effect_type, color = type, pattern_colour = type),
                   position = "stack",
                   alpha = 0.2,
                   pattern_size = 0.35,
                   pattern_angle = 45,
                   pattern_density = 0.01,
                   pattern_spacing = 0.05,
                   pattern_key_scale_factor = .1) +
  geom_point(data = sem_eff %>% filter(effect_type %in% c("total")),
             aes(x = effect, alpha = signif_cat),
             size = 3, position = position_dodge(width = 0.9)) +
  geom_errorbar(data = sem_eff %>% filter(effect_type %in% c("total")),
                aes(xmin = lower_ci, xmax = upper_ci, alpha = signif_cat),
                width = 0, linewidth = 1.2, position = position_dodge(width = 0.9)) +
  geom_text(data = sem_eff %>% filter(effect_type == "total"),
            aes(x = star_x, y = response, label = stars),
            vjust = 0.75,
            position = position_dodge(width = 0.9),
            size = 6,
            show.legend = FALSE) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey40") +
  facet_grid(predictor ~ factor(type, levels = c("Bipartite", "Multitrophic", "Combined")), axes = "all_x", axis.labels = "margins", scales = "free_x") +
  scale_pattern_manual(values = c(direct = "none", indirect = "stripe")) +
  scale_pattern_color_manual(values = type_colors, guide = "none") +
  scale_fill_manual(values = type_colors, guide = "none") +
  scale_color_manual(values = type_colors, guide = "none") +
  scale_alpha_manual(values = c(significant = 1, marginal = 1, ns = 0.3), guide = "none") +
  labs(x = "Standardized coefficients and 95% CI", y = "", alpha = "", color = "", fill = "", pattern = "Effect type") +
  theme_pubr(base_size = 15) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(5, "mm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        plot.margin = margin(0, 0, 0, 0),
        legend.position = "top") +
  guides(pattern = guide_legend(override.aes = list(fill = "white",
                                                    colour = "black",
                                                    pattern_fill = "black",
                                                    pattern_colour = "black",
                                                    pattern_density = 0.01,
                                                    pattern_spacing = 0.08)))
