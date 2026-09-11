library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(ggpubr)
library(broom)

### 1. Load species occurrences, covariates, and the focal metaweb
# fauna_matrix: species per site matrix for the observed animal taxa.
# flora_matrix: species per site matrix for the observed plant taxa.
# metaweb: the "metaweb", i.e. all known/reported trophic (and other) interactions, with each Source/Target given at a taxonomic rank (Species, Genus, Family, Order).
# covariates: per-site covariates.
fauna_matrix <- read.csv("data/01_raw/fauna_matrix.csv", check.names = FALSE)
flora_matrix <- read.csv("data/01_raw/flora_matrix.csv", check.names = FALSE)
species_matrix <- merge(fauna_matrix, flora_matrix, by = "ugs")
metaweb <- read.csv("data/02_processed/focal_metaweb.csv")
covariates <- read.csv("data/01_raw/covariates.csv")

# Plant taxa and detritus form the basal resource groups. Exclude the `ugs`
# site-identifier column, which is not itself a taxon.
resources <- c(setdiff(colnames(flora_matrix), "ugs"), "Detritus")
species_matrix$Detritus <- 1

fauna_taxa <- setdiff(colnames(fauna_matrix), "ugs")

### 2. Assign animal taxa to trophic guilds
# Primary consumers feed only on basal resources, omnivores feed on both basal resources and animal taxa, and predators feed on animal taxa but not basal resources. 
# Pieris is excluded from animal prey taxa as in the original model.
resource_consumers <- unique(metaweb$Source_Name[metaweb$Target_Name %in% resources])
animal_consumers <- unique(metaweb$Source_Name[metaweb$Target_Name %in% setdiff(fauna_taxa, "Pieris")])

omnivores <- intersect(resource_consumers, animal_consumers)
primary_consumers <- setdiff(resource_consumers, animal_consumers)
predators <- setdiff(unique(metaweb$Source_Name[metaweb$Target_Name %in% fauna_taxa]),
                     resource_consumers)

guild_colors <- c("Primary consumers" = "#E69F00",
                  "Omnivores" = "#0471a6",
                  "Predators" = "#df2935")

covariate_labels <- c(impervious_1000 = "% impervious\n(1000 m)", 
                      no.plants = "Plant\nrichness",
                      prop.nonnative = "% non-native\nplants",
                      prop.wood = "% woody\nplants",
                      maint_disturbance = "Disturbance\nintensity",
                      maint_time = "Time\ninvestment")

# 3. Calculate trophic-guild proportions at each site
# Guild richness is divided by consumer richness (all non-basal taxa) to obtain the proportion of consumer nodes belonging to each trophic guild.
primary_consumer_richness <- specnumber(species_matrix %>% select(any_of(primary_consumers)))
omnivore_richness <- specnumber(species_matrix %>% select(any_of(omnivores)))
predator_richness <- specnumber(species_matrix %>% select(any_of(predators)))
basal_richness <- specnumber(species_matrix %>% select(any_of(resources)))
total_richness <- specnumber(species_matrix %>% select(any_of(c(primary_consumers, omnivores, predators, resources))))
consumer_richness <- total_richness - basal_richness

guild_composition <- bind_rows(data.frame(richness = primary_consumer_richness,
                                          fraction = primary_consumer_richness / consumer_richness,
                                          ugs = species_matrix$ugs,
                                          guild = "Primary consumers"),
                               data.frame(richness = omnivore_richness,
                                          fraction = omnivore_richness / consumer_richness,
                                          ugs = species_matrix$ugs,
                                          guild = "Omnivores"),
                               data.frame(richness = predator_richness,
                                          fraction = predator_richness / consumer_richness,
                                          ugs = species_matrix$ugs,
                                          guild = "Predators"))

### 4. Reshape guild proportions and test covariate relationships
guild_composition_long <- guild_composition %>%
  left_join(covariates, by = "ugs") %>%
  select(ugs, guild, all_of(names(covariate_labels)), fraction) %>%
  pivot_longer(cols = all_of(names(covariate_labels)),
               names_to = "covariate",
               values_to = "covariate_value") %>%
  mutate(covariate = factor(covariate, levels = names(covariate_labels)),
         covariate = recode(covariate, !!!covariate_labels))

# Fit one linear model for every guild x covariate combination. 
# The p-value determines whether its fitted line is displayed as solid or dashed.
significance_data <- guild_composition_long %>%
  group_by(covariate, guild) %>%
  do(tidy(lm(fraction ~ covariate_value, data = .))) %>%
  filter(term == "covariate_value") %>%
  transmute(covariate, 
            guild,
            p_value = p.value,
            significance = if_else(p_value < 0.05, "Significant", "NS"))

guild_composition_long <- guild_composition_long %>%
  left_join(significance_data, by = c("guild", "covariate"))

### 5. Figure S11: trophic-guild composition along environmental gradients
ggplot(guild_composition_long, aes(x = covariate_value, y = fraction, color = guild, linetype = significance)) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.1) +
  facet_wrap(~ covariate, scales = "free", nrow = 1, ncol = 6) +
  scale_linetype_manual(values = c("Significant" = "solid", "NS" = "dashed")) +
  scale_color_manual(values = guild_colors) +
  labs(x = NULL,
       y = "% of consumer nodes",
       color = NULL,
       linetype = NULL) +
  theme_pubr(base_size = 16) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 15),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14))
