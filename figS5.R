library(dplyr)
library(tidyverse)
library(stringr)
library(eulerr)
library(ggforce)

### 0. Load raw data
# metaweb       : the "metaweb", i.e. all known/reported trophic (and other) interactions, with each Source/Target given at a taxonomic rank (Species, Genus, Family, Order).
# fauna_taxonomy  : taxonomy table (genus/family/order) for the observed animal taxa.
# flora_taxonomy  : taxonomy table (genus/family/order) for the observed plant taxa.
metaweb <- read.csv("data/01_raw/metaweb.csv")
fauna_taxonomy <- read.csv("data/01_raw/fauna_taxo.csv")
flora_taxonomy <- read.csv("data/01_raw/flora_taxo.csv")

# Manually add two order-level detritivory links that are not in the source metaweb:
# Collembola and Acari (mites) both feed on Detritus.
metaweb <- rbind(metaweb,
                 data.frame(Source_Name = c("Collembola", "Acari"),
                            Source_Rank = rep("Order", 2),
                            Target_Name = rep("Detritus", 2),
                            Target_Rank = rep("Species", 2),
                            Interaction_Type = c("Detritivore")))

# For parasitism records, the metaweb lists Source = host, Target = parasite.
# Flip Source/Target for those rows so that Source is always the "actor" (the parasite) and 
# Target is always the "recipient" (the host) -- consistent with predation direction.
# Then drop the rows that involve three problematic taxa where the source data is unreliable or ambiguous at species level (hybrids / misidentifications): 
# Prunella (accentors), Linaria (linnets/twite), and Arenaria interpres (turnstone).
metaweb <- metaweb %>%
  mutate(old_source = Source_Name,
         old_target = Target_Name,
         Source_Name = if_else(Interaction_Type == "Parasitism", old_target, old_source, missing = old_source),
         Target_Name = if_else(Interaction_Type == "Parasitism", old_source, old_target, missing = old_target)) %>%
  select(-old_source, -old_target) %>%
  filter(!(Source_Name %in% c("Prunella collaris", "Prunella modularis") | Target_Name %in% c("Prunella collaris", "Prunella modularis") | (Target_Name == "Prunella" & Interaction_Type == "Predation") |
             Source_Name %in% c("Linaria cannabina", "Linaria flavirostris") | Target_Name %in% c("Linaria cannabina", "Linaria flavirostris") | (Target_Name == "Linaria" & Interaction_Type == "Predation") |
             Source_Name == "Arenaria interpres" | Target_Name == "Arenaria interpres" | (Target_Name == "Arenaria" & Interaction_Type == "Predation")))

# Reduce any species-level names to just their genus (first word), so everything that isn't already Family/Order rank ends up at Genus rank. Interaction_Type is no longer needed after the parasite/host flip above, so it's dropped.
metaweb <- metaweb %>%
  mutate(Source_Name = if_else(Source_Rank == "Species", word(Source_Name, 1), Source_Name),
         Target_Name = if_else(Target_Rank == "Species", word(Target_Name, 1), Target_Name)) %>%
  select(-Interaction_Type)

### 1. Prepare taxonomy & observed genera
# Combine fauna and flora taxonomies into one lookup table, and get the list of genera that were 
# actually observed in the study (this is the target set we want the final metaweb to be restricted to).
all_tax <- bind_rows(fauna_taxonomy, flora_taxonomy)
observed_genera <- unique(all_tax$genus)

### 2. Identify observed genera that are missing from the metaweb
# For each observed genus, check whether it already appears in the metaweb at genus level (as a Source for fauna, as a Target for flora). 
# Genera that don't appear will need their interactions inferred from a coarser taxonomic level (family or order) below.
metaweb_unknown_genera <- c(setdiff(fauna_taxonomy$genus, metaweb$Source_Name),
                            setdiff(flora_taxonomy$genus, metaweb$Target_Name))
metaweb_unknown_families <- all_tax$family[all_tax$genus %in% metaweb_unknown_genera]

### 3. Family-level expansion
# Keep only the family-level metaweb rows that involve a family containing at least one "missing" genus from step 2 
# These are the rows we need to expand down to genus level.
metaweb_family <- metaweb %>%
  filter(Source_Rank == "Family" | Target_Rank == "Family") %>%
  filter(Source_Name %in% metaweb_unknown_families | Target_Name %in% metaweb_unknown_families)

# If the Source side is a family, replace it with every genus in that family.
family_source_exp <- metaweb_family %>%
  filter(Source_Rank == "Family") %>%
  left_join(all_tax, by = c("Source_Name" = "family")) %>%
  mutate(Source_Name = genus, Source_Rank = "Genus") %>%
  select(-genus, -order)

# Same expansion, but for the Target side.
family_target_exp <- metaweb_family %>%
  filter(Target_Rank == "Family") %>%
  left_join(all_tax, by = c("Target_Name" = "family")) %>%
  mutate(Target_Name = genus, Target_Rank = "Genus") %>%
  select(-genus, -order)

# Combine both directions of family expansion into one table.
family_expanded <- bind_rows(family_source_exp, family_target_exp) %>% distinct()

# Safety net: if any Species-rank names slipped through the join, reduce them to genus too.
family_expanded <- family_expanded %>%
  mutate(Source_Name = if_else(Source_Rank == "Species", word(Source_Name, 1), Source_Name, missing = Source_Name),
         Target_Name = if_else(Target_Rank == "Species", word(Target_Name, 1), Target_Name, missing = Target_Name))


### 4. Order-level expansion (Collembola & Acari only)
# These two taxa are only known at Order rank in the metaweb (see step 0), so their interactions need to 
# be expanded down to genus level using the taxonomy table.
metaweb_order <- metaweb %>%
  filter(Source_Name %in% c("Collembola", "Acari") | Target_Name %in% c("Collembola", "Acari"))

# Expand the Source order into all matching genera, then -- for the Target side -- also expand any Family-rank target into its constituent genera in the same step.
order_exp <- metaweb_order %>%
  left_join(all_tax, by = c("Source_Name" = "order")) %>%
  mutate(Source_Name = genus, Source_Rank = "Genus") %>%
  select(-genus, -family) %>%
  left_join(all_tax, by = c("Target_Name" = "family")) %>%
  mutate(Target_Name = if_else(Target_Rank == "Family", genus, Target_Name, missing = Target_Name),
         Target_Rank = if_else(Target_Rank == "Family", "Genus", Target_Rank, missing = Target_Rank)) %>%
  select(colnames(metaweb))

### 5. Combine the original metaweb with both expansions
metaweb_augmented <- bind_rows(metaweb,
                               family_expanded,
                               order_exp) %>%
  distinct()

### 6. Build the final genus-level metaweb, restricted to observed taxa
# Force everything to Genus rank (reducing any remaining Species names to genus first), keep only links where both ends 
# are either an observed genus or "Detritus", drop self-loops (a genus interacting with itself), and keep only unique Source-Target pairs.
metaweb_cleaned <- metaweb_augmented %>%
  mutate(Source_Name = if_else(Source_Rank == "Species", word(Source_Name, 1), Source_Name, missing = Source_Name),
         Target_Name = if_else(Target_Rank == "Species", word(Target_Name, 1), Target_Name, missing = Target_Name),
         Source_Rank = "Genus",
         Target_Rank = "Genus") %>%
  filter(Source_Name %in% observed_genera | Source_Name == "Detritus",
         Target_Name %in% observed_genera | Target_Name == "Detritus",
         Source_Name != Target_Name) %>%
  select(Source_Name, Target_Name) %>%
  distinct()

### 7. Compare predicted and observed pollinator interactions
pollinator_interactions <- read.csv("data/01_raw/pollinator_interactions.csv") %>%
  distinct()

# Observed interactions that are present in the focal metaweb.
realized_interactions <- pollinator_interactions %>%
  semi_join(metaweb_cleaned, by = c("Source_Name", "Target_Name")) %>%
  distinct()

# Observed interactions not predicted by the focal metaweb.
unexpected_interactions <- pollinator_interactions %>%
  anti_join(metaweb_cleaned, by = c("Source_Name", "Target_Name")) %>%
  distinct()

# Candidate focal-metaweb interactions among the observed pollinator taxa.
expected_interactions <- metaweb_cleaned %>%
  filter(Source_Name %in% pollinator_interactions$Source_Name,
         Target_Name %in% pollinator_interactions$Target_Name) %>%
  distinct()

# Predicted interactions that were not observed.
unobserved_interactions <- expected_interactions %>%
  anti_join(pollinator_interactions,
            by = c("Source_Name", "Target_Name")) %>%
  distinct()

### 8. Figure S5: expected versus observed pollinator interactions
euler_fit <- euler(c(Expected = nrow(unobserved_interactions),
                     Observed = nrow(unexpected_interactions),
                     "Expected&Observed" = nrow(realized_interactions)))

circle_data <- as.data.frame(euler_fit$ellipses)
circle_data$set <- rownames(circle_data)

region_labels <- tibble(x = c(circle_data$h[1] - circle_data$a[1] / 2,
                              circle_data$h[2] + circle_data$a[2] * 1.5,
                              mean(circle_data$h) + 10),
                        y = c(circle_data$k[1],
                              circle_data$k[2],
                              mean(circle_data$k)),
                        label = c(paste0("Unobserved\n\n",
                                         nrow(unobserved_interactions), "\n",
                                         round(nrow(unobserved_interactions) / nrow(expected_interactions) * 100, 2), "%"),
                                  paste0("Unexpected\n\n", 
                                         nrow(unexpected_interactions), "\n",
                                         round(nrow(unexpected_interactions) / nrow(pollinator_interactions) * 100, 2), "%"),
                                  paste0("Observed\n\n", nrow(realized_interactions))))

ggplot() +
  geom_circle(data = circle_data, 
              aes(x0 = h, y0 = k, r = a, fill = set),
              alpha = 0.4,
              color = "black") +
  geom_text(data = region_labels,
            aes(x = x, y = y, label = label),
            size = 5) +
  coord_fixed(clip = "off") +
  labs(fill = "") +
  theme_void() + 
  theme(legend.position = "top", 
        margins = margin(0, 100, 0, 0))
