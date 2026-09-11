library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(legendry)

### 1. Load fauna occurrences, taxonomy, and site covariates
# fauna_matrix: species per site matrix for the observed animal taxa.
# fauna_taxonomy: taxonomy table (genus/family/order) for the observed animal taxa.
# covariates: per-site covariates.
fauna_matrix <- read.csv("data/01_raw/fauna_matrix.csv", check.names = FALSE)
fauna_taxonomy <- read.csv("data/01_raw/fauna_taxo.csv", check.names = FALSE)
covariates <- read.csv("data/01_raw/covariates.csv")

# Define the taxonomic groups displayed in the community-composition figure.
# Some groups include taxa reported at family rank in the taxonomy table.
taxonomic_groups <- list(Hemiptera = c("Hemiptera", "Heteroptera", "Cicadellidae", "Aphididae", "Miridae", "Aphrophoridae"),
                         Hymenoptera = c("Hymenoptera", "Eulophidae"),
                         Diptera = c("Diptera", "Chironomidae"),
                         Coleoptera = c("Coleoptera", "Carabidae", "Coccinellidae"),
                         Arachnida = c("Araneae", "Opiliones", "Arachnida", "Cheiracanthiidae"),
                         Gastropoda = c("Gastropoda", "Stylommatophora"),
                         Collembola = "Collembola",
                         Orthoptera = "Orthoptera",
                         Acari = "Acari",
                         "Other taxa" = c("Thysanoptera", "Psocoptera", "Neuroptera", "Thripidae", "Isopoda", "Dermaptera", "Blattodea", "Insecta",
                                          "Lepidoptera", "Ephemeroptera", "Plecoptera", "Trichoptera", "Polyxenida", "Lithobiomorpha", "Julida"))

### 2. Calculate taxonomic composition for every site
# For each site, calculate the fraction of observed fauna genera belonging to each focal taxonomic group. Fractions are relative to total fauna richness.
calculate_community_composition <- function(species_matrix,
                                            taxonomy,
                                            covariates,
                                            taxonomic_groups) {
  
  composition_results <- vector("list", nrow(species_matrix))
  
  for (site_index in seq_len(nrow(species_matrix))) {
    
    site_id <- species_matrix$ugs[site_index]
    message("Processing site: ", site_id)
    
    present_taxa <- species_matrix[site_index, ] %>%
      select(-ugs) %>%
      select(where(~ .x >= 1)) %>%
      colnames()
    
    species_richness <- length(present_taxa)
    
    group_counts <- sapply(taxonomic_groups, function(group_orders) {
      group_genera <- taxonomy$genus[taxonomy$order %in% group_orders]
      sum(present_taxa %in% group_genera)
    })
    
    composition_results[[site_index]] <- data.frame(ugs = site_id,
                                                    setNames(as.list(group_counts / species_richness),
                                                             paste0("fraction_", make.names(names(group_counts), unique = TRUE))))
  }
  
  composition_data <- bind_rows(composition_results) %>%
    merge(covariates, by = "ugs") 
  
  composition_data %>%
    pivot_longer(cols = starts_with("fraction_"),
                 names_to = "taxonomic_group",
                 values_to = "proportion") %>%
    
    arrange(ugs_type) %>%
    mutate(ugs = stringr::str_split_fixed(ugs, "_", 3)[, 2], 
           ugs_type = interaction(ugs, ugs_type, sep = "_"), 
           taxonomic_group = recode(taxonomic_group,
                                    fraction_Hemiptera = "Hemiptera",
                                    fraction_Hymenoptera = "Hymenoptera",
                                    fraction_Diptera = "Diptera",
                                    fraction_Coleoptera = "Coleoptera",
                                    fraction_Arachnida = "Arachnida",
                                    fraction_Gastropoda = "Gastropoda",
                                    fraction_Collembola = "Collembola",
                                    fraction_Orthoptera = "Orthoptera",
                                    fraction_Acari = "Acari",
                                    fraction_Other.taxa = "Other taxa")) %>%
    arrange(ugs_type)
}

community_composition <- calculate_community_composition(
  species_matrix = fauna_matrix,
  taxonomy = fauna_taxonomy,
  covariates = covariates,
  taxonomic_groups = taxonomic_groups
)

### 3. Figure S10: taxonomic composition of fauna communities
ggplot(community_composition, aes(x = ugs_type, y = proportion, fill = taxonomic_group)) + 
  geom_bar(stat = "identity") + 
  labs(x = "", y = "% of community") + 
  scale_x_discrete(guide = guide_axis_nested(angle = 45, key = key_range_auto(sep = "_"))) + 
  scale_fill_brewer(palette = "Set3", name = "") + 
  theme_pubr(base_size = 16) 
