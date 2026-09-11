library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(patchwork)

### 1. Helper functions
# Build a square, binary (presence/absence) taxon x taxon adjacency matrix from a long-format metaweb (one row per Source -> Target interaction). 
create_adjacency_matrix <- function(metaweb) {
  # Rows are resources/prey and columns are consumers/predators.
  taxa <- sort(unique(c(metaweb$Target_Name, metaweb$Source_Name)))
  
  adjacency_matrix <- matrix(0, nrow = length(taxa), ncol = length(taxa), dimnames = list(taxa, taxa))
  
  interaction_counts <- table(metaweb$Target_Name, metaweb$Source_Name)
  
  adjacency_matrix[rownames(interaction_counts),
                   colnames(interaction_counts)] <- interaction_counts
  
  # Retain only interaction presence/absence, not counts.
  (adjacency_matrix > 0) * 1
}

# For each site, work out what proportion of its observed links fall into each of the 7 possible directed trophic-guild pairings 
# (basal -> primary consumer, basal -> omnivore, primary consumer -> omnivore, primary consumer -> predator, omnivore -> omnivore, omnivore -> predator, predator -> predator), 
# then compare those proportions to a null model that distributes the same number of links uniformly at random among all biologically allowed guild pairings.
# Guild abbreviations used throughout: B = basal, PC = primary consumer, OM = omnivore, PR = predator.
calculate_link_probabilities <- function(metaweb_matrix,
                                         species_matrix,
                                         primary_consumers,
                                         omnivores,
                                         predators,
                                         basal_taxa = NULL,
                                         remove_isolated_nodes = TRUE,
                                         n_randomisations = 100) {
  
  results <- vector("list", nrow(species_matrix))
  
  for (site_index in seq_len(nrow(species_matrix))) {
    
    site_id <- species_matrix$ugs[site_index]
    message("Processing site: ", site_id)
    
    ### Taxa present at this site that are also in the metaweb.
    present_taxa <- species_matrix[site_index, ] %>%
      select(-ugs) %>%
      select(any_of(intersect(colnames(metaweb_matrix), rownames(metaweb_matrix)))) %>%
      select(where(~ .x >= 1)) %>%
      colnames()
    
    local_matrix <- metaweb_matrix[present_taxa, present_taxa, drop = FALSE]
    
    ## Drop taxa with no links at all at this site.
    isolated_nodes <- which(rowSums(local_matrix) == 0 & colSums(local_matrix) == 0)
    if (remove_isolated_nodes && length(isolated_nodes) > 0) {
      local_matrix <- local_matrix[-isolated_nodes, -isolated_nodes, drop = FALSE]
    }
    
    species_richness <- ncol(local_matrix)
    taxa_names <- colnames(local_matrix)
    
    ## Assign each taxon to a trophic guild.
    # Basal taxa: either supplied explicitly via `basal_taxa`, or (by default) inferred as taxa that are never a link source, i.e. colSums == 0.
    basal_idx <- if (!is.null(basal_taxa)) which(taxa_names %in% basal_taxa) else which(colSums(local_matrix) == 0)
    pcons_idx <- which(taxa_names %in% primary_consumers)
    omni_idx  <- which(taxa_names %in% omnivores)
    pred_idx  <- which(taxa_names %in% predators)
    
    total_links <- sum(local_matrix)
    
    ### Allowed dyays per interaction type.
    # Enumerate every individual (from, to) taxon pair that is biologically allowed for each of the 7 guild-pairing types. 
    dyads_B_PC  <- expand.grid(from = basal_idx,  to = pcons_idx)
    dyads_B_OM  <- expand.grid(from = basal_idx,  to = omni_idx)
    dyads_PC_OM <- expand.grid(from = pcons_idx,  to = omni_idx)
    dyads_PC_PR <- expand.grid(from = pcons_idx,  to = pred_idx)
    dyads_OM_OM <- expand.grid(from = omni_idx,   to = omni_idx)
    dyads_OM_PR <- expand.grid(from = omni_idx,   to = pred_idx)
    dyads_PR_PR <- expand.grid(from = pred_idx,   to = pred_idx)
    
    allowed_dyads <- rbind(data.frame(dyads_B_PC,  type = "B_PC"),
                           data.frame(dyads_B_OM,  type = "B_OM"),
                           data.frame(dyads_PC_OM, type = "PC_OM"),
                           data.frame(dyads_PC_PR, type = "PC_PR"),
                           data.frame(dyads_OM_OM, type = "OM_OM"),
                           data.frame(dyads_OM_PR, type = "OM_PR"),
                           data.frame(dyads_PR_PR, type = "PR_PR"))
    
    # Every individual allowed dyad is equally likely to be drawn in the null model below (so guild-pairing types with more possible taxon pairs get proportionally more weight, not each type being equally likely).
    weights <- rep(1, nrow(allowed_dyads)) / nrow(allowed_dyads)
    
    ### Observed link counts and proportions.
    safe_sum <- function(mat, rows, cols) {
      if (length(rows) == 0 || length(cols) == 0) return(0)
      sum(mat[rows, cols, drop = FALSE])
    }
    
    links_B_PC  <- safe_sum(local_matrix, basal_idx, pcons_idx)
    links_B_OM  <- safe_sum(local_matrix, basal_idx, omni_idx)
    links_PC_OM <- safe_sum(local_matrix, pcons_idx, omni_idx)
    links_PC_PR <- safe_sum(local_matrix, pcons_idx, pred_idx)
    links_OM_OM <- safe_sum(local_matrix, omni_idx,  omni_idx)
    links_OM_PR <- safe_sum(local_matrix, omni_idx,  pred_idx)
    links_PR_PR <- safe_sum(local_matrix, pred_idx,  pred_idx)
    
    P_B_PC_obs  <- links_B_PC  / total_links
    P_B_OM_obs  <- links_B_OM  / total_links
    P_PC_OM_obs <- links_PC_OM / total_links
    P_PC_PR_obs <- links_PC_PR / total_links
    P_OM_OM_obs <- links_OM_OM / total_links
    P_OM_PR_obs <- links_OM_PR / total_links
    P_PR_PR_obs <- links_PR_PR / total_links
    

    ### Null models
    # Redistribute the same number of observed links (`total_links`) uniformly at random among all allowed dyads, 
    # `n_randomisations` times, and record what proportion of links fall into each guild-pairing type each time.
    null_P_B_PC  <- numeric(n_randomisations)
    null_P_B_OM  <- numeric(n_randomisations)
    null_P_PC_OM <- numeric(n_randomisations)
    null_P_PC_PR <- numeric(n_randomisations)
    null_P_OM_OM <- numeric(n_randomisations)
    null_P_OM_PR <- numeric(n_randomisations)
    null_P_PR_PR <- numeric(n_randomisations)
    
    for (k in seq_len(n_randomisations)) {
      sampled_types <- allowed_dyads$type[sample(nrow(allowed_dyads), size = total_links, prob = weights, replace = TRUE)]
      
      null_P_B_PC[k]  <- sum(sampled_types == "B_PC")  / total_links
      null_P_B_OM[k]  <- sum(sampled_types == "B_OM")  / total_links
      null_P_PC_OM[k] <- sum(sampled_types == "PC_OM") / total_links
      null_P_PC_PR[k] <- sum(sampled_types == "PC_PR") / total_links
      null_P_OM_OM[k] <- sum(sampled_types == "OM_OM") / total_links
      null_P_OM_PR[k] <- sum(sampled_types == "OM_PR") / total_links
      null_P_PR_PR[k] <- sum(sampled_types == "PR_PR") / total_links
    }
    
    ### Calculate z-scores
    safe_z <- function(obs, null_vec) {
      s <- sd(null_vec)
      if (is.na(s) || s == 0) return(NA)
      (obs - mean(null_vec)) / s
    }
    
    z_B_PC  <- safe_z(P_B_PC_obs,  null_P_B_PC)
    z_B_OM  <- safe_z(P_B_OM_obs,  null_P_B_OM)
    z_PC_OM <- safe_z(P_PC_OM_obs, null_P_PC_OM)
    z_PC_PR <- safe_z(P_PC_PR_obs, null_P_PC_PR)
    z_OM_OM <- safe_z(P_OM_OM_obs, null_P_OM_OM)
    z_OM_PR <- safe_z(P_OM_PR_obs, null_P_OM_PR)
    z_PR_PR <- safe_z(P_PR_PR_obs, null_P_PR_PR)
    
    ### Store results 
    results[[site_index]] <- data.frame(ugs = site_id,
                                        sp_richness = species_richness,
                                        total_links = total_links,
                                        n_basal = length(basal_idx),
                                        n_pcons = length(pcons_idx),
                                        n_omni = length(omni_idx),
                                        n_pred = length(pred_idx),
                                        L_basal_to_pcons = links_B_PC,
                                        L_basal_to_omni = links_B_OM,
                                        L_pcons_to_omni = links_PC_OM,
                                        L_pcons_to_pred = links_PC_PR,
                                        L_omni_to_omni = links_OM_OM,
                                        L_omni_to_pred = links_OM_PR,
                                        L_pred_to_pred = links_PR_PR,
                                        z_basal_to_pcons  = z_B_PC,
                                        z_basal_to_omni = z_B_OM,
                                        z_pcons_to_omni = z_PC_OM,
                                        z_pcons_to_pred = z_PC_PR,
                                        z_omni_to_omni = z_OM_OM,
                                        z_omni_to_pred = z_OM_PR,
                                        z_pred_to_pred = z_PR_PR)
  }
  bind_rows(results)
}

### 2. Load data and classify taxa into trophic guilds
# mw.df: the "metaweb", i.e. all known/reported trophic (and other) interactions, with each Source/Target given at a taxonomic rank (Species, Genus, Family, Order).
# fauna_matrix: species per site matrix for the observed animal taxa.
# flora_matrix: species per site matrix for the observed plant taxa.
# covariates: per-site covariates.
fauna_matrix <- read.csv("data/01_raw/fauna_matrix.csv", check.names = FALSE)
flora_matrix <- read.csv("data/01_raw/flora_matrix.csv", check.names = FALSE)
species_matrix <- merge(fauna_matrix, flora_matrix, by = "ugs")
metaweb <- read.csv("data/02_processed/focal_metaweb.csv")
covariates <- read.csv("data/01_raw/covariates.csv")

# Plant taxa and detritus are the basal resource groups. 
resources <- c(setdiff(colnames(flora_matrix), "ugs"), "Detritus")
species_matrix$Detritus <- 1

combined_matrix <- create_adjacency_matrix(metaweb)

# Classify every taxon that appears as a link Source into one of three guilds, based on what it is recorded eating:
# - omnivores: eats both a basal resource AND a fauna taxon 
# - primary consumers: eats a basal resource, but never a fauna taxon
# - predators: eats a fauna taxon, but never a basal resource
# We exclude "Pieris", as this names matches both a butterfly and a plant and is thus unreliable in determining which kind of consumer each species is
omnivores <- unique(metaweb$Source_Name[metaweb$Source_Name %in% metaweb$Source_Name[metaweb$Target_Name %in% resources] &
                                          metaweb$Source_Name %in% metaweb$Source_Name[metaweb$Target_Name %in% setdiff(colnames(fauna_matrix), "Pieris")]])
primary_consumers <- unique(metaweb$Source_Name[metaweb$Target_Name %in% resources &
                                                  !metaweb$Source_Name %in% metaweb$Source_Name[metaweb$Target_Name %in% setdiff(colnames(fauna_matrix), "Pieris")]])
predators <- unique(metaweb$Source_Name[metaweb$Target_Name %in% colnames(fauna_matrix) &
                                          !metaweb$Source_Name %in% metaweb$Source_Name[metaweb$Target_Name %in% resources]])

set.seed(123)

### 3. Compute per-site link probabilities and derived aggregate metrics.
combined_link_probs <- calculate_link_probabilities(species_matrix = species_matrix,
                                                    metaweb_matrix = combined_matrix,
                                                    primary_consumers = primary_consumers,
                                                    omnivores = omnivores,
                                                    predators = predators,
                                                    n_randomisations = 100)

# Merge with covariates.
link_probs <- merge(combined_link_probs, covariates, by = "ugs")

# Pool the three guild-pairing types that all feed into "omnivores" or into "predators" respectively, for a coarser 2-category view of consumer -> higher consumer links.
link_probs$L_all_consumer_to_omnivores <- link_probs$L_pcons_to_omni + link_probs$L_omni_to_omni + link_probs$L_omni_to_pred
link_probs$z_all_consumer_to_omnivores <- link_probs$z_pcons_to_omni + link_probs$z_omni_to_omni + link_probs$z_omni_to_pred
link_probs$L_all_consumer_to_predator <- link_probs$L_pcons_to_pred + link_probs$L_omni_to_pred + link_probs$L_pred_to_pred
link_probs$z_all_consumer_to_predator <- link_probs$z_pcons_to_pred + link_probs$z_omni_to_pred + link_probs$z_pred_to_pred

### 4. Prepare data formats for plotting. 


covariate_labels <- c(impervious_1000 = "% impervious (1000 m)", 
                      impervious_500 = "% impervious (500 m)", 
                      impervious_150 = "% impervious (150 m)")

# For each category x covariate combination, fit value ~ cov_value (linear model, or negative-binomial GLM for count data) and 
# flag whether the slope is significant at p < 0.05.
fit_sig <- function(df, model = c("lm", "nb")) {
  model <- match.arg(model)
  df %>%
    group_by(cons_type, covariate) %>%
    summarise(p_value = if (model == "lm") {
      summary(lm(value ~ cov_value))$coefficients["cov_value", "Pr(>|t|)"]
      } else {
        summary(MASS::glm.nb(value ~ cov_value))$coefficients["cov_value", "Pr(>|z|)"]
      },
      sig = ifelse(p_value < 0.05, "Significant", "NS"),
      .groups = "drop"
    )
}

# Build the long-format table for both the plotted values and the significance test.
rich_effects <- link_probs %>%
  select(ugs, all_of(names(covariate_labels)), n_pcons, n_omni, n_pred) %>%
  pivot_longer(cols = c(n_pcons, n_omni, n_pred), names_to = "cons_type", values_to = "value") %>%
  pivot_longer(cols = all_of(names(covariate_labels)), names_to = "covariate", values_to = "cov_value") %>%
  mutate(covariate = factor(covariate, levels = names(covariate_labels)),
         cons_type = recode(cons_type, 
                            "n_pcons" = "Primary\nconsumers",
                            "n_omni" = "Omnivores",
                            "n_pred" = "Predators")) %>%
  mutate(covariate = recode(covariate, 
                            impervious_1000 = "% impervious\n(1000 m)", 
                            impervious_500 = "% impervious\n(500 m)", 
                            impervious_150 = "% impervious\n(150 m)"))

rich_effects <- rich_effects %>%
  left_join(fit_sig(rich_effects, model = "nb"), by = c("cons_type", "covariate")) %>%
  mutate(consumer = as.character(cons_type))

# For each covariate x interaction-type combination, fit response ~ scaled covariate and return the standardised slope (beta) and its p-value. 
# Used to build the effect-size "bubble plots" (p2, p3). 
link_effects <- link_probs %>%
  pivot_longer(cols = c(L_basal_to_pcons, L_basal_to_omni, L_all_consumer_to_omnivores, L_all_consumer_to_predator),
               names_to = "interaction",
               values_to = "response") %>%
  pivot_longer(cols = all_of(names(covariate_labels)),
               names_to = "covariate",
               values_to = "cov_value") %>%
  group_by(covariate, interaction) %>%
  summarise(beta = coef(lm(response ~ scale(cov_value)))[2],
            p = summary(lm(response ~ scale(cov_value)))$coefficients[2, 4]) %>%
  ungroup() %>%
  mutate(effect_scaled = bruceR::scaler(beta, min = -1, max = 1), 
         covariate = factor(covariate, levels = rev(names(covariate_labels))),
         interaction = factor(interaction, levels = c("L_basal_to_pcons", "L_basal_to_omni", "L_all_consumer_to_omnivores", "L_all_consumer_to_predator"))) %>%
  mutate(interaction = recode(interaction, 
                              "L_basal_to_pcons" = "Plants →\nprimary consumers",
                              "L_basal_to_omni" = "Plants →\nomnivores",
                              "L_all_consumer_to_omnivores" = "Consumers →\nomnivores",
                              "L_all_consumer_to_predator" = "Consumers →\npredators"),
         covariate = recode(covariate, !!!covariate_labels)) 

z_effects <- link_probs %>%
  pivot_longer(cols = c(z_basal_to_pcons, z_basal_to_omni, z_all_consumer_to_omnivores, z_all_consumer_to_predator),
               names_to = "interaction",
               values_to = "response") %>%
  pivot_longer(cols = all_of(names(covariate_labels)),
               names_to = "covariate",
               values_to = "cov_value") %>%
  group_by(covariate, interaction) %>%
  summarise(beta = coef(lm(response ~ scale(cov_value)))[2],
            p = summary(lm(response ~ scale(cov_value)))$coefficients[2, 4]) %>%
  ungroup() %>%
  mutate(effect_scaled = bruceR::scaler(beta, min = -1, max = 1), 
         covariate = factor(covariate, levels = rev(names(covariate_labels))),
         interaction = factor(interaction, levels = c("z_basal_to_pcons", "z_basal_to_omni", "z_all_consumer_to_omnivores", "z_all_consumer_to_predator"))) %>%
  mutate(interaction = recode(interaction, 
                              "z_basal_to_pcons" = "Plants →\nprimary consumers",
                              "z_basal_to_omni" = "Plants →\nomnivores",
                              "z_all_consumer_to_omnivores" = "Consumers →\nomnivores",
                              "z_all_consumer_to_predator" = "Consumers →\npredators"), 
         covariate = recode(covariate, !!!covariate_labels))

### 5. Figure 2 panels
# Custom geom_smooth with alpha control on both line + ribbon (as ggplot doesn't natively map alpha separately to both).
smooth_alpha_layers <- function(df, sig_level, alpha_sig = c("Significant" = 1.0, "NS" = 0.2)) {
  
  df_sub <- df %>% filter(sig == sig_level)
  a <- alpha_sig[sig_level]
  
  list(geom_line(data = df_sub,
                 aes(x = cov_value, y = value, linetype = consumer),
                 stat      = "smooth",
                 method    = "lm",
                 linewidth = 1.1,
                 alpha     = a),
       geom_ribbon(data = df_sub,
                   aes(x = cov_value, y = value, linetype = consumer),
                   stat  = "smooth",
                   method = "lm",
                   color  = NA,
                   alpha  = a * 0.2))
}

p1 <- ggplot(mapping = aes(x = cov_value, y = value, linetype = consumer)) +
  smooth_alpha_layers(rich_effects, "Significant") +
  smooth_alpha_layers(rich_effects, "NS") +
  facet_wrap(~ covariate, scales = "free_x", nrow = 1) +
  scale_linetype_manual(values = c("Omnivores" = "dashed", "Primary\nconsumers" = "solid", "Predators" = "dotted"),
                        breaks = c("Primary\nconsumers", "Omnivores", "Predators"),
                        name = "Consumer type") +
  labs(x = "Covariate value", y = "Species richness", linetype = "Consumer\ntype") +
  theme_pubr(base_size = 16) +
  theme(strip.background = element_blank(),
        strip.text  = element_text(size = 16, hjust = 0),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16),
        plot.margin = margin(0, 0, 20, 0),
        axis.title.y = element_text(vjust = -35),
        legend.position = "right",
        plot.tag.position = c(0.01, 0.98))

p2 <- ggplot(link_effects, aes(x = interaction, y = covariate)) +
  geom_point(aes(size = abs(effect_scaled), fill = effect_scaled, alpha = ifelse(p < 0.05, 1, 0.2)), stroke = 1,
             shape = 21, color = "black") +
  scale_fill_gradient2(low = "#3D0066",
                       mid = "#EEC7FC",
                       high = "#FFCE1F",
                       name = "Effect\nsize", guide = guide_legend(override.aes = list(shape = 21, size = c(14, 8, 2, 8, 14)), reverse = TRUE)) +
  scale_size_continuous(range = c(2, 14),
                        name = "Effect\nsize",
                        guide = "none") +
  scale_alpha_identity() +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        plot.margin = margin(0, 0, 0, 0),
        plot.tag.position = c(0.01, 0.98))

p3 <- ggplot(z_effects, aes(x = interaction, y = covariate)) +
  geom_point(aes(size = abs(effect_scaled), fill = effect_scaled, alpha = ifelse(p < 0.05, 1, 0.2)), stroke = 1,
             shape = 21, color = "black") +
  scale_fill_gradient2(low = "#3D0066",
                       mid = "#EEC7FC",
                       high = "#FFCE1F",
                       name = "Effect\nsize", guide = guide_legend(override.aes = list(shape = 21, size = c(14, 8, 2, 8, 14)), reverse = TRUE)) +
  scale_size_continuous(range = c(2, 14),
                        name = "Effect\nsize",
                        guide = "none") +
  scale_alpha_identity() +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_blank(),
        plot.tag.position = c(-0.035, 0.98),
        plot.margin = margin(0, 0, 0, 30),
        legend.position = "top")

### 6. Combine panels into the final figure.
layout <- "
AA
BC
"

(p1 + p2 + p3 + theme(legend.position = "none")) +
  plot_layout(design = layout,
              heights = c(1, 0.8, 0.8),
              guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.key.spacing.y = unit(5, "pt"),
        legend.spacing.y = unit(30, "pt"),
        legend.box.margin = margin(30, 0, 30, 0),
        plot.tag = element_text(size = 20, face = "bold"))
