library(dplyr)
library(bipartite)
library(igraph)
library(moments)
library(future.apply)
library(parallel)
library(UNODF)

### 1. Robustness: proportion of primary removals needed to lose 50% of all taxa
# Repeatedly removes one non-basal taxon at random, cascades any consumers that lose all of their resources as a result ("secondary extinctions"), 
# and records the fraction of the network removed (in terms of primary removals only) at the point where >= 50% of all taxa are gone. 
# Averaged over `n_replicates` random removal orders to get a stable R50 estimate.
network_robustness <- function(diet_matrix,
                               n_replicates = 100,
                               remove_basal = FALSE,
                               basal_taxa) {
  
  # Drop non-basal taxa that already have no resources (colSums == 0), and repeat since removing one such taxon can strand another consumer with zero resources.
  isolated_taxa <- setdiff(names(which(colSums(diet_matrix) == 0)), basal_taxa)
  
  while (length(isolated_taxa) > 0) {
    diet_matrix <- diet_matrix[!rownames(diet_matrix) %in% isolated_taxa,
                               !colnames(diet_matrix) %in% isolated_taxa]
    
    isolated_taxa <- setdiff(names(which(colSums(diet_matrix) == 0)), basal_taxa)
  }
  
  original_richness <- nrow(diet_matrix)

  # Run n_replicates independent random-removal sequences in parallel.
  r50_values <- future_lapply(seq_len(n_replicates), function(i) {
    current_matrix <- diet_matrix
    primary_removals <- 0
    total_losses <- 0
    r50_step <- NA_real_
    
    while (nrow(current_matrix) > 0) {
      basal_nodes <- which(rownames(current_matrix) %in% basal_taxa)
      available_nodes <- seq_len(nrow(current_matrix))
      
      # Preserve basal taxa unless `remove_basal` is explicitly requested.
      removable_nodes <- if (!remove_basal && length(basal_nodes) < length(available_nodes)) {
        setdiff(available_nodes, basal_nodes)
      } else {
        available_nodes
      }
      
      removed_node <- sample(removable_nodes, 1)
      primary_removals <- primary_removals + 1
      total_losses <- total_losses + 1
      
      # Remove one taxon directly (both its row and column, since the matrix is square).
      current_matrix <- current_matrix[-removed_node, -removed_node, drop = FALSE]
      
      # Cascade: remove any consumers that have lost all of their resources.
      while (nrow(current_matrix) > 0) {
        basal_nodes <- which(rownames(current_matrix) %in% basal_taxa)
        
        secondary_extinctions <- setdiff(which(colSums(current_matrix) == 0), basal_nodes)
        
        if (length(secondary_extinctions) == 0) {
          break
        }
        
        total_losses <- total_losses + length(secondary_extinctions)
        
        current_matrix <- current_matrix[-secondary_extinctions,
                                         -secondary_extinctions, 
                                         drop = FALSE]
      }
      
      # Record the first primary-removal step at which half the network is lost.
      if (is.na(r50_step) && total_losses >= 0.5 * original_richness) {
        r50_step <- primary_removals
      }
    }
    
    r50_step / original_richness
  }, future.seed = TRUE)
  
  mean(unlist(r50_values), na.rm = TRUE)
}

### 2. Helper functions
# Build a square, binary (presence/absence) taxon x taxon adjacency matrix from a long-format metaweb (one row per Source -> Target interaction).
create_adjacency_matrix <- function(metaweb) {
  # Create adjacency matrix of metaweb (species x species)
  adjacency_matrix <- matrix(0,
                             nrow = length(unique(c(metaweb$Target_Name, metaweb$Source_Name))),
                             ncol = length(unique(c(metaweb$Target_Name, metaweb$Source_Name))),
                             dimnames = list(sort(unique(c(metaweb$Target_Name, metaweb$Source_Name))), 
                                             sort(unique(c(metaweb$Target_Name, metaweb$Source_Name)))))
  
  # Fill adjacency matrix using interaction counts (target = prey, source = predator)
  interaction_counts <- table(metaweb$Target_Name, metaweb$Source_Name)
  
  adjacency_matrix[rownames(interaction_counts),
                   colnames(interaction_counts)] <- interaction_counts
  
  # Binarize the matrix (presence/absence of interaction)
  adjacency_matrix <- ifelse(adjacency_matrix > 0, 1, 0)
  
  return(adjacency_matrix)
}

### 3. Calculate observed metrics and null-model z-scores for every site
# For each site: build its local food web (the subset of the metaweb spanned by taxa actually observed there), compute observed structural metrics, 
# then compare them against `n_randomisations` degree-preserving null networks to get a z-score for each metric.
compute_foodweb_metrics <- function(species_matrix,
                                    metaweb_matrix,
                                    resources,
                                    n_randomisations = 100,
                                    remove_isolated_nodes = TRUE,
                                    type = "unipartite") {
  
  results <- vector("list", nrow(species_matrix))
  
  for (site_index in seq_len(nrow(species_matrix))) {
    site_id <- species_matrix$ugs[site_index]
    message("Processing site: ", site_id)
    
    # Keep taxa observed at this site that are also present in the metaweb.
    present_taxa <- species_matrix[site_index, ] %>%
      select(-ugs) %>%
      dplyr::select_if(colnames(.) %in% c(colnames(metaweb_matrix), rownames(metaweb_matrix))) %>%
      dplyr::select_if(colSums(.) >= 1) %>%
      colnames()
    
    local_matrix <- metaweb_matrix[rownames(metaweb_matrix) %in% present_taxa,
                                   colnames(metaweb_matrix) %in% present_taxa]
    
    # Remove taxa with neither incoming nor outgoing interactions, except
    # designated resource taxa such as plants and detritus (which are allowed
    # to have zero incoming links since they are basal).
    isolated_taxa <- setdiff(names(which(rowSums(local_matrix) == 0 & colSums(local_matrix) == 0)), resources)
    
    if (remove_isolated_nodes && length(isolated_taxa) > 0) {
      present_taxa <- setdiff(present_taxa, isolated_taxa)
      local_matrix <- metaweb_matrix[rownames(metaweb_matrix) %in% present_taxa,
                                     colnames(metaweb_matrix) %in% present_taxa]
    }
    
    local_network <- graph_from_adjacency_matrix(local_matrix, mode = "directed")
    basal_nodes <- which(colSums(local_matrix) == 0)
    
    # Basic network structure.
    species_richness <- vcount(local_network)
    number_of_links <- ecount(local_network)
    connectance <- number_of_links / species_richness^2
    link_density <- number_of_links / species_richness
    
    degree_values <- igraph::degree(local_network, mode = "all")
    degree_values <- degree_values[degree_values > 0]
    
    # Observed topology and robustness metrics.
    modularity <- modularity(cluster_edge_betweenness(local_network))
    
    if(type == "bipartite") {
      nestedness <- unname(nested(local_matrix[basal_nodes, setdiff(present_taxa, basal_nodes)], method = "NODF"))
    } else {
      nestedness <- unname(unlist(unodf(local_matrix, selfloop = FALSE)[1]))
    }
    
    r50 <- network_robustness(local_matrix,
                             n_replicates = 1000,
                             remove_basal = TRUE,
                             basal_taxa = colnames(local_matrix)[basal_nodes])

    observed_metrics <- c(modularity = modularity,
                          nestedness = nestedness,
                          r50 = r50,
                          median_degree = median(degree_values),
                          skew_degree = skewness(degree_values),
                          cv_degree = sd(degree_values) / mean(degree_values))
    
    consumer_nodes <- setdiff(seq_len(ncol(local_matrix)), basal_nodes)
    
    # Null networks: same number of links as the observed web, placed at random among consumer columns only (so basal taxa can never become consumers).
    null_metrics <- future_lapply(seq_len(n_randomisations), function(i) {
      random_matrix <- matrix(0, nrow = nrow(local_matrix), ncol = ncol(local_matrix), dimnames = dimnames(local_matrix))
      
      eligible_cells <- which(col(random_matrix) %in% consumer_nodes)
      
      random_matrix[sample(eligible_cells, sum(local_matrix), replace = FALSE)] <- 1
      
      random_network <- graph_from_adjacency_matrix(random_matrix, mode = "directed")
      
      random_degree_values <- igraph::degree(random_network, mode = "all")
      random_degree_values <- random_degree_values[random_degree_values > 0]
      
      random_modularity <- modularity(cluster_edge_betweenness(random_network))
      
      if(type == "bipartite") {
        random_nestedness <- unname(nested(random_matrix[apply(random_matrix, 1, "sum") > 0, 
                                                         apply(random_matrix, 2, "sum") > 0], method = "NODF"))
      } else {
        random_nestedness <- unname(unlist(unodf(random_matrix, selfloop = FALSE)[1]))
      }
      
      random_r50 <- network_robustness(random_matrix,
                                       n_replicates = 1,
                                       remove_basal = TRUE,
                                       basal_taxa = colnames(local_matrix)[basal_nodes])
      c(modularity = random_modularity,
        nestedness = random_nestedness,
        r50 = random_r50,
        median_degree = median(random_degree_values),
        skew_degree = skewness(random_degree_values),
        cv_degree = sd(random_degree_values) / mean(random_degree_values))
    }, future.seed = TRUE)
    
    null_matrix <- do.call(rbind, null_metrics)
    null_mean <- colMeans(null_matrix, na.rm = TRUE)
    null_sd <- apply(null_matrix, 2, sd, na.rm = TRUE)
    
    # Compare observed values with the null-model distribution.
    z_scores <- (observed_metrics - null_mean) / null_sd
    names(z_scores) <- paste0(names(z_scores), "_z")
    
    results[[site_index]] <- as.data.frame(t(c(ugs = site_id,
                                               no_isolated_nodes = length(isolated_taxa),
                                               species_richness = species_richness,
                                               number_of_links = number_of_links,
                                               connectance = connectance,
                                               link_density = link_density,
                                               observed_metrics,
                                               setNames(null_mean, paste0(names(null_mean), "_null_mean")),
                                               z_scores)))
  }
  
  bind_rows(results)
}

### 4. Load data and construct the three focal metawebs
fauna_matrix <- read.csv("data/01_raw/fauna_matrix.csv", check.names = FALSE)
flora_matrix <- read.csv("data/01_raw/flora_matrix.csv", check.names = FALSE)
species_matrix <- merge(fauna_matrix, flora_matrix, by = "ugs")

# Plant taxa and detritus are the basal resource groups. "ugs" (the site ID
# column) is excluded here since it's not a taxon.
resources <- c(setdiff(colnames(flora_matrix), "ugs"), "Detritus")
species_matrix$Detritus <- 1

metaweb <- read.csv("data/02_processed/focal_metaweb.csv")

combined_matrix <- create_adjacency_matrix(metaweb)

# Split the full metaweb into a strictly bipartite "who eats a basal resource"
# subweb and a unipartite "who eats a consumer" (omnivore/predator) subweb.
primary_consumer_matrix <- create_adjacency_matrix(metaweb %>% filter(Target_Name %in% resources))
omnivore_predator_matrix <- create_adjacency_matrix(metaweb %>% filter(!Target_Name %in% resources))

### 5. Run the food-web metrics for each metaweb type and combine the results
# Use all available cores except one, while retaining reproducible randomisation.
# max(1, ...) guards against detectCores() - 1 == 0 on single-core machines.
plan(multisession, workers = max(1, detectCores() - 1))
RNGkind("L'Ecuyer-CMRG")
set.seed(1)

combined_foodwebs <- compute_foodweb_metrics(species_matrix, combined_matrix, resources, n_randomisations = 100, type = "unipartite")
bipartite_foodwebs <- compute_foodweb_metrics(species_matrix, primary_consumer_matrix, resources, n_randomisations = 100, type = "bipartite")
multitrophic_foodwebs <- compute_foodweb_metrics(species_matrix, omnivore_predator_matrix, resources, n_randomisations = 100, type = "unipartite")

foodweb_properties <- bind_rows(mutate(bipartite_foodwebs, type = "Bipartite"),
                                mutate(multitrophic_foodwebs, type = "Multitrophic"),
                                mutate(combined_foodwebs, type = "Combined"),)

write.csv(foodweb_properties, "data/02_processed/fw_properties.csv", row.names = FALSE)

### 6. Repeat the analysis using the separately observed (non-inferred) metaweb
observed_metaweb <- read.csv("data/01_raw/observed_metaweb.csv")

observed_matrix <- create_adjacency_matrix(observed_metaweb)

observed_foodwebs <- compute_foodweb_metrics(species_matrix, observed_matrix, resources, n_randomisations = 100, type = "bipartite")

write.csv(observed_foodwebs, "data/02_processed/observed_fw_properties.csv", row.names = FALSE)

# Restore sequential execution once the analysis is finished.
plan(sequential)