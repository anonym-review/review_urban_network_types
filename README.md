Data and code for the article titled "Urban densification and vegetation composition reshape ecological networks across trophic layers"

This repository contains the raw data, processed data, and R scripts used to build site-level food webs from a regional metaweb, calculate their structural properties, and relate those properties to environmental covariates.


*********************************


Workflow summary

Place raw files in data/01_raw/.
Run 01_prepare_mw.R to build focal_metaweb.csv.
Run 02_extract_foodweb_metrics.R to build fw_properties.csv and observed_fw_properties.csv.
Run fig2.R–fig4.R and figS1.R–figS12.R to reproduce the corresponding figures, each reading from data/01_raw/ and data/02_processed/ as needed.


*********************************


Raw data (data/01_raw/)

metaweb.csv — The regional metaweb: all known/reported trophic (and other) interactions among taxa, compiled from the literature, before any site-specific filtering. Each row is one interaction.

Source_Name / Target_Name — the two taxa involved in the interaction, at whatever taxonomic rank they were reported at. Source_Rank / Target_Rank — taxonomic rank of Source/Target (Species, Genus, Family, Order). Interaction_Type — type of interaction (e.g. Predation, Parasitism, Detritivore). For parasitism records, Source is the host and Target is the parasite in the raw file; this is corrected in 01_prepare_mw.R. observerd_metaweb.csv (sic — filename as provided) — Interactions that were directly observed in the field (rather than inferred from the literature-based metaweb above). Columns: row index, Target_Name, Source_Name. Used in 02_extract_foodweb_metrics.R as an alternative, non-inferred network for a sensitivity comparison against the main (inferred) food webs.

pollinator_interactions.csv — Observed plant–pollinator visitation records, one row per site x pollinator x plant observation.

ugs — site ID. Source_Name — pollinator taxon (genus). Target_Name — plant taxon (genus) visited. These are folded into the final metaweb in 01_prepare_mw.R, since pollination is not captured in the literature-based metaweb.

fauna_matrix.csv / flora_matrix.csv — Site x taxon occurrence matrices for animals and plants respectively. Rows = sites (ugs), columns = genera; a non-zero value indicates the genus was recorded at that site. flora_matrix.csv covers the basal (plant) resources; fauna_matrix.csv covers all consumer taxa.

fauna_taxo.csv / flora_taxo.csv — Taxonomy lookup tables for the animal and plant genera, respectively, with columns genus, family, order. Used to expand family- or order-level metaweb records down to the genus level for genera that aren't already resolved in the metaweb.

covariates.csv — Site-level environmental covariates, one row per site (ugs). See the data dictionary below.


*********************************


Processed data (data/02_processed/)

focal_metaweb.csv — The final, cleaned, genus-level metaweb used in all downstream analyses, produced by 01_prepare_mw.R. Columns: Source_Name, Target_Name (both genera, or Detritus). Restricted to genera actually observed in the study, with parasite/host direction corrected, problematic taxa dropped, family/order-level records expanded to genus level, and observed pollinator interactions added.

fw_properties.csv — Per-site food-web structural metrics, produced by 02_extract_foodweb_metrics.R, calculated separately for three network subsets (type column: Bipartite, Multitrophic, Combined — see below) and merged into one long table. See the data dictionary below for column definitions.

observed_fw_properties.csv — The same structural metrics as fw_properties.csv, but calculated from observerd_metaweb.csv (the directly observed, non-inferred interactions only) rather than the full inferred metaweb. Used as a robustness check. Same columns as fw_properties.csv, minus type (only one network is used here).


*********************************


Analysis scripts:

01_prepare_mw.R — Builds the final, analysis-ready focal_metaweb.csv from the raw metaweb.csv.

Identifies genera observed in the study that are absent from the metaweb at genus level, and expands the relevant family-level metaweb records down to genus level. Adds the observed plant–pollinator interactions. Writes focal_metaweb.csv.

02_extract_foodweb_metrics.R — Calculates structural properties of each site's local food web.

For each site, subsets the metaweb to the taxa observed at that site and calculates: species richness, number of links, connectance, link density, modularity, nestedness, r50, and the median/skewness/coefficient of variation of node degree. Each metric is also compared to a null model (same number of links, randomly placed among possible consumer links, preserving basal taxa as basal) to obtain a z-score. Splits the full metaweb into three subsets and repeats the site-level calculation for each: Bipartite (basal resource -> primary consumer links only), Multitrophic (consumer -> consumer/predator links only), and Combined (all links). Results are written to fw_properties.csv. Repeats the same calculation using the directly observed (non-inferred) metaweb, written to observed_fw_properties.csv.


*********************************


Figure scripts:

fig2.R — Produces Fig. 2. Calculates the observed proportion of links falling into each of the 7 possible directed guild-pairings (e.g. basal -> primary consumer, omnivore -> predator, etc.) and compares these to a null model.

fig3.R — Produces Fig. 3. Fits standardized linear models of degree-distribution metrics (median, skewness, coefficient of variation, all as z-scores) against each environmental covariate, separately for the Combined, Bipartite, and Multitrophic networks.

fig4.R — Produces Fig. 4. Fits a structural equation model (SEM) for each covariate and network type, following a fixed causal chain: covariate -> species richness -> connectance -> modularity/nestedness -> robustness (r50), with the covariate also allowed a direct path to every downstream variable. Bootstraps (1000 replicates) the standardized direct, indirect, and total effects of each covariate.

figS1.R — Produces Fig. S1. Calculates pairwise Pearson correlations (and p-values) among the standardized environmental covariates.

figS2.R — Produces Fig. S2. Same trophic-guild link-probability analysis as fig2.R, but relating link probabilities to percent impervious cover measured at three spatial scales (150 m, 500 m, 1000 m) instead of the full covariate set, as a robustness check for spatial scale.

figS3.R — Produces Fig. S3. Same degree-distribution model as fig3.R, but using percent impervious cover at the three spatial scales (150 m, 500 m, 1000 m) as the predictor, as a robustness check for spatial scale.

figS4.R — Produces Fig. S4. Same structural-equation-model (SEM) analysis as fig4.R, but using percent impervious cover at the three spatial scales (150 m, 500 m, 1000 m) as the predictor instead of the full covariate set — the SEM counterpart to figS2.R/figS3.R.

figS5.R — Produces Fig. S5. Rebuilds the focal metaweb from the raw data, then checks how well it predicts the independently observed plant–pollinator interactions in pollinator_interactions.csv. Splits interactions into three categories: predicted by the metaweb and observed ("Observed"), predicted but never observed ("Unobserved"), and observed but not predicted by the metaweb ("Unexpected").

figS6.R — Produces Fig. S6. Compares the inferred food web (fw_properties.csv) against the independently, directly observed food web (observed_fw_properties.csv) at the same sites. For each structural metric (connectance, nestedness, modularity, robustness, and the z-scored degree median/skewness/CV), calculates the Pearson correlation between the inferred and observed value across sites.

figS7.R — Produces Fig. S7. For both the inferred and independently observed food webs, fits linear models of each structural metric against each environmental covariate, then classifies each metric–covariate relationship by whether the two data sources agree in effect direction and statistical significance (e.g. "Same trend, both significant", "Different trend, observed-only significant").

figS8.R — Produces Fig. S8. For each network type (Bipartite/Multitrophic/Combined), calculates pairwise Pearson correlations among the three z-scored degree-distribution metrics (median, skewness, coefficient of variation).

figS9.R — Produces Fig. S9. For every structural metric and network type in fw_properties.csv, compares the distribution of observed values across sites to the distribution of the corresponding null-model mean values.

figS10.R — Produces Fig. S10. For each site, groups the observed fauna genera into broad taxonomic groups (Hemiptera, Hymenoptera, Diptera, Coleoptera, Arachnida, Gastropoda, Collembola, Orthoptera, Acari, and "Other taxa"), and calculates the fraction of the site's fauna richness contributed by each group.

figS11.R — Produces Fig. S11. Classifies fauna taxa into the same three trophic guilds used in fig2.R/figS2.R (primary consumer, omnivore, predator) based on the focal metaweb, then calculates the fraction of each site's consumer richness belonging to each guild. Fits a linear model of each guild's proportion against each environmental covariate.

figS12.R — Produces Fig. S12. Fits the same 18 predictor x network-type structural equation models (SEMs) as fig4.R/figS4.R (predictor -> species richness -> connectance -> modularity/nestedness -> robustness), but instead of a bootstrapped effect-size plot, draws a path diagram for each fitted SEM: edge width reflects the standardized coefficient, edge color indicates effect direction (blue = positive, red = negative, grey = non-significant), and line style indicates significance (solid p < 0.05, dashed p < 0.1).


*********************************


Data dictionary — covariates.csv

Variable Description

ugs Site ID

impervious_150 Percentage impervious surface in a 150 m buffer around the site

impervious_500 Percentage impervious surface in a 500 m buffer around the site

impervious_1000 Percentage impervious surface in a 1000 m buffer around the site

no.plants Number of plant species/genera recorded at the site (plant richness)

prop.wood Proportion of woody plant taxa at the site

prop.nonnative Proportion of non-native plant taxa at the site

maint_disturbance Index of disturbance intensity from site maintenance

maint_time Index of time investment in site maintenance


*********************************


Data dictionary — fw_properties.csv / observed_fw_properties.csv

Variable Description

ugs Site ID

type Network subset: Bipartite (basal resource -> primary consumer links), Multitrophic (consumer -> consumer/predator links), or Combined (all links). (fw_properties.csv only — observed_fw_properties.csv covers a single network.)

no_isolated_nodes Number of taxa removed for having no incoming or outgoing links (excluding basal resources)

species_richness Number of taxa (nodes) in the site's local food web

number_of_links Number of interactions (edges) in the site's local food web

connectance number_of_links / species_richness^2

link_density number_of_links / species_richness

modularity Observed network modularity (edge-betweenness community detection)

nestedness Observed network nestedness (NODF/UNODF)

r50 Observed robustness: mean proportion of taxa that must be sequentially removed (via random primary removals, allowing secondary extinctions) before 50% of all taxa are lost

median_degree, skew_degree, cv_degree Median, skewness, and coefficient of variation of node degree, observed

*_null_mean Mean of the corresponding metric across degree/link-count-preserving null networks

*_z Z-score of the observed metric relative to its null distribution: (observed - null_mean) / null_sd
