# Analysis functions and pipelines

This directory contains scripts for analysis and supplementary pipelines: 

* [Various analysis related functions](analysis_functions.R) specifically related to identifying association statistics for certain genome regions and phenotypes
* [Determine analysis results for gene region and phenotype](determine_sig_snps_stats_by_gene.R)
* [Count the number of associations by phenotype](identify_genome_wide_association_count.R)
* Determine genomic inflation value--lambda for [gene-conditioned](calculate_lambdas_gene_conditioned.R) and [non-gene-conditioned](calculate_lambdas.R) cases
* [Execute random bootstraps](bootstrap_analysis.R) to prepare for gene-conditioned lambda calculations

The following directories contain supplementary pipelines for analysis:
* [Conditional analysis for DNTT and Artemis region](conditional_analysis_scripts) to identify independent snp associations within specified gene regions and phenotypes
* [Validation cohort analyses](validation_cohort_analysis) for several SNPs (identical pipeline to the main genome-wide analysis with a few data specific adjustments)
