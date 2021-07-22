# Manuscript plots 

The purpose of these scripts are to (1) visualize genome-wide and gene-level associations between SNPs and TCR repertoire features and (2) visualize other various figures for the manuscript.
The following scripts will create figures of interest: 

* [Plot SNP "dosage" effect](artemis_trimming_dosage.R) for the top trimming associated snp within artemis
* [Compare TRBD2 allele genotype and SNP genotype](d_allele_genotype_by_snp_genotype.R) for the top insertion association snp within dntt
* [Feature variations by gene](feature_by_gene_dist.R): this script creates plots which explore trimming/insertion distributions by gene
* [Plot genome-wide manhattan plot for gene usage](gene_usage_manhattan.R)
* [Plot distributions of genes which contain TRB associations for gene usage](gene_usage_tcrb_distributions.R) 
* [Plot N-insertion associations within the DNTT region](insertion_dntt_loci.R)
* [Plot genome-wide manhattan plot for N-insertion](insertion_j1_subset_manhattan.R) for only TCRs which contain TRBJ1
* [Plot genome-wide (and DNTT region) manhattan plot for N-insertion](insertion_multipanel_plot.R)
* [Compare extent of N-insertions by ancestry group](insertions_by_race.R)
* [Compare MAF by ancestry group](maf_by_race_boxplot.R) for the significantly associated SNPs within the DNTT region for N-insertion
* [Plot genome-wide manhattan plot for trimming, without gene-conditioning](naive_trimming_manhattan.R)
* [Plot genome-wide manhattan plot for P-nucleotide count](p_nuc_count_manhattan.R)
• [Plot genome-wide manhattan plot for the fraction of non-gene-trimmed TCRs which contain P-nucleotides](p_nuc_fraction_manhattan.R)
* [Create PCA parallel coordinates and scree plots](pca_plots.R): this script plots the parallel coordinate plots for the ancestry group PCAs (plotting PCA values versus PC number and colored by ancestry group) and plots the scree plot (proportion variance explained by PC number)
* [Plot genome-wide manhattan plot for trimming while conditioning on TRBD2 genotype](trimming_allele_status_correction_manhattan.R)
* [Plot trimming associations within the Artemis region](trimming_artemis_loci.R)
* [Plot genome-wide manhattan plot for trimming](trimming_j1_subset_manhattan.R) for only TCRs which contain TRBJ1
* [Plot genome-wide manhattan plot for trimming](trimming_manhattan.R)
* Plot SNP "dosage" effect for the validation cohort:
    * [Alpha chain insertions](validation_insertion_alpha_boxplots.R)
    * [Beta chain insertions](validation_insertion_beta_boxplots.R)
    * [Beta chain insertions--Discovery cohort](validation_insertion_snp_discovery_cohort.R)
    * [Alpha chain trimming](validation_trimming_alpha_boxplots.R)
    * [Beta chain trimming](validation_trimming_beta_boxplots.R)
    * [Beta chain trimming--Discovery cohort](validation_trimming_snp_discovery_cohort.R)

