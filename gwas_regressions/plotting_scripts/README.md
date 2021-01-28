# Plotting scripts

The purpose of these scripts are to (1) visualize genome-wide and gene-level associations between SNPs and TCR repertoire features and (2) visualize other various figures for the manuscript.
The following scripts will create figures of interest: 

* [Feature variations by gene](plotting_scripts/feature_by_gene_dist.R): this script creates plots which explore trimming/insertion distributions by gene
* [Insertions by ancestry group](plotting_scripts/insertions_by_race.R): this script creates a boxplot of average insertions per individual repertoire by productivity and ancestry group
    * This script also computes the pvalues presented in the paper (associated with insertions by ancestry group)
* [Minor allele frequency by ancestry group](plotting_scripts/maf_by_race.R): this script computes minor allele frequency by population and by ancestry group for SNPs within the DNTT locus with significant associations for the insertion analysis (without population structure correction) and plots the ancestry group MAFs versus population MAF
* [Minor allele frequency by ancestry group boxplot](plotting_scripts/maf_by_race_boxplot.R): again, this script computes minor allele frequency by population and by ancestry group for SNPs within the DNTT locus with significant associations for the insertion analysis (without population structure correction) and plots the ancestry group MAFs as a boxplot. 
    * This script also computes the pvalues presented in the paper comparing ancestry group MAF to population MAF
* [PCA parallel coordinates and scree plots](plotting_scripts/pca_plots.R): this script plots the parallel coordinate plots for the ancestry group PCAs (plotting PCA values versus PC number and colored by ancestry group) and plots the scree plot (proportion variance explained by PC number)
* [Plot analyses genome-wide via manhattan plot](plotting_scripts/plot_genome_wide.R): this script plots manhattan plots for all of the possible phenotypes. Use the equivalent [rmarkdown](plotting_scripts/plot_genome_wide.Rmd) document to plot individual manhattan plots (if desired).
* [Plot analyses for gene-of-interest via manhattan plot](plotting_scripts/plot_gene.R): this script plots manhattan plots for various genes/phenotypes of interest. Use the equivalent [rmarkdown](plotting_scripts/plot_gene.Rmd) document to plot individual manhattan plots (if desired).
* [Plot significant associations for loci of interest](plotting_scripts/plot_loci.R): this script plots various significant associations for genes/phenotypes of interest. Use the equivalent [rmarkdown](plotting_scripts/plot_loci.Rmd) document to plot individual plots (if desired).



