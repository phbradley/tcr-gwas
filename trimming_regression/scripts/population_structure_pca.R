library("tidyverse")
library("SNPRelate")
library(data.table)
setDTthreads(threads = 1)

snps_gds = snpgdsOpen("/fh/fast/matsen_e/shared/tcr-gwas/snp_data/HSCT_comb_geno_combined_v03_tcr.gds")

# Before running PCA, we use LD pruning to select a set of independent SNPs for analysis
ld_pruned_snps = snpgdsLDpruning(snps_gds)

pruned_snp_list = unlist(unname(ld_pruned_snps))

# Generate pca from pruned snp list
pca <- snpgdsPCA(snps_gds, snp.id = pruned_snp_list, num.thread=2)

# The code below shows how to calculate the percent of variation is accounted for by the top principal components.
pc.percent <- pca$varprop*100

# Create dataframe with pca information
pca_dt = data.table(sample_id = pca$sample.id, EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], EV3 = pca$eigenvect[,3], EV4 = pca$eigenvect[,4], EV5 = pca$eigenvect[,5], EV6 = pca$eigenvect[,6], EV7 = pca$eigenvect[,7], EV8 = pca$eigenvect[,8], EV9 = pca$eigenvect[,9], EV10 = pca$eigenvect[,10],stringsAsFactors = FALSE)

snpgdsClose(snps_gds)

write.table(pca_dt, file='/home/mrussel2/tcr-gwas/_ignore/snp_data/population_structure_pca_by_LD_snps.tsv', quote=FALSE, sep='\t', col.names = NA)

