library("tidyverse")
library("gdsfmt")
library("SNPRelate")
library("GWASTools")
library(data.table)
setDTthreads(threads = 1)
PROJECT_PATH = '/home/mrussel2'
source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts/compile_data_functions.R"))

dntt = find_snp_start_by_position(chromosome = 10, position1 = 98064085, position2 =98098321)
snp_data = snp_file_by_snp_start(dntt[1], dntt[2])
genotype_data = compile_all_genotypes(dntt[1], dntt[2])
genotype_data_asian = genotype_data[rownames(genotype_data) %in% ethnicity_data[race.g == 'Asian']$localID,]

ethnicity_data = fread(file = '/home/mrussel2/tcr-gwas/_ignore/race_pcs_18Nov2020.txt')   
ethnicity_data = ethnicity_data[, c('localID', 'scanID', 'race.g')]

list_of_snps = filtered_snps_by_maf(0.05, genotype_data_asian)
subjects = ethnicity_data[race.g == 'Asian']$scanID
snps = list_of_snps

snps_gds = snpgdsOpen(paste0(PROJECT_PATH, "/tcr-gwas/_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds"))

ld_dntt_asian = snpgdsLDMat(snps_gds, snp.id = snps, sample.id = subjects, slide = -1, with.id = TRUE)

plot_ld_dntt_asian = reshape2::melt(ld_dntt_asian$LD)

plot = ggplot(plot_ld_dntt_asian) + geom_tile(aes(x = Var1, y = Var2, fill = value^2))

ggsave(paste0('figures/ld/ld_dntt_asian.png'), plot = plot, height=10, width=10, units="in", dpi = 500)

