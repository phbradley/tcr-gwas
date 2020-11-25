library(GGally)
library(reshape)
library(data.table)
setDTthreads(1)
library(ggplot2)
library(ggpubr)
PROJECT_PATH = '/home/mrussel2'
OUTPUT_PATH = '/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output'
TRIM_TYPE = 'vd_insert'
# import functions
source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts/config.R"))
source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts/run_bootstrap_regression_all_snps_functions_cluster.R"))
source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts/compile_data_functions.R"))
source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts/compile_regression_functions.R"))
source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts/manha_visualization.R"))


find_significant_snps_by_gene <- function(trim_type, maf_cutoff, signficance_cutoff){
    maf_data = fread(paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/maf_all_snps.tsv'), sep = "\t", fill=TRUE, header = TRUE)[,-c(1)]
    sig_data = compile_manhattan_plot_data(trim_type, maf_data, maf_cutoff, pca_structure_correction = 'False', pca_type = 'none')
    return(sig_data[pvalue < signficance_cutoff])
}

ethnicity = fread(file = paste0(PROJECT_PATH, "/tcr-gwas/_ignore/race_pcs_18Nov2020.txt"))

# look at MAFs
dntt = find_snp_start_by_position(chromosome = 10, position1 = 98064085, position2 =98098321)
snp_start = dntt[1]
count = dntt[2]
snp_data = snp_file_by_snp_start(snp_start = snp_start, count)
genotype_data = compile_all_genotypes(snp_start = snp_start, count)
genotype_data_filtered = as.data.frame(remove_matrix_column_by_genotype(genotype_data))
genotype_data_filtered$localID = rownames(genotype_data_filtered)
genotype_dt = as.data.table(genotype_data_filtered)

# subset snps by those that pass the 0.05 MAF cutoff for total population
list_of_snps = filtered_snps_by_maf(MAF_CUTOFF = 0.05, genotype_list = genotype_data_filtered)

maf_satisfactory_snps = c(as.vector(as.character(list_of_snps)), 'localID')
genotype_dt = genotype_dt[,..maf_satisfactory_snps]

genotype_dt = merge(ethnicity[,c('localID', 'race.g')], genotype_dt, by='localID')

genotype_sums = genotype_dt[,-c('localID')][,lapply(.SD, sum, na.rm=TRUE), by = .(race.g)]
genotype_nonNA_counts = genotype_dt[,-c('localID')][, lapply(.SD, function(x) 2*sum(!is.na(x))), by = race.g]
mafs = cbind(genotype_sums[,1], genotype_sums[,-1]/genotype_nonNA_counts[,-1])

reshaped_mafs = as.data.frame(t(as.data.frame(mafs)))
colnames(reshaped_mafs) = reshaped_mafs[1,]
reshaped_mafs = reshaped_mafs[-1,]
reshaped_mafs$snp = rownames(reshaped_mafs)
reshaped_mafs = sapply(reshaped_mafs, as.numeric)

# pairs_plot = ggpairs(as.data.frame(reshaped_mafs), columns = 1:6, aes(alpha = 0.5))
# ggsave('figures/maf_by_ethnicity.png', plot = pairs_plot, width = 15, height = 15, units = 'in', dpi = 500)

# calculate population maf
total_genotype_sums = genotype_dt[,-c('localID', 'race.g')][,lapply(.SD, sum, na.rm=TRUE)]
total_genotype_nonNA_counts = genotype_dt[,-c('localID', 'race.g')][, lapply(.SD, function(x) 2*sum(!is.na(x)))]
total_mafs = cbind(total_genotype_sums/total_genotype_nonNA_counts)
total_mafs = t(total_mafs)
total_mafs = as.data.frame(total_mafs)
total_mafs$snp = rownames(total_mafs)
colnames(total_mafs) = c('all_races', 'snp')
together_mafs = merge(total_mafs, reshape2::melt(as.data.frame(reshaped_mafs), id = 'snp'))

bonferroni = 0.05/35481497
sig_snps = rbind(find_significant_snps_by_gene(trim_type = 'dj_insert', maf_cutoff = 'False', signficance_cutoff = bonferroni), find_significant_snps_by_gene(trim_type = 'vd_insert', maf_cutoff = 'False', signficance_cutoff = bonferroni))
sig_snps_list = sub('snp', '', unlist(unique(sig_snps$snp)))

sig_snps_together_mafs = as.data.table(together_mafs)[snp %in% sig_snps_list]

ggplot(sig_snps_together_mafs) + 
    geom_point(aes(x = all_races, y = value, color = variable), alpha = 0.5, size = 5) + 
    geom_abline(slope = 1, intercept = 0,  size = 2) +
    ggtitle('MAF by racial group versus by population \nfor DNTT significant snps') +
    theme_classic() +
    theme(text = element_text(size = 18)) + 
    xlab('MAF computed across all individuals') +
    ylab('MAF computed by racial group') +
    labs(fill = "Racial Group")

ggsave('figures/maf_by_ethnicity_population.pdf', plot = last_plot(), width = 8, height = 6, units = 'in', dpi = 500)


