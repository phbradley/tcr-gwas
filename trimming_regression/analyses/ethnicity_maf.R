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
source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts/compile_data_functions.R"))


find_significant_snps_by_gene <- function(trim_type, maf_cutoff, signficance_cutoff){
    maf_data = fread(paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/maf_all_snps.tsv'), sep = "\t", fill=TRUE, header = TRUE)[,-c(1)]
    sig_data = compile_manhattan_plot_data(trim_type, maf_data, maf_cutoff, pca_structure_correction = 'False', pca_type = 'none')
    return(sig_data[pvalue < signficance_cutoff])
}

determine_true_minor_allele <- function(snp, trim_genotype_dt){
    columns = c('total_insert', paste(snp))
    simplified_dt = trim_genotype_dt[productive == 'TRUE'][,..columns]
    regression = lm(simplified_dt$total_insert ~ simplified_dt[[paste(snp)]])
    slope = coef(regression)['simplified_dt[[paste(snp)]]']
    minor_allele = ifelse(slope < 0, 2, 0)
    return(minor_allele)
}

calculate_maf <- function(snp, minor_allele, race, genotype_dt){
    columns = c('localID', 'race.g', paste(snp))
    if (race != 'all'){
        data = genotype_dt[race.g == race,..columns]
    } else {
        data = genotype_dt[,..columns]
    }

    if (minor_allele != 2){
        major_hom = data[get(snp) == 2][, (snp) := 0]
        minor_het = data[get(snp) == 1]
        minor_hom = data[get(snp) == 0][, (snp) := 2]
        data = rbind(minor_hom, minor_het, major_hom)
    }
    total_possible_alleles = 2*sum(!is.na(data[[snp]]))
    observed_alleles = sum(data[[snp]], na.rm = TRUE)
    maf = observed_alleles/total_possible_alleles
    return(maf)
} 

# import race data
ethnicity = fread(file = paste0(PROJECT_PATH, "/tcr-gwas/_ignore/race_pcs_18Nov2020.txt"))

# find snps and genotype data 
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

#compile trimming data
trimming = compile_condensed_trimming_data('vd_insert', 'by_patient')
total_trims = as.data.table(trimming)[,total_insert := dj_insert + vd_insert, by=.(localID, productive)]
trim_genotype = merge(genotype_dt, total_trims[, c('localID', 'productive', 'total_insert')])

maf_dt = data.table()
for (snp in maf_satisfactory_snps[maf_satisfactory_snps != 'localID']){
    minor_allele = determine_true_minor_allele(snp, trim_genotype)
    for (race in c(unique(trim_genotype$race.g), 'all')){
        maf = calculate_maf(snp, minor_allele, race, trim_genotype)
        temp = data.table(snp = snp, race = race, maf = maf, minor_allele = minor_allele)
        maf_dt = rbind(maf_dt, temp)
    }
}

# find significant snps from insertion analysis
bonferroni = 0.05/35481497
sig_snps = rbind(find_significant_snps_by_gene(trim_type = 'dj_insert', maf_cutoff = 'False', signficance_cutoff = bonferroni), find_significant_snps_by_gene(trim_type = 'vd_insert', maf_cutoff = 'False', signficance_cutoff = bonferroni))
sig_snps_list = sub('snp', '', unlist(unique(sig_snps$snp)))

sig_snps_together_mafs = maf_dt[snp %in% sig_snps_list]

# reshape data 
all_race_mafs = sig_snps_together_mafs[race == 'all']
by_race_mafs = sig_snps_together_mafs[race != 'all']
sig_snps_together_reshaped = merge(by_race_mafs, all_race_mafs, by = c('snp', 'minor_allele'))[,c('snp', 'race.x', 'maf.x', 'maf.y')]
colnames(sig_snps_together_reshaped) = c('snp', 'race', 'maf_by_race', 'maf_all_races')


ggplot(sig_snps_together_reshaped) + 
    geom_point(aes(x = maf_all_races, y = maf_by_race, color = race), alpha = 0.5, size = 5) + 
    geom_abline(slope = 1, intercept = 0,  size = 2) +
    ggtitle('MAF by racial group versus by population \nfor DNTT significant snps') +
    theme_classic() +
    theme(text = element_text(size = 40)) + 
    xlab('MAF computed across all individuals') +
    ylab('MAF computed by racial group') +
    labs(fill = "Racial Group")

ggsave(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/figures/maf_by_race.pdf', plot = last_plot(), width = 12, height = 12, units = 'in', dpi = 750, device = 'pdf')


