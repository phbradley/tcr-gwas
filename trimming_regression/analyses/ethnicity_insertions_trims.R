library(GGally)
library(reshape)
library(data.table)
setDTthreads(1)
library(ggplot2)
library(ggpubr)
PROJECT_PATH = '/home/mrussel2'
OUTPUT_PATH = '/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output'
# import functions
source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts/run_bootstrap_regression_all_snps_functions_cluster.R"))
source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts/compile_data_functions.R"))


ethnicity = fread(file = '/home/mrussel2/tcr-gwas/_ignore/race_pcs_18Nov2020.txt')[, 1:3]
insertions = fread(file = '/fh/fast/matsen_e/shared/tcr-gwas/insertion_data/insertions_by_patient.tsv')[,-1]
trims = fread(file = '/fh/fast/matsen_e/shared/tcr-gwas/insertion_data/trim_by_patient.tsv')[,-1]

together_inserts = merge(ethnicity, insertions, by.y = 'patient_id', by.x = 'localID')
together_trims = merge(ethnicity, trims, by.y = 'patient_id', by.x = 'localID')

together_trims = together_trims[,total_trim := v_trim + d0_trim + d1_trim + j_trim]
together_inserts = together_inserts[,total_trim := vd_insert + dj_insert + vj_insert]

together_trims = together_trims[,-c('scanID', 'v_gene', 'd_gene', 'j_gene')][,lapply(.SD, mean), by = .(localID, productive, race.g)]
together_inserts = together_inserts[,-c('scanID', 'v_gene', 'd_gene', 'j_gene')][,lapply(.SD, mean), by = .(localID, productive, race.g)]

together_trims = together_trims[productive == 'TRUE', productivity_status := 'productive']
together_trims = together_trims[productive == 'FALSE', productivity_status := 'NOT_productive']

together_inserts = together_inserts[productive == 'TRUE', productivity_status := 'productive']
together_inserts = together_inserts[productive == 'FALSE', productivity_status := 'NOT_productive']

#look at mean inserts by racial group
for (data in c('together_trims', 'together_inserts')){
    df = get(data)
    title = ifelse(data == 'together_trims', 'Total Trimming Length Average by Ethnicity', 'Total Insertion Length Average by Ethnicity')
    filename = ifelse(data == 'together_trims', 'figures/trim_average_by_ethnicity.png', 'figures/insert_average_by_ethnicity.png')
    average_all = df[, mean(total_trim), by = .(productivity_status)]
    colnames(average_all) = c('productivity_status', 'total_avg')

    ggplot(df, aes(x=race.g, y=total_trim, fill = race.g)) +
        facet_grid(cols = vars(productivity_status)) +
        geom_boxplot() +
        geom_jitter(shape=16, position=position_jitter(0.2), size = 3.5, alpha = 0.5) +
        theme_classic() + 
        theme(text = element_text(size = 30), axis.text.x=element_text(angle = 45, vjust = 0.5), legend.position = "none") +
        ggtitle(title) +
        geom_hline(data = average_all, aes(yintercept = total_avg), size = 2, color = 'green', linetype = 2) +
        stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.")

    ggsave(filename, plot = last_plot(), width = 25, height = 10, units = 'in', dpi = 500)
}

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

genotype_dt = merge(ethnicity[,-c('scanID', 'race.s')], genotype_dt, by='localID')

genotype_sums = genotype_dt[,-c('localID')][,lapply(.SD, sum, na.rm=TRUE), by = .(race.g)]
genotype_nonNA_counts = genotype_dt[,-c('localID')][, lapply(.SD, function(x) 2*sum(!is.na(x))), by = race.g]
mafs = cbind(genotype_sums[,1], genotype_sums[,-1]/genotype_nonNA_counts[,-1])

reshaped_mafs = as.data.frame(t(as.data.frame(mafs)))
colnames(reshaped_mafs) = reshaped_mafs[1,]
reshaped_mafs = reshaped_mafs[-1,]
reshaped_mafs$snp = rownames(reshaped_mafs)
reshaped_mafs = sapply(reshaped_mafs, as.numeric)

pairs_plot = ggpairs(as.data.frame(reshaped_mafs), columns = 1:6, aes(alpha = 0.5))
ggsave('figures/maf_by_ethnicity.png', plot = pairs_plot, width = 15, height = 15, units = 'in', dpi = 500)

# calculate population maf
total_genotype_sums = genotype_dt[,-c('localID', 'race.g')][,lapply(.SD, sum, na.rm=TRUE)]
total_genotype_nonNA_counts = genotype_dt[,-c('localID', 'race.g')][, lapply(.SD, function(x) 2*sum(!is.na(x)))]
total_mafs = cbind(total_genotype_sums/total_genotype_nonNA_counts)
total_mafs = t(total_mafs)
total_mafs = as.data.frame(total_mafs)
total_mafs$snp = rownames(total_mafs)
colnames(total_mafs) = c('all_races', 'snp')
together_mafs = merge(total_mafs, reshape2::melt(as.data.frame(reshaped_mafs), id = 'snp'))

ggplot(together_mafs) + 
    geom_point(aes(x = all_races, y = value, color = variable), alpha = 0.5, size = 5) + 
    geom_abline(slope = 1, intercept = 0,  size = 2) +
    ggtitle('MAF by racial group versus by population') +
    theme_classic() +
    theme(text = element_text(size = 30)) + 
    xlab('MAF computed across all individuals') +
    ylab('MAF computed by racial group') +
    labs(fill = "Racial Group")

ggsave('figures/maf_by_ethnicity_population.pdf', plot = last_plot(), width = 15, height = 15, units = 'in', dpi = 500)


