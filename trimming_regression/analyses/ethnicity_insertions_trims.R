library(GGally)
library(reshape)
library(data.table)
setDTthreads(1)
library(ggplot2)
library(ggpubr)
PROJECT_PATH = '/home/mrussel2'
OUTPUT_PATH = '/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output'
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

insertions = fread(file = '/fh/fast/matsen_e/shared/tcr-gwas/insertion_data/insertions_by_patient.tsv')[,-1]
trims = fread(file = '/fh/fast/matsen_e/shared/tcr-gwas/insertion_data/trim_by_patient.tsv')[,-1]
ethnicity = fread(file = paste0(PROJECT_PATH, "/tcr-gwas/_ignore/race_pcs_18Nov2020.txt"))
together_inserts = merge(ethnicity[,c('localID', 'race.g')], insertions, by.y = 'patient_id', by.x = 'localID')
together_trims = merge(ethnicity[,c('localID', 'race.g')], trims, by.y = 'patient_id', by.x = 'localID')

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
    yname = ifelse(data == 'together_trims', 'Mean total trimming length', 'Mean total insertion length')
    title = ifelse(data == 'together_trims', 'Mean total trimming length by racial group', 'Mean total insertion length by racial group')
    filename = ifelse(data == 'together_trims', 'figures/trim_average_by_ethnicity.pdf', 'figures/insert_average_by_ethnicity.pdf')
    average_all = df[, mean(total_trim), by = .(productivity_status)]
    colnames(average_all) = c('productivity_status', 'total_avg')

    ggplot(df, aes(x=race.g, y=total_trim, fill = race.g)) +
        facet_grid(cols = vars(productivity_status)) +
        geom_boxplot() +
        geom_jitter(shape=16, position=position_jitter(0.2), size = 3.5, alpha = 0.5) +
        stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.") +
        theme_classic() + 
        theme(text = element_text(size = 18), axis.text.x=element_text(angle = 45, vjust = 0.5), legend.position = "none") +
        ggtitle(title) +
        geom_hline(data = average_all, aes(yintercept = total_avg), size = 2, color = 'green', linetype = 2) +
        xlab('Racial Group') +
        ylab(yname) 


    ggsave(filename, plot = last_plot(), width = 18, height = 6, units = 'in', dpi = 500)
}


