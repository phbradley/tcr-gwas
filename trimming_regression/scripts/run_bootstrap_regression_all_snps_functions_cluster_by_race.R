library('plyr')
library("tidyverse")
library('foreach')
library('doParallel')
library('RhpcBLASctl')
library('stringr')
library('data.table')
setDTthreads(threads = 1)
omp_set_num_threads(1)
blas_set_num_threads(1)


source(paste0(project_path, "/tcr-gwas/trimming_regression/scripts/regression_functions.R"))
source(paste0(project_path, "/tcr-gwas/trimming_regression/scripts/bootstrap_functions.R"))
source(paste0(project_path, "/tcr-gwas/trimming_regression/scripts/compile_data_functions.R"))
source(paste0(project_path, "/tcr-gwas/trimming_regression/scripts/execute_regression_function.R"))

# This function runs regressions on the cluster
run_snps_trimming_snp_list_cluster_by_race <- function(snp_list, genotype_list, trim_type, pca_structure_correction, pca_type, write_table, ncpus, maf_cutoff, random_effects){
 
    # parse type (either 'insert' or 'trim')
    stopifnot(trim_type %in% c('v_trim', 'd0_trim', 'd1_trim', 'j_trim', 'vd_insert', 'dj_insert', 'vj_insert'))

    set_regression_parameters(trim_type)
    
    regression_dataframe = data.frame()

    # import condensed trimming file
    condensed_trimming_data = compile_condensed_trimming_data(trim_type, d_infer, condensing)
    # filter out small repertoires
    trimming_data = remove_small_repertoire_observations(condensed_trimming_data, productive_log10_count_cutoff= 4.25, NOT_productive_log10_count_cutoff=3.5, trim_type, d_infer)
    ethnicity_data = fread(file = '/home/mrussel2/tcr-gwas/_ignore/snp_data/ethnicity_data.csv')   

    together = merge(trimming_data, ethnicity_data, by.x = 'localID', by.y = 'id')

    #filter out snps with maf below cutoff
    if (maf_cutoff == "False"){
        list_of_snps = as.numeric(colnames(genotype_list)[colnames(genotype_list) != 'localID'])
    } else {
        maf_file_name = paste0(output_path, '/maf_all_snps.tsv')
        if (!file.exists(maf_file_name)){
            system(command = paste0("Rscript ", project_path, "/tcr-gwas/trimming_regression/scripts/calculate_maf.R ", ncpus))
        }
        maf_data = fread(maf_file_name, sep = "\t", fill=TRUE, header = TRUE)[,-c(1)]
        maf_data_filtered = maf_data %>% filter(maf >= maf_cutoff)
        maf_data_snps = gsub('snp', '', maf_data_filtered$snp)
        list_of_snps = as.numeric(intersect(colnames(genotype_list), maf_data_snps))
    }
    
    print('All data compiled. Starting regressions.')

    count = 0
    together_regressions = data.table()

    for (race_index in unique(together$race)){
        results = NULL
        if (nrow(together[race == race_index])>2){
            print(paste0('processing for ', race_index))
            registerDoParallel(cores=ncpus)
            results = foreach(snp=list_of_snps, .combine='rbind') %dopar% {
                RcppParallel::setThreadOptions(1L) 
                execute_regression(snp, list_of_snps, snp_list, genotype_list = genotype_list[row.names(genotype_list) %in% together[race == race_index]$localID,], trim_type, condensing, trimming_data = together[race == race_index], repetitions,pca_structure_correction, pca_type, boot_cutoff, weighting, gene_conditioning, random_effects, regression_dataframe)
            }
            stopImplicitCluster()
            results$race = race_index
            together_regressions = rbind(together_regressions, results)
        } else {
            print(paste0('skipping ', race_index))
        }
    }
    return(together_regressions)
}


