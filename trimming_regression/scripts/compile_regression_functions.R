library(data.table)
library(ggplot2)
library(RColorBrewer)
library(GWASTools)
library(tidyverse)
#omp_set_num_threads(1)
#blas_set_num_threads(1)
setDTthreads(threads = 1)

source(paste0(PROJECT_PATH, "/tcr-gwas/trimming_regression/scripts/config.R"))
source(paste0(PROJECT_PATH, '/tcr-gwas/trimming_regression/scripts/regression_functions.R'))
source(paste0(PROJECT_PATH, '/tcr-gwas/trimming_regression/scripts/compile_data_functions.R'))


find_files <- function(trim_type, pca_structure_correction){
    #make arbitrary snp list
    snp_list = data.frame(snp = seq(1,100))
    file_name = make_regression_file_name(snp_list, trim_type, pca_structure_correction)
    path_name = paste0('/', paste(strsplit(file_name, '/')[[1]][-c(1,12)], collapse = '/'))
    file_pattern = paste0('*_', paste(strsplit(strsplit(file_name, '/')[[1]][12], '_')[[1]][-c(1:3)], collapse = '_'))
        
    data_files = list.files(path=path_name, pattern=file_pattern, full.names=TRUE)   
    return(data_files)
} 


# This script compiles all regression data from cluster across the entire genome
compile_all_data_from_cluster_sequential <- function(trim_type, pca_structure_correction, pca_type){
    bonferroni = 0.05/35481497

    # parse type (either 'insert' or 'trim')
    stopifnot(trim_type %in% c('v_trim', 'd0_trim', 'd1_trim', 'j_trim', 'vd_insert', 'dj_insert', 'vj_insert'))
    set_regression_parameters(trim_type)
 
    data_files = find_files(trim_type, pca_structure_correction)
  
    print(paste0('compile file list for ', trim_type))
    assign(paste0('together_list_', trim_type), NULL)
    count = 0
    for (file in data_files){
        # read file...
        if (file.size(file) == 1 | file.size(file) == 0){
            next
        }
        temp_file = fread(file, sep = "\t", fill=TRUE, header = TRUE)
        #together = rbind(together, temp_file)
        if (ncol(temp_file) > 2){
            assign(paste0('together_list_', trim_type), 
                   rbindlist(list(get(paste0('together_list_', trim_type)), temp_file)))
        }
        count = count + 1
        print(paste0(count, ' of ', length(data_files), ' completed for ', trim_type))
    }
    assign('together', get(paste0('together_list_', trim_type)))
    subset_and_write_df(together, trim_type, pca_structure_correction, pca_type, ncpus = 1)
}


make_compiled_regression_file_name <- function(productivity, trim_type, pca_structure_correction, pca_type){
    random_effects_name = ifelse(RANDOM_EFFECTS == 'True', '_random_effects', '_no_random_effects')
    d_infer_name = ifelse(D_INFER == 'True', '_d_infer', '_no_d_infer')
    pca_name = ifelse(pca_structure_correction == 'True', '_pca_correction', '_no_pca_correction')
    boots = paste0('_', REPETITIONS, '_bootstraps')
    weight_name = ifelse(WEIGHTING=='True', '_with_weighting', '')

    file_name = paste0(OUTPUT_PATH, 
                       '/results/', 
                       productivity, '_', trim_type, '_snps_regression', weight_name, '_condensing_', CONDENSING, random_effects_name, d_infer_name, pca_name, boots,'_', pca_type, '.tsv')
    if (!dir.exists(paste0(OUTPUT_PATH, '/results'))){
        system(paste0('mkdir ', OUTPUT_PATH, '/results'))
    }
    return(file_name)
}



subset_and_write_df <- function(together, trim_type, pca_structure_correction, pca_type, ncpus){

    setDTthreads(threads = ncpus)

    assign(paste0('together_productive_', trim_type), together[productivity == 'productive'])
    print('finished productive')
    assign(paste0('together_NOT_productive_', trim_type), together[productivity == 'NOT_productive'])
    print('finished not productive')
    print(paste0('processed data for ', trim_type))

    together_productive = get(paste0('together_productive_', trim_type))
    together_NOT_productive = get(paste0('together_NOT_productive_', trim_type))

    write.table(together_productive, file = make_compiled_regression_file_name(productivity = 'productive', 
                                                                               trim_type, 
                                                                               pca_structure_correction, 
                                                                               pca_type) , quote=FALSE, sep='\t', col.names = NA)
    write.table(together_NOT_productive, file = make_compiled_regression_file_name(productivity = 'NOT_productive', 
                                                                                   trim_type, 
                                                                                   pca_structure_correction, 
                                                                                   pca_type), quote=FALSE, sep='\t', col.names = NA)
}

