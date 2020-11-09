library(data.table)
library(ggplot2)
library(RColorBrewer)
library(GWASTools)
library(tidyverse)
#omp_set_num_threads(1)
#blas_set_num_threads(1)
setDTthreads(threads = 1)

source('/home/mrussel2/tcr-gwas/trimming_regression/scripts/regression_functions.R')
source('/home/mrussel2/tcr-gwas/trimming_regression/scripts/compile_data_functions.R')


find_files <- function(trim_type, random_effects, condensing, d_infer, repetitions, pca_structure_correction){
    #make arbitrary snp list
    snp_list = data.frame(snp = seq(1,100))
    file_name = make_regression_file_name(snp_list, trim_type, condensing, random_effects, d_infer, repetitions, pca_structure_correction)
    path_name = paste0('/', paste(strsplit(file_name, '/')[[1]][-c(1,12)], collapse = '/'))
    file_pattern = paste0('*_', paste(strsplit(strsplit(file_name, '/')[[1]][12], '_')[[1]][-c(1:3)], collapse = '_'))
        
    data_files = list.files(path=path_name, pattern=file_pattern, full.names=TRUE)   
    return(data_files)
} 


# This script compiles all regression data from cluster across the entire genome
compile_all_data_from_cluster_sequential <- function(trim_type, pca_structure_correction){
    bonferroni = 0.05/35481497

    # parse type (either 'insert' or 'trim')
    stopifnot(trim_type %in% c('v_trim', 'd0_trim', 'd1_trim', 'j_trim', 'vd_insert', 'dj_insert', 'vj_insert'))
    type = strsplit(trim_type, '_')[[1]][2]
 
    # Set regression parameters
    condensing = ifelse(type == 'insert', 'by_patient', 'by_gene')
    random_effects = ifelse(type == 'insert', 'False', 'True')
    d_infer = ifelse(type == 'insert', 'False', 'True')
    repetitions = ifelse(type == 'insert', 0, 100)
 
    data_files = find_files(trim_type, random_effects, condensing, d_infer, repetitions, pca_structure_correction)
  
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
            assign(paste0('together_list_', trim_type), rbindlist(list(get(paste0('together_list_', trim_type)), temp_file)))
        }
        count = count + 1
        print(paste0(count, ' of ', length(data_files), ' completed for ', trim_type))
    }
    assign('together', get(paste0('together_list_', trim_type)))
    subset_and_write_df(together, trim_type, random_effects, condensing, d_infer, repetitions, pca_structure_correction, ncpus = 1)
}

# This script compiles all regression data from cluster across the entire genome
compile_all_data_from_cluster_for_parallel <- function(file, data_files, trim_type, random_effects, condensing, d_infer, repetitions, pca_structure_correction){
    index = which(file == data_files)

    # read file...
    if (file.size(file) > 1){
        temp_file = fread(file, sep = "\t", fill=TRUE, header = TRUE)
        print(paste0('finished ', index, ' of ', length(data_files)))
    
        return(temp_file)
    }
}


subset_and_write_df <- function(together, trim_type, random_effects, condensing, d_infer, repetitions, pca_structure_correction, ncpus){

    setDTthreads(threads = ncpus)

    assign(paste0('together_productive_', trim_type), together[productivity == 'productive'])
    print('finished productive')
    assign(paste0('together_NOT_productive_', trim_type), together[productivity == 'NOT_productive'])
    print('finished not productive')
    print(paste0('processed data for ', trim_type))

    together_productive = get(paste0('together_productive_', trim_type))
    together_NOT_productive = get(paste0('together_NOT_productive_', trim_type))

    if (random_effects == 'True'){
         file_name = paste0(trim_type, '_snps_regression_with_weighting_condensing_', condensing,'_with_random_effects_', repetitions, '_bootstraps')
    } else {
         file_name = paste0(trim_type, '_snps_regression_with_weighting_condensing_', condensing, '_NO_random_effects_', repetitions, '_bootstraps')
    }

    if (d_infer == 'False'){
        file_name = paste0(file_name, '_NO_d_infer')
    } 

    if (pca_structure_correction == 'True'){
        file_name = paste0(file_name, '_WITH_pca_structure_correction.tsv')
    } else{
        file_name = paste0(file_name, '.tsv')
    }
        

    write.table(together_productive, file=paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/results/productive_', file_name), quote=FALSE, sep='\t', col.names = NA)
    write.table(together_NOT_productive, file=paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/results/NOT_productive_', file_name), quote=FALSE, sep='\t', col.names = NA)
}

