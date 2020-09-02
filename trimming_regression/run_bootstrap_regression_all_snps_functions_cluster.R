library("tidyverse")
library("gdsfmt")
library("SNPRelate")
library("GWASTools")
library("parallel")


source("simple_trimming_regression_functions.R")
source("lmer_trimming_regression_functions.R")
source("bootstrap_functions.R")
source("compile_regression_data_functions.R")
source("execute_regression_function.R")


run_snps_trimming_snp_list_cluster <- function(snp_id_list, trim_type, gene_type, condensing, gene_conditioning, weighting, random_effects, repetitions, write_table, ncpus){
    # Get snp meta data
    snps_gds = snpgdsOpen("../_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds", allow.fork=TRUE)
    
    regression_dataframe = data.frame()

    # import condensed trimming file
    if (condensing == 'by_patient'){
        assign('trimming_data', as.data.frame(read.table(paste0("../_ignore/by_patient_condensed_data_all_patients.tsv"), sep = "\t", fill=TRUE, header = TRUE)[-1]))
        names(trimming_data)[names(trimming_data) == "patient_id"] <- "localID"
    } else if (condensing == 'gene_cross'){
        trimming_data = compile_trimming_data_cross()
        trimming_data = trimming_data %>% filter(gene_class == gene_type)
        names(trimming_data)[names(trimming_data) == 'weighted_gene_count'] <- paste0('weighted_', gene_type, '_count')
        names(trimming_data)[names(trimming_data) == 'gene'] <- paste0(gene_type)
    } else if (condensing == 'phil'){
        assign('trimming_data', as.data.frame(read.table(paste0("../_ignore/by_patient_condensed_data_all_patients_phil.tsv"), sep = "\t", fill=TRUE, header = TRUE)[-1]))
    } else {
        assign('trimming_data', as.data.frame(read.table(paste0("../_ignore/condensed_", trim_type, "_data_all_patients.tsv"), sep = "\t", fill=TRUE, header = TRUE)[-1]))
        names(trimming_data)[names(trimming_data) == "patient_id"] <- "localID"
    }

    results = do.call(rbind, mclapply(snp_id_list, execute_regression, snps_gds, snp_id_list, trim_type, gene_type, condensing, trimming_data, repetitions, weighting, gene_conditioning, random_effects, regression_dataframe, mc.cores = ncpus))
    
    
    if (random_effects == 'True'){
         file_name = paste0('cluster_job_results/',trim_type, '/', trim_type, '_',snp_id_list[1], '-', snp_id_list[length(snp_id_list)],'_snps_regression_with_weighting_condensing_by_gene_with_random_effects_', repetitions, '_bootstraps.tsv')
    } else {
         file_name = paste0('cluster_job_results/',trim_type, '/', trim_type, '_',snp_id_list[1], '-', snp_id_list[length(snp_id_list)],'_snps_regression_with_weighting_condensing_by_gene_NO_random_effects_', repetitions, '_bootstraps.tsv')
    }
            
    closefn.gds(snps_gds)
    if (write_table == "True"){
        # Write tables
        write.table(as.data.frame(results), file= file_name, quote=FALSE, sep='\t', col.names = NA)
    } else if (write_table != "True"){
        return(results)
    }
}


