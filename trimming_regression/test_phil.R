source("run_bootstrap_regression_all_snps_functions.R")
source("manha_visualization.R")

snp_list = as.data.table(read.table("../_ignore/phil_sig_snps/top100_trimming_pvals_for_maggie.tsv", sep = "\t", fill=TRUE, header = TRUE))[-c(1:4)]

test_simple_reg <- function(snp_list, repetitions){
    snp_list_final = snp_list
    snp_list_final$snpid = paste0("snp", snp_list_final$snpid)
    results = data.table()
    for (trim_type in c('v_trim', 'd1_trim','j_trim')){
        snp_list_temp = snp_list[feature == paste0(trim_type, "_IF") | feature == paste0(trim_type, "_OF")]
        snp_id_list = unique(snp_list_temp$snpid)
        condensing = 'by_patient'

        run_snps_trimming_snp_list(snp_id_list = unique(snp_list_temp$snpid), trim_type, condensing = 'by_patient', gene_conditioning = 'false', weighting = 'false', repetitions = 5)
    
        productive_test =  as.data.table(read.table(paste0('regression_bootstrap_results/productive/', trim_type, '/', trim_type, '_productive_snplist_', snp_id_list[1], "-",snp_id_list[length(snp_id_list)], '_snps_simple_condensing_', condensing, '.tsv') , sep = "\t", fill=TRUE, header = TRUE))
    
        not_productive_test =  as.data.table(read.table(paste0('regression_bootstrap_results/NOT_productive/', trim_type, '/', trim_type, '_NOT_productive_snplist_', snp_id_list[1], "-",snp_id_list[length(snp_id_list)],'_snps_simple_condensing_', condensing, '.tsv'), sep = "\t", fill=TRUE, header = TRUE))

        if (nrow(productive_test) != 0){
            results = rbind(merge(snp_list_final[feature == paste0(trim_type, "_IF")], productive_test, by.x = "snpid", by.y = "snp"), results)
        } 
        
        if (nrow(not_productive_test) != 0){
            results = rbind(merge(snp_list_final[feature == paste0(trim_type, "_OF")], not_productive_test, by.x = "snpid", by.y = "snp"), results)
        }
    }
    return(results)
}

results = test_simple_reg(snp_list = snp_list, repetitions = 5)

write.table(results, file= 'regression_bootstrap_results/phil_pvalue_troubleshooting.tsv', quote=FALSE, sep='\t', col.names = NA)

# plot!
compare_phil_me()
