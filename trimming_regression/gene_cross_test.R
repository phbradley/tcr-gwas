source("run_bootstrap_regression_all_snps_functions.R")

snp_list = as.data.table(read.table("../_ignore/phil_sig_snps/top100_trimming_pvals_for_maggie.tsv", sep = "\t", fill=TRUE, header = TRUE))[-c(1:4)]


gene_cross_regression_all <- function(snp_list){
    snp_id_list = unique(snp_list$snpid)
    gene_cross_regression_dt = data.table(snp = rep(paste0('snp', snp_list$snpid), 2), productivity_status = c(rep("productive", nrow(snp_list)), rep("NOT_productive", nrow(snp_list))))
    for (trim in c('v_trim', 'd0_trim', 'd1_trim', 'j_trim')){
        for (gene in c('v_gene', 'd_gene', 'j_gene')){
            temp = run_snps_trimming_snp_list(snp_id_list = unique(snp_list$snpid), trim_type = trim, gene_type = gene, condensing = 'gene_cross', gene_conditioning = 'True', weighting = 'True', repetitions = 0, write_table = 'False')[,c('snp', 'pvalue', 'productivity_status')]
            setnames(temp, "pvalue", paste0("pvalue_",trim, '_', gene))
            gene_cross_regression_dt = merge(gene_cross_regression_dt, temp)
        }
    }
    filename = paste0('regression_bootstrap_results/crosses/', 'crossing_all_snplist_', snp_id_list[1], "-",snp_id_list[length(snp_id_list)], '_snps_lmer_with_gene_with_weighting_condensing_by_gene_0_bootstraps.tsv') 
    write.table(gene_cross_regression_dt, file= filename, quote=FALSE, sep='\t', col.names = NA)
}

gene_cross_regression_all(snp_list)