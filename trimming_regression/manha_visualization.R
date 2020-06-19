library(plotly)
library(manhattanly)
library(data.table)


chr7_data = as.data.table(read.table("../_ignore/tcrb_vgene_oof_hg19_pvals.tsv", header = TRUE))
chr7_j_trim_regression_productive = read.table("regression_bootstrap_results/productive/j_gene/j_gene_productive_snplist_1391_snps.tsv", header = TRUE)

together = merge(chr7_j_trim_regression_productive, chr7_snp_meta_data, all.x = TRUE)

compile_data_manhattan <- function(snp_meta_data, regression_snp_list){
    if (snp_meta_data == "chr7_data"){
        snp_meta_data = snp_meta_data[,.N, by = .(snpnum, chromosome, hg19_pos)]
        snp_meta_data$snpnum = paste0("snp", snp_meta_data$snpnum)
        colnames(snp_meta_data) = c("snp", "chr", "hg19_pos", "count")
        snp_meta_data = snp_meta_data[,1:3]
    } else {
        snp_meta_data = snp_meta_data[,1:3]
        snp_meta_data$snpid = paste0("snp", snp_meta_data$snpid)
        colnames(snp_meta_data) = c("snp", "hg19_pos", "chr")
    }
    together = merge(regression_snp_list, snp_meta_data, all.x = TRUE, by = snp)
    ## NEED TO MAKE NAMING COMPATIBLE WITH MANHATTANLY!!!!!
    manhattanly(together, snp = "snp")
}