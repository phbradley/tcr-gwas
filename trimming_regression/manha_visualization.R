##library(plotly)
##library(manhattanly)
library(data.table)
library(ggplot2)


compile_data_manhattan <- function(snp_meta_data, productivity, trim_type){
    bonferroni = 0.05/35481497
    if (nrow(snp_meta_data) == 11267){
        snp_meta_data = snp_meta_data[,mean(mwu_pval), by = .(snpnum, chromosome, hg19_pos)]
        snp_meta_data$snpnum = paste0("snp", snp_meta_data$snpnum)
        colnames(snp_meta_data) = c("snp", "chr", "hg19_pos", "mwu_P")
    } else {
        snp_meta_data = snp_meta_data[,1:3]
        snp_meta_data$snpid = paste0("snp", snp_meta_data$snpid)
        colnames(snp_meta_data) = c("snp", "hg19_pos", "chr")
    }
    regression_snp_list = as.data.table(read.table(paste0("regression_bootstrap_results/", productivity, "/", trim_type, "/", trim_type, "_", productivity, "_snplist_1391_snps.tsv"), header = TRUE))

    together = merge(regression_snp_list, snp_meta_data, all.x = TRUE)

    productivities = c("productive", "NOT_productive")
    prod = productivities[!productivities %in% c(productivity)]
    regression_snp_list_other = as.data.table(read.table(paste0("regression_bootstrap_results/", prod, "/", trim_type, "/", trim_type, "_", prod, "_snplist_1391_snps.tsv"), header = TRUE))[,c(1, 6)]
    colnames(regression_snp_list_other) = c("snp", "P_other")

    together_other = merge(together, regression_snp_list_other, all.x = TRUE)
    

    if (ncol(snp_meta_data) ==4){
        colnames(together) = c("SNP", "SE", "INTERCEPT", "SLOPE", "ZSCORE", "P", "CHR", "BP", "mwu_P")
        colnames(together_other) = c("SNP", "SE", "INTERCEPT", "SLOPE", "ZSCORE", "P", "CHR", "BP", "mwu_P", "P_other")
        together_other = together_other[P_other < 0.05]
    } else {
        colnames(together) = c("SNP", "SE", "INTERCEPT", "SLOPE", "ZSCORE", "P", "CHR", "BP")
        colnames(together_other) = c("SNP", "SE", "INTERCEPT", "SLOPE", "ZSCORE", "P", "CHR", "BP", "P_other")
    }
    
    par(mar=c(5,6,4,4)+.1)
    plot(together$BP, -1*log(together$P, base =10), col = alpha("black", 0.25), pch=19, cex = 1.25, xlab = 'Chromosome 7 position', ylab = '-log10(p value)', main = paste0('P-value of SNP effect on ', trim_type, ' for ', productivity, " TCRs"), panel.first = grid(), cex.main=1.5, cex.lab=1.5, cex.axis=1)
    points(together_other$BP, -1*log(together_other$P, base =10), col = alpha("blue", 0.35), pch=19, cex = 1.25,)
    abline(h = -1*log(bonferroni, base = 10), col = "red", lwd = 4, lty = 2)
    legend("topleft", box.lty=0, legend=c("-log10(0.05)", paste0("significant for ", prod, " TCRs")), col=c("red", alpha("blue", 0.35)), lty=c(2, NA), lwd = c(3, NA), pch = c(NA, 19), cex = 1.25)
}