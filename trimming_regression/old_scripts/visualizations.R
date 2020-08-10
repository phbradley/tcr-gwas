compile_data_manhattan <- function(snp_meta_data, snp_id_list, productivity, trim_type, chromosome){
    bonferroni = 0.05/35481497

    # sort the snp meta data from phil (i.e. assign common column names)
    if (nrow(snp_meta_data) == 11267){
        snp_meta_data = snp_meta_data[,mean(mwu_pval), by = .(snpnum, chromosome, hg19_pos)]
        snp_meta_data$snpnum = paste0("snp", snp_meta_data$snpnum)
        colnames(snp_meta_data) = c("snp", "chr", "hg19_pos", "phil_pval")
    } else if (nrow(snp_meta_data) == 96){
        snp_meta_data$snpid = paste0("snp", snp_meta_data$snpid)
        colnames(snp_meta_data) = c("chr", "feature", "hg19_pos", "phil_pval", "snp")
    } else {
        snp_meta_data = snp_meta_data[,1:3]
        snp_meta_data$snpid = paste0("snp", snp_meta_data$snpid)
        colnames(snp_meta_data) = c("snp", "hg19_pos", "chr")
    }

    # load regression / pvalue results
    regression_snp_list = as.data.table(read.table(paste0("regression_bootstrap_results/", productivity, "/", trim_type, "/", trim_type, "_", productivity, "_snplist_", snp_id_list[1], "-",snp_id_list[length(snp_id_list)], "_snps_varying_intercepts_by_subject.tsv"), header = TRUE))

    # load correlated snps list (if it exists)
    correlated_snps_filename = paste0('regression_bootstrap_results/', productivity, '/', trim_type, '/', trim_type, '_correlated_snps_from_snplist_', snp_id_list[1], "-",snp_id_list[length(snp_id_list)], '_snps_varying_intercepts_by_subject.tsv')

    if (file.exists(correlated_snps_filename)){
        regression_snp_list_correlated = as.data.table(read.table(correlated_snps_filename, header = TRUE))[, -c(7:9)]
        colnames(regression_snp_list_correlated) = c("index", colnames(regression_snp_list))
        regression_snp_list_no_correlated = subset(regression_snp_list, !(snp %in% regression_snp_list_correlated$snp))
        correlated_snps = paste0("Correlated snps shown by color")
    } else {
        regression_snp_list_no_correlated = data.table()
        correlated_snps = paste0("No correlated snps present!")
    }

    # define dataframe containing correlated snps and meta data
    together_correlated = merge(regression_snp_list_correlated, snp_meta_data, by = "snp", all.x = TRUE)

    # define dataframe containing NOT correlated snps and meta data
    together_not_correlated  = merge(regression_snp_list_no_correlated, snp_meta_data, by = "snp", all.x = TRUE)

    # find snps from other productivity group
    productivities = c("productive", "NOT_productive")
    prod = productivities[!productivities %in% c(productivity)]
    regression_snp_list_other = as.data.table(read.table(paste0("regression_bootstrap_results/", prod, "/", trim_type, "/", trim_type, "_", prod, "_snplist_", snp_id_list[1], "-",snp_id_list[length(snp_id_list)], "_snps_varying_intercepts_by_subject.tsv"), header = TRUE))[,c(1, 5)]
    colnames(regression_snp_list_other) = c("snp", "P_other")

    # combind correlated and not correlated snps
    together_temp = rbind(together_correlated, together_not_correlated, fill=TRUE)

    # find snps significant in other productivity group
    together_other= merge(together_temp, regression_snp_list_other, all.x = TRUE)
    together_other = together_other[P_other < bonferroni]

    # Use only data from indicated chromosome
    together_correlated = together_correlated[CHR == chromosome]
    together_not_correlated = together_not_correlated[CHR == chromosome]
    together_other = together_other[CHR == chromosome]

    par(mar=c(5,6,4,4)+.1)

    plot(together_not_correlated$hg19_pos, -1*log(together_not_correlated$pvalue, base =10), col = alpha("black", 0.4), xlab = paste0('Chromosome ', chromosome,' position'), ylab = '-log10(p value)', main = paste0('P-value of SNP effect on ', trim_type, ' for ', productivity, " TCRs \n on Chromosome ", chromosome, "\n", correlated_snps), panel.first = grid(), cex.main=1.5, cex.lab=1.5, cex.axis=1, pch = 19)
}
    palette(brewer.pal(n = length(unique(together_correlated$index)), name = 'Set2'))
    col = setNames(palette(), levels(together_correlated$index))


compile_data_manhattan <- function(snp_meta_data, snp_id_list, productivity, trim_type, chromosome){
    

    together_temp = rbind(together_correlated, together_not_correlated, fill=TRUE)
    together_other= merge(together_temp, regression_snp_list_other, all.x = TRUE)

    if (ncol(snp_meta_data) ==4){
        colnames(together) = c("SNP", "SE", "INTERCEPT", "SLOPE", "ZSCORE", "P", "CHR", "BP", "phil_pval")
        colnames(together_other) = c("SNP", "SE", "INTERCEPT", "SLOPE", "ZSCORE", "P", "CHR", "BP", "phil_pval", "P_other")
        together_other = together_other[P_other < bonferroni]
    } else if (ncol(snp_meta_data) ==5){
        colnames(together_correlated) = c("SNP", "CHR", "INDEX", "INTERCEPT", "SLOPE", "SE", "P", "START", "group_P", "FEATURE", "BP", "phil_pval")
        colnames(together_not_correlated) = c("SNP", "INTERCEPT", "SLOPE", "SE", "P", "CHR", "FEATURE", "BP", "phil_pval")
        colnames(together_other) = c("SNP", "CHR", "INDEX", "INTERCEPT", "SLOPE", "SE", "P", "START", "group_P", "FEATURE", "BP", "phil_pval", "P_other")
        together_other = together_other[P_other < bonferroni]
    } else {
        colnames(together) = c("SNP", "SE", "INTERCEPT", "SLOPE", "ZSCORE", "P", "CHR", "BP")
        colnames(together_other) = c("SNP", "SE", "INTERCEPT", "SLOPE", "ZSCORE", "P", "CHR", "BP", "P_other")
    }
    
    together_correlated = together_correlated[CHR == chromosome]
    together_not_correlated = together_not_correlated[CHR == chromosome]
    together_other = together_other[CHR == chromosome] 

    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = length(unique(together_correlated$INDEX)), name = 'Set2'))
    col = setNames(palette(), levels(together_correlated$INDEX))

    plot(c(0,1,5), c(10,10,15))
    #plot(together_not_correlated$BP, -1*log(together_not_correlated$P, base =10), col = alpha("black", 0.4), xlab = paste0('Chromosome ', chromosome,' position'), ylab = '-log10(p value)', main = paste0('P-value of SNP effect on ', trim_type, ' for ', productivity, " TCRs \n on Chromosome ", chromosome, "\n", correlated_snps), panel.first = grid(), cex.main=1.5, cex.lab=1.5, cex.axis=1, pch = 19)

    #points(together_correlated$BP, -1*log(together_correlated$P, base =10), col = c(factor(together_correlated$INDEX), alpha = 0.4), pch=19, cex = 1.5)
    
    #points(as.numeric(together_other$BP), -1*log(together_other$P, base =10), col = alpha("red", 0.9), pch=1, cex = 1.5)
    abline(h = -1*log(bonferroni, base = 10), col = "red", lwd = 4, lty = 2)
    legend("topleft", box.lty=0, legend=c("-log10(bonferroni)", paste0("significant for ", prod, " TCRs"), "not in a cluster"), col=c("red", alpha("red", 0.9), alpha("black", 0.4)), lty=c(2, NA, NA), lwd = c(3, NA, NA), pch = c(NA, 1, 19), cex = 1.5)
}

compile_data_manhattan <- function(snp_meta_data, snp_id_list, productivity, trim_type, chromosome){
    bonferroni = 0.05/35481497
    if (nrow(snp_meta_data) == 11267){
        snp_meta_data = snp_meta_data[,mean(mwu_pval), by = .(snpnum, chromosome, hg19_pos)]
        snp_meta_data$snpnum = paste0("snp", snp_meta_data$snpnum)
        colnames(snp_meta_data) = c("snp", "chr", "hg19_pos", "phil_pval")
    } else if (nrow(snp_meta_data) == 96){
        snp_meta_data$snpid = paste0("snp", snp_meta_data$snpid)
        colnames(snp_meta_data) = c("chr", "feature", "hg19_pos", "phil_pval", "snp")
    } else {
        snp_meta_data = snp_meta_data[,1:3]
        snp_meta_data$snpid = paste0("snp", snp_meta_data$snpid)
        colnames(snp_meta_data) = c("snp", "hg19_pos", "chr")
    }

    regression_snp_list = as.data.table(read.table(paste0("regression_bootstrap_results/", productivity, "/", trim_type, "/", trim_type, "_", productivity, "_snplist_", snp_id_list[1], "-",snp_id_list[length(snp_id_list)], "_snps_varying_intercepts_by_subject.tsv"), header = TRUE))

    correlated_snps_filename = paste0('regression_bootstrap_results/', productivity, '/', trim_type, '/', trim_type, '_correlated_snps_from_snplist_', snp_id_list[1], "-",snp_id_list[length(snp_id_list)], '_snps_varying_intercepts_by_subject.tsv')
    
    if (file.exists(correlated_snps_filename)){
        regression_snp_list_correlated = as.data.table(read.table(correlated_snps_filename, header = TRUE))
        regression_snp_list_no_correlated = subset(regression_snp_list, !(snp %in% regression_snp_list_correlated$snp))
        correlated_snps = paste0("Correlated snps shown by color")
    } else {
        regression_snp_list_no_correlated = data.table()
        correlated_snps = paste0("No correlated snps present!")
    }
    

    together_correlated = merge(regression_snp_list_correlated, snp_meta_data, by = "snp", all.x = TRUE)

    together_not_correlated  = merge(regression_snp_list_no_correlated, snp_meta_data, by = "snp", all.x = TRUE)
    colnames(together_not_correlated) = c("snp", "intercept","slope","standard_error","p","chr","feature","hg19_pos", "phil_pval")

    productivities = c("productive", "NOT_productive")
    prod = productivities[!productivities %in% c(productivity)]
    regression_snp_list_other = as.data.table(read.table(paste0("regression_bootstrap_results/", prod, "/", trim_type, "/", trim_type, "_", prod, "_snplist_", snp_id_list[1], "-",snp_id_list[length(snp_id_list)], "_snps_varying_intercepts_by_subject.tsv"), header = TRUE))[,c(1, 5)]
    colnames(regression_snp_list_other) = c("snp", "P_other")

    together_temp = rbind(together_correlated, together_not_correlated, fill=TRUE)
    together_other= merge(together_temp, regression_snp_list_other, all.x = TRUE)

    if (ncol(snp_meta_data) ==4){
        colnames(together) = c("SNP", "SE", "INTERCEPT", "SLOPE", "ZSCORE", "P", "CHR", "BP", "phil_pval")
        colnames(together_other) = c("SNP", "SE", "INTERCEPT", "SLOPE", "ZSCORE", "P", "CHR", "BP", "phil_pval", "P_other")
        together_other = together_other[P_other < bonferroni]
    } else if (ncol(snp_meta_data) ==5){
        colnames(together_correlated) = c("SNP", "CHR", "INDEX", "INTERCEPT", "SLOPE", "SE", "P", "START", "group_P", "FEATURE", "BP", "phil_pval")
        colnames(together_not_correlated) = c("SNP", "INTERCEPT", "SLOPE", "SE", "P", "CHR", "FEATURE", "BP", "phil_pval")
        colnames(together_other) = c("SNP", "CHR", "INDEX", "INTERCEPT", "SLOPE", "SE", "P", "START", "group_P", "FEATURE", "BP", "phil_pval", "P_other")
        together_other = together_other[P_other < bonferroni]
    } else {
        colnames(together) = c("SNP", "SE", "INTERCEPT", "SLOPE", "ZSCORE", "P", "CHR", "BP")
        colnames(together_other) = c("SNP", "SE", "INTERCEPT", "SLOPE", "ZSCORE", "P", "CHR", "BP", "P_other")
    }
    
    together_correlated = together_correlated[CHR == chromosome]
    together_not_correlated = together_not_correlated[CHR == chromosome]
    together_other = together_other[CHR == chromosome] 

    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = length(unique(together_correlated$INDEX)), name = 'Set2'))
    col = setNames(palette(), levels(together_correlated$INDEX))

    plot(c(0,1,5), c(10,10,15))
    #plot(together_not_correlated$BP, -1*log(together_not_correlated$P, base =10), col = alpha("black", 0.4), xlab = paste0('Chromosome ', chromosome,' position'), ylab = '-log10(p value)', main = paste0('P-value of SNP effect on ', trim_type, ' for ', productivity, " TCRs \n on Chromosome ", chromosome, "\n", correlated_snps), panel.first = grid(), cex.main=1.5, cex.lab=1.5, cex.axis=1, pch = 19)

    #points(together_correlated$BP, -1*log(together_correlated$P, base =10), col = c(factor(together_correlated$INDEX), alpha = 0.4), pch=19, cex = 1.5)
    
    #points(as.numeric(together_other$BP), -1*log(together_other$P, base =10), col = alpha("red", 0.9), pch=1, cex = 1.5)
    abline(h = -1*log(bonferroni, base = 10), col = "red", lwd = 4, lty = 2)
    legend("topleft", box.lty=0, legend=c("-log10(bonferroni)", paste0("significant for ", prod, " TCRs"), "not in a cluster"), col=c("red", alpha("red", 0.9), alpha("black", 0.4)), lty=c(2, NA, NA), lwd = c(3, NA, NA), pch = c(NA, 1, 19), cex = 1.5)
}

