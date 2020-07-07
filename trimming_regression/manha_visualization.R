##library(plotly)
##library(manhattanly)
library(data.table)
library(ggplot2)
library(RColorBrewer)


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
        # define dataframe containing correlated snps and meta data
        together_correlated = merge(regression_snp_list_correlated, snp_meta_data, by = "snp", all.x = TRUE)
    } else {
        regression_snp_list_correlated = data.table()
        regression_snp_list_no_correlated = regression_snp_list
        correlated_snps = paste0("No correlated snps present!")
        # define dataframe containing correlated snps and meta data
        together_correlated = data.table()
    }

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
    if (nrow(together_correlated) > 0){
        together_correlated = together_correlated[chr == chromosome]
    }
    together_not_correlated = together_not_correlated[chr == chromosome]
    together_other = together_other[chr == chromosome]

    par(mar=c(5,6,4,4)+.1)
    if (nrow(together_not_correlated) > 0){
        if (nrow(together_correlated) > 0){
            xlimit = c(min(together_not_correlated$hg19_pos, together_correlated$hg19_pos)-1000, max(together_not_correlated$hg19_pos, together_correlated$hg19_pos)+1000)
            ylimit = c(0, max(-1*log(together_not_correlated$pvalue, base =10), -1*log(together_correlated$pvalue, base =10))+5)
        } else {
            xlimit = c(min(together_not_correlated$hg19_pos)-1000, max(together_not_correlated$hg19_pos)+1000)
            ylimit = c(0, max(-1*log(together_not_correlated$pvalue, base =10))+5)
        }
        plot(as.numeric(together_not_correlated$hg19_pos), as.numeric(-1*log(together_not_correlated$pvalue, base =10)), bg = alpha("black", 0.4), col = alpha("black", 0.9), xlab = paste0('Chromosome ', chromosome,' position'), ylab = '-log10(p value)', main = paste0('P-value of SNP effect on ', trim_type, ' for ', productivity, " TCRs \n on Chromosome ", chromosome, "\n", correlated_snps), panel.first = grid(), cex.main=1.5, cex.lab=1.5, cex.axis=1, pch = 21, cex = 1.5, xlim = xlimit, ylim = ylimit)

        palette(brewer.pal(n = length(unique(together_correlated$index)), name = 'Set3'))
        col = setNames(palette(), levels(together_correlated$index))
        if (nrow(together_correlated) > 0){
            points(together_correlated$hg19_pos, -1*log(together_correlated$pvalue, base =10), bg = c(factor(together_correlated$index), alpha = 0.4), col = alpha("black", 0.9), pch=21, cex = 1.5)
        }
        points(as.numeric(together_other$hg19_pos), -1*log(together_other$pvalue, base =10), col = alpha("red", 0.9), pch=1, cex = 1.5)
    } else {
        if (nrow(together_correlated) > 0){
            xlimit = c(min(together_correlated$hg19_pos)-1000, max(together_correlated$hg19_pos)+1000)
            ylimit = c(0, max(-1*log(together_correlated$pvalue, base =10))+5)
        }
        palette(brewer.pal(n = length(unique(together_correlated$index)), name = 'Set3'))
        col = setNames(palette(), levels(together_correlated$index))

        plot(together_correlated$hg19_pos, -1*log(together_correlated$pvalue, base =10), bg = c(factor(together_correlated$index), alpha = 0.4), col = alpha("black", 0.9), xlab = paste0('Chromosome ', chromosome,' position'), ylab = '-log10(p value)', main = paste0('P-value of SNP effect on ', trim_type, ' for ', productivity, " TCRs \n on Chromosome ", chromosome, "\n", correlated_snps), panel.first = grid(), cex.main=1.5, cex.lab=1.5, cex.axis=1, pch = 21, cex = 1.5, xlim = c(min(together_not_correlated$hg19_pos, together_correlated$hg19_pos)-1000, max(together_not_correlated$hg19_pos, together_correlated$hg19_pos)+1000), ylim = c(0, max(-1*log(together_not_correlated$pvalue, base =10), -1*log(together_correlated$pvalue, base =10))+5))

        points(as.numeric(together_other$hg19_pos), -1*log(together_other$pvalue, base =10), col = alpha("red", 0.9), pch=1, cex = 1.5)
    }

    abline(h = -1*log(bonferroni, base = 10), col = "red", lwd = 4, lty = 2)
    legend("topleft", box.lty=0, legend=c("-log10(bonferroni)", paste0("significant for ", prod, " TCRs"), "not in a cluster"), col=c("red", alpha("red", 0.9), alpha("black", 0.4)), lty=c(2, NA, NA), lwd = c(3, NA, NA), pch = c(NA, 1, 19), cex = 1.5)
}


compare_phil <- function(phil_p_vals, snp_id_list, productivity, trim_type){
    regression_snp_list = as.data.table(read.table(paste0("regression_bootstrap_results/", productivity, "/", trim_type, "/", trim_type, "_", productivity, "_snplist_", snp_id_list[1], "-",snp_id_list[length(snp_id_list)], "_snps_varying_intercepts_by_subject.tsv"), header = TRUE))

    # Filter phil's p-values to include only those that correspond to the same trimming region as this analysis
    if (productivity == "productive"){
        productivity_frame = "IF"
    } else if (productivity == "NOT_productive"){
        productivity_frame = "OF"
    }
    phil_p_vals$snpid = paste0("snp", phil_p_vals$snpid)
    trimming_feature = paste0(trim_type, "_", productivity_frame)
    phil_p_vals_filtered = phil_p_vals[feature == trimming_feature][,-c(1:3)]
    phil_meta_data = unique(phil_p_vals[,-c(2,4)])
    together = merge(regression_snp_list, phil_p_vals_filtered, by.x = "snp", by.y = "snpid", all.x = TRUE)
    together_meta = merge(together, phil_meta_data, by.x = "snp", by.y = "snpid", all.x = TRUE)

    
    together_meta[is.na(linreg_pval)]$linreg_pval<-1
    together_meta[is.na(pvalue)]$pvalue<-1

    plot(-1*log(together_meta$linreg_pval, base =10), -1*log(together_meta$pvalue, base =10), xlim = c(0, 60), ylim = c(0, 60), col = alpha("black", 0.35), pch=19, cex = 1.5, ylab = 'P values from random effects model', xlab = 'P values from Phils Analysis', main = paste0('P value comparison for ', trim_type, ' for ', productivity, " TCRs"), cex.main=1.5, cex.lab=1.5, cex.axis=1, panel.first = grid())
    abline(a = 0, b = 1, lty = 2, lwd = 4)
}