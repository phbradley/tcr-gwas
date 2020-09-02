##library(plotly)
##library(manhattanly)
library(data.table)
library(ggplot2)
library(RColorBrewer)
setDTthreads(threads = 1, restore_after_fork=FALSE)

source('lmer_trimming_regression_functions.R')


compile_data_manhattan <- function(snp_meta_data, snp_id_list, correlated_snps, paired_productive_snps, productivity, varying_int, trim_type, chromosome, correlate_snps, gene_conditioning, weighting, gene_region){
    gene_region = gene_region
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
        snp_meta_data$snp = paste0("snp", snp_meta_data$snp)
        colnames(snp_meta_data) = c("snp", "chr", "hg19_pos")
    }

    if (productivity == 'productive'){
        productivity_status = "True"
    } else if (productivity == 'NOT_productive'){
        productivity_status = "False"
    }

    # load regression / pvalue results
    if (varying_int == "True"){
        assign('regression_snp_list', as.data.table(read.table(generate_file_name(snp_id_list, trim_type, gene_type = paste0(substr(trim_type, 1, 1), '_gene'), productivity = productivity_status, gene_conditioning = gene_conditioning, weighting = weighting, condensing = 'by_gene', repetitions = 100), header = TRUE)))
        gene_cond = ifelse(gene_conditioning == 'True', '_gene_conditioning', '')
        weight_cond = ifelse(weighting == 'True', '_weight', '')
        regression_type = paste0('lmer', gene_cond, weight_cond)
    } else {
        regression_snp_list = as.data.table(read.table(paste0("regression_bootstrap_results/", productivity, "/", trim_type, "/", trim_type, "_", productivity, "_snplist_", snp_id_list[1], "-",snp_id_list[length(snp_id_list)], "_snps_simple_condensing_by_patient.tsv"), header = TRUE))
        regression_type = 'Simple regression: no conditioning on gene'
    }

    if (correlate_snps == 'True'){
        regression_snp_list_correlated = correlated_snps
        regression_snp_list_correlated$snp = paste0('snp', regression_snp_list_correlated$snp)
        regression_snp_list_correlated = merge(regression_snp_list_correlated, regression_snp_list, all.x = TRUE)
        regression_snp_list_no_correlated = subset(regression_snp_list, !(snp %in% regression_snp_list_correlated$snp))
        correlated_snps = paste0("Correlated snps shown by color")
        # define dataframe containing correlated snps and meta data
        together_correlated = merge(regression_snp_list_correlated, snp_meta_data, by = "snp", all.x = TRUE)
    } else {
        regression_snp_list_correlated = data.table()
        regression_snp_list_no_correlated = regression_snp_list
        correlated_snps = paste0("Correlated snps not shown")
        # define dataframe containing correlated snps and meta data
        together_correlated = data.table()
    }

    # define dataframe containing NOT correlated snps and meta data
    together_not_correlated  = merge(regression_snp_list_no_correlated, snp_meta_data, by = "snp", all.x = TRUE)

    # find snps from other productivity group
    productivities = c("productive", "NOT_productive")
    prod = productivities[!productivities %in% c(productivity)]
    # load regression / pvalue results
    if (varying_int == "True"){
        product = ifelse(prod == 'productive', 'True', 'False')
        assign('regression_snp_list_other', as.data.table(read.table(generate_file_name(snp_id_list, trim_type, gene_type = paste0(substr(trim_type, 1, 1), '_gene'), productivity = product, gene_conditioning, weighting, condensing = 'by_gene', repetitions = 100), header = TRUE))[,c(1, 5)])
    } else {
        regression_snp_list_other = as.data.table(read.table(paste0("regression_bootstrap_results/", prod, "/", trim_type, "/", trim_type, "_", prod, "_snplist_", snp_id_list[1], "-",snp_id_list[length(snp_id_list)], "_snps_simple_condensing_by_patient.tsv"), header = TRUE))[,c(1, 5)]
    }

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

    par(mar=c(5,5,7,5)+.1)
    if (nrow(together_not_correlated) > 0){
        if (nrow(together_correlated) > 0){
            xlimit = c(min(together_not_correlated$hg19_pos, together_correlated$hg19_pos)-1000, max(together_not_correlated$hg19_pos, together_correlated$hg19_pos)+1000)
            ylimit = c(0, max(-1*log(together_not_correlated$pvalue, base =10), -1*log(together_correlated$pvalue, base =10))+5)
        } else {
            xlimit = c(min(together_not_correlated$hg19_pos)-1000, max(together_not_correlated$hg19_pos)+1000)
            ylimit = c(0, max(-1*log(together_not_correlated$pvalue, base =10))+5)
        }
        plot(as.numeric(together_not_correlated$hg19_pos), as.numeric(-1*log(together_not_correlated$pvalue, base =10)), bg = alpha("black", 0.4), col = alpha("black", 0.9), xlab = paste0('Chromosome ', chromosome,' position'), ylab = '-log10(p value)', main = paste0('P-value of SNP effect on ', trim_type, ' for ', productivity, " TCRs on chr", chromosome, "\n", gene_region), panel.first = grid(), cex.main=1.75, cex.lab=1.75, cex.axis=1.5, pch = 21, cex = 1.5, xlim = xlimit, ylim = ylimit)

        palette(brewer.pal(n = length(unique(together_correlated$cluster)), name = 'Set2'))
        col = setNames(palette(), levels(as.factor(together_correlated$cluster)))
        if (nrow(together_correlated) > 0){
            points(together_correlated$hg19_pos, -1*log(together_correlated$pvalue, base =10), bg = alpha(col[together_correlated$cluster], 0.8), col = alpha("black", 0.9), pch=21, cex = 1.5)
        }
        if (paired_productive_snps == 'True'){
            points(as.numeric(together_other$hg19_pos), -1*log(together_other$pvalue, base =10), col = alpha("red", 0.9), pch=1, cex = 1.5)
        }
    } else {
        if (nrow(together_correlated) > 0){
            xlimit = c(min(together_correlated$hg19_pos)-1000, max(together_correlated$hg19_pos)+1000)
            ylimit = c(0, max(-1*log(together_correlated$pvalue, base =10))+5)
        }
        palette(brewer.pal(n = length(unique(together_correlated$cluster)), name = 'Set2'))
        col = setNames(palette(), levels(as.factor(together_correlated$cluster)))

        plot(together_correlated$hg19_pos, -1*log(together_correlated$pvalue, base =10), bg = alpha(col[together_correlated$cluster], 0.8), col = alpha("black", 0.9), xlab = paste0('Chromosome ', chromosome,' position'), ylab = '-log10(p value)', main = paste0('P-value of SNP effect on ', trim_type, ' for ', productivity, " TCRs on ch", chromosome, "\n", gene_region), panel.first = grid(), cex.main=1.75, cex.lab=1.75, cex.axis=1.5, pch = 21, cex = 1.5, xlim = c(min(together_not_correlated$hg19_pos, together_correlated$hg19_pos)-1000, max(together_not_correlated$hg19_pos, together_correlated$hg19_pos)+1000), ylim = c(0, max(-1*log(together_not_correlated$pvalue, base =10), -1*log(together_correlated$pvalue, base =10))+5))
        
        if (paired_productive_snps == 'True'){
            points(as.numeric(together_other$hg19_pos), -1*log(together_other$pvalue, base =10), col = alpha("red", 0.9), pch=1, cex = 1.5)
        }
    }

    abline(h = -1*log(bonferroni, base = 10), col = "chartreuse3", lwd = 4, lty = 2)
    if (correlate_snps == 'True'){
        legend("topleft", box.lty=0, legend=c("-log10(bonferroni)", paste0("significant for ", prod, " TCRs"), "not in a cluster", rep("cluster", length(unique(together_correlated$cluster)))), col=c("chartreuse3", alpha("red", 0.9), alpha("black", 0.4), alpha(col, 0.8)), lty=c(2, NA, NA, rep(NA, length(unique(together_correlated$cluster)))), lwd = c(3, NA, NA,rep(NA, length(unique(together_correlated$cluster)))), pch = c(NA, 1, 19,rep(19, length(unique(together_correlated$cluster)))), cex = 1.5)
    } else if (correlate_snps == 'True'){
        legend("topleft", box.lty=0, legend=c("-log10(bonferroni)", paste0("significant for ", prod, " TCRs")), col=c("chartreuse3", alpha("red", 0.9)), lty=c(2, NA), lwd = c(3, NA), pch = c(NA, 1), cex = 1.5)
    } 
}

combine_trimming_regression_dfs <- function(pre_regression_snp_list, repetitions=100, condensing, trimming_type_list){
    pre_regression_snp_list_final = pre_regression_snp_list
    pre_regression_snp_list_final$snpid = paste0("snp", pre_regression_snp_list_final$snpid)
    snp_id_list = unique(pre_regression_snp_list$snpid)
    post_regression_snp_list = data.table()
    for (trim in trimming_type_list){
        productive_test =  as.data.table(read.table(paste0('regression_bootstrap_results/productive/', trim, '/', trim, '_productive_snplist_', snp_id_list[1], "-",snp_id_list[length(snp_id_list)], '_snps_lmer_with_gene_with_weighting_condensing_', condensing, '_', repetitions, '_bootstraps.tsv') , sep = "\t", fill=TRUE, header = TRUE))
    
        not_productive_test =  as.data.table(read.table(paste0('regression_bootstrap_results/NOT_productive/', trim, '/', trim, '_NOT_productive_snplist_', snp_id_list[1], "-",snp_id_list[length(snp_id_list)],'_snps_lmer_with_gene_with_weighting_condensing_', condensing, '_', repetitions, '_bootstraps.tsv'), sep = "\t", fill=TRUE, header = TRUE))
        

        if (nrow(productive_test) != 0){
            post_regression_snp_list = rbind(merge(pre_regression_snp_list_final[feature == paste0(trim, "_IF")], productive_test, by.x = "snpid", by.y = "snp"), post_regression_snp_list)
        } 
        
        if (nrow(not_productive_test) != 0){
            post_regression_snp_list = rbind(merge(pre_regression_snp_list_final[feature == paste0(trim, "_OF")], not_productive_test, by.x = "snpid", by.y = "snp"), post_regression_snp_list)
        }
    }
    return(post_regression_snp_list)
}



compare_phil <- function(phil_p_vals, snp_id_list, productivity, varying_int, trim_type, gene_conditioning, weighting){
    # load regression / pvalue results
    if (varying_int == "True"){
        regression_snp_list = as.data.table(read.table(generate_file_name(snp_id_list, trim_type, productivity = 'True', gene_conditioning, weighting, condensing, repetitions = 100), header = TRUE))
        regression_type = 'lmer regression: conditioning on gene'
    } else {
        regression_snp_list = as.data.table(read.table(paste0("regression_bootstrap_results/", productivity, "/", trim_type, "/", trim_type, "_", productivity, "_snplist_", snp_id_list[1], "-",snp_id_list[length(snp_id_list)], "_snps_simple.tsv"), header = TRUE))
        regression_type = 'Simple regression: no conditioning on gene'
    }
    
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

    plot(-1*log(together_meta$linreg_pval, base =10), -1*log(together_meta$pvalue, base =10), xlim = c(0, 60), ylim = c(0, 60), col = alpha("black", 0.35), pch=19, cex = 1.5, ylab = paste0('P values from ', regression_type), xlab = 'P values from Phils Analysis', main = paste0('P value comparison for ', trim_type, ' for ', productivity, " TCRs"), cex.main=1.5, cex.lab=1.5, cex.axis=1, panel.first = grid())
    abline(a = 0, b = 1, lty = 2, lwd = 4)
}

compare_phil_regression_list <- function(regression_list){
    
    regression_list[is.na(linreg_pval)]$linreg_pval<-1
    regression_list[is.na(pvalue)]$pvalue<-1

    palette(brewer.pal(n = length(unique(regression_list$feature)), name = 'Set2'))
    col = setNames(palette(), levels(as.factor(regression_list$feature)))

    plot(-1*log(regression_list$linreg_pval, base =10), -1*log(regression_list$pvalue, base =10), xlim = c(0, 55), ylim = c(0, 55), col = as.factor(regression_list$feature), pch=19, cex = 1.5, ylab = paste0('P values from lmer: gene conditioning, weighting'), xlab = 'P values from Phils Analysis', main = paste0('P value comparison'), cex.main=1.5, cex.lab=1.5, cex.axis=1, panel.first = grid())
    abline(a = 0, b = 1, lty = 2, lwd = 4)

    legend("topleft", box.lty=0, legend=levels(as.factor(regression_list$feature)), col=col, pch = c(rep(19, length(unique(regression_list$feature)))), cex = 1.5)
}

visualize_crosses <- function(cross_regression_list){
    cross_regression_list = unique(cross_regression_list)
    reshaped = data.table()
    for (trim in c('v_trim', 'd0_trim', 'd1_trim', 'j_trim')){
        for (gene in c('v_gene', 'j_gene', 'd_gene')){
            for (prod in c('productive', 'NOT_productive')){
                simple_regression = as.data.table(read.table(paste0("regression_bootstrap_results/", prod, "/", trim, "/", trim, '_',prod,  '_snplist_16814629-16814111_snps_simple_with_weighting_condensing_by_patient_0_bootstraps.tsv')), header = TRUE)[,c('snp', 'pvalue')]
                desired_cols = c('snp', 'productivity_status', paste0('pvalue_', trim, '_', gene))
                assign(paste0(trim, '_', gene,'_dt'), cross_regression_list[productivity_status == prod,..desired_cols])
                new_name_column = rep(paste0(trim, '_', gene), nrow(get(paste0(trim, '_', gene,'_dt'))))
                new_df = cbind(get(paste0(trim, '_', gene,'_dt')), type = new_name_column)
                new_df_with_simple = merge(new_df, simple_regression)
                colnames(new_df_with_simple) = c('snp', 'productivity_status', 'gene_conditioning_regression', 'type', 'simple_regression')
                reshaped = rbind(reshaped, new_df_with_simple)
            }   
        }
    }

    palette(brewer.pal(n = length(unique(reshaped$type)), name = 'Set3'))
    col = setNames(palette(), levels(as.factor(reshaped$type)))

    plot(-1*log(reshaped$simple_regression, base =10), -1*log(reshaped$gene_conditioning_regression, base =10), col = as.factor(reshaped$type), pch=19, cex = 1.5, ylab = paste0('P values from gene conditioning regressions'), xlab = 'P values from simple regressions', main = paste0('P value comparison'), cex.main=1.5, cex.lab=1.5, cex.axis=1, panel.first = grid())
    points(-1*log(reshaped[snp == 'snp20717772' | snp == 'snp20717781']$simple_regression, base= 10), -1*log(reshaped[snp == 'snp20717772' | snp == 'snp20717781']$gene_conditioning_regression, base =10), col = alpha("red", 0.9), pch=1, cex = 1.5)
    abline(a = 0, b = 1, lty = 2, lwd = 4)
    legend("bottomright", box.lty=0, legend=levels(as.factor(reshaped$type)), col=col, pch = c(rep(19, length(unique(reshaped$type)))), cex = 1.1)
}