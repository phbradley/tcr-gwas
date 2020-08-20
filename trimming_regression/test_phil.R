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



compare_phil <- function(phil_p_vals, snp_id_list, productivity, varying_int, trim_type, gene_conditioning, weighting){
    # load regression / pvalue results
    if (varying_int == "True"){
        regression_snp_list = as.data.table(read.table(generate_file_name(snp_id_list, trim_type, productivity = 'True', gene_conditioning, weighting, condensing), header = TRUE))
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


compare_phil_ALL <- function(phil_p_vals, snp_id_list, productivity, varying_int, trim_type){
    # Filter phil's p-values to include only those that correspond to the same trimming region as this analysis
    if (productivity == 'productive'){
        productivity_frame = "IF"
        productivity_status = "True"
    } else if (productivity == 'NOT_productive'){
        productivity_frame = "OF"
        productivity_status = "False"
    }

    # load regression / pvalue results
    regression_snp_list_all = data.table()
    for (gene in c('True', 'False')){
        for (weight in c('True', 'False')){
            assign("data_temp", as.data.table(read.table(generate_file_name(snp_id_list, trim_type, productivity = productivity_status, gene_conditioning = gene, weighting = weight, condensing), header = TRUE, col.names = c('snp', 'intercept', 'slope', 'standard_error', 'pvalue')))[,c(1,5)])
            gene_cond = ifelse(gene == 'True', '_gene_conditioning', '')
            weight_cond = ifelse(weight == 'True', '_weight', '')
            data_temp$regression_type = paste0('lmer', gene_cond, weight_cond)
            regression_snp_list_all = rbind(regression_snp_list_all, data_temp)
        }
    }

    simple_regression_snp_list = as.data.table(read.table(paste0("regression_bootstrap_results/", productivity, "/", trim_type, "/", trim_type, "_", productivity, "_snplist_", snp_id_list[1], "-",snp_id_list[length(snp_id_list)], "_snps_simple.tsv"), header = TRUE))[,c(1,5)]
    simple_regression_snp_list$regression_type = 'glm'

    regression_snp_list_all = rbind(regression_snp_list_all, simple_regression_snp_list)

    phil_p_vals$snpid = paste0("snp", phil_p_vals$snpid)
    trimming_feature = paste0(trim_type, "_", productivity_frame)
    phil_p_vals_filtered = phil_p_vals[feature == trimming_feature][,-c(1:3)]
    phil_meta_data = unique(phil_p_vals[,-c(2,4)])

    together = merge(regression_snp_list_all, phil_p_vals_filtered, by.x = "snp", by.y = "snpid", all.x = TRUE)

    
    together[is.na(linreg_pval)]$linreg_pval<-1
    together[is.na(pvalue)]$pvalue<-1

    palette(brewer.pal(n = length(unique(together$regression_type)), name = 'Set2'))
    col = setNames(palette(), levels(as.factor(together$regression_type)))

    plot(jitter(-1*log(together$linreg_pval, base =10), 0.5), -1*log(together$pvalue, base =10), xlim = c(0, 90), ylim = c(0, 90), col = alpha(col[together$regression_type], 0.6), pch=19, cex = 1.5, ylab = paste0('P values from lmer'), xlab = 'P values from Phils Analysis', main = paste0('P value comparison for ', trim_type, ' for ', productivity, " TCRs"), cex.main=1.5, cex.lab=1.5, cex.axis=1, panel.first = grid())
    abline(a = 0, b = 1, lty = 2, lwd = 4)
    legend("bottomright", box.lty=0, legend=levels(as.factor(together$regression_type)), col= alpha(col, 0.6), pch = c(rep(19, 4)), cex = 1.5)
}

compare_phil_me <- function(){
    png(file='figures/phil_pvalue_troublshoot_plot.png', width=600, height=600)

    snps_pvals = as.data.table(read.table('regression_bootstrap_results/phil_pvalue_troubleshooting.tsv', header = TRUE))

    pvals_maggie_lm_pvalue = snps_pvals[,c(1:4, 10)]
    colnames(pvals_maggie_lm_pvalue) = c('snp', 'chromosome', 'feature', 'hg19_pos', 'maggie_lm_pvalue')

    pvals_maggie_python_pvalue = snps_pvals[,c(1:4, 11)]
    colnames(pvals_maggie_python_pvalue) = c('snp', 'chromosome', 'feature', 'hg19_pos', 'maggie_python_pvalue')

    pvals_phil_python_pvalue = snps_pvals[,c(1:5)]
    colnames(pvals_phil_python_pvalue) = c('snp', 'chromosome', 'feature', 'hg19_pos', 'phil_python_pvalue')

    snps_pvals_plotting = cbind(pvals_maggie_lm_pvalue, pvals_maggie_python_pvalue, pvals_phil_python_pvalue)

    palette(brewer.pal(n = length(unique(snps_pvals_plotting$feature)), name = 'Set2'))
    col = setNames(palette(), levels(as.factor(snps_pvals_plotting$feature)))

    plot(jitter(-1*log(snps_pvals_plotting$phil_python_pvalue, base =10), 0.5), -1*log(snps_pvals_plotting$maggie_python_pvalue, base =10), xlim = c(0, 40), ylim = c(0, 40), col = alpha(col[snps_pvals_plotting$feature], 0.6), pch=17, cex = 1.5, ylab = 'P values from Maggie Analysis', xlab = 'P values from Phil Analysis', main = paste0('P value comparison'), cex.main=1.5, cex.lab=1.5, cex.axis=1, panel.first = grid())

    points(jitter(-1*log(snps_pvals_plotting$phil_python_pvalue, base =10), 0.5), -1*log(snps_pvals_plotting$maggie_lm_pvalue, base =10), pch=15, cex = 1.5, col = alpha(col[snps_pvals_plotting$feature], 0.6))

    abline(a = 0, b = 1, lty = 2, lwd = 4)
    legend("topleft", box.lty=0, legend=c('Maggie python regression', 'Maggie R lm regression', levels(as.factor(snps_pvals_plotting$feature))), col= c('black', 'black', alpha(col, 0.6)), pch = c(17, 15, rep(19, length(unique(snps_pvals_plotting$feature)))), cex = 1.5)
    dev.off()
}

# plot!
compare_phil_me()
