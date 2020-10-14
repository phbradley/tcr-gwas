##library(plotly)
##library(manhattanly)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(GWASTools)
library(tidyverse)
#omp_set_num_threads(1)
#blas_set_num_threads(1)
setDTthreads(threads = 1)

source('/home/mrussel2/tcr-gwas/trimming_regression/scripts/lmer_trimming_regression_functions.R')
source('/home/mrussel2/tcr-gwas/trimming_regression/scripts/genomic_control_calculations.R')

# This script makes a manhattan plot for a specific region
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

# This script compiles all regression data from cluster across the entire genome

compile_all_data_from_cluster <- function(trim_type, random_effects, bootstrap_count){
    bonferroni = 0.05/35481497

    if (random_effects == 'True'){
         file_pattern = paste0('*_snps_regression_with_weighting_condensing_by_gene_with_random_effects_', bootstrap_count, '_bootstraps.tsv')
    } else {
         file_pattern = paste0('no_random/*_snps_regression_with_weighting_condensing_by_gene_NO_random_effects_', bootstrap_count, '_bootstraps.tsv')
    }

    data_files = list.files(path=paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/cluster_job_results/', trim_type, '/', bootstrap_count, '_bootstraps'), pattern=file_pattern, full.names=TRUE)   
    #chrom_sizes = as.data.frame(read.table('/home/mrussel2/tcr-gwas/_ignore/chrom_sizes_hg19.tsv', sep = "\t", header = FALSE))
    #colnames(chrom_sizes) = c('chr', 'length')

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
    together = together[order(together$chr, together$hg19_pos),]
    assign(paste0('together_productive_', trim_type), together %>% filter(productivity == 'productive'))
    assign(paste0('together_NOT_productive_', trim_type), together %>% filter(productivity == 'NOT_productive'))
    print(paste0('processed data for ', trim_type))

    together_productive = get(paste0('together_productive_', trim_type))
    together_NOT_productive = get(paste0('together_NOT_productive_', trim_type))

    if (random_effects == 'True'){
         file_name = paste0(trim_type, '_snps_regression_with_weighting_condensing_by_gene_with_random_effects_', bootstrap_count, '_bootstraps.tsv')
    } else {
         file_name = paste0(trim_type, '_snps_regression_with_weighting_condensing_by_gene_NO_random_effects_', bootstrap_count, '_bootstraps.tsv')
    }

    write.table(together_productive, file=paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/results/productive_', file_name), quote=FALSE, sep='\t', col.names = NA)
    write.table(together_NOT_productive, file=paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/results/NOT_productive_', file_name), quote=FALSE, sep='\t', col.names = NA)
}

# This script compiles all minor allele fraction data (to be used in plotting cutoff)
compile_all_maf_data <- function(){
    data_files = list.files(path=paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/maf_results/'), pattern='*', full.names=TRUE)   

    print(paste0('compile file list'))
    assign(paste0('together'), NULL)
    count = 0
    for (file in data_files){
        # read file...
        if (file.size(file) == 1 | file.size(file) == 0){
            next
        }
        temp_file = fread(file, sep = "\t", fill=TRUE, header = TRUE)
        #together = rbind(together, temp_file)
        if (ncol(temp_file) > 2){
            assign(paste0('together'), rbindlist(list(get(paste0('together')), temp_file)))
        }
        count = count + 1
        print(paste0(count, ' of ', length(data_files), ' completed'))
    }
    
    file_name = paste0('maf_all_snps.tsv')

    write.table(together, file=paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/', file_name), quote=FALSE, sep='\t', col.names = NA)
}

# This script makes a manhattan plot for the entire genome!

manhattan_plot_cluster <- function(trim_type, random_effects, bootstrap_count, plotting_cutoff, gene_annotations, maf_cutoff, bootstrap_rerun_count, genomic_control){
    bonferroni = 0.05/35481497
    if (random_effects == 'True'){
         file_name = paste0('_', trim_type, '_snps_regression_with_weighting_condensing_by_gene_with_random_effects_', bootstrap_count, '_bootstraps.tsv')
    } else {
         file_name = paste0('_',trim_type, '_snps_regression_with_weighting_condensing_by_gene_NO_random_effects_', bootstrap_count, '_bootstraps.tsv')
    }
    productive_data = fread(paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/results/productive', file_name), sep = "\t", fill=TRUE, header = TRUE)
    not_productive_data = fread(paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/results/NOT_productive', file_name), sep = "\t", fill=TRUE, header = TRUE)
    maf_data = fread(paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/maf_all_snps.tsv'), sep = "\t", fill=TRUE, header = TRUE)[,-c(1)]

    if (bootstrap_rerun_count != 'False'){
        if (random_effects == 'True'){
            file_name = paste0('_', trim_type, '_snps_regression_with_weighting_condensing_by_gene_with_random_effects_', bootstrap_rerun_count, '_bootstraps.tsv')
        } else {
            file_name = paste0('_',trim_type, '_snps_regression_with_weighting_condensing_by_gene_NO_random_effects_', bootstrap_rerun_count, '_bootstraps.tsv')
        }
        productive_boot_data = fread(paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/results/productive', file_name), sep = "\t", fill=TRUE, header = TRUE)
        productive_boot_data = productive_boot_data[!duplicated(productive_boot_data$snp)]
        not_productive_boot_data = fread(paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/results/NOT_productive', file_name), sep = "\t", fill=TRUE, header = TRUE)
        not_productive_boot_data = not_productive_boot_data[!duplicated(not_productive_boot_data$snp)]
        data_boot = rbind(productive_boot_data, not_productive_boot_data)[,-c(1,2)]
        data_boot$bootstraps = '100_bootstraps'

        # remove snps from non-bootstrapped dataframe if they have been boostrapped
        productive_data = productive_data[!productive_data$snp %in% productive_boot_data$snp,]
        not_productive_data = not_productive_data[!not_productive_data$snp %in% not_productive_boot_data$snp,]
        data = rbind(productive_data, not_productive_data)[,-c(1,2)]
        data$bootstraps = '0_bootstraps'

        data = rbind(data, data_boot)
        data = merge(data, maf_data, by = 'snp')
    } else {
        data = rbind(productive_data, not_productive_data)[,-c(1,2)]
        data = merge(data, maf_data, by = 'snp')
    }

    if (maf_cutoff != 'False'){
        data = data[maf > maf_cutoff]
        title = paste0(trim_type, ' plot with ', maf_cutoff, ' MAF cutoff')
    } else {
        title = paste0(trim_type, ' plot')
    }

    if (genomic_control == 'True'){
        data = genomic_control_calculation(data)
        data$pvalue = data$pvalue_genomic_control_correction
        title = paste0(title, ' with genomic control pvalue correction')
    }

    if (plotting_cutoff != 'False'){
        data = data[-log10(pvalue) > plotting_cutoff]
        if (nrow(data)==0){
            data = rbind(productive_data, not_productive_data)[,-c(1,2)]
        }
    }

    genes = c('artemis', 'mhc', 'dntt', 'rag', 'tcrb', 'tcra')
    chr = c(10, 6, 10, 11, 7, 14)
    pos1 = c(14939358, 25912984, 98064085, 36510709, 141998851, 22090057)
    pos2 = c(14996431, 33290793, 98098321, 36593156, 142510972, 23021075)

    gene_annotations = data.frame(genes = genes, chr = chr, pos1 = pos1, pos2 = pos2)

    if (gene_annotations == 'False'){
        if (bootstrap_rerun_count != 'False'){
            ggplot(data) + geom_point(aes(x = hg19_pos, y = -log10(pvalue), color=productivity, shape = bootstraps), alpha = 0.5, size = 3) + facet_grid(.~chr, switch="both", space='free_x', scales = "free_x") + theme_classic() + theme(panel.spacing.x=unit(0, "lines"), text = element_text(size = 30), axis.text.x = element_text(size = 8, angle = 90))+ labs(y="-log10(p-value)", x="Chromosome Position") + geom_hline(yintercept=-log10(bonferroni), linetype="dashed", color = "green4", size=1.5) + ggtitle(title) + scale_x_continuous(breaks=seq(0, 2.5e8, 0.75e8))
        } else {
            ggplot(data) + geom_point(aes(x = hg19_pos, y = -log10(pvalue), color = productivity), alpha = 0.5, size = 3) + facet_grid(.~chr, switch="both", space='free_x', scales = "free_x") + theme_classic() + theme(panel.spacing.x=unit(0, "lines"), text = element_text(size = 30), axis.text.x = element_text(size = 8, angle = 90))+ labs(y="-log10(p-value)", x="Chromosome Position") + geom_hline(yintercept=-log10(bonferroni), linetype="dashed", color = "green4", size=1.5) + ggtitle(title) + scale_x_continuous(breaks=seq(0, 2.5e8, 0.75e8))
        }
    } else {
        if (bootstrap_rerun_count != 'False'){
            ggplot(data) + geom_point(aes(x = hg19_pos, y = -log10(pvalue), color=productivity, shape = bootstraps), alpha = 0.5, size = 3) + geom_rect(data = gene_annotations, aes(xmin = pos1, xmax = pos2, ymin = -Inf, ymax = Inf, fill = genes), alpha = 0.5) + facet_grid(.~chr, switch="both", space='free_x', scales = "free_x") + theme_classic() + theme(panel.spacing.x=unit(0, "lines"), text = element_text(size = 30), axis.text.x = element_text(size = 8, angle = 90))+ labs(y="-log10(p-value)", x="Chromosome Position") + geom_hline(yintercept=-log10(bonferroni), linetype="dashed", color = "green4", size=1.5) + ggtitle(title) + scale_x_continuous(breaks=seq(0, 2.5e8, 0.75e8))
        } else {
            ggplot(data) + geom_point(aes(x = hg19_pos, y = -log10(pvalue), color=productivity), alpha = 0.5, size = 3) + geom_rect(data = gene_annotations, aes(xmin = pos1, xmax = pos2, ymin = -Inf, ymax = Inf, fill = genes), alpha = 0.5) + facet_grid(.~chr, switch="both", space='free_x', scales = "free_x") + theme_classic() + theme(panel.spacing.x=unit(0, "lines"), text = element_text(size = 30), axis.text.x = element_text(size = 8, angle = 90))+ labs(y="-log10(p-value)", x="Chromosome Position") + geom_hline(yintercept=-log10(bonferroni), linetype="dashed", color = "green4", size=1.5) + ggtitle(title) + scale_x_continuous(breaks=seq(0, 2.5e8, 0.75e8))  
        }  
    }
}

# This script makes a manhattan plot for a specific regression file

manhattan_plot_specific_file <- function(filename, trim_type, add_gene_annotations, maf_cutoff, bootstrap_count){
    bonferroni = 0.05/35481497
    data = fread(filename, sep = "\t", fill=TRUE, header = TRUE)
    maf_data = fread(paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/maf_all_snps.tsv'), sep = "\t", fill=TRUE, header = TRUE)[,-c(1)]

    data = merge(data, maf_data, by = 'snp')

    if (maf_cutoff != 'False'){
        data = data[maf > maf_cutoff]
        title = paste0(trim_type, ' plot with ', bootstrap_count, ' bootstraps and ', maf_cutoff, ' MAF cutoff')
    } else {
        title = paste0(trim_type, ' plot with ', bootstrap_count, ' bootstraps')
    }

    genes = c('artemis', 'mhc', 'dntt', 'rag', 'tcrb', 'tcra')
    chr = c(10, 6, 10, 11, 7, 14)
    pos1 = c(14939358, 25912984, 98064085, 36510709, 141998851, 22090057)
    pos2 = c(14996431, 33290793, 98098321, 36593156, 142510972, 23021075)

    gene_annotations = data.frame(genes = genes, chr = chr, pos1 = pos1, pos2 = pos2)

    if (add_gene_annotations == 'False'){
        ggplot(data) + geom_point(aes(x = hg19_pos, y = -log10(pvalue), color=productivity), alpha = 0.5, size = 2.5) + facet_grid(.~chr, switch="both", space='free_x', scales = "free_x") + theme_classic() + theme(panel.spacing.x=unit(0, "lines"), text = element_text(size = 30), axis.text.x = element_text(size = 8, angle = 90))+ labs(y="-log10(p-value)", x="Chromosome Position") + geom_hline(yintercept=-log10(bonferroni), linetype="dashed", color = "green4", size=1.5) + ggtitle(title) + scale_x_continuous(breaks=seq(0, 2.5e8, 0.75e8))
    } else if (add_gene_annotations == 'True'){
        ggplot(data) + geom_point(aes(x = hg19_pos, y = -log10(pvalue), color=productivity), alpha = 0.5, size = 2.5) + geom_rect(data = gene_annotations, aes(xmin = pos1, xmax = pos2, ymin = -Inf, ymax = Inf, fill = genes), alpha = 0.5) + facet_grid(.~chr, switch="both", space='free_x', scales = "free_x") + theme_classic() + theme(panel.spacing.x=unit(0, "lines"), text = element_text(size = 30), axis.text.x = element_text(size = 8, angle = 90))+ labs(y="-log10(p-value)", x="Chromosome Position") + geom_hline(yintercept=-log10(bonferroni), linetype="dashed", color = "green4", size=1.5) + ggtitle(title) + scale_x_continuous(breaks=seq(0, 2.5e8, 0.75e8))
    }
}
