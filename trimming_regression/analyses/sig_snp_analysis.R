library(data.table)
setDTthreads(1)

source('/home/mrussel2/tcr-gwas/trimming_regression/scripts/manha_visualization.R')

trim_type = 'all_trim'
random_effects = 'True'
bootstrap_count = 0
condensing = 'by_gene'
d_infer = 'True'
maf_cutoff = 0.05
bootstrap_rerun_count = 100
pca_structure_correction = 'True'



bonferroni = 0.05/35481497
maf_data = fread(paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/maf_all_snps.tsv'), sep = "\t", fill=TRUE, header = TRUE)[,-c(1)]

if (trim_type == 'all_trim'){
    data = data.frame()
    for (trim in c('v_trim', 'd0_trim', 'd1_trim', 'j_trim')){
       data = rbind(data, compile_manhattan_plot_data(trim_type = trim, maf_data, random_effects, condensing, d_infer, bootstrap_count, maf_cutoff, bootstrap_rerun_count, pca_structure_correction))
    }
} else if (trim_type == 'all_insert'){
    data = data.frame()
    for (trim in c('vj_insert', 'vd_insert', 'dj_insert')){
       data = rbind(data, compile_manhattan_plot_data(trim_type = trim, maf_data, random_effects, condensing, d_infer, bootstrap_count, maf_cutoff, bootstrap_rerun_count, pca_structure_correction))
    }
}
    
genes = c('artemis', 'mhc', 'dntt', 'rag', 'tcrb', 'tcra')
chr = c(10, 6, 10, 11, 7, 14)
pos1 = c(14939358, 25912984, 98064085, 36510709, 141998851, 22090057)
pos2 = c(14996431, 33290793, 98098321, 36593156, 142510972, 23021075)

gene_annotations = data.frame(genes = genes, chr = chr, pos1 = pos1, pos2 = pos2)

for (index in seq(nrow(gene_annotations))){
   data[chr == gene_annotations[index,]$chr & hg19_pos  > gene_annotations[index,]$pos1 & hg19_pos  < gene_annotations[index,]$pos2, gene := gene_annotations[index,]$genes]
}

sig_data = data[pvalue < bonferroni]

print(sig_data[,.N, by = .(productivity, trim_type, gene)])

artemis_sigs = sig_data[gene == 'artemis']
    
print(head(artemis_sigs[order(pvalue)][1:25]))

tcrb_sigs = sig_data[gene == 'tcrb']

print(head(tcrb_sigs[order(pvalue)][1:25]))


