library(ggplot2)
library(data.table)
setDTthreads(1)
library(GWASTools)

args = commandArgs(trailingOnly=TRUE)

type = args[1]
stopifnot(type %in% c('insertions', 'trimming'))

if (type =='insertions'){
    types = c('vd_insert', 'dj_insert')
} else {
    types = c('v_trim', 'j_trim', 'd1_trim', 'd1_trim')
}
pca_structure_correction = args[2]
pca_type = args[3]
project_path = args[4]
output_path = args[5]

source(paste0(project_path, '/tcr-gwas/trimming_regression/scripts/compile_regression_functions.R'))
source(paste0(project_path, '/tcr-gwas/trimming_regression/scripts/compile_data_functions.R'))


calculate_lambda <- function(genome_wide_pvalue_dt){
    chisq = (genome_wide_pvalue_dt$slope/genome_wide_pvalue_dt$standard_error)**2

    lambda = median(chisq, na.rm=TRUE)/qchisq(0.5, 1)
    return(lambda)
}

make_qq_plot <- function(productivity, trim_type, weighting, condensing, random_effects, d_infer, repetitions, pca_structure_correction, pca_type){
    genome_wide_pvalue_file = make_compiled_regression_file_name(productivity, trim_type, weighting, condensing, random_effects, d_infer, repetitions, pca_structure_correction, pca_type)
    if(file.exists(genome_wide_pvalue_file)){
        genome_wide_pvalue_dt = fread(file = paste(genome_wide_pvalue_file))
        lambda = calculate_lambda(genome_wide_pvalue_dt) 
        file_name = strsplit(strsplit(genome_wide_pvalue_file, '/')[[1]][9], '.tsv')[[1]][1]
        png(paste0('figures/', file_name, '.png'), units="px", width=1600, height=1600, res=300, pointsize = 8)
        title = paste0(trim_type, ' for ', productivity, ' TCRs with ', pca_type, ' pca correction')
        subtitle = paste0("lambda = ", format(lambda, digits=5))
        qqPlot(genome_wide_pvalue_dt$pvalue, main = title, sub = subtitle, ylim = c(0,max(-log10(genome_wide_pvalue_dt$pvalue))), xlim = c(0, max(-log10(genome_wide_pvalue_dt$pvalue))))
        dev.off()
    } else {
        print(paste0(genome_wide_pvalue_file, ' file does not exist, so no plot'))
    }
}

print(types)
for (trim in types){
    set_regression_parameters(trim)
    for (prod in c('productive', 'NOT_productive')){
        make_qq_plot(productivity = prod, trim_type = trim, weighting, condensing, random_effects, d_infer, repetitions, pca_structure_correction, pca_type)
        print(paste0('finished plotting for ', prod, ' ', trim, ' ', pca_type))
    }
}
