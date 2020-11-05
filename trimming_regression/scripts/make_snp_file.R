library("tidyverse")
library("gdsfmt")
library("SNPRelate")
library("GWASTools")
library(data.table)
setDTthreads(threads = 1)

# This function makes a snp file from chromosme and position data using the gds file
snp_file <- function(chromosome, position1, position2){
    snps_gds = snpgdsOpen("/home/mrussel2/tcr-gwas/_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds")
    snp_chrom <- read.gdsn(index.gdsn(snps_gds, "snp.chromosome"))

    chr_index_start = match(chromosome,snp_chrom)
    chr_index_end = match(chromosome+1,snp_chrom)-1

    snp_pos <- read.gdsn(index.gdsn(snps_gds, "snp.position"), start = chr_index_start, count=(chr_index_end-chr_index_start)+1)
    snp_id <- read.gdsn(index.gdsn(snps_gds, "snp.id"), start = chr_index_start, count=(chr_index_end-chr_index_start)+1)

    snps_chromosome = data.frame(snp = snp_id, chr = snp_chrom[chr_index_start:chr_index_end], hg19_pos = snp_pos)

    if (position1 != 'all'){
        snps_chromosome = snps_chromosome[hg19_pos < position2 & hg19_pos > position1]
    }
    closefn.gds(snps_gds)
    return(snps_chromosome)
}

# This function makes a snp file from snp index start and snp count using the gds file
snp_file_by_snp_start <- function(snp_start, count){
    snps_gds = openfn.gds("/home/mrussel2/tcr-gwas/_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds")
    snp_id <- read.gdsn(index.gdsn(snps_gds, "snp.id"))

    if ((snp_start+count) > length(snp_id)){
        count = length(snp_id) - snp_start
    }

    snp_chrom <- read.gdsn(index.gdsn(snps_gds, "snp.chromosome"), start = snp_start, count=count)
    snp_pos <- read.gdsn(index.gdsn(snps_gds, "snp.position"), start = snp_start, count=count)
    
    snps = data.frame(snp = snp_id[snp_start:(snp_start+count-1)], chr = snp_chrom, hg19_pos = snp_pos)

    closefn.gds(snps_gds)
    return(snps)
}

# This function makes a genotype file from snp index start and snp count using the gds file
compile_all_genotypes <- function(snp_start, count) {
    gfile = openfn.gds("/home/mrussel2/tcr-gwas/_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds")
    bigsize <- 35481497
    #start <- (i-1)*sz + 1
    numrows <- min( count, bigsize-snp_start+1 )
    
    genotype_matrix <- read.gdsn(index.gdsn(gfile, "genotype"), start=c(1,snp_start), count=c(398, numrows))
    sample_ids <- read.gdsn(index.gdsn(gfile, "sample.id"), start=1, count=398)
    snp_ids <- read.gdsn(index.gdsn(gfile, "snp.id"), start=snp_start, count=numrows)
    closefn.gds(gfile)

    subject_id_mapping = as.data.frame(read.table('/home/mrussel2/tcr-gwas/_ignore/snp_data/gwas_id_mapping.tsv', sep = "\t", fill=TRUE, header = TRUE, check.names = FALSE))
    sample_ids_converted = plyr::mapvalues(sample_ids, subject_id_mapping$scanID, subject_id_mapping$localID)
    rownames(genotype_matrix) = sample_ids_converted
    colnames(genotype_matrix) = snp_ids
    return(genotype_matrix)
}

# This function compiles all subjects
compile_subjects <- function() {
    gfile = openfn.gds("/home/mrussel2/tcr-gwas/_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds")
    subjects <- read.gdsn(index.gdsn(gfile, "sample.id"))
    closefn.gds(gfile)
    return(subjects)
}

# This function finds the snp start index given chromosome and positions
find_snp_start_by_position <- function(chromosome, position1, position2){
    library(data.table)
    snp_meta_data = fread('/home/mrussel2/tcr-gwas/_ignore/snps_meta_data.tsv', sep = "\t", fill=TRUE, header = TRUE, check.names = FALSE)
    colnames(snp_meta_data) =c('snpindex', 'snpid', 'snppos', 'snpchrome', 'snpallele')
    filtered = snp_meta_data[snpchrome == chromosome & snppos > position1 & snppos < position2]
    filtered_ordered = filtered[order(filtered$snpindex),]
    return(c(filtered_ordered$snpindex[1], filtered_ordered$snpindex[nrow(filtered_ordered)]-filtered_ordered$snpindex[1]))
}

# This function finds a regression file (from cluster) and subsets based on signficance cutoff (for use in bootstraps)
open_regressed_file_and_subset_by_pval <- function(significance_cutoff, trim_type, random_effects, condensing, d_infer, maf_cutoff){
    if (random_effects == 'True'){
         file_name = paste0('_', trim_type, '_snps_regression_with_weighting_condensing_', condensing, '_with_random_effects_0_bootstraps')
    } else {
         file_name = paste0('_',trim_type, '_snps_regression_with_weighting_condensing_', condensing, '_NO_random_effects_0_bootstraps')
    }

    if (d_infer == 'False'){
        file_name = paste0(file_name, '_NO_d_infer.tsv')
    } else {
        file_name = paste0(file_name, '.tsv')
    }

    productive_data = fread(paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/results/productive', file_name), sep = "\t", fill=TRUE, header = TRUE)
    not_productive_data = fread(paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/results/NOT_productive', file_name), sep = "\t", fill=TRUE, header = TRUE)
    maf_data = fread(paste0('/fh/fast/matsen_e/shared/tcr-gwas/trimming_regression_output/maf_all_snps.tsv'), sep = "\t", fill=TRUE, header = TRUE)[,-c(1)]

    data = rbind(productive_data, not_productive_data)[,-c(1,2)]
    data = merge(data, maf_data, by = 'snp')
    data = data[maf > maf_cutoff]
    data = data[pvalue < significance_cutoff]
    return(data[,c(1:3)])
}

# This function takes a snps dataframe and makes a coninciding snp file (with chromosome and position columns)
make_snp_file_subset_by_count_and_index <- function(snp_dataframe, count, index){
    total = nrow(snp_dataframe)
    converted_index = index - 1
    start = converted_index * count + 1
    if (start > nrow(snp_dataframe)){
        return(data.frame())
    } else {
        end = min(start-1 + count, total)
        return(unique(snp_dataframe[start:end,]))
    }  
}

# This function makes a genotype file given a random snp file
make_genotype_file_given_random_snps <- function(snp_file){
    snps_gds = snpgdsOpen("/home/mrussel2/tcr-gwas/_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds")
    snp_id_list = as.numeric(gsub('snp', '', snp_file$snp))
    subject_id_mapping = as.data.frame(read.table('/home/mrussel2/tcr-gwas/_ignore/snp_data/gwas_id_mapping.tsv', sep = "\t", fill=TRUE, header = TRUE, check.names = FALSE))
    genotype_list = data.frame(localID = subject_id_mapping$localID)
    for (snp in snp_id_list){
        genotype = snpgdsGetGeno(snps_gds, snp.id=snp, with.id = TRUE)
        genotype_matrix = data.frame(genotype$genotype, scanID = as.numeric(as.character(genotype$sample.id)))
        subject_id_HIP_genotypes = merge( genotype_matrix, subject_id_mapping, by = "scanID")
        colnames(subject_id_HIP_genotypes) = c('scanID', paste(snp), 'localID')
        genotype_list = merge(genotype_list, subject_id_HIP_genotypes[,c('localID', paste(snp))], by = 'localID')
    }
    rownames(genotype_list) = genotype_list$localID
    closefn.gds(snps_gds)
    return(genotype_list)
}

