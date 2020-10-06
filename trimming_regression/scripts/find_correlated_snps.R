library("COMBAT")
library("tidyverse")
library("gdsfmt")
library("SNPRelate")
library("GWASTools")

setDTthreads(threads = 1, restore_after_fork=FALSE)


correlate_snps_ld <- function(chrom, snp_list, cutoff){
    snps_gds = snpgdsOpen("/home/mrussel2/tcr-gwas/_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds")
    snp_chrom <- read.gdsn(index.gdsn(snps_gds, "snp.chromosome"))
    sampleid <- read.gdsn(index.gdsn(snps_gds, "sample.id"))

    chr_index_start = match(chrom,snp_chrom)
    chr_index_end = match(chrom+1,snp_chrom)-1

    snp_pos <- read.gdsn(index.gdsn(snps_gds, "snp.position"), start = chr_index_start, count=(chr_index_end-chr_index_start)+1)
    snp_id <- read.gdsn(index.gdsn(snps_gds, "snp.id"), start = chr_index_start, count=(chr_index_end-chr_index_start)+1)
    genotypes <- read.gdsn(index.gdsn(snps_gds, "genotype"), start = c(1, chr_index_start), count=c(398, (chr_index_end-chr_index_start)+1))
    snp_id <- read.gdsn(index.gdsn(snps_gds, "snp.id"), start = chr_index_start, count=(chr_index_end-chr_index_start)+1)

    genotypes_matrix = as.matrix(genotypes)
    colnames(genotypes_matrix) = snp_id
    rownames(genotypes_matrix) = sampleid

    snp_list_chrom = snp_list[chromosome == paste(chrom)]
    snp_list_chrom = snp_list_chrom[order(hg19_pos)]

    genotypes_snp_list = matrix(nrow = 398, ncol = length(unique(snp_list_chrom$snp)))
    colnames(genotypes_snp_list) = unique(snp_list_chrom$snp)
    count = 1
    for (snp in unique(snp_list_chrom$snp)){
        index = match(snp,snp_id)
        genotypes_snp_list[,count] = genotypes_matrix[,index]
        count = count + 1
    }
    genotypes_snp_list[genotypes_snp_list==3]<-NA
    
    # get rid of columns containing ALL NA
    genotypes_snp_list_no_NA = genotypes_snp_list[, colMeans(is.na(genotypes_snp_list)) != 1]
    LD_matrix = ld.Rsquare(genotypes_snp_list_no_NA)

    # Now find LD blocks
    correlated_snps = data.table()
    start = 1
    cluster_id = 0
    for (end in 1:ncol(LD_matrix)){
        correlated = all(LD_matrix[start:end, start:end] > cutoff | LD_matrix[start:end, start:end] < -1*cutoff)
        if (start != end & end == ncol(LD_matrix)){
                cluster_id = cluster_id + 1
                correlated_snps = rbind(correlated_snps, data.table(snp = c(colnames(LD_matrix[start:end, start:end])), cluster = rep(cluster_id, length(colnames(LD_matrix[start:end, start:end])))))
                print(paste0('block', start, '_', (end)))
                next
            }
        if (correlated == TRUE){
            next
        } else {
            if (start != end-1){
                cluster_id = cluster_id + 1
                correlated_snps = rbind(correlated_snps, data.table(snp = c(colnames(LD_matrix[start:end-1, start:end-1])), cluster = rep(cluster_id, length(colnames(LD_matrix[start:end-1, start:end-1])))))
                print(paste0('block', start, '_', (end-1)))
            }
            start = end
        }
    }
    closefn.gds(snps_gds)
    return(correlated_snps)
}