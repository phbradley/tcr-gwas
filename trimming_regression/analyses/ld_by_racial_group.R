library("tidyverse")
library("gdsfmt")
library("SNPRelate")
library("GWASTools")
library(data.table)
setDTthreads(threads = 1)

 This function makes a snp file from chromosme and position data using the gds file
snp_file <- function(chromosome, position1, position2){
    snps_gds = snpgdsOpen(paste0(PROJECT_PATH, "/tcr-gwas/_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds"))
    snp_chrom <- read.gdsn(index.gdsn(snps_gds, "snp.chromosome"))

    chr_index_start = match(chromosome,snp_chrom)
    chr_index_end = match(chromosome+1,snp_chrom)-1

    snp_pos <- read.gdsn(index.gdsn(snps_gds, "snp.position"), 
                         start = chr_index_start, 
                         count=(chr_index_end-chr_index_start)+1)
    snp_id <- read.gdsn(index.gdsn(snps_gds, "snp.id"), 
                        start = chr_index_start, 
                        count=(chr_index_end-chr_index_start)+1)

    snps_chromosome = data.frame(snp = snp_id, 
                                 chr = snp_chrom[chr_index_start:chr_index_end], 
                                 hg19_pos = snp_pos)

    if (position1 != 'all'){
        snps_chromosome = snps_chromosome[hg19_pos < position2 & hg19_pos > position1]
    }
    closefn.gds(snps_gds)
    return(snps_chromosome)
}


