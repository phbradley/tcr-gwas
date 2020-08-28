library("data.table")
library("plyr")
library("readr")
library("stringr")
library("gdsfmt")
library("SNPRelate")
library("GWASTools")

snp_file <- function(chromosome, position1, position2){
    snps_gds = snpgdsOpen("../_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds")
    snp_chrom <- read.gdsn(index.gdsn(snps_gds, "snp.chromosome"))

    chr_index_start = match(chromosome,snp_chrom)
    chr_index_end = match(chromosome+1,snp_chrom)-1

    snp_pos <- read.gdsn(index.gdsn(snps_gds, "snp.position"), start = chr_index_start, count=(chr_index_end-chr_index_start)+1)
    snp_id <- read.gdsn(index.gdsn(snps_gds, "snp.id"), start = chr_index_start, count=(chr_index_end-chr_index_start)+1)

    snps_chromosome = data.table(snp = snp_id, chr = snp_chrom[chr_index_start:chr_index_end], hg19_pos = snp_pos)

    if (position1 != 'all'){
        snps_chromosome = snps_chromosome[hg19_pos < position2 & hg19_pos > position1]
    }
    closefn.gds(snps_gds)
    return(snps_chromosome)
}

snp_file_by_snp_start <- function(snp_start, count){
    snps_gds = snpgdsOpen("../_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds")
    snp_id <- read.gdsn(index.gdsn(snps_gds, "snp.id"))

    if ((snp_start+count) > length(snp_id)){
        count = length(snp_id) - snp_start
    }

    snp_chrom <- read.gdsn(index.gdsn(snps_gds, "snp.chromosome"), start = snp_start, count=count)
    snp_pos <- read.gdsn(index.gdsn(snps_gds, "snp.position"), start = snp_start, count=count)
    
    snps = data.table(snp = snp_id[snp_start:(snp_start+count-1)], chr = snp_chrom, hg19_pos = snp_pos)

    closefn.gds(snps_gds)
    return(snps)
}
