library(data.table)
library(plyr)
library(readr)
library(stringr)
library(gdsfmt)
library(SNPRelate)
library(GWASTools)
library(SeqArray) # May not need this

setwd("../_ignore/")

# ALL SNP DATA
# Get snp meta data
snps_gds = openfn.gds("snp_data/HSCT_comb_geno_combined_v03_tcr.gds")

n <- index.gdsn(snps_gds, "sample.id")
sampleid <- read.gdsn(n) 


# Take a subset of the genotype data into two individuals and 2500 snps(This is where I will get all of the data)
genotypes <- read.gdsn(index.gdsn(snps_gds, "genotype"), start=c(1,1))

row.names(genotypes) = c(sampleid)
genotypes_df = data.frame(scanID = row.names(genotypes), genotypes)
colnames(genotypes_df) = c("scanID", paste0("snp", 1:(ncol(genotypes_df)-1)))
genotypes_dt = as.data.table(genotypes_df)
genotypes_dt$scanID = as.numeric(as.character(genotypes_dt$scanID))

closefn.gds(snps_gds)


# Convert subject names and compile condensed data: 

subject_id_mapping = as.data.table(read.table('snp_data/gwas_id_mapping.tsv', sep = "\t", fill=TRUE, header = TRUE, check.names = FALSE))

genotypes_dt_bothID = merge(genotypes_dt, subject_id_mapping, by = "scanID")


    

