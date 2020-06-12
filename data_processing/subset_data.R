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

n <- index.gdsn(snps_gds, "snp.id")
snpid <- read.gdsn(n)

n <- index.gdsn(snps_gds, "snp.position")
snppos <- read.gdsn(n) #, start=c(1,start), count=c(398, numrows))

n <- index.gdsn(snps_gds, "snp.chromosome")
snpchrome <- read.gdsn(n)

n <- index.gdsn(snps_gds, "snp.allele")
snpallele <- read.gdsn(n)

snp_metadata <- data.frame(snpid, snppos, snpchrome, snpallele)
snp_metadata = as.data.table(snp_metadata)

# write table for all filtered snp data
write.table(snp_metadata, file='snps_meta_data.tsv', quote=FALSE, sep='\t', col.names = NA)

n <- index.gdsn(snps_gds, "sample.id")
sampleid <- read.gdsn(n) 

write.table(sampleid, file='snps_subject_data.tsv', quote=FALSE, sep='\t', col.names = NA)

# Take a subset of the genotype data into two individuals and 2500 snps(This is where I will get all of the data)
genotypes <- read.gdsn(index.gdsn(snps_gds, "genotype"), start=c(1,1), count=c(398,2500))

row.names(genotypes) = c(sampleid)
genotypes_df = data.frame(scanID = row.names(genotypes), genotypes)
colnames(genotypes_df) = c("scanID", paste0("snp", 1:(ncol(genotypes_df)-1)))
genotypes_dt = as.data.table(genotypes_df)
genotypes_dt$scanID = as.numeric(as.character(genotypes_dt$scanID))

closefn.gds(snps_gds)


# Convert subject names and compile condensed data: 

subject_id_mapping = as.data.table(read.table('snp_data/gwas_id_mapping.tsv', sep = "\t", fill=TRUE, header = TRUE, check.names = FALSE))

genotypes_dt_bothID = merge(genotypes_dt, subject_id_mapping, by = "scanID")
write.table(genotypes_dt_bothID, file='subset_snps_data_2500_snps.tsv', quote=FALSE, sep='\t')



# ALL VGENE TRIMMING DATA--condensed by patient
condensed_trimming_data = as.data.table(read.table('condensed_trim_data_all_patients.tsv', sep = "\t", fill=TRUE, header = TRUE, check.names = FALSE)[-1])
colnames(condensed_trimming_data) = c("localID", "v_gene", "productive", "avg_v_trim")
# subset data and write table for two patients
condensed_trimming_data_subset = condensed_trimming_data[localID %in% genotypes_dt$localID]
write.table(condensed_trimming_data_subset, file='subset_condensed_trimming_data_2_patients.tsv', quote=FALSE, sep='\t', col.names = NA)

    

