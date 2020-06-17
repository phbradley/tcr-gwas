library(data.table)
library(plyr)
library(readr)
library(stringr)
library(gdsfmt)
library(SNPRelate)
library(GWASTools)

source("trimming_basic_regression_functions.R")
source("trimming_bootstrap_functions.R")

# Read in trimming data
v_trimming = as.data.table(read.table("../_ignore/condensed_trim_data_all_patients.tsv", sep = "\t", fill=TRUE, header = TRUE)[-1])
colnames(v_trimming) = c("localID", "v_gene", "productive", "v_trim", "v_gene_count", "weighted_v_gene_count")

#d_trimming = as.data.table(read.table("../_ignore/condensed_d_trim_data_all_patients.tsv", sep = "\t", fill=TRUE, header = TRUE)[-1])
#colnames(d_trimming) = c("localID", "d_gene", "productive", "d0_trim", "d1_trim", "d_gene_count", "weighted_d_gene_count")

#j_trimming = as.data.table(read.table("../_ignore/condensed_j_trim_data_all_patients.tsv", sep = "\t", fill=TRUE, header = TRUE)[-1])
#colnames(j_trimming) = c("localID", "j_gene", "productive", "j_trim", "j_gene_count", "weighted_j_gene_count")



# ALL SNP DATA
# Get snp meta data
snps_gds = openfn.gds("../_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds")
n <- index.gdsn(snps_gds, "sample.id")
sampleid <- read.gdsn(n) 
#sz <- 10
sz <- 100000
#nloop <- 10
nloop <- 355
bigsize <- 35481497
for ( i in 1:nloop ) {
    bootstrap_regression_results_v_gene_productive = data.frame()
    bootstrap_regression_results_v_gene_NOT_productive = data.frame()
    start <- (i-1)*sz + 1
    numrows <- min( sz, bigsize-start+1 ) 
    genotypes <- read.gdsn(index.gdsn(snps_gds, "genotype"), start=c(1,start), count = c(398, numrows))
    genotypes_df = data.frame(scanID = c(sampleid), genotypes)
    snpid <- read.gdsn(index.gdsn(snps_gds, "snp.id"), start = start, count=numrows)
    colnames(genotypes_df) = c("scanID", paste0("snp",snpid))
    genotypes_dt = as.data.table(genotypes_df)
    genotypes_dt$scanID = as.numeric(as.character(genotypes_dt$scanID))

    # Convert subject names and compile condensed data: 

    subject_id_mapping = as.data.table(read.table('../_ignore/snp_data/gwas_id_mapping.tsv', sep = "\t", fill=TRUE, header = TRUE, check.names = FALSE))

    snps = merge(genotypes_dt, subject_id_mapping, by = "scanID")
    print(paste0("snp data ", i, " of ", nloop, " compiled"))

    # Replace '3' genotypes with NA (missing genotype)
    snps[snps==3]<-NA
    snps = as.data.table(snps)

    # Get ride of snp columns which have NA for all entries (missing for all individuals)
    snps_no_NA = snps[,which(unlist(lapply(snps, function(x)!all(is.na(x))))),with=F]


    weighted_regression_patient_vgene_results_productive = trimming_snp_regression_weighted(snps_no_NA, v_trimming, productive = "True", gene_type = 'v_gene')
    print("finished weighted_regression_patient_vgene_results_productive")

    weighted_regression_patient_vgene_results_NOT_productive = trimming_snp_regression_weighted(snps_no_NA, v_trimming, productive = "False", gene_type = 'v_gene')
    print("finished weighted_regression_patient_vgene_results_NOT_productive")

    bootstrap_weighted_v_gene_productive = regression_weighted_bootstrap_se(snps_no_NA, v_trimming, productive = "True", repetitions = 500, gene_type = 'v_gene')
    print("finished bootstrap_weighted_productive")

    bootstrap_weighted_v_gene_NOT_productive = regression_weighted_bootstrap_se(snps_no_NA, v_trimming, productive = "False", repetitions = 500, gene_type = 'v_gene')
    print("finished bootstrap_weighted_NOT_productive")

    if (nrow(bootstrap_weighted_v_gene_productive) != 0){
        bootstrap_regression_results_v_gene_productive = rbind(bootstrap_regression_results_v_gene_productive, bootstrap_regression_combine(bootstrap_weighted_v_gene_productive, weighted_regression_patient_vgene_results_productive))
    }
    print("finished bootstrap_regression_weighted_productive")

    if (nrow(bootstrap_weighted_v_gene_NOT_productive) != 0){
        bootstrap_regression_results_v_gene_NOT_productive = rbind(bootstrap_regression_results_v_gene_NOT_productive, bootstrap_regression_combine(bootstrap_weighted_v_gene_NOT_productive, weighted_regression_patient_vgene_results_NOT_productive))
    }
    print("finished bootstrap_regression_weighted_NOT_productive")

    write.table(bootstrap_regression_results_v_gene_NOT_productive, file=paste0('regression_bootstrap_results/productive/v_gene/v_gene_NOT_productive_', i, 'tsv'), quote=FALSE, sep='\t', col.names = NA)
    write.table(bootstrap_regression_results_v_gene_productive, file=paste0('regression_bootstrap_results/NOT_productive/v_gene/v_gene_productive_', i, '.tsv'), quote=FALSE, sep='\t', col.names = NA)
    
    print(paste0("finished bootstrap for snp data", i, " of ", nloop))
}

closefn.gds(snps_gds)



    

