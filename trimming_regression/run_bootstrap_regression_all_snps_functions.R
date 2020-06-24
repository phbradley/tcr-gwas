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
v_trimming = as.data.table(read.table("../_ignore/condensed_v_trim_data_all_patients.tsv", sep = "\t", fill=TRUE, header = TRUE)[-1])
colnames(v_trimming) = c("localID", "v_gene", "productive", "v_trim", "v_gene_count", "weighted_v_gene_count")

d_trimming = as.data.table(read.table("../_ignore/condensed_d_trim_data_all_patients.tsv", sep = "\t", fill=TRUE, header = TRUE)[-1])
colnames(d_trimming) = c("localID", "d_gene", "productive", "d0_trim", "d1_trim", "d_gene_count", "weighted_d_gene_count")

# For now, exclude entries with missing d_gene
d_trimming2 = d_trimming[d_gene != "-"]

j_trimming = as.data.table(read.table("../_ignore/condensed_j_trim_data_all_patients.tsv", sep = "\t", fill=TRUE, header = TRUE)[-1])
colnames(j_trimming) = c("localID", "j_gene", "productive", "j_trim", "j_gene_count", "weighted_j_gene_count")

vj_insert = as.data.table(read.table("../_ignore/condensed_vj_insert_data_all_patients.tsv", sep = "\t", fill=TRUE, header = TRUE)[-1])
colnames(vj_insert) = c("localID", "v_gene", "j_gene", "productive", "vj_insert", "vj_gene_count", "weighted_vj_gene_count")

dj_insert = as.data.table(read.table("../_ignore/condensed_dj_insert_data_all_patients.tsv", sep = "\t", fill=TRUE, header = TRUE)[-1])
colnames(dj_insert) = c("localID", "d_gene", "j_gene", "productive", "dj_insert", "dj_gene_count", "weighted_dj_gene_count")

# For now, exclude entries with missing d_gene
dj_insert2 = dj_insert[d_gene != "-"]

vd_insert = as.data.table(read.table("../_ignore/condensed_vd_insert_data_all_patients.tsv", sep = "\t", fill=TRUE, header = TRUE)[-1])
colnames(vd_insert) = c("localID", "v_gene", "d_gene", "productive", "vd_insert", "vd_gene_count", "weighted_vd_gene_count")

# For now, exclude entries with missing d_gene
vd_insert2 = vd_insert[d_gene != "-"]




run_snps_trimming <- function(snps_per_run, number_of_runs, trim_type, varying_int){
    # ALL SNP DATA
    # Get snp meta data
    snps_gds = openfn.gds("../_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds")
    n <- index.gdsn(snps_gds, "sample.id")
    sampleid <- read.gdsn(n) 
    sz <- snps_per_run
    #sz <- 100000
    nloop <- number_of_runs
    #nloop <- 355
    bigsize <- 35481497
    for ( i in 1:nloop ) {
        bootstrap_regression_results_v_gene_productive = data.table()
        bootstrap_regression_results_v_gene_NOT_productive = data.table()
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
        print(paste0("snp data ", i, " of ", nloop, " compiled for ", trim_type))

        # Replace '3' genotypes with NA (missing genotype)
        snps[snps==3]<-NA
        snps = as.data.table(snps)

        # Get rid of snp columns which have the same entry for all entries and get ride of snp columns which have NA for all entries (missing for all individuals)
        snps_no_NA = Filter(function(x) length(unique(x))!=1, snps)
        snps_no_N2 = Filter(function(x) length(unique(snps_no_NA[x != "NA"])) != 1, snps_no_NA)
        snps_no_NA2 = data.frame(localID = snps_no_NA$localID)
        for (column in names(snps_no_NA)[-c(1, ncol(snps_no_NA))]){
            if (length(unique(na.omit(snps_no_NA[[column]]))) > 1){
                snps_no_NA2 = cbind(snps_no_NA2, snps_no_NA[[column]])
                names(snps_no_NA2)[names(snps_no_NA2) == 'snps_no_NA[[column]]'] <- paste0(column)
            }
        }
        if (trim_type == "v_trim"){
            trimming_data = v_trimming
        } else if (trim_type == "d0_trim" | trim_type == "d1_trim"){
            trimming_data = d_trimming2
        } else if (trim_type == "j_trim"){
            trimming_data = j_trimming
        } else if (trim_type == "vj_insert"){
            trimming_data = vj_insert
        } else if (trim_type == "dj_insert"){
            trimming_data = dj_insert2
        } else if (trim_type == "vd_insert"){
            trimming_data = vd_insert2
        }
        if (varying_int == "FALSE"){
            weighted_regression_bootstrap_patient_vgene_results_productive = trimming_snp_regression_weighted(snps_no_NA2, trimming_data, productive = "True", trim_type = trim_type, repetitions = 300)
            print("finished weighted_regression_patient_vgene_results_productive")

            weighted_regression_bootstrap_patient_vgene_results_NOT_productive = trimming_snp_regression_weighted(snps_no_NA2, trimming_data, productive = "False", trim_type = trim_type,  repetitions = 300)
            print("finished weighted_regression_patient_vgene_results_NOT_productive")
        } else if (varying_int == "TRUE"){
            weighted_regression_bootstrap_patient_vgene_results_productive = trimming_snp_regression_weighted_varying_int_subject(snps_no_NA2, trimming_data, productive = "True", trim_type = trim_type, repetitions = 300)
            print("finished weighted_regression_patient_vgene_results_productive")

            weighted_regression_bootstrap_patient_vgene_results_NOT_productive = trimming_snp_regression_weighted_varying_int_subject(snps_no_NA2, trimming_data, productive = "False", trim_type = trim_type,  repetitions = 300)
            print("finished weighted_regression_patient_vgene_results_NOT_productive")
        }   
        
        if (nrow(weighted_regression_bootstrap_patient_vgene_results_productive) != 0){
            bootstrap_regression_results_v_gene_productive = rbind(bootstrap_regression_results_v_gene_productive, weighted_regression_bootstrap_patient_vgene_results_productive)
            write.table(bootstrap_regression_results_v_gene_productive, file=paste0('regression_bootstrap_results/productive/', trim_type, '/', trim_type, '_productive_', i, '.tsv'), quote=FALSE, sep='\t', col.names = NA)
        }
        print("finished bootstrap_regression_weighted_productive")
        
        if (nrow(weighted_regression_bootstrap_patient_vgene_results_NOT_productive) != 0){
            bootstrap_regression_results_v_gene_NOT_productive = rbind(bootstrap_regression_results_v_gene_NOT_productive, weighted_regression_bootstrap_patient_vgene_results_NOT_productive)
            write.table(bootstrap_regression_results_v_gene_NOT_productive, file=paste0('regression_bootstrap_results/NOT_productive/', trim_type, '/', trim_type, '_NOT_productive_', i, '.tsv'), quote=FALSE, sep='\t', col.names = NA)
        }
        print("finished bootstrap_regression_weighted_NOT_productive")
      

        print(paste0("finished bootstrap for snp data ", i, " of ", nloop, " for ", trim_type))
    }

    closefn.gds(snps_gds)
}


run_snps_trimming_start_end <- function(snp_start, snp_end, trim_type, varying_int){
    # ALL SNP DATA
    # Get snp meta data
    snps_gds = openfn.gds("../_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds")
    n <- index.gdsn(snps_gds, "sample.id")
    sampleid <- read.gdsn(n) 
    bigsize <- 35481497
    bootstrap_regression_results_v_gene_productive = data.table()
    bootstrap_regression_results_v_gene_NOT_productive = data.table()
    numrows <- min( snp_end-snp_start, bigsize-snp_start+1 ) 
    genotypes <- read.gdsn(index.gdsn(snps_gds, "genotype"), start=c(1,snp_start), count = c(398, numrows))
    genotypes_df = data.frame(scanID = c(sampleid), genotypes)
    snpid <- read.gdsn(index.gdsn(snps_gds, "snp.id"), start = snp_start, count=numrows)
    colnames(genotypes_df) = c("scanID", paste0("snp",snpid))
    genotypes_dt = as.data.table(genotypes_df)
    genotypes_dt$scanID = as.numeric(as.character(genotypes_dt$scanID))
    # Convert subject names and compile condensed data: 
    subject_id_mapping = as.data.table(read.table('../_ignore/snp_data/gwas_id_mapping.tsv', sep = "\t", fill=TRUE, header = TRUE, check.names = FALSE))

    snps = merge(genotypes_dt, subject_id_mapping, by = "scanID")
    print(paste0("snp data compiled for ", snp_start, " to ", snp_end, " for ", trim_type))
    # Replace '3' genotypes with NA (missing genotype)
    snps[snps==3]<-NA
    snps = as.data.table(snps)
    # Get rid of snp columns which have the same entry for all entries and get ride of snp columns which have NA for all entries (missing for all individuals)
    snps_no_NA = Filter(function(x) length(unique(x))!=1, snps)
    snps_no_N2 = Filter(function(x) length(unique(snps_no_NA[x != "NA"])) != 1, snps_no_NA)
    snps_no_NA2 = data.frame(localID = snps_no_NA$localID)
    for (column in names(snps_no_NA)[-c(1, ncol(snps_no_NA))]){
        if (length(unique(na.omit(snps_no_NA[[column]]))) > 1){
            snps_no_NA2 = cbind(snps_no_NA2, snps_no_NA[[column]])
            names(snps_no_NA2)[names(snps_no_NA2) == 'snps_no_NA[[column]]'] <- paste0(column)
        }
    }
    if (trim_type == "v_trim"){
        trimming_data = v_trimming
    } else if (trim_type == "d0_trim" | trim_type == "d1_trim"){
        trimming_data = d_trimming2
    } else if (trim_type == "j_trim"){
        trimming_data = j_trimming
    } else if (trim_type == "vj_insert"){
        trimming_data = vj_insert
    } else if (trim_type == "dj_insert"){
        trimming_data = dj_insert2
    } else if (trim_type == "vd_insert"){
        trimming_data = vd_insert2
    }
    if (varying_int == "FALSE"){
        weighted_regression_bootstrap_patient_vgene_results_productive = trimming_snp_regression_weighted(snps_no_NA2, trimming_data, productive = "True", trim_type = trim_type, repetitions = 300)
        print("finished weighted_regression_patient_vgene_results_productive")

        weighted_regression_bootstrap_patient_vgene_results_NOT_productive = trimming_snp_regression_weighted(snps_no_NA2, trimming_data, productive = "False", trim_type = trim_type,  repetitions = 300)
        print("finished weighted_regression_patient_vgene_results_NOT_productive")
    } else if (varying_int == "TRUE"){
        weighted_regression_bootstrap_patient_vgene_results_productive = trimming_snp_regression_weighted_varying_int_subject(snps_no_NA2, trimming_data, productive = "True", trim_type = trim_type, repetitions = 300)
        print("finished weighted_regression_patient_vgene_results_productive")

        weighted_regression_bootstrap_patient_vgene_results_NOT_productive = trimming_snp_regression_weighted_varying_int_subject(snps_no_NA2, trimming_data, productive = "False", trim_type = trim_type,  repetitions = 300)
        print("finished weighted_regression_patient_vgene_results_NOT_productive")
    }   
        
    
    if (nrow(weighted_regression_bootstrap_patient_vgene_results_productive) != 0){
        bootstrap_regression_results_v_gene_productive = rbind(bootstrap_regression_results_v_gene_productive, weighted_regression_bootstrap_patient_vgene_results_productive)
        write.table(bootstrap_regression_results_v_gene_productive, file=paste0('regression_bootstrap_results/productive/', trim_type, '/', trim_type, '_productive_',snp_start,'_', snp_end,'.tsv'), quote=FALSE, sep='\t', col.names = NA)
    }
    print("finished bootstrap_regression_weighted_productive")

    if (nrow(weighted_regression_bootstrap_patient_vgene_results_NOT_productive) != 0){
        bootstrap_regression_results_v_gene_NOT_productive = rbind(bootstrap_regression_results_v_gene_NOT_productive, weighted_regression_bootstrap_patient_vgene_results_NOT_productive)
        write.table(bootstrap_regression_results_v_gene_NOT_productive, file=paste0('regression_bootstrap_results/NOT_productive/', trim_type, '/', trim_type, '_NOT_productive_', snp_start,'_', snp_end, '.tsv'), quote=FALSE, sep='\t', col.names = NA)
    }
    print("finished bootstrap_regression_weighted_NOT_productive")

    print(paste0("finished bootstrap for snp data ", snp_start, " to ", snp_end, " for ", trim_type))

    closefn.gds(snps_gds)
}
    

run_snps_trimming_snp_list <- function(snp_id_list, trim_type, varying_int, repetitions){
    # ALL SNP DATA
    # Get snp meta data
    snps_gds = snpgdsOpen("../_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds")
    n <- index.gdsn(snps_gds, "sample.id")
    sampleid <- read.gdsn(n) 
    bigsize <- 35481497
    bootstrap_regression_results_v_gene_productive = data.table()
    bootstrap_regression_results_v_gene_NOT_productive = data.table()
    i = 0
    for (snp in snp_id_list){
        i = i + 1
        genotype = snpgdsGetGeno(snps_gds, snp.id=snp)
        genotypes_df = data.frame(scanID = c(sampleid), genotype)
        colnames(genotypes_df) = c("scanID", paste0("snp",snp))
        genotypes_dt = as.data.table(genotypes_df)
        genotypes_dt$scanID = as.numeric(as.character(genotypes_dt$scanID))
        # Convert subject names and compile condensed data: 
        subject_id_mapping = as.data.table(read.table('../_ignore/snp_data/gwas_id_mapping.tsv', sep = "\t", fill=TRUE, header = TRUE, check.names = FALSE))

        snps = merge(genotypes_dt, subject_id_mapping, by = "scanID")
        print(paste0("snp data compiled for ", i, " of ", length(snp_id_list), " for ", trim_type))
        # Replace '3' genotypes with NA (missing genotype)
        snps[snps==3]<-NA
        snps = as.data.table(snps)
        # Get rid of snp columns which have the same entry for all entries and get ride of snp columns which have NA for all entries (missing for all individuals)
        snps_no_NA = Filter(function(x) length(unique(x))!=1, snps)
        snps_no_NA = Filter(function(x) length(unique(snps_no_NA[x != "NA"])) != 1, snps_no_NA)
        snps_no_NA2 = data.frame(localID = snps_no_NA$localID)
        for (column in names(snps_no_NA)[-c(1, ncol(snps_no_NA))]){
            if (length(unique(na.omit(snps_no_NA[[column]]))) > 1){
                snps_no_NA2 = cbind(snps_no_NA2, snps_no_NA[[column]])
                names(snps_no_NA2)[names(snps_no_NA2) == 'snps_no_NA[[column]]'] <- paste0(column)
            }
        }
        if (trim_type == "v_trim"){
            trimming_data = v_trimming
        } else if (trim_type == "d0_trim" | trim_type == "d1_trim"){
            trimming_data = d_trimming2
        } else if (trim_type == "j_trim"){
            trimming_data = j_trimming
        } else if (trim_type == "vj_insert"){
            trimming_data = vj_insert
        } else if (trim_type == "dj_insert"){
            trimming_data = dj_insert2
        } else if (trim_type == "vd_insert"){
            trimming_data = vd_insert2
        }

        if (varying_int == "FALSE"){
            weighted_regression_bootstrap_patient_vgene_results_productive = trimming_snp_regression_weighted(snps_no_NA2, trimming_data, productive = "True", trim_type = trim_type, repetitions =  repetitions)
            print("finished weighted_regression_patient_vgene_results_productive")

            weighted_regression_bootstrap_patient_vgene_results_NOT_productive = trimming_snp_regression_weighted(snps_no_NA2, trimming_data, productive = "False", trim_type = trim_type,  repetitions =  repetitions)
            print("finished weighted_regression_patient_vgene_results_NOT_productive")
            prod_name = paste0('regression_bootstrap_results/productive/', trim_type, '/', trim_type, '_productive_snplist_', length(snp_id_list),'_snps.tsv')
            not_prod_name = paste0('regression_bootstrap_results/NOT_productive/', trim_type, '/', trim_type, '_NOT_productive_snplist_', length(snp_id_list),'_snps.tsv')
        } else if (varying_int == "TRUE"){
            weighted_regression_bootstrap_patient_vgene_results_productive = trimming_snp_regression_weighted_varying_int_subject(snps_no_NA2, trimming_data, productive = "True", trim_type = trim_type, repetitions =  repetitions)
            print("finished weighted_regression_patient_vgene_results_productive")

            weighted_regression_bootstrap_patient_vgene_results_NOT_productive = trimming_snp_regression_weighted_varying_int_subject(snps_no_NA2, trimming_data, productive = "False", trim_type = trim_type,  repetitions =  repetitions)
            print("finished weighted_regression_patient_vgene_results_NOT_productive")

            prod_name = paste0('regression_bootstrap_results/productive/', trim_type, '/', trim_type, '_productive_snplist_', length(snp_id_list),'_snps_varying_intercepts_by_subject.tsv')
            not_prod_name = paste0('regression_bootstrap_results/NOT_productive/', trim_type, '/', trim_type, '_NOT_productive_snplist_', length(snp_id_list),'_snps_varying_intercepts_by_subject.tsv')
        }   
        

        if (nrow(weighted_regression_bootstrap_patient_vgene_results_productive) != 0){
            bootstrap_regression_results_v_gene_productive = rbind(bootstrap_regression_results_v_gene_productive, weighted_regression_bootstrap_patient_vgene_results_productive)
        }
        print("finished bootstrap_regression_weighted_productive")

        if (nrow(weighted_regression_bootstrap_patient_vgene_results_NOT_productive) != 0){
            bootstrap_regression_results_v_gene_NOT_productive = rbind(bootstrap_regression_results_v_gene_NOT_productive, weighted_regression_bootstrap_patient_vgene_results_NOT_productive)
        }
        print("finished bootstrap_regression_weighted_NOT_productive")

        print(paste0("finished bootstrap for snp data for ", i, " of ", length(snp_id_list), " for ", trim_type))
    }
    write.table(bootstrap_regression_results_v_gene_productive, file= prod_name, quote=FALSE, sep='\t', col.names = NA)
    write.table(bootstrap_regression_results_v_gene_NOT_productive, file= not_prod_name, quote=FALSE, sep='\t', col.names = NA)
    closefn.gds(snps_gds)
}
