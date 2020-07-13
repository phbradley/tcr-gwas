library("data.table")
library("plyr")
library("readr")
library("stringr")
library("gdsfmt")
library("SNPRelate")
library("GWASTools")

source("trimming_basic_regression_functions.R")
source("trimming_bootstrap_functions.R")


v_trimming = as.data.table(read.table("../_ignore/condensed_v_trim_data_all_patients.tsv", sep = "\t", fill=TRUE, header = TRUE)[-1])
colnames(v_trimming) = c("localID", "v_gene", "productive", "v_trim", "v_gene_count", "weighted_v_gene_count")

d_trimming = as.data.table(read.table("../_ignore/condensed_d0_trim_data_all_patients.tsv", sep = "\t", fill=TRUE, header = TRUE)[-1])
colnames(d_trimming) = c("localID", "d_gene", "productive", "d0_trim", "d1_trim", "d_gene_count", "weighted_d_gene_count")

j_trimming = as.data.table(read.table("../_ignore/condensed_j_trim_data_all_patients.tsv", sep = "\t", fill=TRUE, header = TRUE)[-1])
colnames(j_trimming) = c("localID", "j_gene", "productive", "j_trim", "j_gene_count", "weighted_j_gene_count")

vj_insert = as.data.table(read.table("../_ignore/condensed_vj_insert_data_all_patients.tsv", sep = "\t", fill=TRUE, header = TRUE)[-1])
colnames(vj_insert) = c("localID", "v_gene", "j_gene", "productive", "vj_insert", "vj_gene_count", "weighted_vj_gene_count")

dj_insert = as.data.table(read.table("../_ignore/condensed_dj_insert_data_all_patients.tsv", sep = "\t", fill=TRUE, header = TRUE)[-1])
colnames(dj_insert) = c("localID", "d_gene", "j_gene", "productive", "dj_insert", "dj_gene_count", "weighted_dj_gene_count")

vd_insert = as.data.table(read.table("../_ignore/condensed_vd_insert_data_all_patients.tsv", sep = "\t", fill=TRUE, header = TRUE)[-1])
colnames(vd_insert) = c("localID", "v_gene", "d_gene", "productive", "vd_insert", "vd_gene_count", "weighted_vd_gene_count")



run_snps_trimming_snp_list <- function(snp_id_list, trim_type, varying_int, repetitions){
    # Get snp meta data
    snps_gds = snpgdsOpen("../_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds")

    # get sample ids
    n <- index.gdsn(snps_gds, "sample.id")
    sampleid <- read.gdsn(n) 
    bigsize <- 35481497
    boot_regression_results_productive = data.table()
    boot_regression_results_NOT_productive = data.table()
    i = 0

    for (snp in snp_id_list){
        i = i + 1
        # get snp genotype
        genotype = snpgdsGetGeno(snps_gds, snp.id=snp)
        # combine genotype, sample ids
        genotypes_df = data.frame(scanID = c(sampleid), genotype)
        colnames(genotypes_df) = c("scanID", paste0("snp",snp))
        genotypes_dt = as.data.table(genotypes_df)
        genotypes_dt$scanID = as.numeric(as.character(genotypes_dt$scanID))

        # Convert subject names : 
        subject_id_mapping = as.data.table(read.table('../_ignore/snp_data/gwas_id_mapping.tsv', sep = "\t", fill=TRUE, header = TRUE, check.names = FALSE))
        snps = merge(genotypes_dt, subject_id_mapping, by = "scanID")

        print(paste0("snp data compiled for ", i, " of ", length(snp_id_list), " for ", trim_type))

        # Replace '3' genotypes with NA (missing genotype)
        snps[snps==3]<-NA
        snps = as.data.table(snps)

        # Get rid of snp columns which have the same entry for all entries and get rid of snp columns which have NA for all entries (missing for all individuals)
        snps_no_NA = Filter(function(x) length(unique(x))!=1, snps)
        snps_no_NA = Filter(function(x) length(unique(snps_no_NA[x != "NA"])) != 1, snps_no_NA)
        snps_no_NA2 = data.frame(snps_no_NA$localID, snps_no_NA[[paste0("snp",snp)]])
        colnames(snps_no_NA2) = c("localID", paste0("snp",snp))

        # import condensed trimming data
        if (trim_type == "v_trim"){
            trimming_data = v_trimming
        } else if (trim_type == "d0_trim" | trim_type == "d1_trim"){
            trimming_data = d_trimming
        } else if (trim_type == "j_trim"){
            trimming_data = j_trimming
        } else if (trim_type == "vj_insert"){
            trimming_data = vj_insert
        } else if (trim_type == "dj_insert"){
            trimming_data = dj_insert
        } else if (trim_type == "vd_insert"){
            trimming_data = vd_insert
        }

        # do regression, bootstrap
        if (varying_int == "True"){
            regression_productive = trimming_snp_regression_weighted_varying_int_subject(snps_no_NA2, trimming_data, productive = "True", trim_type =trim_type, repetitions =  repetitions)
            print("finished regression_productive")

            regression_NOT_productive = trimming_snp_regression_weighted_varying_int_subject(snps_no_NA2, trimming_data, productive = "False", trim_type = trim_type,  repetitions =  repetitions)
            print("finished regression_NOT_productive")

             # set path name
            prod_name = paste0('regression_bootstrap_results/productive/', trim_type, '/', trim_type, '_productive_snplist_', snp_id_list[1], "-",snp_id_list[length(snp_id_list)], '_snps_varying_intercepts_by_subject_with_gene.tsv')
            not_prod_name = paste0('regression_bootstrap_results/NOT_productive/', trim_type, '/', trim_type, '_NOT_productive_snplist_', snp_id_list[1], "-",snp_id_list[length(snp_id_list)],'_snps_varying_intercepts_by_subject_with_gene.tsv')  

        } else {
            regression_productive = simple_trimming_snp_regression(snps_no_NA2, trimming_data, productive = "True", trim_type =trim_type, repetitions =  repetitions)
            print("finished simple_regression_productive")

            regression_NOT_productive = simple_trimming_snp_regression(snps_no_NA2, trimming_data, productive = "False", trim_type = trim_type,  repetitions =  repetitions)
            print("finished simple_regression_NOT_productive")

             # set path name
            prod_name = paste0('regression_bootstrap_results/productive/', trim_type, '/', trim_type, '_productive_snplist_', snp_id_list[1], "-",snp_id_list[length(snp_id_list)], '_snps_simple.tsv')
            not_prod_name = paste0('regression_bootstrap_results/NOT_productive/', trim_type, '/', trim_type, '_NOT_productive_snplist_', snp_id_list[1], "-",snp_id_list[length(snp_id_list)],'_snps_simple.tsv')  
        }
        
        
        # if regression, bootstrap conducted, add them to the results tables
        if (nrow(regression_productive) != 0){
            boot_regression_results_productive = rbind(boot_regression_results_productive, regression_productive)
        }
        print("finished bootstrap_regression_productive")

        if (nrow(regression_NOT_productive) != 0){
            boot_regression_results_NOT_productive = rbind(boot_regression_results_NOT_productive, regression_NOT_productive)
        }
        print("finished bootstrap_regression_NOT_productive")

        print(paste0("finished bootstrap for snp data for ", i, " of ", length(snp_id_list), " for ", trim_type))
    } 

    # Write tables
    write.table(boot_regression_results_productive, file= prod_name, quote=FALSE, sep='\t', col.names = NA)
    write.table(boot_regression_results_NOT_productive, file= not_prod_name, quote=FALSE, sep='\t', col.names = NA)

    closefn.gds(snps_gds)
}
