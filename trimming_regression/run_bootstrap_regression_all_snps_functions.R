library("data.table")
library("plyr")
library("readr")
library("stringr")
library("gdsfmt")
library("SNPRelate")
library("GWASTools")

source("trimming_basic_regression_functions.R")
source("trimming_regression_functions.R")
source("trimming_bootstrap_functions.R")


run_snps_trimming_snp_list <- function(snp_id_list, trim_type, condensing, gene_conditioning, weighting, repetitions, write_table){
    # Get snp meta data
    snps_gds = snpgdsOpen("../_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds")

    # get sample ids
    n <- index.gdsn(snps_gds, "sample.id")
    sampleid <- read.gdsn(n) 
    bigsize <- 35481497
    boot_regression_results_productive = data.table()
    boot_regression_results_NOT_productive = data.table()
    i = 0

    # import condensed trimming file
    if (condensing == 'by_patient'){
        assign('trimming_data', as.data.table(read.table(paste0("../_ignore/by_patient_condensed_data_all_patients.tsv"), sep = "\t", fill=TRUE, header = TRUE)[-1]))
        setnames(trimming_data, "patient_id", "localID")
    } else {
        assign('trimming_data', as.data.table(read.table(paste0("../_ignore/condensed_", trim_type, "_data_all_patients.tsv"), sep = "\t", fill=TRUE, header = TRUE)[-1]))
        setnames(trimming_data, "patient_id", "localID")
    }

    for (snp in snp_id_list){
        i = i + 1
        # get snp genotype
        genotype = snpgdsGetGeno(snps_gds, snp.id=snp, with.id = TRUE)
        # combine genotype, sample ids
        genotypes_df = data.frame(scanID = genotype$sample.id, genotype$genotype)
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
        snps_no_NA2 = data.frame(snps$localID, snps[[paste0("snp",snp)]])
        colnames(snps_no_NA2) = c("localID", paste0("snp",snp))

        # skip iteration if the genotypes are all the same...
        if (length(unique(snps_no_NA2$snp)[!is.na(unique(snps_no_NA2$snp))]) <= 1){
            next
        }

        # do regression, bootstrap
        if (condensing == "by_gene"){
            regression_productive = trimming_regression(snps_dataframe = snps_no_NA2, condensed_trimming_dataframe = trimming_data, productive = "True",trim_type = trim_type, bootstrap_repetitions = repetitions, gene_conditioning, weighting)
            print("finished regression_productive")
            regression_NOT_productive = trimming_regression(snps_dataframe = snps_no_NA2, condensed_trimming_dataframe = trimming_data, productive = "False", trim_type = trim_type, bootstrap_repetitions = repetitions, gene_conditioning, weighting)
            print("finished regression_NOT_productive")
             # set path name
            prod_name = generate_file_name(snp_id_list, trim_type, productivity = 'True', gene_conditioning, weighting, condensing, repetitions)
            not_prod_name = generate_file_name(snp_id_list, trim_type, productivity = 'False', gene_conditioning, weighting, condensing, repetitions)
        } else if (condensing == 'by_patient'){
            regression_productive = simple_trimming_snp_regression(snps_no_NA2, trimming_data, productive = "True", trim_type =trim_type, repetitions =  repetitions, python_test = 'True')
            print("finished simple_regression_productive")
            regression_NOT_productive = simple_trimming_snp_regression(snps_no_NA2, trimming_data, productive = "False", trim_type = trim_type,  repetitions =  repetitions, python_test = 'True')
            print("finished simple_regression_NOT_productive")
             # set path name
            prod_name = paste0('regression_bootstrap_results/productive/', trim_type, '/', trim_type, '_productive_snplist_', snp_id_list[1], "-",snp_id_list[length(snp_id_list)], '_snps_simple_condensing_', condensing, '.tsv')
            not_prod_name = paste0('regression_bootstrap_results/NOT_productive/', trim_type, '/', trim_type, '_NOT_productive_snplist_', snp_id_list[1], "-",snp_id_list[length(snp_id_list)],'_snps_simple_condensing_', condensing, '.tsv') 
        } else if (condensing == 'none'){
            regression_productive = simple_trimming_snp_regression_no_condensing(snps_no_NA2, productive = "True", trim_type = trim_type)
            print("finished simple_regression_productive")
            regression_NOT_productive = simple_trimming_snp_regression_no_condensing(snps_no_NA2, productive = "False", trim_type = trim_type)
            print("finished simple_regression_NOT_productive")
             # set path name
            prod_name = paste0('regression_bootstrap_results/productive/', trim_type, '/', trim_type, '_productive_snplist_', snp_id_list[1], "-",snp_id_list[length(snp_id_list)], '_snps_simple_condensing_', condensing, '.tsv')
            not_prod_name = paste0('regression_bootstrap_results/NOT_productive/', trim_type, '/', trim_type, '_NOT_productive_snplist_', snp_id_list[1], "-",snp_id_list[length(snp_id_list)],'_snps_simple_condensing_', condensing, '.tsv')  
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

        if (write_table == "True"){
            # Write tables
            write.table(boot_regression_results_productive, file= prod_name, quote=FALSE, sep='\t', col.names = NA)
            write.table(boot_regression_results_NOT_productive, file= not_prod_name, quote=FALSE, sep='\t', col.names = NA)
        } 
        print(paste0("finished bootstrap for snp data for ", i, " of ", length(snp_id_list), " for ", trim_type))
    } 

    if (write_table != "True"){
        boot_regression_results_productive$productivity_status = "productive"
        boot_regression_results_NOT_productive$productivity_status = "NOT_productive"
        results = rbind(boot_regression_results_productive, boot_regression_results_NOT_productive)
        return(results)
    }
    
    closefn.gds(snps_gds)
}

test_gene_conditioning <- function(snp_id_list, trim_type, condensing, gene_conditioning, weighting, repetitions, write_table){
    # Get snp meta data
    snps_gds = snpgdsOpen("../_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds")

    # get sample ids
    n <- index.gdsn(snps_gds, "sample.id")
    sampleid <- read.gdsn(n) 
    bigsize <- 35481497
    boot_regression_results_productive = data.table()
    boot_regression_results_NOT_productive = data.table()
    i = 0

    # combine all trimming data (from all gene types)
    trimming_data_all = data.table()
    for (trim in c('v_trim', 'd0_trim', 'j_trim')){
        assign(paste0(trim, '_trimming_data'), as.data.table(read.table(paste0("../_ignore/condensed_", trim, "_data_all_patients.tsv"), sep = "\t", fill=TRUE, header = TRUE)[-1]))
        setnames(get(paste0(trim, '_trimming_data')), "patient_id", "localID")
        setnames(get(paste0(trim, '_trimming_data')), paste0(substring(trim, 1, 1), '_gene'), "gene")
        setnames(get(paste0(trim, '_trimming_data')), paste0(substring(trim, 1, 1), '_gene_count'), "gene_count")
        setnames(get(paste0(trim, '_trimming_data')), paste0('weighted_', substring(trim, 1, 1), '_gene_count'), "weighted_gene_count")
        gene_type = data.table(gene_type = rep(trim, nrow(get(paste0(trim, '_trimming_data')))))
        assign(paste0(trim, '_trimming_data'), cbind(get(paste0(trim, '_trimming_data')), gene_type))
        trimming_data_all = rbind(trimming_data_all, get(paste0(trim, '_trimming_data')))
    }
    
}
