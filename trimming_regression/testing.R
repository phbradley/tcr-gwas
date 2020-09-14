source("/home/mrussel2/tcr-gwas/trimming_regression/run_bootstrap_regression_all_snps_functions_cluster.R")
#source("find_correlated_snps.R")
source("/home/mrussel2/tcr-gwas/trimming_regression/import_snp_file.R")


args = c(1001, 'v_trim', 4)


count = 1000
# Read in snp list
snp_data = get_snp_file(snp_start = as.numeric(args[1]), count)
genotype_data = get_genotype_file(snp_start = as.numeric(args[1]), count)
snp_data = snp_data[1:100,]
genotype_data = genotype_data[,1:100]
genotype_data_filtered = remove_matrix_column_by_genotype(genotype_data)


snp_list = snp_data
genotype_list = genotype_data_filtered
trim_type = args[2]
 gene_type = 'same'
 condensing = 'by_gene'
  gene_conditioning = 'True'
   weighting = 'True'
   random_effects = 'True'
    repetitions = 4
    write_table = 'True'
    ncpus = as.numeric(args[3])

    snp = 1003

regression_dataframe = data.frame()

    # import condensed trimming file
    if (condensing == 'by_patient'){
        assign('trimming_data', as.data.frame(read.table(paste0("/home/mrussel2/tcr-gwas/_ignore/by_patient_condensed_data_all_patients.tsv"), sep = "\t", fill=TRUE, header = TRUE)[-1]))
        names(trimming_data)[names(trimming_data) == "patient_id"] <- "localID"
    } else if (condensing == 'gene_cross'){
        trimming_data = compile_trimming_data_cross()
        trimming_data = trimming_data %>% filter(gene_class == gene_type)
        names(trimming_data)[names(trimming_data) == 'weighted_gene_count'] <- paste0('weighted_', gene_type, '_count')
        names(trimming_data)[names(trimming_data) == 'gene'] <- paste0(gene_type)
    } else if (condensing == 'phil'){
        assign('trimming_data', as.data.frame(read.table(paste0("/home/mrussel2/tcr-gwas/_ignore/by_patient_condensed_data_all_patients_phil.tsv"), sep = "\t", fill=TRUE, header = TRUE)[-1]))
    } else {
        assign('trimming_data', as.data.frame(read.table(paste0("/home/mrussel2/tcr-gwas/_ignore/condensed_", trim_type, "_data_all_patients.tsv"), sep = "\t", fill=TRUE, header = TRUE)[-1]))
        names(trimming_data)[names(trimming_data) == "patient_id"] <- "localID"
    }

genotypes_temp = as.data.frame(genotype_list[, as.character(snp)])
    snp_genotypes = data.frame(rownames(genotypes_temp),genotypes_temp)
    colnames(snp_genotypes) = c('localID', paste0('snp',snp))
    rownames(snp_genotypes) = NULL

    index = which(as.numeric(colnames(genotype_list)) == snp)

    # skip iteration if the genotypes are all the same...
    if (length(unique(snp_genotypes$snp)[!is.na(unique(snp_genotypes$snp))]) <= 1){
        temp_regression_dataframe = merge(snp_list, data.frame(snp = snp, intercept = 'NA', slope = 'NA', standard_error = 'NA', pvalue = 'NA', productivity = 'NA'), by = 'snp')
        temp_regression_dataframe$snp = paste0('snp', temp_regression_dataframe$snp)
        regression_dataframe = rbind(regression_dataframe, temp_regression_dataframe)
        print(paste0("no regression needed for snp data for ", index, " of ", ncol(genotype_list), " for ", trim_type))
    }
    
    snps_dataframe = snp_genotypes
    condensed_trimming_dataframe = trimming_data
    productive = "True"
    trim_type = trim_type
    gene_type = gene_type
    bootstrap_repetitions = repetitions
    gene_conditioning

   
   data = snps_trimming_data[snp != "NA"]
   cluster_variable = snps_trimming_data[snp != "NA"]$localID
  varying_int = 'True'
  weighting
  repetitions = bootstrap_repetitions