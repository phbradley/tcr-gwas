source('config/config.R')
###############################################
# Functions to read in SNP and genotype files #
###############################################

map_scanID_to_localID <- function(scanIDs_to_convert){
    ID_map_file = fread(ID_MAPPING_FILE)
    converted_IDs = plyr::mapvalues(scanIDs_to_convert, 
                                    ID_map_file$scanID, 
                                    ID_map_file$localID)
    return(converted_IDs)
}

create_maf_file <- function(){
    snp_gds_file = openfn.gds(SNP_GDS_FILE, readonly = TRUE, allow.fork = TRUE)
    bigsize = 35481497
    snp_starts = seq(1, bigsize, by = 10000)
    mafs = data.table()
    registerDoParallel(cores=NCPU)
    mafs = foreach(start = snp_starts, .combine = 'rbind') %dopar% {
        numrows = min(10000, bigsize-start+1)
        snp_ids = read.gdsn(index.gdsn(snp_gds_file, "snp.id"),
                                       start=start,
                                       count=numrows)
        genotype_matrix = read.gdsn(index.gdsn(snp_gds_file, "genotype"),
                                    start=c(1,start),
                                    count=c(398, numrows))
        colnames(genotype_matrix) = snp_ids
        genotype_matrix[genotype_matrix == 3] <- NA 
        genotype_dt = as.data.table(genotype_matrix)
        subject_counts = colSums(!is.na(genotype_dt))
        nonNA_subject_counts = subject_counts[subject_counts != 0]
        cols = names(nonNA_subject_counts)
        allele_counts = colSums(genotype_dt[, ..cols], na.rm = TRUE)
        temp = data.table(snp = names(nonNA_subject_counts), maf = allele_counts/(2*nonNA_subject_counts))
        temp[maf >= 0.5, maf_flipped := TRUE]
        temp[maf < 0.5,  maf_flipped := FALSE]
        temp[maf >= 0.5, maf := 1-maf]
        print(paste0('finished processing mafs for snps ', start, ' to ', start + 10000))
        temp
    }
    stopImplicitCluster()
    closefn.gds(snp_gds_file)
    fwrite(mafs, paste0(OUTPUT_PATH, '/maf_all_snps.tsv'), sep = '\t')
    return(mafs)
}


filter_snps_by_maf <- function(genotype_matrix){
    file_name = paste0(OUTPUT_PATH, '/maf_all_snps.tsv')
    if (!file.exists(file_name)){
        create_maf_file()
    }
    maf_data = fread(file_name)
    maf_data_filtered = maf_data[maf >= MAF_CUTOFF]
    list_of_snps = as.numeric(intersect(colnames(genotype_matrix), maf_data_filtered$snp))
    return(list_of_snps)
}
    
snp_file_by_snp_start <- function(snp_start, count){
    snp_gds_file = openfn.gds(SNP_GDS_FILE)
    snp_id = read.gdsn(index.gdsn(snp_gds_file, "snp.id"))

    count = ifelse((snp_start+count) > length(snp_id), length(snp_id) - snp_start, count)
    snp_chrom = read.gdsn(index.gdsn(snp_gds_file, "snp.chromosome"),
                           start = snp_start,
                           count=count)
    snp_pos = read.gdsn(index.gdsn(snp_gds_file, "snp.position"),
                        start = snp_start,
                        count=count)
    snps = data.frame(snp = snp_id[snp_start:(snp_start+count-1)],
                      chr = snp_chrom, 
                      hg19_pos = snp_pos)
    closefn.gds(snp_gds_file)
    return(snps)
}


remove_matrix_column_by_genotype <- function(genotype_matrix){
    snps_passing_maf_cutoff = filter_snps_by_maf(genotype_matrix)
    if (length(snps_passing_maf_cutoff) == 0){
        return(data.table())
    } else {
        for (snp in colnames(genotype_matrix)){
            genotypes = unique(genotype_matrix[,snp])
            nonNA_genotypes = genotypes[genotypes != 3]
            if (length(nonNA_genotypes) <= 1 | !(snp %in% snps_passing_maf_cutoff)){
                genotype_matrix = genotype_matrix[, colnames(genotype_matrix) != snp, drop = FALSE]
            }
        }
        genotype_matrix[genotype_matrix == 3] <- NA
        return(genotype_matrix)
    }
}

compile_all_genotypes <- function(snp_start, count){
    snp_gds_file = openfn.gds(SNP_GDS_FILE)
    bigsize = 35481497
    numrows = min(count, bigsize-snp_start+1)
    
    genotype_matrix = read.gdsn(index.gdsn(snp_gds_file, "genotype"),
                                 start=c(1,snp_start),
                                 count=c(398, numrows))
    sample_ids = read.gdsn(index.gdsn(snp_gds_file, "sample.id"),
                            start=1,
                            count=398)
    snp_ids = read.gdsn(index.gdsn(snp_gds_file, "snp.id"),
                         start=snp_start,
                         count=numrows)
    closefn.gds(snp_gds_file)
    
    rownames(genotype_matrix) = map_scanID_to_localID(sample_ids)
    colnames(genotype_matrix) = snp_ids
    genotype_matrix = remove_matrix_column_by_genotype(genotype_matrix)
    colnames(genotype_matrix) = as.character(colnames(genotype_matrix))
    genotype_dt = data.table(localID = row.names(genotype_matrix), genotype_matrix)
    return(genotype_dt)
}

#############################################
# Functions to condense TCR repertoire data #
#############################################

extract_subject_ID <- function(tcr_repertoire_file_path){
    file_name = str_split(tcr_repertoire_file_path, "/")[[1]][8]
    file_root_name = str_split(file_name, ".tsv")[[1]][1]
    localID = str_split(file_root_name, "_")[[1]][3]
    return(localID)
}

combine_genes_by_common_cdr3 <- function(){
    cdr3 = fread(CDR3_GENE_ASSIGNMENT_FILE)
    cdr3$gene = str_split(cdr3$id, fixed('*'), simplify = TRUE)[,1]

    if (!('D' %in% substr(cdr3$gene, 4, 4))){
        cdr3 = rbind(cdr3, data.frame(cdr3_nucseq = c('NA', 'NA', 'NA'), id = c("TRBD1*01", "TRBD2*01", "TRBD2*02"), gene = c("TRBD1*01", "TRBD2*01", "TRBD2*02")))
    }
    
    cdr3 = cdr3[,cdr3_gene_group := .GRP, by = .(cdr3_nucseq, gene)]
    cdr3_groups = cdr3[,c('id', 'cdr3_gene_group')]
    return(cdr3_groups)
}     

infer_d_gene <- function(tcr_repertoire_data){
    tcrs_with_d_gene = tcr_repertoire_data[d_gene != '-']
    tcrs_missing_d_gene = tcr_repertoire_data[d_gene == '-']

    triplet_counts = tcrs_with_d_gene[,.N, by = .(v_gene, d_gene, j_gene)]
    doublet_dj = tcrs_with_d_gene[,.N, by = .(j_gene, d_gene)]
    doublet_vd = tcrs_with_d_gene[,.N, by = .(v_gene, d_gene)]

    for (row in 1:nrow(tcrs_missing_d_gene)){
        observed_v = tcrs_missing_d_gene[row]$v_gene
        observed_j = tcrs_missing_d_gene[row]$j_gene

        empirical_vj = triplet_counts[v_gene == observed_v & j_gene == observed_j]
        empirical_j = doublet_dj[j_gene == observed_j]
        empirical_v = doublet_vd[v_gene == observed_v]

        if(nrow(empirical_vj) != 0){
            empirical = empirical_vj
        } else if (nrow(empirical_j) != 0){
            empirical = empirical_j
        } else {
            empirical = empirical_v
        }

        empirical$probability = empirical$N/sum(empirical$N)
        sample_from_empirical = rmultinom(1,1,empirical$probability)
        index = which(sample_from_empirical == 1)
        tcrs_missing_d_gene[row]$d_gene = empirical[index]$d_gene

        for (variable in c('d0_trim', 'd1_trim', 'vd_insert', 'dj_insert')){
            if (tcrs_missing_d_gene[row]$d_gene == "TRBD2*01" | tcrs_missing_d_gene[row]$d_gene == "TRBD2*02"){
                tcrs_missing_d_gene[row][[variable]] = 4
            } else if (tcrs_missing_d_gene[row]$d_gene == "TRBD1*01"){
                tcrs_missing_d_gene[row][[variable]] = 3
            }
        }
    }
    inferred_tcrs = rbind(tcrs_missing_d_gene, tcrs_with_d_gene)
    return(inferred_tcrs)
}

#######################################
# Import phenotype specific functions #
#######################################

source(paste0(PROJECT_PATH, '/tcr-gwas/analysis_scripts/conditional_analysis_scripts/phenotype_functions/', PHENOTYPE, '.R'))


generate_condensed_tcr_repertoire_file_name <- function(){
    inferred_d_gene_end = ifelse(KEEP_MISSING_D_GENE == 'True', '_with_inferred_d_gene.tsv', '_NO_inferred_d_gene.tsv')
    condensed_tcr_repertoires_file_name = paste0(OUTPUT_PATH, '/condensed_tcr_repertoire_data/', CONDENSING_VARIABLE, '_by_', GENE_TYPE, '_condensed_tcr_repertoire_data_all_subjects_', PHENOTYPE_CLASS, '_for_', PHENOTYPE, inferred_d_gene_end)
    return(condensed_tcr_repertoires_file_name)
}

compile_condensed_tcr_repertoire_data <- function(){
    file_name = generate_condensed_tcr_repertoire_file_name() 

    if (!file.exists(file_name)){
        print('Data CONDENSING required. Computing now.')
        condense_all_tcr_repertoire_data()
    }

    tcr_repertoire_data = fread(file = file_name)
    
    return(tcr_repertoire_data)
}

filter_tcr_repertoire_data_by_productivity <- function(condensed_tcr_repertoire_data, productivity){
    if (productivity == 'productive'){
        condensed_tcr_repertoire_data = condensed_tcr_repertoire_data[productive == 'TRUE']
    } else if (productivity == 'not_productive'){
        condensed_tcr_repertoire_data = condensed_tcr_repertoire_data[productive == 'FALSE']
    }
    return(condensed_tcr_repertoire_data)
}

remove_small_repertoire_observations <- function(tcr_repertoire_data, productive_log10_count_cutoff = 4.25, NOT_productive_log10_count_cutoff = 3.5){
    filtered_tcr_repertoire_data = tcr_repertoire_data[(productive == 'TRUE' & log10(productivity_tcr_count)>productive_log10_count_cutoff) | (productive == 'FALSE' & log10(productivity_tcr_count)>NOT_productive_log10_count_cutoff)]
    return(filtered_tcr_repertoire_data)
}

compile_population_structure_pca_data <- function(){
    pca_file = fread(PCA_FILE)
    pca_file$localID = map_scanID_to_localID(pca_file$scanID)
    pca_columns = c('localID', paste0('EV', seq(1,PCA_COUNT)))
    pca_subset = pca_file[,..pca_columns]
    return(pca_subset)
}

compile_d_allele_status_correction_data <- function(){
    allele_statuses = fread(D_ALLELES)
    allele_statuses[allele_0 == 'TRBD2*02', alt_allele_genotype := 1]
    allele_statuses[allele_0 != 'TRBD2*02', alt_allele_genotype := 0]
    allele_statuses[allele_1 == 'TRBD2*02', alt_allele_genotype := alt_allele_genotype + 1]
    colnames(allele_statuses) = c('allele_0', 'allele_1', 'localID', 'TRBD2_alt_allele_genotype')
    return(allele_statuses)
}

compile_phenotype_data <- function(){
    tcr_repertoire_data = compile_condensed_tcr_repertoire_data()
    if (PHENOTYPE != 'tcr_div'){
        if (!('productivity_tcr_count' %in% colnames(tcr_repertoire_data))){
            tcr_repertoire_data$productivity_tcr_count = tcr_repertoire_data$tcr_count
        }
        tcr_repertoire_data = remove_small_repertoire_observations(tcr_repertoire_data)
    }
    if (!is.na(PCA_COUNT)){
        pca_data = compile_population_structure_pca_data()  
        phenotype_data = merge(tcr_repertoire_data, pca_data, by = 'localID')
    } else {
        phenotype_data = tcr_repertoire_data
    } 

    if (!is.na(ALLELE_STATUS_CORRECTION)){
        allele_status_data = compile_d_allele_status_correction_data()
        phenotype_data = merge(phenotype_data, allele_status_data, by = 'localID')
    }
    return(phenotype_data)
}

########################
# regression functions #
########################

make_regression_file_path <- function(){
    path_to_file = paste0(OUTPUT_PATH, '/cluster_job_results/conditional_analysis/', PHENOTYPE, '_', CONDENSING_VARIABLE, '/D_infer_', KEEP_MISSING_D_GENE, '/', PCA_COUNT, '_PCAir_PCs/')

    if (!dir.exists(path_to_file)){
        system(paste0('mkdir -p ', path_to_file))
    }

    return(path_to_file)
}

find_regression_file_path_for_shell <- function(){
    cat(make_regression_file_path())
}

make_regression_file_name <- function(productivity, round_number){
    path = make_regression_file_path()
    file_name = paste0(PHENOTYPE, '_', productivity, '_regressions_condensing_', CONDENSING_VARIABLE, '_', GENE, '_conditional_analysis_ROUND_', round_number, '.tsv')
    return(paste0(path, file_name)) 
}

make_regression_file_name_final <- function(){
    path = paste0(OUTPUT_PATH, '/results/conditional_results/') 
    file_name = paste0(PHENOTYPE, '_regressions_condensing_', CONDENSING_VARIABLE, '_', GENE, '_conditional_analysis.tsv')
    return(paste0(path, file_name)) 
}

execute_conditional_regressions <- function(genotypes, phenotypes, conditional_snp_start, productivity, significance_threshold){
    regression_data = merge(genotypes, phenotypes, by = 'localID')
    conditional_snp_list = conditional_snp_start
    still_finding_snps = TRUE
    round_number = 0
    while (still_finding_snps == TRUE){
        round_number = round_number + 1
        snps = colnames(genotypes)[colnames(genotypes)!='localID'] 
        snps = snps[!(snps %in% conditional_snp_list)]
        registerDoParallel(cores=NCPU)
        results = foreach(snp = snps, .combine='rbind') %dopar% {
            conditional_regress(snp, snps, regression_data, conditional_snp_list, productivity)
        }
        stopImplicitCluster()
        significant_signals = results[pvalue < significance_threshold]
        for (cond_snp in seq(1, length(conditional_snp_list))){
            significant_signals = significant_signals[get(paste0('conditional_snp_', cond_snp, '_pvalue')) < significance_threshold]
        }

        file_name = make_regression_file_name(productivity, round_number)
        write.table(as.data.frame(results), file = file_name, quote=FALSE, sep='\t', col.names = NA)

        if (nrow(significant_signals) == 0){
            still_finding_snps = FALSE
        } else {
            significant_signal_snp = significant_signals[order(pvalue)][1]$snp
            conditional_snp_list = c(conditional_snp_list, significant_signal_snp)
        }
    }
}

write_final_compiled_file <- function(significance_threshold, original_analysis_data){
    files = fs::dir_ls(path = make_regression_file_path())
    together = data.table()
    for (file in files){
        temp_file = fread(file)
        together = rbind(together, temp_file, fill = TRUE)
    }
    file_name = make_regression_file_name_final()
    sig = together[pvalue < significance_threshold]

    if (nrow(sig) != 0){
        original = original_analysis_data[snp %in% unique(sig$conditional_snps)]
        sig = rbind(sig, original, fill = TRUE)
    }

    fwrite(sig, file_name, sep = '\t')
}
