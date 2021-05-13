source('config/config.R')

compile_all_genotypes <- function(snp_file){
    file = fread(snp_file)
    file[rs12768894 == 3, rs12768894 := NA]
    file[rs3762093 == 3, rs3762093 := NA]
    if ('id' %in% colnames(file)){
        setnames(file, 'id', 'localID')
    }
    colnames(file) = c('localID', '12768894', '3762093')
    return(file)
}

extract_subject_ID <- function(tcr_repertoire_file_path){
    file_name = str_split(tcr_repertoire_file_path, "/")[[1]][length(str_split(tcr_repertoire_file_path, "/")[[1]])]
    file_root_name = str_split(file_name, ".tsv")[[1]][1]
    localID = str_split(file_root_name, "_run4_umi5")[[1]][1]
    return(localID)
}

combine_genes_by_common_cdr3 <- function(){
    if (CHAIN == 'beta'){
        cdr3 = fread(CDR3_GENE_ASSIGNMENT_FILE)
        cdr3$gene = str_split(cdr3$id, fixed('*'), simplify = TRUE)[,1]

        if (!('D' %in% substr(cdr3$gene, 4, 4))){
            cdr3 = rbind(cdr3, data.frame(cdr3_nucseq = c('NA', 'NA', 'NA'), id = c("TRBD1*01", "TRBD2*01", "TRBD2*02"), gene = c("TRBD1*01", "TRBD2*01", "TRBD2*02")))
        }
    
        cdr3 = cdr3[,cdr3_gene_group := .GRP, by = .(cdr3_nucseq, gene)]
        cdr3_groups = cdr3[,c('id', 'cdr3_gene_group')]
    } else if (CHAIN == 'alpha'){
        cdr3 = fread(CDR3_GENE_ASSIGNMENT_FILE)
        cdr3$id = cdr3$cdr3_gene_group
        cdr3_groups = cdr3
    }
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

source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/scripts/phenotype_functions/', PHENOTYPE, '.R'))

PCA_COUNT <<- 3

generate_condensed_tcr_repertoire_file_name <- function(){
    inferred_d_gene_end = ifelse(KEEP_MISSING_D_GENE == 'True', '_with_inferred_d_gene.tsv', '_NO_inferred_d_gene.tsv')
    condensed_tcr_repertoires_file_name = paste0(OUTPUT_PATH, '/condensed_nicaraguan_tcr_repertoire_data/', CHAIN, '_chain_', CONDENSING_VARIABLE, '_by_', GENE_TYPE, '_condensed_tcr_repertoire_data_all_subjects_', PHENOTYPE_CLASS, '_for_', PHENOTYPE, inferred_d_gene_end)
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

compile_population_structure_pca_data <- function(){
    pca_file = fread(PCA_FILE)
    pca_subset = pca_file[,1:4]
    colnames(pca_subset) = c('localID', paste0('EV', seq(1,PCA_COUNT)))
    return(pca_subset)
}


compile_phenotype_data <- function(){
    tcr_repertoire_data = compile_condensed_tcr_repertoire_data()
    if (!is.na(PCA_COUNT)){
        pca_data = compile_population_structure_pca_data()  
        phenotype_data = merge(tcr_repertoire_data, pca_data, by = 'localID')
    } else {
        phenotype_data = tcr_repertoire_data
    } 
    return(phenotype_data)
}

########################
# regression functions #
########################

make_regression_file_path <- function(){
    path_to_file = paste0(OUTPUT_PATH, '/results/validation_cohort/', PHENOTYPE, '_', CONDENSING_VARIABLE, '/D_infer_', KEEP_MISSING_D_GENE, '/', PCA_COUNT, '_PCAir_PCs/')

    if (!dir.exists(path_to_file)){
        system(paste0('mkdir -p ', path_to_file))
    }

    return(path_to_file)
}

find_regression_file_path_for_shell <- function(){
    cat(make_regression_file_path())
}

make_regression_file_name <- function(){
    path = make_regression_file_path()
    file_name = paste0(PHENOTYPE, '_regressions_condensing_', CONDENSING_VARIABLE, '_', CHAIN, '.tsv')
    return(paste0(path, file_name)) 
}

execute_regressions <- function(genotypes, phenotypes, write.table){
    regression_data = merge(genotypes, phenotypes, by = 'localID')
    if (nrow(genotypes) == 0){
        results = data.table(snp = NA, slope = NA, standard_error = NA, pvalue = NA, parameter = NA, phenotype = NA, bootstraps = NA, productive = NA)
    } else {
        snps = colnames(genotypes)[colnames(genotypes)!='localID'] 
        count = 0 
        registerDoParallel(cores=NCPU)
        results = foreach(snp = snps, .combine='rbind') %dopar% {
            regress(snp, snps, regression_data)
        }
        stopImplicitCluster()
    }

    results$snp = paste0('rs', results$snp)

    if (write.table == TRUE){
        file_name = make_regression_file_name()
        write.table(as.data.frame(results), file = file_name, quote=FALSE, sep='\t', col.names = NA)
    } else if (write.table == FALSE){
        return(results)
    }   
}


source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/scripts/phenotype_functions/phenotype_classes/', PHENOTYPE_CLASS, '_class_functions.R'))
