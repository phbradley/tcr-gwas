library('data.table')
setDTthreads(1)
library('tidyverse')

tcr_repertoire_data_directory = paste0(PROJECT_PATH, '/tcr-gwas/_ignore/emerson_stats/')
cdr3_gene_file_path = paste0(PROJECT_PATH, '/tcr-gwas/_ignore/human_vj_allele_cdr3_nucseqs.tsv')

extract_subject_ID <- function(tcr_repertoire_file_path){
    file_name = str_split(tcr_repertoire_file_path, "/")[[1]][8]
    file_root_name = str_split(file_name, ".tsv")[[1]][1]
    localID = str_split(file_root_name, "_")[[1]][3]
    return(localID)
}

combine_genes_by_common_cdr3 <- function(cdr3_gene_file_path){
    cdr3 = fread(cdr3_gene_file_path)
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

condense_tcr_repertoire_data <- function(condensing_method, gene_type = 'NA', tcr_repertoire_data, cdr3_gene_file_path = 'NA'){
    if (condensing_method == 'by_gene_cdr3' | condensing_method == 'by_gene_allele'){
        stopifnot(gene_type != 'NA')
        if (condensing_method == 'by_gene_allele'){
            conditioning_variable = gene_type
        } else if (condensing_method == 'by_gene_cdr3'){
            stopifnot(cdr3_gene_file_path != 'NA')
            cdr3 = combine_genes_by_common_cdr3(cdr3_gene_file_path)
            tcr_repertoire_data = merge(tcr_repertoire_data, 
                                        cdr3, 
                                        by.x = gene_type, 
                                        by.y = 'id')
            conditioning_variable = 'cdr3_gene_group'
        }

        names = c('localID', paste(conditioning_variable), 'productive')
        condensed_tcr_repertoire_data = tcr_repertoire_data[, paste0(gene_type, '_count') := .N, by = mget(names)][, lapply(.SD, mean), by = mget(names), .SDcols = sapply(tcr_repertoire_data, is.numeric)]
        condensed_tcr_repertoire_data[[paste0('weighted_', gene_type, '_count')]] = condensed_tcr_repertoire_data[[paste0(gene_type, '_count')]]/nrow(tcr_repertoire_data)
    } else if (condensing_method =='by_subject'){
        names = c('localID', 'productive')
        condensed_tcr_repertoire_data = tcr_repertoire_data[, tcr_count := .N, by =mget(names)][, lapply(.SD, mean), by = mget(names), .SDcols = sapply(tcr_repertoire_data, is.numeric)]
    }
    return(condensed_tcr_repertoire_data)
}   

generate_tcr_repertoire_file_name <- function(condensing_method, infer_missing_d_gene, gene_type='NA'){
    if (condensing_method != 'by_subject'){
        stopifnot(gene_type != 'NA')
    }

    root = paste0(PROJECT_PATH, '/tcr-gwas/_ignore/condensed_tcr_repertoire_data/', condensing_method)
    middle = ifelse(condensing_method == 'by_subject', '_condensed_tcr_repertoire_data_all_subjects', paste0('_by_', gene_type,'_condensed_tcr_repertoire_data_all_subjects'))
    inferred_d_gene_end = ifelse(infer_missing_d_gene == 'True', '_with_inferred_d_gene.tsv', '_NO_inferred_d_gene.tsv')

    file_name = paste0(root, middle, inferred_d_gene_end)
    return(file_name)
}


compile_tcr_repertoire_data <- function(condensing_method, infer_missing_d_gene, tcr_repertoire_data_directory){
    stopifnot(condensing_method %in% c('by_gene_cdr3', 'by_gene_allele', 'by_subject'))
    stopifnot(infer_missing_d_gene %in% c('True', 'False'))
    stopifnot(!file.exists)
    files = list.files(path=tcr_repertoire_data_directory, 
                       pattern="*.tsv", 
                       full.names=TRUE)

    # make empty directories
    if (condensing_method == 'by_gene_cdr3' | condensing_method == 'by_gene_allele'){ 
        gene_types = c('v_gene', 'j_gene', 'd_gene') 
    } else {             
        gene_types = 'NA'
    }
    
    for (gene_type in gene_types){
        assign(paste0(gene_type, '_tcr_rep_data'), data.table())
    }

    count = 0
    for (file in files){
        count = count + 1
        file_data = fread(file)
        file_data$localID = extract_subject_ID(file)
        
        if (infer_missing_d_gene == 'True'){
            file_data = infer_d_gene(file_data)
        } else {
            file_data = file_data[d_gene != '-']
        }

        for (gene_type in gene_types){
            assign(paste0(gene_type, '_tcr_rep_data'), 
                           rbind(get(paste0(gene_type, '_tcr_rep_data')), 
                                 condense_tcr_repertoire_data(condensing_method, gene_type, file_data, cdr3_gene_file_path)))
            filename = generate_tcr_repertoire_file_name(condensing_method, infer_missing_d_gene, gene_type)
            write.table(get(paste0(gene_type, '_tcr_rep_data')), file = filename, quote=FALSE, sep='\t', col.names = NA)
        }
        print(paste0(count," of ", length(files), " processed"))
    }
    print(paste0("compiling tcr repertoire data!"))
}

