library("data.table")
library("gdsfmt")
library("plyr")
library("readr")
library("stringr")

#need to comment out the next two lines
#args = commandArgs(trailingOnly=TRUE)
#project_path = '/home/mrussel2'

infer_d_gene <- function(patient_trimming_table){
    # filter between entries that have missing d_gene and those that do not have a missing dgene
    patient_trimming_table_NO_missing_d = patient_trimming_table[d_gene != '-']
    patient_trimming_table_missing_d = patient_trimming_table[d_gene == '-']

    # condense and count triplet observations (i.e. VDJ combos)
    triplet_counts = patient_trimming_table_NO_missing_d[,.N, by = .(v_gene, d_gene, j_gene)]
    doublet_dj = patient_trimming_table_NO_missing_d[,.N, by = .(j_gene, d_gene)]
    doublet_vd = patient_trimming_table_NO_missing_d[,.N, by = .(v_gene, d_gene)]

    # for each row with missing d gene...
    for (row in 1:nrow(patient_trimming_table_missing_d)){
        # Find observed v and j gene
        v_gene_observed = patient_trimming_table_missing_d[row]$v_gene
        j_gene_observed = patient_trimming_table_missing_d[row]$j_gene

        # Find observed triplets containing the observed v,j (or both) gene
        empirical_dist_vj = triplet_counts[v_gene == v_gene_observed & j_gene == j_gene_observed]
        empirical_dist_j = doublet_dj[j_gene == j_gene_observed]
        empirical_dist_v = doublet_vd[v_gene == v_gene_observed]

        # If both the v and j gene are observed empirically, sample from the empirical distribution to infer missing d gene
        if (nrow(empirical_dist_vj) != 0) {
            empirical_dist_vj$probability = empirical_dist_vj$N/sum(empirical_dist_vj$N)
            sample = rmultinom(1,1,empirical_dist_vj$probability)
            index = which(sample == 1)
            patient_trimming_table_missing_d[row]$d_gene = empirical_dist_vj[index]$d_gene
        # If both the v and j gene are NOT observed empirically, sample from the empirical distribution of just the j gene observations to infer missing d gene
        } else if (nrow(empirical_dist_j) != 0) {
            empirical_dist_j$probability = empirical_dist_j$N/sum(empirical_dist_j$N)
            sample = rmultinom(1,1,empirical_dist_j$probability)
            index = which(sample == 1)
            patient_trimming_table_missing_d[row]$d_gene = empirical_dist_j[index]$d_gene
        } else {
            empirical_dist_v$probability = empirical_dist_v$N/sum(empirical_dist_v$N)
            sample = rmultinom(1,1,empirical_dist_v$probability)
            index = which(sample == 1)
            patient_trimming_table_missing_d[row]$d_gene = empirical_dist_v[index]$d_gene
        } 
        # Now, trim or insert half of the dgene (length is different depending on dgene)
        if (patient_trimming_table_missing_d[row]$d_gene == "TRBD2*01" | patient_trimming_table_missing_d[row]$d_gene == "TRBD2*02"){
            patient_trimming_table_missing_d[row]$d0_trim = 4
            patient_trimming_table_missing_d[row]$d1_trim = 4
            patient_trimming_table_missing_d[row]$dj_insert = 4
            patient_trimming_table_missing_d[row]$vd_insert = 4
        } else if (patient_trimming_table_missing_d[row]$d_gene == "TRBD1*01"){
            patient_trimming_table_missing_d[row]$d0_trim = 3
            patient_trimming_table_missing_d[row]$d1_trim = 3
            patient_trimming_table_missing_d[row]$dj_insert = 3
            patient_trimming_table_missing_d[row]$vd_insert = 3
        }
    }
    inferred = rbind(patient_trimming_table_missing_d, patient_trimming_table_NO_missing_d)
    return(inferred)
}

combine_genes_by_common_cdr3 <- function(gene_type){
    cdr3 = fread(file = paste0(project_path, '/tcr-gwas/_ignore/human_vj_allele_cdr3_nucseqs.tsv'))
    cdr3$gene = str_split(cdr3$id, fixed('*'), simplify = TRUE)[,1]
    # add d genes
    cdr3 = rbind(cdr3, data.frame(cdr3_nucseq = c('NA', 'NA', 'NA'), id = c("TRBD1*01", "TRBD2*01", "TRBD2*02"), gene = c("TRBD1*01", "TRBD2*01", "TRBD2*02")))
    cdr3 = cdr3[,cdr3_gene_group := .GRP, by = .(cdr3_nucseq, gene)]
    #cdr3$gene_type = paste0(tolower(str_sub(cdr3$id, 4,4)), '_gene')
    #cdr3 = cdr3[gene_type == gene_type]
    cdr3_small = cdr3[,c('id', 'cdr3_gene_group')]
    return(cdr3_small)
}


condense_trim_data <- function(trim_types, infer_d, condensing){
    dir = paste0(project_path, '/tcr-gwas/_ignore/emerson_stats/')
    files = list.files(path=dir, pattern="*.tsv", full.names=TRUE)
    
    for (trim in trim_types){
        assign(paste0('condensed_', trim, '_data'), data.table())
    }
    all_patients = data.table()
    count = 0
    
    # for each patient...
    for (file in files){
        count = count + 1
        # read file...
        temp_file_original = data.table()
        temp_file_condensed = data.table()
        temp_file_original = fread(file, sep = "\t", fill=TRUE, header = TRUE)
        temp_file_original = as.data.table(temp_file_original)
        file_name = str_split(file, "/")[[1]][8]
        file_root_name = str_split(file_name, ".tsv")[[1]][1]
        patient_id = str_split(file_root_name, "_")[[1]][3]
        temp_file_original$patient_id = patient_id

        # Infer d genes (when d gene is missing)
        if (infer_d == 'True'){
            temp_file_original = infer_d_gene(temp_file_original)
        } else {
            temp_file_original = temp_file_original[d_gene != '-']
        }
        if (condensing == 'by_gene' | condensing == 'by_cdr3'){        
            for (trim in trim_types){
                together = data.table()

                if ((str_split(trim, "_")[[1]][2])== "trim"){
                    # extract gene type (separate from '_trim')
                    gene_type = paste0(substr(trim, 1, 1), '_gene')
                    if (condensing == 'by_cdr3'){
                        cdr3 = combine_genes_by_common_cdr3(gene_type)
                        temp_file = merge(temp_file_original, cdr3, by.x = gene_type, by.y = 'id')
                        conditioning_variable = 'cdr3_gene_group' 
                    } else if (condensing == 'by_gene'){
                        conditioning_variable = gene_type
                    }

                    # condense patient file by gene type (and take the mean trimming for each gene type)
                    together = temp_file[,.(mean(v_trim), mean(d0_trim), mean(d1_trim), mean(j_trim), mean(vj_insert), mean(dj_insert), mean(vd_insert), .N), by = .(patient_id, eval(parse(text=conditioning_variable)), productive)]
                    colnames(together) = c("patient_id", paste(gene_type), "productive", "v_trim", "d0_trim", "d1_trim", "j_trim", "vj_insert","dj_insert", "vd_insert", paste0(gene_type, '_count'))
                } else if ((str_split(trim, "_")[[1]][2])== "insert"){
                    # extract gene type (separate from '_trim')
                    gene_type1 = paste0(substr(trim, 1, 1), '_gene')
                    gene_type2 = paste0(substr(trim, 2, 2), '_gene')
                    gene_type = paste0(substr(trim, 1, 2), '_gene')

                    # condense patient file by gene type (and take the mean trimming for each gene type)
                    together = temp_file[,.(mean(v_trim), mean(d0_trim), mean(d1_trim), mean(j_trim), mean(vj_insert), mean(dj_insert), mean(vd_insert), .N), by = .(patient_id, eval(parse(text=gene_type1)), eval(parse(text=gene_type2)), productive)]

                    colnames(together) = c("patient_id", paste(gene_type1), paste(gene_type2), "productive", "v_trim", "d0_trim", "d1_trim", "j_trim", "vj_insert","dj_insert", "vd_insert", paste0(gene_type, '_count'))
                }
            
                together[[paste0('weighted_', gene_type, '_count')]] = together[[paste0(gene_type, '_count')]]/nrow(temp_file)
                assign(paste0('condensed_', trim, '_data'), rbind(get(paste0('condensed_', trim, '_data')), together))
                print(paste0(count," of ", length(files), " processed for ", trim))
            
                filename = ifelse(infer_d == 'True', paste0(project_path, '/tcr-gwas/_ignore/', condensing, '_condensed_', trim,'_data_all_patients.tsv'), paste0(project_path, '/tcr-gwas/_ignore/', condensing, '_condensed_', trim,'_data_all_patients_NO_d_gene_infer.tsv'))
        
                write.table(get(paste0('condensed_', trim, '_data')), file=filename, quote=FALSE, sep='\t', col.names = NA)
            }
        } else if (condensing == 'by_patient'){
            filename = ifelse(infer_d == 'True', paste0(project_path, '/tcr-gwas/_ignore/by_patient_condensed_data_all_patients.tsv'), paste0(project_path, '/tcr-gwas/_ignore/by_patient_condensed_data_all_patients_NO_d_infer.tsv'))
	    
            together = temp_file[,.(mean(v_trim), mean(d0_trim), mean(d1_trim), mean(j_trim), mean(vj_insert), mean(dj_insert), mean(vd_insert), .N), by = .(patient_id, productive)]

	        colnames(together) = c("patient_id", "productive", "v_trim" , "d0_trim", "d1_trim", "j_trim", "vj_insert","dj_insert", "vd_insert", "tcr_count")        
	        all_patients = rbind(all_patients, together)

            write.table(all_patients, file=filename, quote=FALSE, sep='\t', col.names = NA) 

            print(paste0(count," of ", length(files), " processed"))
        }
    }
    print(paste0("finished processing!"))
}


#stopifnot(args[1] %in% c('trimming', 'insertion'))
#if (args[1] == 'trimming'){
#    types = c('v_trim', 'd0_trim', 'd1_trim', 'j_trim')
#} else {
#    types = c('vd_insert', 'vj_insert', 'dj_insert')
#}

#stopifnot(args[2] %in% c('True', 'False'))
#stopifnot(args[3] %in% c('by_cdr3', 'by_gene', 'by_patient'))

#condense_trim_data(trim_types = types, infer_d = args[2], condensing = args[3])

