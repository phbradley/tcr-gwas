library("tidyverse")
library("gdsfmt")
library("SNPRelate")
library("GWASTools")
library(data.table)
setDTthreads(threads = 1)

source(paste0(project_path, '/tcr-gwas/trimming_regression/scripts/condense_trimming_data.R'))

compile_condensed_trimming_data <- function(trim_type, d_infer, condensing){
    print('compiling TCR trimming/insertion data')
    # import condensed trimming file
    if (condensing == 'by_patient'){
        if (d_infer == 'False'){
            filename = paste0(project_path, "/tcr-gwas/_ignore/by_patient_condensed_data_all_patients_NO_d_infer.tsv")
        } else {
            filename = paste0(project_path, "/tcr-gwas/_ignore/by_patient_condensed_data_all_patients.tsv")
        }
    } else {
        if (d_infer == 'False'){
            filename = paste0(project_path, "/tcr-gwas/_ignore/condensed_", trim_type, "_data_all_patients_NO_d_gene_infer.tsv")
        } else {
            filename = paste0(project_path, "/tcr-gwas/_ignore/condensed_", trim_type, "_data_all_patients.tsv")
        }
    }

    if (!file.exists(filename)){
        print('Data condensing required. Computing now.')
        condense_trim_data(trim_types = trim_type, infer_d = d_infer, condensing)
    } 

    assign('trimming_data', as.data.frame(read.table(filename, sep = "\t", fill=TRUE, header = TRUE)[-1]))
    names(trimming_data)[names(trimming_data) == "patient_id"] <- "localID"

    return(trimming_data)
}

remove_small_repertoire_observations <- function(condensed_trimming_data, productive_log10_count_cutoff, NOT_productive_log10_count_cutoff, trim_type, d_infer){
    tcr_count_dt = as.data.table(compile_condensed_trimming_data(trim_type, d_infer, condensing = 'by_patient'))[,c('localID', 'productive', 'tcr_count')]
    condensed_trim_dt = as.data.table(condensed_trimming_data)    
    condensed_trim_dt = merge(condensed_trim_dt, tcr_count_dt)
    condensed_trim_dt_filtered = condensed_trim_dt[(productive == 'TRUE' & log10(tcr_count)>productive_log10_count_cutoff) | (productive == 'FALSE' & log10(tcr_count)>NOT_productive_log10_count_cutoff)]
    return(condensed_trim_dt_filtered)
}

# This function makes a snp file from chromosme and position data using the gds file
snp_file <- function(chromosome, position1, position2){
    snps_gds = snpgdsOpen(paste0(project_path, "/tcr-gwas/_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds"))
    snp_chrom <- read.gdsn(index.gdsn(snps_gds, "snp.chromosome"))

    chr_index_start = match(chromosome,snp_chrom)
    chr_index_end = match(chromosome+1,snp_chrom)-1

    snp_pos <- read.gdsn(index.gdsn(snps_gds, "snp.position"), start = chr_index_start, count=(chr_index_end-chr_index_start)+1)
    snp_id <- read.gdsn(index.gdsn(snps_gds, "snp.id"), start = chr_index_start, count=(chr_index_end-chr_index_start)+1)

    snps_chromosome = data.frame(snp = snp_id, chr = snp_chrom[chr_index_start:chr_index_end], hg19_pos = snp_pos)

    if (position1 != 'all'){
        snps_chromosome = snps_chromosome[hg19_pos < position2 & hg19_pos > position1]
    }
    closefn.gds(snps_gds)
    return(snps_chromosome)
}

# This function makes a snp file from snp index start and snp count using the gds file
snp_file_by_snp_start <- function(snp_start, count){
    snps_gds = openfn.gds(paste0(project_path, "/tcr-gwas/_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds"))
    snp_id <- read.gdsn(index.gdsn(snps_gds, "snp.id"))

    if ((snp_start+count) > length(snp_id)){
        count = length(snp_id) - snp_start
    }

    snp_chrom <- read.gdsn(index.gdsn(snps_gds, "snp.chromosome"), start = snp_start, count=count)
    snp_pos <- read.gdsn(index.gdsn(snps_gds, "snp.position"), start = snp_start, count=count)
    
    snps = data.frame(snp = snp_id[snp_start:(snp_start+count-1)], chr = snp_chrom, hg19_pos = snp_pos)

    closefn.gds(snps_gds)
    return(snps)
}

# This function makes a genotype file from snp index start and snp count using the gds file
compile_all_genotypes <- function(snp_start, count) {
    gfile = openfn.gds(paste0(project_path, "/tcr-gwas/_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds"))
    bigsize <- 35481497
    #start <- (i-1)*sz + 1
    numrows <- min( count, bigsize-snp_start+1 )
    
    genotype_matrix <- read.gdsn(index.gdsn(gfile, "genotype"), start=c(1,snp_start), count=c(398, numrows))
    sample_ids <- read.gdsn(index.gdsn(gfile, "sample.id"), start=1, count=398)
    snp_ids <- read.gdsn(index.gdsn(gfile, "snp.id"), start=snp_start, count=numrows)
    closefn.gds(gfile)

    subject_id_mapping = as.data.frame(read.table(paste0(project_path, '/tcr-gwas/_ignore/snp_data/gwas_id_mapping.tsv'), sep = "\t", fill=TRUE, header = TRUE, check.names = FALSE))
    sample_ids_converted = plyr::mapvalues(sample_ids, subject_id_mapping$scanID, subject_id_mapping$localID)
    rownames(genotype_matrix) = sample_ids_converted
    colnames(genotype_matrix) = snp_ids
    return(genotype_matrix)
}

# This function compiles all subjects
compile_subjects <- function() {
    gfile = openfn.gds(paste0(project_path, "/tcr-gwas/_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds"))
    subjects <- read.gdsn(index.gdsn(gfile, "sample.id"))
    closefn.gds(gfile)
    return(subjects)
}

# This function finds the snp start index given chromosome and positions
find_snp_start_by_position <- function(chromosome, position1, position2){
    library(data.table)
    snp_meta_data = fread(paste0(project_path, '/tcr-gwas/_ignore/snps_meta_data.tsv'), sep = "\t", fill=TRUE, header = TRUE, check.names = FALSE)
    colnames(snp_meta_data) =c('snpindex', 'snpid', 'snppos', 'snpchrome', 'snpallele')
    filtered = snp_meta_data[snpchrome == chromosome & snppos > position1 & snppos < position2]
    filtered_ordered = filtered[order(filtered$snpindex),]
    return(c(filtered_ordered$snpindex[1], filtered_ordered$snpindex[nrow(filtered_ordered)]-filtered_ordered$snpindex[1]))
}

# This function takes a snps dataframe and makes a coninciding snp file (with chromosome and position columns)
make_snp_file_subset_by_count_and_index <- function(snp_dataframe, count, index){
    total = nrow(snp_dataframe)
    converted_index = index - 1
    start = converted_index * count + 1
    if (start > nrow(snp_dataframe)){
        return(data.frame())
    } else {
        end = min(start-1 + count, total)
        return(unique(snp_dataframe[start:end,]))
    }  
}


# This function filters trimming data based on productivity status
filter_by_productivity <- function(condensed_trimming_dataframe, productive){
    # subset trimming data to include only productive or not productive entires
    if (productive == "True"){
        condensed_trimming_dataframe = condensed_trimming_dataframe %>% filter(productive == "TRUE")
    } else if (productive == "False"){
        condensed_trimming_dataframe = condensed_trimming_dataframe %>% filter(productive == "FALSE")
    } 
    return(condensed_trimming_dataframe)
}

# remove snp_genotype columns that are either all NA, or only have one genotype (for everyone)
remove_matrix_column_by_genotype <- function(genotype_matrix){
    for (snp in colnames(genotype_matrix)){
        genotypes = unique(genotype_matrix[,snp])
        nonNA_genotypes = genotypes[genotypes != 3]
        if (length(nonNA_genotypes) <= 1){
            genotype_matrix = genotype_matrix[, colnames(genotype_matrix) != snp]
        }
    }
    # replace 3 with NA
    genotype_matrix[genotype_matrix == 3] <- NA
    return(genotype_matrix)
}

read_genotype_pca <- function(type){
    stopifnot(type %in% c('none', 'pc_air', 'pc_fixed_maggie', 'pc_fixed_dave'))
    subject_id_conversion = read.table(paste0(project_path, '/tcr-gwas/_ignore/snp_data/gwas_id_mapping.tsv'), sep = "\t", fill=TRUE, header = TRUE)
    if (type == 'pc_fixed_maggie'){
        pca_file_name = paste0(project_path, '/tcr-gwas/_ignore/snp_data/population_structure_pca_by_LD_snps.tsv')
        if (!file.exists(pca_file_name)){
            system(command = paste0("Rscript ", project_path, "/tcr-gwas/trimming_regression/scripts/population_structure_pca.R "))
        }
    } else if (type == 'pc_fixed_dave'){
        pca_file_name = paste0(project_path, '/tcr-gwas/_ignore/snp_data/pc_fixed_08Nov2020.txt')
    } else if (type == 'pc_air'){
        pca_file_name = paste0(project_path, '/tcr-gwas/_ignore/snp_data/pc_pcair_08Nov2020.txt')
    } 
    if (type != 'none'){
        pca = read.table(pca_file_name, sep = '\t', fill = TRUE, header = TRUE)
        if ('sample_id' %in% colnames(pca)){
            pca = merge(subject_id_conversion, pca, by.x = 'scanID', by.y = 'sample_id')[,-c(1,3)]
        }
    } else {
        pca = NULL
    }
    return(pca)
}
    
make_regression_file_name <- function(snp_list, trim_type, condensing, random_effects, d_infer, repetitions, pca_structure_correction){
    random_effects_name = ifelse(random_effects == 'True', 'random_effects', 'no_random_effects')
    d_infer_name = ifelse(d_infer == 'True', '_d_infer', '_no_d_infer')
    pca_name = ifelse(pca_structure_correction == 'True', '_pca_correction', '_no_pca_correction')
   boots = paste0('_', repetitions, '_bootstraps')

    file_name = paste0(output_path, '/cluster_job_results/', trim_type, '/', random_effects_name, d_infer_name, pca_name, boots, '/', condensing, '/', trim_type, '_',snp_list$snp[1], '-', snp_list$snp[nrow(snp_list)],'_snps_regression_with_weighting_condensing_', condensing, '.tsv')
    if (!dir.exists(paste0(output_path, '/cluster_job_results'))){
        system(paste0('mkdir ', output_path, '/cluster_job_results'))
    }
 
    if (!dir.exists(paste0(output_path, '/cluster_job_results/', trim_type))){
        system(paste0('mkdir ', output_path, '/cluster_job_results/', trim_type))
    }
    
    if (!dir.exists(paste0(output_path, '/cluster_job_results/', trim_type, '/', random_effects_name, d_infer_name, pca_name, boots))){
        system(paste0('mkdir ', output_path, '/cluster_job_results/', trim_type, '/', random_effects_name, d_infer_name, pca_name, boots))
    }

    if (!dir.exists(paste0(output_path, '/cluster_job_results/', trim_type, '/', random_effects_name, d_infer_name, pca_name, boots, '/', condensing))){
        system(paste0('mkdir ', output_path, '/cluster_job_results/', trim_type, '/', random_effects_name, d_infer_name, pca_name, boots, '/', condensing))
    }

    return(file_name)
}

make_regression_file_path <- function(trim_type, condensing, random_effects, d_infer, repetitions, pca_structure_correction, output_path){
    random_effects_name = ifelse(random_effects == 'True', 'random_effects', 'no_random_effects')
    d_infer_name = ifelse(d_infer == 'True', '_d_infer', '_no_d_infer')
    pca_name = ifelse(pca_structure_correction == 'True', '_pca_correction', '_no_pca_correction')
    boots = paste0('_', repetitions, '_bootstraps')

    path_name = paste0(output_path, '/cluster_job_results/', trim_type, '/', random_effects_name, d_infer_name, pca_name, boots, '/', condensing)
    if (!dir.exists(paste0(output_path, '/cluster_job_results'))){
        system(paste0('mkdir ', output_path, '/cluster_job_results'))
    }
 
    if (!dir.exists(paste0(output_path, '/cluster_job_results/', trim_type))){
        system(paste0('mkdir ', output_path, '/cluster_job_results/', trim_type))
    }
    
    if (!dir.exists(paste0(output_path, '/cluster_job_results/', trim_type, '/', random_effects_name, d_infer_name, pca_name, boots))){
        system(paste0('mkdir ', output_path, '/cluster_job_results/', trim_type, '/', random_effects_name, d_infer_name, pca_name, boots))
    }

    if (!dir.exists(paste0(output_path, '/cluster_job_results/', trim_type, '/', random_effects_name, d_infer_name, pca_name, boots, '/', condensing))){
        system(paste0('mkdir ', output_path, '/cluster_job_results/', trim_type, '/', random_effects_name, d_infer_name, pca_name, boots, '/', condensing))
    }

    cat(path_name)
}


# This script compiles all minor allele fraction data (to be used in plotting cutoff)
compile_all_maf_data <- function(){
    data_files = list.files(path=paste0(output_path, '/maf_results/'), pattern='*', full.names=TRUE)   

    print(paste0('compile file list'))
    assign(paste0('together'), NULL)
    count = 0
    for (file in data_files){
        # read file...
        if (file.size(file) == 1 | file.size(file) == 0){
            next
        }
        temp_file = fread(file, sep = "\t", fill=TRUE, header = TRUE)
        #together = rbind(together, temp_file)
        if (ncol(temp_file) > 2){
            assign(paste0('together'), rbindlist(list(get(paste0('together')), temp_file)))
        }
        count = count + 1
        print(paste0(count, ' of ', length(data_files), ' completed'))
    }
    
    file_name = paste0('maf_all_snps.tsv')

    write.table(together, file=paste0(output_path, '/', file_name), quote=FALSE, sep='\t', col.names = NA)
}


