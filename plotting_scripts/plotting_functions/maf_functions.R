find_significant_snps <- function(phenotype_list, signficance_cutoff){
    sig_data = compile_manhattan_plot_data(phenotype_list) 
    return(sig_data[pvalue < signficance_cutoff])
}

determine_true_minor_allele <- function(snp, phenotype_genotype_dt){
    columns = c('feature_of_interest', paste(snp))
    simplified_dt = phenotype_genotype_dt[productivity == 'productive'][,..columns]
    regression = lm(simplified_dt$feature_of_interest ~ simplified_dt[[paste(snp)]])
    slope = coef(regression)['simplified_dt[[paste(snp)]]']
    minor_allele = ifelse(slope < 0, 2, 0)
    return(minor_allele)
}

calculate_maf <- function(snp, minor_allele, race, genotype_dt){
    columns = c('localID', 'race.g', paste(snp))
    if (race != 'all'){
        data = genotype_dt[race.g == race,..columns]
    } else {
        data = genotype_dt[,..columns]
    }

    if (minor_allele != 2){
        major_hom = data[get(snp) == 2][, (snp) := 0]
        minor_het = data[get(snp) == 1]
        minor_hom = data[get(snp) == 0][, (snp) := 2]
        data = rbind(minor_hom, minor_het, major_hom)
    }
    total_possible_alleles = 2*sum(!is.na(data[[snp]]))
    observed_alleles = sum(data[[snp]], na.rm = TRUE)
    maf = observed_alleles/total_possible_alleles
    return(maf)
}

find_snp_start_by_position <- function(chromosome, position1, position2){
    snp_meta_data = fread(SNP_META_DATA_FILE)
    colnames(snp_meta_data) =c('snpid','snppos', 'snpchrome', 'snpallele', 'rsid')
    snp_meta_data$snpindex = seq(1, nrow(snp_meta_data))
    filtered = snp_meta_data[snpchrome == chromosome & snppos > position1 & snppos < position2]
    filtered_ordered = filtered[order(filtered$snpindex),]
    return(c(filtered_ordered$snpindex[1],
            filtered_ordered$snpindex[nrow(filtered_ordered)]-filtered_ordered$snpindex[1]))
}

map_scanID_to_localID <- function(scanIDs_to_convert){
    ID_map_file = fread(ID_MAPPING_FILE)
    converted_IDs = plyr::mapvalues(scanIDs_to_convert, 
                                    ID_map_file$scanID, 
                                    ID_map_file$localID)
    return(converted_IDs)
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


