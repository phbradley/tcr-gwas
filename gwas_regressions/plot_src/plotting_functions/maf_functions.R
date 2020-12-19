find_significant_snps <- function(phenotype_list, signficance_cutoff){
    sig_data = compile_manhattan_plot_data(phenotype_list) 
    return(sig_data[pvalue < signficance_cutoff])
}

determine_true_minor_allele <- function(snp, phenotype_genotype_dt){
    columns = c('total_inserts', paste(snp))
    simplified_dt = phenotype_genotype_dt[productive == 'TRUE'][,..columns]
    regression = lm(simplified_dt$total_inserts ~ simplified_dt[[paste(snp)]])
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
    snp_meta_data = fread(paste0(PROJECT_PATH, '/tcr-gwas/_ignore/snps_meta_data.tsv'))
    colnames(snp_meta_data) =c('snpindex', 'snpid', 'snppos', 'snpchrome', 'snpallele')
    filtered = snp_meta_data[snpchrome == chromosome & snppos > position1 & snppos < position2]
    filtered_ordered = filtered[order(filtered$snpindex),]
    return(c(filtered_ordered$snpindex[1],
             filtered_ordered$snpindex[nrow(filtered_ordered)]-filtered_ordered$snpindex[1]))
}
