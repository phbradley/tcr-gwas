library("data.table")
library("COMBAT")
library("plyr")
library("readr")
library("stringr")
library("gdsfmt")
library("SNPRelate")
library("GWASTools")


# Need to first define a cutoff for correlating snps with one another. 
# Let's define a proximity cutoff and a p_value difference cutoff

bonferroni = 0.05/35481497
proximity_cutoff = 1000
pvalue_cutoff = 3.594536e-05 #0.05/nrow(together)

# This method finds a bonferonni significant snp, then searches forward and backward in space until the snps are either 1. too far away (with proximity_cutoff) or 2. not significant enough (with pvalue_cutoff)

correlate_snps <- function(snp_id_list, snp_meta_data, productivity, trim_type, bonferroni, proximity_cutoff, pvalue_cutoff){
    if (nrow(snp_meta_data) == 11267){
        snp_meta_data = snp_meta_data[,mean(mwu_pval), by = .(snpnum, chromosome, hg19_pos)]
        snp_meta_data$snpnum = paste0("snp", snp_meta_data$snpnum)
        colnames(snp_meta_data) = c("snp", "chr", "hg19_pos", "phile_pval")
    } else if (colnames(snp_meta_data)[1] == "chromosome"){
        snp_meta_data$snpid = paste0("snp", snp_meta_data$snpid)
        colnames(snp_meta_data) = c("chr", "feature", "hg19_pos", "phil_pval", "snp")
    } else {
        snp_meta_data = snp_meta_data[,1:3]
        snp_meta_data$snpid = paste0("snp", snp_meta_data$snpid)
        colnames(snp_meta_data) = c("snp", "hg19_pos", "chr")
    }

    regression_snp_list = as.data.table(read.table(paste0("regression_bootstrap_results/", productivity, "/", trim_type, "/", trim_type, "_", productivity, "_snplist_", snp_id_list[1], "-",snp_id_list[length(snp_id_list)], "_snps_varying_intercepts_by_subject.tsv"), header = TRUE))
    together = merge(regression_snp_list, snp_meta_data, all.x = TRUE, by = "snp")

    correlated_snps = data.table()
    j = 1
    together = together[order(hg19_pos)]

    while (j < nrow(together)){
        if (together[j]$pvalue < bonferroni){
            sig_snp_temp = data.table()
            index = j
            for (i in seq(j-1, 1)){
                if (index == 1){ break }
                else if ((abs(together[i]$hg19_pos - together[index]$hg19_pos )<= proximity_cutoff) && (together[i]$pvalue < pvalue_cutoff) && (together[i]$chr == together[index]$chr)){
                    sig_snp_temp = rbind(sig_snp_temp, together[i])
                    index = index - 1 
                } else { break }  
            }

            sig_snp_temp = rbind(sig_snp_temp, together[j])

            index = j
            for(k in seq(j+1, nrow(together))){
                if ((abs(together[k]$hg19_pos - together[index]$hg19_pos)<= proximity_cutoff) && (together[k]$pvalue < pvalue_cutoff)){
                    sig_snp_temp = rbind(sig_snp_temp, together[k]) 
                    index = index + 1 
                } else { break }
            }

            correlated_snps = rbind(correlated_snps, data.table(index = j, sig_snp_temp))  
            j = j + nrow(sig_snp_temp)
        } else {
            j = j + 1
        }
    }
    print("finished correlating snps")
    if (nrow(correlated_snps) >= 1){
        if (ncol(correlated_snps) == 10){
            colnames(correlated_snps) = c("index","snp", "intercept", "slope", "standard_error", "p", "chr", "feature", "start", "phil_pval")
        }
    }
    return(correlated_snps)
}


group_p_values <- function(snp_id_list, snp_meta_data, productivity, trim_type, bonferroni, proximity_cutoff, pvalue_cutoff){
    correlated_snps = correlate_snps(snp_id_list, snp_meta_data, productivity, trim_type, bonferroni, proximity_cutoff, pvalue_cutoff)
    if (nrow(correlated_snps) >= 1){
        correlated_snps_with_group_p_val = data.table()
        snps_gds = snpgdsOpen("../_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds")
        n <- index.gdsn(snps_gds, "sample.id")
        sampleid <- read.gdsn(n)
        for (i in unique(correlated_snps$index)){
            if(length(correlated_snps$feature) > 0){
                snps = unique(correlated_snps[index == i][,-c(8,10)])
            } else {
                snps = unique(correlated_snps[index == i])
            }
            if (nrow(snps)>1){
                genotype_list = data.table(scanID = c(sampleid))
                colnames(genotype_list) = c("scanID")
                genotype_list$scanID = as.numeric(as.character(genotype_list$scanID))
                subject_id_mapping = as.data.table(read.table('../_ignore/snp_data/gwas_id_mapping.tsv', sep = "\t", fill=TRUE, header = TRUE, check.names = FALSE))
                genotypes_subjects = merge(genotype_list, subject_id_mapping, by = "scanID")
                for (snp in snps$snp){
                    snp_id = as.numeric(substring(snp, 4))
                    genotype = snpgdsGetGeno(snps_gds, snp.id=snp_id)
                    genotypes_df = data.table(genotype)
                    colnames(genotypes_df) = c(snp)
                    # Convert subject names and compile condensed data: 
                    genotypes_subjects = cbind(genotypes_subjects, genotypes_df)
                }
                #This may throw an error if there are NA entries: 
                snps$group_p_val = COMBAT(x = snps$p, snp.ref = genotypes_subjects[,-c(1,2)])[1]
            } else {
                snps$group_p_val = snps$p
            }
            correlated_snps_with_group_p_val= rbind(correlated_snps_with_group_p_val, snps)
        }
        closefn.gds(snps_gds)
        print("finished p-value grouping")

        name = paste0('regression_bootstrap_results/', productivity, '/', trim_type, '/', trim_type, '_correlated_snps_from_snplist_', snp_id_list[1], "-",snp_id_list[length   (snp_id_list)], '_snps_varying_intercepts_by_subject.tsv')

        write.table(correlated_snps_with_group_p_val, file= name, quote=FALSE, sep='\t', col.names = NA)
    } else {
        print("no correlated snps!!!")
    } 
}

#condense_correlated_snps_pvals <- function(correlated_snps){

#}