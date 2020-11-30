# Miscellaneous functions!

#library("SNPRelate")
# This function makes a genotype file given a random snp file
make_genotype_file_given_random_snps <- function(snp_file){
    snps_gds = snpgdsOpen(paste0(project_path, "/tcr-gwas/_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds"))
    snp_id_list = as.numeric(gsub('snp', '', snp_file$snp))
    subject_id_mapping = as.data.frame(read.table(paste0(project_path, '/tcr-gwas/_ignore/snp_data/gwas_id_mapping.tsv'), sep = "\t", fill=TRUE, header = TRUE, check.names = FALSE))
    genotype_list = data.frame(localID = subject_id_mapping$localID)
    for (snp in snp_id_list){
        genotype = snpgdsGetGeno(snps_gds, snp.id=snp, with.id = TRUE)
        genotype_matrix = data.frame(genotype$genotype, scanID = as.numeric(as.character(genotype$sample.id)))
        subject_id_HIP_genotypes = merge( genotype_matrix, subject_id_mapping, by = "scanID")
        colnames(subject_id_HIP_genotypes) = c('scanID', paste(snp), 'localID')
        genotype_list = merge(genotype_list, subject_id_HIP_genotypes[,c('localID', paste(snp))], by = 'localID')
    }
    rownames(genotype_list) = genotype_list$localID
    closefn.gds(snps_gds)
    return(genotype_list)
}


# This function compiles trimming data crosses
compile_trimming_data_cross <- function(){
    trimming_data_by_gene_all = data.frame()
    for (trim in c('v_trim', 'd0_trim', 'j_trim')){
        assign(paste0(trim, '_trimming_data'), as.data.frame(read.table(paste0(project_path, "/tcr-gwas/_ignore/condensed_", trim, "_data_all_patients.tsv"), sep = "\t", fill=TRUE, header = TRUE)[-1]))
        setnames(get(paste0(trim, '_trimming_data')), "patient_id", "localID")
        setnames(get(paste0(trim, '_trimming_data')), paste0(substring(trim, 1, 1), '_gene'), "gene")
        setnames(get(paste0(trim, '_trimming_data')), paste0(substring(trim, 1, 1), '_gene_count'), "gene_count")
        setnames(get(paste0(trim, '_trimming_data')), paste0('weighted_', substring(trim, 1, 1), '_gene_count'), "weighted_gene_count")
        gene_type = data.frame(gene_class = rep(paste0(substr(trim, 1, 1), '_gene'), nrow(get(paste0(trim, '_trimming_data')))))
        assign(paste0(trim, '_trimming_data'), cbind(get(paste0(trim, '_trimming_data')), gene_type))
        trimming_data_by_gene_all = rbind(trimming_data_by_gene_all, get(paste0(trim, '_trimming_data')))
    }
    return(trimming_data_by_gene_all)
}


# This function finds a regression file (from cluster) and subsets based on signficance cutoff (for use in bootstraps)
open_regressed_file_and_subset_by_pval <- function(significance_cutoff, trim_type, random_effects, condensing, d_infer, maf_cutoff){
    if (random_effects == 'True'){
         file_name = paste0('_', trim_type, '_snps_regression_with_weighting_condensing_', condensing, '_with_random_effects_0_bootstraps')
    } else {
         file_name = paste0('_',trim_type, '_snps_regression_with_weighting_condensing_', condensing, '_NO_random_effects_0_bootstraps')
    }

    if (d_infer == 'False'){
        file_name = paste0(file_name, '_NO_d_infer.tsv')
    } else {
        file_name = paste0(file_name, '.tsv')
    }

    productive_data = fread(paste0(output_path, '/results/productive', file_name), sep = "\t", fill=TRUE, header = TRUE)
    not_productive_data = fread(paste0(output_path, '/results/NOT_productive', file_name), sep = "\t", fill=TRUE, header = TRUE)
    maf_data = fread(paste0(output_path, '/maf_all_snps.tsv'), sep = "\t", fill=TRUE, header = TRUE)[,-c(1)]

    data = rbind(productive_data, not_productive_data)[,-c(1,2)]
    data = merge(data, maf_data, by = 'snp')
    data = data[maf > maf_cutoff]
    data = data[pvalue < significance_cutoff]
    return(data[,c(1:3)])
}


# This script makes a manhattan plot for a specific region
compile_data_manhattan <- function(snp_meta_data, snp_id_list, correlated_snps, paired_productive_snps, productivity, varying_int, trim_type, chromosome, correlate_snps, gene_conditioning, weighting, gene_region){
    gene_region = gene_region
    bonferroni = 0.05/35481497

    # sort the snp meta data from phil (i.e. assign common column names)
    if (nrow(snp_meta_data) == 11267){
        snp_meta_data = snp_meta_data[,mean(mwu_pval), by = .(snpnum, chromosome, hg19_pos)]
        snp_meta_data$snpnum = paste0("snp", snp_meta_data$snpnum)
        colnames(snp_meta_data) = c("snp", "chr", "hg19_pos", "phil_pval")
    } else if (nrow(snp_meta_data) == 96){
        snp_meta_data$snpid = paste0("snp", snp_meta_data$snpid)
        colnames(snp_meta_data) = c("chr", "feature", "hg19_pos", "phil_pval", "snp")
    } else {
        snp_meta_data = snp_meta_data[,1:3]
        snp_meta_data$snp = paste0("snp", snp_meta_data$snp)
        colnames(snp_meta_data) = c("snp", "chr", "hg19_pos")
    }

    if (productivity == 'productive'){
        productivity_status = "True"
    } else if (productivity == 'NOT_productive'){
        productivity_status = "False"
    }

    # load regression / pvalue results
    if (varying_int == "True"){
        assign('regression_snp_list', as.data.table(read.table(generate_file_name(snp_id_list, trim_type, gene_type = paste0(substr(trim_type, 1, 1), '_gene'), productivity = productivity_status, gene_conditioning = gene_conditioning, weighting = weighting, condensing = 'by_gene', repetitions = 100), header = TRUE)))
        gene_cond = ifelse(gene_conditioning == 'True', '_gene_conditioning', '')
        weight_cond = ifelse(weighting == 'True', '_weight', '')
        regression_type = paste0('lmer', gene_cond, weight_cond)
    } else {
        regression_snp_list = as.data.table(read.table(paste0("regression_bootstrap_results/", productivity, "/", trim_type, "/", trim_type, "_", productivity, "_snplist_", snp_id_list[1], "-",snp_id_list[length(snp_id_list)], "_snps_simple_condensing_by_patient.tsv"), header = TRUE))
        regression_type = 'Simple regression: no conditioning on gene'
    }

    if (correlate_snps == 'True'){
        regression_snp_list_correlated = correlated_snps
        regression_snp_list_correlated$snp = paste0('snp', regression_snp_list_correlated$snp)
        regression_snp_list_correlated = merge(regression_snp_list_correlated, regression_snp_list, all.x = TRUE)
        regression_snp_list_no_correlated = subset(regression_snp_list, !(snp %in% regression_snp_list_correlated$snp))
        correlated_snps = paste0("Correlated snps shown by color")
        # define dataframe containing correlated snps and meta data
        together_correlated = merge(regression_snp_list_correlated, snp_meta_data, by = "snp", all.x = TRUE)
    } else {
        regression_snp_list_correlated = data.table()
        regression_snp_list_no_correlated = regression_snp_list
        correlated_snps = paste0("Correlated snps not shown")
        # define dataframe containing correlated snps and meta data
        together_correlated = data.table()
    }

    # define dataframe containing NOT correlated snps and meta data
    together_not_correlated  = merge(regression_snp_list_no_correlated, snp_meta_data, by = "snp", all.x = TRUE)

    # find snps from other productivity group
    productivities = c("productive", "NOT_productive")
    prod = productivities[!productivities %in% c(productivity)]
    # load regression / pvalue results
    if (varying_int == "True"){
        product = ifelse(prod == 'productive', 'True', 'False')
        assign('regression_snp_list_other', as.data.table(read.table(generate_file_name(snp_id_list, trim_type, gene_type = paste0(substr(trim_type, 1, 1), '_gene'), productivity = product, gene_conditioning, weighting, condensing = 'by_gene', repetitions = 100), header = TRUE))[,c(1, 5)])
    } else {
        regression_snp_list_other = as.data.table(read.table(paste0("regression_bootstrap_results/", prod, "/", trim_type, "/", trim_type, "_", prod, "_snplist_", snp_id_list[1], "-",snp_id_list[length(snp_id_list)], "_snps_simple_condensing_by_patient.tsv"), header = TRUE))[,c(1, 5)]
    }

    colnames(regression_snp_list_other) = c("snp", "P_other")

    # combind correlated and not correlated snps
    together_temp = rbind(together_correlated, together_not_correlated, fill=TRUE)

    # find snps significant in other productivity group
    together_other= merge(together_temp, regression_snp_list_other, all.x = TRUE)
    together_other = together_other[P_other < bonferroni]

    # Use only data from indicated chromosome
    if (nrow(together_correlated) > 0){
        together_correlated = together_correlated[chr == chromosome]
    }
    together_not_correlated = together_not_correlated[chr == chromosome]
    together_other = together_other[chr == chromosome]

    par(mar=c(5,5,7,5)+.1)
    if (nrow(together_not_correlated) > 0){
        if (nrow(together_correlated) > 0){
            xlimit = c(min(together_not_correlated$hg19_pos, together_correlated$hg19_pos)-1000, max(together_not_correlated$hg19_pos, together_correlated$hg19_pos)+1000)
            ylimit = c(0, max(-1*log(together_not_correlated$pvalue, base =10), -1*log(together_correlated$pvalue, base =10))+5)
        } else {
            xlimit = c(min(together_not_correlated$hg19_pos)-1000, max(together_not_correlated$hg19_pos)+1000)
            ylimit = c(0, max(-1*log(together_not_correlated$pvalue, base =10))+5)
        }
        plot(as.numeric(together_not_correlated$hg19_pos), as.numeric(-1*log(together_not_correlated$pvalue, base =10)), bg = alpha("black", 0.4), col = alpha("black", 0.9), xlab = paste0('Chromosome ', chromosome,' position'), ylab = '-log10(p value)', main = paste0('P-value of SNP effect on ', trim_type, ' for ', productivity, " TCRs on chr", chromosome, "\n", gene_region), panel.first = grid(), cex.main=1.75, cex.lab=1.75, cex.axis=1.5, pch = 21, cex = 1.5, xlim = xlimit, ylim = ylimit)

        palette(brewer.pal(n = length(unique(together_correlated$cluster)), name = 'Set2'))
        col = setNames(palette(), levels(as.factor(together_correlated$cluster)))
        if (nrow(together_correlated) > 0){
            points(together_correlated$hg19_pos, -1*log(together_correlated$pvalue, base =10), bg = alpha(col[together_correlated$cluster], 0.8), col = alpha("black", 0.9), pch=21, cex = 1.5)
        }
        if (paired_productive_snps == 'True'){
            points(as.numeric(together_other$hg19_pos), -1*log(together_other$pvalue, base =10), col = alpha("red", 0.9), pch=1, cex = 1.5)
        }
    } else {
        if (nrow(together_correlated) > 0){
            xlimit = c(min(together_correlated$hg19_pos)-1000, max(together_correlated$hg19_pos)+1000)
            ylimit = c(0, max(-1*log(together_correlated$pvalue, base =10))+5)
        }
        palette(brewer.pal(n = length(unique(together_correlated$cluster)), name = 'Set2'))
        col = setNames(palette(), levels(as.factor(together_correlated$cluster)))

        plot(together_correlated$hg19_pos, -1*log(together_correlated$pvalue, base =10), bg = alpha(col[together_correlated$cluster], 0.8), col = alpha("black", 0.9), xlab = paste0('Chromosome ', chromosome,' position'), ylab = '-log10(p value)', main = paste0('P-value of SNP effect on ', trim_type, ' for ', productivity, " TCRs on ch", chromosome, "\n", gene_region), panel.first = grid(), cex.main=1.75, cex.lab=1.75, cex.axis=1.5, pch = 21, cex = 1.5, xlim = c(min(together_not_correlated$hg19_pos, together_correlated$hg19_pos)-1000, max(together_not_correlated$hg19_pos, together_correlated$hg19_pos)+1000), ylim = c(0, max(-1*log(together_not_correlated$pvalue, base =10), -1*log(together_correlated$pvalue, base =10))+5))
        
        if (paired_productive_snps == 'True'){
            points(as.numeric(together_other$hg19_pos), -1*log(together_other$pvalue, base =10), col = alpha("red", 0.9), pch=1, cex = 1.5)
        }
    }

    abline(h = -1*log(bonferroni, base = 10), col = "chartreuse3", lwd = 4, lty = 2)
    if (correlate_snps == 'True'){
        legend("topleft", box.lty=0, legend=c("-log10(bonferroni)", paste0("significant for ", prod, " TCRs"), "not in a cluster", rep("cluster", length(unique(together_correlated$cluster)))), col=c("chartreuse3", alpha("red", 0.9), alpha("black", 0.4), alpha(col, 0.8)), lty=c(2, NA, NA, rep(NA, length(unique(together_correlated$cluster)))), lwd = c(3, NA, NA,rep(NA, length(unique(together_correlated$cluster)))), pch = c(NA, 1, 19,rep(19, length(unique(together_correlated$cluster)))), cex = 1.5)
    } else if (correlate_snps == 'True'){
        legend("topleft", box.lty=0, legend=c("-log10(bonferroni)", paste0("significant for ", prod, " TCRs")), col=c("chartreuse3", alpha("red", 0.9)), lty=c(2, NA), lwd = c(3, NA), pch = c(NA, 1), cex = 1.5)
    } 
}



#compile_genotype_data <- function(snps_gds_file, snp_id){
#    # get snp genotype
#    genotype = snpgdsGetGeno(snps_gds_file, snp.id=snp_id, with.id = TRUE)
#    
#    # combine genotype, sample ids
#    genotypes_dt = data.frame(as.numeric(as.character(genotype$sample.id)), genotype$genotype)
#    colnames(genotypes_dt) = c("scanID", paste0("snp",snp_id))
#    # Convert subject names : 
#    subject_id_mapping = as.data.frame(read.table('/home/mrussel2/tcr-gwas/_ignore/snp_data/gwas_id_mapping.tsv', sep = "\t", fill=TRUE, header = TRUE, check.names = FALSE))
#    snps_genotypes = merge(genotypes_dt, subject_id_mapping, by = "scanID")
#    return(snps_genotypes[,c('localID', paste0("snp",snp_id))])
#}

#OLD VERSION:
trimming_snp_regression_weighted_varying_int_subject <- function(snps_dataframe, condensed_trimming_dataframe, productive, trim_type, repetitions){
    # remove warning messages (about singularity)
    control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))
    # set bonferroni correction to us the full group of snps from the gwas (regardless of how many we want to analyze)
    bonferroni = 0.05/35481497

    condensed_trimming_dataframe = as.data.table(condensed_trimming_dataframe)
    simple_regression_results = data.table()
    bootstrap_results = data.table()

    # subset trimming data to include only productive or not productive entires
    if (productive == "True"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "True"]
    } else if (productive == "False"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "False"]
    }
 
    # For each snpID:
    for (snpID in names(snps_dataframe)[-c(1)]){
        sub = NULL
        sub2 = NULL
        regression = NULL

        # merge snp data and trimming data
        sub = data.table(localID = snps_dataframe$localID, snp = snps_dataframe[[snpID]])
        sub2 = merge(sub, condensed_trimming_dataframe, by = "localID")
    
        # set regression formula given 
        form = formula(get(paste0(trim_type)) ~ snp + (1|localID))
        if (trim_type =='v_trim'){
            regression = lmer(formula = v_trim ~ snp + (1|localID), data = sub2[snp != "NA"], weights = weighted_v_gene_count, control=control)
        } else if (trim_type =='d0_trim'){
            regression = lmer(formula = d0_trim ~ snp + (1|localID), data = sub2[snp != "NA"], weights = weighted_d_gene_count, control=control)
        } else if (trim_type =='d1_trim'){
            regression = lmer(formula = d1_trim ~ snp + (1|localID), data = sub2[snp != "NA"], weights = weighted_d_gene_count, control=control)
        } else if (trim_type =='j_trim'){
            regression = lmer(formula = j_trim ~ snp + (1|localID), data = sub2[snp != "NA"], weights = weighted_j_gene_count, control=control)
        } else if (trim_type =='vj_insert'){
            regression = lmer(formula = vj_insert ~ snp + (1|localID), data = sub2[snp != "NA"], weights = weighted_vj_gene_count, control=control)
        } else if (trim_type =='dj_insert'){
            regression = lmer(formula = dj_insert ~ snp + (1|localID), data = sub2[snp != "NA"], weights = weighted_dj_gene_count, control=control)
        } else if (trim_type =='vd_insert'){
            regression = lmer(formula = vd_insert ~ snp + (1|localID), data = sub2[snp != "NA"], weights = weighted_vd_gene_count, control=control)
        }
        simple_regression_results = rbind(simple_regression_results, data.table(snp = snpID, intercept = mean(as.numeric(as.character(unlist(coefficients(regression)[[1]][1])))), slope = mean(as.numeric(as.character(unlist(coefficients(regression)[[1]][2]))))))
        
        se = clusboot_lmer(regression, data = sub2[snp != "NA"], cluster_variable = sub2[snp != "NA"]$localID, repetitions)[2,2]
        
        bootstrap_results = rbind(bootstrap_results, data.table(snp = snpID, standard_error = se))
    }
    together = bootstrap_regression_combine(bootstrap_results, simple_regression_results, bonferroni)
    return(together)
}

trimming_snp_regression_weighted <- function(snps_dataframe, condensed_trimming_dataframe, productive, trim_type, repetitions){
    condensed_trimming_dataframe = as.data.table(condensed_trimming_dataframe)
    simple_regression_results = data.table()
    bootstrap_results = data.table()
    if (productive == "True"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "True"]
    } else if (productive == "False"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "False"]
    }
    for (snpID in names(snps_dataframe)[-c(1)]){
        sub = NULL
        sub2 = NULL
        regression = NULL
        sub = data.table(localID = snps_dataframe$localID, snp = snps_dataframe[[snpID]])
        sub2 = merge(sub, condensed_trimming_dataframe, by = "localID")
    
        if (trim_type =='v_trim'){
            regression = glm(formula = v_trim ~ snp, data = sub2[snp != "NA"], weights = weighted_v_gene_count)
        } else if (trim_type =='d0_trim'){
            regression = glm(formula = d0_trim ~ snp, data = sub2[snp != "NA"], weights = weighted_d_gene_count)
        } else if (trim_type =='d1_trim'){
            regression = glm(formula = d1_trim ~ snp, data = sub2[snp != "NA"], weights = weighted_d_gene_count)
        } else if (trim_type =='j_trim'){
            regression = glm(formula = j_trim ~ snp, data = sub2[snp != "NA"], weights = weighted_j_gene_count)
        } else if (trim_type =='vj_insert'){
            regression = glm(formula = vj_insert ~ snp, data = sub2[snp != "NA"], weights = weighted_vj_gene_count)
        } else if (trim_type =='dj_insert'){
            regression = glm(formula = dj_insert ~ snp, data = sub2[snp != "NA"], weights = weighted_dj_gene_count)
        } else if (trim_type =='vd_insert'){
            regression = glm(formula = vd_insert ~ snp, data = sub2[snp != "NA"], weights = weighted_vd_gene_count)
        }
        simple_regression_results = rbind(simple_regression_results, data.table(snp = snpID, intercept = as.numeric(coef(regression)[1]), slope = as.numeric(coef(regression)[2])))
        se = bootstrap_cluster(sub2[snp != "NA"], repetitions, trim_type)
        bootstrap_results = rbind(bootstrap_results, data.table(snp = snpID, standard_error = se))
    }
    together = bootstrap_regression_combine(bootstrap_results, simple_regression_results)
    return(together)
}


#### BOOTSTRAP

subset_data_snp <- function(snpID, snps_dataframe, condensed_trimming_dataframe, productive){
    condensed_trimming_dataframe = as.data.table(condensed_trimming_dataframe)
    if (productive == "True"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "True"]
    } else if (productive == "False"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "False"]
    }
    sub = data.table(localID = snps_dataframe$localID, snp = snps_dataframe[[snpID]])
    sub2 = merge(sub, condensed_trimming_dataframe, by = "localID")
    return(sub2[snp != "NA"])
}

#model_coef <- function(data, index){
#    coef(glm(formula = v_trim ~ snp, data = data, weights = weighted_v_gene_count, subset = index))
#}
#
#bootstrap_coef <- function(data, repetitions){
#    boot = boot(data, model_coef, repetitions, weights = data$weighted_v_gene_count)
#    standard_error = sd(boot$t[,2])
#    return(standard_error)
#}

model_se <- function(data, index, gene_type){
    if (gene_type =='v_gene'){
        object = glm(formula = v_trim ~ snp, data = data, weights = weighted_v_gene_count, subset = index)
    } else if (gene_type =='d0_gene'){
        object = glm(formula = d0_trim ~ snp, data = data, weights = weighted_d_gene_count, subset = index)
    } else if (gene_type =='d1_gene'){
        object = glm(formula = d1_trim ~ snp, data = data, weights = weighted_d_gene_count, subset = index)
    } else if (gene_type =='j_gene'){
        object = glm(formula = j_trim ~ snp, data = data, weights = weighted_j_gene_count, subset = index)
    }
    coefs = summary(object)$coefficients
    coefs[,'Std. Error']
}

bootstrap_se <- function(data, repetitions, gene_type){
    if (gene_type =='v_gene'){
        weight = data$weighted_v_gene_count
    } else if (gene_type =='d0_gene' | gene_type =='d1_gene'){
        weight = data$weighted_d_gene_count
    } else if (gene_type =='j_gene'){
        weight = data$weighted_j_gene_count
    }
    boot = boot(data, model_se, repetitions, gene_type = gene_type, weights = weight, parallel="multicore", ncpus = 20)
    standard_error = colMeans(boot$t)[2]
    return(standard_error)
}


bootstrap_cluster <- function(data, repetitions, trim_type){
    if (trim_type =='v_trim'){
        boot = clusbootglm(v_trim ~ snp, data, clusterid = localID, B = repetitions, n.cores = 6)
    } else if (trim_type =='d0_trim'){
        boot = clusbootglm(d0_trim ~ snp, data, clusterid = localID, B = repetitions, n.cores = 6)
    } else if (trim_type =='d1_trim'){
        boot = clusbootglm(d1_trim ~ snp, data, clusterid = localID, B = repetitions, n.cores = 6)
    } else if (trim_type =='j_trim'){
        boot = clusbootglm(j_trim ~ snp, data, clusterid = localID, B = repetitions, n.cores = 6)
    } else if (trim_type =='vj_insert'){
        boot = clusbootglm(vj_insert ~ snp, data, clusterid = localID, B = repetitions, n.cores = 6)
    } else if (trim_type =='dj_insert'){
        boot = clusbootglm(dj_insert ~ snp, data, clusterid = localID, B = repetitions, n.cores = 6)
    } else if (trim_type =='vd_insert'){
        boot = clusbootglm(vd_insert ~ snp, data, clusterid = localID, B = repetitions, n.cores = 6)
    }
    standard_error = boot$boot.sds[2]
    return(standard_error)
}

bootstrap_cluster_lmer <- function(regression, repetitions){
    # preserves grouping from lmer regression
    subject_slope = function(.) {coef(.)$localID[,"snp"]}
    merBoot <- bootMer(regression, subject_slope, nsim = repetitions, re.form = NA, ncpus = 10, parallel = "multicore")
    standard_error <- mean(apply(merBoot$t, 2, sd))
    return(standard_error)
}


regression_weighted_bootstrap_se <- function(snps_dataframe, condensed_trimming_dataframe, productive, repetitions, gene_type){
    bootstrap_results = data.frame()
    for (snpID in names(snps_dataframe)[-c(1,ncol(snps_dataframe))]){
        data = subset_data_snp(snpID, snps_dataframe, condensed_trimming_dataframe, productive)
        se = bootstrap_cluster(data, repetitions, gene_type)
        bootstrap_results = rbind(bootstrap_results, data.frame(snp = snpID, standard_error = se))
    }
    return(bootstrap_results)
}


### bootstrap/regression on genotype data functions; 


run_snps_trimming <- function(snps_per_run, number_of_runs, trim_type, varying_int){
    # ALL SNP DATA
    # Get snp meta data
    snps_gds = openfn.gds("../_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds")
    n <- index.gdsn(snps_gds, "sample.id")
    sampleid <- read.gdsn(n) 
    sz <- snps_per_run
    #sz <- 100000
    nloop <- number_of_runs
    #nloop <- 355
    bigsize <- 35481497
    for ( i in 1:nloop ) {
        bootstrap_regression_results_v_gene_productive = data.table()
        bootstrap_regression_results_v_gene_NOT_productive = data.table()
        start <- (i-1)*sz + 1
        numrows <- min( sz, bigsize-start+1 ) 
        genotypes <- read.gdsn(index.gdsn(snps_gds, "genotype"), start=c(1,start), count = c(398, numrows))
        genotypes_df = data.frame(scanID = c(sampleid), genotypes)
        snpid <- read.gdsn(index.gdsn(snps_gds, "snp.id"), start = start, count=numrows)
        colnames(genotypes_df) = c("scanID", paste0("snp",snpid))
        genotypes_dt = as.data.table(genotypes_df)
        genotypes_dt$scanID = as.numeric(as.character(genotypes_dt$scanID))

        # Convert subject names and compile condensed data: 

        subject_id_mapping = as.data.table(read.table('../_ignore/snp_data/gwas_id_mapping.tsv', sep = "\t", fill=TRUE, header = TRUE, check.names = FALSE))

        snps = merge(genotypes_dt, subject_id_mapping, by = "scanID")
        print(paste0("snp data ", i, " of ", nloop, " compiled for ", trim_type))

        # Replace '3' genotypes with NA (missing genotype)
        snps[snps==3]<-NA
        snps = as.data.table(snps)

        # Get rid of snp columns which have the same entry for all entries and get ride of snp columns which have NA for all entries (missing for all individuals)
        snps_no_NA = Filter(function(x) length(unique(x))!=1, snps)
        snps_no_N2 = Filter(function(x) length(unique(snps_no_NA[x != "NA"])) != 1, snps_no_NA)
        snps_no_NA2 = data.frame(localID = snps_no_NA$localID)
        for (column in names(snps_no_NA)[-c(1, ncol(snps_no_NA))]){
            if (length(unique(na.omit(snps_no_NA[[column]]))) > 1){
                snps_no_NA2 = cbind(snps_no_NA2, snps_no_NA[[column]])
                names(snps_no_NA2)[names(snps_no_NA2) == 'snps_no_NA[[column]]'] <- paste0(column)
            }
        }
        if (trim_type == "v_trim"){
            trimming_data = v_trimming
        } else if (trim_type == "d0_trim" | trim_type == "d1_trim"){
            trimming_data = d_trimming2
        } else if (trim_type == "j_trim"){
            trimming_data = j_trimming
        } else if (trim_type == "vj_insert"){
            trimming_data = vj_insert
        } else if (trim_type == "dj_insert"){
            trimming_data = dj_insert2
        } else if (trim_type == "vd_insert"){
            trimming_data = vd_insert2
        }
        if (varying_int == "FALSE"){
            weighted_regression_bootstrap_patient_vgene_results_productive = trimming_snp_regression_weighted(snps_no_NA2, trimming_data, productive = "True", trim_type = trim_type, repetitions = 300)
            print("finished weighted_regression_patient_vgene_results_productive")

            weighted_regression_bootstrap_patient_vgene_results_NOT_productive = trimming_snp_regression_weighted(snps_no_NA2, trimming_data, productive = "False", trim_type = trim_type,  repetitions = 300)
            print("finished weighted_regression_patient_vgene_results_NOT_productive")
        } else if (varying_int == "TRUE"){
            weighted_regression_bootstrap_patient_vgene_results_productive = trimming_snp_regression_weighted_varying_int_subject(snps_no_NA2, trimming_data, productive = "True", trim_type = trim_type, repetitions = 300)
            print("finished weighted_regression_patient_vgene_results_productive")

            weighted_regression_bootstrap_patient_vgene_results_NOT_productive = trimming_snp_regression_weighted_varying_int_subject(snps_no_NA2, trimming_data, productive = "False", trim_type = trim_type,  repetitions = 300)
            print("finished weighted_regression_patient_vgene_results_NOT_productive")
        }   
        
        if (nrow(weighted_regression_bootstrap_patient_vgene_results_productive) != 0){
            bootstrap_regression_results_v_gene_productive = rbind(bootstrap_regression_results_v_gene_productive, weighted_regression_bootstrap_patient_vgene_results_productive)
            write.table(bootstrap_regression_results_v_gene_productive, file=paste0('regression_bootstrap_results/productive/', trim_type, '/', trim_type, '_productive_', i, '.tsv'), quote=FALSE, sep='\t', col.names = NA)
        }
        print("finished bootstrap_regression_weighted_productive")
        
        if (nrow(weighted_regression_bootstrap_patient_vgene_results_NOT_productive) != 0){
            bootstrap_regression_results_v_gene_NOT_productive = rbind(bootstrap_regression_results_v_gene_NOT_productive, weighted_regression_bootstrap_patient_vgene_results_NOT_productive)
            write.table(bootstrap_regression_results_v_gene_NOT_productive, file=paste0('regression_bootstrap_results/NOT_productive/', trim_type, '/', trim_type, '_NOT_productive_', i, '.tsv'), quote=FALSE, sep='\t', col.names = NA)
        }
        print("finished bootstrap_regression_weighted_NOT_productive")
      

        print(paste0("finished bootstrap for snp data ", i, " of ", nloop, " for ", trim_type))
    }

    closefn.gds(snps_gds)
}


run_snps_trimming_start_end <- function(snp_start, snp_end, trim_type, varying_int){
    # ALL SNP DATA
    # Get snp meta data
    snps_gds = openfn.gds("../_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds")
    n <- index.gdsn(snps_gds, "sample.id")
    sampleid <- read.gdsn(n) 
    bigsize <- 35481497
    bootstrap_regression_results_v_gene_productive = data.table()
    bootstrap_regression_results_v_gene_NOT_productive = data.table()
    numrows <- min( snp_end-snp_start, bigsize-snp_start+1 ) 
    genotypes <- read.gdsn(index.gdsn(snps_gds, "genotype"), start=c(1,snp_start), count = c(398, numrows))
    genotypes_df = data.frame(scanID = c(sampleid), genotypes)
    snpid <- read.gdsn(index.gdsn(snps_gds, "snp.id"), start = snp_start, count=numrows)
    colnames(genotypes_df) = c("scanID", paste0("snp",snpid))
    genotypes_dt = as.data.table(genotypes_df)
    genotypes_dt$scanID = as.numeric(as.character(genotypes_dt$scanID))
    # Convert subject names and compile condensed data: 
    subject_id_mapping = as.data.table(read.table('../_ignore/snp_data/gwas_id_mapping.tsv', sep = "\t", fill=TRUE, header = TRUE, check.names = FALSE))

    snps = merge(genotypes_dt, subject_id_mapping, by = "scanID")
    print(paste0("snp data compiled for ", snp_start, " to ", snp_end, " for ", trim_type))
    # Replace '3' genotypes with NA (missing genotype)
    snps[snps==3]<-NA
    snps = as.data.table(snps)
    # Get rid of snp columns which have the same entry for all entries and get ride of snp columns which have NA for all entries (missing for all individuals)
    snps_no_NA = Filter(function(x) length(unique(x))!=1, snps)
    snps_no_NA2 = Filter(function(x) length(unique(snps_no_NA[x != "NA"])) != 1, snps_no_NA)
    snps_no_NA2 = data.frame(localID = snps_no_NA$localID)
    for (column in names(snps_no_NA)[-c(1, ncol(snps_no_NA))]){
        if (length(unique(na.omit(snps_no_NA[[column]]))) > 1){
            snps_no_NA2 = cbind(snps_no_NA2, snps_no_NA[[column]])
            names(snps_no_NA2)[names(snps_no_NA2) == 'snps_no_NA[[column]]'] <- paste0(column)
        }
    }
    if (trim_type == "v_trim"){
        trimming_data = v_trimming
    } else if (trim_type == "d0_trim" | trim_type == "d1_trim"){
        trimming_data = d_trimming2
    } else if (trim_type == "j_trim"){
        trimming_data = j_trimming
    } else if (trim_type == "vj_insert"){
        trimming_data = vj_insert
    } else if (trim_type == "dj_insert"){
        trimming_data = dj_insert2
    } else if (trim_type == "vd_insert"){
        trimming_data = vd_insert2
    }
    if (varying_int == "FALSE"){
        weighted_regression_bootstrap_patient_vgene_results_productive = trimming_snp_regression_weighted(snps_no_NA2, trimming_data, productive = "True", trim_type = trim_type, repetitions = 300)
        print("finished weighted_regression_patient_vgene_results_productive")

        weighted_regression_bootstrap_patient_vgene_results_NOT_productive = trimming_snp_regression_weighted(snps_no_NA2, trimming_data, productive = "False", trim_type = trim_type,  repetitions = 300)
        print("finished weighted_regression_patient_vgene_results_NOT_productive")
    } else if (varying_int == "TRUE"){
        weighted_regression_bootstrap_patient_vgene_results_productive = trimming_snp_regression_weighted_varying_int_subject(snps_no_NA2, trimming_data, productive = "True", trim_type = trim_type, repetitions = 300)
        print("finished weighted_regression_patient_vgene_results_productive")

        weighted_regression_bootstrap_patient_vgene_results_NOT_productive = trimming_snp_regression_weighted_varying_int_subject(snps_no_NA2, trimming_data, productive = "False", trim_type = trim_type,  repetitions = 300)
        print("finished weighted_regression_patient_vgene_results_NOT_productive")
    }   
        
    
    if (nrow(weighted_regression_bootstrap_patient_vgene_results_productive) != 0){
        bootstrap_regression_results_v_gene_productive = rbind(bootstrap_regression_results_v_gene_productive, weighted_regression_bootstrap_patient_vgene_results_productive)
        write.table(bootstrap_regression_results_v_gene_productive, file=paste0('regression_bootstrap_results/productive/', trim_type, '/', trim_type, '_productive_',snp_start,'_', snp_end,'.tsv'), quote=FALSE, sep='\t', col.names = NA)
    }
    print("finished bootstrap_regression_weighted_productive")

    if (nrow(weighted_regression_bootstrap_patient_vgene_results_NOT_productive) != 0){
        bootstrap_regression_results_v_gene_NOT_productive = rbind(bootstrap_regression_results_v_gene_NOT_productive, weighted_regression_bootstrap_patient_vgene_results_NOT_productive)
        write.table(bootstrap_regression_results_v_gene_NOT_productive, file=paste0('regression_bootstrap_results/NOT_productive/', trim_type, '/', trim_type, '_NOT_productive_', snp_start,'_', snp_end, '.tsv'), quote=FALSE, sep='\t', col.names = NA)
    }
    print("finished bootstrap_regression_weighted_NOT_productive")

    print(paste0("finished bootstrap for snp data ", snp_start, " to ", snp_end, " for ", trim_type))

    closefn.gds(snps_gds)
}

# old version: 

run_snps_trimming_snp_list <- function(snp_id_list, trim_type, varying_int, repetitions){
    # Get snp meta data
    snps_gds = snpgdsOpen("../_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds")

    # get sample ids
    n <- index.gdsn(snps_gds, "sample.id")
    sampleid <- read.gdsn(n) 
    bigsize <- 35481497
    bootstrap_regression_results_v_gene_productive = data.table()
    bootstrap_regression_results_v_gene_NOT_productive = data.table()
    i = 0

    for (snp in snp_id_list){
        i = i + 1
        # get snp genotype
        genotype = snpgdsGetGeno(snps_gds, snp.id=snp)
        # combine genotype, sample ids
        genotypes_df = data.frame(scanID = c(sampleid), genotype)
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
        snps_no_NA = Filter(function(x) length(unique(x))!=1, snps)
        snps_no_NA = Filter(function(x) length(unique(snps_no_NA[x != "NA"])) != 1, snps_no_NA)
        snps_no_NA2 = data.frame(snps_no_NA$localID, snps_no_NA[[paste0("snp",snp)]])
        colnames(snps_no_NA2) = c("localID", paste0("snp",snp))

        # import condensed trimming data
        if (trim_type == "v_trim"){
            trimming_data = v_trimming
        } else if (trim_type == "d0_trim" | trim_type == "d1_trim"){
            trimming_data = d_trimming
        } else if (trim_type == "j_trim"){
            trimming_data = j_trimming
        } else if (trim_type == "vj_insert"){
            trimming_data = vj_insert
        } else if (trim_type == "dj_insert"){
            trimming_data = dj_insert
        } else if (trim_type == "vd_insert"){
            trimming_data = vd_insert
        }

        if (varying_int == "FALSE"){
            weighted_regression_bootstrap_patient_vgene_results_productive = trimming_snp_regression_weighted(snps_no_NA2, trimming_data, productive = "True", trim_type = trim_type, repetitions =  repetitions)
            print("finished weighted_regression_patient_vgene_results_productive")

            weighted_regression_bootstrap_patient_vgene_results_NOT_productive = trimming_snp_regression_weighted(snps_no_NA2, trimming_data, productive = "False", trim_type = trim_type,  repetitions =  repetitions)
            print("finished weighted_regression_patient_vgene_results_NOT_productive")
            prod_name = paste0('regression_bootstrap_results/productive/', trim_type, '/', trim_type, '_productive_snplist_', length(snp_id_list),'_snps.tsv')
            not_prod_name = paste0('regression_bootstrap_results/NOT_productive/', trim_type, '/', trim_type, '_NOT_productive_snplist_', length(snp_id_list),'_snps.tsv')
        } else if (varying_int == "TRUE"){
            weighted_regression_bootstrap_patient_vgene_results_productive = trimming_snp_regression_weighted_varying_int_subject(snps_no_NA2, trimming_data, productive = "True", trim_type = trim_type, repetitions =  repetitions)
            print("finished weighted_regression_patient_vgene_results_productive")

            weighted_regression_bootstrap_patient_vgene_results_NOT_productive = trimming_snp_regression_weighted_varying_int_subject(snps_no_NA2, trimming_data, productive = "False", trim_type = trim_type,  repetitions =  repetitions)
            print("finished weighted_regression_patient_vgene_results_NOT_productive")

            prod_name = paste0('regression_bootstrap_results/productive/', trim_type, '/', trim_type, '_productive_snplist_', length(snp_id_list),'_snps_varying_intercepts_by_subject.tsv')
            not_prod_name = paste0('regression_bootstrap_results/NOT_productive/', trim_type, '/', trim_type, '_NOT_productive_snplist_', length(snp_id_list),'_snps_varying_intercepts_by_subject.tsv')
        }   
        

        if (nrow(weighted_regression_bootstrap_patient_vgene_results_productive) != 0){
            bootstrap_regression_results_v_gene_productive = rbind(bootstrap_regression_results_v_gene_productive, weighted_regression_bootstrap_patient_vgene_results_productive)
        }
        print("finished bootstrap_regression_weighted_productive")

        if (nrow(weighted_regression_bootstrap_patient_vgene_results_NOT_productive) != 0){
            bootstrap_regression_results_v_gene_NOT_productive = rbind(bootstrap_regression_results_v_gene_NOT_productive, weighted_regression_bootstrap_patient_vgene_results_NOT_productive)
        }
        print("finished bootstrap_regression_weighted_NOT_productive")

        print(paste0("finished bootstrap for snp data for ", i, " of ", length(snp_id_list), " for ", trim_type))
    }
    write.table(bootstrap_regression_results_v_gene_productive, file= prod_name, quote=FALSE, sep='\t', col.names = NA)
    write.table(bootstrap_regression_results_v_gene_NOT_productive, file= not_prod_name, quote=FALSE, sep='\t', col.names = NA)
    closefn.gds(snps_gds)
}


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


# Weighted mixed model function

trimming_snp_regression_weighted_varying_int_subject <- function(snps_dataframe, condensed_trimming_dataframe, productive, trim_type, repetitions, bonferroni){
    varying_int = "True"
    # remove warning messages (about singularity)
    control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))
    # set bonferroni correction to us the full group of snps from the gwas (regardless of how many we want to analyze)
    bonferroni = 0.05/35481497

    condensed_trimming_dataframe = as.data.table(condensed_trimming_dataframe)
    simple_regression_results = data.table()
    bootstrap_results = data.table()

    # subset trimming data to include only productive or not productive entires
    if (productive == "True"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "TRUE"]
    } else if (productive == "False"){
        condensed_trimming_dataframe = condensed_trimming_dataframe[productive == "FALSE"]
    }
    
    # For each snpID:
    for (snpID in names(snps_dataframe)[-c(1)]){
        sub = NULL
        sub2 = NULL
        regression = NULL

     # set weights and gene type for regression
        if (trim_type =='d1_trim' | trim_type =='d0_trim'){
            weight = paste0("weighted_d_gene_count")
            gene_type = paste0("d_gene")
        } else {
            weight = paste0("weighted_", str_split(trim_type, "_")[[1]][1], "_gene_count")
            gene_type = paste0(str_split(trim_type, "_")[[1]][1], "_gene")
            if (str_split(trim_type, "_")[[1]][2] == 'insert'){
                gene_type1 = paste0(substr(trim_type, 1, 1), '_gene')
                gene_type2 = paste0(substr(trim_type, 2, 2), '_gene')
            }
        }

        # merge snp data and trimming data
        sub = data.table(localID = snps_dataframe$localID, snp = snps_dataframe[[snpID]])
        sub2 = merge(sub, condensed_trimming_dataframe, by = "localID")
    
        # set regression formula given 
        if (str_split(trim_type, "_")[[1]][2] == 'insert'){
            form = formula(get(paste0(trim_type)) ~ snp + get(paste0(gene_type1)) + get(paste0(gene_type2)) + (1|localID))
        } else {
            form = formula(get(paste0(trim_type)) ~ snp + get(paste0(gene_type))+ (1|localID))
        }
        

        # REGRESSION!
        regression = lmer(formula = form, data = sub2[snp != "NA"], weights = get(weight), control=control)
        
        # Calculate slope, intercept 
        # Add the Intercept term with a mean of the gene specific intercept
        intercept = summary(regression)$coefficients[,'Estimate']['(Intercept)'] + mean(summary(regression)$coefficients[,'Estimate'][-c(1,2)])
        slope = summary(regression)$coefficients[,'Estimate']['snp']

        if (slope != "NA"){
            # Pvalue screen before doing bootstrap (so that we only bootstrap things that may be significant)
            boot_screen = bootstrap_screen(regression)
            if (boot_screen[2]< (bonferroni*10)){
                bootstrap_results = calculate_pvalue(regression, data = sub2[snp != "NA"], cluster_variable = sub2[snp != "NA"]$localID, varying_int, repetitions)
                if (bootstrap_results[2]<bonferroni){
                    bootstrap_results = calculate_pvalue(regression, data = sub2[snp != "NA"], cluster_variable = sub2[snp != "NA"]$localID, varying_int, repetitions=1000)
                } else {
                    bootstrap_results = bootstrap_results
                }
            } else {
                bootstrap_results = boot_screen
            }
        } else {
            bootstrap_results = data.table()
        }
        together = cbind(data.table(snp = snpID, intercept = intercept, slope = slope), bootstrap_results)
        
        # Combine snpID, intercept, slope, etc.
        results = rbind(simple_regression_results, together)
    }
    return(results)
}

# combine bootstrap, results, and calculate p value
bootstrap_regression_combine <- function(bootstrap_results, regression_results, bonferroni, regression, data, cluster, repetitions){
    # If there was a valid regression, bootstrap
    if (nrow(bootstrap_results) != 0 & nrow(regression_results) != 0){
        # merge bootstrap and regression results
        together = merge(bootstrap_results, regression_results, by = "snp")
        # calculate zscore
        together$zscore = together$slope/together$standard_error
        # calculate two sided pvalue
        together$pvalue = 2*pnorm(-abs(together$zscore))
        # if p-value is significant, repeat bootstrap for 1000 repetitions, recalculate zscore, pvalue, etc.
        if(together$pvalue < bonferroni){
            se = clusboot_lmer(regression, data, cluster, trim_type, varying_int, weighting, repetitions)[2,2]
            together$standard_error = se
            together$zscore = together$slope/together$standard_error
            together$pvalue = 2*pnorm(-abs(together$zscore))
        }
    } else {
        together = data.table()
    }
    return(together)
}


# no data condensing at all...this function eneters
simple_trimming_snp_regression_no_condensing <- function(snps_dataframe, productive, trim_type){
    varying_int = "False"

    files = list.files(path="../_ignore/emerson_stats", pattern="*.tsv", full.names=TRUE)

    # For each snpID:
    for (snpID in names(snps_dataframe)[-c(1)]){
        sub = NULL
        sub2 = NULL
        regression = NULL

        # load snp data
        sub = data.table(localID = snps_dataframe$localID, snp = snps_dataframe[[snpID]])

        # skip patients until we get to one that has snp data associated...
        i = 0
        empty = 'True'
        while (empty == 'True'){
            i = i + 1
            file_name = str_split(files[i], "/")[[1]][4]
            file_root_name = str_split(file_name, ".tsv")[[1]][1]
            patient_id = str_split(file_root_name, "_")[[1]][3]
            if (patient_id %in% sub$localID == TRUE){
                empty = 'False'
            }
        }

        first_trimming_file = fread(files[i], sep = "\t", fill=TRUE, header = TRUE)
        first_trimming_file$localID = patient_id

        # subset trimming data to include only productive or not productive entires
        if (productive == "True"){
            first_trimming_file = first_trimming_file[productive == "TRUE"]
        } else if (productive == "False"){
            first_trimming_file = first_trimming_file[productive == "FALSE"]
        }

        sub2 = merge(sub, first_trimming_file, by = "localID")
    
        # set regression formula given 
        form = formula(get(paste0(trim_type)) ~ snp )
        
        # REGRESSION!
        regression = biglm(formula = form, data = sub2[snp != "NA"])

        count = 0
        for (file in files[-c(1:i)]){
            sub2 = NULL
            # read file...
            temp_file = fread(file, sep = "\t", fill=TRUE, header = TRUE)

            file_name = str_split(file, "/")[[1]][4]
            file_root_name = str_split(file_name, ".tsv")[[1]][1]
            patient_id = str_split(file_root_name, "_")[[1]][3]
            temp_file$localID = patient_id

            # skip iteration if patient does not have snp data
            if (patient_id %in% sub$localID == FALSE){
                next
            }

            # subset trimming data to include only productive or not productive entires
            if (productive == "True"){
                temp_file = temp_file[productive == "TRUE"]
            } else if (productive == "False"){
                temp_file = temp_file[productive == "FALSE"]
            }

            sub2 = merge(sub, temp_file, by = "localID")

            regression = update(object = regression, moredata = sub2[snp != "NA"])
            count = count + 1
        }

        # Calculate slope, intercept 
        # Add the Intercept term with a mean of the gene specific intercept
        intercept = coef(regression)['(Intercept)']
        slope = coef(regression)['snp']

        simple_regression_results = data.table()
        bootstrap_results = data.table()

        if (slope != "NA"){
            #NO BOOTSTRAP FOR THIS ANALYSIS
            # Pvalue screen before doing bootstrap (so that we only bootstrap things that may be significant)
            se = summary(regression)$mat[2,4]
            zscore = slope/se
            # calculate two sided pvalue
            pvalue = 2*pnorm(-abs(zscore))
            bootstrap_results = data.frame(standard_error = se, pvalue = pvalue)
        } else {
            bootstrap_results = data.table()
        }

        together = cbind(data.table(snp = snpID, intercept = intercept, slope = slope), bootstrap_results)
        
        # Combine snpID, intercept, slope, etc.
        results = rbind(simple_regression_results, together)
    }
    return(results)
}
