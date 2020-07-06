# Miscellaneous functions!

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






