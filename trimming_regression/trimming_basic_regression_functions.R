library("lme4")
library("data.table")

source("trimming_bootstrap_functions.R")

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


simple_trimming_snp_regression <- function(snps_dataframe, condensed_trimming_dataframe, productive, trim_type, repetitions, bonferroni){
    varying_int = "False"
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

    # get an individual mean for trimming....
    condensed_trimming_dataframe = condensed_trimming_dataframe[, .(mean(v_trim), mean(d0_trim), mean(d1_trim), mean(j_trim), mean(vj_insert), mean(dj_insert), mean(vd_insert)), by = .(localID, productive)]
    colnames(condensed_trimming_dataframe) = c('localID', 'productive', 'v_trim', 'd0_trim', 'd1_trim', 'j_trim', 'vj_insert', 'dj_insert', 'vd_insert')
 
    # For each snpID:
    for (snpID in names(snps_dataframe)[-c(1)]){
        sub = NULL
        sub2 = NULL
        regression = NULL

        # merge snp data and trimming data
        sub = data.table(localID = snps_dataframe$localID, snp = snps_dataframe[[snpID]])
        sub2 = merge(sub, condensed_trimming_dataframe, by = "localID")
    
        # set regression formula given 
        form = formula(get(paste0(trim_type)) ~ snp )

        # REGRESSION!
        regression = glm(formula = form, data = sub2[snp != "NA"])
        
        # Calculate slope, intercept 
        # Add the Intercept term with a mean of the gene specific intercept
        intercept = summary(regression)$coefficients[,'Estimate']['(Intercept)']
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