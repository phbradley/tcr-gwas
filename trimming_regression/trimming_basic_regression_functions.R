library("lme4")
library("data.table")
library('biglm')
library(reticulate)
use_python('/home/mrussel2/miniconda3/envs/py/bin/python')

source("trimming_bootstrap_functions.R")

# fully condensed data (mean by patient)
simple_trimming_snp_regression <- function(snps_dataframe, condensed_trimming_dataframe, productive, trim_type, repetitions, bonferroni, python_test){
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
            #NO BOOTSTRAP FOR THIS ANALYSIS
            # Pvalue screen before doing bootstrap (so that we only bootstrap things that may be significant)
            se = summary(regression)$coefficients[,'Std. Error']['snp']
            zscore = slope/se
            # calculate two sided pvalue
            pvalue = 2*pnorm(-abs(zscore))
            bootstrap_results = data.frame(standard_error = se, pvalue = pvalue)
        } else {
            bootstrap_results = data.table()
        }

        if (python_test == 'True'){
            source_python('test_phil.py')
            pval_py = linear_reg_phil(sub2[snp != "NA"]$snp, sub2[snp != "NA"][[paste0(trim_type)]])
            bootstrap_results$pval_py = pval_py
        }

        #if (slope != "NA"){
        #    # Pvalue screen before doing bootstrap (so that we only bootstrap things that may be significant)
        #    boot_screen = bootstrap_screen(regression)
        #    if (boot_screen[2]< (bonferroni*10)){
        #        bootstrap_results = calculate_pvalue(regression, data = sub2[snp != "NA"], cluster_variable = sub2[snp != "NA"]$localID, varying_int, repetitions)
        #        if (bootstrap_results[2]<bonferroni){
        #            bootstrap_results = calculate_pvalue(regression, data = sub2[snp != "NA"], cluster_variable = sub2[snp != "NA"]$localID, varying_int, #repetitions=1000)
        #        } else {
        #            bootstrap_results = bootstrap_results
        #        }
        #    } else {
        #        bootstrap_results = boot_screen
        #    }
        #} else {
        #    bootstrap_results = data.table()
        #}
        together = cbind(data.table(snp = snpID, intercept = intercept, slope = slope), bootstrap_results)
        
        # Combine snpID, intercept, slope, etc.
        results = rbind(simple_regression_results, together)
    }
    return(results)
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

