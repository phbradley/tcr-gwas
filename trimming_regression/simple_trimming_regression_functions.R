library("lme4")
library('reticulate')
use_python('/home/mrussel2/miniconda3/envs/py/bin/python')

source("bootstrap_functions.R")

# fully condensed data (mean by patient)
simple_trimming_snp_regression <- function(snps_dataframe, condensed_trimming_dataframe, productive, trim_type, gene_type, weighting, gene_conditioning, python_test, snp_list){
    # subset trimming data to include only productive or not productive entires
    condensed_trimming_dataframe = filter_by_productivity(as.data.frame(condensed_trimming_dataframe), productive)

    if (ncol(condensed_trimming_dataframe) == 10){
        colnames(condensed_trimming_dataframe) = c('localID', 'productive', 'v_trim', 'd0_trim', 'd1_trim', 'j_trim', 'vj_insert', 'dj_insert', 'vd_insert', 'total_tcr')
    } else if (ncol(condensed_trimming_dataframe) == 7){
        colnames(condensed_trimming_dataframe) = c('localID', 'productive', 'v_trim', 'd0_trim', 'd1_trim', 'j_trim', 'total_tcr')
    } 
    
    # For each snpID:
    snpID=names(snps_dataframe)[-c(1)]
    
    # merge snp data and trimming data
    sub = data.frame(localID = snps_dataframe$localID, snp = snps_dataframe[[snpID]])
    sub2 = merge(sub, condensed_trimming_dataframe, by = "localID")
    snps_trimming_data = sub2 %>% filter(!is.na(snp))

    # define trim type and gene type
    if (gene_type == 'same'){
        gene_type = paste0(substr(trim_type, 1, 1), '_gene')
        weight = paste0("weighted_", substr(trim_type, 1, 1), '_gene_count')
    }
    gene_type = paste0(gene_type)
    if (str_split(trim_type, "_")[[1]][2] == 'insert'){
        gene_type1 = paste0(substr(trim_type, 1, 1), '_gene')
        gene_type2 = paste0(substr(trim_type, 2, 2), '_gene')
        weight = paste0("weighted_", substr(trim_type, 1, 2), '_gene_count')
    }

    # set regression formula given 
    form = formula(get(paste0(trim_type)) ~ snp )

    if (gene_conditioning == 'True'){
        if (str_split(trim_type, "_")[[1]][2] == 'insert'){
            form = update(form, ~ . + get(paste0(gene_type1)) + get(paste0(gene_type2)))
        } else {
            form = update(form, ~ . + get(paste0(gene_type)))
        }
    } else {
        weight = 'total_tcr'
    }

    # REGRESSION!
    if (weighting == 'True'){
        regression = lm(formula = form, data = snps_trimming_data, weights = get(weight))
    } else {
        regression = lm(formula = form, data = snps_trimming_data)
    }
    
        
    # Calculate slope, intercept 
    # Add the Intercept term with a mean of the gene specific intercept
    intercept = summary(regression)$coefficients[,'Estimate']['(Intercept)']
    slope = summary(regression)$coefficients[,'Estimate']['snp']
    se = summary(regression)$coefficients[,'Std. Error']['snp']
    pvalue = summary(regression)$coefficients[,'Pr(>|t|)']['snp']
    bootstrap_results = data.frame(standard_error = se, pvalue = pvalue)


    if (python_test == 'True'){
        source_python('test_phil.py')
        pval_py = linear_reg_phil(sub2[snp != "NA"]$snp, sub2[snp != "NA"][[paste0(trim_type)]])
        bootstrap_results$pval_py = pval_py
    }

    snp_list$snp = paste0('snp', snp_list$snp)
    results_temp = merge(snp_list, data.frame(snp = snpID, intercept = intercept, slope = slope), by = 'snp')

    together = cbind(results_temp, bootstrap_results)
    return(together)
}




