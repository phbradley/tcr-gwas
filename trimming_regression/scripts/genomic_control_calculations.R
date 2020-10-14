library(data.table)
library(tidyverse)
setDTthreads(threads = 1)

genomic_control_calculation <- function(data){

    productive_pval_median = median(data[productivity == 'productive']$pvalue)
    NOTproductive_pval_median = median(data[productivity == 'NOT_productive']$pvalue)

    # Here, we can use the pvalue median to calculate a correction factor where
    # we first convert the median pvalue to a chi-squared value (with 1 degree
    # of freedom). To calculate the correction factor, we divide the median
    # chi-squared value by the chi-squared value of a pvalue of 0.5 (which is
    # what the median should be if there is no population structure)
    
    correction_factor_productive = qchisq(productive_pval_median, df = 1, lower.tail = FALSE)/qchisq(0.5, df = 1, lower.tail = FALSE) 

    correction_factor_NOT_productive = qchisq(NOTproductive_pval_median, df = 1, lower.tail = FALSE)/qchisq(0.5, df = 1, lower.tail = FALSE) 

    # Now, we can calculate a chi-squared value for each pvalue in the dataset
    data$chisq = qchisq(data$pvalue, df = 1, lower.tail = FALSE)
    
    # Using the chi-squared values and the correction factor (which is on the
    # chi-squared scale), we can re-scale the chi-squared values by dividing by
    # the correction factor
    data[productivity == 'productive', chisq_corrected := chisq/correction_factor_productive]

    data[productivity == 'NOT_productive', chisq_corrected := chisq/correction_factor_NOT_productive]

    # Now, with the corrected chi-squared values for each row, we can
    # convert the chi-squared values back to pvalues!

    data$pvalue_genomic_control_correction = pchisq(data$chisq_corrected, df = 1, lower.tail = FALSE) 
    
    return(data)
}
    

