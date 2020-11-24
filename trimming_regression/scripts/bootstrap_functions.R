library("lme4")
library('RhpcBLASctl')
omp_set_num_threads(1)
blas_set_num_threads(1)

# This file houses functions related to the boostrapping protocol. 

# This function calculates the bootstrap-free p-value 
# (so the pvalue from the initial regression is calculated using the variance, covariance matrix from the regression)
bootstrap_screen <- function(regression, RANDOM_EFFECTS){
    # extract standard error from regression object for snp slope
    se = sqrt(diag(vcov(regression)))[2]
    slope = ifelse(RANDOM_EFFECTS == 'True', fixef(regression)['snp'], coef(regression)['snp'])
    # zscore calculation
    zscore = slope/se
    # calculate two sided pvalue
    pvalue = 2*pnorm(-abs(zscore))
    return(data.frame(standard_error = se, pvalue = pvalue, bootstraps = 0))
}

# This one definitely does clustering (which is what we want I think)--this incorporates fixed and random effects!
clusboot_lmer <- function(regression, data, cluster_variable, trim_type, RANDOM_EFFECTS, GENE_CONDITIONING, WEIGHTING, REPETITIONS){
    # forms cluster variable (here, extract patient names)
    clusters <- names(table(cluster_variable))
    # Set empty matrix to hold results
    standard_errors <- matrix(NA, nrow=REPETITIONS, ncol=2)

   # set regression weights 
    if (WEIGHTING == 'True'){
        stopifnot(str_split(trim_type, "_")[[1]][2] == 'trim')
        weight = paste0("weighted_", substr(trim_type, 1, 1), '_gene_count')
    } 

    for(i in 1:REPETITIONS){
        # Sample clusters (patients) with replacement
        index <- sample(1:length(clusters), length(clusters), replace=TRUE)
        # Assign clusters to index
        clusters_by_index <- clusters[index]
        # Create a contingency table for the clusters
        contingency_table_clusters <- table(clusters_by_index)
        bootdat <- NULL

        # For each cluster group (based on how many repeats in contingency table)
        for(j in 1:max(contingency_table_clusters)){
            # subset data to include only the cluster variables indicated in j
            data_subset <- data[cluster_variable %in% names(contingency_table_clusters[contingency_table_clusters %in% j]),]
            # Add data to bootdata table
            for(k in 1:j){
                bootdat <- rbind(bootdat, data_subset)
            }
        }
        if (length(unique(bootdat$snp)) <= 1){
            next
        }
        # Repeat regression for this data subset, calculate standard errors
        formula = formula(regression) 
        if (RANDOM_EFFECTS == 'True'){
            # Hide warnings for singularity
            control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))

            if (WEIGHTING == 'True'){
                regression_temp = eval(bquote(lmer(formula, bootdat, weights = .(as.name(weight)), control = control)))
                standard_errors[i,] <- c(colMeans(coef((regression_temp))$localID[1]), colMeans(coef(regression_temp)$localID[2]))
            } else {
                regression_temp =lmer(formula, bootdat, control = control)
                standard_errors[i,] <- c(colMeans(coef(regression_temp)$localID[1]), colMeans(coef(regression_temp)$localID[2]))
            }
        } else {
            if (WEIGHTING == 'True'){
                regression_temp = eval(bquote(lm(formula = formula, data = bootdat, weights = .(as.name(weight)))))
                standard_errors[i,] <- c(coef(regression_temp)[1], coef(regression_temp)[2])
            } else {
                regression_temp = lm(formula = formula, data = bootdat)   
                standard_errors[i,] <- c(coef(regression_temp)[1], coef(regression_temp)[2])
            }
        }
    }
    if (RANDOM_EFFECTS == 'True'){
        # combine results into output
        # Note...the standard error of the intercept does NOT include variation by gene choice...but the snp standard error is CORRECT
        coefficients_se <- cbind(c(fixef(regression)['(Intercept)'], fixef(regression)['snp']),apply(na.omit(standard_errors),2,sd))
        colnames(coefficients_se) <- c("Estimate","Std. Error")
    } else {
        coefficients_se <- cbind(c(coef(summary(regression))[1], coef(summary(regression))[2]),apply(na.omit(standard_errors),2,sd))
        colnames(coefficients_se) <- c("Estimate","Std. Error")
    }
    return(coefficients_se)
}

# calculates pvalue from bootstrap
calculate_pvalue <- function(regression, data, cluster_variable, trim_type, RANDOM_EFFECTS, GENE_CONDITIONING, WEIGHTING, REPETITIONS){
    # cluster bootstrap (clustered by individual)
    se = clusboot_lmer(regression, data, cluster_variable, trim_type, RANDOM_EFFECTS, GENE_CONDITIONING, WEIGHTING, REPETITIONS)[2,2]
    slope = ifelse(RANDOM_EFFECTS == 'True', fixef(regression)['snp'], coef(regression)['snp'])

    # zscore calculation
    zscore = slope/se
    # calculate two sided pvalue
    pvalue = 2*pnorm(-abs(zscore))
    return(data.frame(standard_error = se, pvalue = pvalue, bootstraps = REPETITIONS))
}


