library("data.table")
library("lme4")
library("boot")
library("ClusterBootstrap")

bootstrap_screen <- function(regression){
    # extract standard error from regression object for snp slope
    se = summary(regression)$coefficients[,'Std. Error']['snp']
    slope = summary(regression)$coefficients[,'Estimate']['snp']

    # zscore calculation
    zscore = slope/se
    # calculate two sided pvalue
    pvalue = 2*pnorm(-abs(zscore))
    return(data.frame(standard_error = se, pvalue = pvalue))
}

# This one definitely does clustering (which is what we want I think)--this incorporates fixed and random effects!
clusboot_lmer <- function(regression, data, cluster_variable, repetitions){
    # Hide warnings for singularity
    control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))
    # forms cluster variable (here, extract patient names)
    clusters <- names(table(cluster_variable))
    # Set empty matrix to hold results
    standard_errors <- matrix(NA, nrow=repetitions, ncol=2)

    for(i in 1:repetitions){
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
        # Repeat regression for this data subset, calculate standard errors
        formula = formula(regression)   
        standard_errors[i,] <- c(colMeans(coef((lmer(formula, bootdat, control = control)))$localID[1]), colMeans(coef(lmer(formula, bootdat, control = control))$localID[2]))
    }

    # combine results into output
    # Note...the standard error of the intercept does NOT include variation by gene choice...but the snp standard error is CORRECT
    coefficients_se <- cbind(c(coef(summary(regression))[1], coef(summary(regression))[2]),apply(standard_errors,2,sd))
    colnames(coefficients_se) <- c("Estimate","Std. Error")
    return(coefficients_se)
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
            se = clusboot_lmer(regression, data, cluster, repetitions)[2,2]
            together$standard_error = se
            together$zscore = together$slope/together$standard_error
            together$pvalue = 2*pnorm(-abs(together$zscore))
        }
    } else {
        together = data.table()
    }
    return(together)
}


calculate_pvalue <- function(regression, data, cluster_variable, repetitions){
    # cluster bootstrap (clustered by individual)
    se = clusboot_lmer(regression, data, cluster_variable, repetitions)[2,2]
    slope = mean(as.numeric(as.character(unlist(coefficients(regression)[[1]][2]))))
    # zscore calculation
    zscore = slope/se
    # calculate two sided pvalue
    pvalue = 2*pnorm(-abs(zscore))
    return(data.frame(standard_error = se, pvalue = pvalue))
}


