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
clusboot_lmer <- function(regression, data, cluster_variable, trim_type, varying_int, weighting, repetitions){
    # forms cluster variable (here, extract patient names)
    clusters <- names(table(cluster_variable))
    # Set empty matrix to hold results
    standard_errors <- matrix(NA, nrow=repetitions, ncol=2)

    # set regression weights
    if (trim_type =='d1_trim' | trim_type =='d0_trim'){
        weight = paste0("weighted_d_gene_count")
    } else {
        weight = paste0("weighted_", str_split(trim_type, "_")[[1]][1], "_gene_count")
    }

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
        if (varying_int == 'True'){
            # Hide warnings for singularity
            control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))

            if (weighting == 'True'){
                regression_temp = lmer(formula, bootdat, weights = get(weight), control = control)
                standard_errors[i,] <- c(colMeans(coef((regression_temp))$localID[1]), colMeans(coef(regression_temp)$localID[2]))
            } else {
                regression_temp =lmer(formula, bootdat, control = control)
                standard_errors[i,] <- c(colMeans(coef(regression_temp)$localID[1]), colMeans(coef(regression_temp)$localID[2]))
            }
            # combine results into output
            # Note...the standard error of the intercept does NOT include variation by gene choice...but the snp standard error is CORRECT
            coefficients_se <- cbind(c(coef(summary(regression))[1], coef(summary(regression))[2]),apply(standard_errors,2,sd))
            colnames(coefficients_se) <- c("Estimate","Std. Error")
        } else {
            standard_errors[i,] <- c(coef(glm(formula = formula, data = bootdat))[1], coef(glm(formula = formula, data = bootdat))[2])
            coefficients_se <- cbind(c(coef(summary(regression))[1], coef(summary(regression))[2]),apply(standard_errors,2,sd))
            colnames(coefficients_se) <- c("Estimate","Std. Error")
        }
    }
    return(coefficients_se)
}


calculate_pvalue <- function(regression, data, cluster_variable, trim_type, varying_int, weighting, repetitions){
    # cluster bootstrap (clustered by individual)
    se = clusboot_lmer(regression, data, cluster_variable, trim_type, varying_int, weighting, repetitions)[2,2]
    if (varying_int == 'True'){
        slope = mean(as.numeric(as.character(unlist(coefficients(regression)[[1]][2]))))
    } else {
        slope = coefficients(regression)['snp']
    }
    # zscore calculation
    zscore = slope/se
    # calculate two sided pvalue
    pvalue = 2*pnorm(-abs(zscore))
    return(data.frame(standard_error = se, pvalue = pvalue))
}


