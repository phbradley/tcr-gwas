library(data.table)
library(lme4)
library(boot)

source("trimming_bootstrap_functions.R")
#source("trimming_regressions.R")

bootstrap_weighted_productive = regression_weighted_bootstrap_se(subset_snps_no_NA, trimming, productive = "True", repetitions = 500)

print("finished bootstrap_regression_weighted_productive")

bootstrap_weighted_NOT_productive = regression_weighted_bootstrap_se(subset_snps_no_NA, trimming, productive = "False", repetitions = 500)

print("finished bootstrap_regression_weighted_NOT_productive")

bootstrap_regression_weighted_productive = bootstrap_regression_combine(bootstrap_weighted_productive, weighted_regression_patient_vgene_results_productive)
print("combined_1")

bootstrap_regression_weighted_NOT_productive = bootstrap_regression_combine(bootstrap_weighted_NOT_productive, weighted_regression_patient_vgene_results_NOT_productive)
print("combined_2")