library(data.table)
library(lme4)
library(boot)

source("trimming_bootstrap_functions.R")

print("loaded functions")

source("trimming_all_regressions.R")

print("finished processing regressions")

bootstrap_weighted_v_gene_productive = regression_weighted_bootstrap_se(snps_no_NA, trimming, productive = "True", repetitions = 500, gene_type = 'v_gene')
print("finished bootstrap_weighted_productive")

bootstrap_weighted_v_gene_NOT_productive = regression_weighted_bootstrap_se(snps_no_NA, trimming, productive = "False", repetitions = 500, gene_type = 'v_gene')
print("finished bootstrap_weighted_NOT_productive")

bootstrap_regression_weighted_productive = bootstrap_regression_combine(bootstrap_weighted_v_gene_productive, weighted_regression_patient_vgene_results_productive)
print("finished bootstrap_regression_weighted_productive")

bootstrap_regression_weighted_NOT_productive = bootstrap_regression_combine(bootstrap_weighted_v_gene_NOT_productive, weighted_regression_patient_vgene_results_NOT_productive)
print("finished bootstrap_regression_weighted_NOT_productive")

print("finished all bootstraps")