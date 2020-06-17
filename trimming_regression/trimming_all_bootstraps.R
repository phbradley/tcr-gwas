library(data.table)
library(lme4)
library(boot)

source("trimming_bootstrap_functions.R")
source("trimming_all_regressions.R")

print("finished processing regressions")

bootstrap_weighted_productive = regression_weighted_bootstrap_se(snps_no_NA, trimming, productive = "True", repetitions = 500)
print("finished bootstrap_weighted_productive")

bootstrap_weighted_NOT_productive = regression_weighted_bootstrap_se(snps_no_NA, trimming, productive = "False", repetitions = 500)
print("finished bootstrap_weighted_NOT_productive")

bootstrap_regression_weighted_productive = bootstrap_regression_combine(bootstrap_weighted_productive, weighted_regression_patient_vgene_results_productive)
write.table(bootstrap_regression_weighted_productive, file='snps_subset_productive_v_gene_regression_bootstrap.tsv', quote=FALSE, sep='\t', col.names = NA)
print("finished bootstrap_regression_weighted_productive")

bootstrap_regression_weighted_NOT_productive = bootstrap_regression_combine(bootstrap_weighted_NOT_productive, weighted_regression_patient_vgene_results_NOT_productive)
write.table(bootstrap_regression_weighted_NOT_productive, file='snps_subset_NOT_productive_v_gene_regression_bootstrap.tsv', quote=FALSE, sep='\t', col.names = NA)
print("finished bootstrap_regression_weighted_NOT_productive")

print("finished all bootstraps")