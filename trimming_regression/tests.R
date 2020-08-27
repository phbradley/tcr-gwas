source("run_bootstrap_regression_all_snps_functions.R")


# Read in snp list
snp_list = as.data.table(read.table("../_ignore/phil_sig_snps/top100_trimming_pvals_for_maggie.tsv", sep = "\t", fill=TRUE, header = TRUE))[-c(1:4)]

snp_data = snp_list[11]

# Run regression/bootstrap
run_snps_trimming_snp_list(snp_id_list = unique(snp_data$snp), trim_type = 'v_trim', condensing = 'by_gene', gene_conditioning = 'True', weighting = 'True', repetitions = 0, write_table = 'False')

# for run_snps_trimming_snp_list:
snp_id_list = unique(snp_data$snp)
trim_type = 'v_trim'
condensing = 'gene_cross'
gene_type = 'j_gene'
gene_conditioning = 'True'
weighting = 'True'
repetitions = 0
write_table = 'False'

# for trimming_regression

snps_dataframe = snp_genotypes
condensed_trimming_dataframe = trimming_data
productive = "True"
bootstrap_repetitions = repetitions

run_snps_trimming_snp_list(snp_id_list = unique(snp_data$snp), trim_type = 'v_trim', gene_type = 'j_gene', condensing = 'gene_cross', gene_conditioning = 'True', weighting = 'True', repetitions = 0, write_table = 'False')


source("manha_visualization.R")
source("find_correlated_snps.R")
source("make_snp_file.R")

pre_regression_snp_list = as.data.table(read.table("../_ignore/phil_sig_snps/top100_trimming_pvals_for_maggie.tsv", sep = "\t", fill=TRUE, header = TRUE))[-c(1:4)]

phil_snp_list = combine_trimming_regression_dfs(pre_regression_snp_list, condensing = 'by_gene', trimming_type_list = c('v_trim', 'd0_trim', 'd1_trim', 'j_trim', "vd_insert", "vj_insert", "dj_insert"))

artemis = snp_file(chromosome = 10, position1= 14397359, position2= 15454432)
tdt = snp_file(chromosome = 10, position1= 95804409, position2= 96838564)
mhc = snp_file(chromosome = 6, position1= 0, position2= 4050069)
rag = snp_file(chromosome = 11, position1= 36010709, position2= 37093156)

product = "productive"
gene = "artemis"

#correlated_snps_list = correlate_snps_ld(chromosome, artemis, cutoff = 0.85)

trim = "v_trim"

snp_meta_data = get(paste(gene))
snp_id_list = unique(get(gene)$snp)
correlated_snps = 'False'
paired_productive_snps = 'False'
productivity = product
varying_int = 'True'
trim_type = trim
chromosome = unique(get(gene)$chr)
correlate_snps = "False"
gene_conditioning = 'True'
weighting = 'True'
gene_region = gene

