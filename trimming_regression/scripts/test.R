args = c(1000, 'vd_insert', 1, 0)

count = 1000

# Read in snp list
snp_data = snp_file_by_snp_start(snp_start = as.numeric(args[1]), count)
genotype_data = compile_all_genotypes(snp_start = as.numeric(args[1]), count)
genotype_data_filtered = remove_matrix_column_by_genotype(genotype_data) 

snp_list = snp_data
genotype_list = genotype_data_filtered
trim_type = args[2]

# type is 'insert' or 'trim'
type = strsplit(args[2], '_')[[1]][2]

condensing = ifelse(type == 'insert', 'by_patient', 'by_gene')
gene_conditioning = ifelse(type == 'insert', 'False', 'True')
random_effects = ifelse(type == 'insert', 'False', 'True')
d_infer = ifelse(type == 'insert', 'False', 'True')


# set pca_structure_correction variable
pca_structure_correction = ifelse(as.numeric(args[4]) == 0, 'False', 'True')

gene_type = 'same'
weighting = 'True'
repetitions = as.numeric(args[4])
write_table = 'False'
ncpus = as.numeric(args[3])
maf_cutoff = 0.05
