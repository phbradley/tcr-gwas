trim_type = 'all_trim'

# type is 'insert' or 'trim'
type = strsplit(trim_type, '_')[[1]][2]
bootstrap_count = 0
condensing = ifelse(type == 'insert', 'by_patient', 'by_gene')
gene_conditioning = ifelse(type == 'insert', 'False', 'True')
random_effects = ifelse(type == 'insert', 'False', 'True')
bootstrap_rerun_count = 100
d_infer = ifelse(type == 'insert', 'False', 'True')


# set pca_structure_correction variable
pca_structure_correction = 'True'

maf_cutoff = 0.05
