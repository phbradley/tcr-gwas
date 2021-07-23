# SNP genotype file
SNP_GDS_FILE <<- paste0(PROJECT_PATH, "/tcr-gwas/_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds")

# File to map scanIDs to localID
ID_MAPPING_FILE <<- paste0(PROJECT_PATH, '/tcr-gwas/_ignore/snp_data/gwas_id_mapping.tsv')

# PCA file
PCA_FILE <<- paste0(PROJECT_PATH, '/tcr-gwas/_ignore/snp_data/pc_pcair_08Nov2020.txt')

# TCR repertoire directory path
TCR_REPERTOIRE_DATA_DIRECTORY = paste0(PROJECT_PATH, '/tcr-gwas/_ignore/emerson_stats/')

# CDR3 sequence to gene file
CDR3_GENE_ASSIGNMENT_FILE = paste0(PROJECT_PATH, '/tcr-gwas/_ignore/human_vj_allele_cdr3_nucseqs.tsv')

# TCR diversity measure file path (optional)
TCR_DIV_FILE = paste0(PROJECT_PATH, '/tcr-gwas/_ignore/emerson_tcrdiv_scores.tsv')

# snp meta data file
SNP_META_DATA_FILE = paste0('/fh/fast/matsen_e/shared/tcr-gwas/emerson_snp_data.csv')

#rsid meta data file 
RSIDS = paste0('/fh/fast/matsen_e/shared/tcr-gwas/emerson_snp_rs_data.tsv')
