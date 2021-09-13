# SNP genotype file (download https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001918.v1.p1)
SNP_GDS_FILE <<- paste0(PROJECT_PATH, "/tcr-gwas/data/downloaded_data/HSCT_comb_geno_combined_v03_tcr.gds")

# File to map scanIDs to localID (download from zenodo--TODO, put link)
ID_MAPPING_FILE <<- paste0(PROJECT_PATH, '/tcr-gwas/data/downloaded_data/gwas_id_mapping.tsv')

# PCA file (download from zenodo--TODO, put link)
PCA_FILE <<- paste0(PROJECT_PATH, '/tcr-gwas/data/downloaded_data/all_pc_air.txt')

# TCR repertoire directory path (download https://doi.org/10.21417/B7001Z)
TCR_REPERTOIRE_DATA_DIRECTORY = paste0(PROJECT_PATH, '/tcr-gwas/data/downloaded_data/emerson_data')

# CDR3 sequence to gene file
CDR3_GENE_ASSIGNMENT_FILE = paste0(PROJECT_PATH, '/tcr-gwas/data/human_vj_allele_cdr3_nucseqs.tsv')

# snp meta data file (TODO write function to generate this)
SNP_META_DATA_FILE = paste0('/fh/fast/matsen_e/shared/tcr-gwas/emerson_snp_rs_data.tsv')

# D gene allele genotypes (download from zenodo--TODO, put link)
D_ALLELES = paste0(PROJECT_PATH, '/tcr-gwas/data/downloaded_data/emerson_trbd2_alleles.tsv')
