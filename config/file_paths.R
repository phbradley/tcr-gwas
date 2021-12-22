# SNP genotype file (download https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001918.v1.p1)
SNP_GDS_FILE <<- paste0(PROJECT_PATH, "/tcr-gwas/data/downloaded_data/HSCT_comb_geno_combined_v03_tcr.gds")

# File to map scanIDs to localID (data located within the `gwas_id_mapping.tsv` file available at https://doi.org/10.5281/zenodo.5719520)
ID_MAPPING_FILE <<- paste0(PROJECT_PATH, '/tcr-gwas/data/downloaded_data/gwas_id_mapping.tsv')

# PCA file (data located within the `all_pc_air.txt` file available at https://doi.org/10.5281/zenodo.5719520)
PCA_FILE <<- paste0(PROJECT_PATH, '/tcr-gwas/data/downloaded_data/all_pc_air.txt')
# PCA variance file (data located within the `all_pc_air_variance.txt` file available at https://doi.org/10.5281/zenodo.5719520)
PCA_VARIANCE_FILE <<- paste0(PROJECT_PATH, '/tcr-gwas/data/downloaded_data/all_pc_air_variance.txt')

# Parsed TCR repertoire data (data located within the `emerson_parsed_TCRB.tgz` file available at https://doi.org/10.5281/zenodo.5719520)
TCR_REPERTOIRE_DATA_DIRECTORY = paste0(PROJECT_PATH, '/tcr-gwas/data/downloaded_data/emerson_data')

# TCRB CDR3 gene file (data located within the `human_vj_allele_cdr3_nucseqs.tsv` file available at https://doi.org/10.5281/zenodo.5719520)
CDR3_GENE_ASSIGNMENT_FILE = paste0(PROJECT_PATH, '/tcr-gwas/data/downloaded_data/human_vj_allele_cdr3_nucseqs.tsv')

# snp meta data file (data located within the `emerson_snp_rs_data.tsv` file available at https://doi.org/10.5281/zenodo.5719520) 
SNP_META_DATA_FILE = paste0(PROJECT_PATH, '/tcr-gwas/data/downloaded_data/emerson_snp_rs_data.tsv')

# D gene allele genotypes (data located within the `emerson_trbd2_alleles.tsv` file available at https://doi.org/10.5281/zenodo.5719520)
D_ALLELES = paste0(PROJECT_PATH, '/tcr-gwas/data/downloaded_data/emerson_trbd2_alleles.tsv')
