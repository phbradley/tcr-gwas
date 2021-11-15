library(data.table)
setDTthreads(1)
# hg19 gene annotations from UCSC genome browswer
source(paste0(PROJECT_PATH, '/tcr-gwas/plotting_scripts/plotting_functions/gene_annotation_functions.R'))

GENE_ANNOTATIONS = get_gene_locations()

dntt_features = get_gene_features('dntt', c('exon', 'promoter'))
artemis_features = get_gene_features('artemis', c('exon', 'promoter'))
tcrb_features = get_genes_in_tcrb(simplify = FALSE)
#add trbd2 and trbj1-x genes from liftover igmt ncbi36 to hg19...

GENE_FEATURES = rbind(dntt_features, artemis_features, tcrb_features, fill = TRUE)
GENE_FEATURES[, pos1 := as.numeric(pos1)]
GENE_FEATURES[, pos2 := as.numeric(pos2)]

# from IMGT and Adaptive (https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-015-0238-z)
TRB_NONPROD_ALLELE_GENES = c('TRBV6-4', 'TRBV12-5', 'TRBV7-3', 'TRBV11-1', 'TRBV11-3', 'TRBV10-1', 'TRBV30') 
