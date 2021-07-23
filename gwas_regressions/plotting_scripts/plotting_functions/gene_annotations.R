library(data.table)
setDTthreads(1)
# hg19 gene annotations from UCSC genome browswer
source(paste0(PROJECT_PATH, '/tcr-gwas/gwas_regressions/plotting_scripts/plotting_functions/gene_annotation_functions.R'))

GENE_ANNOTATIONS = get_gene_locations()

dntt_features = get_gene_features('dntt', c('exon', 'promoter'))
artemis_features = get_gene_features('artemis', c('exon', 'promoter'))
tcrb_features = get_genes_in_tcrb(simplify = FALSE)
#add trbd2 and trbj1-x genes from liftover igmt ncbi36 to hg19...

GENE_FEATURES = rbind(dntt_features, artemis_features, tcrb_features, fill = TRUE)
GENE_FEATURES[, pos1 := as.numeric(pos1)]
GENE_FEATURES[, pos2 := as.numeric(pos2)]
