
get_gene_locations <- function(){
    # genes_common_names = c('artemis', 'mhc', 'dntt', 'rag', 'tcrb', 'tcra', 'tcrb_v_dj_zoom', 'tcrb_dj_zoom')
    genes_common_names = c('artemis', 'mhc', 'dntt', 'rag', 'tcrb', 'tcra')

    # genes = c('DCLRE1C', 'MHC', 'DNTT', 'RAG', 'TCRB', 'TCRA', 'TCRB_zoom', 'TCRB_zoom')
    genes = c('DCLRE1C', 'MHC', 'DNTT', 'RAG', 'TCRB', 'TCRA') 

    # chr = c(10, 6, 10, 11, 7, 14, 7, 7)
    chr = c(10, 6, 10, 11, 7, 14)

    # pos1 = c(14948871, 25912984, 98064085, 36589563, 141998851, 22090057, 142336000, 142448742)
    pos1 = c(14948871, 25912984, 98064085, 36589563, 141998851, 22090057)

    # pos2 = c(14996094, 33290793, 98098321, 36619829, 142510972, 23021075, 142504000, 142504000)
    pos2 = c(14996094, 33290793, 98098321, 36619829, 142510972, 23021075)
    
    annotations = data.table(gene_common_name = genes_common_names, gene_locus = genes, chr = chr, pos1 = pos1, pos2 = pos2)

    return(annotations)
}

get_gene_features <- function(gene_name, features){
    stopifnot(gene_name %in% c('dntt', 'artemis'))
    features_dt = data.table()
    gene_location = get_gene_locations()
    gene = gene_location[gene_common_name == gene_name]

    for (feature in features){
        stopifnot(feature %in% c('promoter', 'intron', 'exon'))
        temp_feature_dt = fread(paste0(PROJECT_PATH, '/tcr-gwas/data/', gene_name, '_', feature, '.txt'))[, 1:3]
        colnames(temp_feature_dt) = c('chr', 'pos1', 'pos2')
        temp_feature_dt$name = feature
        temp_feature_dt$chr = substring(temp_feature_dt$chr, 4)
        if (feature == 'intron'){
            temp_feature_dt = temp_feature_dt[pos1 >= gene$pos1 & pos2 <= gene$pos2]
        }
        features_dt = rbind(features_dt, temp_feature_dt, fill = TRUE)
    }
    features_dt$gene_locus = gene$gene_locus
    features_dt$gene_common_name = gene$gene_common_name
    return(features_dt)
}

get_genes_in_tcrb <- function(simplify){
    tcrb_features = fread(paste0(PROJECT_PATH, '/tcr-gwas/data/tcrb_features_gencode_24_37.txt'), fill = TRUE)[, c('chrom', 'txStart', 'txEnd', 'name2')]
    colnames(tcrb_features) = c('chr', 'pos1', 'pos2', 'gene')
    tcrb_features$chr = as.numeric(substring(tcrb_features$chr, 4))
    tcrb_features = tcrb_features[substring(gene, 1, 3) == 'TRB']
    tcrb_features = tcrb_features[gene != 'TRBC1' & gene != 'TRBC2']

    if (simplify == TRUE){
        tcrb_features = tcrb_features[gene != 'TRBV30']
        tcrb_features[substring(gene, 4, 4) == 'V', name := 'total span of V-gene (except TRBV30)']
        tcrb_features[substring(gene, 4, 4) == 'D', name := gene]
        tcrb_features[substring(gene, 4, 4) == 'J', name := paste0('total span of ', substring(gene, 1, 5), '-genes')]
        simple_tcrb = tcrb_features[, min(pos1), by = .(chr,name)]
        colnames(simple_tcrb) = c('chr', 'name', 'pos1')
        tcrb_features = merge(simple_tcrb, tcrb_features[, max(pos2), by = .(chr, name)])
        colnames(tcrb_features) = c('chr', 'name', 'pos1', 'pos2')
    } else {
        tcrb_features[substring(gene, 4, 4) == 'V', name := 'V-genes']
        tcrb_features[substring(gene, 4, 4) == 'D', name := 'D-genes']
        tcrb_features[substring(gene, 4, 4) == 'J', name := 'J-genes']
    }
    tcrb_features$gene_common_name = 'tcrb'
    tcrb_features$gene_locus = 'TCRB'
    return(tcrb_features)
}
