library(data.table)
library(ggplot2)

calculate_gene_frequency <- function(file){
    totals = file[,.N, by = patient_id]
    colnames(totals) = c('patient_id', 'total')
    for (gene in c('v_gene', 'd_gene', 'j_gene')){
        cols = c('patient_id', gene)
        assign(paste0('dt_', gene), file[,.N, by = cols])
        assign('temp', merge(get(paste0('dt_', gene)), totals, by = 'patient_id'))
        temp[[paste0(gene, '_frequency')]] = temp$N/temp$total
        temp = temp[,c(1,2,5)]
        file = merge(file, temp, by = c('patient_id', gene))
    }
    return(file)
} 


plot_insertions_distribution_by_gene <- function(gene_type, insertion_type, by_allele, frequency_filter){
    if (insertion_type %in% c('vj_insert', 'vd_insert', 'dj_insert')){
        file = fread(file = '/fh/fast/matsen_e/shared/tcr-gwas/insertion_data/insertions_by_patient.tsv')
    } else if (insertion_type %in% c('v_trim', 'd0_trim', 'd1_trim', 'j_trim')){
        file = fread(file = '/fh/fast/matsen_e/shared/tcr-gwas/insertion_data/trim_by_patient.tsv')
    }
    file[productive == 'TRUE', productivity_status := 'productive']
    file[productive == 'FALSE', productivity_status := 'NOT_productive']
    # leave out unassigned d gene cases
    file = file[d_gene != '-']
    
    if (frequency_filter > 0){
        file = calculate_gene_frequency(file)
    }  

    for (gene in gene_type){
        for (insert in insertion_type){
            if (by_allele == 'False'){
                file$gene_no_allele = substr(file[[gene]],1,nchar(file[[gene]])-3)
                gene_type_plot = 'gene_no_allele'
            } else {
                gene_type_plot = gene
            }

            file = file[get(paste0(gene,'_frequency')) > frequency_filter]
    
            ggplot(file, aes_string(x = paste(insert), color = paste(gene_type_plot))) + facet_wrap(~productivity_status, ncol = 1, nrow = 2) + geom_density(aes(y = ..scaled..), size = 0.5, adjust = 6) + ggtitle(paste0(insert, ' distribution by ', gene)) + theme_bw() + theme(title=element_text(size=20, hjust=0.5), axis.title=element_text(size=18), axis.text = element_text(size=18), legend.text = element_text(size=16), strip.text = element_text(face="bold", size=18)) 

            ggsave(paste0('figures/insertion_distribution_', insert, '_by_', gene, '.png'), plot = last_plot(), height=8, width=16, units="in", dpi = 500) 
        }
    }
}
