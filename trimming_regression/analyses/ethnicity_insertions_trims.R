library(data.table)
setDTthreads(1)
library(ggplot2)

ethnicity = fread(file = '/home/mrussel2/tcr-gwas/_ignore/snp_data/ethnicity_data.csv')
insertions = fread(file = '/fh/fast/matsen_e/shared/tcr-gwas/insertion_data/insertions_by_patient.tsv')
trims = fread(file = '/fh/fast/matsen_e/shared/tcr-gwas/insertion_data/trim_by_patient.tsv')

together_inserts = merge(insertions[-1], ethnicity, by.x = 'patient_id', by.y = 'id')[,-c(2:5)]
together_trims = merge(trims[-1], ethnicity, by.x = 'patient_id', by.y = 'id')[,-c(2:5)]

together_trims = together_trims[,total_trim := v_trim + d0_trim + d1_trim + j_trim]
together_inserts = together_inserts[,total_trim := vd_insert + dj_insert + vj_insert]

together_trims = together_trims[,lapply(.SD, mean), by = .(patient_id, productive, race)]
together_inserts = together_inserts[,lapply(.SD, mean), by = .(patient_id, productive, race)]

together_trims = together_trims[productive == 'TRUE', productivity_status := 'productive']
together_trims = together_trims[productive == 'FALSE', productivity_status := 'NOT_productive']

together_inserts = together_inserts[productive == 'TRUE', productivity_status := 'productive']
together_inserts = together_inserts[productive == 'FALSE', productivity_status := 'NOT_productive']



for (data in c('together_trims', 'together_inserts')){
    df = get(data)
    title = ifelse(data == 'together_trims', 'Total Trimming Length Average by Ethnicity', 'Total Insertion Length Average by Ethnicity')
    filename = ifelse(data == 'together_trims', 'figures/trim_average_by_ethnicity.png', 'figures/insert_average_by_ethnicity.png')
    average_all = df[, mean(total_trim), by = .(productivity_status)]
    colnames(average_all) = c('productivity_status', 'total_avg')

    ggplot(df, aes(x=race, y=total_trim)) +
        facet_grid(cols = vars(productivity_status)) +
        geom_boxplot() +
        geom_jitter(shape=16, position=position_jitter(0.2), size = 3, alpha = 0.5) +
        theme_classic() + 
        theme(text = element_text(size = 30)) +
        ggtitle(title) +
        geom_hline(data = average_all, aes(yintercept = total_avg))

    ggsave(filename, plot = last_plot(), width = 20, height = 8, units = 'in', dpi = 500)
}


