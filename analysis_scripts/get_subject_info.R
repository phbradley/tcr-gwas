# This script will run regressions for all snps in the indicated snp-chunk
source('config/config.R')
source(paste0(PROJECT_PATH, '/tcr-gwas/config/file_paths.R'))

library(data.table)
setDTthreads(1)
library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)
library(foreach)
library(doParallel)
library(tidyverse)

subjects = fread(ID_MAPPING_FILE) 
meta = fread(paste0(PROJECT_PATH, '/tcr-gwas/data/downloaded_data/emerson_subject_meta_data.tsv'))

together = merge(subjects, meta, by.x = 'localID', by.y = 'sample_name')

together_good = together %>% mutate(meta = str_split(as.character(sample_tags), ",")) %>% unnest(meta) %>% as.data.table()

# collapse HLA
together_good[meta %like% 'HLA', HLA := paste(meta, collapse = ','), by = localID]
together_good = together_good %>% group_by(localID) %>% fill(HLA, .direction='up') %>% fill(HLA) %>% as.data.table()

together_good[meta %like% 'Years', age := meta]
together_good = together_good %>% group_by(localID) %>% fill(age, .direction='up') %>% fill(age) %>% as.data.table()

together_good[meta %like% 'Cyto', CMV := meta]
together_good = together_good %>% group_by(localID) %>% fill(CMV, .direction='up') %>% fill(CMV) %>% as.data.table()

together_good[meta %like% 'ale', sex := meta]
together_good = together_good %>% group_by(localID) %>% fill(sex, .direction='up') %>% fill(sex) %>% as.data.table()

together_good[(meta %like% 'Hispanic') | (meta %like% 'Ethnicity'), ethnicity := meta]
together_good = together_good %>% group_by(localID) %>% fill(ethnicity, .direction='up') %>% fill(ethnicity) %>% as.data.table()

together_good[meta %in% c(" Unknown racial group", " Caucasian", " Native American or Alaska Native", "Caucasian", " Asian or Pacific Islander", " African Race"), race := meta]
together_good = together_good %>% group_by(localID) %>% fill(race, .direction='up') %>% fill(race) %>% as.data.table()

# together_good[meta %like% 'Childhood', age := meta]

cleaned_meta = unique(together_good[, -c('meta')])

# remove grouped ages
cleaned = cleaned_meta[!(age %like% "\\(") & !(age %like% "\\-") & !(age %like% "\\+")]

# clean ages
cleaned$age = trimws(cleaned$age)
cleaned$age = str_remove_all(cleaned$age, ' Years')
cleaned$age = as.numeric(cleaned$age)

#clean race and ethnicity
cleaned$race = trimws(cleaned$race)
cleaned$ethnicity = trimws(cleaned$ethnicity)

# group_ages 
# less than 10
cleaned[age <= 10, grouped_age := '< 10 Years']
# 11-20
cleaned[age > 10 & age <= 20, grouped_age := '11-20 Years']

# 21-30
cleaned[age > 20 & age <= 30, grouped_age := '21-30 Years']

# 31-40
cleaned[age > 30 & age <= 40, grouped_age := '31-40 Years']

# 41-50
cleaned[age > 40 & age <= 50, grouped_age := '41-50 Years']

# 51-60
cleaned[age > 50 & age <= 60, grouped_age := '51-60 Years']

# 60+
cleaned[age > 60, grouped_age := '> 60 Years']


