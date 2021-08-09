# TCR-GWAS 
The goal of this project is to identify genetic biases within T-cell receptor repertoires by testing whether any genome-wide genetic variations are associated with TCR repertoire feature biases (i.e. increased trimming, N-insertion, P-addition, etc.). 

# Install
Everything R 4.0.3 based. R packages that are required can be installed via [miniconda](https://docs.conda.io/en/latest/miniconda.html): 

```bash 
conda env create -f environment.yml
conda activate tcr-gwas 
```

# Requirements: 
Due to the computationally intensive nature of running these analyses genome-wide, this pipeline requires a cluster to run. 
These scripts are written specifically for a cluster set up to use the Slurm job scheduler. 
(Some minor modifications to the [single job submission script](submission_scripts/run_regressions_genome_wide.sh) and the [genome-wide submission script](submission_scripts/submit_cluster_jobs_continuously.sh) could allow this pipeline to be run using a different cluster workload manager.) 

# Analysis outline: 
1. Download data into the directory `tcr-gwas/data/downloaded_data` 
2. Edit [config](config/) files to be project and/or computer specific
2. Run the genome-wide GWAS analysis for a phenotype or phenotype class of interest
    * You can run [submit_phenotype.sh](submit_phenotype.sh) or [submit_phenotype_class.sh](submit_phenotype_class.sh). Both scripts take three arguments: 
        1. phenotype (i.e. v_trim, vd_insert, etc.) or phenotype class (i.e. trimming, insertion, etc.) 
        2. cluster partition (these scripts are set up to run with a slurm cluster scheduler)
        3. number of cluster cores per job
    * i.e. `bash submit_phenotype.sh v_trim partition-name 2` will run the pipeline for v_trim by submitting 2 core jobs to the cluster partition named "partition-name"
    * A full list of possible phenotypes and phenotype classes is listed below
3. Plot [figures](plotting_scripts/final_plots) from the manuscript.
    * Also, have a look at the plotting [README](plotting_scripts/README.md) for more details.
4. Run conditional analysis to identify independent snp associations within specified gene regions and phenotypes
    * To do this, you can run [run_conditional_analysis.sh](run_conditional_analysis.sh) either locally or by submitting to a cluster. 
    * This script takes the following arguments:
        1. gene (i.e. artemis, dntt, etc.)
        2. phenotype (i.e. v_trim, vd_insert, etc.)
        3. number of cluster cores
5. Run analysis for Nicaraguan validation cohort (two snps). Since this analysis is only for two snps, it can easily be run locally. To run this, you can run `Rscript analysis_scripts/validation_cohort_analysis/execute_validation_data_regressions.R [phenotype] [number of cores] [TCR chain]` which can take the arguments:
    1. phenotype (i.e. v_trim, vd_insert, etc.)
    2. number of cluster cores
    3. TCR chain (either `beta` or `alpha`)

**Note: all output files will be located at the indicated `OUTPUT_PATH` as specified in the [config](config) files**

# About the genome-wide analysis

With this model, we want to measure the effect of genome-wide SNPs on certain TCR repertoire features. 
See the manuscript for specific model details. (TODO add ms details) 

Specifically, you can run the analysis one of two ways: 
1. Using a specific phenotype (using [submit_phenotype.sh](submit_phenotype.sh) as described above) from the options in the "TCR Phenotypes in Class" column:
2. Using a group of these phenotypes called a "phenotype class" (using [submit_phenotype_class.sh](submit_phenotype_class.sh) as described above) from the "Phenotype Class" column:
 
   
| Phenotype Class                     | TCR Phenotypes in Class                | Phenotype Description                                                                                          |
|-------------------------------------|----------------------------------------|----------------------------------------------------------------------------------------------------------------|
| trimming                            | v_trim                                 | Amount of V-gene trimming (with population structure correction)                                               |
|                                     | d0_trim                                | Amount of 5'-end-D-gene trimming (with population structure correction)                                        |
|                                     | d1_trim                                | Amount of 3'-end-D-gene trimming (with population structure correction)                                        |
|                                     | j_trim                                 | Amount of 3'-end-D-gene trimming (with population structure correction)                                        |
| trimming_naive                      | v_trim_naive                           | Amount of V-gene trimming (without population structure correction)                                            |
|                                     | d0_trim_naive                          | Amount of 5'-end-D-gene trimming (without population structure correction)                                     |
|                                     | d1_trim_naive                          | Amount of 3'-end-D-gene trimming (without population structure correction)                                     |
|                                     | j_trim_naive                           | Amount of 3'-end-D-gene trimming (without population structure correction)                                     |
| p-addition_count                    | v_pnucs_count                          | Number of V-gene P-nucleotides (with population structure correction)                                          |
|                                     | d0_pnucs_count                         | Number of 5'-end-D-gene P-nucleotides (with population structure correction)                                   |
|                                     | d1_pnucs_count                         | Number of 3'-end-D-gene P-nucleotides (with population structure correction)                                   |
|                                     | j_pnucs_count                          | Number of J-gene P-nucleotides (with population structure correction)                                          |
| p-addition_fraction_trimming_subset | v_pnucs_fraction_zero_trimming_subset  | Fraction of non-V-gene-trimmed TCRs with V-gene P-nucleotides (with population structure correction)           |
|                                     | d0_pnucs_fraction_zero_trimming_subset | Fraction of non-5'-D-gene-trimmed TCRs with 5'-end-D-gene P-nucleotides (with population structure correction) |
|                                     | d1_pnucs_fraction_zero_trimming_subset | Fraction of non-3'-D-gene-trimmed TCRs with 3'-end-D-gene P-nucleotides (with population structure correction) |
|                                     | j_pnucs_fraction_zero_trimming_subset  | Fraction of non-J-gene-trimmed TCRs with J-gene P-nucleotides (with population structure correction)           |
| insertion                           | dj_insert                              | Number of D-J nucleotide insertions (with population structure correction)                                     |
|                                     | vd_insert                              | Number of V-D nucleotide insertions (with population structure correction)                                     |
|                                     | total_insert                           | Number of V-D and D-J nucleotide insertions (with population structure correction)                             |




