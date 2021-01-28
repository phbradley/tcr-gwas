# TCR-GWAS 
The goal of this project is to identify genetic biases within T-cell receptor repertoires by testing whether any genetic variations genome-wide are associated with repertoire feature biases (i.e. increased trimming, N-insertion, P-addition, etc.). 

# Install
Everything R 4.0.2 based. R packages that are required can be installed via [miniconda](https://docs.conda.io/en/latest/miniconda.html): 

```bash 
cd gwas-regressions/
conda env create -f environment.yml
conda activate r
```

# Analysis outline: 
1. Download data into the directory `tcr-gwas/_ignore/` (TODO add sources and/or write script to do this)
2. Edit [config](config/) files to be project and/or computer specific
2. Run the genome-wide GWAS analysis for a phenotype or phenotype class of interest
    * You can run [submit_phenotype.sh](submit_phenotype.sh) or [submit_phenotype_class.sh](submit_phenotype_class.sh). Both scripts take three arguments: 
        1. phenotype (i.e. v_trim, vd_insert, etc.) or phenotype class (i.e. trimming, insertion, etc.) 
        2. cluster partition (these scripts are set up to run with a slurm cluster scheduler)
        3. number of cluster cores per job
    * i.e. `bash submit_phenotype.sh v_trim partition-name 2` will run the pipeline for v_trim by submitting 2 core jobs to the cluster partition named "partition-name"
    * A full list of possible phenotypes and phenotype classes is listed below
3. Plot results across the [genome](plotting_scripts/plot_genome_wide.Rmd) or for a specific [gene](plotting_scripts/plot_gene.Rmd) using the Rmarkdown notebook.
    * Also, have a look at the plotting [README](plotting_scripts/README.md) for more details.

# About the genome-wide analysis

With this model, we want to measure the effect of genome wide SNPs on certain TCR repertoire features. 
See the manuscript for specific model details. (TODO add ms details) 

Specifically, you can run the analysis one of two ways: 
1. Using a specific phenotype (using [submit_phenotype.sh](submit_phenotype.sh) as described above) from the options in the "TCR Phenotypes in Class" column:
2. Using a group of these phenotypes called a "phenotype class" (using [submit_phenotype_class.sh](submit_phenotype_class.sh) as described above) from the "Phenotype Class":
 
   
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
| zero_trimming_fraction              | v_trim_zero_trimming_fraction          | Fraction of TCRs with zero V-gene trimming (with population structure correction)                              |
|                                     | d0_trim_zero_trimming_fraction         | Fraction of TCRs with zero 5'-end-D-gene trimming (with population structure correction)                       |
|                                     | d1_trim_zero_trimming_fraction         | Fraction of TCRs with zero 3'-end-D-gene trimming (with population structure correction)                       |
|                                     | j_trim_zero_trimming_fraction          | Fraction of TCRs with zero J-gene trimming (with population structure correction)                              |
| p-addition_count                    | v_pnucs_count                          | Number of V-gene P-nucleotides (with population structure correction)                                          |
|                                     | d0_pnucs_count                         | Number of 5'-end-D-gene P-nucleotides (with population structure correction)                                   |
|                                     | d1_pnucs_count                         | Number of 3'-end-D-gene P-nucleotides (with population structure correction)                                   |
|                                     | j_pnucs_count                          | Number of J-gene P-nucleotides (with population structure correction)                                          |
| p-addition_fraction                 | v_pnucs_fraction                       | Fraction of TCRs with V-gene P-nucleotides (with population structure correction)                              |
|                                     | d0_pnucs_fraction                      | Fraction of TCRs with 5'-end-D-gene P-nucleotides (with population structure correction)                       |
|                                     | d1_pnucs_fraction                      | Fraction of TCRs with 3'-end-D-gene P-nucleotides (with population structure correction)                       |
|                                     | j_pnucs_fraction                       | Fraction of TCRs with J-gene P-nucleotides (with population structure correction)                              |
| p-addition_fraction_trimming_subset | v_pnucs_fraction_zero_trimming_subset  | Fraction of non-V-gene-trimmed TCRs with V-gene P-nucleotides (with population structure correction)           |
|                                     | d0_pnucs_fraction_zero_trimming_subset | Fraction of non-5'-D-gene-trimmed TCRs with 5'-end-D-gene P-nucleotides (with population structure correction) |
|                                     | d1_pnucs_fraction_zero_trimming_subset | Fraction of non-3'-D-gene-trimmed TCRs with 3'-end-D-gene P-nucleotides (with population structure correction) |
|                                     | j_pnucs_fraction_zero_trimming_subset  | Fraction of non-J-gene-trimmed TCRs with J-gene P-nucleotides (with population structure correction)           |
| insertion                           | dj_insert                              | Number of D-J nucleotide insertions (with population structure correction)                                     |
|                                     | vd_insert                              | Number of V-D nucleotide insertions (with population structure correction)                                     |
|                                     | total_insert                           | Number of V-D and D-J nucleotide insertions (with population structure correction)                             |




