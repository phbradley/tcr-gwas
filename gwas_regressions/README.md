# TCR-GWAS 
TODO add description

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
3. Plot results [genome-wide](plotting_scripts/plot_genome_wide.Rmd) or for a [specific gene](plotting_scripts/plot_gene.Rmd) using the Rmarkdown notebooks

# About the genome-wide analysis

(TODO expand)
With this model, we want to measure the effect of genome wide SNPs on certain TCR repertoire features. 

Specifically, you can run the analysis one of two ways: 
1. Using a specific phenotype (using [submit_phenotype.sh](submit_phenotype.sh) as described above) from the options below: 
    * v_trim
    * v_trim_naive
    * v_trim_zero_trimming_fraction
    * v_pnucs_count
    * v_pnucs_fraction
    * v_pnucs_fraction_zero_trimming_subset
    * d0_trim
    * d0_trim_naive
    * d0_trim_zero_trimming_fraction
    * d0_pnucs_count
    * d0_pnucs_fraction
    * d0_pnucs_fraction_zero_trimming_subset
    * d1_trim
    * d1_trim_naive
    * d1_trim_zero_trimming_fraction
    * d1_pnucs_count
    * d1_pnucs_fraction
    * d1_pnucs_fraction_zero_trimming_subset
    * j_trim
    * j_trim_naive
    * j_trim_zero_trimming_fraction
    * j_pnucs_count
    * j_pnucs_fraction
    * j_pnucs_fraction_zero_trimming_subset
    * dj_insert
    * vd_insert
    * vj_insert
    * total_insert
   
2. Using a group of these phenotypes called a "phenotype class" (using [submit_phenotype_class.sh](submit_phenotype_class.sh) as described above) from the options below:
    * trimming
        * v_trim
        * d0_trim
        * d1_trim
        * j_trim
    * trimming_naive
        * v_trim_naive
        * d0_trim_naive
        * d1_trim_naive
        * j_trim_naive
    * insertion
        * dj_insert
        * vd_insert
        * vj_insert
        * total_insert
    * p-addition_count
        * v_pnucs_count
        * d0_pnucs_count
        * d1_pnucs_count
        * j_pnucs_count
    * p-addition_fraction
        * v_pnucs_fraction
        * d0_pnucs_fraction
        * d1_pnucs_fraction
        * j_pnucs_fraction
    * p-addition_fraction_trimming_subset
        * v_pnucs_fraction_zero_trimming_subset
        * d0_pnucs_fraction_zero_trimming_subset
        * d1_pnucs_fraction_zero_trimming_subset
        * j_pnucs_fraction_zero_trimming_subset
    * zero_trimming_fraction
        * v_trim_zero_trimming_fraction
        * d0_trim_zero_trimming_fraction
        * d1_trim_zero_trimming_fraction
        * j_trim_zero_trimming_fraction




