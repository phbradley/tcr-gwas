# TCR-GWAS 

# Install
Everything R 4.0.2 based. R packages that are required can be installed via [miniconda](https://docs.conda.io/en/latest/miniconda.html): 

```bash 
cd gwas-regressions/
conda env create -f environment.yml
conda activate r
```

# Analysis outline: 
1. Download data into the directory `tcr-gwas/_ignore/` (add sources and/or write script to do this)
2. Edit [config](config/) files to be project and/or computer specific
2. Run the genome-wide GWAS analysis for a phenotype or phenotype class of interest
    * You can run [submit_phenotype.sh](submit_phenotype.sh) or [submit_phenotype_class.sh](submit_phenotype_class.sh). Both scripts take three arguments: 
        1. phenotype (i.e. v_trim, vd_insert, etc.) or phenotype class (i.e. trimming, insertion, etc.) 
        2. cluster partition (these scripts are set up to run with a slurm cluster scheduler)
        3. number of cluster cores per job
    * i.e. `bash submit_phenotype.sh v_trim partition-name 2` will run the pipeline for v_trim by submitting 2 core jobs to the cluster partition named "partition-name"
    * A full list of possible phenotypes and phenotype classes is listed below
3. Plot results [genome-wide](plotting_scripts/plot_genome_wide.Rmd) or for a [specific gene](plotting_scripts/plot_gene.Rmd) using the Rmarkdown notebooks







