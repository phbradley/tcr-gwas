visualizations.html:visualizations.Rmd 
	R -e "rmarkdown::render('visualizations.Rmd')"

visualize_genome.html:visualize_genome.Rmd 
	R -e "rmarkdown::render('visualize_genome.Rmd')"

visualize_genome_boot_no_pca.html:visualize_genome.Rmd
   	R -e "rmarkdown::render('visualize_genome.Rmd', params = list(bootstrap_rerun = 100, pca = 'False'), output_file = 'visualize_genome_boot_no_pca.html')"
visualize_genome_boot_pca.html:visualize_genome.Rmd
   	R -e "rmarkdown::render('visualize_genome.Rmd', params = list(bootstrap_rerun = 100, pca = 'True'), output_file = 'visualize_genome_boot_pca.html')"

visualize_genome_no_boot-pca.html:visualize_genome.Rmd
   	R -e "rmarkdown::render('visualize_genome.Rmd', params = list(bootstrap_rerun = 'False', pca = 'False'), output_file = 'visualize_genome_no_boot-pca.html')"

visualize_specific_file.html:visualize_specific_file.Rmd 
	R -e "rmarkdown::render('visualize_specific_file.Rmd')"
