all:
	R -e "rmarkdown::render('pca.Rmd')"
	#R -e "rmarkdown::render('pca.Rmd', output_format = 'pdf_document')"

continuous:
	while :; do inotifywait -e modify *.Rmd; make; done

.PHONY: continuous
