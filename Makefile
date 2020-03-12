all:
	R -e "rmarkdown::render('pca.Rmd')"

continuous:
	while :; do inotifywait -e modify *.Rmd; make; done

.PHONY: continuous
