library(data.table)
library(plyr)
library(readr)
library(stringr)
library(gdsfmt)
library(SNPRelate)
library(GWASTools)

setwd("../_ignore/")

# Get snp meta data
gfile = openfn.gds("snp_data/HSCT_comb_geno_combined_v03_tcr.gds")

n <- index.gdsn(gfile, "sample.id")
sampleid <- read.gdsn(n) 

write.table(sampleid, file='big_subject_data.tsv', quote=FALSE, sep='\t', col.names = NA)

n <- index.gdsn(gfile, "snp.id")
snpid <- read.gdsn(n) #, start=c(1,start), count=c(398, numrows))
n <- index.gdsn(gfile, "snp.position")
snppos <- read.gdsn(n) #, start=c(1,start), count=c(398, numrows))

n <- index.gdsn(gfile, "snp.chromosome")
snpchrome <- read.gdsn(n) #, start=c(1,start), count=c(398, numrows))

n <- index.gdsn(gfile, "snp.allele")
snpallele <- read.gdsn(n) #, start=c(1,start), count=c(398, numrows))

snp_info <- data.frame(snpid, snppos, snpchrome, snpallele)
write.table(snp_info, file='big_snp_data.tsv', quote=FALSE, sep='\t', col.names = NA)

n <- index.gdsn(gfile, "genotype")
​
sz <- 10
#sz <- 10000
​
nloop <- 10
#nloop <- 3549
​
bigsize <- 35481497
​
for ( i in 1:nloop ) {
    start <- (i-1)*sz + 1
    numrows <- min( sz, bigsize-start+1 ) 
    snp.matrix <- read.gdsn(n, start=c(1,start), count=c(398, numrows))
    filename<-paste0("snp_matrix_",start,"_",numrows,".tsv");
    print(filename)
    write.table(snp.matrix, file=filename, quote=FALSE, sep='\t', col.names = NA)
}

closefn.gds(gfile)
