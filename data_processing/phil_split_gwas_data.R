library(gdsfmt)
​
gfile = openfn.gds("../_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds")
​
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
    filename<-paste0("../_ignore/snp_matrix_",start,"_",numrows,".csv");
    print(filename)
    write.csv( snp.matrix, filename )
}

