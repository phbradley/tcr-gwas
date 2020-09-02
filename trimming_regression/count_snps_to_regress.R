
library("gdsfmt")
library(data.table)
setDTthreads(threads = 1, restore_after_fork=FALSE)

gfile = openfn.gds("../_ignore/snp_data/HSCT_comb_geno_combined_v03_tcr.gds")
n <- index.gdsn(gfile, "genotype")
#sz <- 10
#sz <- 10000
#sz <- 50000
sz <- 5000
#nloop <- 10
#nloop <- 3549
#nloop <- 71
nloop <- 7098
bigsize <- 35481497
count_dataframe = data.table()
for ( i in 1:nloop ) {
    count = 0
    start <- (i-1)*sz + 1
    numrows <- min( sz, bigsize-start+1 ) 
    snp.matrix <- read.gdsn(n, start=c(1,start), count=c(398, numrows))
    for (k in 1:ncol(snp.matrix)){
        if (length(unique(snp.matrix[,k])[unique(snp.matrix[,k]) != 3]) > 1){
            count = count + 1
        } else {
            count = count
        }
    }
    count_dataframe = rbind(count_dataframe, data.table(set = i, count = count))
    print(paste0('finished loop ', i, ' with count of ', count))
}
print(paste0('finished counting. COUNT IS ', count))




    