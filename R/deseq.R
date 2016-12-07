# A deseq gene expression difference analyze pipeline
# reference: https://dwheelerau.com/2013/04/15/how-to-use-deseq-to-analyse-rnaseq-data/

library(DESeq)
# load count table
countsTable <- read.delim('../merged.csv',header=TRUE)
rownames(countsTable) <- countsTable$gene
countsTable <- countsTable[,-1]
# conditions
conds <- factor( c( "t1", "t2", "t3", "t4", "t5", "c1", "c5"  ) )
cds <- newCountDataSet(countsTable, conds)
cds <- estimateSizeFactors( cds )
# set size factor manually
#sizefactor <- c()
#sizeFactors( cds ) <- sizefactor

cds <- estimateDispersions(cds)
# single sample
#cds <- estimateDispersions(cds, method='blind', sharingMode="fit-only") 

res_c1_t1 = nbinomTest(cds, 'c1', 't1')
res_c1_t5 = nbinomTest(cds, 'c1', 't5')
res_c1_c5 = nbinomTest(cds, 'c1', 'c5')

# write result to csv file
write.csv(res_c1_t1, file='./c1_t1')
write.csv(res_c1_t5, file='./c1_t5')
write.csv(res_c1_c5, file='./c1_c5')
