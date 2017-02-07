
# laungh edgeR pakcage
library(edgeR)
# modify the scirpt from the sequencing pipeline, majorly for the data import. export list, key and count file from the pipeline
args = c("ahr-yating","ahr-yating.csv")
# Get counts file from analysis/fname/fname.T.csv
bname = basename(args[1])
fname = paste(bname,"-count.T.csv",sep='')
dat = read.csv(file.path("analysis",bname,fname), header = TRUE, row.names=1)
# Filter low count reads
keep = rowSums(cpm(dat) > 3) >= 3
counts = dat[keep,]
## Read in key file
key = read.csv(args[2], header=TRUE, row.names=1)
# create groups or model comparision matrix with sample key
sel = grepl("MT-.*", rownames(counts)) + grepl("ERCC-.*", rownames(counts)) + grepl("mt-.*", rownames(counts))
counts = counts[!sel,]
factors = key[order(rownames(key)), c("type")]
groups = factors
counts = counts[,order(names(counts))]
#run pairwise analysis
y = DGEList(counts=counts, group=factors)
y = calcNormFactors(y)
y = estimateCommonDisp(y)
y = estimateTagwiseDisp(y)

## get normalized counts for each group for outputting to summary spreadsheet
scaled.counts = data.frame(mapply(`*`, counts, y$samples$lib.size *
                                    y$samples$norm.factors/mean(y$samples$lib.size)))
rownames(scaled.counts) = rownames(counts)
dfs = split.data.frame(t(scaled.counts), groups)
dfss = sapply(dfs, colMeans)


# write results
#pair-wise comparison group1 and group2
lrt = exactTest(y, pair = c(1,2))
ot1 = topTags(lrt,n=nrow(counts),sort.by="PValue")$table
ot1 = merge(ot1, dfss, by=0)
ot1 = ot1[order(ot1$FDR),] # Sort by ascending FDR
write.csv(ot1,file="Indoel_DMSO.csv",row.names=FALSE)

#pair-wise comparison group1 and group3
lrt = exactTest(y, pair = c(1,3))
ot1 = topTags(lrt,n=nrow(counts),sort.by="PValue")$table
ot1 = merge(ot1, dfss, by=0)
ot1 = ot1[order(ot1$FDR),] # Sort by ascending FDR
write.csv(ot1,file="TCDD_DMSO.csv",row.names=FALSE)


