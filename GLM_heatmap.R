
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


counts = counts[,order(names(counts))]

# run GLM analysis 
y = DGEList(counts=counts)
y = calcNormFactors(y)
design = model.matrix(~0+factors, data=y$samples)
colnames(design) <- levels(factors)
y = estimateGLMCommonDisp(y, design)
y = estimateGLMTrendedDisp(y, design)
y = estimateGLMTagwiseDisp(y, design)
fit = glmFit(y, design)


## get normalized counts for each group for outputting to summary spreadsheet
scaled.counts = data.frame(mapply(`*`, counts, y$samples$lib.size *
                                    y$samples$norm.factors/mean(y$samples$lib.size)))
rownames(scaled.counts) = rownames(counts)
dfs = split.data.frame(t(scaled.counts), groups)
dfss = sapply(dfs, colMeans)


# write results
Indole_dmso <- glmLRT(fit, contrast=c(-1,1,0))
TCDD_dmso <- glmLRT(fit, contrast=c(-1,0,1))

##########################plotSmear##############################################
topTags(TCDD_dmso,n=10)

dt_significant
vctr_names_sig

dt_significant <- decideTestsDGE( TCDD_dmso, adjust.method="BH", p.value=0.05)
vctr_names_sig <- rownames(y)[ as.logical( dt_significant )]
plotSmear( TCDD_dmso, de.tags = vctr_names_sig )
abline( h = c( -2, 2 ), col = "blue")

##########################heatmap##############################################
vctr_names_top <- rownames( topTags(TCDD_dmso, n = 30 ) )
vctr_names_top <- c( vctr_names_top, rownames( topTags( TCDD_dmso, n = 30 ) ) )
vctr_names_top <- unique( c( vctr_names_top, rownames( topTags( TCDD_dmso, n = 30 ) ) ) )

vctr_sig <- as.logical( decideTestsDGE( TCDD_dmso, adjust.method="BH", p.value=0.000005) )
vctr_sig <- vctr_sig | as.logical( decideTestsDGE( Indole_dmso, adjust.method="BH", p.value=0.000005) )
vctr_sig <- vctr_sig | as.logical( decideTestsDGE( Indole_dmso, adjust.method="BH", p.value=0.000005) )
vctr_names_hcl <- rownames(y)[ vctr_sig ]
length( vctr_names_hcl )
head( vctr_names_hcl )

mtrx_significant <- y$counts[ vctr_names_top, ]
vctr_colors = as.factor( c( "black", "red", "green" ))
vctr_colors
vctr_sample_colors <- as.character( vctr_colors[ as.numeric(factors) ] )
vctr_sample_colors
heatmap.2( log2( mtrx_significant + 1 ), ColSideColors=vctr_sample_colors, key=TRUE, trace="none", col=heat.colors(200), scale="row" )
