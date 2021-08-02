library(fishpond)

library('tximport')
library(DESeq2)
library(limma)
library(Glimma)
library(edgeR)
library(tximeta)
library(GenomicFeatures)
library(DRIMSeq)
library(stageR)
library(tidyverse)

# To import salmon results into tximeta, we need the following files: quant.sf, cmd_info.json, lib_format_counts.json,
#  aux_info folder having meta file, libParams folder having flenDist file

dir <- system.file("extdata", package="macrophage")
list.files(dir)
library(readr)
library(dplyr)
coldata <- read_csv(file.path(dir,"coldata.csv"))
lvls <- c("naive","IFNg","SL1344","IFNg_SL1344")
coldata <- coldata %>% 
  dplyr::select(names, id=sample_id, line=line_id, 
                condition=condition_name) %>%
  mutate(line=factor(line),
         condition=factor(condition, levels=lvls),
         files=file.path(dir, "quants", names, "quant.sf.gz"))
coldata <- coldata %>% filter(condition %in% c("naive","IFNg"))
coldata$condition <- droplevels(coldata$condition)
library(tximeta)
suppressPackageStartupMessages(library(SummarizedExperiment))
y <- tximeta(coldata)


# create SummarizedExperiment object

se <- tximeta(coldata, dropInfReps=TRUE)

assayNames(se)

rowRanges(se)

# Because the tximeta function has identified the correct reference transcripts that were used for 
# quantification, and stored this as metadata in the SummarizedExperiment object, we can obtain additional 
# information, such as the correspondence of transcripts to genes. we can then perform a summarization task, 
# without having to look up this information manually.

gse <- summarizeToGene(se)

# We can easily subset the genes (rows) based on overlaps with a given genomic range using standard square bracket indexing

rowRanges(gse)

x <- GRanges("chr1", IRanges(10e6, 11e6))
gse[gse %over% x, ]

gse[, gse$line == "OCT4"]

#Finally, I demonstrate that it is easy to add alternative identifiers, making use of the organism 
# annotation packages in Bioconductor. For example, to add gene symbols, I can use tximeta’s addIds 
# function, and specify the SYMBOL column.

library(org.Mm.eg.db)
gse <- addIds(gse, column="SYMBOL")

# To check all the available columns for an organism package, use the columns function:

columns(org.Mm.eg.db)

# DGE analysis

# Count observation

cs <- colSums(assay(gse, "counts"))
hist(cs/1e6, col="grey", border="white",
     main="", xlab="column sums (per million)")
cts <- assay(gse, "counts")[,c(1,4)]
idx <- rowSums(cts >= 5) == 2
cts <- cts[idx,]
p <- sweep(cts, 2, colSums(cts), "/")

library(rafalib)
maplot(log10(p[,1]), log10(p[,2]), n=nrow(p), 
       cex=.3, col=rgb(0,0,0,.2),
       xlab="log10 proportion (geometric mean)", 
       ylab="log10 fold change")
abline(h=0, col=rgb(0,0,1,.5))

mean.log.p <- rowMeans(log10(p))
hist(mean.log.p,  col="grey", border="white",
     main="", xlab="log10 proportion (geometric mean)")


library(rafalib)
maplot(log10(p[,1]), log10(p[,2]), n=nrow(p),
       xlim=c(-6.5, -3.5), ylim=c(-2,2),
       cex=.3, col=rgb(0,0,0,.2),
       xlab="log10 proportion (geometric mean)", 
       ylab="log10 fold change")
abline(h=0, col=rgb(0,0,1,.5))

library(DESeq2)
dds <- DESeqDataSet(gse, design=~line + line:condition)
keep <- rowSums(counts(dds) >= 10) >= 3
table(keep)
dds <- dds[keep,]
dds <- DESeq(dds)
plotDispEsts(dds, ylim=c(1e-3, .5), xlim=c(5,1e5))
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("line","condition"))
resultsNames(dds)
res <- results(dds, name="lineOCT4.conditiontrt")
summary(res)
head(res, 3)

# We  can order adjusted p-values from small to large (but then remember that this object is no long 
# aligned with rows of dds for example).

res.ord <- res[order(res$padj),]

# I will show how to examine the differences across all genes in an MA-plot. But first, I will compute a new 
# estimate of fold change. The estimate in the results table above is the MLE, or maximum likelihood estimate,
# which is highly variable for low count genes (as we saw in the simulated example). Here I compute a Bayesian
# estimate for the fold change using methods in the apeglm package (27). The apeglm functions are wrapped up in 
# a DESeq2 function called lfcShrink, which produces a table similar to results but with shrunken LFC instead of 
# MLE.

library(apeglm)
lfc <- lfcShrink(dds, coef="lineOCT4.conditiontrt", type="apeglm")

# DTE analysis
y <- tximeta(coldata)
y <- scaleInfReps(y, quiet=TRUE)
y <- labelKeep(y)
y <- y[mcols(y)$keep,]

# Because the method makes use of permutations, it is required to set a seed for computational 
# reproducibility. We specify to test across the condition variable, while controlling for a pairing 
# variable line. The line variable indicates which donor the cell line came from.

set.seed(1)
y <- swish(y, x="condition", pair="line", quiet=TRUE)

# After running swish, all of the results are stored in the metadata columns (mcols) of the object y. 
# We look to see how many transcripts have a small q-value (analogous to an adjusted p-value, this 
# should provide a set with a nominal FDR control).

names(mcols(y))

table(mcols(y)$qvalue < .05)

# One important aspect in testing across many features, in particular where the uncertainty level is so 
# heterogeneous, is to consider if the p-value distribution is roughly uniform, with the exception of the 
# rejected tests. Here Swish provides a roughly uniform distribution, with a spike on the left side representing the rejections of the null hypothesis.

hist(mcols(y)$pvalue, col="grey",
     main="", xlab="p-values")

# We  can make an MA-plot, with the differential transcripts highlighted in blue (here at 5% FDR).

plotMASwish(y, alpha=.05)

# We can also examine individual transcripts with evidence of differential expression. As each sample is 
# represented by a distribution of possible estimated counts from Salmon, Swish uses boxplots to represent 
# the differences in expression across samples:

idx <- with(mcols(y), which(pvalue < .05 & log2FC > 4))
plotInfReps(y, idx[1], x="condition", cov="line", xaxis=FALSE)
