library(DESeq2)
library(data.table)
library(limma)
library(calibrate)
setwd("/media/hp/04c65089-71ff-4b33-9a30-c21b8c77eda2/li/uci")
DF<-read.table("uci.txt", head=TRUE)
b<-aggregate(DF[, -c(1:2)], by=list(DF$EntrezID, DF$Name), mean)
write.table(b, file="uci_averageif.txt")

setwd("/media/hp/04c65089-71ff-4b33-9a30-c21b8c77eda2/li/uci/deg")
#delta11_vs_mock
cts<-read.table("delta11_vs_mock.txt",head=TRUE) 
#setup DESeqDataSetFromMatrix
cts_round<-round(cts)
countData<-cts_round
condition <- factor(c(rep("delta11",11),rep("mock",9)), levels = c("delta11","mock"))
coldata<-data.frame(row.names=colnames(countData),condition)
condition
coldata
dds<-DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=~condition)
head(dds)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10#keep sum counts >=10
dds <- dds[keep,]
dim(dds)
dds$condition <- factor(dds$condition, levels = c("delta11","mock"))
dds$condition <- relevel(dds$condition, ref = "mock")
dds$condition <- droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
dds_norm<-counts(dds,normalized=TRUE)
head(dds_norm)
write.csv(dds_norm,"delta11_vs_mock_norm.csv")
dds2<-DESeq(dds)
resultsNames(dds2)
res <- results(dds2)
res <- results(dds2, contrast=c("condition","delta11","mock"))
resultsNames(dds2)
write.csv(res,"delta11_vs_mock_deg.csv")

#delta11_vs_wt
cts<-read.table("delta11_vs_wt.txt",head=TRUE)
#setup DESeqDataSetFromMatrix
cts_round<-round(cts)
countData<-cts_round
condition <- factor(c(rep("delta11",11),rep("wt",11)), levels = c("delta11","wt"))
coldata<-data.frame(row.names=colnames(countData),condition)
condition
coldata
dds<-DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=~condition)
head(dds)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10#keep sum counts >=10
dds <- dds[keep,]
dim(dds)
dds$condition <- factor(dds$condition, levels = c("delta11","wt"))
dds$condition <- relevel(dds$condition, ref = "wt")
dds$condition <- droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
dds_norm<-counts(dds,normalized=TRUE)
head(dds_norm)
write.csv(dds_norm,"delta11_vs_wt_norm.csv")

dds2<-DESeq(dds)
resultsNames(dds2)
res <- results(dds2)
res <- results(dds2, contrast=c("condition","delta11","wt"))
resultsNames(dds2)
write.csv(res,"delta11_vs_wt_deg.csv")

#mock_vs_wt
cts<-read.table("mock_vs_wt.txt",head=TRUE)
#setup DESeqDataSetFromMatrix
cts_round<-round(cts)
countData<-cts_round
condition <- factor(c(rep("mock",9),rep("wt",11)), levels = c("mock","wt"))
coldata<-data.frame(row.names=colnames(countData),condition)
condition
coldata
dds<-DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=~condition)
head(dds)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10#keep sum counts >=10
dds <- dds[keep,]
dim(dds)
dds$condition <- factor(dds$condition, levels = c("mock","wt"))
dds$condition <- relevel(dds$condition, ref = "wt")
dds$condition <- droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
dds_norm<-counts(dds,normalized=TRUE)
head(dds_norm)
write.csv(dds_norm,"mock_vs_wt_norm.csv")

dds2<-DESeq(dds)
resultsNames(dds2)
res <- results(dds2)
res
resultsNames(dds2)
write.csv(res,"mock_vs_wt_deg.csv")

#volcano plot delta11_vs_mock_DEG
res <- read.csv("delta11_vs_mock_deg.csv", header=TRUE)
head(res)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-5,7)))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
with(subset(res, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, pvalue<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)
#dev.off()

#volcano plot delta11_vs_wt_DEG
res <- read.csv("delta11_vs_wt_deg.csv", header=TRUE)
head(res)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-5,6)))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
with(subset(res, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, pvalue<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)
#dev.off()

#volcano plot mock_vs_wt_DEG
res <- read.csv("mock_vs_wt_deg.csv", header=TRUE)
head(res)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-6,7)))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
with(subset(res, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, pvalue<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)
#dev.off()


library(fgsea)
library(tidyverse)
library(data.table)
library(ggplot2)
setwd("/media/hp/04c65089-71ff-4b33-9a30-c21b8c77eda2/li/uci/deg")

#delta11_vs_mock
res<-read.table("delta11_vs_mock_fgsea.txt")
ranks <- deframe(res)
pathways.kegg<- gmtPathways("c2.cp.kegg.v7.2.symbols.gmt")
pathways<-gmtPathways("c2.cp.v7.2.symbols.gtm")
fgseaRes_kegg<- fgsea(pathways = pathways.kegg, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple = 10000)
fgseaRes<- fgsea(pathways = pathways, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple = 10000)
fwrite(fgseaRes_kegg, file="delta11_vs_mock_kegg.txt", sep="\t", sep2=c("", " ", ""))
fwrite(fgseaRes, file="delta11_vs_mock_all.txt", sep="\t", sep2=c("", " ", ""))

#delta11_vs_wt
res<-read.table("delta11_vs_wt_fgsea.txt")
ranks <- deframe(res)
pathways.kegg<- gmtPathways("c2.cp.kegg.v7.2.symbols.gmt")
pathways<-gmtPathways("c2.cp.v7.2.symbols.gtm")
fgseaRes_kegg<- fgsea(pathways = pathways.kegg, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple=10000)
fgseaRes<- fgsea(pathways = pathways, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple=10000)
fwrite(fgseaRes_kegg, file="delta11_vs_wt_kegg.txt", sep="\t", sep2=c("", " ", ""))
fwrite(fgseaRes, file="delta11_vs_wt_all.txt", sep="\t", sep2=c("", " ", ""))

plotEnrichment(pathways.kegg[["KEGG_HUNTINGTONS_DISEASE"]],
               ranks) + labs(title="HUNTINGTONS_DISEASE")

plotEnrichment(pathways.kegg[["KEGG_OXIDATIVE_PHOSPHORYLATION"]],
               ranks) + labs(title="OXIDATIVE_PHOSPHORYLATION")

#mock_vs_wt
res<-read.table("mock_vs_wt_fgsea.txt")
ranks <- deframe(res)
pathways.kegg<- gmtPathways("c2.cp.kegg.v7.2.symbols.gmt")
pathways<-gmtPathways("c2.cp.v7.2.symbols.gtm")
fgseaRes_kegg<- fgsea(pathways = pathways.kegg, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple=10000)
fgseaRes<- fgsea(pathways = pathways, stats  = ranks,  minSize = 1,maxSize  = Inf, nPermSimple=10000)
fwrite(fgseaRes_kegg, file="mock_vs_wt_kegg.txt", sep="\t", sep2=c("", " ", ""))
fwrite(fgseaRes, file="mock_vs_wt_all.txt", sep="\t", sep2=c("", " ", ""))


