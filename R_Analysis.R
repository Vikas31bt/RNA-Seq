library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)
# 8. load the sample information
pheno_data <- read.csv("chrX_data/geuvadis_phenodata.csv")

# 9. Read in the expression data that were calculated by StringTie and compute the FPKM.
bg_chrX <- ballgown(dataDir = "ballgown",samplePattern = "ERR",pData = pheno_data)
class(bg_chrX)
bg_chrX

# Export Gene expression - FPMK
gene_expression <- gexpr(bg_chrX)
head(gexpr(bg_chrX), 1)


# Export Transcript expression - FPMK
transcript_expression <- texpr(bg_chrX)
head(texpr(bg_chrX), 1)

transcript_expression = texpr(bg_chrX, 'FPKM')
transcript_cov = texpr(bg_chrX, 'cov')
whole_tx_table = texpr(bg_chrX, 'all')
exon_mcov = eexpr(bg_chrX, 'mcov')
junction_rcount = iexpr(bg_chrX)
whole_intron_table = iexpr(bg_chrX, 'all')
gene_expression = gexpr(bg_chrX)

# 10.Filter to remove low-abundance genes.The genes with a low read counts have to be discarded. It is also possible to apply a low variance filter. Here, we remove all the transcripts having a variance lower than 1.
bg_chrX_filt <- subset(bg_chrX, "rowVars(texpr(bg_chrX)) >1", genomesubset=TRUE)
bg_chrX_filt
head(pData(bg_chrX_filt), 3)

# 11.Identify transcripts that show statistically significant differences between groups.
results_transcripts <- stattest(bg_chrX_filt,feature="transcript",covariate="sex",adjustvars = c("population"),getFC=TRUE, meas="FPKM")
dim(results_transcripts)
table(results_transcripts$qval < 0.05)

# 12.Identify genes that show statistically significant differences between groups.
results_genes <- stattest(bg_chrX_filt,feature="gene",covariate="sex",adjustvars = c("population"),getFC=TRUE, meas="FPKM")
class(results_genes)
dim(results_genes)
table(results_genes$qval<0.05)

# 13.Add gene names and gene IDs to the results_transcripts data frame:
results_transcripts <- data.frame(geneNames = geneNames(bg_chrX_filt), geneIDs = geneIDs(bg_chrX_filt),results_transcripts)

# 14.Sort the results from the smallest P value to the largest:
results_transcripts = arrange(results_transcripts,pval)
results_genes = arrange(results_genes,pval)
head(results_transcripts)

# 15.Write the results to a csv file that can be shared and distributed:
write.csv(results_transcripts, "chrX_transcript_results.csv", row.names=FALSE)
write.csv(results_genes, "chrX_gene_results.csv", row.names=FALSE)
write.csv(gene_expression, "chrX_gene_expression_results.csv", row.names=FALSE)
write.csv(transcript_expression, "chrX_transcript_expression_results.csv", row.names=FALSE)

# 16.Identify transcripts and genes with a q value <0.05:
subset(results_transcripts,results_transcripts$qval<0.05)
subset(results_genes,results_genes$qval<0.05)

library(gplots)
library(RColorBrewer)
library(ggplot2)
library(cowplot)

# 17.Make the plots
tropical= c('darkorange', 'dodgerblue','hotpink', 'limegreen', 'yellow')
palette(tropical)

# 18. Show the distribution of gene abundances (measured as FPKM values) across samples, colored by sex
fpkm = texpr(bg_chrX,meas="FPKM")
fpkm = log2(fpkm+1)
boxplot(fpkm,col=as.numeric(pheno_data$sex),las=2,ylab='log2(FPKM+1)')

# 19. Make plots of individual transcripts across samples.
ballgown::transcriptNames(bg_chrX)[12]
ballgown::geneNames(bg_chrX)[12]
plot(fpkm[12,] ~ pheno_data$sex, border=c(1,2), main=paste(ballgown::geneNames(bg_chrX)[12],' : ',ballgown::transcriptNames(bg_chrX)[12]),pch=19, xlab="Sex",
ylab='log2(FPKM+1)')

# 20.Plot the structure and expression levels in a sample of all transcripts that share the same gene locus.
plotTranscripts(ballgown::geneIDs(bg_chrX)[1729], bg_chrX, main=c('Gene XIST in sample ERR188234'), sample=c('ERR188234'))
plotTranscripts(ballgown::geneIDs(bg_chrX)[ballgown::geneNames(bg_chrX) == "XIST"], bg_chrX, main=c('Gene XIST in sample ERR188234'), sample=c('ERR188234'))

# 21.average expression levels for all transcripts of a gene within different groups
plotMeans(ballgown::geneIDs(bg_chrX)[203], bg_chrX_filt, groupvar="sex", legend=FALSE)

# 22.MA plot 
results_transcripts$mean <- rowMeans(texpr(bg_chrX_filt))
ggplot(results_transcripts, aes(log2(mean), log2(fc), colour = qval<0.05)) + scale_color_manual(values=c("#999999", "#FF0000")) + geom_point() +
geom_hline(yintercept=0)

# 23. The log fold change (logFC) can be computed and Visualisation- Volcano plot
results_transcripts$logFC <- log2(results_transcripts$fc)
results_genes$logFC <- log2(results_genes$fc)
volcano_plot <- function(data){
  logfc.threshold <- 1
  with(data, plot(logFC, -log10(qval), pch=20, main="Volcano plot"))
  with(subset(data, qval<.05 ), points(logFC, -log10(qval), pch=20, col="red"))
  with(subset(data, abs(logFC)>logfc.threshold), points(logFC, -log10(qval), pch=20, col="orange"))
  with(subset(data, qval<.05 & abs(logFC)>logfc.threshold), points(logFC, -log10(qval), pch=20, col="green"))
}
volcano_plot(results_genes)

# 24. Heatmap
my_palette <- colorRampPalette(c("green", "yellow", "red"))(n = 299)
mat_data <- data.matrix(fpkm)
heatmap.2(mat_data, main = "FPKM", notecol="black", trace="none", margins =c(12,9), col=my_palette, dendrogram="both")
         
          

results_transcripts%>%
  filter(feature == 'transcript') %>%
  ggplot(., aes(x = geneNames, y = logFC, fill = id)) + geom_col()
 
