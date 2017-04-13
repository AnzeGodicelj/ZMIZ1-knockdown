
setwd("/Users/AnzeGodicelj/Dropbox/T47D siZMIZ1")
load("/Users/AnzeGodicelj/Dropbox/rnaseq/named_newcounts.rda")
cellLine <- gsub(' .+', '', colnames(named_newcounts))
treatment <- gsub('.+ ', '', colnames(named_newcounts))
passage <- gsub("^.*? ","", colnames(named_newcounts))
passage <- gsub(' .*?$','', passage)
experimental_data <- cbind(cellLine, treatment, passage)

# Extracting counts only for MCF7 cell line from the matric named_newcounts
T47D_named_newcounts <- as.matrix(named_newcounts[ , grep("T47D", colnames(named_newcounts))])
save(T47D_named_newcounts, file = "T47D_named_newcounts.rda")

# Extracting the experimental data for MCF7 cell line from the matrix experimental_data
T47D_experimental_data <- as.matrix(experimental_data[grep("T47D", experimental_data), ])


# Load counts and experiment data into a DESeq object named MCF7_deSeqData
rownames(T47D_experimental_data)<-colnames(T47D_named_newcounts)
View(T47D_experimental_data)
T47D_deSeqData <- DESeqDataSetFromMatrix(# counts are in the matrix T47D_names_newcounts.
  countData = T47D_named_newcounts, 
  # experimental data for each column of T47D_named_newcounts 
  # are in matrix T47D_experimental_data.
  colData = T47D_experimental_data, 
  # when analysing the data, take into account treatment and passage od the sample
  # (stored in T47D_experimental_data)
  design = ~treatment + passage)
T47D_deSeqData

# Let's filter out non- or very low expressed genes which are be detected in only one or zero samples,
# after all , we dont want to model our differential expression model on the genes that do not change at all
# (removes 6678 non-expressed genes) and stroes it in a vector named expressed_genes.
T47D_expressed_genes <- T47D_deSeqData[rowSums(counts(T47D_deSeqData))>1]

summary(rowSums(counts(T47D_deSeqData)))
nrow(T47D_deSeqData)

summary(rowSums(counts(T47D_expressed_genes)))
nrow(T47D_expressed_genes)

# Now it is time to normalise our counts accross the samples by estimating size factors. 
# This function corrects the read counts for a gene in each sample 
# based on the distribution of the number of reads in that sample.
T47D_normalized_genes <- estimateSizeFactors(T47D_expressed_genes)
counts(T47D_normalized_genes, normalized = TRUE) ["ZMIZ1", ] # See that siZMIZ1 samples show lower normalised reads values than siCTRL or RNAiMAX samples.

#Let's save MCF7 normalised counts in a matrix named MCF7_normalised_matrix.
T47D_normalized_matrix <- counts(T47D_normalized_genes, normalized = TRUE)
View(T47D_normalized_matrix)
colData(T47D_normalized_genes)
save(T47D_normalized_matrix, file = "T47D normalised matrix.rda")
write.csv(T47D_normalized_matrix, file = "T47D normalised matrix.csv")

# Estimating dispersions to add to the data on which will be used in the negative binomial test.
T47D_gene_dispersion <- estimateDispersions(T47D_normalized_genes)
mcols(T47D_gene_dispersion )
png("Dispersion estimates of the T47D dataset.png",w=8000,h=5000,p=100)
plotDispEsts(T47D_gene_dispersion, 
             xlab = "Mean of normalised counts", 
             ylab = "Dispersion",
             ymin = 1e-5,
             cex.lab=1.5,
             cex.axis = 1,
             main = "Dispersion estimates of the T47D dataset",
             cex.main=2,
             legend = FALSE)
legend("bottomright", 
       pch = 19, 
       legend = c("Gene-wise dispersion estimates", "Fitted curve", "Final dispersion estimates"), 
       col = c("black", "red", "dodgerblue"), 
       pt.cex=1, 
       cex=1)
grid()
dev.off()


# SHOWTIME! Time to do a statistical test on normalised counts and a dispersion estimate per each gene.
# Genes that have a similar dispersion and expression will be used to create a distribution of the gene counts.
# Ultimatelly these distributions will be compared in the negative binomial test. 

T47D_MA_data <- nbinomWaldTest(T47D_gene_dispersion)
T47D_MA_results <- results(T47D_MA_data, contrast = c('treatment', 'siZMIZ1', 'siCTRL')) ##Choosing the right contrast
MA_results
View(T47D_MA_results)
summary(T47D_MA_results)

save(T47D_MA_results, file = "T47D MA results.rda")

# Plotting the results of negative binomial test.
par(las=1)
png("DESeq2 T47D MA Plot.png",w=8000,h=5000,p=100)
plotMA(nbinomWaldTest(T47D_gene_dispersion), 
       main = "DESeq2 T47D MA Plot", 
       xlab = "Mean of normalised counts",
       cex.main=2,
       alpha = 0.01)

# Summing up genes that are significantly differentially expressed (p<0.01).
number_of_significant_genes <- sum(T47D_MA_results$padj<0.01, na.rm = TRUE)
mtext(paste0("Significant genes: ", number_of_significant_genes))
dev.off()



