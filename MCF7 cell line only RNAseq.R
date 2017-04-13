#### RNAseq
library(Rsubread)
library(DESeq2)
library(devtools)
library(Biobase)
library(preprocessCore)
library("gplots")

setwd("/Users/godice01/Documents/rnaseq")
#List BAM files

FilesToCount <- dir("/Users/godice01/Documents/rnaseq/SLX-13365", pattern = '.bam$', full.names = TRUE)
FilesToCount 
# pattern = '.bam$' means that the string (file name) ends with .bam and not by .bami etc.
# pattern = '.bam' means what?

# Counting reads 
filename <- 'rawcounts.rda' #filename is an object that is in file rawcount.rda
if(!file.exists(filename)) # if file rawdata.rda does not exist run the next bit of code 
{
  rawcounts <- featureCounts(FilesToCount, annot.inbuilt = "hg38", ignoreDup = FALSE, countMultiMappingReads = FALSE, minOverlap = 1)
  rawcounts$counts 
  save(rawcounts, file = filename)
  
}else{
  load(filename)
}
#Merge replicates three by three
newcounts<-matrix(0,nrow=nrow(rawcounts$counts),ncol=ncol(rawcounts$counts)/3)
rownames(newcounts)<-rownames(rawcounts$counts)
samplenames<-gsub("X.Users.godice01.Documents.rnaseq.SLX.13365.SLX.13365.","",colnames(rawcounts$counts))
samplenames<-gsub(".HHGWGBBXX.s_[0-9].tophat.homo_sapiens.bam","",samplenames)
samplenames<-unique(samplenames)
colnames(newcounts)<-samplenames
newcounts
# And finally we merge all three into a single sample
i <- 1
ii <- 1
while(i<ncol(rawcounts$counts)){
  replicates<-rawcounts$counts[,i:(i+2)]
  merged<-apply(replicates,1,sum)
  newcounts[,ii]<-merged
  i<-i+3
  ii<-ii+1
}
newcounts

# Great! Now we can do stuff like converting entrez ids to gene symbols 

library("org.Hs.eg.db")
mapped_genes <- mappedkeys(org.Hs.egSYMBOL)
entrez2symbol <- as.list(org.Hs.egSYMBOL[mapped_genes])
entrez2symbol[["7157"]]

# Convert all gene names in the newcounts
genesymbols<-unlist(entrez2symbol[rownames(newcounts)])
geneswithsymbols<-names(genesymbols)
newcounts<-newcounts[geneswithsymbols,]
rownames(newcounts)<-genesymbols

######Rename columns with meaningful names
##Importing the SLX_ID file (from sample submission form)
SLX_ID <- read.csv("/Users/godice01/Documents/rnaseq/RNAseq data analysis/SLX_ID.csv", as.is = TRUE)
SLX_ID
SLX_ID[,2] <- gsub('-', '_', SLX_ID$SLX_ID)
View(SLX_ID)

#creating the transforming vector
experiment_names <- SLX_ID$Sample_ID
names(experiment_names) <- SLX_ID$SLX_ID
experiment_names

#Changing the names of the columens from SLX_IDs to experiment names

filename2 <- 'named_newcounts.rda'

named_newcounts <- newcounts

colnames(named_newcounts)<-as.vector(experiment_names[colnames(newcounts)])
View(named_newcounts)

save(named_newcounts, file = filename2)

colnames(newcounts)
View(newcounts)
colnames(named_newcounts)
View(named_newcounts)


#Creating DESeqDataSet
load('named_newcounts.rda')

#Delete names of rownames of tanle replicate_newcounts

names(rownames(named_newcounts)) <- NULL

#Load experiment data
#experiment_data <- read.csv('/Users/AnzeGodicelj/Library/Mobile Documents/com~apple~CloudDocs/Documents/rnaseq/experiment_data.csv')
#experimental_data <- experiment_data[,c(1,3)]
#experimental_data

cellLine <- gsub(' .+', '', colnames(named_newcounts))
treatment <- gsub('.+ ', '', colnames(named_newcounts))
passage <- gsub("^.*? ","", colnames(named_newcounts))
passage <- gsub(' .*?$','', passage)
experimental_data <- cbind(cellLine, treatment, passage)

# Extracting counts only for MCF7 cell line from the matric named_newcounts
MCF7_named_newcounts <- as.matrix(named_newcounts[ , grep("MCF7", colnames(named_newcounts))])

# Extracting the experimental data for MCF7 cell line from the matrix experimental_data
MCF7_experimental_data <- as.matrix(experimental_data[grep("MCF7", experimental_data), ])


# Load counts and experiment data into a DESeq object named MCF7_deSeqData
rownames(MCF7_experimental_data)<-colnames(MCF7_named_newcounts)
View(MCF7_experimental_data)
MCF7_deSeqData <- DESeqDataSetFromMatrix(# counts are in the matrix MCF7_names_newcounts.
  countData = MCF7_named_newcounts, 
  # experimental data for each column of MCF7_named_newcounts 
  # are in matrix MCF7_experimental_data.
  colData = MCF7_experimental_data, 
  # when analysing the data, take into account treatment and passage od the sample
  # (stored in MCF7_experimental_data)
  design = ~treatment + passage)
MCF7_deSeqData

# Let's filter out non- or very low expressed genes which are be detected in only one or zero samples,
# after all , we dont want to model our differential expression model on the genes that do not change at all
# (removes 6678 non-expressed genes) and stroes it in a vector named expressed_genes.
MCF7_expressed_genes <- MCF7_deSeqData[rowSums(counts(MCF7_deSeqData))>1]

summary(rowSums(counts(MCF7_deSeqData)))
nrow(MCF7_deSeqData)

summary(rowSums(counts(MCF7_expressed_genes)))
nrow(MCF7_expressed_genes)

# Now it is time to normalise our counts accross the samples by estimating size factors. 
# This function corrects the read counts for a gene in each sample 
# based on the distribution of the number of reads in that sample.
MCF7_normalized_genes <- estimateSizeFactors(MCF7_expressed_genes)
counts(MCF7_normalized_genes, normalized = TRUE) ["ZMIZ1", ] # See that siZMIZ1 samples show lower normalised reads values than siCTRL or RNAiMAX samples.

#Let's save MCF7 normalised counts in a matrix named MCF7_normalised_matrix.
MCF7_normalized_matrix <- counts(MCF7_normalized_genes, normalized = TRUE)
View(MCF7_normalized_matrix)
colData(MCF7_normalized_genes)
write.csv(MCF7_normalized_matrix, file = "MCF7 normalised matrix.csv")

# Estimating dispersions to add to the data on which will be used in the negative binomial test.
MCF7_gene_dispersion <- estimateDispersions(MCF7_normalized_genes)
mcols(MCF7_gene_dispersion )
png("Dispersion estimates of the MCF7 dataset.png",w=8000,h=5000,p=100)
plotDispEsts(MCF7_gene_dispersion, 
             xlab = "Mean of normalised counts", 
             ylab = "Dispersion",
             ymin = 1e-5,
             cex.lab=1.5,
             cex.axis = 1,
             main = "Dispersion estimates of the MCF7 dataset",
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

MCF7_MA_data <- nbinomWaldTest(MCF7_gene_dispersion)
MCF7_MA_results <- results(MCF7_MA_data, contrast = c('treatment', 'siZMIZ1', 'siCTRL')) ##Choosing the right contrast
MA_results
View(MCF7_MA_results)
summary(MCF7_MA_results)

# Saving the results of the negative binomial test on the MCF7_gene_dispersion.
MCF7_MA_results_table <- as.data.frame(MCF7_MA_results) 
write.csv(MCF7_MA_results_table, file = "MCF7 MA results table.csv") 

# Plotting the results of negative binomial test.
par(las=1)
png("DESeq2 MCF7 MA Plot.png",w=8000,h=5000,p=100)
plotMA(MCF7_MA_data, 
      main = "DESeq2 MCF7 MA Plot", 
      alpha=0.01,
      xlab = "Mean of normalised counts", 
      ylim = c(-1,1),      
      cex.main=2)
# Summing up genes that are significantly differentially expressed (p<0.01).
number_of_significant_genes <- sum(MCF7_MA_results_table$padj<0.01, na.rm = TRUE)
mtext(paste0("Significant genes: ", number_of_significant_genes))
dev.off()

# Significant genes and heatmap
#extracting a list of significant genes
MCF7_significant_genes_padj <- MCF7_MA_results_table$padj[which(MCF7_MA_results_table$padj<0.01)]
MCF7_significant_genes_names <- rownames(MCF7_MA_results_table)[which(MCF7_MA_results_table$padj<0.01)]

MCF7_significant_genes <- as.data.frame(cbind(MCF7_significant_genes_names, MCF7_significant_genes_padj))
rownames(MCF7_significant_genes) <- rownames(MCF7_MA_results_table)[which(MCF7_MA_results_table$padj<0.01)]
MCF7_significant_genes$MCF7_significant_genes_names <- NULL

save(MCF7_significant_genes, file = "MCF7 sifnificant genes.rda")
write.csv(MCF7_significant_genes, file = "MCF7 significant genes.csv")
# Taken from above normalisation: 

# Applying the list of significant genes from MA analysis on the matrix of normalised counts
# (normalized_matrix <- counts(normalized_genes, normalized = TRUE) - creates a  matrix of 
# all 6758 differentially expressed genes and their normalised counts in each experiment.

significant_counts<-normalized_matrix[significant_genes,]
colors<-colorpanel(1000,"navy","white","red3")
png("all_genes_heatmap.png",w=8000,h=8000,p=100)
heatmap.2(significant_counts,trace="none",scale="row",col=colors, margins = c(13,8), labRow = FALSE, 
          ylab = "All differentially expressed genes (N=6758, p<0.01)", dendrogram="col")
dev.off()

# Select top variance genes
vars<-apply(normalized_matrix,1,var)
topvars<-names(sort(vars,decreasing=TRUE))[1:1000]
submatrix<-normalized_matrix[topvars,]
# You can now calculate PCA with the top varying genes, this saves time and memory, and it gives more weight
# to genes that actually change, rather than giving equal weight also to noise



test<-rbind(c(1,2,3),c(3,0,5),c(3,2,9),c(-1,1,1))
apply(test,1,sum)
apply(test,2,sum)


###PCA plots for MCF7 cell line

# Proportion of variance explained barplot
MCF7_pca <- prcomp(MCF7_normalized_matrix, scale. = FALSE)

# Assess how much variance is explained by each of the top n component
MCF7_totvar<-sum(MCF7_pca$sdev^2)
MCF7_pcavar<-((MCF7_pca$sdev^2)/MCF7_totvar)*100
n <- 10
png("MCF7_pca_varianceExplained.png",width=10000,height=8000,pointsize=60, res = 300)
bp<-barplot(
  MCF7_pcavar[1:n],
  ylab="% Variance Explained",
  xlab=paste0("Principal Component (1 to ",n,")"),
  main="Variance Explained by each Component",
  col="firebrick2",
  ylim = c(0,100)
)
axis(1,at=bp[,1],labels=1:n)
mtext("MCF7 cell line",cex=1)
dev.off()

# Lets start plotting the actual PCA plots

cellLine <- gsub(' .+', '', colnames(named_newcounts))
treatment <- gsub('.+ ', '', colnames(named_newcounts))
passage <- gsub("^.*? ","", colnames(named_newcounts))
passage <- gsub(' .*?$','', passage)
experimental_data <- cbind(cellLine, treatment, passage)

MCF7_experimental_data <- as.matrix(experimental_data[grep("MCF7", experimental_data), ])
MCF7_cellLine <- MCF7_experimental_data[ ,"cellLine"] ## important for colours
MCF7_cell_line <- as.factor(MCF7_experimental_data[ ,"cellLine"])
MCF7_treatment <- MCF7_experimental_data[ ,"treatment"]
MCF7_treat <- as.factor(MCF7_treatment)

pchs<-setNames(rep(19,times = ncol(MCF7_normalized_matrix)),
               colnames(MCF7_normalized_matrix)) #### setNames(function, names to be assigned to the function)
pchs[grep("RNAiMAX",names(pchs))] <- 0
pchs[grep("siCTRL",names(pchs))] <- 2
pchs[grep("siZMIZ1",names(pchs))] <- 19

png("PCA on MCF7 cell line (non-scaled MCF7 data only).png",width=11000,height=8000,pointsize=60, res = 300)
plot(MCF7_pca$rotation[,1],MCF7_pca$rotation[,2],
     col=c("#009E73"),
     pch = pchs, cex = 1, 
     lty = "solid", 
     lwd = 2,
     main = "Principal component analysis on MCF7 cell line", 
     xlab = paste0("PC1 (",signif(MCF7_pcavar[1],3),"%)"), 
     ylab = paste0("PC2 (",signif(MCF7_pcavar[2],3),"%)"))

legend("bottomright", title="Treatment:", 
       pch = pchs, 
       legend = c("siZMIZ1", "siCTRL", "RNAiMAX"), 
       col = c("#009E73"), 
       pt.cex=1, 
       cex=0.75)
grid()
#text(MCF7_pca$rotation[,1],MCF7_pca$rotation[,2],substr(rownames(MCF7_pca$rotation),6,8),cex=0.7, pos=4)
dev.off()


png("PC2 vs PC3 on MCF7 cell line (non-scaled MCF7 data only).png",width=11000,height=8000,pointsize=60, res = 300)
plot(MCF7_pca$rotation[,2],MCF7_pca$rotation[,3],
     col=c("#009E73"),
     pch = pchs, cex = 1, 
     lty = "solid", 
     lwd = 2,
     main = "Principal component analysis on MCF7 cell line", 
     xlab = paste0("PC2 (",signif(MCF7_pcavar[2],3),"%)"), 
     ylab = paste0("PC3 (",signif(MCF7_pcavar[3],3),"%)"))
legend("topright", title="Treatment:", 
       pch = pchs, 
       legend = c("siZMIZ1", "siCTRL", "RNAiMAX"), 
       col = c("#009E73"), 
       pt.cex=1, 
       cex=0.75)
grid()
#text(MCF7_pca$rotation[,2],MCF7_pca$rotation[,3],substr(rownames(MCF7_pca$rotation),6,8),cex=0.7, pos=1)
dev.off() 


##### RUNNING GSEA ala Federico Giorgi 
# Load Federico's GSEA functions (freely developed from Subramainan et al., 200x)
source("gseaFunctions.R")

# Obtain our signature
View(MCF7_MA_results_table)
signature <- MCF7_MA_results_table$stat
names(signature) <- rownames(MCF7_MA_results_table)
signature

# Load pathways from Msigdb Broad v 5.0
whatis<-load("MSigDB_v5.0_human.rda")
allpathways<-msigDBsymbol$c2.all.v5.0.symbols.gmt
keggpathways<-msigDBsymbol$c2.cp.kegg.v5.0.symbols.gmt



# Repeat this block of code for every pathway in the allpathways list
GSEAresults<-matrix(0,nrow=0,ncol=3)
colnames(GSEAresults)<-c("Pathway Name","NES","pvalue")
collection<-allpathways
pb<-txtProgressBar(0,length(collection),style=3)
for(pathwayname in names(collection)){
  #message("Hello Anze, I am doing ", pathwayname)
  
  # Extract pathway corresponding to pathawayname
  pathway<-collection[[pathwayname]]
  set.seed(1) # Control random permutation
  
  # Calculate GSEA
  obj<-gsea(reflist=signature,set=pathway,method="pareto")
  nes<-obj$nes
  pvalue<-obj$p.value
  
  entry<-c(pathwayname,nes,pvalue)
  GSEAresults<-rbind(GSEAresults,entry)
  
  setTxtProgressBar(pb,nrow(GSEAresults))
}


save(GSEAresults, file = "MCF7 GSEA all pahways.rda")
write.csv(GSEAresults, file = "MCF7 GSEA all pathways.csv")



# Adjust your pvalues
adjusted<-p.adjust(GSEAresults[,"pvalue"],method="bonferroni")
adjustedGSEAresults<-cbind(GSEAresults,adjusted)
write.table(adjustedGSEAresults,file="adjustedGSEAresults")
save(adjustedGSEAresults, file = "MCF7_adjustedGSEAresults.rda")
write.csv(adjustedGSEAresults, file = "MCF7adjustedGSEAresults.csv")