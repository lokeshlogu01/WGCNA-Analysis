# Load necessary libraries
library(WGCNA)
library(GEOquery)
library(tidyverse)
library(gridExtra)
library(corrplot)
library(limma)
library(genefilter)
library(preprocessCore)

# Enable multi-threading for WGCNA
allowWGCNAThreads()

# Load GEO dataset
geo_data <- getGEO("GSE41011")
datExpr <- exprs(geo_data[[1]])
expression_data <- datExpr

# Extract phenotype data
phenoData <- pData(geo_data[[1]])

# Check for good samples and genes
gsg <- goodSamplesGenes(t(expression_data), verbose = 3)
summary(gsg)

# Hierarchical clustering of samples
htree <- hclust(dist(t(datExpr)), method = "average")
plot(htree, main = "Sample Clustering Dendrogram")

# Transpose expression data for WGCNA
exp.data <- t(expression_data)

# Pick a suitable soft-thresholding power
power <- c(1:10, seq(12, 20, by = 2))
sft <- pickSoftThreshold(exp.data, powerVector = power, verbose = 5)

# Plot scale-free topology model fit
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale-free Topology Model Fit, signed R^2",
     main = "Scale Independence", col = "red")
abline(h = 0.90, col = "blue")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (Power)",
     ylab = "Mean Connectivity",
     main = "Mean Connectivity", col = "red")

# Choose the appropriate soft power
soft_power <- 20

# Build the network
bwnet <- blockwiseModules(exp.data,
                          power = soft_power,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          networkType = "signed",
                          deepSplit = 2,
                          pamRespectsDendro = FALSE,
                          mergeCutHeight = 0.25,
                          numericLabels = TRUE,
                          saveTOMs = TRUE,
                          saveTOMFileBase = "WGCNA",
                          verbose = 3)

# Convert module labels to colors
mergedcolors <- labels2colors(bwnet$colors)

# Plot dendrogram with module colors
plotDendroAndColors(
  bwnet$dendrograms[[1]], 
  mergedcolors[bwnet$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)

# Save module assignments
module_df <- data.frame(gene_id = names(bwnet$colors), colors = mergedcolors)
write.csv(module_df, file = "WGCNA_modules.csv", row.names = FALSE)

# Extract module eigengenes
MElist <- moduleEigengenes(exp.data, colors = mergedcolors)
MEs <- MElist$eigengenes

# Cluster module eigengenes
ME.dissimilarity <- 1 - cor(MEs, use = "complete")
METree <- hclust(as.dist(ME.dissimilarity), method = "average")
plot(METree, main = "Clustering of Module Eigengenes")
abline(h = 0.4, col = "red")

# Merge similar modules
merge <- mergeCloseModules(exp.data, mergedcolors, cutHeight = 0.3)
mergedColors <- merge$colors
writeLines(mergedColors, "mergedModules.txt")

# Prepare trait data
traitData <- phenoData
traitdata <- traitData %>%
  mutate(stage1 = as.numeric(grepl("paired_tumor_stage: 1", characteristics_ch1.4)),
         stage2 = as.numeric(grepl("paired_tumor_stage: 2", characteristics_ch1.4)),
         stage3 = as.numeric(grepl("paired_tumor_stage: 3", characteristics_ch1.4)),
         stage4 = as.numeric(grepl("paired_tumor_stage: 4", characteristics_ch1.4)))

datTraits <- traitdata[, c("stage1", "stage2", "stage3", "stage4")]

# Calculate module-trait relationships
module.trait.correlation <- cor(MEs, datTraits, use = "p")
module.trait.Pvalue <- corPvalueStudent(module.trait.correlation, nrow(exp.data))

# Create labeled heatmap
textMatrix <- paste(signif(module.trait.correlation, 2), "\n(",
                   signif(module.trait.Pvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(module.trait.correlation)

labeledHeatmap(Matrix = module.trait.correlation,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.4,
               zlim = c(-1,1),
               main = "Module-Trait Relationships")
