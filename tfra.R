
# DESCRIPTION: 

# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------
library(data.table)
library(limma)
library(biomaRt)
library(tidyverse)
library(readxl)
library(Mfuzz)
library(stringr)
library(Biobase)
library(fame)
library(ggfortify)
library(ggplot2)
library(Rtsne)
library(umap)
library(GSA)
library(dorothea)


# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------
kmeansBIC <- function(fit){
  # Source: http://sherrytowers.com/2013/10/24/k-means-clustering/
  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(D + log(n)*m*k)
}


# ------------------------------------------------------------------------------
# Prepare the data
# ------------------------------------------------------------------------------
setwd('/home/mario/Projects/az_project/')

# Load annotation data
anno.data <- fread("./data/processed_annotation.tsv")

# Load expression data
expr.data <- read.csv('./data/Normalised_bone_marrow_data_Ensembl_18122020.csv')

# Translate to gene symbols
mart <- useEnsembl(biomart='ensembl', dataset='hsapiens_gene_ensembl')
translation <-  getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                      filters = 'ensembl_gene_id',
                      values = expr.data$Ensembl,
                      mart=mart)
colnames(translation) <- c("Ensembl", "hgnc_symbol")
genes <- merge(x = expr.data, y = translation, by = 'Ensembl')
genes <- genes[, c("Ensembl", "hgnc_symbol")]


# ------------------------------------------------------------------------------
# Create the gene clusters for all the lineages
# ------------------------------------------------------------------------------
lineages <- c('BFUE', 'GM', 'Mk')
first.lineage <- T
for (selected.lineage in lineages) {
  filtered.anno <- anno.data %>% dplyr::filter(lineage == selected.lineage,
                                               sample_ID %in% colnames(expr.data))
  filtered.anno <- filtered.anno[order(filtered.anno$day), ]
  selected.samples <- filtered.anno$sample_ID
  lineage.expr.data <- expr.data[, selected.samples]
  
  # Compute log2FC with respect day 0
  filtered.anno.0 <- filtered.anno %>% dplyr::filter(day == 0)
  first <- T
  for (local.donor in unique(filtered.anno$donor)) {
    # Iterate for every donor to divide by the reference for that donor
    donor.samples.0 <- filtered.anno.0[filtered.anno.0$donor == local.donor, ]$sample_ID
    donor.expr.0 <- expr.data[, unlist(donor.samples.0)]
    donor.reference <- data.frame(rowMeans(donor.expr.0))
    donor.samples <- filtered.anno[filtered.anno$donor == local.donor, ]$sample_ID
    donor.expr <- expr.data[, unlist(donor.samples)]
    donor.reference <- donor.reference[, rep(1, times=dim(donor.expr)[2])]
    colnames(donor.reference) <- colnames(donor.expr)
    log.donor.expr <- log(donor.expr / donor.reference, base = 2)
    if (first) {
      log.lineage.expr.data <- log.donor.expr
      first <- F
    } else {
      log.lineage.expr.data <- cbind(log.lineage.expr.data, log.donor.expr)
    }
  }
  
  # Average donors
  first <- T
  for (local.day in unique(filtered.anno$day)) {
    samples <- filtered.anno[filtered.anno$day == local.day, ]$sample_ID
    columns <- log.lineage.expr.data[, samples]
    day.means <- rowMeans(columns)
    if (first) {
      avg.log.lineage.expr.data <- data.frame("0" = day.means)
      first <- F
    } else {
      avg.log.lineage.expr.data <- cbind(avg.log.lineage.expr.data, day.means)
    }
  }
  colnames(avg.log.lineage.expr.data) <- unique(filtered.anno$day)
  
  # Compute the z-score for every column
  z.score.data <- apply(avg.log.lineage.expr.data, 2, function(x) {(x - mean(x))/sd(x)})
  
  # Compute the k-means with the chosen parameters
  selected.n.centres <- 7
  selected.nstart <- 50
  iter.max <- 1000
  computed.kmeans <- kmeans(z.score.data, selected.n.centres, iter.max, selected.nstart)
  
  # Prepare data for dorothea
  dorothea.data <- avg.log.lineage.expr.data
  dorothea.data$Ensembl <- expr.data$Ensembl
  dorothea.data$cluster <- computed.kmeans$cluster
  dorothea.data <- dorothea.data %>% melt(id = c("cluster", "Ensembl"))
  colnames(dorothea.data) <- c("cluster", "gene", "day", "intensity")
  
  # DOROTHEA
  filtered.dorothea_hs <- dorothea_hs[!(dorothea_hs$confidence %in% c("D", "E")), ] 
  tfactors <- unique(filtered.dorothea_hs$tf)
  # Translate to hgnc symbols
  doro.data <- merge(x = dorothea.data, y = translation, by.x = "gene", by.y = "Ensembl")
  doro.data <- doro.data[order(doro.data$cluster, doro.data$day), ]
  # Compute contingency table. Note that there is a test for every donor and 
  first <- T
  contingency.table <- data.frame(tfactor = "", cluster = 0, N = 0, k = 0, m = 0, n = 0, x = 0, p.value = 0)
  # Number of genes among all the regulons
  N <- length(filtered.dorothea_hs$target)
  for (n.cluster in 1:length(unique(doro.data$cluster))) {
    # Genes in the cluster
    cluster.genes <- doro.data[doro.data$cluster == n.cluster, ]$hgnc_symbol
    # Number of genes in the cluster
    k <- length(unique(doro.data[doro.data$cluster == n.cluster, ]$gene))
    for (tfactor in tfactors) {
      # Number of genes in the regulon
      m <- length(filtered.dorothea_hs[filtered.dorothea_hs$tf == tfactor, ]$target)
      # Number of genes not associated with the regulon
      n <- N - m
      # Regulon genes
      tfactor.genes <- filtered.dorothea_hs[filtered.dorothea_hs$tf == tfactor, ]$target
      # Number of regulon genes that are in the cluster
      x <- sum(tfactor.genes %in% cluster.genes)
      # Hypergeometric test: p-value
      p_value <-  phyper(q=x - 1, m=m, n=n, k=k, lower.tail=FALSE)
      # Make contingency table
      contingency.table <- rbind(contingency.table, c(tfactor, n.cluster, N, k, m, n, x, p_value))
    }
  }
  contingency.table <- contingency.table[2:dim(contingency.table)[1], ]
  # Correct for multiple testing
  contingency.table$adjusted.p.value <- p.adjust(contingency.table$p.value, method = "BH", n = length(contingency.table$p.value))
  
  # Obtain and store the relevant results
  cutoff <- 0.25
  lineage.relevant <- contingency.table[contingency.table$adjusted.p.value <= cutoff, ]
  lineage.relevant$lineage <- rep(paste(unlist(str_split(selected.lineage, "-")), collapse = ""), times=dim(lineage.relevant)[1])
  
  # Store the transcription factor results
  if (selected.lineage == "BFUE") {
    relevant <- lineage.relevant
  } else {
    relevant <- rbind(relevant, lineage.relevant)
  }
  
  # Store gene clustering
  dorothea.data$lineage <- rep(selected.lineage, times=dim(dorothea.data)[1])
  if (first.lineage) {
    first.lineage <- F
    gene.data <- dorothea.data
  } else {
    gene.data <- rbind(gene.data, dorothea.data)
  }
}


# ------------------------------------------------------------------------------
# Print the clusters for the selected lineage
# ------------------------------------------------------------------------------
selected.lineage <- "BFUE"
print(paste("Lineage:", selected.lineage))
for (selected.cluster in 1:7) {
  print(paste("Cluster:", selected.cluster))
  filtered.gene.data <- gene.data %>% 
    dplyr::filter(lineage == selected.lineage, cluster == selected.cluster)
  trend.data <- filtered.gene.data %>%
    dplyr::filter(lineage == selected.lineage, cluster == selected.cluster) %>%
    dplyr::group_by(day) %>%
    dplyr::summarise(trend = median(intensity)) %>% data.frame
  plot <- filtered.gene.data %>% 
    ggplot + 
    geom_line(mapping = aes(x = day, y = intensity, group = gene), alpha = 0.5, color = "lightblue") +
    geom_line(data = trend.data, mapping = aes(x = day, y = trend, group = 1), color = "red") +
    labs(x = "Time (days)", y = "log2 fold-change") + 
    theme_bw() +
    theme(text = element_text(size=14),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14)) 
  print(plot)
}

