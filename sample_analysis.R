
# DESCRIPTION: 
# This script show the computation of the charts introduced in the sample analysis
# of the report. All the charts were computed twice. The first time with raw
# intensities, the second with logFC. The report only shows the plots computed with
# raw intensity values.

# Author: Mario Rubio Chavarría


# --------------------------------------------------------------
# Libraries
# --------------------------------------------------------------
library(data.table)
library(limma)
library(biomaRt)
library(tidyverse)
library(reshape2)
library(ggfortify)
library(umap)
library(Rtsne)
library(gridExtra)


# --------------------------------------------------------------
# Functions
# --------------------------------------------------------------

prepare_log_data <- function(merged.data, base.day) {
  # Create merged.data with logFC expression
  # Take all the samples from day 0
  base <- merged.data[merged.data$day == base.day, ]
  base <- base %>%
    dplyr::select(starts_with("ENSG")) %>%
    aggregate(list(Sample = base$Sample, limma_group = base$limma_group), mean)
  # Create the data to divide
  n.numbers <- merged.data %>% count(Sample)
  for (i in 1:dim(n.numbers)[1]) {
    base.row <- base[i, ] %>% data.frame
    n.copies <- n.numbers[n.numbers$Sample == base.row$Sample, ]$n
    repeated.rows <- base.row[rep(seq_len(nrow(base.row)), each = n.copies), ]
    if (i == 1) {
      divide.data <- repeated.rows
    } else {
      divide.data <- rbind(divide.data, repeated.rows)
    }
  }
  # Divide the original expression and log it
  log.merged.data <- merged.data
  log.merged.data[ , grepl("ENSG" , colnames(log.merged.data))] <- 
    log(dplyr::select(merged.data, starts_with("ENSG")) / dplyr::select(divide.data, starts_with("ENSG")), 
        base = 2)
  return(log.merged.data)
}


# --------------------------------------------------------------
# Prepare the data
# --------------------------------------------------------------

setwd('/home/mario/Projects/az_project/')

# Load expression data
expr_data = read.csv('./data/Normalised_bone_marrow_data_Ensembl_18122020.csv')
expr_data$probe_id = NULL
rownames(expr_data) <- expr_data$Ensembl
expr_data$Ensembl <- NULL

# Load annotation data
anno.data <- fread('./data/processed_annotation.tsv')

# Lodel model results
BFUE_res <- fread('./data/BFUE_res.tsv')
GM_res <- fread('./data/GM_res.tsv')
Mk_res <- fread('./data/Mk_res.tsv')

# Translate to gene symbols
mart <- useEnsembl(biomart='ensembl', dataset='hsapiens_gene_ensembl') # Some times the Ensembl site is unresponsive.
translation <-  getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description'),
                      filters = 'ensembl_gene_id',
                      values = rownames(expr_data),
                      mart=mart)

# Combine expression and annotation data to plot charts
expr.data <- t(as.matrix(expr_data))
expr.data <- data.frame(expr.data)
expr.data$sample_ID <- rownames(expr.data)
rownames(expr.data) <- NULL
merged.data <- merge(x = anno.data, y = expr.data, by = "sample_ID")
merged.data$sample_ID <- factor(merged.data$sample_ID)
ordered.days <- 0:14
merged.data$day <- factor(merged.data$day, levels = ordered.days)
ordered.donors <- 1:4
merged.data$donor <- factor(merged.data$donor, levels = ordered.donors)
ordered.lineages <- c("BFUE", "Mk", "GM")
merged.data$lineage <- factor(merged.data$lineage, levels = ordered.lineages)
merged.data$Sample <- factor(substr(merged.data$sample_ID, 4, 5),
                             levels = c("A1","A2", "B1", "B2", "C1", "C2"))
merged.data <- merged.data %>% data.frame

# Prepare expression data in logFC form
base.day <- 0
log.merged.data <- prepare_log_data(merged.data, base.day)

# Working data with raw expression
working.data <- merged.data[ , grepl("ENSG" , colnames(merged.data))]

# Working data with logFC expression
log.working.data <- log.merged.data[ , grepl("ENSG" , colnames(log.merged.data))]


# --------------------------------------------------------------
# Dimensionality reduction techniques
# --------------------------------------------------------------

# UMAP
set.seed(7) # To make UMAP reproducible (selected value)
data.umap <- umap(working.data)
log.data.umap <- umap(log.working.data)

# PCA
data.pca <- prcomp(working.data, scale = T)
log.data.pca <- prcomp(log.working.data, scale = T)

# t-SNE
set.seed(2) # To make t-SNE reproducible (selected value)
perplexity.value <- 110
data.tsne <- Rtsne(data.pca$x, perplexity=perplexity.value)
log.data.tsne <- Rtsne(log.data.pca$x, perplexity=perplexity.value)

# --------------------------------------------------------------
# Plot samples by lineage, day, and test
# --------------------------------------------------------------

# NORMAL intensities
# Plot by lineage
labels <- merged.data$lineage
plot.data <- cbind(labels, data.frame(data.tsne$Y))
colnames(plot.data) <- c("labels", "component1", "component2")
p.lineage <- plot.data %>%
  ggplot(aes(x = component1, y = component2, colour = labels)) +
  geom_point(size = 4) +
  theme_bw() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16)) +
  scale_color_brewer(palette="Dark2") +
  labs(x = "t-SNE 1", y = "t-SNE 2", colour = "Lineage")
# Plot by day
days <- merged.data$day
lineage <- merged.data$lineage
plot.data <- cbind(lineage, cbind(as.numeric(days), data.frame(data.tsne$Y)))
colnames(plot.data) <- c("lineage", "labels", "component1", "component2")
p.day <- plot.data %>%
  ggplot(aes(x = component1, y = component2, colour = labels)) +
  geom_point(size = 4) +
  annotate(geom="text", size = 6, x=1.25, y=-5, label="Mk", color="black") +
  annotate(geom="text", size = 6, x=-4.5, y=-1.25, label="BFUE", color="black") +
  annotate(geom="text", size = 6, x=5, y=-1.8, label="GM", color="black") +
  scale_fill_brewer(palette="Paired") +
  theme_bw() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16)) +
  labs(x = "t-SNE 1", y = "t-SNE 2", colour = "Time (days)")
# Plot by test
labels <- merged.data$Sample
plot.data <- cbind(labels, data.frame(data.tsne$Y))
colnames(plot.data) <- c("labels", "component1", "component2")
p.test <- plot.data %>%
  ggplot(aes(x = component1, y = component2, colour = labels)) +
  geom_point(size = 4) +
  annotate(geom="text", size = 6, x=1.25, y=-5, label="Mk", color="black") +
  annotate(geom="text", size = 6, x=-4.5, y=-1.25, label="BFUE", color="black") +
  annotate(geom="text", size = 6, x=5, y=-1.8, label="GM", color="black") +
  theme_bw() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16)) +
  scale_color_brewer(palette="Paired") +
  labs(x = "t-SNE 1", y = "t-SNE 2", colour = "Test")
# Show plots
grid.arrange(p.lineage, p.day, p.test, ncol=3)

# LOGFC intensities
# Plot by lineage
labels <- merged.data$lineage
plot.data <- cbind(labels, data.frame(log.data.tsne$Y))
colnames(plot.data) <- c("labels", "component1", "component2")
log.p.lineage <- plot.data %>%
  ggplot(aes(x = component1, y = component2, colour = labels)) +
  geom_point(size = 4) +
  theme_bw() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16)) +
  scale_color_brewer(palette="Dark2") +
  labs(x = "t-SNE 1", y = "t-SNE 2", colour = "Lineage")
# Plot by day
days <- merged.data$day
lineage <- merged.data$lineage
plot.data <- cbind(lineage, cbind(as.numeric(days), data.frame(log.data.tsne$Y)))
colnames(plot.data) <- c("lineage", "labels", "component1", "component2")
log.p.day <- plot.data %>%
  ggplot(aes(x = component1, y = component2, colour = labels)) +
  geom_point(size = 4) +
  annotate(geom="text", size = 6, x=1.25, y=-5, label="Mk", color="black") +
  annotate(geom="text", size = 6, x=-4.5, y=-1.25, label="BFUE", color="black") +
  annotate(geom="text", size = 6, x=5, y=-1.8, label="GM", color="black") +
  scale_fill_brewer(palette="Paired") +
  theme_bw() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16)) +
  labs(x = "t-SNE 1", y = "t-SNE 2", colour = "Time (days)")
# Plot by test
labels <- merged.data$Sample
plot.data <- cbind(labels, data.frame(log.data.tsne$Y))
colnames(plot.data) <- c("labels", "component1", "component2")
log.p.test <- plot.data %>%
  ggplot(aes(x = component1, y = component2, colour = labels)) +
  geom_point(size = 4) +
  annotate(geom="text", size = 6, x=1.25, y=-5, label="Mk", color="black") +
  annotate(geom="text", size = 6, x=-4.5, y=-1.25, label="BFUE", color="black") +
  annotate(geom="text", size = 6, x=5, y=-1.8, label="GM", color="black") +
  theme_bw() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16)) +
  scale_color_brewer(palette="Paired") +
  labs(x = "t-SNE 1", y = "t-SNE 2", colour = "Test")
# Show plots
grid.arrange(log.p.lineage, log.p.day, log.p.test, ncol=3)

# --------------------------------------------------------------
# Plot biology of the donors
# --------------------------------------------------------------

# NORMAL intensities
# PCA
p.donor.pca <- autoplot(data.pca, data = merged.data, colour = 'donor', shape='donor', frame=T) +
  #labs(x = "PCA 1", y = "PCA 2") +
  theme_bw() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16)) +
  labs(colour = 'Donor', shape = 'Donor', fill = 'Donor')
# PCA 1:  14.51%
# PCA 2:  12.3%
# t-SNE
# Plot by donor
labels <- merged.data$donor
plot.data <- cbind(merged.data, data.frame(data.tsne$Y))
plot.data <- plot.data[, c("donor","X1", "X2")]
colnames(plot.data) <- c("labels", "component1", "component2")
p.donor.tsne <- plot.data %>%
  ggplot(aes(x = component1, y = component2, colour = labels)) +
  geom_point(size = 4) +
  annotate(geom="text", size = 6, x=1.25, y=-5, label="Mk", color="black") +
  annotate(geom="text", size = 6, x=-4.5, y=-1.25, label="BFUE", color="black") +
  annotate(geom="text", size = 6, x=5, y=-1.8, label="GM", color="black") +
  theme_bw() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16)) +
  scale_color_brewer(palette="Dark2") +
  labs(x = "t-SNE 1", y = "t-SNE 2", colour = "Donor")
# UMAP
labels <- merged.data$donor
plot.data <- cbind(labels, data.frame(data.umap$layout))
colnames(plot.data) <- c("labels", "component1", "component2")
p.donor.umap <- plot.data[merged.data$day == 0, ] %>%
  ggplot(aes(x = component1, y = component2, colour = labels)) +
  geom_point(size = 2) +
  geom_density2d(contour_var = "density", alpha = 0.5, show.legend = F) +
  # annotate(geom="text", size = 8, x=-1.1, y=7.5, label="Mk", color="black") +
  # annotate(geom="text", size = 8, x=-3.5, y=-5.5, label="BFUE", color="black") +
  # annotate(geom="text", size = 8, x=2.5, y=-4.4, label="GM", color="black") +
  theme_bw() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16)) +
  labs(x = "UMAP 1", y = "UMAP 2", colour = "Donor")
# Show plots
grid.arrange(p.donor.pca, p.donor.tsne, p.donor.umap, ncol=3)

# LOGFC intensities
# PCA
log.p.donor.pca <- autoplot(log.data.pca, data = merged.data, colour = 'donor', shape='donor', frame=T) +
  labs(x = "PCA 1", y = "PCA 2") +
  theme_bw() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16)) +
  labs(colour = 'Donor', shape = 'Donor', fill = 'Donor')
# PCA 1:  14.51%
# PCA 2:  12.06%
# t-SNE
# Plot by donor
labels <- merged.data$donor
plot.data <- cbind(merged.data, data.frame(log.data.tsne$Y))
plot.data <- plot.data[, c("donor","X1", "X2")]
colnames(plot.data) <- c("labels", "component1", "component2")
log.p.donor.tsne <- plot.data %>%
  ggplot(aes(x = component1, y = component2, colour = labels)) +
  geom_point(size = 4) +
  annotate(geom="text", size = 6, x=1.25, y=-5, label="Mk", color="black") +
  annotate(geom="text", size = 6, x=-4.5, y=-1.25, label="BFUE", color="black") +
  annotate(geom="text", size = 6, x=5, y=-1.8, label="GM", color="black") +
  theme_bw() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16)) +
  scale_color_brewer(palette="Dark2") +
  labs(x = "t-SNE 1", y = "t-SNE 2", colour = "Donor")
# UMAP
labels <- merged.data$donor
plot.data <- cbind(labels, data.frame(log.data.umap$layout))
colnames(plot.data) <- c("labels", "component1", "component2")
log.p.donor.umap <- plot.data[merged.data$day == 0, ] %>%
  ggplot(aes(x = component1, y = component2, colour = labels)) +
  geom_point(size = 2) +
  geom_density2d(contour_var = "density", alpha = 0.5, show.legend = F) +
  # annotate(geom="text", size = 8, x=-1.1, y=7.5, label="Mk", color="black") +
  # annotate(geom="text", size = 8, x=-3.5, y=-5.5, label="BFUE", color="black") +
  # annotate(geom="text", size = 8, x=2.5, y=-4.4, label="GM", color="black") +
  theme_bw() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16)) +
  labs(x = "UMAP 1", y = "UMAP 2", colour = "Donor")
# Show plots
grid.arrange(log.p.donor.pca, log.p.donor.tsne, log.p.donor.umap, ncol=3)


