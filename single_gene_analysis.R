
# --------------------------------------------------------------
# Libraries
# --------------------------------------------------------------
library(data.table)
library(limma)
library(biomaRt)
library(tidyverse)
library(gridExtra)
library(readxl)
library(reshape2)
library(ggfortify)
library(writexl)
library(umap)
library(reshape)
library(Rtsne)
library(ggrepel)
library(gridExtra)
library(fame)
library(wesanderson)


# --------------------------------------------------------------
# Functions
# --------------------------------------------------------------

prepare_data <- function(lineage, anno.data, merged.data, days) {
  selected <- merged.data %>%
    dplyr::filter(Lineage == lineage) %>%
    dplyr::filter(Day %in% days)
  selected <- selected %>%
    dplyr::select(starts_with("ENSG")) %>%
    aggregate(list(selected$limma_group), mean)
  avg.selected <- t(selected)
  avg.selected <- data.frame(avg.selected)
  avg.selected$Gene <- rownames(avg.selected)
  colnames(avg.selected) <- c(avg.selected[1, 1], avg.selected[1, 2], "ensembl_gene_id")
  rownames(avg.selected) <- NULL
  avg.selected <- avg.selected[2:dim(avg.selected)[1], ]
  avg.selected[, 1] <- as.numeric(avg.selected[, 1])
  avg.selected[, 2] <- as.numeric(avg.selected[, 2])
  avg.selected$M <- log(avg.selected[, 2],2) - log(avg.selected[, 1],2)
  avg.selected$A <- 0.5 * (log(avg.selected[, 2],2) + log(avg.selected[, 1],2))
  avg.selected$ensembl_gene_id <- substr(avg.selected$ensembl_gene_id, 1, 15)
  merged.selected <- merge(avg.selected, anno.data, by = "ensembl_gene_id")
}

print_MA_plot <- function(plot.data, threshold.adj.P.Val, n.relevants, no.legend, no.y) {
  plot.data <- plot.data[order(-abs(plot.data$logFC)), ]
  # Obtain the relevant dots to plot
  relevants <- plot.data[1:n.relevants, ]
  # Plot data
  if (no.legend) {
    p <- plot.data %>%
      ggplot(aes(x = A, y = M)) +
      geom_point(aes(colour = adj.P.Val < threshold.adj.P.Val, alpha = 0.5), size = 4) + 
      scale_colour_manual(values=c("#5bb9f3", "#0d3147")) +
      guides(size = F, alpha = F, colour = F) +
      labs(colour = paste("Adjusted P-Value <", threshold.adj.P.Val, sep = " ")) + 
      theme_bw() +
      theme(text = element_text(size=20),
            axis.text.x = element_text(size=20),
            axis.text.y = element_text(size=20)) 
  } else {
    p <- plot.data %>%
      ggplot(aes(x = A, y = M)) +
      geom_point(aes(colour = adj.P.Val < threshold.adj.P.Val, alpha = 0.5), size = 4) + 
      scale_colour_manual(values=c("#5bb9f3", "#0d3147")) +
      guides(size = F, alpha = F) +
      labs(colour = paste("Adjusted P-Value <", threshold.adj.P.Val, sep = " ")) + 
      theme_bw() +
      theme(text = element_text(size=20),
            axis.text.x = element_text(size=20),
            axis.text.y = element_text(size=20)) 
  }
  if (no.y) {
    p <- p + ylab("")
  }
  print(p)
  return(p)
}

print_volcano_plot <- function(plot.data, threshold.adj.P.Val, n.relevants) {
  plot.data$mlogAdj <- -log(plot.data$adj.P.Val, 10)
  plot.data <- plot.data[order(-abs(plot.data$logFC)), ]
  # Obtain the relevant dots to plot
  relevants <- plot.data %>% dplyr::filter(adj.P.Val < 5E-5)
  relevants <- plot.data[1:n.relevants, ]
  # Plot data
  p <- plot.data %>%
    ggplot(aes(x = logFC, y = mlogAdj)) +
    # geom_point(aes(colour = adj.P.Val < threshold.adj.P.Val, alpha = 0.5), size = 0.75) + 
    geom_point(alpha = 1, size = 4, color = "lightblue") + 
    geom_text_repel(data = relevants,
                    size = 5,
                    segment.alpha = 1,
                    segment.size = 0.1,
                    max.iter = 7e7,
                    box.padding = unit(0.3, "lines"),
                    point.padding = unit(0.0, "lines"),
                    arrow = arrow(length = unit(0.01, "inches"), type = "closed", ends = "first"),
                    aes(x = logFC, y = mlogAdj, label = hgnc_symbol)) +
    # scale_size_manual(values=c(0.01, 0.8)) +
    # scale_alpha_manual(values=c(0.25, 0.5)) +
    # scale_colour_manual(values=c("#5bb9f3", "#0d3147")) +
    guides(size = F, alpha = F) +
    # labs(x = "log2 fold-change", y = "-log10 adjusted p-value", colour = paste("Adjusted P-Value <", threshold.adj.P.Val, sep = " ")) + 
    labs(x = "log2 fold-change", y = "-log10 adjusted p-value") + 
    theme_bw() + 
    theme(text = element_text(size=20),
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20)) 
  print(p)
  return(p)
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
BFUE_res <- fread('/home/mario/Projects/az_project/data/BFUE_res.tsv')
GM_res <- fread('/home/mario/Projects/az_project/data/GM_res.tsv')
Mk_res <- fread('/home/mario/Projects/az_project/data/Mk_res.tsv')

# Translate to gene symbols
mart <- useEnsembl(biomart='ensembl', dataset='hsapiens_gene_ensembl') # Some times the Ensembl site is unresponsive. 
translation <-  getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description'),
                      filters = 'ensembl_gene_id',
                      values = rownames(expr_data),
                      mart=mart)

BFUE_hgnc <- merge(x = BFUE_res, y = translation, by = "ensembl_gene_id")
BFUE_hgnc <- BFUE_hgnc[order(-abs(BFUE_hgnc$logFC)), ]
Mk_hgnc <- merge(x = Mk_res, y = translation, by = "ensembl_gene_id")
Mk_hgnc <- Mk_hgnc[order(-abs(Mk_hgnc$logFC)), ]
GM_hgnc <- merge(x = GM_res, y = translation, by = "ensembl_gene_id")
GM_hgnc <- GM_hgnc[order(-abs(GM_hgnc$logFC)), ]

# Combine expression and annotation data to plot charts
expr.data <- t(as.matrix(expr_data))
expr.data <- data.frame(expr.data)
expr.data$Sample_ID <- rownames(expr.data)
rownames(expr.data) <- NULL
merged.data <- merge(x = anno.data, y = expr.data, by = "Sample_ID")
merged.data$Sample_ID <- factor(merged.data$Sample_ID)
ordered.days <- 0:14
merged.data$Differentation_day <- substr(merged.data$Differentation_day, 4, 5)
merged.data$Day <- factor(merged.data$Differentation_day, levels = ordered.days)
merged.data$Differentation_day <- NULL
merged.data$Donor <- substr(merged.data$Donor, 6, 6)
ordered.donors <- 1:4
merged.data$Donor <- factor(merged.data$Donor, levels = ordered.donors)
ordered.lineages <- c("BFUE", "Mk", "GM")
merged.data$Lineage <- factor(merged.data$Lineage, levels = ordered.lineages)
merged.data$Sample <- factor(substr(merged.data$Sample_ID, 4, 5),
                             levels = c("A1","A2", "B1", "B2", "C1", "C2"))
merged.data <- merged.data %>% data.frame

# Reshape the data for the MA and volcano plots. It takes some time
days <- c(0, 10)
lineage <- "BFUE"
anno.data <- BFUE_hgnc
bfue.plot.data <- prepare_data(lineage, anno.data, merged.data, days) 
days <- c(0, 14)
lineage <- "GM"
anno.data <- GM_hgnc
gm.plot.data <- prepare_data(lineage, anno.data, merged.data, days)
days <- c(0, 14)
lineage <- "Mk"
anno.data <- Mk_hgnc
mk.plot.data <- prepare_data(lineage, anno.data, merged.data, days)

# --------------------------------------------------------------
# MA plots by lineage
# --------------------------------------------------------------

# Common parameters
threshold <- 1e-5
n.relevants <- 25

# Erythrocytes
no.legend <- F
no.y <- F
p.ma.bfue <- print_MA_plot(bfue.plot.data, threshold, n.relevants, no.legend, no.y)

# Myeloid
no.legend <- F
no.y <- F
p.ma.gm <- print_MA_plot(gm.plot.data, threshold, n.relevants, no.legend, no.y)

# Megakaryocytes
no.legend <- F
no.y <- F
p.ma.mk <- print_MA_plot(mk.plot.data, threshold, n.relevants, no.legend, no.y)


# --------------------------------------------------------------
# Volcano plots by lineage
# --------------------------------------------------------------

# Common parameters
threshold <- 1e-5
n.relevants <- 25

# Erythrocytes
p.volcano.bfue <- print_volcano_plot(bfue.plot.data, threshold, n.relevants)

# Myeloid
p.volcano.gm <- print_volcano_plot(gm.plot.data, threshold, n.relevants)

# Megakaryocytes
p.volcano.mk <- print_volcano_plot(mk.plot.data, threshold, n.relevants)


# # --------------------------------------------------------------
# # Print the plots by lineage
# # --------------------------------------------------------------
# # NOTE: there are no overlaps when the plots are printed independently
# 
# # Erythrocytes
# grid.arrange(p.ma.bfue, p.volcano.bfue, ncol=2)
# 
# # Myeloid
# grid.arrange(p.ma.gm, p.volcano.gm, ncol=2)
# 
# # Megakaryocytes
# grid.arrange(p.ma.mk, p.volcano.mk, ncol=2)


