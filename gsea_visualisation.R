
# DESCRIPTION:
# This scripts shows the evolution of the genes in the results obtained from GSEA.
# It is shown the temporal evolution of all the genes in the gene set. These gene
# sets are shown lineage.

# Author: Mario Rubio Chavarría


# Libraries
library(data.table)
library(limma)
library(biomaRt)
library(tidyverse)
library(reshape2)
library(ggfortify)
library(umap)
library(gridExtra)


# --------------------------------------------------------------
# Functions
# --------------------------------------------------------------
g_legend<-function(a.gplot){
  # DESCRIPTION: Function to extract a legend from a plot.
  # Source: https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
  # :param a.gplot: [ggplot] the plot from which the legend is to be taken.
  # :return: [gDesc] plot legend.
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

prepare_log_data <- function(merged.data, base.day) {
  # DESCRIPTION: changes the data provided in raw intensities to log2 fold-change.
  # In essence, for each donor this function computes the average of gene 
  # expression for the day base.day. Afterwards, all the intensity dataset is
  # divided by the average (by donor). Finally, the function computes the log2 
  # of the ratio.
  # :param merged.data: [data.frame] intensity dataset with the annotations merged.
  # :param base.day: [numeric] day used to compute the average of every donor.
  # :return: [data.frame] merged.data with log2 fold-change intensities. 
  
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

print_main_pathways <- function(selected.lineage, pathways.folder, pathways.file, divided, n.genes.pathway, p.value.limit, translation, days, no.title, base.data, log.base.data) {
  # Read pathways data
  y.label <- "Intensity"
  log.y.label <- "log2 fold-change"
  initial.day <- days[1]
  last.day <- days[2]
  pathways.data <- read.csv(paste(pathways.folder, pathways.file, sep=""), sep="\t")
  # Filter pathways (size and p-value) and order by normalised NES and take the first ten
  pathways.data <- pathways.data %>% filter(FDR.q.val < p.value.limit)
  pathways.data <- pathways.data[order(-pathways.data$NES), ]
  # Select the n most relevant pathways (we plot all the pathways)
  relevant.pathways <- pathways.data
  # Select the genes for every pathway
  print(paste("Pathways to plot: ", dim(relevant.pathways)[1], sep=""))
  print(paste("day: ", unlist(str_split(pathways.file, "_"))[5]))
  print(paste("Lineage: ", selected.lineage))
  print("The list is printed in inverse order")
  # The printing is in reversed order to plot the most relevant chart the last, 
  # what will make the reading easier.
  for (i in dim(relevant.pathways)[1]:1) {
    
    # Obtain the genes of the leading edge of the gene set
    pathway.name <- pathways.data$NAME[i]
    if (is.na(pathway.name)) {
      print("No pathways found")
      break
    }
    pathway.data <- read.csv(paste(pathways.folder, pathway.name, ".tsv", sep=""), sep="\t", fill = FALSE, header = FALSE)
    name.cols <- pathway.data[1, !is.na(pathway.data[1, ])]
    pathway.data <- pathway.data[2:dim(pathway.data)[1], ]
    colnames(pathway.data) <- name.cols
    pathway.data <- pathway.data[, c(2, 5, 6, 7)]
    colnames(pathway.data)[2] <- "RANK_METRIC_SCORE"
    colnames(pathway.data)[4] <- "CORE_ENRICHMENT"
    pathway.data <- pathway.data %>% filter(CORE_ENRICHMENT == "Yes")
    pathway.data <- pathway.data[order(-abs(as.numeric(pathway.data$RANK_METRIC_SCORE))), ]
    pathway.data$IMPORTANT <- 0
    pathway.data[1:5, ]$IMPORTANT <- seq(1, 5)
    pathway.data$IMPORTANT <- factor(pathway.data$IMPORTANT, levels = c(0, 1, 2, 3, 4, 5))
    # Select all the leading edge or only a subset
    pathway.data <- if (n.genes.pathway != 'all') pathway.data[1:n.genes.pathway, ] else pathway.data
    remark.data <- merge(x = pathway.data[, c(1, 5)], y = translation,
                         by.x = "SYMBOL", by.y = "hgnc_symbol")
    remark.data <- remark.data[, c(1, 2, 3)]
    # Obtain the expression of those genes 
    genes.names <- merge(x = pathway.data,
                         y = translation,
                         by.x = "SYMBOL",
                         by.y = "hgnc_symbol")$ensembl_gene_id
    genes.data <- cbind(base.data[, !grepl("ENSG" , colnames(log.base.data))],
                        base.data[, genes.names]) %>% dplyr::filter(lineage == selected.lineage)
    log.genes.data <- cbind(log.base.data[, !grepl("ENSG" , colnames(log.base.data))],
                            log.base.data[, genes.names]) %>% dplyr::filter(lineage == selected.lineage)
    genes.means <- aggregate(genes.data[, genes.names], list(genes.data$day), mean)
    log.genes.means <- aggregate(log.genes.data[, genes.names], list(genes.data$day), mean)
    colnames(genes.means)[1] <- "day"
    colnames(log.genes.means)[1] <- "day"
    plot.data <- melt(genes.means, id=c("day")) %>% merge(y = remark.data,
                                                          by.x = "variable",
                                                          by.y = "ensembl_gene_id")
    log.plot.data <- melt(log.genes.means, id=c("day")) %>% merge(y = remark.data,
                                                                  by.x = "variable",
                                                                  by.y = "ensembl_gene_id")
    plot.data$day <- as.integer(as.integer(plot.data$day) - 1)
    remark.data <- plot.data %>% filter(IMPORTANT != 0)
    plot.data <- plot.data %>% filter(IMPORTANT == 0)
    log.plot.data$day <- as.integer(as.integer(log.plot.data$day) - 1)
    log.remark.data <- log.plot.data %>% filter(IMPORTANT != 0)
    log.plot.data <- log.plot.data %>% filter(IMPORTANT == 0)
    # Plot the evolution for all the genes (raw data)
    remark.data <- remark.data %>% dplyr::filter(day <= last.day)
    p.raw <- plot.data %>%
      dplyr::filter(day <= last.day) %>%
      ggplot() + 
      geom_line(mapping = aes(x = day, y = value, group = variable), color = "lightblue", alpha = 0.65) +
      geom_line(data = remark.data, mapping = aes(x = day, y = value, group = variable, color = SYMBOL)) +
      scale_color_brewer(palette="Dark2") +
      theme_bw() + 
      theme(text = element_text(size=14),
            axis.text.x = element_text(size=14),
            axis.text.y = element_text(size=14)) +
      labs(x = "Time(days)", y = y.label, colour = "Most enriched") +
      # guides(color = F) +
      xlim(initial.day, last.day)
    # Plot the evolution for all the genes (logFC)
    log.remark.data <- log.remark.data %>% dplyr::filter(day <= last.day)
    p.log <-log.plot.data %>%
      dplyr::filter(day <= last.day) %>%
      ggplot() + 
      geom_line(mapping = aes(x = day, y = value, group = variable), color = "lightblue", alpha = 0.65) +
      geom_line(data = log.remark.data, mapping = aes(x = day, y = value, group = variable, color = SYMBOL)) +
      scale_color_brewer(palette="Dark2") +
      theme_bw() + 
      theme(text = element_text(size=14),
            axis.text.x = element_text(size=14),
            axis.text.y = element_text(size=14)) +
      labs(x = "Time(days)", y = log.y.label, colour = "Main genes") +
      xlim(initial.day, last.day)
    # Show the pathway printed
    print(pathway.name)
    # Build and print the final chart
    if (no.title) {pathway.name <- ""}
    if (divided) {
      print(p.raw + ggtitle(pathway.name) + guides(colour = F))
      print(p.log + ggtitle(pathway.name))
    } else {
      legend <- g_legend(p.log)
      grid.arrange(
        arrangeGrob(p.raw + theme(legend.position="none"),
                    p.log + theme(legend.position="none"),
                    nrow=1),
        legend, ncol=2, widths=c(10, 1), heights=c(10, 1), top = pathway.name)
    }
  }
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

# Translate to gene symbols
mart <- useEnsembl(biomart='ensembl', dataset='hsapiens_gene_ensembl') # Some times the Ensembl site is unresponsive. 
translation <-  getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description'),
                      filters = 'ensembl_gene_id',
                      values = rownames(expr_data),
                      mart=mart)

# Load annotation data
anno.data <- fread('./data/processed_annotation.tsv')

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


# --------------------------------------------------------------
# Plots by lineage of all the pathways (GSEA)
# --------------------------------------------------------------

# Erythrocytes
lineage <- "BFUE"
pathways.folder.bfue <- "./gsea_results/bfue_gsea_results/"
n.genes.pathway <- 'all'
p.value.limit <- 0.15
translation <- translation
days <- c(0, 10)
# More enriched in Day 0
pathways.file.bfue.0 <- "gsea_report_for_Day_0_BFUE_1612444397464.tsv"
no.title <- F
divided <- F
print_main_pathways(lineage, pathways.folder.bfue, pathways.file.bfue.0, divided,
                    n.genes.pathway, p.value.limit, translation, days, no.title,
                    merged.data, log.merged.data)
# More enriched in Day 10
pathways.file.bfue.10 <- "gsea_report_for_Day_10_BFUE_1612444397464.tsv"
no.title <- F
divided <- F
print_main_pathways(lineage, pathways.folder.bfue, pathways.file.bfue.10, divided,
                    n.genes.pathway, p.value.limit, translation, days, no.title,
                    merged.data, log.merged.data)

# Myeloid
lineage <- "GM"
pathways.folder.gm <- "./gsea_results/gm_gsea_results/"
n.genes.pathway <- 'all'
p.value.limit <- 0.15
translation <- translation
days <- c(0, 14)
# More enriched in Day 0
pathways.file.gm.0 <- "gsea_report_for_Day_0_GM_1612444535888.tsv"
no.title <- F
divided <- F
print_main_pathways(lineage, pathways.folder.gm, pathways.file.gm.0, divided,
                    n.genes.pathway, p.value.limit, translation, days, no.title,
                    merged.data, log.merged.data)
# More enriched in Day 14
pathways.file.gm.14 <- "gsea_report_for_Day_14_GM_1612444535888.tsv"
no.title <- F
divided <- F
print_main_pathways(lineage, pathways.folder.gm, pathways.file.gm.14, divided,
                    n.genes.pathway, p.value.limit, translation, days, no.title,
                    merged.data, log.merged.data)

# Megakaryocytes
lineage <- "Mk"
pathways.folder.mk <- "./gsea_results/mk_gsea_results/"
n.genes.pathway <- 'all'
p.value.limit <- 0.25
translation <- translation
days <- c(0, 14)
# More enriched in Day 0
pathways.file.mk.0 <- "gsea_report_for_Day_0_Mk_1612444630696.tsv"
no.title <- F
divided <- F
print_main_pathways(lineage, pathways.folder.mk, pathways.file.mk.0, divided,
                    n.genes.pathway, p.value.limit, translation, days, no.title,
                    merged.data, log.merged.data)
# More enriched in Day 14
pathways.file.mk.14 <- "gsea_report_for_Day_14_Mk_1612444630696.tsv"
no.title <- F
divided <- F
print_main_pathways(lineage, pathways.folder.mk, pathways.file.mk.14, divided,
                    n.genes.pathway, p.value.limit, translation, days, no.title,
                    merged.data, log.merged.data)

