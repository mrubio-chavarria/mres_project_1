
# DESCRIPTION: 
# This file gathers the code used for the computation of the KEGG  analysis.
# Additionally, a similar analysis was performed with GO, although it is not
# finally introduced in the analysis. The code is in this file too.

# Author: Mario Rubio Chavarría


# Libraries
library(data.table)
library(limma)
library(biomaRt)
library(tidyverse)
library(ggfortify)
library(stringr)
library(DescTools)


# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

findMid <- function(dtf) {
  # DESCRIPTION:
  # This function computes the mean of the ratios of the given dataset, the data
  # used in kegga_analysis (below). This value is used to compute the dotted red
  # line in the summary chart.
  # :param dtf: [data.frame] the dataset computed by kegga. Essentially, any 
  # data.frame with a numeric ratio column.
  # :return: [numeric] the mean between the maximum and minimum value of the
  # ratios.
  n_results <- dim(dtf)[1]
  rownames(dtf) <- seq(1, n_results, 1)
  ratios <- dtf$ratio
  m1 <- max(ratios)
  m2 <- min(ratios)
  ratio_level <- mean(c(m1, m2))
  r1 <- rownames(dtf[dtf$ratio >= ratio_level, ])[1]
  r2 <- tail(rownames(dtf[dtf$ratio < ratio_level, ]), 1)
  height <- mean(as.integer(c(r1, r2)))
  return(height)
}

kegga_analysis <- function(entrezgene_ids, species, lower_limit, upper_limit, lineage, n_results = NULL, fdr.threshold = 0.25, no.title = F, no.legend = F) {
  # DESCRIPTION:
  # Wrapper of the test for over representation of KEGG pathways (kegga).
  # :param entrezgene_ids: [vector] the list of genes that should be compared against KEGG.
  # :param species: [string] the code for the species whose pathways should be compared.
  # :param lower_limit: [integer] minimum size of the pathways to compare.
  # :param upper_limit: [integer] maximum size of the pathways to compare.
  # :param lineage: [string] code for the lineage to analyse.
  # :param fdr.threshold: [boolean] FDR threshold to filter the sets. The default
  # value is 0.25 to be consistent with GSEA, whose default value is 0.25.
  # :param n_results: [integer] n results to select with the highest 100 * DE / N.
  # If it is null, all the results are returned.
  # :param no.title: [boolean] code to control the title in the printed image.
  # :param no.legend: [boolean] code to control the legend in the printed image.
  # :return: [data.frame] result of the test for over representation.
  
  # Compute the test
  kegga_response <- kegga(entrezgene_ids, species = species, FDR = fdr.threshold)
  analysis_results <- kegga_response[kegga_response$N <= upper_limit, ] 
  analysis_results <- analysis_results[analysis_results$N >= lower_limit, ] 
  results <- analysis_results %>% transmute(p.value = P.DE,
                                            term = Pathway,
                                            term.id = rownames(analysis_results),
                                            adj.p.value = p.adjust(P.DE, method = "BH"),
                                            ratio = 100 * DE / N, de = DE, n = N)
  results <- results[order(-results$ratio), ]
  # The threshold for the analysis correction is set to 0.05. This is an ORA 
  # method like TFRA. The same adjusted p-value is used in both.
  results <- results[results$adj.p.value <= 0.05, ]
  main_results <- if (is.null(n_results)) {results} else {head(results, n_results)}
  main_results$term <- StrCap(main_results$term, method="first")
  height <- findMid(main_results)
  plot <- main_results %>%
    ggplot(aes(y = reorder(term, ratio), x = de, fill = ratio)) + 
    labs(fill = "Annotation\npercentage", y = "KEGG term description",
         x = "Annotated genes of the set") +
    geom_bar(stat = "identity") + 
    geom_hline(yintercept = height , colour = "red", linetype = "dashed") +
    theme_bw() + 
    theme(text = element_text(size=14),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14))
  # Layout parameters
  if (!no.title) {plot + ggtitle(paste("KEGGA: Pathways ", "(", lineage, ")", sep = ""))}
  if (no.legend) {plot <- plot + guides(color = F, fill = F)}
  # Indications on terminal
  print("KEGG Analysis")
  print(paste("Lineage:", lineage, sep = " "))
  print("Plotted terms:")
  for (term in main_results$term) {
    print(term)
  }
  # Print final plot
  print(plot)
  # Return the summary table
  summary <- data.frame(term = sapply(main_results$term.id, function(x) {str_split(x, ":")[[1]][2]}),
                        description = main_results$term,
                        annotation = main_results$ratio,
                        size = main_results$n)
  rownames(summary) <- NULL
  summary <- summary[, c(1, 2, 4, 3)]
  return(summary)
}


# ------------------------------------------------------------------------------
# Read data
# ------------------------------------------------------------------------------

setwd('/home/mario/Projects/az_project/')

# Load data from the linear models
BFU_res <- fread('./data/BFUE_res.tsv')
GM_res <- fread('./data/GM_res.tsv')
MK_res <- fread('./data/Mk_res.tsv')

# Load gene expression data
expr_data = read.csv('./data/Normalised_bone_marrow_data_Ensembl_18122020.csv')
expr_data$probe_id <- NULL
rownames(expr_data) <- expr_data$Ensembl
expr_data$Ensembl <- NULL

# Translate to entrezgene_ids
mart <- useEnsembl(biomart='ensembl', dataset='hsapiens_gene_ensembl') # Some times the Ensembl site is unresponsive.
translation <-  getBM(attributes = c('ensembl_gene_id', 'entrezgene_id'),
                      filters = 'ensembl_gene_id',
                      values = rownames(expr_data),
                      mart=mart)

# Translate to HGNC symbols
species <- "Hs"
BFU_go <- merge(x = BFU_res, y = translation, merge = "ensembl_gene_id")
BFU_go <- BFU_go[order(BFU_go$P.Value), ]
GM_go <- merge(x = GM_res, y = translation, merge = "ensembl_gene_id")
GM_go <- GM_go[order(GM_go$P.Value), ]
MK_go <- merge(x = MK_res, y = translation, merge = "ensembl_gene_id")
MK_go <- MK_go[order(MK_go$P.Value), ]


# ------------------------------------------------------------------------------
# Analyses
# ------------------------------------------------------------------------------

# KEGG
n_results <- NULL # Show the n_results most relevant terms
n_top_genes <- 3000 # Select the n most significant genes
folder <- "/home/mario/Projects/az_project/data/" # Folder to store the results
# Limits in gene set size to filter
lower_limit <- 10
upper_limit <- 200
# BFUE
bfue.kegg.results <- kegga_analysis(BFU_go[1:n_top_genes, ]$entrezgene_id, species, lower_limit, upper_limit, "BFUE", n_results)
filename <- "KEGG_BFUE_raw_summary.tsv"
write.table(bfue.kegg.results, paste(folder, filename, sep=""), row.names = F, col.names = F, quote = F, sep = "\t")
# GM
gm.kegg.results <- kegga_analysis(GM_go[1:n_top_genes, ]$entrezgene_id, species, lower_limit, upper_limit, "GM", n_results)
filename <- "KEGG_GM_raw_summary.tsv"
write.table(gm.kegg.results, paste(folder, filename, sep=""), row.names = F, col.names = F, quote = F, sep = "\t")
# MK
mk.kegg.results <- kegga_analysis(MK_go[1:n_top_genes, ]$entrezgene_id, species, lower_limit, upper_limit, "MK", n_results)
filename <- "KEGG_MK_raw_summary.tsv"
write.table(mk.kegg.results, paste(folder, filename, sep=""), row.names = F, col.names = F, quote = F, sep = "\t")

