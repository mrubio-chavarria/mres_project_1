
# Libraries
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
library(stringr)


# --------------------------------------------------------------
# Read data
# --------------------------------------------------------------

setwd('/home/mario/Projects/az_project/')
anno.data = fread('./data/A4240_sample_annotation_and_QC.csv')

# Remove unused samples
anno.data = anno.data[anno.data$`Additional sample description 1` !='HumanRef',]

# Keep only chips passing QC
anno.data = anno.data[anno.data$Comments =='',]

# Remove experiments unrelated to our question
anno.data = anno.data[!anno.data$`Additional sample description 3` %in% c('Glucose', 'Galactose'),]

# Remove columns which we dont need
anno.data = anno.data[,-c(5,6)]

# Rename the columns to be more informative
colnames(anno.data) = c('sample_ID', 'day', 'donor', 'lineage')

# Remove whitespace from entries
anno.data$day = gsub(' ', '', anno.data$day)
anno.data$donor = gsub(' ', '', anno.data$donor)

# Change BFU-E, otherwise it messes up contrast notation
anno.data$lineage = gsub('BFU-E', 'BFUE', anno.data$lineage)

# Create specific groups
anno.data$limma_group = paste(anno.data$lineage, anno.data$day, sep='_')

# Load expression data
expr_data = read.csv('./data/Normalised_bone_marrow_data_Ensembl_18122020.csv')
expr_data$probe_id = NULL
rownames(expr_data) <- expr_data$Ensembl
expr_data$Ensembl <- NULL


# --------------------------------------------------------------
# Compute linear models
# --------------------------------------------------------------

# Create a design matrix for
anno.data = anno.data[anno.data$sample_ID %in% colnames(expr_data),]
main_design = model.matrix(~0 + limma_group, data = anno.data)

# Contrasts to run first
# Main comparison last day vs day 0 for each lineage
# BFU-E_Day10-BFU-E_Day0
# Mk_Day14-Mk_Day0
# GM_Day14-GM_Day0
main_fit = lmFit(expr_data, main_design)

BFU_fit = contrasts.fit(main_fit, contrasts=makeContrasts(limma_groupBFUE_Day10-limma_groupBFUE_Day0, levels = main_design) )
MK_fit = contrasts.fit(main_fit, contrasts=makeContrasts(limma_groupMk_Day14-limma_groupMk_Day0, levels = main_design))
GM_fit = contrasts.fit(main_fit, contrasts=makeContrasts(limma_groupGM_Day14-limma_groupGM_Day0, levels = main_design))

BFU_fit = eBayes(BFU_fit, trend=T, robust=T)
MK_fit = eBayes(MK_fit, trend=T, robust=T)
GM_fit = eBayes(GM_fit, trend=T, robust=T)

BFU_res = topTable(BFU_fit, sort.by = 'p', number=Inf)
MK_res = topTable(MK_fit, sort.by = 'p', number=Inf)
GM_res = topTable(GM_fit, sort.by = 'p', number=Inf)

BFU_res$ensembl_gene_id <- rownames(BFU_res)
rownames(BFU_res) <- NULL
MK_res$ensembl_gene_id <- rownames(MK_res)
rownames(MK_res) <- NULL
GM_res$ensembl_gene_id <- rownames(GM_res)
rownames(GM_res) <- NULL

# --------------------------------------------------------------
# Store data
# --------------------------------------------------------------

# Modified by Mario
anno.data$day = as.integer(str_replace(anno.data$day, 'Day', ''))
anno.data$donor = as.integer(str_replace(anno.data$donor, 'Donor', ''))

# Store annotation data
write.table(anno.data, file = "./data/processed_annotation.tsv", append = F, quote = F, sep="\t", row.names = F)

# Store model results
write.table(BFU_res, file = "./data/BFUE_res.tsv", append = F, quote = F, sep="\t", row.names = F)
write.table(GM_res, file = "./data/GM_res.tsv", append = F, quote = F, sep="\t", row.names = F)
write.table(MK_res, file = "./data/Mk_res.tsv", append = F, quote = F, sep="\t", row.names = F)

