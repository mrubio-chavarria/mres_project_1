
# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------
library(biomaRt)
library(dplyr)
library(qusage)
library(gdata)


# ------------------------------------------------------------------------------
# Author: Piotr Grabowski
# ------------------------------------------------------------------------------

setwd('/home/mario/Projects/az_project')

# Load expression data
bm_data = read.csv('/home/mario/Projects/az_project/data/Normalised_bone_marrow_data_Ensembl_18122020.csv', sep=',', header=T)
bm_data$probe_id = NULL
row.names(bm_data) = bm_data$Ensembl
bm_data$Ensembl = NULL

# Translate to Entrez genes for GSEA
mart = useEnsembl(biomart='ensembl', dataset = 'hsapiens_gene_ensembl')
translation = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                    filters = 'ensembl_gene_id',
                    values = row.names(bm_data),
                    mart = mart)
translation = translation[!duplicated(translation$ensembl_gene_id),]
translation = translation[!duplicated(translation$hgnc_symbol),]

# Remap to hgnc
bm_data_hgnc = bm_data
bm_data_hgnc = bm_data_hgnc[row.names(bm_data_hgnc) %in% translation$ensembl_gene_id, ]
row.names(bm_data_hgnc) = translation$hgnc_symbol[match(row.names(bm_data_hgnc),translation$ensembl_gene_id)]
row.names(bm_data_hgnc) = toupper(row.names(bm_data_hgnc))

# Load metadata
bm_meta = read.csv('/home/mario/Projects/az_project/data/A4240_sample_annotation_and_QC.csv', sep=',', header=T)
bm_meta$Experimental_group = paste(bm_meta$Additional.sample.description.1, bm_meta$Additional.sample.description.3)
bm_meta$Experimental_group = gsub(' ', '_', bm_meta$Experimental_group)
bm_meta$Experimental_group = gsub('-', '', bm_meta$Experimental_group)
# ------------------------------------------------------------------------------
# Adapted by Mario
temp <- bm_meta %>% filter(Comments != 'FAILED',
                           Comments != 'Marginal') %>% select("Unique.sample.ID")
bm_meta = bm_meta[bm_meta$Unique.sample.ID %in% temp$Unique.sample.ID, ]
# bm_meta = bm_meta[bm_meta$Unique.sample.ID %in% colnames(temp), ]
# ------------------------------------------------------------------------------

# Required input format for phenotype (description of what is what)
# More help here: http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#CLS:_Categorical_.28e.g_tumor_vs_normal.29_class_file_format_.28.2A.cls.29

# 1st line: X Y Z, where X is number of samples, Y is total num of classes, Z is always 1 (??)
# 2nd line: # A B ..., where A and B will be names of classes
# 3rd line: 0 0 1 1..., where 0 and 1 representing classes, could be names e.g. 'aml'

################################################################################
# 1 Create annotation + subset data for BFUE, only day 10 and day 0 initially
################################################################################

bm_meta_bfue = bm_meta[bm_meta$Additional.sample.description.3 == 'BFU-E', ]
bm_meta_bfue = bm_meta_bfue[bm_meta_bfue$Additional.sample.description.1 %in% c('Day 0', 'Day 10'),]
bm_meta_bfue = bm_meta_bfue[order(bm_meta_bfue$Experimental_group),]

# subset expression data and order in the same way (important!)
bfue_ordered = bm_data_hgnc[ ,bm_meta_bfue$Unique.sample.ID]

# export a phenotype file
create_pheno_description = function(metadata){
  num_samples = nrow(metadata)
  num_classes = length(unique(metadata$Experimental_group))
  
  first_line = as.character(c(num_samples, num_classes, 1))
  second_line = paste('#', paste(unique(metadata$Experimental_group), collapse=' '), collapse='')
  third_line = metadata$Experimental_group
  
  return( list(first_line, second_line, third_line))
}

# get the list with lines needed by GSEA as sample annotation
BFUE_pheno = create_pheno_description(bm_meta_bfue)

# Write out the CLS file with sample description
bfue_pheno_file = file('./gsea_input_files/bfue_pheno.cls')
open(bfue_pheno_file, 'w')
for(line in BFUE_pheno){
  line = paste(line, collapse=' ')
  writeLines(line, con = bfue_pheno_file)
}
close(bfue_pheno_file)

# Save microarray data in tab-separated txt format for GSEA
# It has to have two ID columns
bfue_ordered$NAME = row.names(bfue_ordered)
row.names(bfue_ordered) = NULL
bfue_ordered$DESCRIPTION = NA
bfue_ordered = bfue_ordered[,c('NAME','DESCRIPTION',setdiff(names(bfue_ordered),c('NAME', 'DESCRIPTION')))] 
bfue_ordered = bfue_ordered[bfue_ordered$NAME != '',]

# Write out
write.table(bfue_ordered, './gsea_input_files/bfue_gsea_input.txt',
            sep='\t', quote=F, row.names = F)

# TODO: 
# Mario: repeat for two other lineages (re-use the code from above)
# Use function create_pheno_description() I created above and repeat the same formatting I did for BFU-E

################################################################################
# 2 Create annotation + subset data for MK, only day 14 and day 0 initially
################################################################################

bm_meta_mk = bm_meta[bm_meta$Additional.sample.description.3 == 'Mk', ]
bm_meta_mk = bm_meta_mk[bm_meta_mk$Additional.sample.description.1 %in% c('Day 0', 'Day 14'),]
bm_meta_mk = bm_meta_mk[order(bm_meta_mk$Experimental_group),]

# subset expression data and order in the same way (important!)
mk_ordered = bm_data_hgnc[ ,bm_meta_mk$Unique.sample.ID]

# get the list with lines needed by GSEA as sample annotation
MK_pheno = create_pheno_description(bm_meta_mk)

# Write out the CLS file with sample description
mk_pheno_file = file('./gsea_input_files/mk_pheno.cls')
open(mk_pheno_file, 'w')
for(line in MK_pheno){
  line = paste(line, collapse=' ')
  writeLines(line, con = mk_pheno_file)
}
close(mk_pheno_file)

# Save microarray data in tab-separated txt format for GSEA
# It has to have two ID columns
mk_ordered$NAME = row.names(mk_ordered)
row.names(mk_ordered) = NULL
mk_ordered$DESCRIPTION = NA
mk_ordered = mk_ordered[,c('NAME','DESCRIPTION',setdiff(names(mk_ordered),c('NAME', 'DESCRIPTION')))] 
mk_ordered = mk_ordered[mk_ordered$NAME != '',]

# Write out
write.table(mk_ordered, './gsea_input_files/mk_gsea_input.txt',
            sep='\t', quote=F, row.names = F)

################################################################################
# 3 Create annotation + subset data for GM, only day 14 and day 0 initially
################################################################################

bm_meta_gm = bm_meta[bm_meta$Additional.sample.description.3 == 'GM', ]
bm_meta_gm = bm_meta_gm[bm_meta_gm$Additional.sample.description.1 %in% c('Day 0', 'Day 14'),]
bm_meta_gm = bm_meta_gm[order(bm_meta_gm$Experimental_group),]

# subset expression data and order in the same way (important!)
gm_ordered = bm_data_hgnc[ ,bm_meta_gm$Unique.sample.ID]

# get the list with lines needed by GSEA as sample annotation
GM_pheno = create_pheno_description(bm_meta_gm)

# Write out the CLS file with sample description
gm_pheno_file = file('./gsea_input_files/gm_pheno.cls')
open(gm_pheno_file, 'w')
for(line in GM_pheno){
  line = paste(line, collapse=' ')
  writeLines(line, con = gm_pheno_file)
}
close(gm_pheno_file)

# Save microarray data in tab-separated txt format for GSEA
# It has to have two ID columns
gm_ordered$NAME = row.names(gm_ordered)
row.names(gm_ordered) = NULL
gm_ordered$DESCRIPTION = NA
gm_ordered = gm_ordered[,c('NAME','DESCRIPTION',setdiff(names(gm_ordered),c('NAME', 'DESCRIPTION')))] 
gm_ordered = gm_ordered[gm_ordered$NAME != '',]

# Write out
write.table(gm_ordered, './gsea_input_files/gm_gsea_input.txt',
            sep='\t', quote=F, row.names = F)


# ------------------------------------------------------------------------------
# Author: Mario Rubio
# ------------------------------------------------------------------------------

# Read data
file <- './data/c8.all.v7.2.symbols.gmt'
old_sets <- read.gmt(file)

# Delete gene sets that are not of interest
sets <- old_sets
gene_sets <- names(sets)
names(sets) <- NULL
rows <- c()

# Build the matrix with the gene sets of interest
for (i in 1:length(sets)) {
  name <- gene_sets[i]
  if (startsWith(name, 'HAY_BONE_MARROW') | startsWith(name, 'ZHENG_CORD_BLOOD')) {
    print(i)
    rows <- c(rows, paste(gene_sets[i], "https://www.google.com", paste(unlist(sets[i]), collapse = "\t"), sep = "\t"))
    
  }
}

# Write the new file
fileConn<-file("./gsea_input_files/custom_c8.gmt")
writeLines(rows, fileConn)
close(fileConn)

