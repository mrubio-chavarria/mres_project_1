
# MRes in Biomedical Research
## Project 1
### Defining the transcriptional profile of in vitro differentiated haematopoietic cells.
##  Description
The analysis desribed in the report has been computed with the coded gathered in this repository.

## Files
The files follow the order described in the report. Particularly, data_processing should be executed the first script in being executed. The rest are independent from each other.
                
1. data_processing: this script should be executed the first. It reads the data, performs the linear models (t-tests), and store the results.
2. sample_analysis: this script computes and visualises the dimensionality reduction techniques described in the report: PCA, t-SNE, and UMAP.
3. expression_analysis: the study of the gene expression through MA and volcano plots. The result of the linear models.
4. kegg_analysis: gene set analysis with the KEGG terms.
5. gsea_files: this script generates the input files for GSEA. 
6. gsea_visualisation: the code to plot the evolution of the enriched genes resulting from GSEA.
7. tfra: the code to execute the transcription factor regulon analysis.
                
## Libraries
* data.table == 1.14.0
* limma == 3.46.0
*  biomaRt == 2.46.3
* stringr == 1.4.0
* tidyverse == 1.3.0
* reshape == 1.4.4
*  ggfortify == 0.4.11
* umap == 0.2.7.0
* Rtsne == 0.15
* gridExtra == 2.3
* dorothea == 1.2.1
* GSA == 1.3.1
* Biobase == 2.50.0
* DescTools == 0.99.40
* qusage == 2.24.0
* gdata == 2.18.0

This code uses R version 4.0.3.