# WGCNA-Analysis
WGCNA Analysis for Gene Co-expression Networks   - Correlates modules with clinical traits or experimental conditions   

This project focuses on **Weighted Gene Co-expression Network Analysis (WGCNA)** to identify gene modules associated with biological traits. WGCNA helps uncover highly connected gene clusters (modules) and their correlation with external phenotypic traits, aiding in biomarker discovery.  

## Features  
- Constructs a **gene co-expression network** from MICRO ARRAY AND RNA-seq data  
- Identifies **co-expressed gene modules**  
- Correlates modules with clinical traits or experimental conditions  
- Visualizes network structures and hub genes  

## Prerequisites  
Ensure the following R packages are installed:  
```r
install.packages("BiocManager")
BiocManager::install(c("WGCNA", "Biobase", "limma"))
