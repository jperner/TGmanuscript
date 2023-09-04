# TGmanuscript

##INTRODUCTION

This repository contains the R code for generating RNAseq analysis related figures and tables in the manuscript J. Beil et al. "Unaltered hepatic wound-healing response in male rats with ancestral liver injury".

##DATA

RNA sequencing data originating from the study have been deposited in NCBI GEO under the accession code: GSE229524. The GEO repository includes also the gene expression quantification and average transcript length matrices generated via Salmon (version 0.14.0). 

Both matrices need to be downloaded from the GEO repository and are the starting point for applying the script in src/. For simplicity the sample metadata has been summarized under data/.

##RUNNING THE CODE

All exploratory and differential gene expression analyses were performed in R (version 4.1.1) and the code for all steps is provided under src/. Please refer to the rmarkdown pdf for versions of individual packages used.


