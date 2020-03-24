# global_ancestry
Global ancestry pipeline  

Panel for PGC-PTSD global ancestry estimation. Designed to run with most Illumina panels and Affymetrix UKBB arrays. Based on built in cutoffs, classifies subject as belonging to a given ancestry group. Useful for GWAS which is typically run within ancestry group.

## Installation Requirements:
R (version 3 or above)  
SNPweights https://www.hsph.harvard.edu/alkes-price/software/  
Eigensoft https://www.hsph.harvard.edu/alkes-price/software/  
PLINK 1.9 beta https://www.cog-genomics.org/plink/  
PRE-quality control genotype data. Absolutely must be pre-QC genotype data!  

## Installation
Download and unzip the repository. Maintain the directory structure of the github (i.e. scripts must remain in the scripts directory).

## Operation
Open ancestry_pipeline.sh. Modify the script where directed (where it says "User:"), giving the paths to relevant files.

## Output
the .predpc_oneweek.header file contains all relevant information (subject IDs, projections onto reference population PCs, ancestry proportion estimates, and overall ancestry group classification). The pdf file contains a plot of subjects projected onto reference population PCs.
