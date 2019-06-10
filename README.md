# Blood Atlas

## Description
Blood is the predominant source for molecular analyses in humans, both in clinical and research settings, and it is the target for many therapeutic strategies, emphasizing the need for comprehensive molecular maps of the cells constituting human blood. Here, we have performed a genome-wide transcriptomic analysis of protein-coding genes in 18 sorted blood immune cell populations to characterize the expression levels of each individual gene across the blood cell types. All data is presented in an interactive, open access Blood Atlas as part of the Human Protein Atlas (https://www.proteinatlas.org/) and is integrated with expression profiles across 37 tissues to provide a new classification of all protein-coding genes. This allows for a genome-wide exploration of the spatial expression profiles across human immune cell populations and all major human tissues and organs.


## Normalization pipeline
![image](https://github.com/human-protein-atlas/BloodAtlas/blob/master/Fig.S1.A%20schematic%20view%20of%20the%20normalization%20strategy.png)


## Code explanation
1.	The R script normalize.R inputs transcript expression data as transcript per million (TPM) values for tissues from HPA and GTEx as well as Tags Per Million values from FANTOM5 and outputs normalized expression values for each gene and source. The normalization is performed using TMM scaling, pareto scaling per gene and source and LIMMA to correct for batch effects between different sources. The full pipeline is described in the paper and as an overview in the figure S1A provided in the paper as well as in the repository.

2.	The R script function.R contains the following functions utilized by the script normalize.R to perform the normalization of all transcript expression values: ‘under_limit’, ‘impute_expression’,  ‘tmm_method_normalization’, ‘pareto_scale_method_gene’, and ‘limma_method_correction.

3.	The php script categorize.php is used to classify all Ensembl genes into the specificity and distribution categories based on the limit of detection of 1 NX and the normalized expression values (NX) for each gene and tissue. The functions used in the script depend on fetching and inserting data from and into the Human Protein Atlas internal MySQL database.
