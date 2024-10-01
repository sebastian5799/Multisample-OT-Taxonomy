# Multisample-OT-Taxonomy
![taxonomy_plos](https://github.com/user-attachments/assets/3b94b0ba-b83b-46db-89bb-c56cabc2a5b2)
Single-cell RNA sequencing has signficantly advanced the study of cellular diversity and behavior, greatly enhancing our understanding of various diseases and their progression. A key part of this research is comparing cell subpopulations (clusters) across different samples. However, several challenges complicate this comparison. Often, cell clusters are unlabeled in some or all samples, the inherent variability and noise in scRNA-seq data are significant, and batch effects—systematic differences in gene expression across different batches—are common. To tackle these issues, we have developed a new system called Multisample OT Taxonomy (MOTT). MOTT automatically constructs a hierarchical taxonomy of cell clusters from multiple samples, aiding in the identification of clusters of the same cell type without manual labeling. This system integrates the advanced technique of optimal transport with relaxed marginal constraints (OT-RMC) and hierarchical clustering. Our experiments demonstrate that MOTT excels in annotating cell clusters and generating useful features for further analysis, such as classifying samples into different disease types.

The Datasets folder contains:
* the R code for loading and formatting the scRNA-seq datasets.
* Two examples od formatted datasets ready to be used in Matlab.

The OTcode folder contains:
* Functions for Optimal Transport and Optimal Transport with Marginal Constraints.  
* A function to compute the cost matrix between two samples.
* Sample formatted dataset. See Data folder for how to generate data. 
* Two Matlab notebooks that implement OT and OT-RMC functions to generate the taxonomies
* Matlab notebook to implement random forest to predict sample class using the built taxonomy
* Matlab script that implements Reference Alignment. This is used as a comparison to our method.
