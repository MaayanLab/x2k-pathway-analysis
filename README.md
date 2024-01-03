# X2K pathway analysis
This repository implements the eXpression2Kinases pipeline to infer upstream TFs, intermediates, and kinases from user-inputted lists of genes.

## Steps
1. [run_x2k](https://github.com/MaayanLab/x2k-pathway-analysis/blob/main/run_x2k.ipynb): Input up- and down- gene set files (containing the up- and down- gene sets for each sample as columns) to the notebook. The notebook then produces files with the predicted TFs, intermediates, and kinases for each sample. If phosphoproteomics are available, you can find the overlap of kinases enriched in both transcriptomics and phosphoproteomics.

2. [visualize_x2k_results](https://github.com/MaayanLab/x2k-pathway-analysis/blob/main/visualize_x2k_results.ipynb): After running the X2K pipeline, the pathways commonly enriched in a cohort of interest can be visualized as networks, where nodes are the inferred proteins and edges are inferred pathway connections. The weight of an edge is given by its frequency in the cohort of interest, and the size of the network can be adjusted by pruning the less common edges with smaller weights. 
