# Overview

Features drawn from the topology of the human protein-protein interaction (PPI) network have previously been used to successfully predict essential genes in a range of organisms including humans. Whilst these predictions have been valuable these studies have generally focused on genes commonly which are essential in almost all cells and as such offer limited opportunity for drug targeting and tailored cancer therapy.
Here we explore how well this approach can be extended to predict genes dependencies in specific cancer cell lines that have been acquired following disruptions in the PPI network due to mutations and changes in gene expression.  The overriding motivation for predicting cancer cell line specific gene dependencies is to provide a low-cost approach to identifying personalised cancer drug targets without the cost of exhaustive loss of function screening.

Using recent experimental cancer cell line gene dependency data along with PPI network and genomic alteration data we have built a range of machine learning classification models to predict cell line specific acquired gene dependencies. Genetic alterations found in each individual cell line were modeled using a novel approach. This involved removing protein nodes to reflect loss of function mutations and changing the weights of edges in each PPI to reflect gain of function mutations and gene expression changes.

We have found that PPI networks can be used to successfully classify human cell line specific gene dependencies within individual cell lines and between cell lines, even across tissue types with AUC ROC scores of between 0.75 and 0.85.  Our perturbed PPI network models further improved prediction power compared to the base PPI model and are shown to be more sensitive to genes on which the cell becomes dependent as a result of other changes: an improvement which offers opportunities for personalised therapy.  

# Setting up source data

The Dependant source code is dependnent on data from a range of sources such as the orignal depmap data and string ptotein interacton data. Below is a list of the files you will need along with links and the directory structure they take.

* Note, links below might need trailing ] characters removed

Download id map (For matching gene and protein ids / names etc) via 
[https://bitbucket.org/bioinformatics_lab_sussex/dependant2/downloads/human_id_map.txt] and save to
[data_dir]/supporting_files/human_id_map.txt

Download dependency data via 
[https://ndownloader.figshare.com/files/12704096]
Or find this file manually:

* Visit: https://depmap.org/portal/download/

* Select: DepMap Public 18Q3

* Download: gene_dependency.csv

and add data to [data_dir]/dependency/gene_dependency.csv

Download weights file via [https://bitbucket.org/bioinformatics_lab_sussex/dependant2/downloads/all_preds_raw_0.csv]
[https://bitbucket.org/bioinformatics_lab_sussex/dependant2/downloads/base_weights.csv]
and save both to [data_dir]/weights/...
Alternative you can use the python scripts availbale in the python directory to create custom weights files. Documentation to do this is inlcuded in the jupyter notebook files

Protein interaction data via string
https://stringdb-static.org/download/protein.links.detailed.v11.0.txt.gz unzipped and 
to be saved in [data_dir]/supporting_files/9606.protein.links.detailed.v11.txt

# Running the programme

Once these files are in place you should be able to run the mutlipipe.R file to run all preprocessing, feature generation, model training and computational validation.

Please note that this script was originally produced to run each script seperately and the pipeline has been added for ease of use. If you have any issues tryin running each script seperately for better error messages.

Once the multipipe.R file as been run successfully you should be able to run prediction.R to predict your own cell lines as required.


