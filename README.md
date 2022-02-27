## News

Our publication introducing PlaqView has been accepted at Atherosclerosis (link coming soon!)

PlaqView for scATAC-seq is coming in 2022!

## About *PlaqView*
*PlaqView* is a standalone, interactive, and reproducible Shiny and R-based tool to explore atherosclerosis-related single-cell RNA-sequencing data. Our goal is to make these valuable data analytic tools publicly available to non-bioinformaticians. This is an open-sourced tool that is freely available and can be modified for other scRNA-seq data. Our preprint is available [HERE](https://www.biorxiv.org/content/10.1101/2020.10.27.357715v2). 

PlaqView is LIVE! [Check it out here](plaqview.com)

### Data Source
The data used in this interactive session were first published by Wirka et al. (Nature Medicine, 2019). The data were later extensively analyzed and made available here by Ma et al. (2021). Overtime, we have included additional datasets, please see the [application webpage](plaqview.com) or the 'available_dataset' file for the latest available datasets!

#### Submitting Your Data
Currently, PlaqView supports scRNA-seq type data. We are working on expanding to ATAC-seq and spatial transcriptomics. We welcome submission of all kinds of single-cell data related to cardiovascular genetics.

To learn about how data is processed or to **submit your own data**, [Click Here](plaqview.com) and click Contact Us.

For scRNA-seq data, we require the following:
1) *Seurat* object containing your count matrix and cell labels used in your publication. (Cell labels are often stored in the 'Metadata' tab in the *Seurat* obj.).

If this is not available, we will also accept a complete count matrix (or matrices) and, if available, the differential gene list used to annotate your cell type. Please note that we are unable to accept FASTQ files, please first process your raw data through cellranger before submitting.

2) a short description of your dataset, including number of cells, patients (or biological replicates), platform used (10x or other), publication DOI, types of tissue and species.

To transfer your data, we will work with you if your institution require a data transfer agreement. If you are within UVA, we prefer the use of Box.

### Reproducibility and Code Availability
*PlaqView* is designed to be open-source and publicly available and serve as a template for other scRNA-seq analysis pipelines. The code scripts are available on our [GitHub Page](https://github.com/MillerLab-CPHG/PlaqView). This application calls for .rds datafiles that are too large to host on Github, but you may request the latest files at weima@virginia.edu or via a GitHub request. We will respond to your requests within 24hrs. To Run this application successfully, place the .rds object files in a folder/subdiretory called 'data.'


### Funding and Collaborations
[UVA MSTP](https://mstp.med.virginia.edu/)
[PlaqOmics (Leducq Foundation)](https://jeanette-erdmann.jimdo.com/)

<img src="www/millerlablogo.png" alt="millerlab" width="200"/> <img src="www/PlaqOmics.png" alt="PlaqOmics" width="200"/> <img src="www/Leducq.png" alt="leducq" width="200"/><img src="www/umc.png" alt="leducq" width="200"/><img src="www/MSTPlogo.png" alt="MSTP" width="200"/>
