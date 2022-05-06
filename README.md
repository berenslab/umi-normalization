# Analytic Pearson residuals for normalization of single-cell RNA-seq UMI data
###### Jan Lause, Philipp Berens & Dmitry Kobak


### How to use this repository

Version `3.0` of this repository contains the code to reproduce the analysis presented in our *Genome Biology* paper on UMI data normalization ([Lause, Berens & Kobak, 2021](https://doi.org/10.1186/s13059-021-02451-7)) and the corresponding [preprint (v3)](https://doi.org/10.1101/2020.12.01.405886 ). The code used for versions v1 and v2 of the paper is available under the tags `1.0` and `2.0` in this repository.

To start, follow these steps:

 - install the required software listed below
 - clone this repository to your system
 - go to `tools.py` and adapt the three import paths as needed
 - follow the dataset download instructions below
 
Then, you can step through our full analysis by simply following the sequence of the notebooks. If you want to reproduce only parts of our analysis, there are six independent analysis pipelines that you can run individually:

 - Reproduction of the NB regression model by Hafemeister & Satija (2019) and investigation of alternative models (Notebookes `01` & `02`, producing Figure 1 from our paper)
 - Estimation of technical overdispersion from negative control datasets (Notebooks `01` & `03`, producing Figure S1)
 - Benchmarking normalization by Analytical Pearson residuals vs. GLM-PCA vs. standard methods:
     - on the 33k PBMC dataset (Notebooks `01`, `041`, `042`, `05`, producing Figures 2, S2, S4, S5, and additional figures)
     - on different retinal datasets (Notebooks `06`, `07`, `081`, producing Figures 3, S3, and additional figures)
     - on the ground-truth dataset created from FACS-sorted PBMCs (Notebook `101`, `102`, producing Figures 5 and S7)
 - Analysis of the 2-million cell mouse organogenesis dataset (Notebook `091`, producing Figures 4 and S6, and additional figures)
 - Comparison to [Sanity](https://github.com/jmbreda/Sanity) (Notebooks `06`, `07`, `081` and `082` for retina datasets and `091` and `092` for the organogenesis dataset, producing additional figures). These pipelines will require you to run Sanity from the command line; see notebooks `082` and `092` for instructions.

Note that `041` and `101` are R notebooks, the remaining are Python notebooks.

Each of the analyses will first preprocess and filter the datasets. Next, computationally expensive tasks are done (NB regression fits, GLM-PCA, t-SNE, simulations of negative control data, ..) and the results are saved as files. For some analyses, this is done in separate notebooks. Finally, the results files are loaded for plotting (again in separate notebooks for some analyses).

We recommend to run the code on a powerful machine with at least 250 GB RAM.

For questions or feedback, feel free to use the issue system or email us.

### Pre-requisites

We used the following software environments:

##### Python

- Python         `3.8.0`
- IPython        `7.21.0`
- `numpy`        `1.20.1`
- `pandas`       `1.2.0`
- `scipy`        `1.6.0`
- `seaborn`      `0.11.1`
- `matplotlib`   `3.3.3`
- `statsmodels`  `0.12.2`
- `sklearn`      `0.24.0`
- `anndata`      `0.7.5` from https://anndata.readthedocs.io
- `scanpy`       `1.7.1` from https://scanpy.readthedocs.io
- `FIt-SNE` by George C. Linderman from https://github.com/KlugerLab/FIt-SNE, using this version https://github.com/KlugerLab/FIt-SNE/tree/47ff14f1defc1ff3a8065c8b7baaf45c33e7e0b2
- `FFTW`         `3.3.8` from http://www.fftw.org
- `rpy2`          `3.4.2` from https://rpy2.github.io/ 
- `glmpca-py` by Will Townes from https://github.com/willtownes/glmpca-py/, using this version:
 https://github.com/willtownes/glmpca-py/tree/a6fc417b08ab5bc21d8ac9e197f4f5518d093385
- `rna-seq-tsne` by Dmitry Kobak from https://github.com/berenslab/rna-seq-tsne, using this version:  https://github.com/berenslab/rna-seq-tsne/tree/21e3601782d37dd3f0c8e02ed9f239b005c4100f

##### R

- R                       `3.6.3 (2020-02-29)`
- `glmpca`                `0.2.0` by Will Townes from https://github.com/willtownes/glmpca, https://github.com/willtownes/glmpca/releases/tag/v0.2.0. If not present in the local R version, our code will try to install the most recent version from CRAN.
- `SingleCellExperiment`  `1.8.0`
- `MASS`                  `7.3-53.1`
- `sctransform`           `0.3.2` by Christoph Hafemeister (https://github.com/ChristophH/sctransform, 

The full R environment used was

```
attached base packages:
parallel  stats4    stats     graphics  grDevices utils     datasets     methods   base     

other attached packages:
MASS_7.3-53.1               sctransform_0.3.2          SingleCellExperiment_1.8.0  SummarizedExperiment_1.16.1   
DelayedArray_0.12.3         BiocParallel_1.20.1        matrixStats_0.58.0          Biobase_2.46.0             
GenomicRanges_1.38.0        GenomeInfoDb_1.22.1        IRanges_2.20.2              S4Vectors_0.24.4           
BiocGenerics_0.32.0         glmpca_0.2.0               

loaded via a namespace (and not attached):
tidyselect_1.1.0       listenv_0.8.0          purrr_0.3.4           reshape2_1.4.4         lattice_0.20-41        colorspace_2.0-0      
vctrs_0.3.7            generics_0.1.0         utf8_1.2.1            rlang_0.4.10           pillar_1.6.0           glue_1.4.2            
DBI_1.1.1              GenomeInfoDbData_1.2.2 lifecycle_1.0.0       plyr_1.8.6             stringr_1.4.0          zlibbioc_1.32.0       
munsell_0.5.0          gtable_0.3.0           future_1.21.0         codetools_0.2-18       fansi_0.4.2            Rcpp_1.0.6            
scales_1.1.1           XVector_0.26.0         parallelly_1.24.0     gridExtra_2.3          ggplot2_3.3.3          digest_0.6.27         
stringi_1.5.3          dplyr_1.0.5            grid_3.6.3            tools_3.6.3            bitops_1.0-6           magrittr_2.0.1        
RCurl_1.98-1.3         tibble_3.1.1           crayon_1.4.1          future.apply_1.7.0     pkgconfig_2.0.3        ellipsis_0.3.1        
Matrix_1.3-2           assertthat_0.2.1       R6_2.5.0              globals_0.14.0         compiler_3.6.3  
```



### Download instructions for presented datasets

All accession numbers can also be found in Table S2 of our paper.

##### 33k PBMC dataset

###### Counts & Annotations
 - visit https://support.10xgenomics.com/single-cell-gene-expression/datasets 
 - look for '33k PBMC from a healty donor' under "Chromium Demonstration (v1 Chemistry)"
 - provide contact details to proceed to downloads
 - download 'Gene / cell matrix (filtered)' (79.23 MB) 
 - extract files `genes.tsv` and `matrix.mtx` to `umi-normalization/datasets/33k_pbmc/`
 - download 'Clustering analysis' (23.81 MB) from the same website
 - extract folder `analysis` to `umi-normalization/datasets/33k_pbmc/` as well
 
##### 10x control / Svensson et al. 2017
 - visit https://figshare.com/articles/svensson_chromium_control_h5ad/7860092 
 - download `svensson_chromium_control.h5ad` (18.82 MB)
 - save at `umi-normalization/datasets/10x/`

##### inDrop control / Klein et al. 2015
 - visit https://www.ncbi.nlm.nih.gov/geo/ 
 - search for `GSE65525`
 - download the `*.csv.bz2` file for the sample `GSM1599501` (human K562 pure RNA control, 953 samples, 5.1 MB))
 - extract file `GSM1599501_K562_pure_RNA.csv` to `umi-normalization/datasets/indrop/`

##### MicrowellSeq control / Han et al. 2018
 - visit https://www.ncbi.nlm.nih.gov/geo/ 
 - search for `GSE108097`
 - search for sample `GSM2906413` and download `GSM2906413_EmbryonicStemCell_dge.txt.gz` (EmbryonicStemCell.E14, 7.9 MB)
 - save to `umi-normalization/datasets/microwellseq/`

##### Retina: All cell classes/ Macosko et al. 2015

###### Counts
 - visit https://www.ncbi.nlm.nih.gov/geo/
 - search for `GSE63472`
 - download `GSE63472_P14Retina_merged_digital_expression.txt.gz` (50.7 MB)
 - extract to `GSE63472_P14Retina_merged_digital_expression.txt`
 - save to `umi-normalization/datasets/retina/macosko_all`
###### Cluster annotations
 - download from http://mccarrolllab.org/wp-content/uploads/2015/05/retina_clusteridentities.txt
 - save to `umi-normalization/datasets/retina/macosko_all/`


##### Retina: Bipolar cells / Shekhar et al. 2016

###### Counts
 - visit https://www.ncbi.nlm.nih.gov/geo/
 - search for `GSE81904`
 - download `GSE81904_BipolarUMICounts_Cell2016.txt.gz` (42.9 MB)
 - save to `umi-normalization/datasets/retina/shekhar_bipolar/`
###### Cluster annotations
 - visit https://singlecell.broadinstitute.org/single_cell/study/SCP3/retinal-bipolar-neuron-drop-seq
 - click `Download`
 - register with google account
 - download `clust_retinal_bipolar.txt` (1.5 MB)
 - save to `umi-normalization/datasets/retina/shekhar_bipolar/`


##### Retina: Ganglion cells / Tran et al. 2019

###### Raw counts
 - visit https://www.ncbi.nlm.nih.gov/geo/
 - search for `GSE133382`
 - download `GSE133382_AtlasRGCs_CountMatrix.csv.gz` (129.3 MB)
 - extract to `GSE133382_AtlasRGCs_CountMatrix.csv`
 - save to `umi-normalization/datasets/retina/tran_ganglion/`

###### Annotations and original gene selection
 - visit https://singlecell.broadinstitute.org/single_cell/study/SCP509/mouse-retinal-ganglion-cell-adult-atlas-and-optic-nerve-crush-time-series
 - sign in with google account
 - download `RGC_Atlas.csv` (1.05 GB) and `RGC_Atlas_coordinates.txt` (927 KB)
 - save to `umi-normalization/datasets/retina/tran_ganglion/`
  
 
##### 2-million cells: Mouse Organogenesis / Cao et al. 2019

###### Raw counts and annotations
 - visit https://oncoscape.v3.sttrcancer.org/atlas.gs.washington.edu.mouse.rna/downloads
 - download `gene_count.txt` (18 GB),  `gene_annotation.csv` (1.1 MB) and `cell_annotation.csv` (828 MB)
 - save to `umi-normalization/datasets/cao`


##### FACS-sorted PBMC cells / Zheng et al. (2017) and Duò et al (2018)
 - make sure you have an R environment with the Bioconductor package `SingleCellExperiment` installed
 - visit http://imlspenticton.uzh.ch/robinson_lab/DuoClustering2018/ or https://github.com/markrobinsonuzh/scRNAseq_clustering_comparison
 - download `DuoClustering2018.tar.gz` (4.93 GB)
 - extract into `umi-normalization/datasets`
 - make sure that `sce_full_Zhengmix8eq.rds` exists at `umi-normalization/datasets/DuoClustering2018/sce_full/`
