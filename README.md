# umi-normalization

### Pre-requisites


- `FIt-SNE` by George C. Linderman from https://github.com/KlugerLab/FIt-SNE (version used in the paper: `1.2.1`, https://github.com/KlugerLab/FIt-SNE/releases/tag/v1.2.1)
- `FFTW` from http://www.fftw.org (version used in the paper: `3.3.8`)
- `rpy2` from https://rpy2.github.io/ (version used in the paper: `3.2.6`)
- R package `glmpca` by Will Townes from https://github.com/willtownes/glmpca (version used in the paper: `0.2.0`, https://github.com/willtownes/glmpca/releases/tag/v0.2.0)
- `glmpca-py` by Will Townes from https://github.com/willtownes/glmpca-py/
(version used in the paper: https://github.com/willtownes/glmpca-py/tree/a6fc417b08ab5bc21d8ac9e197f4f5518d093385)
- `rna-seq-tsne` by Dmitry Kobak from https://github.com/berenslab/rna-seq-tsne (version used in the paper: https://github.com/berenslab/rna-seq-tsne/tree/21e3601782d37dd3f0c8e02ed9f239b005c4100f)

### Download instructions for presented datasets

##### 33k PBMC dataset

###### Counts & Annotations
 - visit https://support.10xgenomics.com/single-cell-gene-expression/datasets 
 - look for '33k PBMC from a healty donor' under "Chromium Demonstration (v1 Chemistry)"
 - provide contact details to proceed to downloads
 - download 'Gene / cell matrix (filtered)' (79.23 MB) 
 - extract files `genes.tsv` and `matrix.mtx` to `umi-normalization/datasets/33k_pbmc/`
 - download 'Clustering analysis' (23.81 MB) from the same website
 - extract folder `analysis` to `umi-normalization/datasets/33k_pbmc/` as well




##### 10x control / Svensson 2017
 - visit https://figshare.com/articles/svensson_chromium_control_h5ad/7860092 
 - download `svensson_chromium_control.h5ad` (18.82 MB)
 - save at `umi-normalization/datasets/10x/`

##### inDrop control / Klein 2015
 - visit https://www.ncbi.nlm.nih.gov/geo/ 
 - search for `GSE65525`
 - download the `*.csv.bz2` file for the sample `GSM1599501` (human K562 pure RNA control, 953 samples, 5.1 MB))
 - extract file `GSM1599501_K562_pure_RNA.csv` to `umi-normalization/datasets/indrop/`

##### MicrowellSeq control / Han 2018
 - visit https://www.ncbi.nlm.nih.gov/geo/ 
 - search for `GSE108097`
 - search for sample `GSM2906413` and download `GSM2906413_EmbryonicStemCell_dge.txt.gz` (EmbryonicStemCell.E14, 7.9 MB)
 - save to `umi-normalization/datasets/microwellseq/`

##### Retina: All cell classes/ Macosko 2015

###### Counts
 - visit https://www.ncbi.nlm.nih.gov/geo/
 - search for `GSE63472`
 - download `GSE63472_P14Retina_merged_digital_expression.txt.gz` (50.7 MB)
 - extract to `GSE63472_P14Retina_merged_digital_expression.txt`
 - save to `umi-normalization/datasets/retina/macosko_all`
###### Cluster annotations
 - download from http://mccarrolllab.org/wp-content/uploads/2015/05/retina_clusteridentities.txt
 - save to `umi-normalization/datasets/retina/macosko_all/`


##### Retina: Bipolar cells / Shekhar 2016

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


##### Retina: Ganglion cells / Tran 2019

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

