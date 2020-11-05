# umi-normalization

### Pre-requisites

- `rna-seq-tsne` from https://github.com/berenslab/rna-seq-tsne
- `FIt-SNE` from https://github.com/KlugerLab/FIt-SNE (requires `FFTW`)
- `rpy2`


### Download instructions for presented datasets

##### 33k PBMC dataset
 - visit https://support.10xgenomics.com/single-cell-gene-expression/datasets 
 - look for '33k PBMC from a healty donor' under "Chromium Demonstration (v1 Chemistry)"
 - provide contact details to proceed to downloads
 - download 'Gene / cell matrix (filtered)' (79.23 MB)
 - extract files `genes.tsv` and `matrix.mtx` to `umi-normalization/datasets/33k_pbmc`

##### 10x control / Svensson 2017
 - visit https://figshare.com/articles/svensson_chromium_control_h5ad/7860092 
 - download `svensson_chromium_control.h5ad` (18.82 MB)
 - save at `umi-normalization/datasets/10x`

##### inDrop control / Klein 2015
 - visit https://www.ncbi.nlm.nih.gov/geo/ 
 - search for `GSE65525`
 - download the `*.csv.bz2` file for the sample `GSM1599501` (human K562 pure RNA control, 953 samples, 5.1 MB))
 - extract file `GSM1599501_K562_pure_RNA.csv` to `umi-normalization/datasets/indrop`

##### MicrowellSeq control / Han 2018
 - visit https://www.ncbi.nlm.nih.gov/geo/ 
 - search for `GSE108097`
 - search for sample `GSM2906413` and download `GSM2906413_EmbryonicStemCell_dge.txt.gz` (EmbryonicStemCell.E14, 7.9 MB)
 - save to `umi-normalization/datasets/microwellseq`

##### Retina: All cell classes/ Macosko 2015

###### Counts
 - visit https://www.ncbi.nlm.nih.gov/geo/
 - search for `GSE63472`
 - download `GSE63472_P14Retina_merged_digital_expression.txt.gz` (50.7 MB)
 - extract to `GSE63472_P14Retina_merged_digital_expression.txt`
 - save to `umi-normalization/datasets/retina/macosko_all`
###### Cluster annotations
 - download from http://mccarrolllab.org/wp-content/uploads/2015/05/retina_clusteridentities.txt
 - save to `umi-normalization/datasets/retina/macosko_all`


##### Retina: Bipolar cells / Shekhar 2016

###### Counts
 - visit https://www.ncbi.nlm.nih.gov/geo/
 - search for `GSE81904`
 - download `GSE81904_BipolarUMICounts_Cell2016.txt.gz` (42.9 MB)
 - save to `umi-normalization/datasets/retina/shekhar_bipolar`
###### Cluster annotations
 - visit https://singlecell.broadinstitute.org/single_cell/study/SCP3/retinal-bipolar-neuron-drop-seq
 - click `Download`
 - register with google account
 - download `clust_retinal_bipolar.txt` (1.5 MB)
 - save to `umi-normalization/datasets/retina/shekhar_bipolar`


##### Retina: Ganglion cells / Tran 2019

###### Counts & annotations
 - we obtained `RGC_atlas.rds` from personal communication with the authors.

