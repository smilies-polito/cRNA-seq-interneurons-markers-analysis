# ScRNA-seq mammalian cortex interneurons markers meta-analysis
This repository contains a set of R scripts to perform a meta-analysis on the reliability of selected biomarkers in discriminating among different interneurons subtypes, and all the data and results from the related study.

## "Scripts" folder
In the "Scripts" folder one finds three R scripts:
- Dataset Processing.R
- Metrics Calculation Mouse.R
- Metrics Calculation Human.R

### Dataset Processing
This script allow the processing of all the datasets. One needs to change the address of the raw files to load. The scripts firstly create a Seurat Object, since  it is more flexible on the input data, and then it is transformed in a Cell Data Set(CDS, i.e., a Monocle object) through the package "Seurat Wrappers". Then, it is processed with the typical workflow.

### Metrics Calculation
It is divided in two, one for mouse datasets and one for human ones, since the genes have different notations, but they are totally equivalent.

Firstly, the script needs a processed CDS in input that must be loaded. 
The script is divided in three blocks, which represent the three levels of resolution explained in the paper. 
For each level, it creates a subset of the CDS which includes only the markers of the respective level. 

Then, one needs to rename the interneurons clusters accordingly to the level. This means that, for example, if one consider the first level where we want to analyze all the interneuron cells together, one must to rename the four cluster representing the subgroups with a univocal name. This step must be performed manually by the user. To facilitate it, it is provided the file "Processed datasets cluster numbers.txt" which includes, for each one of the datasets processed, the association of cluster and respective interneuron subgroup.

Then, the scripts calculates all the metrics, and creates a final table with the list of markers and for each one all the related results and saves it on a file.
As last thing it creates the plots:
- Specificity-Fraction Expressing
- Information Gain-Impurity Reduction
- Marker Score-Information Gain
- Marker Score-Mean Expression


## "DATA" folder
In the "DATA" folder one can find all the datasets employed in the study. In particular, one can chose between:
- the original datasets
- the processed objects

### Datasets download links.txt
In this text file one can find the direct download link to the files

### "Processed CDS" folder
In the "CDS file" folder one can find the CDS (i.e., the datasets already processed with Monocle3) for each dataset.

It is also provided the file "Processed datasets cluster numbers.txt" which include the useful list of cluster number for the interneurons subgroups of interest. These can be helpful for the user-needed passage in the metrics calculation script. 

The "Cells Plots" folder includes all the plots of the clustered cells of the processed datasets.

## "Results" folder
In the "Plots" folder there are all the final plots of the markers metrics, that one can find also in the paper.

In the "Tables" folder there are all the tables with the results of the metrics calculations. One can easily load them on R.







