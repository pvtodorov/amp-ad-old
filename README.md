# amp-ad: Analysis of AMP-AD datasets

**Authors:** Petar Todorov<sup>1</sup>, Artem Sokolov<sup>1</sup><br />
<sup>1</sup>Laboratory of Systems Pharmacology, Harvard Medical School

## Preliminaries

Nearly all of the scripts in this repository offer a command-line interface, allowing for modular design and execution of new analyses. However, to ensure that the command-line interfaces behave as expected, it is important to preconfigure the underlying `R` and `Python` environments as follows.

### Set up an R environment
To execute R scripts in this repository, you will need to pre-install several R packages. Since some of these packages come from Bioconductor, the easiest thing to do is to run the following two commands inside R:

    source("https://bioconductor.org/biocLite.R")
    biocLite( c("tidyverse","stringr","synapseClient","EnsDb.Hsapiens.v86","ensembldb") )
    
(Go get a cup of coffee...)

### Set up a python environment
Make sure `virtualenv` is installed. If not, install it via pip.
    
    pip install virtualenv

Make a new virtualenv. I named my ampad-env
    
    virtualenv ampad-env

Activate the env
    
    activate ampad-env

Clone and install the btr package

    git clone https://github.com/pvtodorov/btr.git
    cd btr
    pip install -e .


## Wrangle data
The datasets used by the analyses in this repository can be easily downloaded using command-line scripts. To download Mount Sinai Brain Bank (MSBB) dataset, run the following on the command line:

    Rscript msbb.R <path to data>

If `<path to data>` is left blank, then `/data/AMP-AD/MSBB` will be used by default.

## Obtain background predictions
We will use a random forest regression to predict the Braak score of of samples
based on the gene expression profile in the BM36 region. This region was picked
because it has multiple stages of the disease. The background will be generated
for gene sets of 10 to 1000, moving in steps of 10 randomly selected genes, where models will be
trained on the data and the R2 out-of-bag score will be recorded as the performance.
This will be repeated at least 10,000 times.

To perform this, run the `predict_background.py` script in the command line.
The script will run for the specified number of repetitions and dump a csv file
for each one. The optional `-i` tag tells the script how many times to repeat.

    btr-predict <settings.json> -i <number of repetitions>


## Score background runs
The repeated predictions will need to be aggregated into a single file and their
AUCs will be computed.

    btr-score <settings.json>


## Predict the preformance of gene sets
To predict the performance of a gene set, obtain the necessary .gmt file or
create a folder with .txt files which contain gene sets formatted as HGNC ID on
each line. These can then be supplied to the `btr-predict` as follows:

    btr-predict <settings.json> -g <path to gmt file or folder of txts>

The predictions will need to be scored

    btr-score <settings.json> -g <path to gmt file or folder of txts>

The output of `btr-score` can be supplied to `btr-stats` to compute whether a
particular gene set is producing predictions that are statistically better than
the random background distribution of predictions.

    btr-stats <settings.json> -g <path to gmt file or folder of txts>

The output will contain an identifier for the set, a description, the number of genes
in the set, the number of genes which were actually used (if there is incomplete
overlap between the dataset and the gene set not all genes are used), the
R-squared value, the p_value, and an Benjamini-Hochberg adjusted p_value.

