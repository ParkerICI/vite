**--> WARNING: This package is in Beta and under active development. Things may break or change without notice. Use it at your own risk and makre sure you have a backup copy of your data <--**

Please use github issues to report bugs and for feature requests



## Installation

1. In order to install this package, you must be able to compile C++ source packages on your system:
    - Windows: Install the [Rtools](https://cran.r-project.org/bin/windows/Rtools/) package
    - MacOS:  Install the XCode software from Apple that is freely available on the App Store. Depending on the specific version of XCode you are using you might also need to install the "Command Line Tools" package separately. Please refer to the XCode Documentation
    - Linux: Install GCC. Refer to the documentation of your distribution to find the specific package name

2. Make sure you have a recent (>= 1.9) version of `devtools` installed
```R
install.packages("devtools")
```
3. install the `flowCore` package
```R
source("http://bioconductor.org/biocLite.R")
biocLite("flowCore")
```

4. install `scgraphs` with the following command

```R
devtools::install_github("ParkerICI/scfeatures")
```


## Usage

This package enables the analysis of single-cell data using graphs, both unsupervised graphs as well as scaffold maps. While the package is designed to work with clusters generated from the [scfeatures](https://github.com/ParkerICI/scfeatures), any kind of tabular input data can be used as input

The documentation of each function can be accessed directly within R. The following snippet demonstrate typical usage

```R
# Use as input files that have been generated using scfeatures
input.files <- c("A.clustered.txt", "B.clustered.txt")

# Optional: Define a table of sample-level metadata. All the nodes derived from the corresponding cluster file will
# have vertex properties corresponding to this metadata ("response" and "pfs" in this example)
metadata.tab <- data.frame(filename = input.files, response = c("R", "NR"), pfs = c(12, 7))

# Define which columns contain variables that are going to be used to calculate similarities between the nodes
col.names <- c("foo", "bar", "foobar")


# The clusters in each one of the input files will be pooled together in a single graph
# This function also performs graph clustering by community detection. The community assignments are contained in
# the "community_id" vertex property of the resulting graph
G <- scgraphs::get_unsupervised_graph_from_files(input.files, metadata.tab = metadata.tab, 
            metadata.filename.col = "filename", col.names = col.names, filtering.threshold = 15)

# Write the resulting graph in graphml format. The Gephi software package (https://gephi.org/) is an excellent 
# solution to interactively manipulate the data
igraph::write.graph(G, "unsupervised.graphml", format = "grpahml")

```

