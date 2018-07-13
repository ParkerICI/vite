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

4. install `vite` with the following command

```R
devtools::install_github("ParkerICI/vite")
```


## Usage

This package enables the analysis of single-cell data using graphs, both unsupervised graphs as well as scaffold maps. While the package is designed to work with clusters generated from the [grappolo](https://github.com/ParkerICI/grappolo) package, any kind of tabular input data can be used as input

The documentation of each function can be accessed directly within R. The following snippet demonstrate typical usage. Please refer to the full documentation for a complete breakdown of all the options

### Creating an unsupervised graph

```R
# Use as input files that have been generated using grappolo
input.files <- c("A.clustered.txt", "B.clustered.txt")

# Optional: Define a table of sample-level metadata. All the nodes derived from the corresponding cluster file will
# have vertex properties corresponding to this metadata ("response" and "pfs" in this example)
metadata.tab <- data.frame(filename = input.files, response = c("R", "NR"), pfs = c(12, 7))

# Define which columns contain variables that are going to be used to calculate similarities between the nodes
col.names <- c("foo", "bar", "foobar")


# The clusters in each one of the input files will be pooled together in a single graph
# This function also performs graph clustering by community detection. The community assignments are contained in
# the "community_id" vertex property of the resulting graph
G <- vite::get_unsupervised_graph_from_files(input.files, metadata.tab = metadata.tab, 
            metadata.filename.col = "filename", col.names = col.names, filtering.threshold = 15)

# Write the resulting graph in graphml format. 
igraph::write.graph(G, "unsupervised.graphml", format = "grpahml")

```

### Running a Scaffold analysis

This code snippet demonstrates how to construct scaffold maps. This assumes that the data for the landmark nodes, i.e. the gated populations, is in a subfolder called `gated`. The gated populations have to be provided as single FCS files (one for each population). The software will split the name of the FCS files using `"_"` as separator and the last field will be used as the population name. For instance if a file is called `foo_bar_foobar_Bcells.fcs` the corresponding population will be called `Bcells` in the scaffold analysis.


```R
# Use as input files that have been generated using grappolo
input.files <- c("A.clustered.txt", "B.clustered.txt")

# Define which columns contain variables that are going to be used to calculate similarities between the nodes
col.names <- c("foo", "bar", "foobar")

# Load the data for the landmarks
landmarks.data <- load_landmarks_from_dir("gated/", asinh.cofactor = 5, transform.data = T)
    
# Run the analysis. By default results will be save in a directory called "scaffold_result"
run_scaffold_analysis(input.files, input.files[1], landmarks.data, col.names)
```

#### Scaffold output

By the default the output of the analysis will be saved in a folder called `scaffold_result`. The directory will contain a `graphml` file for each Scaffold map and two sub-folders called `clusters_data` and `landmarks_data`.

These folders contain downsampled single-cell data for the clusters and landmarks, to be used for visualization. The `clusters_data` folder will contain a separate sub-folder for each `graphml` file, containing the data specific to that sample. The data is split in multiple `rds` files, one for each cluster (or landmark in `landmarks_data`). 

If the Scaffold analysis was constructed from data that was pooled before clustering (i.e. using `grappolo::cluster_fcs_files_groups`), the `clusters_data` folder will also contain a subfolder called `pooled`, containing the pooled data, in addition to the sample-specific folders described above.

### Using the GUI

A GUI is available to launch either an unsupervised graph analysis or a Scaffold analysis. The GUI allows you to specify all the input options in a graphical environment, instead of having to write R code.

To launch the GUI type the following in your R console

```R
vite::vite_GUI()
```

When the GUI starts you will be prompted to select a working directory. This directory must contain all the files that you want to include in the analysis. Select any file in that directory, and the directory that contains the file will be selected as working directory.

