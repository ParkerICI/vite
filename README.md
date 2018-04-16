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

Functions documentation can be accessed directly within R

