![breakpointR](breakpointR_logo.png)
====================================

An R-package for breakpoint detection in single cell Strand-seq data.

Collaborators: David Porubsky, Ashley D. Sanders, Aaron Taudt

Installation
------------

### Stable bioconductor version
To install the current stable version from Bioconductor, please visit http://bioconductor.org/packages/breakpointR/

### Development version from Github
To install the development version from Github, follow the steps given below. The installation has only been tested on Ubuntu so far, if you need to install on Windows or Mac additional steps might be necessary (e.g. installation of Rtools from https://cran.r-project.org/bin/windows/Rtools/)

1. Install a recent version of R (>=3.5.0) from https://www.r-project.org/
2. Optional: For ease of use, install Rstudio from https://www.rstudio.com/
3. Open R and install all dependencies. Please ensure that you have writing permissions to install packages. Execute the following lines one by one:

   #### To install required packages  
   if (!requireNamespace("BiocManager", quietly=TRUE))
   install.packages("BiocManager")
   install("GenomicRanges") 
	 install("GenomicAlignments")
	 install.packages("devtools")
	 library(devtools)  

4. To install breakpointR package from github	 
	 #### Option1
	 install_github("daewoooo/breakpointRdata", force=TRUE)  
	 install_github("daewoooo/breakpointR", force=TRUE)  
	 #### Option2 
	 install_git("git://github.com/daewoooo/breakpointRdata.git", branch = "master")  
	 install_git("git://github.com/daewoooo/breakpointR.git", branch = "master")  

How to use breakpointR
----------------------

Please refert to the [vignette](https://github.com/daewoooo/breakpointR/blob/master/vignettes/breakpointR-knitr.pdf) for tutorials on all breakpointR features.

Report Errors
-------------

If you encounter errors of any kind, please report an [issue here](https://github.com/daewoooo/breakpointR/issues/new).

NOTE
----

The breakpointR package is currently under development and contains unpublished work. Any usage for publishing is strictly prohibited without permission.
