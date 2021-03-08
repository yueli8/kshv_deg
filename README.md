# kshv_deg
This repository includes code for processing data of kshv_deg. kshv_deg.Rmd script has code for the results.

## Install the software, download the data and set up the directory
This pipeline is designed to be used in R environment.

1. Install the R statistical package. We used version 4.0.4.

2. Install the following R packages, which can be obtained using either the install.packages function in R or via the Bioconductor framework:

* limma
* DESeq2
* calibrate
* ggplot
* fgsea
* tidyverse
* data.table

### deg
using uci.txt, the average_if value of each genes can be got. Then, using delta11_vs_mock.txt, delta_vs_wt.txt and mock_vs_wt.txt, the differentially expressed genes (degs) can be got.
