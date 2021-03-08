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
using uci.txt, the average_if value of each genes can be got. Then, using delta11_vs_mock.txt, delta_vs_wt.txt and mock_vs_wt.txt, Deseq2_uci.R, the differentially expressed genes (degs) can be got.

###fgsea
Using c2.cp.v7.2.symbols.gtm and c2.cp.kegg.v7.2.symbols.gmt, gene name and the t-test value to generate the kegg and canonical pathway. The input files are delta11_vs_mock_fgsea.txt, delta11_vs_we_fgsea.txt and mock_vs_wt_fgsea.txt, by using fgsea_uci.R, can get the kegg and the whole canonical pathways.

## Contact information

* Moom R. Roosam. [roosan@chapman.edu](mailto:roosan@chapman.edu)
* Yue Li. [yli1@chapman.edu](mailto:yli1@chapman.edu)
