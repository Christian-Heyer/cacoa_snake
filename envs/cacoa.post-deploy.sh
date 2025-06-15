#!/bin/bash

Rscript -e 'install.packages(c("coda.base"), repos ="https://cloud.r-project.org")'
Rscript -e 'devtools::install_github("kharchenkolab/sccore")'
Rscript -e 'devtools::install_github("Christian-Heyer/cacoa")'
