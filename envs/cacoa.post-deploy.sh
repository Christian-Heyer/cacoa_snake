#!/bin/bash

Rscript -e 'install.packages(c("coda.base", "svglite"), repos ="https://cloud.r-project.org")'
Rscript -e 'devtools::install_github("kharchenkolab/sccore", ref="dev")'
Rscript -e 'devtools::install_github("kharchenkolab/cacoa")'

