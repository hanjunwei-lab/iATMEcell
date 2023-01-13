# iATMEcell
A software R package for identification of abnormal tumor microenvironment (TME) cells.

> A computing tool is developed to identify abnormal TME cells based on cell-cell crosstalk network and to predict the clinical outcomes by constructing survival models. The operation modes including: i) calculate the gene differential expression level between labels of two different conditions (e.g. disease and normal, dead and alive). ii) construct the cell-cell crosstalk network used a network propagated algorithm to calculate the centrality scores of cells. iii) Select abnormal cells to construct survival models using regression analysis to verify cells prognostic efficacy.

# Installation
> You can install this package from GitHub.
```
library(devtools); 
install_github("hanjunwei-lab/iATMEcell", build_vignettes = TRUE)

Use：
library(iATMEcell)
```
# Please cite the following article when using iTMEcell:
> Waiting for publication
