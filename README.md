# GOplot

## Dependency of R packages
* library(dplyr)
* library(tibble)
* library(readr)
* library(ggplot2)
* library(ggdendro)
* library(rlang)
* library(msigdb) # https://github.com/mw201608/msigdb
* library(multiplot) # https://github.com/ericaenjoy3/multiplot
* library(RColorBrewer)
* library(methods)
* library(stats)
* library(grDevices)


## Installation
```
devtools::install_github("ericaenjoy3/GOplot")
```
## Load the package
```
library(GOplot)
```
## Use selDB function to create/initialize a new gene set (gset) object
Current available  database:
* "C1.CYTO": positional gene sets;
* "C2.CGP": chemical and genetic perturbations;
* "C2.CP": Canonical pathways (including BioCarta, KEGG and Reactome gene sets);
* "C3.MIR": microRNA targets;
* "C3.TFT": transcription factor targets;
* "C4.CGN": cancer gene neighborhoods;
* "C4.CM": cancer modules;
* "C5.BP": GO biological process;
* "C5.CC": GO cellular component;
* "C5.MF": GO molecular function;
* "C6.ONCOGENE": oncogenic signatures;
* "C7.IMMUNE": immunologic signatures.

For example, to load canonical pathways
```
gset.obj <- selDB(major="C2.CP", minor="Reactome", type="symbols", species="mouse")
```
or to load gene ontology biological process
```
gset.obj <- selDB(major="C5.BP", minor=NA, type="symbols", species="mouse")
```
## Initialize a gclus object (with colnames of gene and clus)
For example, prepare gene signatures to be tested
```
signatures=data.frame(gene = as.character(c("Nanog","Rpl3","Rpl4","Mbl2","Ubr1","Herc2","Asb4","Rnf123","Klf4","Uba5")),
	clus = factor(rep(c('Group1','Group2'),c(6,4))))
gclus.obj <- new("gclus", tbl=tibble:::as_tibble(signatures))
```
## Perform GO analysis
```
res.list <- GO(gclus.obj, gset.obj, filterPADJ=FALSE, filterOR=TRUE)
go_set.obj <- res.list$go_set.obj
go_res.obj <- res.list$go_res.obj
```
## Output GO results to file
* nms is the file prefix for output
* output file named as [nms]_GOres.txt
```
write_GO(go_set.obj, go_res.obj, nms='test')
```
## Plot GO results to png
* output file named as [nms]_GOhclus[clusname].png
```
simi(go_set.obj, go_res.obj, nms='test')
```
