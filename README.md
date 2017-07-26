# GOplot

## Dependency of R packages
*library(tidyverse)
*library(ggplot2)
*library(ggdendro)
*library(RcolorBrewer)
*library(GOtest) # adapted library from MHW
*library(multiplot) # adapted code of [multiple graphs on one page] from (http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/)


## Installation
devtools::install_github("ericaenjoy3/GOplot")

## Use selDB function to create/initialize a new gene set (gset) object
Current available  database:
"C1.CYTO": positional gene sets;
"C2.CGP": chemical and genetic perturbations;
"C2.CP": Canonical pathways (including BioCarta, KEGG and Reactome gene sets);
"C3.MIR": microRNA targets;
"C3.TFT": transcription factor targets;
"C4.CGN": cancer gene neighborhoods;
"C4.CM": cancer modules;
"C5.BP": GO biological process;
"C5.CC": GO cellular component;
"C5.MF": GO molecular function;
"C6.ONCOGENE": oncogenic signatures;
"C7.IMMUNE": immunologic signatures.

For example,

gset.obj <- selDB(major="C2.CP", minor="Reactome", type="symbols", species="mouse")
or
gset.obj <- selDB(major="C5.BP")

## Initialize a gclus object (with colnames of gene and clus)
For example,

gclus.obj <- new("gclus", tbl=as_data_frame(data.frame(
  gene = as.character(c("Nanog","Rpl3","Rpl4","Mbl2","Ubr1","Herc2","Asb4","Rnf123","Klf4","Uba5")),
  clus = factor(c(rep(1,6),rep(2,4)))
  ))
)

## Perform GO analysis
res.list <- GO(gclus.obj, gset.obj, filterPADJ=FALSE, filterOR=TRUE)
go_set.obj <- res.list$go_set.obj
go_res.obj <- res.list$go_res.obj

## Output GO results to file
* nms is a directory with prefix for output
* output file named as [nms]_GOres.txt
write_GO((go_set.obj, go_res.obj, nms)

## Plot GO results to png
* output file named as [nms]_GOhclus[clusname].png
simi(go_set.obj, go_res.obj, nms)
