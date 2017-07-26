# GOplot

## Use selDB to create/initialize a new gene set (gset) object
## "C1.CYTO": positional gene sets;
## "C2.CGP": chemical and genetic perturbations;
## "C2.CP": Canonical pathways (including BioCarta, KEGG and Reactome gene sets);
## "C3.MIR": microRNA targets;
## "C3.TFT": transcription factor targets;
## "C4.CGN": cancer gene neighborhoods;
## "C4.CM": cancer modules;
## "C5.BP": GO biological process;
## "C5.CC": GO cellular component;
## "C5.MF": GO molecular function;
## "C6.ONCOGENE": oncogenic signatures;
## "C7.IMMUNE": immunologic signatures.

For example,

gset.obj <- selDB(major="C2.CP", minor="Reactome")
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
