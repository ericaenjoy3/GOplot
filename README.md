# GOplot

## Use selDB to create/initialize a new gene set (gset) object
For example,
gset.obj<-selDB(major="C2.CP", minor="Reactome")

## Initialize a gclus object (with colnames of gene and clus)
For example,
gclus.obj<-new("gclus", tbl=as_data_frame(data.frame(
  gene=c("Nanog","Sox2","Sox3","Sox15","Sox18","Klf1","Klf2","Klf4","Klf5"),
  clus=factor(c(rep(1,5),rep(2,4)))
  ))
)

## Perform GO analysis
