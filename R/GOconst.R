#' @include GOclass.R
#' select db for GO
#' initate gset class
#' minor must be set to NULL to de-select minor DB.
#' @export selDB
selDB<-function(major=c("C2.CP","C3.MIR","C3.TFT","C5.BP","C5.CC"), minor=c("Reactome"), type=c('symbols', 'entrez'), species=c('human','mouse')) {
  major <- match.arg(major,several.ok=FALSE)
  type <- match.arg(type)
  species <- match.arg(species)
  db <- msigdb.genesets(sets=major, type="symbols", species="mouse")
  if (!missing(minor)) {
    minor <- match.arg(minor,several.ok=FALSE)
    idx <- grep(minor,db[["geneset.names"]],ignore.case=TRUE)
    db <- lapply(db,function(lvec)lvec[idx])
  }
  db[["geneset.names"]] <- gsub(".+:", "", db[["geneset.names"]])
  if (any(grepl("Reactome_*", db[["geneset.names"]], ignore.case=TRUE))) {
    db[["geneset.names"]] <- gsub("Reactome_*", "", db[["geneset.names"]],ignore.case=TRUE)
  }
  sizes <- sapply(db[["genesets"]], length)
  dat <- as_data_frame(data.frame(gene = unlist(db[["genesets"]]), set = rep(db[["geneset.names"]], sizes), stringsAsFactors = TRUE))
  return(new("gset",tbl=dat))
}
