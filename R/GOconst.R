#' @include GOclass.R
#' @title Select db for GO
#' @description Initate gset class by making use of gene sets from package \code{msgidb}
#' @param major, the major category.
#' @param minor, a character string to subset through grep. Set NA to ignore.
#' @param type, either symbols or entrez.
#' @param species, either human or mouse.
#' @details For details about the gene sets, please see package \code{msgidb}.
#' @return
#' A \code{gset} object.
#' @export selDB
selDB<-function(major=c("C2.CP","C3.MIR","C3.TFT","C5.BP","C5.CC"), minor=NA, type=c('symbols', 'entrez'), species=c('human','mouse')) {
  major <- match.arg(major,several.ok=TRUE)
  type <- match.arg(type)
  species <- match.arg(species)
  db <- msigdb.genesets(sets=major, type=type, species=species)
  if (!is.na(minor)) {
    idx <- grep(minor,db[["geneset.names"]],ignore.case=TRUE)
    db <- lapply(db,function(lvec)lvec[idx])
  }
  db[["geneset.names"]] <- gsub(".+:", "", db[["geneset.names"]])
  if (any(grepl("Reactome_", db[["geneset.names"]], ignore.case=TRUE))) {
    db[["geneset.names"]] <- gsub("Reactome_", "", db[["geneset.names"]],ignore.case=TRUE)
  }
  sizes <- sapply(db[["genesets"]], length)
  dat <- as_tibble(data.frame(gene = unlist(db[["genesets"]]), set = rep(db[["geneset.names"]], sizes), stringsAsFactors = TRUE))
  return(new("gset",tbl=dat))
}
