#' @include GOclass.R

#' @title Odds Ratio Analysis
#' @name GO
#' @rdname GO-methods
# @aliases GO
#' @description
#' Odds ratio and P value are computed from \code{\link{fisher.test}}.
#' @param gclus.obj A \code{gclus} object.
#' @param gset.obj A \code{gset} object.
#' @param filterPADJ Logical, whether to filter with a P value threshold set by \code{Padj.cutoff}.
#' @param filterOR Logical, whether to filter out infinite ORs.
#' @param Padj.cutoff Numeric, specififying significance cutoff on adjust P value.
#' @details See \code{\link{fisher.test}}. An alternative and slower approach is to compute odds ratio by a contingency table and estiamte P value from repeated reshuffle analysis.
#' @return
#' A list with two entries \code{go_res.obj} and \code{go_set.obj}.
#' @examples
#' \dontrun{
#' library(GOplot)
#' gset.obj <- selDB(major="C2.CP", minor="Reactome", type="symbols", species="mouse")
#' gclus.obj <- new("gclus", tbl=tibble:::as_tibble(data.frame( gene = c("Nanog","Rpl3","Rpl4","Mbl2","Ubr1","Herc2","Asb4","Rnf123","Klf4","Uba5"), clus = factor(c(rep(1,6),rep(2,4))) )) )
#' res.list <- GO(gclus.obj, gset.obj, filterPADJ=FALSE, filterOR=TRUE)
#' go_set.obj <- res.list$go_set.obj
#' go_res.obj <- res.list$go_res.obj
#' write_GO(go_set.obj, go_res.obj, 'output_txt')
#' simi(go_set.obj, go_res.obj, 'output_fig')
#' }
#' @exportMethod GO
setGeneric(name="GO",
  def=function(gclus.obj, gset.obj, filterPADJ=TRUE, filterOR=TRUE, Padj.cutoff=0.05){
    standardGeneric("GO")
  },signature=c("gclus.obj", "gset.obj")
)


#' @title Output GO Results
#' @name write_GO
#' @rdname write_GO-methods
# @aliases write_GO
#' @description Write GO results into a tsv file
#' @param go_set.obj A go_set object.
#' @param go_res.obj A go_res object.
#' @param nms A directory with prefix for output.
#' @export write_GO
setGeneric(name="write_GO",
  def=function(go_set.obj, go_res.obj, nms) {
    standardGeneric("write_GO")
  }
)

#' @title Presentation of similar gene sets
#' @name simi
#' @rdname simi-methods
# @aliases simi
#' @description Estimate similarity between gene set terms and plot clustering
#' @param go_set.obj A go_set object.
#' @param go_res.obj A go_res object.
#' @param nms A directory with prefix for output.
#' @export simi
setGeneric(name="simi",
  def=function(go_set.obj, go_res.obj, nms) {
    standardGeneric("simi")
  }
)
