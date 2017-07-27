#' @include GOclass.R GOgeneric.R
#' @export GO
setGeneric(name="GO",
  def=function(gclus.obj, gset.obj, filterP=TRUE, filterOR=TRUE){
    standardGeneric("GO")
  }
)

#' @export write_GO
setGeneric(name="write_GO",
  def=function(go_set.obj, go_res.obj, nms) {
    standardGeneric("write_GO")
  }
)

#' @export simi
setGeneric(name="simi",
  def=function(go_set.obj, go_res.obj, nms) {
    standardGeneric("simi")
  }
)
