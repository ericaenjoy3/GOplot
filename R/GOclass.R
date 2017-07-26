#' @import GOtest
#' @import tidyverse
#' @import ggplot2
#' @import ggdendro
#' @import RColorBrewer
#' @import multiplot
#' @importFrom msigdb msigdb.genesets
#' @importClassesFrom tidyverse tbl_df

# genes-clusters 1-to-1 relation
# gene sets: list (gene sets) of vectors (genes)
# per cluster GO terms: GOlabel, OR, pvalue, qvalue
# per cluster for a subset of GO terms: GO relationship matrix
#' @slot tbl tbl_df object  contains = c("tbl_df"),
#' @exportClass gset
gset<-setClass("gset",
  representation(tbl="tbl_df"),
  validity = function(object) {
    if (nrow(object@tbl)<1) {
      return("tbl is empty.")
    }
    if (!identical(colnames(object@tbl),c("gene","set"))) {
      return("colnames must be set to \"gene\" and \"set\".")
    }
    return(TRUE)
  }
)

#' @slot tbl tbl_df object
#' @exportClass gclus
gclus<-setClass("gclus",
  representation(tbl="tbl_df"),
  validity = function(object){
    if (nrow(object@tbl)<1) {
      return("tbl is empty.")
    }
    if (!identical(colnames(object@tbl),c("gene","clus"))){
      return("colnames must be set to \"gene\" and \"clus\".")
    }
    if (sum(object@tbl %>% group_by(gene) %>% summarise(any=n()>1) %>% pull(any))>0){
      return("gene column must be unique.")
    }
    return(TRUE)
  }
)

#' @slot tbl tbl_df object
#' @exportClass go_res
go_res<-setClass("go_res",
  representation(tbl="tbl_df"),
  validity = function(object){
    if (nrow(object@tbl)<1) {
      return("tbl is empty.")
    }
    if (!identical(colnames(object@tbl),c("set", "clus", "or", "pval", "padj"))){
      return("colnames must be set to \"set\", \"clus\", \"or\", \"pval\" and \"padj\".")
    }
    return(TRUE)
  }
)
