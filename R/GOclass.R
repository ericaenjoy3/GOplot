#####rlang .data prevents R CMD check from giving a NOTE about undefined global variables 

#' @import ggplot2
#' @import ggdendro
#' @import RColorBrewer
#' @import methods
#' @import multiplot
#' @importFrom msigdb msigdb.genesets
#' @importFrom dplyr %>% group_by summarise pull filter mutate left_join
#' @importFrom tibble as_tibble
#' @importFrom grDevices dev.off png
#' @importFrom stats as.dist fisher.test hclust p.adjust
#' @importFrom readr write_tsv
#' @importFrom rlang .data


# genes-clusters 1-to-1 relation
# gene sets: list (gene sets) of vectors (genes)
# per cluster GO terms: GOlabel, OR, pvalue, qvalue
# per cluster for a subset of GO terms: GO relationship matrix


#' @title An S4 class to represent gene sets.
#' @name gset-class
#' @rdname gset-class
#' @description Store gene sets in a tbl object.
#' @slot tbl tbl_df object  contains = c("tbl_df"),
#' @exportClass gset
gset<-setClass("gset",
  representation(tbl="tbl"),
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

#' @title An S4 class to represent gene cluster
#'
#' @description Store cluster sets in a tbl object.
#' @slot tbl tbl_df object
#' @exportClass gclus
gclus<-setClass("gclus",
  representation(tbl="tbl"),
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

#' @title An S4 class to represent GO enrichment analysis result
#'
#' @description Store GO results in a tbl object.
#' @slot tbl tbl_df object
#' @exportClass go_res
go_res<-setClass("go_res",
  representation(tbl="tbl"),
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
