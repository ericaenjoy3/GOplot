#' @include GOclass.R GOgeneric.R
#' @rdname GO-methods
setMethod(f = "GO",
  signature=c(gclus.obj = "gclus", gset.obj = "gset"),
  definition = function(gclus.obj, gset.obj, filterPADJ=TRUE, filterOR=TRUE, Padj.cutoff=0.05) {
    # filter gclus.obj based all.gene
    gclus.tbl <- filter(gclus.obj@tbl, .data$gene %in% gset.obj@tbl$gene) %>% droplevels()
    gset.tbl <- filter(gset.obj@tbl, .data$gene %in% gclus.tbl$gene) %>% droplevels()
    res <- expand.grid(levels(gset.tbl$set),levels(gclus.tbl$clus))
	if(any(dim(res)==0)) stop('The groups in gclus.obj and gset.obj must be factors\n')
    pval <- OR <- rep(NA, nlevels(gclus.tbl$clus)*nlevels(gset.tbl$set))
    # per cluster
    iter <- 0
    for (i in seq_along(levels(gclus.tbl$clus))) {
      clus.in <- {filter(gclus.tbl, .data$clus==levels(gclus.tbl$clus)[i]) %>% pull('gene')}
      clus.out <- {filter(gclus.tbl, .data$clus!=levels(gclus.tbl$clus)[i]) %>% pull('gene')}
      # per gene set
      for (j in seq_along(levels(gset.tbl$set))){
        iter <- iter+1
        path.in <- {filter(gset.tbl, .data$set==levels(gset.tbl$set)[j]) %>% pull('gene') %>% unique() }
        # g in clus c and in pathway p
        cp1 <- length(intersect(clus.in, path.in))
        # g in clus c and not in pathway p
        cp2 <- length(setdiff(clus.in, path.in))
        # g not in clus c but in pathway p
        cp3 <- length(setdiff(path.in, clus.in))
        # g not in clus c and not in pathway p
        cp4 <- length(setdiff(clus.out, path.in))
        htest.obj <- fisher.test(matrix(c(cp1,cp2,cp3,cp4),nrow=2,byrow=TRUE), alternative="greater")
        OR[iter] <- htest.obj$estimate
        pval[iter] <- htest.obj$p.value
      }
    }
    res <- as_tibble(cbind(res, OR, pval, padj = p.adjust(pval, method = "fdr")))
    # filter by pval
    if (filterPADJ) res <- filter(res, .data$pval<Padj.cutoff) %>% droplevels()
    if (filterOR) res <- filter(res, is.finite(.data$OR)) %>% droplevels()
    colnames(res) <- c("set", "clus", "or", "pval", "padj")
    go_set.obj <- new("gset", tbl = filter(gset.tbl, .data$set %in% res$set) %>% droplevels())
    go_res.obj <- new("go_res",tbl = res)
    return(list(go_res.obj = go_res.obj, go_set.obj = go_set.obj))
  }
)

#' @rdname write_GO-methods
setMethod(f = "write_GO",
  signature = c(go_set.obj="gset", go_res.obj="go_res", nms="character"),
  definition = function(go_set.obj, go_res.obj, nms) {
    #if(is.na(file.info(nms)$isdir) || !file.info(nms)$isdir) dir.create(nms, showWarnings = FALSE, recursive = TRUE)
    # merge genes for each set into one-line
    gset.dict <- go_set.obj@tbl %>% group_by(.data$set) %>% mutate(gene=paste(.data$gene,collapse=";")) %>% unique()
    dat <- go_res.obj@tbl[,c("clus","set", "or", "pval", "padj")] %>% left_join(gset.dict, by="set")
    write_tsv(dat, paste0(nms,"_GOres.txt"), col_names = TRUE)
    invisible(TRUE)
  }
)

#' @rdname simi-methods
setMethod(f = "simi",
  signature = c(go_set.obj="gset", go_res.obj="go_res", nms="character"),
  definition = function(go_set.obj, go_res.obj, nms) {
    #if(is.na(file.info(nms)$isdir) || !file.info(nms)$isdir) dir.create(nms, showWarnings = FALSE, recursive = TRUE)
    for (i in seq_along(levels(go_res.obj@tbl$clus))) {
      per_clus <- filter(go_res.obj@tbl, .data$clus==levels(go_res.obj@tbl$clus)[i]) %>% droplevels()
      path <- filter(go_set.obj@tbl, .data$set %in% pull(per_clus, "set")) %>% droplevels()
      mat <- matrix(NA, nrow=nlevels(path$set), ncol=nlevels(path$set))
      colnames(mat) <- rownames(mat) <-levels(path$set)
      diag(mat) <- 1
      for (n1 in 1:(nlevels(path$set)-1)) {
        n1.g <- filter(path, .data$set==levels(path$set)[n1]) %>% pull('gene') %>% unique()
        for (n2 in (n1+1):nlevels(path$set)) {
          n2.g <- filter(path, .data$set==levels(path$set)[n2]) %>% pull('gene') %>% unique()
          shared.g <- length(intersect(n1.g,n2.g))
          min.g <- min(length(n1.g), length(n2.g))
          shared.p <- shared.g/min.g
          mat[n1,n2] <- mat[n2,n1] <- shared.p
        }
      }
      # cluster analysis
      hc <- hclust(as.dist(1-mat))
      dendr <- dendro_data(hc, type="rectangle")
      colors <- per_clus[order(match(per_clus$set,colnames(mat))),"or"]
      p1=ggplot() +
      geom_segment(data=segment(dendr), aes_(x=~x, y=~y, xend=~xend, yend=~yend)) +
      geom_text(data=label(dendr), aes_(~x, ~y, label=~label, hjust=0, color=~colors), size=3) +
      labs(y="", color="OR")+
      coord_flip() + scale_y_reverse(expand = c(0, 5), breaks=seq(0,1,0.5), labels=seq(1,0,-0.5))+
      scale_colour_gradientn(colours=rev(brewer.pal(3, 'Dark2')))+
      theme(legend.position="bottom",
        axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())
      png(filename=paste0(nms, "_GOhclus", levels(go_res.obj@tbl$clus)[i], ".png"), height=3000, width=3000, res=300)
      multiplot(p1,cols=1)
      dev.off()
      invisible(TRUE)
    }
  }
)
