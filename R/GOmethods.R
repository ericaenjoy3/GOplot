#' @include GOclass.R
#' calculate OR
#' ((g present in clus c and pathway p)/(g present in clus c but absent in pathway p))((g absent in clus c but present in pathway p)/(g absent in clus c and absent in pathway p))
#' reshuffle gene sets 100,000 times, and calcualte OR
#' calculate Pval
#' calculate Padj
setMethod(f = "GO",
  signature=c(gclus.obj = "gclus", gset.obj = "gset", filterPADJ="logical", filterOR="logical"),
  definition = function(gclus.obj, gset.obj, filterPADJ, filterOR) {
    # filter gclus.obj based all.gene
    gclus.tbl <- filter(gclus.obj@tbl, gene %in% gset.obj@tbl$gene) %>% droplevels()
    gset.tbl <- filter(gset.obj@tbl, gene %in% gclus.tbl$gene) %>% droplevels()
    res <- expand.grid(levels(gset.tbl$set),levels(gclus.tbl$clus))
    pval <- OR <- rep(NA, nlevels(gclus.tbl$clus)*nlevels(gset.tbl$set))
    # per cluster
    iter <- 0
    for (i in seq_along(levels(gclus.tbl$clus))) {
      clus.in <- {filter(gclus.tbl, clus==levels(gclus.tbl$clus)[i]) %>% pull(gene)}
      clus.out <- {filter(gclus.tbl, clus!=levels(gclus.tbl$clus)[i]) %>% pull(gene)}
      # per gene set
      for (j in seq_along(levels(gset.tbl$set))){
        iter <- iter+1
        path.in <- {filter(gset.tbl, set==levels(gset.tbl$set)[j]) %>% pull(gene) %>% unique() }
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
    # filter by frd<0.01
    if (filterPADJ) res <- filter(res, padj<0.01) %>% droplevels()
    if (filterOR) res <- filter(res, is.finite(OR)) %>% droplevels()
    colnames(res) <- c("set", "clus", "or", "pval", "padj")
    go_set.obj <- new("gset", tbl = filter(gset.tbl, set %in% res$set) %>% droplevels())
    go_res.obj <- new("go_res",tbl = res)
    return(list(go_res.obj = go_res.obj, go_set.obj = go_set.obj))
  }
)

#' output GO results
setMethod(f = "write_GO",
  signature = c(go_set.obj="gset", go_res.obj="go_res", nms="character"),
  definition = function(go_set.obj, go_res.obj, nms) {
    if(is.na(file.info(nms)$isdir) || !file.info(nms)$isdir) dir.create(nms, showWarnings = FALSE, recursive = TRUE)
    # merge genes for each set into one-line
    gset.dict <- go_set.obj@tbl %>% group_by(set) %>% mutate(gene=paste(gene,collapse=";")) %>% unique()
    dat <- go_res.obj@tbl[,c("clus","set", "or", "pval", "padj")] %>% left_join(gset.dict, by="set")
    write_tsv(dat, paste0(nms,"_GOres.txt"), col_names = TRUE)
    invisible(TRUE)
  }
)

#' estimate similarity between gene set terms
setMethod(f = "simi",
  signature = c(go_set.obj="gset", go_res.obj="go_res", nms="character"),
  definition = function(go_set.obj, go_res.obj, nms) {
    if(is.na(file.info(nms)$isdir) || !file.info(nms)$isdir) dir.create(nms, showWarnings = FALSE, recursive = TRUE)
    for (i in seq_along(levels(go_res.obj@tbl$clus))) {
      per_clus <- filter(go_res.obj@tbl, clus==levels(go_res.obj@tbl$clus)[i]) %>% droplevels()
      path <- filter(go_set.obj@tbl, set %in% pull(per_clus, "set")) %>% droplevels()
      mat <- matrix(NA, nrow=nlevels(path$set), ncol=nlevels(path$set))
      colnames(mat) <- rownames(mat) <-levels(path$set)
      diag(mat) <- 1
      for (n1 in 1:(nlevels(path$set)-1)) {
        n1.g <- filter(path, set==levels(path$set)[n1]) %>% pull(gene) %>% unique()
        for (n2 in (n1+1):nlevels(path$set)) {
          n2.g <- filter(path, set==levels(path$set)[n2]) %>% pull(gene) %>% unique()
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
      geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) +
      geom_text(data=label(dendr), aes(x, y, label=label, hjust=0, color=colors), size=3) +
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
      png(pngfout=paste0(nms, "_GOhclus", levels(go_res.obj@tbl$clus)[i], ".png"), height=3000, width=3000, res=300)
      multiplot(p1,cols=1)
      dev.off()
      invisible(TRUE)
    }
  }
)
