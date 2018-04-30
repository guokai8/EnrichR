gnet<-function (df, rhs, top = 50, pvalue.cutoff = 0.05, padj.cutoff = NULL,
                weightcut = 0.2, useTerm = TRUE, writeCyt = FALSE, vertex.label.font = 2,
                vertex.label.color = "black", vertex.label.cex = 0.5, layout = layout.fruchterman.reingold,savefig=FALSE,filename="network",
                width=7,height=7,...)
{
  options(stringsAsFactors = F)
  suppressMessages(library(reshape2))
  suppressMessages(library(igraph))
  suppressMessages(library(GGally))
  suppressMessages(library(ggrepel))
  suppressMessages(library(intergraph))
  if (!is.null(padj.cutoff)) {
    rhs <- rhs[rhs$Padj < padj.cutoff, ]
  }
  else {
    rhs <- rhs[rhs$Pvalue < pvalue.cutoff, ]
  }
  if (nrow(rhs) <= top) {
    rhs <- rhs
  }
  else {
    rhs <- rhs[1:top, ]
  }
  if (is.data.frame(df)) {
    gene_p <- -log10(df$padj)
    names(gene_p) <- rownames(df)
  }
  else {
    gene_p <- rep(1, length(df))
    names(gene_p) <- df
  }
  pvalue = rhs$Pvalue
  names(pvalue) <- rownames(rhs)
  go2gen <- strsplit(x = as.vector(rhs$GeneID), split = ",")
  names(go2gen) <- rownames(rhs)
  gen2go <- reverseList(go2gen)
  golen <- rhs$Significant
  names(golen) <- rownames(rhs)
  gen2golen <- lapply(gen2go, function(x) golen[x])
  gen2gosum <- lapply(gen2golen, function(x) sum(x)/x)
  gen2res <- lapply(gen2gosum, function(x) x/sum(x))
  id <- rownames(rhs)
  n = nrow(rhs)
  w <- matrix(NA, nrow = n, ncol = n)
  colnames(w) <- rownames(w) <- rownames(rhs)
  for (i in 1:n) {
    ni <- id[i]
    for (j in i:n) {
      nj <- id[j]
      genein = intersect(go2gen[[ni]], go2gen[[nj]])
      geneup <- sum(gene_p[genein] * unlist(lapply(lapply(gen2res[genein],
                                                          "[", c(ni, nj)), sum)))
      genei <- setdiff(go2gen[[ni]], go2gen[[nj]])
      genej <- setdiff(go2gen[[nj]], go2gen[[ni]])
      geneid <- sum(gene_p[genei] * unlist(lapply(lapply(gen2res[genei],
                                                         "[", ni), sum)))
      genejd <- sum(gene_p[genej] * unlist(lapply(lapply(gen2res[genej],
                                                         "[", nj), sum)))
      gened <- geneup + geneid + genejd
      w[i, j] <- geneup/gened
    }
  }
  if (useTerm == TRUE) {
    colnames(w) <- rownames(w) <- rhs$Term
    names(pvalue) = rhs$Term
  }
  wn <- melt(w, as.is = TRUE)
  wn <- wn[wn[, 1] != wn[, 2], ]
  wn <- wn[!is.na(wn[, 3]), ]
  wn <- wn[wn[, 3] > 0, ]
  g <- igraph::graph.data.frame(wn[, -3], directed = F)
  E(g)$width = sqrt(wn[, 3] * 5)
  pvalue = pvalue[V(g)$name]
  if (useTerm == TRUE) {
    idx <- unlist(sapply(V(g)$name, function(x) which(x ==
                                                        rhs$Term)))
  }
  else {
    idx <- unlist(sapply(V(g)$name, function(x) which(x ==
                                                        rownames(rhs))))
  }
  cols <- .color_scale("red", "orange")
  V(g)$color <- cols[sapply(pvalue, .getIdx, min(pvalue), max(pvalue))]
  g <- igraph::delete.edges(g, E(g)[wn[, 3] < weightcut])
  gs <- rhs$Significant
  if (useTerm == TRUE) {
    names(gs) <- rhs$Term
  }
  else {
    names(gs) <- rownames(rhs)
  }
  V(g)$size <- log(gs[V(g)$name], base = 10) * 10
  p<-ggnet2(g, node.size = V(g)$size, node.color = V(g)$color,
         edge.size = E(g)$width/10) + geom_text_repel(label = V(g)$name,size=vertex.label.cex,segment.size=0.2)+theme(legend.position = "none")
  print(p)
    if(savefig==TRUE){
    ggsave(p,file=paste(filename,"pdf",sep="."),width=width,height = height)
  }
  }
