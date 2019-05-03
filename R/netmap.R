#' Plot network of Terms
#' @param df DGE files (DESeq2 result files) or vector contains gene names
#' @param rhs Enrichment analsyis result
#' @param pvalue.cutoff the cut-off P value for selecting significant Terms
#' @param padj.cutoff the cut-off P adjust value for selecting significant Terms
#' @param weightcut the weight cut value for remove edges
#' @param useTerm use the Term description or not(defalutTRUE)
#' @param writeCyt export file for Cyt software
#' @param vertex.label.color color of label(defaultblack)
#' @param vertex.label.cex size of label(default0.5)
#' @param layout layout format (defultlayout.fruchterman.reingold)
#' @param visNet use VisNetwork method to display network(defaultFALSE)
#' @param smooth use smooth edge for visNetwork
#' @param nodeselect select some interesting node(defaultFALSE)
#' @param editvis choose to edit network(defaultFALSE)
#' @param savevis save visnetwork to html
#' @param savefig save figure to pdf
#' @param filename output filename
#' @param top number of Terms you want to display
#' @export
#' @author Kai Guo
netmap<-function (df, rhs, top = 50, pvalue.cutoff = 0.05, padj.cutoff = NULL,
                  weightcut = 0.2, useTerm = TRUE, writeCyt = FALSE,cytoscapeFile = "network-file-for-cytoscape.txt", vertex.label.font = 2,
                  vertex.label.color = "black", vertex.label.cex = 0.5, layout = layout.fruchterman.reingold,
                  visNet = FALSE,smooth=TRUE,nodeselect=FALSE,editvis=FALSE,savevis=FALSE,savefig=FALSE,filename="network",...)
{
  options(stringsAsFactors = F)
  suppressMessages(library(reshape2))
  suppressMessages(library(igraph))
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
  if (writeCyt == TRUE) {
    write.table(wn, file = cytoscapeFile, sep = "\t",
                row.names = F, quote = F)
  }
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
  if (visNet == TRUE) {
    suppressMessages(library(visNetwork))
    graph <- visIgraph(g, smooth = smooth)
    if(nodeselect==TRUE){
      graph<-graph%>%visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)%>%
        visInteraction(navigationButtons = TRUE)
    }
    if(editvis==TRUE){
      graph<-graph%>%visInteraction(navigationButtons = TRUE)%>%visOptions(manipulation = TRUE)

    }
    if(savevis==TRUE){
      visSave(graph,file = paste(filename,"html",sep="."))
    }
    graph
  }
  else {
    l <- layout_with_fr(g)
    l <- norm_coords(l, ymin = -1, ymax = 1, xmin = -1, xmax = 1)
    par(mar = c(2, 2, 2, 2))
    plot.igraph(g, vertex.label.font = vertex.label.font,
                vertex.label.color = "black", vertex.label.cex = vertex.label.cex,
                vertex.frame.color = V(g)$color, rescale = F, layout = l *
                  0.8)
    if(savefig==TRUE){
      dev.print(pdf,file=paste(filename,"pdf",sep="."))
    }
  }
}
.color_scale <- function(c1="pink", c2="red") { #modified from DOSE
  pal <- colorRampPalette(c(c1, c2))
  colors <- pal(200)
  return(colors)
}
.getIdx <- function(v, MIN, MAX) { #modified from DOSE
 # if ( MIN == MAX ) {
  #  return(200)
 # }
  intervals <- seq(MIN, MAX, length.out=200)
  max(which(intervals <= v))
}
