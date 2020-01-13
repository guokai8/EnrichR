##' richplot
##' @description plot the sigificant terms and shared genes with network format
##' @importFrom GGally ggnet2
##' @importFrom ggrepel geom_text_repel
##' @importFrom igraph graph_from_data_frame
##' @importFrom igraph simplify
##' @importFrom ggrepel geom_text_repel
##' @importFrom igraph V
##' @importFrom ggplot2 geom_text
##' @param resultFis Enrichment results
##' @param top number of terms to show (default: 50)
##' @param pvalue.cutoff cutoff p value for enrichment result
##' @param padj.cutoff cutoff p adjust value for enrichment result
##' @param low color used for small value
##' @param high color used for large value
##' @param useTerm use terms for nodes (default: TRUE)
##' @param writeCyt write out the cytoscape file
##' @param cytoscapeFile output cytoscape File
##' @param label.color label color
##' @param label.size label size
##' @param node.shape vector of shape and names of the vector should be the terms (default: 20)
##' @param layout layout method
##' @param savefig save figures or not
##' @param filename output figure name
##' @param width width for output figure
##' @param height height for output figure
##' @param node.alpha alpha-transparency scales
##' @param node.shape shape of the node
##' @param repel use ggrepel text function or not
##' @param segment.size segment size for ggrepel text
##' @export
richplot <- function(resultFis,top=50, pvalue.cutoff=0.05, padj.cutoff=NULL,
                usePadj =TRUE, useTerm=TRUE,low="orange",high="red",
                   writeCyt=FALSE, cytoscapeFile="network-file-for-cytoscape.txt",
                   label.color = "black", label.size = 2, node.shape=NULL,
                   layout = layout.fruchterman.reingold,savefig=FALSE,filename="network",
                   width=7,height=7,node.alpha=0.7,repel=TRUE,segment.size=0.2){
    suppressMessages(library(igraph))
    if(!is.null(padj.cutoff)){
        resultFis<-resultFis[resultFis$Padj<padj.cutoff,]
    }else{
        resultFis<-resultFis[resultFis$Pvalue<pvalue.cutoff,]
    }
    if(nrow(resultFis)>=top){
        resultFis<-resultFis[1:top,]
    }
    resultFis$GeneID<-as.vector(resultFis$GeneID)
    lhs <- strsplit(resultFis$GeneID,",")
    if(isTRUE(useTerm)){
        rhs <- data.frame(from=unlist(lhs),to=rep(resultFis$Term,lapply(lhs,length)))
        rownames(resultFis) <- resultFis$Term
    }else{
        rhs <- data.frame(from=unlist(lhs),to=rep(resultFis$Annot,lapply(lhs,length)))
        rownames(resultFis) <- resultFis$Annot
    }
    if(isTRUE(usePadj)){
        rhs$value <- -log10(resultFis[rhs$to,"Padj"])

    }else{
        rhs$value <- -log10(resultFis[rhs$to,"Pvalue"])
    }
    if (writeCyt == TRUE) {
      write.table(rhs, file = cytoscapeFile, sep = "\t",
                  row.names = F, quote = F)
    }
    pvalue1 <- rhs$value
    names(pvalue1) <- rhs$from
    pvalue1 <- unlist(lapply(split(pvalue1,names(pvalue1)),mean))
    pvalue2 <- rhs$value
    names(pvalue2)<-rhs$to
    pvalue2 <- unlist(lapply(split(pvalue2,names(pvalue2)),mean))
    pvalue<-c(pvalue1,pvalue2)
    g <- graph_from_data_frame(rhs, directed = F)
    pvalue = pvalue[V(g)$name]
    cols <- .color_scale(high, low)
    V(g)$color <- cols[sapply(pvalue, .getIdx, min(pvalue), max(pvalue))]
    if(!is.null(node.shape)){
      node.shape=rep(20,length(V(g)$name))
      names(node.shape)<-V(g)$name
      node.shape[names(node.shape)]<-node.shape
    }else{
      node.shape=rep(20,length(V(g)$name))
      names(node.shape)<-V(g)$name
    }
    p <- ggnet2(g, node.size = degree(g), node.color = V(g)$color,
                node.shape=node.shape,legend.position = "none",
                node.alpha=node.alpha)
    if(isTRUE(repel)){
      p <- p + geom_text_repel(label = V(g)$name,size=label.size, color=label.color,
                      segment.size=segment.size)
    }else{
      p <- p +geom_text(label = V(g)$name,size=label.size, color=label.color)
    }
    print(p)
    if(savefig==TRUE){
      ggsave(p,file=paste(filename,"pdf",sep="."),width=width,height = height)
    }
}
