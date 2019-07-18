#' Enrichment analysis for any type of annotation data
#' @param x a vector include all log2FC with gene name
#' @param annot annotation file for all genes
#' @param annot.info Term with annotation details
#' @param filename output filename
#' @param padj.method p value adjust method
#' @param nperm Number of permutations to do. Minimial possible nominal p-value is about 1/nperm
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @importFrom fgsea fgsea
#' @export
#' @author Kai Guo
gsea<-function(x,annot,annot.info=NULL,minSize=15,maxSize=500,nperm=5000,filename=NULL,padj.method="BH",table = TRUE){
  x<-sort(x)
  if(!is.null(annot.info)){
    rownames(annot.info)<-annot.info[,1]
    annot.info$annot<-paste(rownames(annot.info),annot.info[,1],sep="_")
    annot[,2]<-annot.info[annot[,2],"Annot"]
  }else{
    if(!is.null(annot$Annot)){
      annot[,2]<-annot$Annot
    }
  }
  annod<-sf(annot);
  res<-fgsea(pathways=annod,stats=x,minSize=minSize,maxSize=maxSize,nperm=nperm)
  if(isTRUE(table)){
    res$leadingEdge<-unlist(lapply(res$leadingEdge, function(x)paste(x,collapse = ",",sep="")))
  }
  if(!is.null(filename)){
    tmp<-res
    tmp$leadingEdge<-unlist(lapply(tmp$leadingEdge, function(x)paste(x,collapse = ",",sep="")))
    write.csv(tmp,filename)
  }
  return(res)
}
#' @name plotgsea
#' @title plot the gsea result
#' @param x a vector include all log2FC with gene name
#' @param annot annotation file for all genes
#' @param term the significant term
#' @param annot.info Term with annotation details
#' @param
#' @importFrom fgsea plotEnrichment plotGseaTable
#' @export
#' @author Kai Guo
plotgsea<-function(x,term,annot,annot.info=NULL,gseaRes=NULL){
  x<-sort(x)
  if(!is.null(annot.info)){
    rownames(annot.info)<-annot.info[,1]
    annot.info$annot<-paste(annot.info[,1],annot.info[,2],sep="_")
    annot[,2]<-annot.info[annot[,2],"Annot"]
  }else{
    if(!is.null(annot$Annot)){
      annot[,2]<-annot$Annot
    }
  }
  annod <- sf(annot)
  if(length(term)>1&!is.null(gseaRes)){
    plotGseaTable(annod[term],stats=x,gseaRes,gseaParam=0.5)
  }else{
    plotEnrichment(annod[[term]],stats=x)
  }
}
