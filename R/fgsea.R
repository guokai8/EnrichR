#' Enrichment analysis for any type of annotation data
#' @param x a vector include all log2FC with gene name
#' @param annot annotation file for all genes
#' @param annot.info Term of all annotation
#' @param filenam output filename
#' @param padj.method p value adjust method
#' @param nperm Number of permutations to do. Minimial possible nominal p-value is about 1/nperm
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @importFrom fgsea fgsea plotEnrichment plotGseaTable
#' @export
#' @author Kai Guo
enrich_fgsea<-function(x,annot,annot.info=NULL,minSize=15,maxSize=500,nperm=5000,filename=NULL,padj.method="BH"){
  suppressMessages(require(fgsea))
  x<-sort(x)
  annod<-sf(annot);
  res<-fgsea(pathways=annod,stats=x,minSize=minSize,maxSize=maxSize,nperm=nperm)
  return(res)
}
#' @name plot_gsea
#' @title plot the gsea result
#' @param x a vector include all log2FC with gene name
#' @param annot annotation file for all genes
#' @param term the significant term
#' @importFrom fgsea fgsea plotEnrichment plotGseaTable
#' @export
#' @author Kai Guo
plot_fgsea<-function(term,x,annot){
  x<-sort(x)
  plotEnrichment(sf(annot)[[term]],x)
}
