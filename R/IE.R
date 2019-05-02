#' InterPro Enrichment analysis function
#' @param df DGE files (DESeq2 result files) or vector contains gene names
#' @param IO_FILE InterPro annotation data
#' @param minSize minimal number of genes included in significant terms
#' @param padj.method p value adjust method (default: BH)
#' @param cutoff cutoff value for filtering significant terms (default: 0.05)
#' @param maxSize maximum number of genes included in significant terms
#' @param keepRich keep terms with rich factor value equal 1 or not (default: TRUE)
#' @param filename output filename
#' @param gene.cutoff the cut-off value for select DEGs (default: 0.01)
#' @export
#' @author Kai Guo
IE.ensemble<-function(df,IO_FILE,gene.cutoff=0.01,minSize=2,maxSize=500,keepRich=TRUE,padj.method="BH",filename=NULL,cutoff=0.05){
  data(interproanno)
  rownames(interano)<-interano[,1]
  res=enrich(df,annot=IO_FILE,annot.info = interano,filename=filename,gene.cutoff = gene.cutoff,padj.method = padj.method,minSize=minSize,maxSize=maxSize,keepRich=keepRich)
  res<-res[res$Pvalue<cutoff,]
  return(res)
}
#' InterPro Enrichment analysis function for plant
#' @param df DGE files (DESeq2 result files) or vector contains gene names
#' @param IO_FILE InterPro annotation data
#' @param minSize minimal number of genes included in significant terms
#' @param padj.method p value adjust method (default: BH)
#' @param cutoff cutoff value for filtering significant terms (default: 0.05)
#' @param maxSize maximum number of genes included in significant terms
#' @param keepRich keep terms with rich factor value equal 1 or not (default: TRUE)
#' @param filename output filename
#' @param gene.cutoff the cut-off value for select DEGs (default: 0.01)
#' @export
#' @author Kai Guo
IE.plant<-function(df,IO_FILE,gene.cutoff=0.01,minSize=2,maxSize=500,keepRich=TRUE,padj.method="BH",filename=NULL,cutoff=0.05){
  data(interproanno)
  rownames(interano)<-interano[,1]
  res=enrich(df,annot=IO_FILE,annot.info = interano,filename=filename,gene.cutoff = gene.cutoff,padj.method = padj.method,minSize=minSize,maxSize=maxSize,keepRich=keepRich)
  res<-res[res$Pvalue<cutoff,]
  return(res)
}
