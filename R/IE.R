#' InterPro Enrichment analysis function
#' @param df: DGE files (DESeq2 result files) or vector contains gene names
#' @param IO_FILE: InterPro annotation data, you can get the data by using makeesanno function
#' @param filename: output filename
#' @param gene.cutoff: DGE singificant cutoff value
#' @export
#' @author Kai Guo
IE.ensemble<-function(df,IO_FILE,gene.cutoff=0.01,minSize=2,maxSize=500,keepRich=TRUE,padj.method="BH",filename=NULL){
  data(interproanno)
  rownames(interano)<-interano[,1]
  res=enrich(df,annot=IO_FILE,annot.info = interano,filename=filename,gene.cutoff = gene.cutoff,padj.method = padj.method,minSize=minSize,maxSize=maxSize,keepRich=keepRich)
  return(res)
}
#' InterPro Enrichment analysis function for plant
#' @param df: DGE files (DESeq2 result files) or vector contains gene names
#' @param IO_FILE: InterPro annotation data, you can get the data by using makeplantann function
#' @param filename: output filename
#' @param gene.cutoff: DGE singificant cutoff value
#' @export
#' @author Kai Guo
IE.plant<-function(df,IO_FILE,gene.cutoff=0.01,minSize=2,maxSize=500,keepRich=TRUE,padj.method="BH",filename=NULL){
  data(interproanno)
  rownames(interano)<-interano[,1]
  res=enrich(df,annot=IO_FILE,annot.info = interano,filename=filename,gene.cutoff = gene.cutoff,padj.method = padj.method,minSize=minSize,maxSize=maxSize,keepRich=keepRich)
  return(res)
}
