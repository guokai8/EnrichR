#' PFAM Enrichment analysis function based on ensemble annotation data
#' @param df: DGE files (DESeq2 result files) or vector contains gene names
#' @param PO_FILE: PFAM annotation data, you can get the data by using makeesanno function
#' @param filename: output filename
#' @param gene.cutoff: DGE singificant cutoff value
#' @export
#' @author Kai Guo
PE.ensemble<-function(df,PO_FILE,gene.cutoff=0.01,padj.method="BH",minSize=2,maxSize=500,keepRich=TRUE,filename=NULL){
  data(pfanno)
  annot.info=pfanno
  res=enrich(df,annot=PO_FILE,annot.info =pfanno,filename=filename,gene.cutoff = gene.cutoff,padj.method = padj.method)
  return(res)
}
#' PFAM Enrichment analysis function based on ensemble plant annotation data
#' @param df: DGE files (DESeq2 result files) or vector contains gene names
#' @param PO_FILE: PFAM annotation data, you can get the data by using makeplantann function
#' @param filename: output filename
#' @param gene.cutoff: DGE singificant cutoff value
#' @export
#' @author Kai Guo
PE.plant<-function(df,PO_FILE,gene.cutoff=0.01,padj.method="BH",minSize=2,maxSize=500,keepRich=TRUE,filename=NULL){
  data(pfanno)
  res=enrich(df,annot=PO_FILE,annot.info = pfanno,filename=filename,gene.cutoff = gene.cutoff,padj.method = padj.method,minSize=minSize,maxSize=maxSize,keepRich=keepRich)
  return(res)
}
