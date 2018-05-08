#' Get detail information from enrichment results
#' @importFrom tidyverse left_join
#' @param rese Functional enrichment results
#' @param resd DEG result or vector of input genes
#' @export
#' @author Kai Guo
getdetail<-function(rese,resd){
  if(!is.data.frame(resd)){
    resd=data.frame(gene=resd)
    }
  gene<-strsplit(rese$GeneID,split=",")
  names(gene)<-rese$Annot
  gened<-data.frame("TERM"=rep(names(gene),times=unlist(lapply(gene,length))),"Annot"=rep(rese$Term,times=unlist(lapply(gene,length))),"GeneID"=unlist(gene),row.names=NULL)
  gened$GeneID<-as.character(gened$GeneID)
  res<-left_join(gened,resd,by=c("GeneID"="gene"))
  return(res)
}
