#' Get detail information from enrichment results
#' @param rese Functional enrichment results
#' @param resd DEG result with one gene name column or vector of input genes 
#' @export
#' @author Kai Guo
getdetail<-function(rese,resd){
  require(tidyverse)
  if(!is.data.frame(resd)){
    resd=data.frame(gene=resd)
    }
  gene<-strsplit(as.vector(rese$GeneID),split="\\,")
  names(gene)<-rese$Annot
  gened<-data.frame("TERM"=rep(names(gene),times=unlist(lapply(gene,length))),"Annot"=rep(rese$Term,times=unlist(lapply(gene,length))),"GeneID"=unlist(gene),row.names=NULL)
  gened$GeneID<-as.character(gened$GeneID)
  res<-left_join(gened,resd,by=c("GeneID"="gene"))
  return(res)
}
