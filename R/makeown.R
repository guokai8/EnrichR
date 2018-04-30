#' Make user define annotation data
#' @param anno: data.frame includes gene_id and GO or Pathway information
#' @param godata: godata inlude BP,MF,CC annotation
#' @param use: check if the user choose to load internal godata
#' @export
#' @author Kai Guo
makeOwnGO<-function(anno,godata=godata,use=TRUE){
  if(use==TRUE){
    data(godata)
    godata<-godata
  }else{
    #id<-keys(GO.db)
    #MF = "GO:0003674", node of MF
    #BP = "GO:0008150", node of BP
    #CC = "GO:0005575", node of CC
    gobp<-sapply(c("GO:0008150",GOBPOFFSPRING[["GO:0008150"]]),function(x)GO_child(x,ontology = "BP"))
    gocc<-sapply(c("GO:0005575",GOBPOFFSPRING[["GO:0005575"]]),function(x)GO_child(x,ontology = "CC"))
    gomf<-sapply(c("GO:0003674",GOBPOFFSPRING[["GO:0003674"]]),function(x)GO_child(x,ontology = "MF"))
    godata<-c(gobp,c(gocc,gomf))
  }
  anno.in<-sf(anno)
  anno.bp<-lapply(godata, function(x)as.vector(unlist(anno.in[x])))
  res<-Filter(Negate(is.null), anno.bp)
  res<-data.frame("Gene"=as.vector(unlist(res)),"GO"=rep(names(res),times=lapply(res, length)))
  return(res)
}
makeOwnKO<-function(anno){
  anno<-na.omit(anno)
  return(annot)
}
