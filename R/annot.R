##' @importFrom dplyr filter_
##' @importFrom AnnotationDbi select
.get_go_dat<-function(ont="BP"){
  require(GO.db)
  key<-keys(GO.db)
  suppressMessages(go_dat<-select(GO.db, keys=key, columns=c("TERM","ONTOLOGY"),keytype="GOID"))
  if(ont=="BP") res<-as.data.frame(filter_(go_dat,~ONTOLOGY=="BP"))
  if(ont=="CC") res<-as.data.frame(filter_(go_dat,~ONTOLOGY=="CC"))
  if(ont=="MF") res<-as.data.frame(filter_(go_dat,~ONTOLOGY=="MF"))
  rownames(res)<-res[,1]
  res<-res[, 2, drop = FALSE]
  colnames(res)<-"annotation"
  return(res)
}
##' @importFrom KEGGREST keggList
.get_kg_dat<-function(builtin=TRUE){
   if(isTRUE(builtin)){
      data(kegg)
      return(kegg.db)
     }else{
     pathway<-cbind(keggList('pathway'))
     rownames(pathway)<-sub('path:map','',rownames(pathway))
     colnames(pathway)<-"annotation"
     pathway<-as.data.frame(pathway)
     pathway$annotation<-as.vector(pathway$annotation)
     return(pathway)
     }
}
##' build annotaion for kegg
##' @param ontype GO or KEGG
##' @examples
##' annot = getann("GO")
##' @export
##' @author Kai Guo
getann<-function(ontype="GO"){
    if(ontype=="GO"){
      res<-rbind(.get_go_dat("BP"),.get_go_dat("MF"),.get_go_dat("CC"))
    }
    if(ontype=="KEGG"){
        res<-.get_kg_dat(builtin=F)
    }
    return(res)
}
