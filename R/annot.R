.get_go_dat<-function(ont="BP"){
  require(GO.db)
  suppressMessages(library(dplyr))
  key<-keys(GO.db)
  suppressMessages(go_dat<-AnnotationDbi::select(GO.db, keys=key, columns=c("TERM","ONTOLOGY"),keytype="GOID"))
  if(ont=="BP") res<-as.data.frame(dplyr::filter(go_dat,ONTOLOGY=="BP"))
  if(ont=="CC") res<-as.data.frame(dplyr::filter(go_dat,ONTOLOGY=="CC"))
  if(ont=="MF") res<-as.data.frame(dplyr::filter(go_dat,ONTOLOGY=="MF"))
  return(res)
}
.get_kg_dat<-function(builtin=TRUE){
   if(isTRUE(builtin)){
      data(kegg)
      return(kegg.db)
     }else{
     require(KEGGREST)
     pathway<-cbind(keggList('pathway'))
     rownames(pathway)<-sub('path:map','',rownames(pathway))
     colnames(pathway)<-"annotation"
     pathway<-as.data.frame(pathway)
     return(pathway)
     }
}
