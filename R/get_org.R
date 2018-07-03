#' make GO annotation data function
#' @param species: you can check the support species by using showData()
#' @param keytype: the gene type
#' @export
#' @author Kai Guo
makeGOdat<-function(species="human",keytype="ENTREZID"){
  dbname<-.getdbname(species);
  suppressMessages(require(AnnotationDbi))
  sel<-AnnotationDbi::select
  if (!require(dbname,character.only=TRUE)){
    source("http://bioconductor.org/biocLite.R")
    biocLite(dbname)
  }else{
    suppressMessages(require(dbname,character.only = T,quietly = T))
  }
  dbname<-eval(parse(text=dbname))
  GO_FILE=sel(dbname,keys=keys(dbname,keytype=keytype),keytype=keytype,columns=c("GOALL","ONTOLOGYALL"))
  return(GO_FILE)
}
#' make KEGG annotation data function
#' @param species: you can check the support species by using showData()
#' @param keytype: the gene type
#' @export
#' @author Kai Guo
makeKOdat<-function(species="human",keytype="ENTREZID",builtin=TRUE){
  dbname<-.getdbname(species=species);
  if(builtin==TRUE){
  suppressMessages(require(AnnotationDbi))
  sel<-AnnotationDbi::select
  if (!require(dbname,character.only=TRUE)){
    source("http://bioconductor.org/biocLite.R")
    biocLite(dbname)
  }else{
    suppressMessages(require(dbname,character.only = T,quietly = T))
  }
  dbname<-eval(parse(text=dbname))
  KO_FILE=sel(dbname,keys=keys(dbname,keytype=keytype),keytype=keytype,columns="PATH")
  KO_FILE<-na.omit(KO_FILE)
  }else{
    suppressMessages(require(KEGGREST))
    spe=.getspeices(species)
    tmp<-keggLink("pathway",spe)
    tmp<-substr(tmp,9,13)
    names(tmp)<-sub('.*:','',names(tmp))
    tmp<-vec_to_df(tmp,name=c(keytype,"PATH"))
    if(keytype!="ENTREZID"){
      tmp[,1]<-idconvert(species,keys=tmp[,1],fkeytype = "ENTREZID",tkeytype = keytype)
      tmp<-na.omit(tmp)
    }
    KO_FILE=tmp
  }
  return(KO_FILE)
}
#' Convert ID between ENTREZID to SYMBOL or other type ID based on bioconductor annotation package
#' @param species: you can check the support species by using showData()
#' @param fkeytype: the gene type you want to convert
#' @param tkeytype: the gene type you want to get
#' @export
#' @author Kai Guo
idconvert<-function(species,keys,fkeytype,tkeytype){
  dbname<-.getdbname(species);
  suppressMessages(require(dbname,character.only = T))
  dbname<-eval(parse(text=dbname))
  mapIds(dbname,keys=as.vector(keys),
                  column=tkeytype,
                  keytype=fkeytype,
                  multiVals="first")
}
.getdbname<-function(species="human"){
    dbname=.getdb(species=species);
  if(is.null(dbname)){
    cat("You must check if your request database is avaliable by using showData,
        If not you could make your database by using makeOwnGO and makeOwnKO
        and give a user defined database\n")
    stop("databse error!")
  }
    return(dbname)
}
.getdb<-function(species=species){
  species=tryCatch(match.arg(species,c("anopheles","arabidopsis","bovine","celegans","canine","fly","zebrafish",
                                       "ecoli","ecsakai","chicken","human","mouse","rhesus","malaria","chipm","rat",
                                       "toxoplasma","streptomyces","pig","yeast","xenopus","warm")),
           error=function(cond){return("unsupported")})
  if (species == "anopheles") {
    dbname <- "org.Ag.eg.db"
  } else if (species == "arabidopsis") {
    dbname <- "org.At.tair.db"
  } else if (species == "bovine") {
    dbname <- "org.Bt.eg.db"
  } else if (species == "canine") {
    dbname <- "org.Cf.eg.db"
  } else if (species == "worm" || species == "celegans") {
    dbname <- "org.Ce.eg.db"
  } else if (species == "chicken") {
    dbname <- "org.Gg.eg.db"
  } else if (species == "ecolik12") {
    dbname <- "org.EcK12.eg.db"
  } else if (species == "ecsakai") {
    dbname <- "org.EcSakai.eg.db"
  } else if (species == "fly") {
    dbname <- "org.Dm.eg.db"
  } else if (species == "human") {
    dbname <- "org.Hs.eg.db"
  } else if (species == "malaria") {
    dbname <- "org.Pf.plasmo.db"
  } else if (species == "chipm") {
    dbname <- "org.Pt.eg.db"
  }else if (species == "mouse") {
    dbname <- "org.Mm.eg.db"
  } else if (species == "pig") {
    dbname <- "org.Ss.eg.db"
  } else if (species == "rat") {
    dbname <- "org.Rn.eg.db"
  } else if (species == "rhesus") {
    dbname <- "org.Mmu.eg.db"
  } else if (species == "xenopus") {
    dbname <- "org.Xl.eg.db"
  } else if (species == "yeast") {
    dbname <- "org.Sc.sgd.db"
  } else if (species == "streptomyces") {
    dbname <- "org.Sco.eg.db"
  } else if (species == "zebrafish") {
    dbname <- "org.Dr.eg.db"
  } else if (species == "toxoplasma"){
    dbname<- "org.Tgondii.eg.db"
  } else {
    dbname <- NULL
  }
  return(dbname)
}
.getspeices<-function(species="human"){
  species=tryCatch(match.arg(species,c("anopheles","arabidopsis","bovine","celegans","canine","fly","zebrafish",
                                       "ecoli","ecsakai","chicken","human","mouse","rhesus","malaria","chipm","rat",
                                       "toxoplasma","sco","pig","yeast","xenopus","warm")),
                   error=function(cond){return("unsupported")})
  if (species == "anopheles") {
    species <- "aga"
  } else if (species == "arabidopsis") {
    species <- "ath"
  } else if (species == "bovine") {
    species <- "bta"
  } else if (species == "canine") {
    species <- "cfa"
  } else if (species == "chicken") {
    species <- "gga"
  } else if (species == "chipm") {
    species <- "ptr"
  } else if (species == "ecolik12") {
    species <- "eco"
  } else if (species == "ecsakai") {
    species <- "ecs"
  } else if (species == "fly") {
    species <- "dme"
  } else if (species == "human") {
    species <- "hsa"
  } else if (species == "malaria") {
    species <- "pfa"
  } else if (species == "mouse") {
    species <- "mmu"
  } else if (species == "pig") {
    species <- "ssc"
  } else if (species == "rat") {
    species <- "rno"
  } else if (species == "rhesus") {
    species <- "mcc"
  } else if (species == "worm" || species == "celegans") {
    species <- "cel"
  } else if (species == "xenopus") {
    species <- "xla"
  } else if (species == "yeast") {
    species <- "sce"
  } else if (species =="streptomyces"){
    species <- "sco"
  } else if (species == "zebrafish") {
    species <- "dre"
  } else {
    species <- NULL
  }
  return(species)
}
#' show avaliable data based on bioconductor annotation package
#' @export
#' @author Kai Guo
showData<-function(){
  species=c("anopheles","arabidopsis","bovine","celegans","canine","fly","zebrafish",
            "ecoli","ecsakai","chicken","human","mouse","rhesus","malaria","chipm","rat",
            "toxoplasma","sco","pig","yeast","xenopus")
  dbname=c("org.Ag.eg.db","org.At.tair.db","org.Bt.eg.db","org.Ce.eg.db","org.Cf.eg.db","org.Dm.eg.db",
           "org.Dr.eg.db","org.EcK12.eg.db","org.EcSakai.eg.db","org.Gg.eg.db","org.Hs.eg.db","org.Mm.eg.db",
           "org.Mmu.eg.db","org.Pf.plasmo.db","org.Pt.eg.db","org.Rn.eg.db","org.Sc.sgd.db","org.Sco.eg.db",
           "org.Ss.eg.db","org.Tgondii.eg.db","org.Xl.eg.db")
  dbdata<-data.frame(dbname=dbname,species=species)
  dbdata
}

vec_to_df<-function(x,name){
  dd<-data.frame(names(x),x)
  colnames(dd)<-name
  return(dd)
}
