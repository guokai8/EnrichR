#' Show avaliable plant annotation based on ENSEMBLE
#' @param host the ensemble API host,for plant you can use plants.ensembl.org and for human and other species you can use uswest.ensembl.org
#' @param species the species you want to search
#' @param ann_type the type of function annotation you want get from ensemble
#' @export
#' @author Kai Guo
showplant<-function(host="plants.ensembl.org"){
  suppressMessages(require(biomaRt))
  mart=useMart("plants_mart",host="plants.ensembl.org")
  res<-listDatasets(mart)
  colnames(res)[2]<-"species"
  res
}
#' Show avaliable annotation except plant based on ENSEMBLE
#' @importFrom biomaRt useMart
#' @importFrom biomaRt listDatasets
#' @param host the ensemble API host,for plant you can use plants.ensembl.org and for human and other species you can use uswest.ensembl.org
#' @param species the sepcies you want to search
#' @export
#' @author Kai Guo
showensemble<-function(host="uswest.ensembl.org"){
  cat("You could choose different host to get high speed!\n")
  cat("Ensembl US West: uswest.ensembl.org\nEnsembl US East: useast.ensembl.org\nEnsembl Asia: asia.ensembl.org\n" )
  mart=useMart("ENSEMBL_MART_ENSEMBL",host=host)
  res<-listDatasets(mart)
  colnames(res)[2]<-"species"
  res
}
#' make plant annotation data based on ENSEMBLE
#' @importFrom biomaRt useMart
#' @importFrom biomaRt useDataset
#' @importFrom biomaRt getBM
#' @param host the ensemble API host,for plant you can use plants.ensembl.org and for human and other species you can use uswest.ensembl.org
#' @param species the sepcies you want to search, you can use showplant to get the species name
#' @param ann_type the type of function annotation(GO,KEGG,PFAM,InterPro) you want get from ensemble
#' @export
#' @author Kai Guo
makeplantann<-function(species="Arabidopsis t",host="plants.ensembl.org",ann_type="GO"){
  mart=useMart("plants_mart",host=host)
  dbinfo<-.getmartdb(species,mart)
  dbname=as.character(dbinfo$dbname)
  dataset<-useDataset(dbname,mart=mart)
  if(ann_type=="GO"){
   res<-getBM(attributes = c("ensembl_gene_id","go_id","name_1006"),filters ="chromosome_name",values = as.vector(dbinfo$chr_info$name),dataset)
  }else if(ann_type=="KEGG"){
   res<-getBM(attributes = c("ensembl_gene_id","kegg_enzyme"),filters ="chromosome_name",values = as.vector(dbinfo$chr_info$name),dataset)
   res[,2]<-sub('\\+.*','',res[,2])
   #res$Annot<-kegg.db[res[,2],]
   #res<-res[,c(1,3,2)]
  }else if(ann_type=="PFAM"){
   res<-getBM(attributes = c("ensembl_gene_id","pfam"),filters ="chromosome_name",values = as.vector(dbinfo$chr_info$name),dataset)
  }else if(ann_type=="InterPro"){
    res<-getBM(attributes = c("ensembl_gene_id","interpro"),filters ="chromosome_name",values = as.vector(dbinfo$chr_info$name),dataset)
  }else if(ann_type=="InterPro"){
    res<-getBM(attributes = c("ensembl_gene_id","gramene_plant_reactome"),filters ="chromosome_name",values = as.vector(dbinfo$chr_info$name),dataset)
  }else{
    stop("You must specify one type of annotation!\n")
  }
  res<-res[nchar(res[,2])>1,]
  if(ann_type=="KEGG"){
    res$Annot<-kegg.db[res[,2],]
   # res<-res[,c(1,3,2)]
  }
  colnames(res)[ncol(res)]<-"Annot"
  rownames(res)<-NULL
  return(res)
}
#' make annotation data except plant based on ENSEMBLE
#' @importFrom biomaRt useMart
#' @importFrom biomaRt useDataset
#' @importFrom biomaRt getBM
#' @param host the ensemble API host,for plant you can use plants.ensembl.org and for human and other species you can use uswest.ensembl.org
#' @param species the species you want to search, you can use showplant to get the species name
#' @param ann_type the type of function annotation(GO,KEGG,PFAM,InterPro) you want get from ensemble
#' @export
#' @author Kai Guo
makeesanno<-function(species="Human",host="uswest.ensembl.org",ann_type="GO"){
  mart=useMart("ENSEMBL_MART_ENSEMBL",host=host)
  dbinfo<-.getmartdb(species,mart)
  dbname=as.character(dbinfo$dbname)
  dataset<-useDataset(dbname,mart=mart)
  if(ann_type=="GO"){
    res<-getBM(attributes = c("ensembl_gene_id","go_id","name_1006"),filters ="chromosome_name",values = as.vector(dbinfo$chr_info$name),dataset)
  }else if(ann_type=="KEGG"){
    res<-getBM(attributes = c("ensembl_gene_id","kegg_enzyme"),filters ="chromosome_name",values = as.vector(dbinfo$chr_info$name),dataset)
    res[,2]<-sub('\\+.*','',res[,2])
  }else if(ann_type=="PFAM"){
    res<-getBM(attributes = c("ensembl_gene_id","pfam"),filters ="chromosome_name",values = as.vector(dbinfo$chr_info$name),dataset)
  }else if(ann_type=="InterPro"){
    res<-getBM(attributes = c("ensembl_gene_id","interpro"),filters ="chromosome_name",values = as.vector(dbinfo$chr_info$name),dataset)
  }else{
    stop("You must specify one type of annotation!\n")
  }
  res<-res[nchar(res[,2])>1,]
  if(ann_type=="KEGG"){
    res$Annot<-kegg.db[res[,2],]
    # res<-res[,c(1,3,2)]
  }
  colnames(res)[ncol(res)]<-"Annot"
  rownames(res)<-NULL
  return(res)
}
##' @importFrom dplyr select_
##' @importFrom dplyr collect
##' @importFrom dplyr pull
##' @importFrom magrittr %>%
##' @importFrom biomaRt listDatasets
##' @importFrom jsonlite fromJSON
.getmartdb<-function(spe,mart){
  lhs<- listDatasets(mart)
  spe=simpleCap(spe);
  spe<-gsub(' ','\\\\s',spe)
  sel<-grepl(spe,lhs$description,ignore.case = T)
  tmp<-lhs[sel,]
  dataset<-tmp%>%select_(~dataset)%>%collect%>%pull(1)
  if((length(dataset)==0)|(length(dataset)>1)){
    stop("Maybe you need first check the avaliable database by using showensemble() or showplant()\n")
  }
  chr<-tmp%>%select_(~description)%>%collect%>%pull(1)
  organism<-gsub(' ','_',sub(' genes.*','',chr))
  if(organism=="Oryza_sativa_Japonica"){
    organism="Oryza_sativa"
  }
  if(mart@biomart=="plants_mart"){
    pre_site="http://rest.ensembl.org/info/assembly/"
  }else{
    pre_site="http://rest.ensembl.org/info/assembly/"
  }
  tryCatch({
    chr_d <-
      fromJSON(
        paste0(
          pre_site,
          organism,
          "?content-type=application/json"
        )
      )
  }, error = function(e)
    stop(
      "The API 'http://rest.ensemblgenomes.org' does not seem to work properly. Are you connected to the internet? Is the homepage 'http://rest.ensemblgenomes.org' currently available?", call. = FALSE
    ))
  chr_info<-chr_d$top_level_region
  chr_version<-chr_d$assembly_name
  chr_assembly_date<-chr_d$assembly_date
  rhs<-list(dbname=dataset,chr_info=chr_info,chr_version=chr_version,chr_assembly_date=chr_assembly_date)
  return(rhs)
}

