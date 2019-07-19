##' Download database from Msigdb and prepare for enrichment analysis
##' @name msigdbr
##' @importFrom msigdbr msigdbr
##' @importFrom dplyr filter_
##' @importFrom dplyr select_
##' @importFrom magrittr %>%
##' @param species the species for query
##' @param keytype the gene ID type
##' @param category Gene set category
##' @param subcategory Gene Set subcategory
##' @param save save the dataset or not
##' @param path path to save the dataset
##' @export
##' @author Kai Guo
makeMSIGDB<-function(species="human",keytype="SYMBOL",subcategory=NULL,
                     category=NULL,path="./",save=FALSE){
  flag=0
  if(!is.null(subcategory)){
  if(subcategory=="CGP"){
    category<-"C2"
  }
  if(subcategory=="CP"){
    category<-"C2"
  }
  if(subcategory=="KEGG"){
    subcategory<-"CP:KEGG"
    category<-"C2"
  }
  if(subcategory=="REACTOME"){
    subcategory<-"CP:REACTOME"
    category<-"C2"
  }
  if(subcategory=="BIOCARTA"){
    subcategory<-"CP:BIOCARTA"
    category<-"C2"
  }
  if(subcategory=="MIR"){
    category<-"C3"
  }
  if(subcategory=="TFT"){
    category<-"C3"
  }
  if(subcategory=="CGN"){
    category=="C4"
  }
  if(subcategory=="CM"){
    category<-"C4"
  }
  if(subcategory=="BP"){
    category<-"C5"
  }
  if(subcategory=="CC"){
    category<-"C5"
  }
  if(subcategory=="MF"){
    category<-"C5"
  }
  }
  mspe<-.getmsig(species)
  if(is.null(mspe)){
    stop(cat("can't find support species!\n"))
  }
  if(keytype=="SYMBOL"){
    key="gene_symbol"
  }else if(keytype=="ENTREZID"){
    key="entrez_gene"
  }else{
    key="entrez_gene"
    flag=1
  }
  cat("Downloading msigdb datasets ...\n")
  res <- msigdbr(species=mspe)
  res <- res%>%filter_(~gs_cat==category)
  if(!is.null(subcategory)){
    res <- res%>%filter_(~gs_subcat==subcategory)
  }
  res<-res%>%select_(~key,~gs_name)
  if(flag==1){
    res[,1]<-idconvert(species,keys=res[,1],fkeytype="ENTREZID", tkeytype=keytype)
    res<-na.omit(res)
  }
  if(isTRUE(save)){
    write.csv(res%>%select_(~key,~gs_name,~gs_cat,~gs_subcat),
              file=paste(path,species,"_msigdb.csv",sep=""))
  }
  return(res)
}


##' msigdb support species
##' @param species with common name
.getmsig<-function(species="human"){
  out<-NULL
  if(species=="human"){
    out<-"Homo sapiens"
  }else if(species=="mouse"){
    out<-"Mus musculus"
  }else if(species=="rat"){
    out<-"Rattus norvegicus"
  }else if(species=="celegans"){
    out<-"Caenorhabditis elegans"
  }else if(species=="fly"){
    out<-"rosophila melanogaster"
  }else if(species=="yeast"){
    out<-"Saccharomyces cerevisiae"
  }else if(species=="bovine"){
    out<-"Bos taurus"
  }else if(species=="canine"){
    out<-"Canis lupus familiaris"
  }else if(species=="pig"){
    out<-"Sus scrofa"
  }else if(species=="chicken"){
    out<-"Gallus gallus"
  }else if(species=="zebrafish"){
    out<-"Danio rerio"
  }else{
    out<-NULL
  }
}
##' Print MSIGDB infomation
##' @export
msigdbinfo <- function() {
  cat("#--------------------------------------------------------------#\n")
  cat("# Molecular Signatures Database                        v6.2.1  #\n")
  cat("#--------------------------------------------------------------#\n")
  cat("# Category | Subcategory # Details ----------------------------#\n")
  cat("# C1               # Positional (326)                          #\n")
  cat("# C2 | CGP         # Chemical and Genetic Perturbations (3433) #\n")
  cat("# C2 | CP          # Canonical Pathways (252)                  #\n")
  cat("# C2 | BIOCARTA # Canonical BIOCARTA (217)                     #\n")
  cat("# C2 | KEGG     # Canonical KEGG (186)                         #\n")
  cat("# C2 | CPREACTOME  # Canonical REACTOME (674)                  #\n")
  cat("# C3 | MIR         # Motif miRNA Targets (221)                 #\n")
  cat("# C3.| TFT         # Motif Transcription Factor Targets (615)  #\n")
  cat("# C4.| CGN         # Cancer Gene Neighborhoods (427)           #\n")
  cat("# C4.| CM          # Cancer Modules (431)                      #\n")
  cat("# C5.| BP          # GO Biological Process (4436)              #\n")
  cat("# C5.| CC          # GO Cellular Component (580)               #\n")
  cat("# C5.| MF          # GO Molecular Function (901)               #\n")
  cat("# C6               # Oncogenic Signatures (189)                #\n")
  cat("# C7               # Immunologic Signatures (4872)             #\n")
  cat("# H                # Hallmark (50)                             #\n")
  cat("#--------------------------------------------------------------#\n")
  cat("# Source: http://software.broadinstitute.org/gsea/msigdb       #\n")
  cat("#--------------------------------------------------------------#\n")
  listspe<-c("human","mouse","rat","celegans","fly","yeast","bovine","canine",
             "pig","chicken","zebrafish")
  cat("# Support species:                                             #\n")
  cat(sort(listspe),"\n")
}
