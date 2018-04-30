#' make Reactome annotation data function
#' @param species: you can check the supported species by using showAvailableRO and showAvailablePlants
#' @export
#' @author Kai Guo
makeROdata<-function(species="Homo_sapiens"){
   suppressMessages(require(AnnotationDbi))
   dbname<-.getrodbname(species=species);
   if (!require("reactome.db",character.only=TRUE)){
     source("http://bioconductor.org/biocLite.R")
     biocLite("reactome.db")
   }else{
     suppressMessages(require("reactome.db",character.only = T,quietly = T))
   }
   dbname=sapply(strsplit(dbname,"_"),'[[',1)
   lhs<-AnnotationDbi::as.list(reactomePATHNAME2ID)
   lhs<-lhs[grep(dbname,names(lhs))]
   roid<-AnnotationDbi::as.list(reactomePATHID2EXTID)[unique(as.vector(unlist(lhs)))]
   roid<-lapply(roid, function(x)unique(x))
   roid<-data.frame("GeneID"=unlist(roid),"RO"=rep(names(roid),times=lapply(roid, length)),row.names=NULL)
   ll<-lapply(lhs,function(x)unique(x))
   roan<-data.frame("RO"=unlist(ll),"Description"=rep(names(ll),times=lapply(ll,length)),row.names=NULL)
   res<-list()
   res$ro<-roid
   res$roan<-roan
   return(res)
}
#' make Plant Reactome annotation data function
#' @param species: you can check the supported species by using showAvailableRO and showAvailablePlants
#' @export
#' @author Kai Guo
makeplantROdat<-function(species="Arabidopsis_thaliana"){
  suppressMessages(library(dplyr));
  data(rodata)
  dbname<-.getplantrodbname(species=species);
  RO_FILE=rodata%>%dplyr::filter(species==dbname)%>%dplyr::select(GeneID,RO,Description,species)
  return(as.data.frame(RO_FILE));
}
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}
.getplantrodbname<-function(species="Arabidopsis_thaliana"){
  species=simpleCap(species);
  dbname=.getplantrodb(species=species);
  if(is.null(dbname)){
    cat("You must check if your request database is avaliable by using showAvailablePlants()\n")
    stop("databse error!")
  }
  return(dbname)
}
.getrodbname<-function(species){
  species=simpleCap(species);
  dbname=.getrodb(species=species);
  if(is.null(dbname)){
    cat("You must check if your request database is avaliable by using showAvailableRO()\n")
    stop("databse error!")
  }
  return(dbname)
}
.getrodb<-function(species){
  dbname=tryCatch(match.arg(species,c("Homo_sapiens","Dictyostelium_discoideum","Plasmodium_falciparum","Schizosaccharomyces_pombe","Saccharomyces_cerevisiae",
                 "Caenorhabditis_elegans","Sus_scrofa","Bos_taurus","Canis_familiaris","Mus_musculus","Rattus_norvegicus","Taeniopygia_guttata",
                 "Xenopus_tropicalis","Danio_rerio","Drosophila_melanogaster","Arabidopsis_thaliana","Oryza_sativa","Gallus_gallus",
                 "Mycobacterium_tuberculosis")),
         error=function(cond){return(NULL)})
  return(dbname)
}
.getplantrodb<-function(species="Arabidopsis_thaliana"){
  speices=simpleCap(species)
  dbname=tryCatch(match.arg(species,c("Arachis_duranensis","Arachis_ipaensis","Arabidopsis_lyrata","Aegilops_tauschii","Arabidopsis_thaliana",
                                       "Amborella_trichopoda","Brachypodium_distachyon","Brassica_napus","Brassica_oleracea","Brassica_rapa",
                                       "Beta_vulgaris","Capsicum_annuum","Cicer_arietinum","Cajanus_cajan","Coffea_canephora","Cyanidioschyzon_merolae",
                                       "Chlamydomonas_reinhardtii","Citrus_sinensis","Eucalyptus_grandis","Fragaria_vesca","Glycine_max","Gossypium_raimondii",
                                       "Hordeum_vulgare","Jatropha_curcas","Leersia_perrieri","Musa_acuminata","Malus_domestica","Manihot_esculenta","Mimulus_guttatus",
                                       "Medicago_truncatula","Oryza_australiensis","Oryza_barthii","Oryza_brachyantha","Oryza_glaberrima","Oryza_granulata",
                                       "Oryza_glumaepatula","Oryza_kasalath","Oryza_longistaminata","Ostreococcus_lucimarinus","Oryza_meridionalis","Oryza_minuta",
                                       "Oryza_nivara","Oryza_officinalis","Oryza_punctata","Oryza_rufipogon","Oryza_sativa","Oryza_sativa_Indica","Picea_abies",
                                       "Phoenix_dactylifera","Physcomitrella_patens","Prunus_persica","Pinus_taeda","Populus_trichocarpa","Phaseolus_vulgaris",
                                       "Sorghum_bicolor","Setaria_italica","Solanum_lycopersicum","Selaginella_moellendorffii","Synechocystis_pcc6803","Solanum_tuberosum",
                                       "Triticum_aestivum","Theobroma_cacao","Trifolium_pratense","Triticum_turgidum","Triticum_urartu","Vitis_vinifera","Zea_mays")),
                   error=function(cond){return(NULL)})
  return(dbname)
}
#' show avaliable plant Reactome annotation data function
#' @export
#' @author Kai Guo
showAvailableplantsRO<-function(){
  cbind(unique(rodata$species))
}
#' show avaliable Reactome annotation data except plant
#' @export
#' @author Kai Guo
showAvailableRO<-function(){
  cbind(c("Homo_sapiens","Dictyostelium_discoideum","Plasmodium_falciparum","Schizosaccharomyces_pombe","Saccharomyces_cerevisiae",
          "Caenorhabditis_elegans","Sus_scrofa","Bos_taurus","Canis_familiaris","Mus_musculus","Rattus_norvegicus","Taeniopygia_guttata",
          "Xenopus_tropicalis","Danio_rerio","Drosophila_melanogaster","Arabidopsis_thaliana","Oryza_sativa","Gallus_gallus",
          "Mycobacterium_tuberculosis"))
}
