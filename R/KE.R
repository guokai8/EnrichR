#' KEGG Pathway Enrichment analysis function
#' @param df: DGE files (DESeq2 result files) or vector contains gene names
#' @param KG_FILE: KEGG Pathway annotation data
#' @param filenam: output filename
#' @param cutoff: DGE singificant cutoff value
#' @export
KE<-function(df,KO_FILE,filename=NULL,gene.cutoff=0.01,minSize=2,maxSize=500,keepRich=TRUE,padj.method="BH"){
  ko2gene<-sf(KO_FILE)
  ko2gene_num<-name_table(ko2gene)
  gene2ko<-sf(KO_FILE[,c(2,1)])
  if(is.data.frame(df)){
    IGE<-rownames(df)[df$padj<gene.cutoff]
  }else{
    IGE=as.vector(df)
  }
  fgene2ko=gene2ko[IGE]
  fko2gene=reverseList(fgene2ko)
  k=name_table(fko2gene)
  n=length(unique(unlist(fko2gene)))
  M=ko2gene_num[names(k)]
  N=length(unique(KO_FILE[,1]))
  rhs<-hyper_bench_vector(k,M,N,n)
  lhs<-p.adjust(rhs,method=padj.method)
  all_ko<-.get_kg_dat()
 # rhs_an<-all_ko[rownames(all_ko)%in%names(rhs),]
  rhs_an<-all_ko[names(rhs),]
  rhs_gene<-unlist(lapply(fko2gene, function(x)paste(unique(x),sep="",collapse = ",")))
  resultFis<-data.frame("Annot"=names(rhs),"Term"=rhs_an,"Annotated"=M[names(rhs)],
                        "Significant"=k[names(rhs)],"Pvalue"=as.vector(rhs),"Padj"=lhs,
                        "GeneID"=rhs_gene[names(rhs)])

  resultFis<-resultFis[order(resultFis$Pvalue),]
  resultFis<-resultFis[resultFis$Pvalue<0.05,]
  resultFis<-resultFis%>%dplyr::filter(Significant<=maxSize)
  if(keepRich==FALSE){
    resultFis<-resultFis%>%dplyr::filter(Significant>=minSize)
  }else{
    resultFis<-resultFis%>%dplyr::filter(Significant>=minSize|(Annotated/Significant)==1)
  }
  if(!is.null(filename)){
    write.table(resultFis,file=paste(filename,".txt",sep=""),sep="\t",quote=F,row.names=F)
  }
  return(resultFis)
}
#' Display KEGG enrichment result
#' @param  resultFis: KEGG ennrichment analysis result data.frame
#' @param  top: Number of Terms you want to display
#' @param  filenam: output filename
#' @param  pvalue.cutoff: the cut-off value for selecting Term
#' @param  padj.cutoff: the padj cut-off value for selecting Term
#' @export
#' @author Kai Guo
KE.plot<-function(resultFis,pvalue.cutoff=0.05,top=50,order=FALSE,fontsize.x=10,fontsize.y=10,padj.cutoff=NULL,usePadj=TRUE,filename=NULL){
  library(ggplot2)
  if(!is.null(padj.cutoff)){
    resultFis<-resultFis[resultFis$Padj<padj.cutoff,]
  }else{
    resultFis<-resultFis[resultFis$Pvalue<pvalue.cutoff,]
  }
  if(nrow(resultFis)>=top){
    dd<-resultFis[1:top,]
  }else{
    dd<-resultFis
  }
  if(nrow(dd)>=1){
  dd[,3]<-dd[,4]/dd[,3]
  colnames(dd)[3]<-"rich";
  if(order==TRUE){
    dd$Term<-factor(dd$Term,levels=dd$Term[order(dd$rich)])
  }
  if(usePadj==FALSE){
    p<-ggplot(dd,aes(x=rich,y=Term))+geom_point(aes(size=Significant,color=-log10(Pvalue)))++theme_minimal()+
      theme(axis.text.y=element_text(face="bold",size=fontsize.y),axis.text.x=element_text(face="bold",color="black",size=fontsize.x))+
    scale_colour_gradient(low="lightpink",high="red")+ylab("Pathway name")+
    xlab("Rich factor")+labs(size="Gene number")
    print(p)
  }else{
    p<-ggplot(dd,aes(x=rich,y=Term))+geom_point(aes(size=Significant,color=-log10(Padj)))+theme_minimal()+
      theme(axis.text.y=element_text(face="bold",size=fontsize.y),axis.text.x=element_text(face="bold",color="black",size=fontsize.x))+
      scale_colour_gradient(low="lightpink",high="red")+ylab("Pathway name")+
      xlab("Rich factor")+labs(size="Gene number")
    print(p)
  }}else{
    cat("No Pathway enrichment results were found!\n")
  }
  if(!is.null(filename)){
  ggsave(p,file=paste(filename,"KEGG.pdf",sep="_"),width=10,height=9)
  }
}
