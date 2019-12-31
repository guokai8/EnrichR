#' KEGG Pathway Enrichment analysis function
#' @importFrom dplyr filter_
#' @param df DGE files (DESeq2 result files) or vector contains gene names
#' @param KO_FILE KEGG annotation data
#' @param minSize minimal number of genes included in significant terms
#' @param padj.method p value adjust method (default: BH)
#' @param cutoff cutoff value for filtering significant terms (default: 0.05)
#' @param maxSize maximum number of genes included in significant terms
#' @param keepRich keep terms with rich factor value equal 1 or not (default: TRUE)
#' @param filename output filename
#' @param gene.cutoff the cut-off value for select DEGs (default: 0.01)
#' @param bulitin use KEGG bulit in KEGG annotation or not(set FALSE if you want use newest KEGG data)
#' @export
KE<-function(df,KO_FILE,filename=NULL,gene.cutoff=0.01,minSize=2,maxSize=500,keepRich=TRUE,padj.method="BH",cutoff=0.05,builtin=TRUE){
  suppressMessages(require(tidyr))
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
  all_ko<-.get_kg_dat(builtin=builtin)
  rhs_an<-all_ko[names(rhs),]
  rhs_gene<-unlist(lapply(fko2gene, function(x)paste(unique(x),sep="",collapse = ",")))
  resultFis<-data.frame("Annot"=names(rhs),"Term"=rhs_an,"Annotated"=M[names(rhs)],
                        "Significant"=k[names(rhs)],"Pvalue"=as.vector(rhs),"Padj"=lhs,
                        "GeneID"=rhs_gene[names(rhs)])

  resultFis<-resultFis[order(resultFis$Pvalue),]
  resultFis<-resultFis[resultFis$Pvalue<cutoff,]
  resultFis<-filter_(resultFis, ~Significant<=maxSize)
  resultFis<-resultFis[order(resultFis$Pvalue),]
  resultFis<-resultFis[resultFis$Pvalue<cutoff,]
  resultFis<-filter_(resultFis, ~Significant<=maxSize)
  if(keepRich==FALSE){
    resultFis<-filter_(resultFis, ~Significant>=minSize)
  }else{
    resultFis<-filter_(resultFis, ~Significant>=minSize|(~Annotated/~Significant)==1)
  }
  if(!is.null(filename)){
    write.table(resultFis,file=paste(filename,".txt",sep=""),sep="\t",quote=F,row.names=F)
  }
  return(resultFis)
}
#' Display KEGG enrichment result
#' @param resultFis: KEGG ennrichment analysis result data.frame
#' @param top: Number of Terms you want to display
#' @param filename: output filename
#' @param pvalue.cutoff: the cut-off value for selecting Term
#' @param padj.cutoff: the padj cut-off value for selecting Term
#' @param usePadj use adjust pvalue or not
#' @param order order bar or not
#' @param low color used for small value
#' @param high color used for large value
#' @param alpha alpha-transparency scales
#' @param horiz use horiz or not
#' @param fontsize.x fontsize for x axis
#' @param fontsize.y fontsize for y axis
#' @param filename output filename
#' @param width width for output file
#' @param height height for output file
#' @export
#' @author Kai Guo
KE.plot<-function(resultFis,pvalue.cutoff=0.05,top=50,order=FALSE,
                  low="lightpink",high="red",alpha=0.7,
                font.x="bold",font.y="bold",fontsize.x=10,fontsize.y=10,
                padj.cutoff=NULL,usePadj=TRUE,filename=NULL,width=10,height=8){
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
    p<-ggplot(dd,aes(x=rich,y=Term))+geom_point(aes(size=Significant,color=-log10(Pvalue)),alpha=alpha)+theme_minimal()+
      theme(axis.text.y=element_text(face=font.y,size=fontsize.y),axis.text.x=element_text(face=font.x,color="black",size=fontsize.x))+
    scale_colour_gradient(low=low,high=high)+ylab("Pathway name")+
    xlab("Rich factor")+labs(size="Gene number")
    print(p)
  }else{
    p<-ggplot(dd,aes(x=rich,y=Term))+geom_point(aes(size=Significant,color=-log10(Padj)),alpha=alpha)+theme_minimal()+
      theme(axis.text.y=element_text(face=font.y,size=fontsize.y),axis.text.x=element_text(face=font.x,color="black",size=fontsize.x))+
      scale_colour_gradient(low=low,high=high)+ylab("Pathway name")+
      xlab("Rich factor")+labs(size="Gene number")+theme(color=guide_colourbar(order = 1),size=guide_legend(order = 2))
    print(p)
  }}else{
    cat("No Pathway enrichment results were found!\n")
  }
  if(!is.null(filename)){
  ggsave(p,file=paste(filename,"KEGG.pdf",sep="_"),width=width,height=height)
  }
}
