#' Enrichment analysis for any type of annotation data
#' @importFrom dplyr filter_
#' @importFrom magrittr %>%
#' @param df DGE files (DESeq2 result files) or vector contains gene names
#' @param annot annotation file for all genes
#' @param annot.info Term of all annotation
#' @param filename output filename
#' @param gene.cutoff DGE singificant cutoff value
#' @param padj.method p value adjust method
#' @param keepRich keep terms with high rich factor
#' @export
#' @author Kai Guo
enrich<-function(df,annot,annot.info=NULL,minSize=1,maxSize=500,keepRich=TRUE,filename=NULL,gene.cutoff=0.01,padj.method="BH"){
  suppressMessages(require(tidyr))
  ao2gene<-sf(annot)
  ao2gene_num<-name_table(ao2gene)
  gene2ao<-sf(annot[,c(2,1)])
  if(is.data.frame(df)){
    IGE<-rownames(df)[df$padj<gene.cutoff]
  }else{
    IGE=as.vector(df)
  }
  fgene2ao=gene2ao[IGE]
  fao2gene=reverseList(fgene2ao)
  k=name_table(fao2gene)
  n=length(unique(unlist(fao2gene)))
  M=ao2gene_num[names(k)]
  N=length(unique(annot[,1]))
  rhs<-hyper_bench_vector(k,M,N,n)
  lhs<-p.adjust(rhs,method=padj.method)
  if(!is.null(annot.info)){
    all_ao<-annot.info
    rhs_an<-as.vector(all_ao[names(rhs),2])
    rhs_an<-sub('.*:','',rhs_an)
    rhs_gene<-unlist(lapply(fao2gene, function(x)paste(unique(x),sep="",collapse = ",")))
    resultFis<-data.frame("Annot"=names(rhs),"Term"=rhs_an,"Annotated"=M[names(rhs)],
                          "Significant"=k[names(rhs)],"Pvalue"=as.vector(rhs),"Padj"=lhs,
                          "GeneID"=rhs_gene[names(rhs)])
  }else{
    rhs_gene<-unlist(lapply(fao2gene, function(x)paste(unique(x),sep="",collapse = ",")))
    resultFis<-data.frame("Annot"=names(rhs),"Term"=names(rhs),"Annotated"=M[names(rhs)],
                          "Significant"=k[names(rhs)],"Pvalue"=as.vector(rhs),"Padj"=lhs,
                          "GeneID"=rhs_gene[names(rhs)])
  }
  resultFis<-resultFis[order(resultFis$Pvalue),]
  colnames(resultFis)[2]="Term"
  resultFis<-resultFis%>%filter_(~Significant<=maxSize)
  if(keepRich==FALSE){
    resultFis<-resultFis%>%filter_(~Significant>=minSize)
  }else{
    resultFis<-resultFis%>%filter_(~Significant>=minSize|(~Annotated/~Significant)==1)
  }
  if(!is.null(filename)){
    write.table(resultFis,file=paste(filename,".txt",sep=""),sep="\t",quote=F,row.names=F)
  }
  return(resultFis)
}
#' Display enrichment result By using barchart
#' @param  resultFis Ennrichment analysis result data.frame
#' @param  top Number of Terms you want to display
#' @param  pvalue.cutoff the cut-off value for selecting Term
#' @param  padj.cutoff the padj cut-off value for selecting Term
#' @param low color used for small value
#' @param high color used for large value
#' @param fontsize.x fontsize for x axis
#' @param fontsize.y fontsize for y axis
#' @param horiz show horiz or not (default: FALSE)
#' @param filename output filename
#' @param width width for output file
#' @param height height for output file
#' @export
#' @author Kai Guo
enrichbar<-function(resultFis,top=50,pvalue.cutoff=0.05,padj.cutoff=NULL,
                    low="lightpink",high="red",
                    order=FALSE,horiz=FALSE,fontsize.x=10,fontsize.y=10,
                    fontsize.text=3,angle=75,usePadj=TRUE,filename=NULL,
                    width=10,height=8){
  require(ggplot2)
  if(!is.null(padj.cutoff)){
    resultFis<-resultFis[resultFis$Padj<padj.cutoff,]
  }else{
    resultFis<-resultFis[resultFis$Pvalue<pvalue.cutoff,]
  }
  if(nrow(resultFis)>=top){
    resultFis<-resultFis[1:top,]
  }
  if(max(resultFis$Significant/(resultFis$Annotated+0.1))<=1){
    yheight=max(resultFis$Significant/resultFis$Annotated)+0.1
  }else{
    yheight=1
  }
  if(order==TRUE){
    resultFis$rich<-as.numeric(resultFis$Significant)/as.numeric(resultFis$Annotated)
    resultFis$Term<-factor(resultFis$Term,levels=resultFis$Term[order(resultFis$rich)])
  }
  if(usePadj==FALSE){
    p<-ggplot(resultFis,aes(x=Term,y=round(as.numeric(Significant/Annotated),2)))+geom_bar(stat="identity",aes(fill=-log10(as.numeric(Pvalue))))
    p<-p+scale_fill_gradient(low=low,high=high)+theme_light()
      if(horiz==TRUE){
      p<-p+theme(axis.text.y=element_text(face="bold",size=fontsize.y),axis.text.x=element_text(face="bold",color="black",size=fontsize.x,angle=angle))+labs(fill="-log10(Pvalue)")
      p<-p+coord_flip()
      p<-p+geom_text(aes(label=Significant),hjust=-0.3,size=fontsize.text)+xlab("Annotation")+ylab("Rich Factor")+ylim(0,yheight)
    }else{
      p<-p+theme(axis.text.y=element_text(face="bold",size=fontsize.y),axis.text.x=element_text(face="bold",color="black",size=fontsize.x,angle=angle,vjust=1,hjust=1))+labs(fill="-log10(Pvalue)")
      p<-p+geom_text(aes(label=Significant),vjust=-0.3,size=fontsize.text)+xlab("Annotation")+ylab("Rich Factor")+ylim(0,yheight)
    }
    print(p)
  }else{
    p<-ggplot(resultFis,aes(x=Term,y=round(as.numeric(Significant/Annotated),2)))+geom_bar(stat="identity",aes(fill=-log10(as.numeric(Padj))))
    p<-p+scale_fill_gradient2(low=low,high=high)+theme_light()
      if(horiz==TRUE){
      p<-p+theme(axis.text.y=element_text(face="bold",size=fontsize.y),axis.text.x=element_text(face="bold",color="black",size=fontsize.x,angle=angle))+labs(fill="-log10(Pvalue)")
      p<-p+coord_flip()
      p<-p+geom_text(aes(label=Significant),hjust=-0.3,size=fontsize.text)+xlab("Annotation")+ylab("Rich Factor")+ylim(0,yheight)
    }else{
      p<-p+theme(axis.text.y=element_text(face="bold",size=fontsize.y),axis.text.x=element_text(face="bold",color="black",size=fontsize.x,angle=angle,vjust=1,hjust=1))+labs(fill="-log10(Pvalue)")
      p<-p+geom_text(aes(label=Significant),vjust=-0.3,size=fontsize.text)+xlab("Annotation")+ylab("Rich Factor")+ylim(0,yheight)
    }
    print(p)
  }
  if(!is.null(filename)){
    ggsave(p,file=paste(filename,"barplot.pdf",sep="_"),width=width,height=height)
  }
}
#' Display enrichment result By using dotchart
#' @param resultFis: Ennrichment analysis result data.frame
#' @param top: Number of Terms you want to display
#' @param pvalue.cutoff: the cut-off value for selecting Term
#' @param padj.cutoff: the padj cut-off value for selecting Term
#' @param low color used for small value
#' @param high color used for large value
#' @param alpha alpha-transparency scales
#' @param fontsize.x fontsize for x axis
#' @param fontsize.y fontsize for y axis
#' @param filename output filename
#' @param width width for output file
#' @param height height for output file
#' @export
#' @author Kai Guo
enrichdot<-function(resultFis,top=50,pvalue.cutoff=0.05,order=FALSE,
                    low="lightpink",high="red",alpha=0.7,
                    padj.cutoff=NULL,fontsize.x=10,fontsize.y=10,
                    usePadj=TRUE,filename=NULL,width=10,height=8){
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
    dd$rich<-dd$Significant/dd$Annotated
    if(order==TRUE){
      dd$Term<-factor(dd$Term,levels=dd$Term[order(dd$rich)])
    }
    if(usePadj==FALSE){
      p<-ggplot(dd,aes(x=rich,y=Term))+geom_point(aes(size=Significant,color=-log10(Pvalue)),alpha=alpha)+
        theme(axis.text.y=element_text(face="bold",size=fontsize.y),axis.text.x=element_text(face="bold",color="black",size=fontsize.x))+
        scale_colour_gradient(low=low,high=high)+theme_minimal()+ylab("Pathway name")+
        xlab("Rich factor")+labs(size="Gene number")
      print(p)
    }else{
      p<-ggplot(dd,aes(x=rich,y=Term))+geom_point(aes(size=Significant,color=-log10(Padj)),alpha=alpha)+
        theme(axis.text.y=element_text(face="bold",size=fontsize.y),axis.text.x=element_text(face="bold",color="black",size=fontsize.x))+
        scale_colour_gradient(low=low,high=high)+theme_minimal()+ylab("Pathway name")+
        xlab("Rich factor")+labs(size="Gene number")
      print(p)
    }
      if(!is.null(filename)){
        ggsave(p,file=paste(filename,"dotplot.pdf",sep="_"),width=width,height=height)
      }
    }else{
      cat("No Pathway enrichment results were found!\n")
    }

}
