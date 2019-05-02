#' plot combine network of Terms
#' @param df DGE files (DESeq2 result files) or vector contains gene names
#' @param gores GO Enrichment analsyis result
#' @param kores KEGG Enrichment analsyis result
#' @param rores Reactome Pathway Enrichment analsyis result
#' @param pfres PFAM Enrichment analysis result result
#' @param itres InterPro Enrichment analysis result
#' @param pvalue.cutoff the cut-off P value for selecting significant Terms
#' @param padj.cutoff the cut-off P adjust value for selecting significant Terms
#' @param weightcut the weight cut value for remove edges
#' @param useTerm use the Term description or not(defalutTRUE)
#' @param writeCyt export file for Cyt software
#' @param vertex.label.color color of label(defaultblack)
#' @param vertex.label.cex size of label(default0.5)
#' @param vertex.node.shape vector of shape and names of the vector should be the terms (default 20)
#' @param layout layout format (defultlayout.fruchterman.reingold)
#' @param visNet use VisNetwork method to display network(defaultFALSE)
#' @param top number of Terms you want to display(defaultthe total number of all Significant number)
#' @export
#' @author Kai Guo
mnetmap<-function(df,gores=NULL,kores=NULL,rores=NULL,pfres=NULL,
                   itres=NULL,top=NULL,top.display=NULL,pvalue.cutoff=0.05,padj.cutoff=NULL,
                   visNet=F,weightcut=0.2,useTerm=TRUE,layout=NULL,vertex.label.cex=0.5,gnet=FALSE,vertex.node.shape=NULL,...){
   if(!is.null(top)){
     rhs<-Reduce(function(x, y) rbind(x, y), list(gores[1:top,], kores[1:top,], rores[1:top,],pfres[1:top,],itres[1:top,]))
   }
   rhs<-Reduce(function(x, y) rbind(x, y), list(gores, kores, rores,pfres,itres))
   rhs<-rhs[order(rhs$Pvalue),]
   rhs<-rhs[rhs$Pvalue<pvalue.cutoff,]
   if(is.null(top.display)){
     top.display=nrow(rhs)
   }
   if(gnet==TRUE){
     gnet(df,rhs,top=top.display,pvalue.cutoff = pvalue.cutoff,padj.cutoff = padj.cutoff,visNet=visNet,weightcut = weightcut,useTerm=TRUE,layout=NULL,vertex.node.shape=vertex.node.shape,...)

   }else{
     netmap(df,rhs,top=top.display,pvalue.cutoff = pvalue.cutoff,padj.cutoff = padj.cutoff,visNet=visNet,weightcut = weightcut,useTerm=TRUE,layout=NULL,...)
   }
}
#' Plot compare heatmap of Enrichment result among DEG groups
#' @param rhslist of enrchment analysis result among DEG groups
#' @param top the number of Terms you want to display
#' @param colnames the compare DEG group names
#' @param xsize cex of group name
#' @param ysize cex of Terms name
#' @export
#' @author Kai Guo
lheatmap<-function (rhs, top = 50, colnames = NULL, xsize = 6, ysize = 6,padj=NULL,horizontal=FALSE,pval=0.05,returnData=FALSE,...)
{
  options(stringsAsFactors = F)
  suppressMessages(library(dplyr))
  suppressMessages(library(reshape2))
  suppressMessages(library(ggplot2))
  if(!is.null(padj)){
    rhs <- lapply(rhs, function(x) x[, c("Term", "Padj")])
    rhs<-lapply(rhs,function(x)x[x$Padj<padj,])
    sel<-unique(unlist(lapply(rhs,function(x)x%>%arrange(Padj)%>%
                                dplyr::select(Term)%>%.[[1]]%>%.[1:top])))
  }else{
    rhs <- lapply(rhs, function(x) x[, c("Term", "Pvalue")])
    rhs<-lapply(rhs,function(x)x[x$Pval<pval,])
    sel<-unique(unlist(lapply(rhs,function(x)x%>%arrange(Pvalue)%>%
                                dplyr::select(Term)%>%.[[1]]%>%.[1:top])))
  }
  res <- Reduce(function(x, y) full_join(x, y, by = "Term"),rhs)
  if (!is.null(colnames)) {
    colnames(res)[2:ncol(res)] <- colnames
  }
  else {
    colnames(res)[2:ncol(res)] <- paste("Group", 1:(ncol(res) -
                                                      1), sep = "_")
  }
  res[is.na(res)] <- 1
  res<-res%>%dplyr::filter(Term%in%sel)
  res<-as.data.frame(res)
  rownames(res)<-res$Term
  cor_mat<-cor(t(res[,2:ncol(res)]))
  dd <- as.dist((1-cor_mat)/2);
  hc <- hclust(dd);
  melted <- melt(res[hc$order,])
  melted$Term<-factor(melted$Term,levels=res$Term[hc$order])
  maxp = max(-log10(melted[, 3])) + 0.5
  if(!is.null(padj)){
    colnames(melted)[3] <- "Padj"
    p<-ggplot(melted, aes(x = variable, y = Term, fill = -log10(Padj))) +coord_equal(ratio = 0.8)+
      geom_tile(color = "white") + scale_fill_gradient2(low = "white",
                                                        high = "red", midpoint = 0, limit = c(0, maxp)) +theme_minimal()+theme(axis.text.y = element_text(size = ysize), axis.text.x = element_text(angle = 70,
                                                                                                                                                                                                    vjust = 1, size = xsize, hjust = 1,face = "bold")) +
      theme(axis.title.x = element_blank(),axis.title.y = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
  }else{
    colnames(melted)[3] <- "Pvalue"
    p<-ggplot(melted, aes(x = variable, y = Term, fill = -log10(Pvalue))) +coord_equal(ratio = 0.8)+
      geom_tile(color = "white") + scale_fill_gradient2(low = "white",
                                                        high = "red", midpoint = 0, limit = c(0, maxp)) +theme_minimal()+theme(axis.text.y = element_text(size = ysize), axis.text.x = element_text(angle = 70,
                                                                                                                                                                                                    vjust = 1, size = xsize, hjust = 1,face = "bold")) +
      theme(axis.title.x = element_blank(),axis.title.y = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
  }
  if(horizontal==TRUE){
    p<-p+coord_flip()
  }
  if(returnData==TRUE){
    return(res)
  }else{
    print(p)
  }
}
