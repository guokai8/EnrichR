#' Disease Ontology Enrichment analysis function
#' @param df DGE files (DESeq2 result files) or vector contains gene names
#' @param keytype Keytype
#' @param filename output filename
#' @export
#' @author Kai Guo
DE<-function(df, keytype = "SYMBOL",padj.method = "BH", filename = NULL) {
  require(DO.db)
  data(hsdo)
  annot<-hsdo
  annotf <- toTable(DOTERM)
  annotf <- annotf[, c("do_id", "Term")]
  annotf <- unique(annotf)
  rownames(annotf)<-annotf[,1]
  annot.info =annotf
  if (keytype != "ENTREZID") {
    id = idconvert(species = "human", keys = annot$gene,
                   fkeytype = "ENTREZID", tkeytype = keytype)
    annot2 = annot[annot$gene %in% names(id), ]
    annot2$gene = id[annot2$gene %in% names(id)]
    annot=annot2
  }
  res = enrich(df, annot = annot, annot.info = annotf,
               filename = filename, padj.method = padj.method)
  res<-res[res$Pvalue<0.05,]
  return(res)
}
