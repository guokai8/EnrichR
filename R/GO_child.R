#' Get all children terms of node
#' @param  node : input node of GO
#' @param  ontology: ontology term of BP
#' @author Kai Guo
GO_child <- function(node = "GO:0008150", ontology = "BP") {
  #MF = "GO:0003674", node of MF
  #BP = "GO:0008150", node of BP
  #CC = "GO:0005575", node of CC
  if (ontology == "BP") res <- c(node,GOBPOFFSPRING[[node]])
  if (ontology == "CC") res <- c(node,GOCCOFFSPRING[[node]])
  if (ontology == "MF") res <- c(node,GOMFOFFSPRING[[node]])
  return(res[!is.na(res)])
}
