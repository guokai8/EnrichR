#' reverse List
#' @param lhs: list with names
#' @export
#' @author Kai Guo
reverseList<-function(lhs){
  lhs_n<-rep(names(lhs),times=lapply(lhs,function(x)length(x)))
  res<-sf(as.data.frame(cbind(lhs_n,unlist(lhs))))
  return(res)
}
