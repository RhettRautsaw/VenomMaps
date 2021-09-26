bind.tip2<-function(tree, tip.label, sp1, sp2=NA, perc=0.5) {
  if(is.na(sp2)){
    pos<-which(tree$tip.label==sp1)
    len<-tree$edge.length[which(tree$edge[,2]==pos)]
    len<-len*perc
    tree2<-bind.tip(tree, tip.label, edge.length=len, position=len, where=pos)
  }
  if(!is.na(sp2)){
    pos<-fastMRCA(tree, sp1, sp2)
    len<-tree$edge.length[which(tree$edge[,2]==pos)]
    len<-len*perc
    len2<-len+(fastDist(tree,sp1,sp2)/2)
    tree2<-bind.tip(tree, tip.label, edge.length=len2, position=len, where=pos)
  }
  return(tree2)
}