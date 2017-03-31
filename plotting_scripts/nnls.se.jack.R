se.jack=function(nnls,nnlsjack,chromovec){
  #nnls=my.newNNLSorig
  #nnlsjack=jacKnife
  #chromovec=chrom.nsnps
  poplist=matrix(NA,nrow=nrow(nnls),ncol=ncol(nnls))
  h.j=sum(chromovec)/chromovec
  jackValue=function(nnlsj.vec){sum((1.0-chromovec/sum(chromovec))*nnlsj.vec)}
  se.value=function(nnls.obs,nnls.jvec,jack.val){sum((1.0/(h.j-1.0))*
                                                       (h.j*nnls.obs-(h.j-1.0)*
                                                          nnls.jvec-length(chromovec)*nnls.obs+jack.val)^2)}
  i=1
  if (unique((rownames(nnls)==rownames(nnlsjack[[1]])))==TRUE){
    for (x in c(1)){
      tmp.nnls=t(sapply(nnlsjack,function(X)as.numeric(X[x,])))
      tmp.nnls <- data.frame(tmp.nnls)
      names(tmp.nnls) <- names(nnlsjack[[1]])
      tmp.nnls <- tmp.nnls[names(nnls)]
      tmp.jack=apply(tmp.nnls,2,jackValue)
      a=list()
      for (y in 1:ncol(nnls)){
        a[[y]]=se.value(nnls[x,y],tmp.nnls[,y],tmp.jack[y])
      }
      a=unlist(a)
      poplist[i,]=a
      i=i+1
      
    } } else {stop ("Rownames of the two inputs are DIFFERENT")}
  colnames(poplist)=colnames(nnls)
  rownames(poplist)=rownames(nnls)
  poplist=sqrt(poplist/22)
  return(poplist)}