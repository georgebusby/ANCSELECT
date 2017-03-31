getfit=function(predmat,fitdata,restrict=1){
  
  temp=matrix(predmat[-restrict,],ncol=dim(predmat)[2])
  for(i in 1:nrow(temp)) temp[i,]=temp[i,]-predmat[restrict,]
  
  fitdata2=fitdata-predmat[restrict,]
  
  v=nnls(t(temp),fitdata2)
  x=v$x
  newx=1:nrow(predmat)
  newx[!((1:nrow(predmat))==restrict)]=x
  newx[restrict]=1-sum(x)
  v$x=newx
  names(v$x)=rownames(predmat)
  
  return(v)
}

#####the following does an initial linear model fit quicker than before...
getoverallfit=function(predmat,fitdata){
  restrict=1
  rep=1
  i=1
  while(rep==1){
    q=getfit(predmat,fitdata,restrict=i)
    
    if(q$x[i]>0) rep=0
    i=i+1
  }
  
  return(q)
}

#####need to require sums to equal 1
getfitfull=function(newpredmat,target,restrict){
  
  temp=newpredmat[,-restrict]
  for(i in 1:(ncol(temp)/2)) temp[,i]=temp[,i]-newpredmat[,restrict[1]]
  for(i in (ncol(temp)/2+1):(ncol(temp))) temp[,i]=temp[,i]-newpredmat[,restrict[2]]
  target2=target-newpredmat[,restrict[1]]-newpredmat[,restrict[2]]
  
  v=nnls(temp,target2)
  x=v$x
  newx=1:ncol(newpredmat)
  newx[!((1:ncol(newpredmat))%in%restrict)]=x
  newx[restrict[1]]=1-sum(x[1:(length(x)/2)])
  newx[restrict[2]]=1-sum(x[(length(x)/2+1):(length(x))])
  v$x=newx
  
  return(v)
}

####two populations, need 
getoverallfitfull=function(newpredmat,target){
  ourguess=nnls(newpredmat,target)$x
  ####use as a ranking
  
  guess=matrix(ourguess,nrow=2,byrow=T)
  orderi=order(guess[1,],decreasing=T)
  orderj=order(guess[2,],decreasing=T)+length(orderi)
  
  maxi=ncol(newpredmat)/2
  minj=ncol(newpredmat)/2+1
  maxj=ncol(newpredmat)
  
  finish=0;
  for(i in orderi){
    for(j in orderj){
      
      q=getfitfull(newpredmat,target,restrict=c(i,j))
      if(q$x[i]>0 & q$x[j]>0) finish=1;
      if(finish==1) break;
    }
    if(finish==1) break;
  }
  return(q)
}

getafit=function(alpha,mu,weights.mat, predmat, intercepts, fitdata){
  
  test=weights.mat%*%t(predmat)
  
  temppredmat=predmat
  vals=colSums(weights.mat)
  ######deal with lack of self copying
  for(i in 1:nrow(predmat)) temppredmat[i,]=temppredmat[i,]*vals
  for(i in 1:nrow(predmat)) temppredmat[i,]=temppredmat[i,]/sum(temppredmat[i,])
  
  newpredmat=cbind(mu*sqrt(alpha)*sqrt(1-alpha)*test,-mu*sqrt(alpha)*sqrt(1-alpha)*test)
  newpredmat=rbind(newpredmat,cbind((1-mu)*alpha*t(predmat), (1-mu)*(1-alpha)*t(predmat)))
  
  target=c(mu*intercepts,(1-mu)*fitdata)
  
  r=getoverallfitfull(newpredmat,target)
  names(r$x)=colnames(newpredmat)
  
  return(r)
}
