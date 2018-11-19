#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats predict aggregate terms
#' @importFrom party cforest cforest_control
.computePDelta=function(formula,data, intervention, forestControl){
  classes=unique(data[[intervention]])
  nclasses=length(classes)
  nvar=length(attr(terms(formula),"term.labels"))
  B=forestControl$nBoots

  n=dim(data)[1]
  P=matrix(0, nrow=n,ncol=nclasses)
  Delta=matrix(0, nrow=n,ncol=nclasses^2)

  dataF=data
  Boots=array(0, dim=c(n,nclasses,B))
  cat("Bootstrap simulations running...\n")
  pb=txtProgressBar(style=3)
  for(b in 1:B){
    setTxtProgressBar(pb, b/B)
    for (k in 1:nclasses){
      ids=which(data[[intervention]]==classes[k])
      bootid=sample(ids, replace = T)
      mod=cforest(formula, data=data[bootid,],
                            controls=party::cforest_control(ntree=forestControl$nTrees,minsplit=forestControl$minsplit,
                                                            mincriterion=forestControl$mincriterion,
                                                           replace=TRUE,mtry=forestControl$mtry,  testtype = "Univariate"))

      out=predict(mod, newdata = data)
      out1=predict(mod, OOB=T)
      Boots[,k,b]=out
      Boots[bootid,k,b]=out1
    }
  }
  close(pb)

  for (i in 1:n){
    Boot=Boots[i,,]
    Maxs=apply(Boot,MARGIN = 2, FUN = which.max)
    Pr=table(Maxs)/B
    P[i,as.numeric(names(Pr))]=Pr
    Boot1=data.frame(Bs=t(Boot[, order(Maxs)]),Ms=sort(Maxs))
    M0=aggregate(Boot1, by=list(Boot1$Ms), FUN=mean)
    Gr=M0$Ms
    M0$Ms=c()
    M0[[1]]=c()
    M01=matrix(0, nrow=nclasses, ncol=nclasses)
    M01[Gr,]=as.matrix(M0)
    M1=t(M01)
    Mp=matrix(rep(diag(M1), times=nclasses), nrow=nclasses, byrow =T)
    Delta[i, ]=as.vector(t(Mp-M1))

  }
  #Loss from setting number i to j is Delta[j,i]

  dataF$P=I(P)
  dataF$Delta=I(Delta)
  return(list(data=dataF, nclasses=nclasses,classes=classes))
}

#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats predict aggregate terms rnorm
#' @importFrom partykit ctree ctree_control
.computePDeltaNormal=function(formula,data, intervention, forestControl){

  classes=unique(data[[intervention]])
  nclasses=length(classes)
  nvar=length(attr(terms(formula),"term.labels"))
  M=forestControl$nTrees
  B=forestControl$nBoots

  n=dim(data)[1]
  P=matrix(0, nrow=n,ncol=nclasses)
  Delta=matrix(0, nrow=n,ncol=nclasses^2)

  dataF=data
  Boots=array(0, dim=c(n,nclasses,M))
  BootsOOB=array(NA, dim=c(n,nclasses,M))
  cat("Computing variances by infinitesimal jacknife with bias correction...\n")
  pb=txtProgressBar(style=3)
  Scal=array(0, dim=c(n,nclasses,M))
  ns=numeric(nclasses)
  for(b in 1:M){
    setTxtProgressBar(pb, b/M)
    for (k in 1:nclasses){
      ids=which(data[[intervention]]==classes[k])
      ns[k]=length(ids)
      bootid=sample(ids,ns[k], replace=TRUE)

     mod=ctree(formula, data=data[bootid,],
                         control=partykit::ctree_control(minsplit=forestControl$minsplit,
                                                        mincriterion=forestControl$mincriterion,
                                                        mtry =forestControl$mtry, testtype = "Univariate" ))
       Nb0=table(bootid)
      Nb1=numeric(n)
      Nb1[as.numeric(names(Nb0))]=Nb0
      Nb1[ids]=Nb1[ids]-1
      Nb=Nb1
      out=predict(mod, newdata = data)
      out2=predict(mod, newdata=data[-bootid,])
      Boots[,k,b]=out

      BootsOOB[-bootid,k,b]=out2

      Scal[,k,b]=Nb
    }
  }
  BootsN=array(0, dim=c(n,nclasses,B))
  close(pb)
  cat("Generating bootstrap Samples...\n")
  pb1=txtProgressBar(style=3)
  BootsMO= apply(BootsOOB, MARGIN = c(1,2), FUN = mean, na.rm=T)
  for (i in 1:n){
    setTxtProgressBar(pb, i/n)
    Booti=Boots[rep(i,n),,]
    BootsM= apply(Booti, MARGIN = c(1,2), FUN = mean)

    BootsMexp=array(BootsM, dim=c(n,nclasses,M))
    Boots1=Booti-BootsMexp
    Covs=apply(Boots1*Scal, MARGIN=c(1,2), FUN=mean)
    T1=apply(Covs^2, MARGIN=c(2), FUN=sum)
    T2=apply(Boots[i,,], MARGIN=c(1), FUN=function(x) sum((x-mean(x))^2)/(M^2))
    Vi=T1-T2*ns+2*T2*sqrt(2*ns)
    if(any(Vi<0)) {
      ii=which(Vi<0)
      Vi[ii]=2*T2[ii]*sqrt(2*ns[ii])
      warning("Estimate may be unreliable, increase value of nTrees parameter!")
    }
    for(j in 1:nclasses){
      BootsN[i,j,]=rnorm(B,BootsMO[i,j], sqrt(Vi[j]))
    }


  }
  Boots=BootsN
  close(pb1)


  for (i in 1:n){
    Boot=Boots[i,,]
    Maxs=apply(Boot,MARGIN = 2, FUN = which.max)
    Pr=table(Maxs)/B
    P[i,as.numeric(names(Pr))]=Pr
    Boot1=data.frame(Bs=t(Boot[, order(Maxs)]),Ms=sort(Maxs))
    M0=aggregate(Boot1, by=list(Boot1$Ms), FUN=mean)
    Gr=M0$Ms
    M0$Ms=c()
    M0[[1]]=c()
    M01=matrix(0, nrow=nclasses, ncol=nclasses)
    M01[Gr,]=as.matrix(M0)
    M1=t(M01)
    Mp=matrix(rep(diag(M1), times=nclasses), nrow=nclasses, byrow =T)
    Delta[i, ]=as.vector(t(Mp-M1))

  }
  #Loss from setting number i to j is Delta[j,i]

  dataF$P=I(P)
  dataF$Delta=I(Delta)
  return(list(data=dataF, nclasses=nclasses,classes=classes))
}

#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats predict aggregate terms
#' @importFrom randomForest randomForest
.computePDeltaRF=function(formula,data, intervention, forestControl){
  classes=unique(data[[intervention]])
  nclasses=length(classes)
  B=forestControl$nBoots

  n=dim(data)[1]
  P=matrix(0, nrow=n,ncol=nclasses)
  Delta=matrix(0, nrow=n,ncol=nclasses^2)

  dataF=data
  Boots=array(0, dim=c(n,nclasses,B))
  cat("Bootstrap simulations running...\n")
  pb=txtProgressBar(style=3)
  for(b in 1:B){
    setTxtProgressBar(pb, b/B)
    for (k in 1:nclasses){
      ids=which(data[[intervention]]==classes[k])
      bootid=sample(ids, replace = T)
      mod=randomForest::randomForest(ntree=forestControl$nTrees, formula, data=data,subset=bootid,
                            nodesize=forestControl$minsplit/2)
      out=predict(mod, newdata = data)
      Boots[,k,b]=out
    }
  }
  close(pb)

  for (i in 1:n){
    Boot=Boots[i,,]
    Maxs=apply(Boot,MARGIN = 2, FUN = which.max)
    Pr=table(Maxs)/B
    P[i,as.numeric(names(Pr))]=Pr
    Boot1=data.frame(Bs=t(Boot[, order(Maxs)]),Ms=sort(Maxs))
    M0=aggregate(Boot1, by=list(Boot1$Ms), FUN=mean)
    Gr=M0$Ms
    M0$Ms=c()
    M0[[1]]=c()
    M01=matrix(0, nrow=nclasses, ncol=nclasses)
    M01[Gr,]=as.matrix(M0)
    M1=t(M01)
    Mp=matrix(rep(diag(M1), times=nclasses), nrow=nclasses, byrow =T)
    Delta[i, ]=as.vector(t(Mp-M1))

  }
  #Loss from setting number i to j is Delta[j,i]

  dataF$P=I(P)
  dataF$Delta=I(Delta)
  return(list(data=dataF, nclasses=nclasses,classes=classes))
}



#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats predict aggregate terms
#' @importFrom BayesTree bart
.computePDeltaBart=function(formula,data, intervention, forestControl){
  classes=unique(data[[intervention]])
  nclasses=length(classes)
  B=forestControl$nBoots

  n=dim(data)[1]
  P=matrix(0, nrow=n,ncol=nclasses)
  Delta=matrix(0, nrow=n,ncol=nclasses^2)

  dataF=data
  Boots=array(0, dim=c(n,nclasses,B))
  vars=rownames(attr(terms(formula),"factors"))
  nvars=length(vars)
  x.vars=vars[2:nvars]
  y.var=vars[1]
  if(nvars==2){
    x.train=data[,x.vars]

  } else{
    x.train=as.data.frame(data[,x.vars])
  }

  y.train=data[,y.var]


  cat("BART simulations running...\n")
  pb=txtProgressBar(style=3)

  for (k in 1:nclasses){
    setTxtProgressBar(pb, k/nclasses)
      ids=which(data[[intervention]]==classes[k])

      if(nvars==2)
        mod=BayesTree::bart(x.train[ids], y.train[ids], x.train, ntree=forestControl$nTrees, ndpost=B, verbose = F, sigest=1)
      else
        mod=BayesTree::bart(x.train[ids,], y.train[ids], x.train, ntree=forestControl$nTrees, ndpost=B, verbose = F, sigest=1)


      Boots[,k,]=t(mod$yhat.test)
  }

  close(pb)

  for (i in 1:n){
    Boot=Boots[i,,]
    Maxs=apply(Boot,MARGIN = 2, FUN = which.max)
    Pr=table(Maxs)/B
    P[i,as.numeric(names(Pr))]=Pr
    Boot1=data.frame(Bs=t(Boot[, order(Maxs)]),Ms=sort(Maxs))
    M0=aggregate(Boot1, by=list(Boot1$Ms), FUN=mean)
    Gr=M0$Ms
    M0$Ms=c()
    M0[[1]]=c()
    M01=matrix(0, nrow=nclasses, ncol=nclasses)
    M01[Gr,]=as.matrix(M0)
    M1=t(M01)
    Mp=matrix(rep(diag(M1), times=nclasses), nrow=nclasses, byrow =T)
    Delta[i, ]=as.vector(t(Mp-M1))

  }

  dataF$P=I(P)
  dataF$Delta=I(Delta)
  return(list(data=dataF, nclasses=nclasses,classes=classes))
}
