
minit=function(y, offset, parms,wt){
  sfun=function(yval,dev,wt,ylevel,digits){
    paste("NAN")
  }
  text=function (yval, dev, wt, ylevel, digits, n, use.n)
  {
  }
  list(y=y, parms=parms, numresp=1, numy=1, summary=sfun, text=text)
}

#' @importFrom stats chisq.test sd
#' @importFrom utils combn
msplit=function(y,wt,x,parms, continuous){
  n=length(x)


  goodnessC=0
  ux=unique(x)
  nu=length(ux)

  l=getLabel(parms$P, parms$Delta, parms$m, parms$classes, y[1:n], parms$cost, parms$confidence)
  if(continuous){
    goodness=numeric(n-1)
    for(i in 1:(n-1)){
      # browser()
      ind1=y[1:i]
      ind2=y[(i+1):n]
      l1=getLabel(parms$P, parms$Delta, parms$m, parms$classes, ind1, parms$cost, parms$confidence)
      l2=getLabel(parms$P, parms$Delta, parms$m, parms$classes, ind2, parms$cost, parms$confidence)
      goodness[i]=l$deviance-l1$deviance-l2$deviance

     if(parms$prune){
        if(i>1) P1=colSums(parms$P[ind1,])*0.28/(mean(apply(parms$P[ind1,], MARGIN = 2,FUN=sd))+0.05) else P1=parms$P[ind1,]
        if (i+1<n) P2=colSums(parms$P[ind2,])*0.28/(mean(apply(parms$P[ind2,], MARGIN = 2,FUN=sd))+0.05) else P2=parms$P[ind2,]
        PP=matrix(c(P1,P2), ncol = 2)
        PP[PP<=0.5]=1
        if(all(PP>6)) pv=suppressWarnings(chisq.test(round(PP))$p.value) else pv=suppressWarnings(chisq.test(round(PP))$p.value)
        goodness[i]=goodness[i]*ifelse(pv<0.05,1,0)

      }


    }
    direction=rep(-1, n-1)

  } else{
    goodness=numeric(nu-1)
    ChoosenOrder=c()
    opti=NA
    for (i in (1: round(nu/2))){
      variants=combn(1:nu, i)
      for(j in 1:ncol(variants)){
        ind=numeric(nu)
        ind[variants[,j]]=1
        indic=ind[as.numeric(x)]
        ind1=which(indic==0)
        ind2=which(indic==1)
        yind1=y[ind1]
        yind2=y[ind2]
        if (length(ind1)>0 & length(ind2)>0){
          l1=getLabel(parms$P, parms$Delta, parms$m, parms$classes, yind1, parms$cost, parms$confidence)
          l2=getLabel(parms$P, parms$Delta, parms$m, parms$classes, yind2, parms$cost, parms$confidence)
          g=l$deviance-l1$deviance-l2$deviance
          if(parms$prune){
            if(length(yind1)>1) P1=colSums(parms$P[yind1,])*0.28/(mean(apply(parms$P[yind1,], MARGIN = 2,FUN=sd))+0.05) else P1=parms$P[yind1,]
            if (length(yind2)>1) P2=colSums(parms$P[yind2,])*0.28/(mean(apply(parms$P[yind2,], MARGIN = 2,FUN=sd))+0.05) else P2=parms$P[yind2,]
            PP=matrix(c(P1,P2), ncol = 2)
            PP[PP<=0.5]=1
            if(all(PP>6)) pv=suppressWarnings(chisq.test(round(PP))$p.value) else pv=suppressWarnings(chisq.test(round(PP))$p.value)
            g=g*ifelse(pv<0.05,1,0)

          }
          if(g>goodnessC){
            goodnessC=g
            ChoosenOrder=variants[,j]
            opti=i
          }
        }

      }
      goodness[opti]=goodnessC
      direction=c(ux[ChoosenOrder], ux[setdiff(1:nu, ChoosenOrder)])
    }

  }

  return(list(goodness=goodness, direction=direction))
}



getLabelP=function(P, Delta, m, classes, ind, confidence){
  if(length(ind)>1) Ps=colSums(P[ind,]) else Ps=P[ind,]
  Ps=Ps/sum(Ps)

  Pso=sort(Ps, decreasing=T)
  Psi=order(Ps, decreasing = T)
  Psoc=cumsum(Pso)
  Psoi=which(Psoc>=confidence)
  ind2=Psoi[1]
  if(is.na(ind2)||is.nan(ind2)) browser()
  S=Psi[1:ind2]
  PsA=Ps
  PsA[-S]=0
  PsA=PsA/sum(PsA)
  Loss=sum(Ps*(1-PsA))*length(ind)


  V=Ps
  St0=paste(paste("p(",classes,")", "=",sep=""),round(V,digits=2), collapse="\n")
  St1=paste("[", paste(classes[S], collapse = ","), "]\n")
  St=paste(St1,St0, sep = "")

  return(list(label=St, l=classes[S],deviance=Loss, Ps=PsA))

}
getLabelDelta=function(P, Delta, m, classes, ind, confidence){
  # browser()
  Prep=matrix(rep(P[ind,],m), ncol=m^2)
  if(length(ind)>1) DeltaP=colSums(Delta[ind,]*Prep) else DeltaP=Delta[ind,]*Prep
  DeltaP1=colSums(matrix(DeltaP,nrow=m))


  if(length(ind)>1) Ps=colSums(P[ind,]) else Ps=P[ind,]
  Ps=Ps/sum(Ps)

  Pso=sort(Ps, decreasing=T)
  Psi=order(Ps, decreasing = T)
  Psoc=cumsum(Pso)
  Psoi=which(Psoc>confidence)
  ind2=Psoi[1]
  if(is.na(ind2)||is.nan(ind2)) browser()
  S=Psi[1:ind2]
  PsA=Ps
  PsA[-S]=0
  PsA=PsA/sum(PsA)
  Loss=sum(PsA*DeltaP1)*length(ind)


  V=Ps
  St0=paste(paste("p(",classes,")", "=",sep=""),round(V,digits=2), collapse="\n")
  St1=paste("[", paste(classes[S], collapse = ","), "]\n")
  St=paste(St1,St0, sep = "")

  return(list(label=St, l=classes[S],deviance=Loss, Ps=PsA))

}

getLabel=function(P, Delta, m, classes, ind, cost, confidence){
  if (cost) return(getLabelDelta(P, Delta, m, classes, ind, confidence))
  else return(getLabelP(P, Delta, m, classes, ind, confidence))
}



meval=function(y,wt,parms){
  res=getLabel(parms$P, parms$Delta,parms$m, parms$classes, y, parms$cost, parms$confidence)
  env1=sys.frame(-2)
  env1$LABELS[[env1$nLabel]]=res$label
  env1$pLabel[[env1$nLabel]]=res$l
  env1$Probs[[env1$nLabel]]=res$Ps
  env1$nLabel=env1$nLabel+1
  return(list(label=paste(env1$nLabel-1), deviance=res$deviance))
}


createpsicaTree=function(formula, data, parms, control, ...){
  LABELS=list()
  Probs=list()
  nLabel=1
  pLabel=list()

  res=rpart::rpart(formula,data,
            method=list(eval=meval,split=msplit, init=minit),
            parms = parms, control=control,...)
  class(res)<-c("psicaTree", class(res))
  res$LABELS=LABELS
  res$nLabel=nLabel
  res$pLabel=pLabel
  res$Probs=Probs

  return(res)
}

#' Create a tree that discovers groups having similar treatment
#' (intervention) effects.
#'
#' @param formula Formula that shows the dependent variable (effect)
#' and independent variables (separated by '+'). The treatment
#' variable should not be present among dependent variables
#' @param data Data frame containing dependent and independent variables
#' and the categorical treatment variable
#' @param method Choose "boot" for computing probabilities by bootstrapping random forests,
#' "normal" for computing probabilities by appoximating random forest variance with
#' infinitesimal jackknife with bias correction.
#' @param intervention The name of the treatment variable
#' @param forestControl parameters of forest growing, a list with parameters
#' \itemize{
#' \item minsplit:  minimum number of observation in the node to
#' be splitted when growing random forest, default 10
#'  \item mincriterion: "mincriterion" setting of the random forest,
#' see ctree in package partykit
#' \item nBoots: number of trees in random forest.
#' \item nTrees: amount of trees in each random forest
#' \item mtry: number of variables to be selected at each split. Choose either 'sqrt(amount_of_inputs)'
#' if amount of input variables is large or 'amount_of_inputs' if there are few input variables.
#' }
#' @param treeControl Parameters for decision tree growing,
#' see rpart.control()
# @param cost Parameter that splits the tree by only taking the probabilities that
# the a treatment is better than other treatments (cost=F) or also how much this
# treatment is better than the other treatments (cost=T)
#' @param confidence Parameter that defines the cut-off probability in the loss function
#' and also which treatments are included in the labels of the PSICA tree. More specifically,
#' labels in the terminal nodes show all treatments except of useless treatments, i.e. the
#' treatments that altogether have a probability to be the best which is smaller than 1-confidence.
#' @param prune should the final tree be pruned or is (possibly) overfitted tree desired?
#' @param ... further argumets passed to rpart object.
#' @return Object of a class psicaTree
#' @description
#' The PSICA method operates by first building regression trees for each treament group
#' and then obtaining the distributions of the effect size for given levels of independent
#' variables by either bootstrap or by means of the bias-corrected infinitesimal jackknife.
#' The obtained distributions are used for computing the probabilities that one treatment
#' is better (effect size is greater) than the other treatments for a given
#' set of input values. These probabilities are then summarised in the form of a decision tree
#' built with a special loss function. The terminal nodes of the resulting tree show
#' the probabilities that one  treatment is better than the other treatments
#' as well as a label containing the possible best treatments.
#' @examples
#' n=100
#' X1=runif(n)
#' X1=sort(X1)
#' f1<- function(x){
#'   2*tanh(4*x-2)+3
#' }
#' X2=runif(n)
#' X2=sort(X2)
#' f2<- function(x){
#'   2*tanh(2*x-1)+2.3 #2.8
#' }
#' plot(X1,f1(X1),ylim=c(0,5), type="l")
#' points(X2,f2(X2), type="l")
#' Y1=f1(X1)+rnorm(n, 0, 0.8)
#' Y2=f2(X2)+rnorm(n,0,0.8)

#' points(X1,Y1, col="blue")
#' points(X2,Y2, col="red")

#' data=data.frame(X=c(X1,X2), Y=c(Y1,Y2), interv=c(rep("treat",n), rep("control",n)))

#' pt=psica(Y~X, data=data, method="normal",intervention = "interv",
#'  forestControl=list(nBoots=200, mtry=1))
#' print(pt)
#' plot(pt)
#'
#' @references
#' \insertRef{psica}{psica}
#' @export
#' @importFrom stats as.formula terms
#' @importFrom Rdpack reprompt
psica=function(formula, data, intervention, method="normal",
               forestControl=list(minsplit=10, mincriterion=0.95, nBoots=500, nTrees=200, mtry=5),
               treeControl=rpart::rpart.control(minsplit=20, minbucket=10, cp=0.003),
                confidence=0.95, prune=TRUE,...){

  #browser()
  cost=F
  forestControl1=list(minsplit=10, mincriterion=0.95, nBoots=500, nTrees=30, mtry=5)
  fNames=names(forestControl1)
  for (nm in fNames){
    if (is.null(forestControl[[nm]])) forestControl[[nm]]=forestControl1[[nm]]
  }

  treeControl1=rpart::rpart.control(minsplit=20, minbucket=10, cp=0.003)
  tNames=names(treeControl1)
  for (nm in tNames){
    if (is.null(treeControl[[nm]])) treeControl[[nm]]=treeControl1[[nm]]
  }


  if (method=="boot") {
    res=.computePDelta(formula,data,
                    intervention, forestControl)
  }
  # else if (method=="bootRF"){
  #   res=.computePDeltaRF(formula,data,
  #                     intervention, forestControl)
  # }

  else if (method=="normal"){
    res=.computePDeltaNormal(formula,data,
                        intervention, forestControl)
  }
  # else{
  #
  #   res=.computePDeltaBart(formula,data,
  #                   intervention,  forestControl)
  # }

  In=paste(all.vars(formula)[1],"1", sep="")
  cat("Computing PSICA tree...\n")

  formula1=as.formula(paste(In, "~",paste(attr(terms(formula), "term.labels"), collapse="+")))
  data1=data
  data1[[all.vars(formula)[1]]]=c()
  data1[[In]]=1:nrow(data)
  parms=list(Delta=res$data$Delta,P=res$data$P, m=res$nclasses, classes=res$classes, cost=cost, confidence=confidence, prune=prune)
  t1=createpsicaTree(
    formula1,data1, parms=parms,
    control=treeControl,...)
  t1$In=In
  t1$parms=parms
  t1$Y=data[[all.vars(formula)[1]]]
  t1$interv=data[[intervention]]

  return(t1)
}
#' Plots a PSICA tree.
#' @param x the PSICA tree to be plotted.
#' @param type If type=1 a default plot showing the probability of each treatment to be the best are
#' displayed, if type=2 boxplots of the effect values are displayed in each terminal node.
#' @param ... Other parameter passed to the generic plot function
#' @description
#' A plot is produced that shows a PSICA tree in which the terminal nodes
#' contain the probabilities that one treatment is better than the other
#' ones as well as a label containing the possible best treatments

#' @export
#' @importFrom graphics par boxplot plot
#' @importFrom stats predict
#' @importFrom gridBase gridFIG
#'
plot.psicaTree=function(x, type=1,...){
#  browser()
  pT<-x
  class(pT)<-"rpart"
  opar <- par(no.readonly = T )
  nm=pT$In
  Preds=predict(pT)
  assign(nm, Preds, envir=environment(pT$terms))
  pT1=partykit::as.party(pT)

  panel=function(node){
    id=partykit::id_node(node)
    vp1 <- grid::viewport( name = "vp1")
    grid::pushViewport(vp1)
    grid::grid.rect(gp = grid::gpar(col="white"))
    grid::grid.rect(width = grid::unit(0.9, "npc"), height = grid::unit(0.9, "npc"))
    grid::grid.text(c(pT$LABELS[[pT$frame$yval[as.numeric(id)]]]))
    grid::popViewport()
  }
  panel2=function(node){

    id=partykit::id_node(node)
    vp1 <- grid::viewport( name = "vp1")
    grid::pushViewport(vp1)
    grid::grid.rect(gp = grid::gpar(col="white"))
    vp2=grid::viewport(width = grid::unit(0.6, "npc"), height=grid::unit(0.8, "npc"))
    grid::pushViewport(vp2)
    par(plt=gridBase::gridFIG(), new=TRUE)
    ind=which(Preds==pT$frame$yval[as.numeric(id)])
    Ys=pT$Y[ind]
    Cats=pT$interv[ind]
    boxplot(Ys~Cats, xlab="", ylab="")
    grid::popViewport()
    grid::popViewport()
  }

  #plotting
  if (type==1)
    plot(pT1, terminal_panel = panel,...)
  else if (type==2)
    plot(pT1, terminal_panel = panel2,...)

  par(opar)

}

#' Prints a PSICA tree.
#' @param x the PSICA tree to be printed
#' @param ... Other parameter passed to the generic print function
#' @export
#' @importFrom stats predict
print.psicaTree=function(x,...){
  pT<-x
  class(pT)<-"rpart"
  nm=pT$In
  pLabel1=lapply(pT$pLabel,FUN = function(x)paste("[",paste(x, collapse=",", sep=""), "]", sep=""))
  pLabel2=unlist(pLabel1)

  assign(nm, as.factor(pLabel2[predict(pT)]), envir=environment(pT$terms))

  pT1=partykit::as.party(pT)


  print(pT1,...)
}

#' Makes predictions for PSICA tree.
#' @param object the PSICA tree to be predicted
#' @param prob should the matrix or class probabilities be predicted (TRUE)
#' or the classes themselves (FALSE)
#' @param ... further parameters to be passed to 'rpart' object.
#' @export
#' @importFrom stats predict
predict.psicaTree=function(object,prob=FALSE,...){
  pT<-object
  if(!prob){
    pLabel1=lapply(pT$pLabel,FUN = function(x)paste("[",paste(x, collapse=",", sep=""), "]", sep=""))
    class(pT)<-"rpart"
    pLabel2=unlist(pLabel1)

    Prediction=as.factor(pLabel2[predict(pT,...)])
  } else{
    class(pT)<-"rpart"
    Prediction= pT$Probs[predict(pT,...)]
  }


  return(Prediction)
}

prune <- function (pT, ...) {
  UseMethod("prune", pT)
}

#' Prunes a PSICA tree.
#' @param pT the PSICA tree to be pruned.
#' @param cp Cost complexity parameter that defines how much the tree
#' must be pruned. See rpart::control() for description of the
#' cost-complexity parameter.
#' @param ... further parameters to be passed to 'rpart' object.
#' @export
prune.psicaTree=function(pT,cp=0.001,...){
  pT1=rpart::prune.rpart(pT, cp,...)
  class(pT1)<-c("psicaTree", class(pT1))
  return(pT1)
}

snip <- function (pT, ...) {
  UseMethod("prune", pT)
}


#' Prunes desired nodes in a PSICA tree.
#' @param pT the PSICA tree to be pruned.
#' @param indices vector of indices of the nodes. All contents below the nodes with these
#' indices will be pruned.
#' @param ... further parameters to be passed to 'rpart' object.
#' @export

snip.psicaTree<-function(pT, indices){
  tree1=rpart::snip.rpart(pT,indices)
  return(tree1)
}
