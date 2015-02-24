

######## general functions for the 'gamMD' class: print summary and plot
#' print.gamMD
#'
#' print method for gamMD
#' 
#' @export
print.gamMD<-function(object){
  if(any(names(object)=="fitG")) {
    if (object$converged) cat("algorithm has converged in ",object$iter," iterations","\n\n")
    if (!object$converged) cat("algorithm has not converged","\n\n")
    # print.default("\n\n")
    if(!is.null(object$fitM$sm.p)) cat("lambdaM is",round(as.vector(object$fitM$sm.p),3),"\n")
      else cat("only parametric components in the mean \n")
    cat("the e.d.f. for the mean is",round(object$fitM$df,3),"\n\n")
    if(!is.null(object$fitG$sm.p)) cat("lambdaG is",round(as.vector(object$fitG$sm.p),3),"\n")
      else cat("only parametric components in the dispersion \n")
    cat("the e.d.f. for the dispersion is",round(object$fitG$df,3),"\n")
  }  
  if(!any(names(object)=="fitG")) {
    if (object$converged) cat("algorithm has converged ",object$iter," iterations","\n\n")
    if (!object$converged) cat("algorithm has not converged","\n\n")
    if(!is.null(object$fitM$sm.p)) cat("lambdaM is",round(as.vector(object$fitM$sm.p),3),"\n")
      else cat("only parametric components in the mean \n")
    cat("the e.d.f. for the mean is",round(object$df,3),"\n\n")
  }
}

#' summary.gamMD
#'
#' summary method for gamMD
#' @export
summary.gamMD<-function(object, ...){
  if(any(names(object)=="fitG")) {
    if (object$converged) cat("algorithm has converged in ",object$iter," iterations","\n\n")
    if (!object$converged) cat("algorithm has not converged","\n\n")
    cat('     mean estimation:',"\n")
    for( i in 1:length(object$fitM$vecdf)) cat('for',paste('x',i,sep=''),"lambdaM is",round(as.vector(object$fitM$sm.p)[i],3),"corresponding to",round(object$fitM$vecdf[i],3),'degrees of freedom',"\n")
    cat('     dispersion estimation:',"\n")
    for( i in 1:length(object$fitG$vecdf)) cat('for',paste('x',i,sep=''),"lambdaG is",round(as.vector(object$fitG$sm.p)[i],3),"corresponding to",round(object$fitG$vecdf[i],3),'degrees of freedom',"\n")
  }
  if(!any(names(object)=="fitG")) {
    if (object$converged) cat("algorithm has converged ",object$iter," iterations","\n\n")
    if (!object$converged) cat("algorithm has not converged","\n\n")
    cat('     mean estimation:',"\n")
    for( i in 1:length(object$vecdf)) cat('for',paste('x',i,sep=''),"lambdaM is",as.vector(object$sm.p)[i],"corresponding to",object$vecdf[i],'degrees of freedom',"\n")
  }
}




### needed for plot
read.formPLOT<-function(ff,data){
  ### same as read form in the beginning, this just takes out the covariates singularly
  int<-1
  ff<-as.formula(ff)
  if (missing(data)) data <- environment(ff)
  else {ddaatt <- NULL; data <- data[,all.vars(ff)]}
  ft<-terms.formula(ff,specials=c("bsp","I"))
  termsl<-attr(ft,"term.labels")
  nt<-length(termsl)
  dimsP<-NULL  #### this should be included in the dims part of the real reading in gam.Intern
  lNAl<-NULL       ### I know this is not a good way to program this, should be improved, but missing values are nasty...
  for(k in 1:length(all.vars(ff))){
    oo<-eval(parse(text=all.vars(ff)[k]),envir=data)
    lNAl<-unique(c(lNAl,seq(1,length(oo))[is.na(oo)]))
  }
  if(length(lNAl)) {
    listOR<-list()
    for(k in 1:length(all.vars(ff))){
      nn<-all.vars(ff)[k]
      oo<-eval(parse(text=all.vars(ff)[k]),envir=data)
      listOR[[k]]<-oo
      if(is.environment(data)) assign(nn,oo[-lNAl],envir=data)
      else  ddaatt<-cbind(ddaatt,oo[-lNAl])
    }
  }
  else{ if(!is.environment(data)) ddaatt<-data}
  if(!is.environment(data)) {
    data<-as.data.frame(ddaatt)
    colnames(data)<-all.vars(ff)
  }
  y<-eval(parse(text=attr(ft,"variables")[2]),envir=data);n<-NROW(y)
  indsmooth<-attr(ft,"specials")$bsp-1;lensmooth<-length(indsmooth)
  indNsmoothI<-attr(ft,"specials")$I-1
  indNsmooth<-seq(1,nt)
  for(i in sort(c(indsmooth,indNsmoothI),decreasing=T)) indNsmooth<-indNsmooth[-(i)]
  lenNsmooth<-length(indNsmooth)
  dataMAT<-names<-NULL
  if(length(indNsmooth)){ ## building the parameteric part
    allParsInd<-parMat<-names<-allfirstOrd<-allhigherOrd<-NULL
    j<-0
    for(i in (indNsmooth)){
      j<-j+1
      firstOrd<-termsl[i]
      obj<-eval(as.expression(parse(text=firstOrd)),envir=data);names<-c(names,firstOrd)
      allfirstOrd<-c(allfirstOrd,firstOrd)
      allParsInd<-c(allParsInd,j)
    }
    allPars<-c(allfirstOrd,allhigherOrd);uniquePars<-unique(c(allfirstOrd,allhigherOrd))
    for(i in 1:length(uniquePars)){dataMAT<-cbind(dataMAT,eval(parse(text=uniquePars[i]),envir=data))}
  }
  j<-1 ## building the NON parameteric part
  for(i in (indsmooth)){
    a<-(strsplit(strsplit(termsl[i],',')[[1]][1],')'))
    b<-substr(a,start=5,stop=nchar(a))
    obj<-eval(as.expression(parse(text=b)),envir=data)
    dataMAT<-cbind(dataMAT,obj);names<-c(names,b)
  }
  if(!attr(ft,"intercept")) {
    int<-0
  }
  if(length(lNAl) & is.environment(data)) { ### putting missing values back, when needed
    for(k in 1:length(all.vars(ff))){
      nn<-all.vars(ff)[k]
      assign(nn,listOR[[k]],envir=data)
    }
  }
  if(!length(lNAl)) lNAl<-NULL # FALSE
  colnames(dataMAT)<-names
  list(dataMAT=dataMAT,int=int)
}



#' @name gamMD
#' @title Double Generalised Additive Models class
#'
#' @description Typically the result of a DoubleGam or DoubleRobGam fit. Method functions as \code{summary}, \code{print} and \code{plot}. 
#' Depending on whether the mean only or both the mean and the dispersion function are estimated the results will be different.
#' In both cases the information on the convergence status og the overall estimation is given by \code{converged}.
#' \code{iter} gives the total number of (outer) iterations; 
#' \code{data} contains the dataset to which the model was fitted, if given.
#' 
#' If both the mean and the dispersion function are estimated the object will contain the following elements: 
#' \code{converged}, \code{convVec}, \code{data}, \code{fitG}, \code{fitM}, \code{iter}, \code{relE}.
#' \code{fitM} contains information on the mean estimation procedure;  
#' \code{fitG} contains information on the dispersion estimation procedure and has the same elements as \code{fitM}. 
#' 
#' If only the mean function is estimated all the values of the \code{fitM} object are given in \code{gamMD} object.  
#' 
#' @param object an object of the gamMD class
#' @param col colour in which the lines in the plot should be drawn
#' @param xlab,ylab optional axis labs
#' @param one logic value to indicate whether mean and disperion figures should be included togther in one plot
#' @param conf.level confidence level at which confidence bands should be drawn (default is 0.05)
#' @param ci.plot logic value to indicate whether confidence bands should be drawn
#' @seealso \code{\link{DoubleRobGam}}, \code{\link{DoubleGam}}
#' @export
plot.gamMD<-function(object, col=1, xlab=NULL, ylab="yt", one=TRUE, conf.level=0.05, ci.plot = TRUE, ...){
  ############ this now works and plots both parametric and nonparametric components
  ##### it has a lot of if conditions and it could be greatly improved, but it does its jobs, and plotting the parametric parts was not easy
  dimsM<-c(c(object$dimsP,object$dims[-1]),c(object$fitM$dimsP,object$fitM$dims[-1]))
  if(dimsM[1]==0) dimsM<-dimsM[-1]; if(dimsM[length(dimsM)]==0) dimsM<-dimsM[-length(dimsM)]
  dimsG<-NULL
  cutoff<-qnorm((1-conf.level/2),0,1)
  if(any(names(object)=="fitG"))  {dimsG<-c(object$fitG$dimsP,object$fitG$dims[-1]) 
                                   if(dimsG[1]==0) dimsG<-dimsG[-1]; if(dimsG[length(dimsG)]==0) dimsG<-dimsG[-length(dimsG)]}
  nM<-length(dimsM);nG<-length(dimsG);nPlots<-max(length(dimsM),length(dimsG))
  nn<-nrow(cbind(object$fitM$desMat,object$desMat))
  if(one) par(mfrow=c(2,trunc((nPlots))))
  if(!one | !any(names(object)=="fitG")) par(mfrow=c(1,trunc(nM)))
  if(dimsM[1]>0){
    for(i in 2:(1+nM)){
      fM<-c(object$fitM$formula,object$formula)[[1]]
      if(NROW(object$data)==1) { ff<-read.formPLOT(fM)}
      if(NROW(object$data)!=1) { ff<-read.formPLOT(fM,data=object$data)}
      #         subseq<-seq(max(1+sum(dimsM[0:(i-1)])),sum(dimsM[0:(i)]))
      ddM<-c(ff$int,dimsM)
      subseq<-((1+(sum(ddM[0:(i-1)]))):((sum(ddM[0:i]))))
      nametemp<-colnames(ff$dataMAT)[(i-1)]
      xtemp<-ff$dataMAT[,(i-1)]
      ytemp<-c(object$fitM$yt,object$yt)
      BB<-as.matrix(cbind(object$fitM$desMat[,subseq],object$desMat[,subseq]))
      aa<-as.vector(c(object$fitM$coef[subseq],object$coef[subseq]))
      ci<-matrix(unlist(list(object$cov.coef[subseq,subseq],object$fitM$cov.coef[subseq,subseq])),ncol=length(subseq))
      lintemp<-BB%*%aa;rangetemp<-range(c(lintemp[is.finite(lintemp)],ytemp[is.finite(ytemp)]));ytemp[!is.finite(ytemp)]<-min(rangetemp)
      plot(xtemp[order(xtemp)],ytemp[order(xtemp)],pch=19,col="grey",bty="l",cex=0.85,xlab=nametemp,ylab=ylab,ylim=rangetemp)
      lines(xtemp[order(xtemp)],(lintemp)[order(xtemp)],col=col,lwd=2)
      if(ci.plot){
        lines(xtemp[order(xtemp)],(lintemp+cutoff*sqrt(diag(BB%*%ci%*%t(BB))))[order(xtemp)],col=col,lwd=2,lty=2)
        lines(xtemp[order(xtemp)],(lintemp-cutoff*sqrt(diag(BB%*%ci%*%t(BB))))[order(xtemp)],col=col,lwd=2,lty=2)
      }
    }
  }
  if(any(names(object)=="fitG")){
    if(nM < nPlots) par(mfg=c(2,1))
    if(!one) {dev.new();par(mfrow=c(1,trunc((nG))))}
    for(i in 2:(1+nG)) {
      yd<-rep(0,length=nn)
      Gform<-formula(object$fitG$formula)
      if (length(Gform)==2){
        Gform[3] <- Gform[2];Gform[2] <- formula(yd~x)[2]
      }
      if(NROW(object$data)==1) {data <- environment(Gform);}
      if(NROW(object$data)!=1) {data <- object$data;}
      if(is.environment(data)) { assign("yd",yd,envir=data);}
      if(!is.environment(data)) {n2<-names(data);data<-as.data.frame(cbind(yd,data));names(data)<-c("yd",n2)}
      #### end of preparation to read the formula
      if(NROW(object$data)==1) { ff<-read.formPLOT(Gform);rm(yd,envir=data)} ### data is not given
      if(NROW(object$data)!=1) { ff<-read.formPLOT(Gform,data=data)}
      ddG<-c(ff$int,dimsG)
      subseq<-((1+(sum(ddG[0:(i-1)]))):((sum(ddG[0:i]))))
      nametemp<-colnames(ff$dataMAT)[(i-1)]
      xtemp<-ff$dataMAT[,(i-1)]
      plot(xtemp[order(xtemp)],object$fitG$yt[order(xtemp)],pch=19,col="grey",bty="l",cex=0.85,xlab=nametemp,ylab="log(dev)")
      BB<-as.matrix(object$fitG$desMat[,subseq])
      aa<-as.vector(object$fitG$coef[subseq])
      ci<-object$fitG$cov.coef[subseq,subseq]
      lines(xtemp[order(xtemp)],(BB%*%aa)[order(xtemp)],col=col,lwd=2)
      if(ci.plot){
        lines(xtemp[order(xtemp)],(BB%*%aa+cutoff*sqrt(diag(BB%*%ci%*%t(BB))))[order(xtemp)],col=col,lwd=2,lty=2)
        lines(xtemp[order(xtemp)],(BB%*%aa-cutoff*sqrt(diag(BB%*%ci%*%t(BB))))[order(xtemp)],col=col,lwd=2,lty=2)
      }
    }
  }
}
