######## The DoubleGam function was build at the same time as the DoubleRobGam, the two functions are built up in the same way
    ####### you need the general functions.R file as well

###### some parts come from the glmrob function in the robustbase library



#### In this file, first you find all the functions to estimate the mean function
      ######### then all the function for the dispersion estimation
      ######### the Double part comes at the end

######## All functions written by Ilaria Prosdocimi (Ilaria.Prosdocimi  AT gmail DOT com)


############### mean estimation


gamMeanGivenSP<-function(sm.p,y,XM,Pstr,dimsS=dimsS,dimsT=dimsT,start=NULL,family="gaussian",trace=FALSE,acc=10^-6,maxit=30,scale=1,conf.type="bayesian"){
### fits the model for the given smoothing parameters
    nobs <- NROW(y);  ncoef <- ncol(XM)
    sm.ps<-sm.p
    weights<-rep(1,l=nobs)
    if (is.character(family))
	family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
	family <- family()
    variance <- family$variance
    linkinv <- family$linkinv
    if (!is.function(variance) || !is.function(linkinv)) stop("illegal 'family' argument")
    mu.eta <- family$mu.eta
    if (is.null(valideta <- family$valideta)) valideta <- function(eta) TRUE
    if (is.null(validmu	 <- family$validmu))  validmu <-  function(mu) TRUE
    etastart<-mustart<-NULL
    eval(family$initialize) ## --> n, mustart, y and weights (=ni)
    if(is.null(start)) start<-c(mean(family$linkfun(mustart)),rep(0.001,l=(ncoef-1))) 
    betaOld <- beta <- as.vector(start) # as.v*(): dropping attributes
    eta <- etastart <-  as.vector(XM %*% beta); mustart <- mu <- linkinv(eta) # mu estimates pi (in [0,1]) at the binomial model
    if (!(validmu(mu) && valideta(eta)))
           stop("Cannot find valid starting values: You need help")
    ni <- as.vector(weights) # dropping attributes for computation
    if(dimsT[1] != 0) supportMat<-matrix(0,nrow=nrow(Pstr),ncol=dimsT[1])
         else supportMat<-NULL
    if(length(sm.ps) >= 1)
    for(i in 1:length(sm.ps)) supportMat<-cbind(supportMat,matNum(sm.ps[i],dim1=nrow(Pstr),dim2=dimsS[i]))
    P<-supportMat*Pstr
    ### Iterations for the mean
    if(trace && ncoef) {
         cat("The iteration begins: \n")
    }
    conv <- FALSE
    if(length(sm.ps)==0) inversion<-function(...) {solve(...)} ##### when the model is only parametric use solve()
    if(length(sm.ps)!=0) inversion<-function(...) {ginv(...)}  ##### otherwise we have results different from glm()
    if(ncoef) for (nit in 1:maxit) {
        W<- diag(variance(mu)/scale); z<- eta+(y-mu)/variance(mu)
        beta<- inversion(t(XM)%*%W%*%XM+P)%*%t(XM)%*%W%*%z
        Dbeta<-betaOld-beta ###  I use the same stopping rule (code) as in robustbase
        eta <- as.vector(XM %*% beta); mu <- linkinv(eta)
          ## Check convergence: relative error < tolerance
          relE <- sqrt(sum(Dbeta^2)/max(1e-20, sum(betaOld^2)))
          conv <- relE <= acc
          if(trace) cat("at iteration ", nit, "the relative change in the parameter is ",relE,'\n')
          if(conv) break
        betaOld <- beta
    } ## end of iteration
    else { ## ncoef == 0
          conv <- TRUE
        nit <- 0
    }
    if (!conv) 	warning("Algorithm did not converge")
    eps <- 10 * .Machine$double.eps
    names(mu) <- names(eta) <- names(beta) # re-add after computation
    Dev<-family$dev.resid(y=y,mu=mu,wt=ni)/scale
    df<- sum(diag(XM%*%inversion(t(XM)%*%W%*%XM+P)%*%t(XM)%*%W))
    GCV<-nobs*sum(Dev)/sum((nobs-df)^2)
    AIC<-sum(Dev)+2*df
    mid<-family$linkfun(y)
    mid<-mean(mid[is.finite(mid) & !is.na(mid)])
    yt<-family$linkfun(y)-mid
#     ym<-y;ym[ym==0]<-10^-5 #### to fix problems when zero are present for poisson and binomial, yt only gets used in the plots...
#     mid<-family$linkfun(ym)
#     mid<-mean(mid[is.finite(mid) & !is.na(mid)])
#     yt<-family$linkfun(ym)-mid
    vecdfM<-NULL
    if(dimsT[2]!= 0){
       for(i in 2:length(dimsT)){
       mm<-XM[,(1+sum(dimsT[0:(i-1)])):(sum(dimsT[0:i]))]
       pm<-P[(1+sum(dimsT[0:(i-1)])):(sum(dimsT[0:i])),(1+sum(dimsT[0:(i-1)])):(sum(dimsT[0:i]))]
       vecdfM<-c(vecdfM,sum(diag(mm%*%inversion(t(mm)%*%W%*%mm+pm)%*%t(mm)%*%W)))
       }
    }
#     scale<-rep(.2,l=nobs);W<-diag(variance(mu))
#     print((diag((inversion(t(XM)%*%(diag(rep(scale,length=nobs))*W)%*%XM+P)%*%t(XM)%*%(diag(rep(scale^2,length=nobs))*W)%*%XM%*%inversion(t(XM)%*%(diag(rep(scale,length=nobs))*W)%*%XM+P))) )/(diag(inversion(t(XM)%*%W%*%XM+P)%*%t(XM)%*%W%*%XM%*%inversion(t(XM)%*%W%*%XM+P))))
    if(conf.type=="naive") asCov<-inversion(t(XM)%*%W%*%XM+P)%*%t(XM)%*%(diag(rep(scale,length=nobs))*W)%*%XM%*%inversion(t(XM)%*%W%*%XM+P)
    if(conf.type=="bayesian") asCov<-inversion(t(XM)%*%W%*%XM+P)
    list(coefficients = as.vector(beta), fitted.values = mu,desMat=XM,dims=dimsT,dispersion = scale, family = family,
          linear.predictors = eta, deviance = Dev, residuals= (y-mu)/sqrt(variance(mu)) ,s.resid= (y-mu)/sqrt(scale*variance(mu)),
          iter=nit,y=y,converged=conv,sm.p=sm.ps,GCV=GCV,AIC=AIC,yt=yt,df=df,vecdf=vecdfM,cov.coef=asCov)
}


gamMeanGivenSPsmall<-function(sm.p,y,XM,Pstr,dimsS=dimsS,dimsT=dimsT,start=NULL,family,trace=FALSE,acc=10^-6,maxit=30,scale=1){
### fits the model for the given smoothing parameters
### smaller output
    y<-y
    nobs <- NROW(y);  ncoef <- ncol(XM)
    sm.ps<-sm.p
    weights<-rep(1,l=nobs)
    if (is.character(family))
	family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
	family <- family()
    variance <- family$variance
    linkinv <- family$linkinv
    if (!is.function(variance) || !is.function(linkinv)) stop("illegal 'family' argument")
    mu.eta <- family$mu.eta
    if (is.null(valideta <- family$valideta)) valideta <- function(eta) TRUE
    if (is.null(validmu	 <- family$validmu))  validmu <-  function(mu) TRUE
    etastart<-mustart<-NULL
    eval(family$initialize) ## --> n, mustart, y and weights (=ni)
    if(is.null(start)) start<-c(mean(family$linkfun(mustart)),rep(0.001,l=(ncoef-1))) 
    betaOld <- beta <- as.vector(start) # as.v*(): dropping attributes
    eta <- etastart <- as.vector(XM %*% beta); mu <- linkinv(eta) # mu estimates pi (in [0,1]) at the binomial model
    if (!(validmu(mu) && valideta(eta)))
	stop("Cannot find valid starting values: You need help")
    ni <- as.vector(weights) # dropping attributes for computation
    if(dimsT[1] != 0) supportMat<-matrix(0,nrow=nrow(Pstr),ncol=dimsT[1])
         else supportMat<-NULL
    if(length(sm.ps) >= 1) 
    for(i in 1:length(sm.ps)) supportMat<-cbind(supportMat,matNum(sm.ps[i],dim1=nrow(Pstr),dim2=dimsS[i]))
    P<-supportMat*Pstr
    ### Iterations for the mean
    if(trace && ncoef) {
         cat("The iteration begins: \n")
    }
    conv <- FALSE
    if(length(sm.ps)==0) inversion<-function(...) {solve(...)} ##### when the modle is only parametric use solve()
    if(length(sm.ps)!=0) inversion<-function(...) {ginv(...)}  ##### otherwise w ehave results different from glm()
    if(ncoef) for (nit in 1:maxit) {
        W<-   diag(variance(mu)/scale); z<- eta+(y-mu)/variance(mu)
        beta<- ginv(t(XM)%*%W%*%XM+P)%*%t(XM)%*%W%*%z
        Dbeta<-betaOld-beta ###  I use the same stopping rule (code) as in robustbase
        eta <- as.vector(XM %*% beta); mu <- linkinv(eta)
	  ## Check convergence: relative error < tolerance
        relE <- sqrt(sum(Dbeta^2)/max(1e-20, sum(betaOld^2)))
        conv <- relE <= acc
        if(trace) cat("at iteration ", nit, "the relative change in the parameter is ",relE,'\n')
	if(conv) break
        betaOld <- beta
    } ## end of iteration
    else { ## ncoef == 0
	conv <- TRUE
	nit <- 0
    }
    eps <- 10 * .Machine$double.eps
    names(mu) <- names(eta) <- names(beta) # re-add after computation
    Dev<-family$dev.resid(y=y,mu=mu,wt=ni)/scale
    df<- sum(diag(XM%*%inversion(t(XM)%*%W%*%XM+P)%*%t(XM)%*%W))
    GCV<-nobs*sum(Dev)/sum((nobs-df)^2)
    AIC<-sum(Dev)+2*df
    mid<-family$linkfun(y)
    mid<-mean(mid[is.finite(mid) & !is.na(mid)])
    yt<-family$linkfun(y)-mid
    list(GCV=GCV,AIC=AIC)
}



gamMeanIC<-function(lsm.p,y,XM,Pstr,dimsS,dimsT,start = NULL,family,acc = 10^-6, maxit = 30, crit="GCV",scale=1){
### computes the information criterion of the mode for the given smoothing parameters
    if(crit != "GCV" & crit!="AIC") {print('can not recognize the criterion, use GCV');crit<-"GCV"}
    gg<-gamMeanGivenSPsmall(sm.p=exp(lsm.p),y=y,XM=XM,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,family=family,
    start=start,trace=FALSE,acc=acc,maxit=maxit,scale=scale)
    if(crit=="GCV") res<-gg$GCV
    if(crit=="AIC") res<-gg$AIC
res
}




######## This one collects everything
gamMeanIntern<-function(sm.p=sm.ps,y,XM,Pstr,dimsS,dimsT,start=NULL,family="gaussian",trace=FALSE,acc=10^-6,maxit=30,selection="GCV",data,scale=1,minlambda=0.001,maxlambda=5000,conf.type){
    sm.pinF<-sm.p
    if(selection!="none" & selection!="GCV" & selection!="AIC" ){
    print('how should I select the smoothing parameters? default is GCV')
    selection<-"GCV"
    }
    if(selection=="none") sm.ps<-sm.pinF
#      else sm.ps<-exp(nlm(p=log(sm.pinF),gamMeanIC,y=y,XM=XM,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,start=start,
#                     family=family,acc=acc,maxit=maxit,scale=scale,crit=selection)$est) )   ###### with nlm I can not give upper and lower limits
     else sm.ps<-exp(optim(par = log(sm.pinF),gamMeanIC,y=y,XM=XM,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,start=start,family=family,acc=acc,
               maxit=maxit,scale=scale,crit=selection,method="L-BFGS-B",lower=log(minlambda),upper=log(maxlambda))$par)
   if(trace==TRUE) cat('the smoothing parameters are', sm.ps, '\n')
   fit<-gamMeanGivenSP(sm.p=sm.ps,y=y,XM=XM,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,start=start,family=family,trace=trace,
         acc=acc,maxit=maxit,scale=scale,conf.type=conf.type)
   fit
}





gamMean<-function(formula,start=NULL,trace=FALSE,acc=10^-6,maxit=30,family,selection="GCV",data,scale=1){
    fread<-read.form(formula,data=data)
    sm.pinF<-fread$sm.p
    XM<-fread$dataMAT
    params<-smooths<-NULL
#     if(length(fread$params)) for( j  in 1:length(fread$params)) paramts<-cbind(params,fread$params[[j]])
    if(length(fread$smooths)) for( j  in 1:length(fread$smooths)) smooths<-cbind(smooths,fread$smooths[[j]])
    intercept<-fread$int
    Pstr<-fread$P;y<-fread$y
    nobs <- NROW(y);  ncoef <- ncol(XM)
    dimsS<-if(!length(fread$smooths)) 0 else sapply(fread$smooths,ncol)
    dimsT<-c(ncol(Pstr)-sum(dimsS),dimsS)
    weights<-rep(1,l=nobs)
    if(selection!="none" & selection!="GCV" & selection!="AIC" ){
    print('how should I select the smoothing parameters? default is GCV')
    selection<-"GCV"
    }
    if(selection=="none") sm.ps<-sm.pinF
     else sm.ps<-exp(nlm(p=log(sm.pinF),gamMeanIC,y=y,XM=XM,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,start=start,
                    acc=acc,maxit=maxit,crit=selection,scale=scale)$est)
   if(trace==TRUE) cat('the smoothing parameters are', sm.ps, '\n')
   fit<-gamMeanGivenSP(sm.p=sm.ps,y=y,XM=XM,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,start=start,
               trace=trace,acc=acc,maxit=maxit,scale=scale)
   fit$seqNA<-fread$seqNA
   fit$dimsP<-fread$dimsP
   fit
}






############## dispersion estimation


gamDispGivenSP<-function(sm.p,y,XG,Pstr,dimsS=dimsS,dimsT=dimsT,start=NULL,trace=FALSE,acc=10^-6,maxit=30,link="log",conf.type){
### fits the model for the given smoothing parameters
    yd<-y
    nobs <- NROW(yd);  ncoef <- ncol(XG)
    sm.ps<-sm.p
    weights<-rep(1,l=nobs)
    family<-Gamma(link="log") #### if I give the option of having differnt links, I should change W as well!!
    if (is.character(family))
	family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
	family <- family()
    variance <- family$variance
    linkinv <- family$linkinv
    if (!is.function(variance) || !is.function(linkinv)) stop("illegal 'family' argument")
    gamma.xi <- family$gamma.xi
    if (is.null(validxi <- family$validxi)) validxi <- function(xi) TRUE
    if (is.null(validgamma <- family$validgamma))  validgamma <-  function(gamma) TRUE
    eval(family$initialize) ## --> n, gammastart, y and weights (=ni)
    if(is.null(start)) start<-rep(family$linkfun(mean(y)),l=ncoef)
    deltaOld <- delta <- as.vector(start) # as.v*(): dropping attributes
    xi <- as.vector(XG %*% delta); gamma <- linkinv(xi) # gamma estimates pi (in [0,1]) at the binomial model
    ni <- as.vector(weights) # dropping attributes for computation
    if (!(validgamma(gamma) && validxi(xi)))
	stop("Cannot find valid starting values: You need help")
     #### for now I have fixed smoothing parameters
     #### I hope/should be able to chooes them optimally and put this in the iteration
    if(dimsT[1] != 0) supportMat<-matrix(0,nrow=nrow(Pstr),ncol=dimsT[1])
         else supportMat<-NULL
    matNum<-function(num,dim1,dim2) matrix(num,ncol=dim2,nrow=dim1)
    if(length(sm.ps) >= 1) 
    for(i in 1:length(sm.ps)) supportMat<-cbind(supportMat,matNum(sm.ps[i],dim1=nrow(Pstr),dim2=dimsS[i]))
    P<-supportMat*Pstr
    ### Iterations for the mean
    if(trace && ncoef) {
         cat("The iteration begins: \n")
    }
    conv <- FALSE
    if(length(sm.ps)==0) inversion<-function(...) {solve(...)} ##### when the model is only parametric use solve()
    if(length(sm.ps)!=0) inversion<-function(...) {ginv(...)}  ##### otherwise we have results different from glm()
    if(ncoef) for (nit in 1:maxit) {
        W<-   diag(rep.int(0.5,nobs)); q<- xi+(y-gamma)/gamma
        delta<- inversion(t(XG)%*%W%*%XG+P)%*%t(XG)%*%W%*%q
        Ddelta<-deltaOld-delta ###  I use the same stopping rule (code) as in robustbase
	xi <- as.vector(XG %*% delta);	gamma <- linkinv(xi)
	## Check convergence: relative error < tolerance
	relE <- sqrt(sum(Ddelta^2)/max(1e-20, sum(deltaOld^2)))
	conv <- relE <= acc
        if(trace) cat("at iteration ", nit, "the relative change in the parameter is ",relE,'\n')
	if(conv) break
	deltaOld <- delta
    } ## end of iteration
    else { ## ncoef == 0
	conv <- TRUE
	nit <- 0
    }
    if (!conv) 	warning("Algorithm did not converge")
    eps <- 10 * .Machine$double.eps
    names(gamma) <- names(xi) <- names(delta) # re-add after computation
    Dev<-family$dev.resid(y=y,mu=gamma,wt=ni)/2
    df<- sum(diag(XG%*%inversion(t(XG)%*%W%*%XG+P)%*%t(XG)%*%W))
    GCV<-nobs*sum(Dev)/sum((nobs-df)^2)
    AIC<-sum(Dev)+2*df
    mid<-family$linkfun(y)
    mid<-mean(mid[is.finite(mid) & !is.na(mid)])
    yt<-family$linkfun(y)-mid
    vecdfG<-NULL
    if(dimsT[2]!= 0){
    for(i in 2:length(dimsT)){
        mm<-XG[,(1+sum(dimsT[0:(i-1)])):(sum(dimsT[0:i]))]
        pm<-P[(1+sum(dimsT[0:(i-1)])):(sum(dimsT[0:i])),(1+sum(dimsT[0:(i-1)])):(sum(dimsT[0:i]))]
        vecdfG<-c(vecdfG,sum(diag(mm%*%inversion(t(mm)%*%W%*%mm+pm)%*%t(mm)%*%W)))
        }
    }
    if(conf.type=="naive") asCov<-inversion(t(XG)%*%diag(0.5,nobs)%*%XG+P)%*%t(XG)%*%diag(0.5,nobs)%*%XG%*%inversion(t(XG)%*%diag(0.5,nobs)%*%XG+P)
    if(conf.type=="bayesian") asCov<-inversion(t(XG)%*%diag(0.5,nobs)%*%XG+P)
    list(coefficients = as.vector(delta), fitted.values = gamma,desMat=XG,dims=dimsT,dispersion = 2, family = family,
           linear.predictors = xi, deviance = Dev, residuals= (yd-gamma)/sqrt(gamma), s.resid= (yd-gamma)/sqrt(2*gamma),
           yd = yd, converged = conv, sm.p=sm.ps,GCV=GCV,AIC=AIC,yt=yt,df=df,vecdf=vecdfG,cov.coef=asCov)
  #,relE=relE
}
# iter = nit, 


gamDispGivenSPsmall<-function(sm.p,y,XG,Pstr,dimsS,dimsT,start=NULL,trace=FALSE,acc=10^-6,maxit=30,link="log"){
### fits the model for the given smoothing parameters
    yd<-y
    nobs <- NROW(yd);  ncoef <- ncol(XG)
    sm.ps<-sm.p
    weights<-rep(1,l=nobs)
    family<-Gamma(link="log") #### if I give the option of having differnt links, I should change W as well!!
    if (is.character(family))
	family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
	family <- family()
    variance <- family$variance
    linkinv <- family$linkinv
    if (!is.function(variance) || !is.function(linkinv)) stop("illegal 'family' argument")
    gamma.xi <- family$gamma.xi
    if (is.null(validxi <- family$validxi)) validxi <- function(xi) TRUE
    if (is.null(validgamma	 <- family$validgamma))  validgamma <-  function(gamma) TRUE
    eval(family$initialize) ## --> n, gammastart, y and weights (=ni)
    if(is.null(start)) start<-rep(family$linkfun(mean(y)),l=ncoef)
    deltaOld <- delta <- as.vector(start) # as.v*(): dropping attributes
    xi <- as.vector(XG %*% delta); gamma <- linkinv(xi) # gamma estimates pi (in [0,1]) at the binomial model
    ni <- as.vector(weights) # dropping attributes for computation
    if (!(validgamma(gamma) && validxi(xi)))
	stop("Cannot find valid starting values: You need help")
    if(dimsT[1] != 0) supportMat<-matrix(0,nrow=nrow(Pstr),ncol=dimsT[1])
         else supportMat<-NULL
    matNum<-function(num,dim1,dim2) matrix(num,ncol=dim2,nrow=dim1)
    if(length(sm.ps) >= 1) 
    for(i in 1:length(sm.ps)) supportMat<-cbind(supportMat,matNum(sm.ps[i],dim1=nrow(Pstr),dim2=dimsS[i]))
    P<-supportMat*Pstr
    if(trace && ncoef) {
         cat("The iteration begins: \n")
    }
    conv <- FALSE
    if(length(sm.ps)==0) inversion<-function(...) {solve(...)} ##### when the model is only parametric use solve()
    if(length(sm.ps)!=0) inversion<-function(...) {ginv(...)}  ##### otherwise we have results different from glm()
    if(ncoef) for (nit in 1:maxit) {
        W<-   diag(rep.int(0.5,nobs)); q<- xi+(y-gamma)/gamma
        delta<- inversion(t(XG)%*%W%*%XG+P)%*%t(XG)%*%W%*%q
        Ddelta<-deltaOld-delta ###  I use the same stopping rule (code) as in robustbase
	xi <- as.vector(XG %*% delta);	gamma <- linkinv(xi)
	## Check convergence: relative error < tolerance
	relE <- sqrt(sum(Ddelta^2)/max(1e-20, sum(deltaOld^2)))
	conv <- relE <= acc
        if(trace) cat("at iteration ", nit, "the relative change in the parameter is ",relE,'\n')
	if(conv) break
	deltaOld <- delta
    } ## end of iteration
    else { ## ncoef == 0
	conv <- TRUE
	nit <- 0
    }
    eps <- 10 * .Machine$double.eps
    names(gamma) <- names(xi) <- names(delta) # re-add after computation
    Dev<-family$dev.resid(y=y,mu=gamma,wt=ni)/2
    df<- sum(diag(XG%*%inversion(t(XG)%*%W%*%XG+P)%*%t(XG)%*%W))
    GCV<-nobs*sum(Dev)/sum((nobs-df)^2)
    AIC<-sum(Dev)+2*df
    mid<-family$linkfun(y)
    mid<-mean(mid[is.finite(mid) & !is.na(mid)])
    yt<-family$linkfun(y)-mid
    list(GCV=GCV,AIC=AIC)
}



gamDispIC<-function(lsm.p,y,XG,Pstr,dimsS,dimsT,start = NULL,acc = 10^-6, maxit = 30, crit="GCV"){
### computes the information criterion of the mode for the given smoothing parameters
    if(crit != "GCV" & crit!="AIC") {print('can not recognize the criterion, use GCV');crit<-"GCV"}
    gg<-gamDispGivenSPsmall(sm.p=exp(lsm.p),y,XG,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,start=start,trace=FALSE,acc=acc,maxit=maxit)
    if(crit=="GCV") res<-gg$GCV
    if(crit=="AIC") res<-gg$AIC
res
}



###### this collects a bit everything 
gamDispIntern<-function(sm.p,yd,XG,Pstr,dimsS,dimsT,start=NULL,trace=FALSE,acc=10^-6,maxit=30,selection="GCV",minlambda=0.001,maxlambda=5000,conf.type){
    sm.pinF<-sm.p
    if(selection!="none" & selection!="GCV" & selection!="AIC" ){
    print('how should I select the smoothing parameters? default is GCV')
    selection<-"GCV"
    }
    if(selection=="none") sm.ps<-sm.pinF
     else 
#      sm.ps<-exp(nlm(p=log(sm.pinF),gamDispIC,y=yd,XG=XG,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,start=start,
#                     acc=acc,maxit=maxit,crit=selection)$est)   ###### with nlm I can not give upper and lower limits
       sm.ps<-exp(optim(par = log(sm.pinF),gamDispIC,y=yd,XG=XG,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,
           start=start,acc=acc,maxit=maxit,crit=selection,method="L-BFGS-B",lower=log(minlambda),upper=log(maxlambda))$par)
   if(trace==TRUE) cat('the smoothing parameters are', sm.ps, '\n')
   fit<-gamDispGivenSP(sm.p=sm.ps,y=yd,XG=XG,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,start=start,trace=trace,acc=acc,maxit=maxit,conf.type=conf.type)
   fit
}



gamDisp<-function(formula,start=NULL,trace=FALSE,acc=10^-6,maxit=30,selection="GCV",data){
    fread<-read.form(formula,data=data)
    sm.pinF<-fread$sm.p
    XG<-fread$dataMAT
    params<-smooths<-NULL
#     if(length(fread$params)) for( j  in 1:length(fread$params)) paramts<-cbind(params,fread$params[[j]])
    if(length(fread$smooths)) for( j  in 1:length(fread$smooths)) smooths<-cbind(smooths,fread$smooths[[j]])
    intercept<-fread$int
    Pstr<-fread$P;yd<-fread$y
    nobs <- NROW(y);  ncoef <- ncol(XG)
    dimsS<-if(!length(fread$smooths)) 0 else sapply(fread$smooths,ncol)
    dimsT<-c(ncol(Pstr)-sum(dimsS),dimsS)
    weights<-rep(1,l=nobs)
    if(selection!="none" & selection!="GCV" & selection!="AIC" ){
    print('how should I select the smoothing parameters? default is GCV')
    selection<-"GCV"
    }
    if(selection=="none") sm.ps<-sm.pinF
     else sm.ps<-exp(nlm(p=log(sm.pinF),gamDispIC,y=yd,XG=XG,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,start=start,
                    acc=acc,maxit=maxit,crit=selection)$est)
   if(trace==TRUE) cat('the smoothing parameters are', sm.ps, '\n')
   fit<-gamDispGivenSP(sm.p=sm.ps,y=yd,XG=XG,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,start=start,trace=trace,acc=acc,maxit=maxit)
   fit$seqNA<-fread$seqNA
   fit$dimsP<-fread$dimsP
   fit
}












######## Double estimation


#' Setting the DoubleGam fitting defaults 
#' 
#' This function is used internally in the DoubleGam fitting procedure. 
#' The function's parameters control the numerical properties of the Double fitting procedure. 
#' Care should be taken when changing the parameters default values, things can go wrong 
#' (but sometimes fitting doesn't work well with the initial parameters, so you'll need to fiddle aroud with this)
#' 
#' @param maxitM the maximum number of itearations allowed for each inner mean estimation procedures
#' @param maxitG the maximum number of itearations allowed for each inner dispersion estimation procedures
#' @param tol  the level of accuarancy required for each inner mean and dispersion estimation procedure to converge
#' @param acc the level of accuracy required for the outer iteration to converge 
#' @param maxitOUT the maximum number of outer iteration allowed in the Double estimation procedure
#' @param lambdaM range of acceptable smoothing parameter values for the mean estiamation procedure 
#' @param lambdaG range of acceptable smoothing parameter values for the dispersion estiamation procedure 
#' @return a list of control parameters
#' @export
DoubleGamControl<-function(maxitM=30,maxitG=30,tol=10^-5,acc=5*10^-3,maxitOUT=55,lambdaM=c(0.0001,35500),lambdaG=c(0.0001,13500)){
list(maxitM=maxitM,maxitG=maxitG,tol=tol,maxitOUT=maxitOUT,acc=acc,lambdaM=lambdaM,lambdaG=lambdaG)
}


#' Double Generalized Additive Models 
#' 
#' A function to estimate both the mean and the dispersion function using B-splines bases. 
#' 
#' @param formulaM the mean formula - typically of the type y ~ bs(x1) + bs(x2). Parametric models or mixtures of 
#' parametric and non parametric bases can be specified
#' @param formulaG the dispersion formula - as  for formulaM, but no response variable to be specified. 
#' \code{formulaG = ~ 1} corresponds to the estimation of the dispersion as a constant (i.e. standard GAM) 
#' @param family as in glm, can be "gaussian", "poisson" or "Binomial"
#' @param data dataset where variables are stored (optional)
#' @param startM starting value for the coefficients needed to estimate the mean function
#' @param startG starting value for the coefficients needed to estimate the dispersion function
#' @param method can be "quasi" (the default) or "pseudo", according to whether deviance or pearson's residuals should be used 
#' as response variables when fitting the dispersion function
#' @param selection the method used to select the smoothing parameter. can be "GCV" (the default), "AIC" or 
#' "none", in which case no selection is done and a smoothing parameter value should be provided in the \code{bs()} call.
#' @param weights as in glm -  to be used with care, since the data get always weigthed by the estimated variance function
#' @param control a list to control some behaviours of the estimation procedure, see {\link{DoubleGamControl}}
#' @param trace should a trace to follow the convergence be printed? Default is FALSE
#' @param scale similar to weigth, to be used if the dispersion function should be kept fixed
#' @param conf.type character describing how to compute confidence intervals for the smooth components; can "bayesian" (default) or "naive" as per Wood (2006)
#' @return An object of the \code{\link{gamMD}} class. 
#' If both the mean and the dispersion function are estimated the object will contain the following elements: 
#' \code{converged}, \code{convVec}, \code{data}, \code{fitG}, \code{fitM}, \code{iter}, \code{relE}.
#' \code{fitM} contains information on the mean estimation procedure;  
#' \code{fitG} contains information on the dispersion estimation procedure and has the same elements as \code{fitM}. 
#' 
#' If only the mean function is estimated all the values of the \code{fitM} object are given in \code{gamMD} object.  
#' 
#' The elements of \code{fitM} and \code{fitG} are: 
#' 
#' \item{coefficients}{the estimated regression coefficients}
#' \item{fitted.values}{the fitted values of the model}
#' \item{desMat}{the full design matrix of the model}
#' \item{dims}{the dimension of the design matrix associated to each component}
#' \item{family}{family used for fitting the model - set to \code{Gamma} for \code{fitG}}
#' \item{linear.predictors}{estimated linear predictors for the model}
#' \item{deviance}{deviance residuals - in \code{fitM} these are be scaled by the estimated dispersion function, if present}
#' \item{residuals}{Pearson residuals}
#' \item{s.resid}{standardised Pearson residuals}
#' \item{y, yd}{response variable used in, respectively, the mean and diseprsion estimation procedure}
#' \item{converged}{logic indicator on whether the last inner iteration has converged}
#' \item{sm.p}{smoothing parameters used in the estimation procedure - either fixed or chosen}
#' \item{GCV}{Generalised Cross Validation value}
#' \item{AIC}{Akaike Information criterion value}
#' \item{yt}{the variable response variable transformed according to the link function and then centered}
#' \item{df}{overall equivalent sdegrees of freedom for the fit}
#' \item{vecdf}{a vector given the giving of freedom used by each covariate in the model}
#' \item{cov.coef}{covariance matrix of the coefficients}
#' \item{formula}{formula used}
#' \item{dimsP}{size of the parametric part of the model}
#' 
#' @references Gijbels, Prosdocimi and Claeskens (2010), Nonparametric estimation of mean and dispersion functions in extended generalized linear models, Test, 19(3), 560-608, doi:10.1007/s11749-010-0187-1,
#' @references Gijbels and Prosdocimi (2011), Smooth estimation of mean and dispersion function in extended generalized additive models with application to Italian induced abortion data, Journal of Applied Statistics, 38(11), 2391-2411, doi:10.1080/02664763.2010.550039
#' @seealso  \code{\link{DoubleGamControl}}, \code{\link{DoubleRobGam}}, \code{\link{gamMD}}
#' @author Ilaria Prosdocimi (ilapro@@ceh.ac.uk)
#' @export
DoubleGam <- function (formulaM,formulaG=NULL,family="gaussian", data, startM = NULL, startG = NULL, method="quasi",selection="GCV", weights=NULL, control=DoubleGamControl(), trace = FALSE,scale=1,conf.type = "bayesian", ...) {
    if(method!="quasi" & method!="pseudo") warning('the given method argument was not recognized: use -quasi- instead')
    if(selection!="none" & selection!="GCV" & selection!="AIC"){
    print('how should I select the smoothing parameters? default is GCV')
    selection<-"GCV"
    }
    isgammaEst<-(!is.null(formulaG))
    missdat<-missing(data)
    freadM<-read.form(formulaM,data=data)#     ; print((freadM$dimsP))
       sm.pM<-freadM$sm.p
       XM<-freadM$dataMAT
       interceptM<-freadM$int
       PstrM<-freadM$P
       y<-freadM$y
       dimsSM<-if(!length(freadM$smooths)) 0 else sapply(freadM$smooths,ncol)
       dimsTM<-c(ncol(PstrM)-sum(dimsSM),dimsSM)
    tt<-(trace & !isgammaEst)
    fitM<-gamMeanIntern(sm.p=sm.pM,y=y,XM=XM,Pstr=PstrM,dimsS=dimsSM,dimsT=dimsTM,start=startM,family=family,trace=tt,acc=control$tol,
       maxit=control$maxitG,selection=selection,scale=scale,minlambda=min(control$lambdaM),maxlambda=max(control$lambdaM),conf.type=conf.type)
    if(length(scale)==1 & family=="gaussian"){
       scale<-sum((y-fitM$fitted)^2)/(length(y)-fitM$df)
       fitM<-gamMeanIntern(sm.p=sm.pM,y=y,XM=XM,Pstr=PstrM,dimsS=dimsSM,dimsT=dimsTM,start=startM,family=family,trace=FALSE,acc=control$tol,
          maxit=control$maxitG,selection=selection,scale=scale,minlambda=min(control$lambdaM),maxlambda=max(control$lambdaM),conf.type=conf.type)
    }
    fitM$formula<-formulaM;fitM$family<-family;fitM$dimsP<-freadM$dimsP
    res<-fitM
    if(isgammaEst){
      if(formulaG == "~1"){
          gammaEST<-sum(scale*fitM$deviance)/(length(fitM$deviance)-fitM$df)
          sm.pM<-fitM$sm.p
          fitM<-gamMeanIntern(sm.p=sm.pM,y=y,XM=XM,Pstr=PstrM,dimsS=dimsSM,dimsT=dimsTM,start=startM,family=family,trace=FALSE,acc=control$tol,
                 maxit=control$maxitM,selection=selection,scale=gammaEST,minlambda=min(control$lambdaM),maxlambda=max(control$lambdaM),conf.type=conf.type)
#           gammaEST<-sum(fitM$deviance)/(length((fitM$deviance))-fitM$df)
          convT<-TRUE
          fitM$convT<-convT
          fitM$formula<-formulaM;fitM$family<-family;fitM$dimsP<-freadM$dimsP;names(fitM)
          res<-fitM
          }
   else {
      #### from here there is a pretty complicated thing to read the formula for the dispersion
      #### it can surely be improved
      if(method=="pseudo") dDev<-fitM$residuals^2
      else dDev<-fitM$dispersion*fitM$deviance
      if(length(freadM$seqNA)) { ### when we have missing data, we need to put them back when reading the formula
      dDev<-rep(0,l=length(fitM$deviance)+length(freadM$seqNA))
      dDev[freadM$seqNA]<-NA
      if(method=="pseudo") dDev[-freadM$seqNA]<-fitM$residuals^2
      else  dDev[-freadM$seqNA]<-fitM$dispersion*fitM$deviance # if(method=="quasi")
      }
      Gform<-formula(formulaG)
          if (length(Gform)==2){
                   Gform[3] <- Gform[2];Gform[2] <- formula(dDev~x)[2]
          }
      if (missdat) data <- environment(Gform)
      if (is.environment(data)) assign("dDev",dDev,envir = data)
         else  {data<-data.frame(dDev,data)}
      #### end of preparation to read the formula
      if (missdat) {freadG<-read.form(Gform)} ### actually reading the formula, when the dataset is missing
         else {freadG<-read.form(Gform,data=data)} ### actually reading the formula, if the dataset is given
           sm.pG<-freadG$sm.p
           XG<-freadG$dataMAT
           interceptG<-freadG$int
           PstrG<-freadG$P
           ddev<-freadG$y
           dimsSG<-if(!length(freadG$smooths)) 0 else sapply(freadG$smooths,ncol)
           dimsTG<-c(ncol(PstrG)-sum(dimsSG),dimsSG)
       fitG<-gamDispIntern(sm.p=sm.pG,yd=ddev,XG=XG,Pstr=PstrG,dimsS=dimsSG,dimsT=dimsTG,start=startG,trace=FALSE,acc=control$tol,
                     maxit=control$maxitG,selection=selection,minlambda=min(control$lambdaG),maxlambda=max(control$lambdaG),conf.type=conf.type)
       if(trace) cat("The outer iteration begins","\n")
       for(nit in 1:control$maxitOUT) {
          startM<-fitM$coef
          startG<-fitG$coef
          meanEST<-fitM$fitted.values
          gammaEST<-fitG$fitted.values
          sM<- startM
          sG<- startG
          if(selection != "none"){
             sm.pM<-fitM$sm.p
             sm.pG<-fitG$sm.p
             }
          fitM<-gamMeanIntern(sm.p=sm.pM,y=y,XM=XM,Pstr=PstrM,dimsS=dimsSM,dimsT=dimsTM,start=sM,family=family,trace=FALSE,acc=control$tol,
                 maxit=control$maxitM,selection=selection,scale=gammaEST,minlambda=min(control$lambdaM),maxlambda=max(control$lambdaM),conf.type=conf.type)
          if(method=="pseudo") resp<-fitM$residuals^2
          else  resp<-fitM$dispersion*fitM$deviance
          fitG<-gamDispIntern(sm.p=sm.pG,yd=resp,XG=XG,Pstr=PstrG,dimsS=dimsSG,dimsT=dimsTG,start=sG,trace=FALSE,acc=control$tol,
                        maxit=control$maxitG,selection=selection,minlambda=min(control$lambdaG),maxlambda=max(control$lambdaG),conf.type=conf.type)
#           relE <- c(max(abs((fitM$desMat%*%fitM$coef)/(fitM$desMat%*%startM)),
#                          abs((fitG$desMat%*%fitG$coef)/(fitG$desMat%*%startG))))
	  relE <- c(sqrt(sum((fitM$fitted-meanEST)^2)/max(1e-20, sum(meanEST^2))),
                         sqrt(sum((fitG$fitted-gammaEST)^2)/max(1e-20, sum(gammaEST^2))))
          if(trace) {
	  cat("iteration", nit,"\n")
	  cat("the chosen smoothing parameters are","\n")
	  cat(sm.pM,"for the mean -",fitM$df,"degrees of freedom","\n")
	  cat(sm.pG,"for the dispersion -",fitG$df,"degrees of freedom","\n")
# 	  cat("the relative change in the coefficients is",relE,"\n","\n")
	  cat("the relative change in the estimates is",relE,"\n","\n")
	  }
      convT <- all(relE <= control$acc)
      if(convT) break
      }
   if(!convT) warning('no convergence for the general algorithm')
   convVec<-c(fitM$conv,fitG$conv,convT);names(convVec)<-c("mean conv","dispersion conv","two-steps conv")
   fitM$formula<-formulaM;fitG$formula<-formulaG
   fitM$dimsP<-freadM$dimsP;fitG$dimsP<-freadG$dimsP
   ### really ugly, but some things should not be in output....
   fitM<- fitM[seq(1:length(fitM))[-which(names(fitM) == "iter")]]
   fitM<- fitM[seq(1:length(fitM))[-which(names(fitM) == "dispersion")]]
   fitG<- fitG[seq(1:length(fitG))[-which(names(fitG) == "dispersion")]]
   res<-list(fitM=fitM,fitG=fitG,convVec=convVec,converged=convT,iter=nit,relE=relE)
   if(missdat) rm(dDev,envir = data)
     }
     }
res$seqNA<-freadM$seqNA
if(missdat) res$data<-1
if(!missdat) res$data<-data
class(res)<-"gamMD"
res
}

