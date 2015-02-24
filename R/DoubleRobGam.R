# 
# 
# ######## The DoubleRobGam function
#     ####### you need the general functions.R and the robFunctions.R files as well
# 
# #### In this file, first you find all the functions to estimate the mean function
#       ######### then all the function for the dispersion estimation
#       ######### the Double part comes at the end
# ######## All functions written by Ilaria Prosdocimi (ilapro@ceh.ac.uk)
# 
# 
# 
# 
# #### mean estimation
# 
# 
# gamRobMeanGivenSP<-function(sm.p,y,XM,Pstr,dimsS,dimsT,start=NULL,family,trace=FALSE,acc=10^-6,scale=1,tcc=1.5,weights.on.x="none",maxit=30,weights=NULL){
#       nobs <- NROW(y);  ncoef <- ncol(XM)
#       variance <- family$variance
#       linkinv <- family$linkinv
#       if (!is.function(variance) || !is.function(linkinv)) stop("illegal 'family' argument")
#       mu.eta <- family$mu.eta
#       if (is.null(valideta <- family$valideta)) valideta <- function(eta) TRUE
#       if (is.null(validmu <- family$validmu))  validmu <-  function(mu) TRUE
#       if (is.null(weights))  weights<-rep(1,l=nobs)
#       ## note that 'weights' are used and set by binomial()$initialize !
#       etastart<-mustart<-NULL
#       residPS <- sV <- dmu.deta <- residP <- NULL ## to pass R CMD CHECK
#       eval(family$initialize) ## --> n, mustart, y and weights (=ni)
#       if(is.null(start)) start<-c(mean(family$linkfun(mustart)),rep(0.001,l=(ncoef-1))) 
#       w.x <- if(ncoef) { ##### straight from robustbase
#          switch(weights.on.x,
#                "none" = rep.int(1, nobs),
#                "hat" = wts_HiiDist(X = XM),
#                "robCov" = wts_RobDist(X=XM, intercept = as.logical(all(XM[,1] < 1.0001 & XM[,1]>0.9999)), covFun = MASS::cov.rob),
#                                         # ARu said  'method="mcd" was worse'
#                "covMcd" = wts_RobDist(X=XM, intercept = as.logical(all(XM[,1] < 1.0001 & XM[,1]>0.9999)), covFun = covMcd),
#                stop("Weighting method", sQuote(weights.on.x)," is not implemented"))
#     }
#     else rep.int(1,nobs) ## ncoef == 0 
#     ### Initializations
#     stopifnot(maxit >= 1, tcc >= 0)
# #       eval(family$initialize) ## --> n, mustart, y and weights (=ni)
#     ni <- as.vector(weights) # dropping attributes for computation
#     betaOld <- beta <- as.vector(start) # as.v*(): dropping attributes
#     etaOld<- eta <- as.vector(XM %*% beta)
#     muOld<- mu <- linkinv(eta) # mu estimates pi (in [0,1]) at the binomial model
#     if (!(validmu(mu) && valideta(eta)))
# 	stop("Cannot find valid starting values: You need help")
#     switch(family$family,  ##### straight from robustbase
# 	   "binomial" = {
# 	       Epsi.init <- EpsiBin.init
# 	       Epsi <- EpsiBin
# 	       EpsiS <- EpsiSBin
# 	       Epsi2 <- Epsi2Bin
#                phiEst <- phiEst.cl <- expression({1})
#                Huber2<-Huber2.Bin
#                DevFun<-DevBinom
# 	   },
# 	   "poisson" = {
#                Epsi.init <- EpsiPois.init
#                Epsi <- EpsiPois
#                EpsiS <- EpsiSPois
#                Epsi2 <- Epsi2Pois
#                phiEst <- phiEst.cl <- expression({1})
#                Huber2<-Huber2.Pois
#                DevFun<-DevPois
#     },
#     "gaussian" = {
#                Epsi.init <- EpsiNorm.init
#                Epsi <- EpsiNorm
#                EpsiS <- EpsiSNorm
#                Epsi2 <- Epsi2Norm
#                phiEst <- phiEst.cl <- phiNormEst.cl
#                Huber2<-Huber2.Norm
#                DevFun<-DevNorm
#     },
#       ## else
#       stop(gettextf("family '%s' not yet implemented", family$family))
#       )
#     comp.V.resid <- expression({
#         Vmu <- variance(mu)
#         if (any(is.na(Vmu)))  stop("NAs in V(mu)")
#         if (any(Vmu == 0))    stop("0s in V(mu)")
#         sVF <- sqrt(Vmu)   # square root of variance function
#         residP <- (y - mu)* sni/sVF  # Pearson residuals
#     })
#     comp.scaling <- expression({
#       sV <- sVF * sqrt(phi)
#       residPS <- residP/sqrt(phi) # scaled Pearson residuals
#     })
#     comp.Epsi.init <- expression({
#         ## d mu / d eta :
#         dmu.deta <- mu.eta(eta)
#         if (any(is.na(dmu.deta))) stop("NAs in d(mu)/d(eta)")
#         ## "Epsi init" :
#         H <- floor(mu*ni - tcc* sni*sV)
#         K <- floor(mu*ni + tcc* sni*sV)
#         eval(Epsi.init)
#    })
# #### This build the penalty matrix for the given smoothing parameters
#      if(dimsT[1] != 0) supportMat<-matrix(0,nrow=nrow(Pstr),ncol=dimsT[1])
#          else supportMat<-NULL
#      if(length(sm.p) >= 1) for(i in 1:length(sm.p))
#      supportMat<-cbind(supportMat,matNum((sm.p[i]/100),dim1=nrow(Pstr),dim2=dimsS[i]))
#      P<-supportMat*Pstr
# ### Iterations for the mean
#      if(trace && ncoef) {
#          cat("The iteration begins: \n")
#     }
#     sni <- sqrt(ni)
# #     if(length(scale)==1) phi<-eval(phiEst.cl)
# #     if(length(scale)!=1) phi <- scale ### eval(phiEst.cl) shall I do this? no I guess
#     phi <- scale
#     conv <- FALSE
#     if(length(sm.p)==0) inversion<-function(...) {solve(...)} ##### when the model is only parametric use solve()
#     if(length(sm.p)!=0) inversion<-function(...) {ginv(...)}  ##### otherwise we have results different from glm()
#     if(ncoef) for (nit in 1:maxit) {
#         sni <- sqrt(ni)
#         eval(comp.V.resid) #-> (Vmu, sVF, residP)
#         eval(comp.scaling) #-> (sV, residPS)
#         eval(comp.Epsi.init)
#         cpsi <- pmax2(-tcc, pmin(residPS,tcc)) - eval(Epsi)
# # 	This is what you have in robustbase:
# # 	EEq <- colMeans(cpsi * w.x * sni/sV * dmu.deta * X)
# #       this is what is in robustbase
# # 	Dtheta <- solve(crossprod(X, DiagB*X)/nobs, EEq)
# # 	print("inversion is");print(solve((t(X)%*%diag(DiagB)%*%X)/nobs))
# #	Dtheta <- solve((t(X)%*%diag(DiagB)%*%X)/nobs)%*%EEq
# #  	theta <- thetaOld + Dtheta
# # 	eta <- as.vector(X %*% theta) + offset
# # 	##
# # 	## Solve  1/n (t(XM) %*% B %*% XM) %*% delta.coef	  = EEq
# # cat('here it is',eval(Epsi2),'\n')
# 	DiagB <- eval(EpsiS) /(sni*sV) * w.x * (ni*dmu.deta)^2
# 	### This is different from robustbase! 
#         BB<-diag(DiagB)/nobs
#         beta<-inversion(t(XM)%*%BB%*%XM+P)%*%t(XM)%*%BB%*%(eta+diag(1/DiagB)%*%as.vector(cpsi*w.x*sni/sV*dmu.deta))
#         Dbeta<-betaOld-beta  ### the difference. I use the same stopping rule (code) as in robustbase
# 	if (any(!is.finite(Dbeta))) {
# 	    warning("Non-finite coefficients at iteration ", nit)
# 	    break
# 	}
# 	eta <- as.vector(XM %*% beta);mu <- linkinv(eta)
#         ## estimation of the dispersion parameter no longer here
#         eval(comp.V.resid)
# #         phi <- scale##*eval(phiEst)
# 	## Check convergence: relative error < tolerance
# 	relE0 <- sqrt(sum(Dbeta^2)/max(1e-20, sum(betaOld^2)))
# 	relE <- sqrt(sum((mu-muOld)^2)/max(1e-20, sum(muOld^2)))
# #         cat('coef is', relE0,'   func is' ,relE,'   and coef < func ?  ', relE0<relE , "\n")
# 	conv <- relE <= acc
#         if(trace) cat("at iteration ",nit,"the relative change in the function estimate is ",relE,'\n')
# 	if(conv) break
# 	betaOld <- beta
# 	muOld<-mu ;etaOld<-eta
#     } ## end of iteration
#     else { ## ncoef == 0
# 	conv <- TRUE
# 	nit <- 0
#     }
#     eps <- 10 * .Machine$double.eps
#     switch(family$family,
# 	   "binomial" = {
# 	       if (any(mu/weights > 1 - eps) || any(mu/weights < eps))
# 		   warning("fitted probabilities numerically 0 or 1 occurred")
# 	   },
# 	   "poisson" = {
# 	       if (any(mu < eps))
# 		   warning("fitted rates numerically 0 occurred")
# 	   })
#     eval(comp.V.resid) #-> (Vmu, sVF, residP)
#     eval(comp.scaling) #-> (sV, residPS)
#     ## Estimated asymptotic covariance of the robust estimator
#          ######### shouldn't be here anymore
#     if(ncoef) {
#         eval(comp.Epsi.init)
#         DiagB <- eval(EpsiS) / (sni*sV)* w.x * (ni*dmu.deta)^2
#         if(any(n0 <- ni == 0)) DiagB[n0] <- 0 # instead of NaN
#     } else { ## ncoef == 0
#         DiagB <- matrix(, 0,0)
#     }
#     beta<-as.vector(beta)
#     w.r <- pmin(1, tcc/abs(residPS))
#     Dev<-DevFun(y=y,mu=mu,wt=ni)/scale ### scaled deviances, so that when I estimate the dispersion they get individual weights
#     df<-sum(diag(inversion(t(XM)%*%BB%*%XM+P)%*%t(XM)%*%BB%*%XM))
#     RGCV<- sum(pmax2(-tcc, pmin(Dev,tcc)))/((1-df/nobs)^2)
#     RAIC<- sum(pmax2(-tcc, pmin(Dev,tcc)))+2*df
#     GCV<- sum(Dev)/((1-df/nobs)^2)
#     AIC<- sum(Dev)+2*df
#     mid<-family$linkfun(y)
#     mid<-mean(mid[is.finite(mid) & !is.na(mid)])
#     yt<-family$linkfun(y)-mid
#     eval(comp.Epsi.init);Epsi2val <- eval(Epsi2)
#     uu<-try(uniroot(Huber2,ns.resid=residP,eta=eta,Vmu=mu,tcc=tcc,Epsi2val=Epsi2val,interval=range(residP^2)),silent=TRUE)
# if(length(uu) == 1) uu<-try(uniroot(Huber2,ns.resid=residP,eta=eta,Vmu=mu,tcc=tcc,Epsi2val=Epsi2val,interval=c(0.1,5)),silent=TRUE)
# if(length(uu) == 1) uu<-try(uniroot(Huber2,ns.resid=residP,eta=eta,Vmu=mu,tcc=tcc,Epsi2val=Epsi2val,interval=c(0.51,2.5)),silent=TRUE)
# if(length(uu) == 1) uu<-try(uniroot(Huber2,ns.resid=residP,eta=eta,Vmu=mu,tcc=tcc,Epsi2val=Epsi2val,interval=c(0.0051,8000012.5)),silent=TRUE)
# if(length(uu) == 1) {uu<-list(root=sum(scale*Dev)/(nobs-df))}
# ######## I know that the robust dispersion parameter estimation is not always working, I need to fix this
# #  warning('dispersion parameter estimation is not working (SHOULD THIS MESSAGE APPEAR?)') ### if you uncomment this it prints it out A LOT
#     if(ncoef) {
# 	eval(comp.Epsi.init)
# 	alpha <- colMeans(eval(Epsi) * w.x * sni/sV * dmu.deta * XM)
# # 	print(alpha)# ;print(w.x);print(sni);print(sV);print(mu)
# 	DiagA <- eval(Epsi2) / (ni*sV^2)* w.x^2* (ni*dmu.deta)^2
# # 	print((P%*%beta%*%t(beta)%*%P)[1:5,1:5])
# 	matQ  <- crossprod(XM, DiagA*XM)/nobs  - tcrossprod(alpha, alpha)# + (P%*%beta%*%t(beta)%*%P)
# 	DiagB <- eval(EpsiS) / (sni*sV)* w.x * (ni*dmu.deta)^2
#         if(any(n0 <- ni == 0)) DiagB[n0] <- 0 # instead of NaN
# 	matM <- (crossprod(XM, DiagB*XM)/nobs+P)
# 	matMinv <- inversion(matM)
# 	asCov <-  matMinv %*% matQ %*% matMinv / nobs
#     } else { ## ncoef == 0
# 	matM <- matQ <- asCov <- matrix(, 0,0)
#     }
#     vecdfM<-NULL
#     if(dimsT[2]!= 0){
#        for(i in 2:length(dimsT)){
#        mm<-XM[,(1+sum(dimsT[0:(i-1)])):(sum(dimsT[0:i]))]
#        pm<-P[(1+sum(dimsT[0:(i-1)])):(sum(dimsT[0:i])),(1+sum(dimsT[0:(i-1)])):(sum(dimsT[0:i]))]
#        vecdfM<-c(vecdfM,sum(diag(mm%*%inversion(t(mm)%*%BB%*%mm+pm)%*%t(mm)%*%BB)))
#        }
#     }
#     list(coefficients=as.vector(beta),fitted.values=mu,desMat=XM,dims=dimsT,dispersion=phi,family=family,
#         linear.predictors=eta,deviance=Dev,residuals=residP,s.resid=residPS,iter=nit,y=y,converged=conv,weights=weights,
#         sm.p=sm.p,GCV=GCV,AIC=AIC,RGCV=RGCV,RAIC=RAIC,w.x=w.x,estRphi=uu$root,tcc=tcc,yt=yt,df=df,vecdf=vecdfM,
#         Mmat=matM,Qmat=matQ,cov.coef=asCov)
# }
# 
# 
# 
# 
# 
# 
# gamRobMeanIC<-function(lsm.p,y,XM,Pstr,dimsS,dimsT,start = NULL,family,acc = 10^-6, maxit = 30, crit="RGCV",scale=1,tcc=1.345,weights.on.x="none",weights=NULL){
# ### computes the information criterion of the mean for the given smoothing parameters
#     if(crit=="GCV" | crit=="AIC" | crit=="RGCV" | crit=="RAIC") {
#         gg<-gamRobMeanGivenSP(sm.p=exp(lsm.p),y=y,XM=XM,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,start=start,
#             family=family,trace=FALSE,acc=acc,maxit=maxit,scale=scale,tcc=tcc,weights.on.x=weights.on.x,weights=weights)
#     if(crit=="GCV") res<-gg$GCV
#     if(crit=="AIC") res<-gg$AIC 
#     if(crit=="RGCV") res<-gg$RGCV
#     if(crit=="RAIC") res<-gg$RAIC
#     }
# res
# }
# 
# 
# ##### this is a function which would serve to only estimate the mean. 
# ##### you can do this using DoubleGam and then formulaG=NULL and it is out of my interest
# # gamRobMean <- function(formula, weights = NULL, start = NULL, family= "gaussian", weights.on.x = "none", acc = 10^-6 , 
# #      maxit = 50, tcc = 1.345, trace = FALSE, scale=1,selection="none",data){
# #        fread<-read.form(formula,data=data)
# #        XM<-fread$dataMAT
# #        Pstr<-fread$P
# #        y<-fread$y
# #        sm.ps<-fread$sm.p
# #        dimsS<-if(!length(fread$smooths)) 0 else sapply(fread$smooths,ncol)
# #        dimsT<-c(ncol(Pstr)-sum(dimsS),dimsS)
# #        if (is.character(family))
# #          family <- get(family, mode = "function", envir = parent.frame())
# #        if (is.function(family))
# #          family <- family()
# #        if (is.null(weights))
# #          weights <- rep.int(1, NROW(y))
# #        else if(any(weights <= 0))
# #          stop("All weights must be positive")
# #        if(selection!="none" & selection!="GCV" & selection!="AIC" & selection!="RGCV" & selection!="RAIC"){
# #           cat("how should I select the mean smoothing parameters? default is RGCV \n")
# #           selection<-"RGCV"
# #        }
# #        if (selection!="none") {
# # #          minim<-nlm(f=gamRobMeanIC,p=log(sm.ps),y=y,XM=XM,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,start=NULL,family=family,
# # #                  acc=acc,scale=scale,tcc=tcc,weights.on.x=weights.on.x,weights=weights,maxit=maxit,crit=selection)
# # #         sm.ps<-exp(minim$est)
# #          sm.ps<-exp(optim(p=log(sm.ps),gamRobMeanIC,yy=y,XM=XM,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,start=NULL,family=family,
# #                  acc=acc,scale=scale,tcc=tcc,weights.on.x=weights.on.x,weights=weights,maxit=maxit,crit=selection,
# #                  method="L-BFGS-B",lower=log(0.001),upper=log(5000))$par)
# # #
# # # #         sm.ps<-exp(optim(p=log(sm.ps),gamRobMeanIC,yy=y,XM=XM,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,start=NULL,family=family,
# # # #                  acc=acc,scale=scale,tcc=tcc,weights.on.x=weights.on.x,weights=weights,maxit=maxit,crit=selection,
# # # #                  method="L-BFGS-B",lower=log(0.001),upper=log(5000))$par)
# #         }
# #      fit<-gamRobMeanGivenSP(sm.p=sm.ps,y=y,XM=XM,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,start=NULL,family=family,
# #                trace=trace,acc=acc,scale=scale,tcc=tcc,weights.on.x=weights.on.x,weights=weights,maxit=maxit,compRdev=compRdev)
# #      fit$seqNA<-fread$seqNA
# #      fit
# # }
# 
# 
# 
# 
# gamRobMeanIntern<-function(sm.p,y,XM,Pstr,dimsS,dimsT,start,family,trace=FALSE,acc=10^-6,maxit=30,tcc,weights.on.x,weights,
#                                     selection,scale=1,minlambda=0.001,maxlambda=5000){
#     sm.pinF<-sm.p
#     if(selection=="none") sm.ps<-sm.pinF
# #     else sm.ps<-exp(nlm(p=log(sm.pinF),gamRobMeanIC,
# #          y=y,XM=XM,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,start=start,family=family,acc=acc,maxit=maxit,
# #          scale=scale,tcc=tcc,weights.on.x=weights.on.x,weights=weights,crit=selection)$est)
# #
#      else sm.ps<-exp(optim(par=log(sm.pinF),fn=gamRobMeanIC,y=y,XM=XM,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,start=start,family=family,
#                  acc=acc,scale=scale,tcc=tcc,weights.on.x=weights.on.x,weights=weights,maxit=maxit,crit=selection,
#                  method="L-BFGS-B",lower=log(minlambda),upper=log(maxlambda))$par)
#    if(trace==TRUE) cat('the smoothing parameters are', sm.ps, '\n')
#    fit<-gamRobMeanGivenSP(sm.p=sm.ps,y=y,XM=XM,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,start=start,family=family,
#                trace=trace,acc=acc,scale=scale,tcc=tcc,weights.on.x=weights.on.x,weights=weights,maxit=maxit)
#    fit
# }
# 
# 
# 
# 
# 
# ########## dispersion estimation
# 
# 
# #### this functions takes a lot from glmrob for a gamma family
# 
# gamRobDispGivenSP<-function(sm.p,yd,XG,Pstr,dimsS,dimsT,start=NULL,trace=FALSE,acc=10^-6,tcc=1.345,weights.on.x="none",maxit=30,weights=NULL,link="log"){
# ### fits the model for the given smoothing parameters
#       scale<-phi<-2
#       nobs <- NROW(yd);ncoef <- ncol(XG)
#       variance <- function(gamma)  gamma^2
#       linkinv  <- function(xi) pmax(exp(xi),.Machine$double.eps)
#       gamma.xi <- function(xi) pmax(exp(xi),.Machine$double.eps)
#       validxi <-  function(xi) TRUE
#       validgamma <- function(gamma)  all(gamma > 0)
#       w.x <- if(ncoef) {
# ##### straight from robustbase
#          switch(weights.on.x,
#                "none" = rep.int(1, nobs),
#                "hat" = wts_HiiDist(X = XG),
#                "robCov" = wts_RobDist(X =XG, intercept =  as.logical(all(XG[,1] < 1.0001 & XG[,1]>0.9999)), covFun = MASS::cov.rob),
#                                         # ARu said  'method="mcd" was worse'
#                "covMcd" = wts_RobDist(X =XG, intercept = as.logical(all(XG[,1] < 1.0001 & XG[,1]>0.9999)), covFun = covMcd),
#                stop("Weighting method",sQuote(weights.on.x)," is not implemented"))
#     }
#     else rep.int(1,nobs) ## ncoef == 0 
#   ### Initializations
#   residPS <- sV <- dgamma.dxi <- residP <- NULL ## to pass R CMD CHECK
#   stopifnot(maxit >= 1, tcc >= 0)
#     initf<-expression({
#     ### this is Gamma(link="log")$initialize , with some necessary modifications
#         if (any(yd <= 0)) stop("non-positive values not allowed for the gamma family")
#         n <- rep.int(1, nobs) ;  gammastart <- yd
#     })
#     eval(initf) ## --> n, gammastart, yd and weights (=ni)
#     ni <- as.vector(weights)
#     if(is.null(start)) start <- c(mean(log(yd[yd != 0]),na.rm=TRUE),rep(0.0001,l=ncoef-1))
#     # rep(mean(log(yd[yd != 0]),na.rm=TRUE),l=ncoef) #### default starting value
#     deltaOld <- delta <- as.vector(start)
#     xiOld<- xi <- as.vector(XG %*% delta)
#     gammaOld<-gamma <- linkinv(xi)
#     if (!(validgamma(gamma) && validxi(xi)))
#        stop("Cannot find valid starting values: You need help")
#     Epsi.init <- expression({
#          nu <- 1/phi      ## form parameter nu
#          snu <- 1/sqrt(phi) ## sqrt (nu)
#          pPtc <- pgamma(snu+c(-tcc,tcc),shape=nu,rate=snu)
#          pMtc <- pPtc[1]
#          pPtc <- pPtc[2]
#          aux2 <- tcc*snu
#          GLtcc <- Gmn(-tcc,nu)
#          GUtcc <- Gmn( tcc,nu)
#     })
#     Epsi <- expression( tcc*(1-pPtc-pMtc) + GLtcc - GUtcc )
#     EpsiS <- expression( ((GLtcc - GUtcc) + snu*(pPtc-pMtc))/gamma )
#     Epsi2 <- expression({
#           (tcc^2*(pMtc+1-pPtc)+ (pPtc-pMtc) +
#           (GLtcc*(1-aux2) - GUtcc*(1+aux2))/snu )
#     })
#     comp.V.resid <- expression({
#         Vgamma <- variance(gamma)
#         if (any(is.na(Vgamma)))  stop("NAs in V(gamma)")
#         if (any(Vgamma == 0))    stop("0s in V(gamma)")
#         sVF <- sqrt(Vgamma)   # square root of variance function
#         residP <- (yd - gamma)* sni/sVF  # Pearson residuals
#     })
#     comp.scaling <- expression({
#       sV <- sVF*sqrt(phi);residPS<-residP/sqrt(phi) # scaled Pearson residuals
#     })
#     comp.Epsi.init <- expression({
#         ## d gamma / d xi :
#         dgamma.dxi <- gamma.xi(xi)
#         if (any(is.na(dgamma.dxi))) stop("NAs in d(gamma)/d(xi)")
#         ## "Epsi init" :
#         H <- floor(gamma*ni - tcc* sni*sV)
#         K <- floor(gamma*ni + tcc* sni*sV)
#         eval(Epsi.init)
#     })
#     #### this build the penalty matrix for smoothing
#     if(dimsT[1] != 0) supportMat<-matrix(0,nrow=nrow(Pstr),ncol=dimsT[1])
#         else supportMat<-NULL
#     matNum<-function(num,dim1,dim2) matrix(num,ncol=dim2,nrow=dim1)
#     if(length(sm.p) >= 1) for(i in 1:length(sm.p))
#     supportMat<-cbind(supportMat,matNum((sm.p[i]/100),dim1=nrow(Pstr),dim2=dimsS[i]))
#     P<-supportMat*Pstr
#     ### Iterations
#     sni <- sqrt(ni) ### sqrt of weigths 
#     eval(comp.V.resid) #-> (Vgamma, sVF, residP)
#     phi <- scale ### that's because we work with a chisq(and scale = 2)
#     conv <- FALSE
#     if(length(sm.p)==0) inversion<-function(...) {solve(...)} ##### when the model is only parametric use solve()
#     if(length(sm.p)!=0) inversion<-function(...) {ginv(...)}  ##### otherwise we have results different from glm()
#     if(ncoef) for (nit in 1:maxit) {
#         if(trace)  cat("The iteration begins: \n")
#         eval(comp.scaling) #-> (sV, residPS)
#         eval(comp.Epsi.init)
#         cpsi <- pmax2(-tcc, pmin(residPS,tcc)) - eval(Epsi)
#         EEq <- colMeans(cpsi * w.x * sni/sV * dgamma.dxi * XG)
#         DiagB <- eval(EpsiS) /(sni*sV) * w.x * (ni*dgamma.dxi)^2
#         if(any(n0 <- ni == 0)) DiagB[n0] <- 0 # instead of NaN
#         BB<-diag(DiagB)/nobs
#         delta<-inversion(t(XG)%*%BB%*%XG+P)%*%t(XG)%*%BB%*%(xi+diag(1/DiagB)%*%as.vector(cpsi*w.x*sni/sV*dgamma.dxi))
#         deltaOld <- delta
#         Ddelta<-deltaOld-delta ### change in the estimate. if I use the same stopping rule (code) as in robustbase 
#         if (any(!is.finite(Ddelta))) {
#         warning("Non-finite coefficients at iteration ", nit)
#             break
#         }
#         xi <- as.vector(XG %*% delta);gamma <- linkinv(xi)
#         ## dispersion parameter is equal to 2
#         eval(comp.V.resid); phi <- 2
#         ## Check convergence: relative error < tolerance
#         relE0 <- sqrt(sum(Ddelta^2)/max(1e-20, sum(deltaOld^2)))
#         relE <- sqrt(sum((gamma-gammaOld)^2)/max(1e-20, sum(gammaOld^2)))
#         conv <- relE <= acc
#         if(trace) cat("at iteration ",nit,"the relative change in the function estimate is ",relE,'\n')
#         if(conv) break
#         deltaOld <- delta;gammaOld<-gamma;xiOld<-xi
#         }
#         else { ## ncoef == 0
#            conv <- TRUE
#            nit <- 0
#         }
#     eps <- 10 * .Machine$double.eps
#     eval(comp.V.resid) #-> (Vmu, sVF, residP)
#     eval(comp.scaling) #-> (sV, residPS)
#     ## Estimated asymptotic covariance of the robust estimator
#     ##### should go
#     if(ncoef) {
#         eval(comp.Epsi.init)
#         DiagB <- eval(EpsiS) / (sni*sV)* w.x * (ni*dgamma.dxi)^2
#     } else { ## ncoef == 0
#         matM <- matQ <- asCov <- matrix(, 0,0)
#     }
#     w.r <- pmin(1, tcc/abs(residPS))
#     names(gamma) <- names(xi) <- names(delta) # re-add after computation
#     Dev<-Gamma.resid(y=yd,mu=gamma,wt=ni)/2
#     df<- sum(diag(inversion(t(XG)%*%BB%*%XG+P)%*%t(XG)%*%BB%*%XG))
#     GCV<- sum(Dev)/((1-df/nobs)^2)
#     AIC<- sum(Dev)+2*df
#     RGCV<- nobs*sum(pmax2(-tcc, pmin(Dev,tcc)))/((1-df/nobs)^2)
#     RAIC<- sum(pmax2(-tcc, pmin(Dev,tcc)))+2*df
#     mid<-log(yd)
#     mid<-mean(mid[is.finite(mid) & !is.na(mid)])
#     yt<-log(yd)-mid
#     delta<-as.vector(delta)
#     vecdfG<-NULL
#     if(ncoef) {
# 	eval(comp.Epsi.init)
# 	alpha <- colMeans(eval(Epsi) * w.x * sni/sV * dgamma.dxi * XG)
# # 	print(alpha)# ;print(w.x);print(sni);print(sV)
# 	DiagA <- eval(Epsi2) / (ni*sV^2)* w.x^2* (ni*dgamma.dxi)^2
# # 	print(eval(Epsi2));print(gamma[1:10])
# 	matQ  <- crossprod(XG, DiagA*XG)/nobs - tcrossprod(alpha, alpha)
# 	DiagB <- eval(EpsiS) / (sni*sV)* w.x * (ni*dgamma.dxi)^2
#         if(any(n0 <- ni == 0)) DiagB[n0] <- 0 # instead of NaN
# 	matM <- (crossprod(XG, DiagB*XG)/nobs+P)
# 	matMinv <- inversion(matM)
# 	asCov <-  matMinv %*% matQ %*% matMinv / nobs
#     } else { ## ncoef == 0
# 	matM <- matQ <- asCov <- matrix(, 0,0)
#     }
#    if(dimsT[2]!= 0){
#    for(i in 2:length(dimsT)){
#         mm<-XG[,(1+sum(dimsT[0:(i-1)])):(sum(dimsT[0:i]))]
#         pm<-P[(1+sum(dimsT[0:(i-1)])):(sum(dimsT[0:i])),(1+sum(dimsT[0:(i-1)])):(sum(dimsT[0:i]))]
#         vecdfG<-c(vecdfG,sum(diag(mm%*%inversion(t(mm)%*%BB%*%mm+pm)%*%t(mm)%*%BB)))
#         }
#     }
#     list(coefficients=as.vector(delta),fitted.values=gamma,desMat=XG,dims=dimsT,dispersion=2,family=Gamma(link="log"),
#             linear.predictors=xi,deviance=Dev,residuals=residP,s.resid=residPS,iter=nit,yd=yd,converged=conv,weights=weights,
#             sm.p=sm.p,GCV=GCV,AIC=AIC,RGCV=RGCV,RAIC=RAIC,w.x=w.x,tcc=tcc,yt=yt,df=df,vecdf=vecdfG,Mmat=matM,Qmat=matQ,cov.coef=asCov,penalty=P)
#             # w.r=w.r,,Pstr=Pstr
# }
# 
# 
# 
# 
# 
# Gamma.resid<- function(y, mu, wt) -2 * wt * (log(ifelse(y == 0, 1, y/mu)) - (y - mu)/mu)
# 
# 
# 
# gamRobDispIC<-function(lsm.p,yd,XG,Pstr,dimsS,dimsT,start=NULL,acc=10^-6,maxit=30, crit="RGCV",tcc=1.345,weights.on.x="none",weights=NULL,link="log"){
# ### computes the information criterion of the mean for the given smoothing parameters
#     gg<-gamRobDispGivenSP(sm.p=exp(lsm.p),yd=yd,XG=XG,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,start=start,
#             trace=FALSE,acc=acc,maxit=maxit,tcc=tcc,weights.on.x=weights.on.x,weights=weights,link=link)
#     if(crit=="GCV") res<-gg$GCV
#     if(crit=="AIC") res<-gg$AIC
#     if(crit=="RGCV") res<-gg$RGCV
#     if(crit=="RAIC") res<-gg$RAIC
# res
# }
# 
# 
# 
# 
# 
# #############  not really needed....
# # gamRobDisp<-
# #     function(formula,weights=NULL,start=NULL,weights.on.x="none",
# #                 acc = 10^-6,maxit=50,tcc=1.345,trace=FALSE,link="log",selection="GCV",data){
# #     fread<-read.form(formula,data=data)
# #       XG<-fread$dataMAT
# #       Pstr<-fread$P
# #       yd<-fread$y
# #       sm.ps<-fread$sm.p
# #       intercept<-fread$int
# #       nobs <- NROW(yd);ncoef<-ncol(XG)
# #       dimsS<-if(!length(fread$smooths)) 0 else sapply(fread$smooths,ncol)
# #       dimsT<-c(ncol(Pstr)-sum(dimsS),dimsS)
# #       if (is.null(weights))
# #         weights <- rep.int(1, NCOL(yd))
# #       else if(any(weights <= 0))
# #         stop("All weights must be positive")
# #      if(selection!="none" & selection!="GCV" & selection!="AIC" & selection!="RGCV" &  selection!="RAIC"){
# #      print('how should I select the smoothing parameters? default is RGCV')
# #      selection<-"RGCV"
# #      }
# # #      if(selection=="none") sm.ps<-sm.ps
# # #       else sm.ps<-exp(nlm(p=log(sm.ps),gamRobDispIC,y=yd,XG=XG,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,start=start,
# # #                      acc=acc,maxit=maxit,crit=selection,tcc=tcc,weights.on.x=weights.on.x,weights=weights,link=link)$est)
# # #       else 
# #      if(selection!="none")
# # #         sm.ps<-exp(nlm(p=log(sm.ps),gamRobDispIC,y=yd,XG=XG,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,start=start,
# # #              acc=acc,maxit=maxit,crit=selection,tcc=tcc,weights.on.x=weights.on.x,weights=weights,link=link)$est)
# #      sm.ps<-exp(optim(p=log(sm.ps),gamRobDispIC,y=yd,XG=XG,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,
# #            start=start,acc=acc,maxit=maxit,crit=selection,tcc=tcc,weights.on.x=weights.on.x,weights=weights,link=link,
# #            method="L-BFGS-B",lower=log(0.05),upper=log(3000))$par)
# #      if(trace==TRUE) cat('the smoothing parameters are', sm.ps, '\n')
# #      fit<-gamRobDispGivenSP(sm.p=sm.ps,yd=yd,XG=XG,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,start=start,trace=FALSE,
# #           acc=acc,tcc=tcc,weights.on.x=weights.on.x,maxit=maxit,weights=weights,link=link,compRdev=compRdev)
# #      fit$seqNA<-fread$seqNA
# #      fit
# # }
# 
# 
# gamRobDispIntern<-function(sm.p,yd,XG,Pstr,dimsS,dimsT,start,trace=FALSE,acc=10^-6,maxit=30,tcc,weights.on.x,weights,
#                                      selection,link="log",minlambda=0.001,maxlambda=5000){
#     sm.pinFd<-sm.p
#     if(selection=="none") sm.ps<-sm.pinFd
# #      else sm.ps<-exp(nlm(p=log(sm.pinF),gamRobDispIC,
# #          yd=yd,XG=XG,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,start=start,acc=acc,maxit=maxit,
# #          tcc=tcc,weights.on.x=weights.on.x,weights=weights,crit=selection,link=link)$est)
#    if(selection!="none")  sm.ps<-exp(optim(par=log(sm.pinFd),fn=gamRobDispIC,y=yd,XG=XG,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,
#            start=start,acc=acc,maxit=maxit,crit=selection,tcc=tcc,weights.on.x=weights.on.x,weights=weights,link=link,
#            method="L-BFGS-B",lower=log(minlambda),upper=log(maxlambda))$par)
#    if(trace==TRUE) cat('the smoothing parameters are',sm.ps,'\n')
#    fit<-gamRobDispGivenSP(sm.p=sm.ps,yd=yd,XG=XG,Pstr=Pstr,dimsS=dimsS,dimsT=dimsT,start=start,
#                acc=acc,tcc=tcc,weights.on.x=weights.on.x,weights=weights,maxit=maxit,link=link)
# fit
# }
# 
# 
# 
# 
# 
# 
# 
# ######## double estimation (robust and smooth)
# 
# 
# 
# 

#' Setting the DoubleRobGam fitting defaults 
#' 
#' This function is used internally in the DoubleRobGam fitting procedure. 
#' The function's parameters control the numerical properties of the Double Robust fitting procedure. 
#' Care should be taken when changing the parameters default values, things can go wrong 
#' (but sometimes fitting doesn't work well with the initial parameters, so you'll need to fiddle around with this)
#' @param maxitM the maximum number of iterations allowed for each inner mean estimation procedures
#' @param maxitG the maximum number of iterations allowed for each inner dispersion estimation procedures
#' @param tol the level of accuracy required for each inner mean and dispersion estimation procedure to converge
#' @param acc the level of accuracy required for the outer iteration to converge 
#' @param maxitOUT the maximum number of outer iteration allowed in the Double estimation procedure
#' @param tccM the tuning constant c used in the robust estimation of the mean function
#' @param tccG the tuning constant c used in the robust estimation of the dispersion function
#' @param lambdaM range of acceptable smoothing parameter values for the mean estimation procedure 
#' @param lambdaG range of acceptable smoothing parameter values for the dispersion estimation procedure 
#' @return a list of control parameters
#' @export
DoubleRobGamControl <- function(maxitM=35,maxitG=35,tol=10^-7,acc=5*10^-3,maxitOUT=60,tccM=1.345,tccG=1.345,lambdaM=c(0.0001,3500),lambdaG=c(0.0001,3500)){
list(maxitM=maxitM,maxitG=maxitG,tol=tol,maxitOUT=maxitOUT,acc=acc,tccM=tccM,tccG=tccG,lambdaM=lambdaM,lambdaG=lambdaG)
}






#' Double Robust Generalized Additive Models 
#' 
#' A function to robustly estimate both the mean and the dispersion function using B-splines bases. 
#' The robustification is done in a similar fashion to the robustbase::glmrob function.  
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
#' @param selection the method used to select the smoothing parameter. 
#' can be "RGCV" (the default), "RAIC", "GCV", "AIC" or "none", in which case no selection is done and a smoothing parameter value should be provided in the \code{bs()} call. 
#' The RGCV and RAIC are robustified version of GCV and AIC (see Croux et al. (2012))
#' @param weights as in glm - to be used with care, since the data get always weigthed by the estimated variance function
#' @param control a list to control some behaviours of the estimation procedure, see {\link{DoubleRobGamControl}} 
#' @param trace should a trace to follow the convergence be printed? Default is FALSE
#' @param scale similar to weigth, to be used if the dispersion function should be kept fixed
#' @param weights.on.x mutated from robustbase::glmrob, a vector of weights for the x-variables. Used to accomodated for leverage points. 
#' 
#' 
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
#' \item{RGCV}{Roubust Generalised Cross Validation value}
#' \item{RAIC}{Roubust Akaike Information criterion value}
#' \item{w.x}{the weights on the covariates used to correct for the effect of leverage points}
#' \item{estRphi}{roubustly estimated dispersion parameter}
#' \item{tcc}{tuning constant c used in the robust estimation}
#' \item{yt}{the variable response variable transformed according to the link function and then centered}
#' \item{df}{overall equivalent degrees of freedom for the fit}
#' \item{vecdf}{a vector giving the degrees of freedom used by each covariate in the model}
#' \item{cov.coef}{covariance matrix of the coefficients}
#' \item{formula}{formula used}
#' \item{dimsP}{size of the parametric part of the model}
#' 
#' 
#' @seealso  \code{\link{DoubleRobGamControl}}, \code{\link{DoubleGam}}, \code{\link{gamMD}}
#' @author Ilaria Prosdocimi (ilapro@@ceh.ac.uk)
#' @references  Croux, C., Gijbels, I. and Prosdocimi, I. (2012), Robust Estimation of Mean and Dispersion Functions in Extended Generalized Additive Models. Biometrics, 68: 31-44. doi: 10.1111/j.1541-0420.2011.01630.x
#' @export 
DoubleRobGam <- function(formulaM,formulaG=NULL,family="gaussian", data, startM = NULL, startG = NULL, selection="RGCV",weights=NULL,control=DoubleRobGamControl(),trace=FALSE,scale=1,weights.on.x="none",method="quasi",...){
    if(method!="quasi" & method!="pseudo") warning('the given method argument was not recognized: use -quasi- instead')
    if(selection!="none" & selection!="GCV" & selection!="AIC" & selection!="RGCV" & selection!="RAIC"){
    print('how should I select the smoothing parameters? default is RGCV')
    selection<-"RGCV"
    }
    Glink="log" ## if the gamma is fitted, the link function is a log
    isgammaEst<-(!is.null(formulaG))
    missdat<-missing(data)
    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family <- family()
    freadM<-read.form(formulaM,data=data)
       sm.pM<-freadM$sm.p
       XM<-freadM$dataMAT
       interceptM<-freadM$int
       PstrM<-freadM$P
       y<-freadM$y
       dimsSM<-if(!length(freadM$smooths)) 0 else sapply(freadM$smooths,ncol)
       dimsTM<-c(ncol(PstrM)-sum(dimsSM),dimsSM)
    if (is.null(weights))
        weights <- rep.int(1, NROW(y))
    else if(any(weights <= 0))
        stop("All weights must be positive")
   tt<-(trace & !isgammaEst)
   fitM<-gamRobMeanIntern(sm.p=sm.pM,y=y,weights=weights,XM=XM,Pstr=PstrM,dimsS=dimsSM,dimsT=dimsTM,
        start=startM,family=family,trace=tt,acc=control$tol,maxit=control$maxitM,selection=selection,
        scale=scale,tcc=control$tccM,weights.on.x=weights.on.x,minlambda=min(control$lambdaM),maxlambda=max(control$lambdaM) )
   if(length(scale)==1 & family$family=="gaussian"){
       scale<-sum((y-fitM$fitted)^2)/(length(y)-fitM$df)
       fitM<-gamRobMeanIntern(sm.p=sm.pM,y=y,weights=weights,XM=XM,Pstr=PstrM,dimsS=dimsSM,dimsT=dimsTM,start=startM,
           family=family,trace=FALSE,acc=control$tol,maxit=control$maxitM,selection=selection,scale=scale,
           tcc=control$tccM,weights.on.x=weights.on.x,minlambda=min(control$lambdaM),maxlambda=max(control$lambdaM) )
    }
    fitM$convT<-fitM$conv
    fitM$formula<-formulaM;fitM$family<-family
    res<-fitM
    if(isgammaEst){
      if(formulaG == "~1"){
         if(control$tccG <= 0) gammaEST<-sum(scale*fitM$deviance)/(length(fitM$deviance)-fitM$df)
         if(control$tccG > 0)  gammaEST<-fitM$estRphi
         fitM<-gamRobMeanIntern(sm.p=sm.pM,y=y,weights=weights,XM=XM,Pstr=PstrM,dimsS=dimsSM,dimsT=dimsTM,
                 start=startM,family=family,trace=FALSE,acc=control$tol,maxit=control$maxitM,selection=selection,
                 scale=gammaEST,tcc=control$tccM,weights.on.x=weights.on.x,minlambda=min(control$lambdaM),maxlambda=max(control$lambdaM) )
         if(control$tccG <= 0) gammaEST<-sum(fitM$deviance)/(length(fitM$deviance)-fitM$df)
         if(control$tccG > 0)  gammaEST<-fitM$estRphi
         convT<-TRUE
         fitM$dispersion<-gammaEST
         fitM$convT<-convT
         fitM$formula<-formulaM;fitM$family<-family;fitM$dimsP<-freadM$dimsP
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
         else {freadG<-read.form(Gform,data=data)}### actually reading the formula, if the dataset is given
           sm.pG<-freadG$sm.p
           XG<-freadG$dataMAT
           interceptG<-freadG$int
           PstrG<-freadG$P
           ddev<-freadG$y
           dimsSG<-if(!length(freadG$smooths)) 0 else sapply(freadG$smooths,ncol)
           dimsTG<-c(ncol(PstrG)-sum(dimsSG),dimsSG)
       fitG<-gamRobDispIntern(sm.p=sm.pG,yd=dDev,XG=XG,Pstr=PstrG,dimsS=dimsSG,dimsT=dimsTG,start=startG,
              trace=FALSE,acc=control$tol,maxit=control$maxitG,selection=selection,link=Glink,
              tcc=control$tccG,weights.on.x=weights.on.x,weights=weights,minlambda=min(control$lambdaG),maxlambda=max(control$lambdaG) )
       if(trace) cat("The outer iteration begins","\n")
       for(nit in 1:control$maxitOUT) {
          oldcoefM<-fitM$coef #  NULL+runif(length(fitM$coef),0,fitM$coef/10)*sign(runif(length(fitM$coef.p))-0.5)
          oldcoefG<-fitG$coef #  NULL+runif(length(fitG$coef),0,fitG$coef/10)*sign(runif(length(fitG$coef.p))-0.5)
          startM<- fitM$coef
          startG<- fitG$coef
          meanEST<- fitM$fitted.values
          gammaEST<-fitG$fitted.values
          if(selection != "none"){
          sm.pM<-fitM$sm.p
          sm.pG<-fitG$sm.p
          }
          gammaEST[gammaEST<10^-5]<-10^-5
          fitM<-gamRobMeanIntern(sm.p=sm.pM,y=y,weights=weights,XM=XM,Pstr=PstrM,dimsS=dimsSM,dimsT=dimsTM,
                 start=startM,family=family,trace=FALSE,acc=control$tol,maxit=control$maxitM,selection=selection,
                 scale=gammaEST,tcc=control$tccM,weights.on.x=weights.on.x,minlambda=min(control$lambdaM),maxlambda=max(control$lambdaM) )
          if(method=="pseudo") resp<-fitM$residuals^2
          else  resp<-fitM$dispersion*fitM$deviance
          fitG<-gamRobDispIntern(sm.p=sm.pG,yd=resp,XG=XG,Pstr=PstrG,dimsS=dimsSG,dimsT=dimsTG,start=startG,
                trace=FALSE,acc=control$tol,maxit=control$maxitG,selection=selection,link=Glink,
                tcc=control$tccG,weights.on.x=weights.on.x,weights=weights,minlambda=min(control$lambdaG),maxlambda=max(control$lambdaG) )
#           relE <- c((sum((fitM$coef-startM)^2)/max(1e-20, sum(startM^2))),
#                          (sum((fitG$coef-startG)^2)/max(1e-20, sum(startG^2))))
	  relE <- c(sqrt(sum((fitM$fitted-meanEST)^2)/max(1e-20, sum(meanEST^2))),
                         sqrt(sum((fitG$fitted-gammaEST)^2)/max(1e-20, sum(gammaEST^2))))
          if(trace) {
	  cat("iteration", nit,"\n")
	  cat("the chosen smoothing parameters are","\n")
	  cat(sm.pM,"for the mean -",fitM$df,"degrees of freedom","\n")
	  cat(sm.pG,"for the dispersion -",fitG$df,"degrees of freedom","\n")
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
   fitM<- fitM[seq(1:length(fitM))[-which(names(fitM) == "Mmat")]]
   fitG<- fitG[seq(1:length(fitG))[-which(names(fitG) == "Mmat")]]
   fitM<- fitM[seq(1:length(fitM))[-which(names(fitM) == "Qmat")]]
   fitG<- fitG[seq(1:length(fitG))[-which(names(fitG) == "Qmat")]]
   res<-list(fitM=fitM,fitG=fitG,convVec=convVec,converged=convT,iter=nit,relE=relE)
   if(missdat) rm(dDev,envir = data)
     }
   }
if(missdat) res$data<-1
if(!missdat) res$data<-data
class(res)<-"gamMD"
res
}



