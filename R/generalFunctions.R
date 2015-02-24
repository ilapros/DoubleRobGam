
### functions used in the DoubleRobGam library

######## All functions written by Ilaria Prosdocimi (Ilaria.Prosdocimi@wis.kuleuven.be)

#' Builder for B-splines basis and penalty matrices typically used in DoubleGam and DoubleRobGam
#' @param x data vector
#' @param nknots number of knots used in building the B-splines
#' @param p polynomial degree of the basis
#' @param center logical indicator on whether the design matrix should be centered (and it should be the case when more than one covariate enter the model) 
#' @param sm.p smoothing parameter used in the penalty matrix build
#' @param order penalty matrix order 
#' 
#' @export
bsp<-function(x,nknots=15,p=3,center=TRUE,sm.p=1,order=2) {
n<-length(x)
## Bsplines built like in Eilers and Marx (1996)
maxx<-(max(x)+0.001);minn<-(min(x)-0.001)
dx<-(maxx-minn)/nknots
kn<-seq(minn-p*dx,maxx+p*dx,by=dx)
BB<-splineDesign(knots=kn,x=x,ord=p+1,0*x,outer.ok=TRUE)
if(center) BB<-(diag(n)-rep(1,l=n)%*%t(rep(1,l=n))/n)%*%BB
### this is needed in GAM, when having only one bspline it could be avoided, but I'd rather not...
p<-ncol(BB)
P<-diag(p)
if(order>0) P<-t(diff(P,diff=order))%*%diff(P,diff=order)
### build the structure of the P matrix
res<-list(B=BB,P=P,sm.p=sm.p)
}


read.form<-function(ff,data){
    int<-1
    if (missing(data)) data <- environment(ff)
     else {ddaatt <- NULL; data <- data[,all.vars(ff)]}
    ft<-terms.formula(ff,specials=c("bsp","I"))
    termsl<-attr(ft,"term.labels")
    nt<-length(termsl)
    lNAl<-NULL       ### I know this is not a good way to program this, should be improved...
    for(k in 1:length(all.vars(ff))){
      oo<-eval(parse(text=all.vars(ff)[k]),envir = data)
      lNAl<-unique(c(lNAl,seq(1,length(oo))[is.na(oo)]))
      }
    if(length(lNAl)) {
      listOR<-list()
       for(k in 1:length(all.vars(ff))){
       nn<-all.vars(ff)[k]
       oo<-eval(parse(text=all.vars(ff)[k]),envir = data)
       listOR[[k]]<-oo
       if(is.environment(data)) assign(nn,oo[-lNAl],envir = data)
          else  ddaatt<-cbind(ddaatt,oo[-lNAl])
       }
    }
    else{ if(!is.environment(data)) ddaatt<-data}
    if(!is.environment(data)) {
       data<-as.data.frame(ddaatt)
       colnames(data)<-all.vars(ff)
    }
    y<-eval(parse(text=attr(ft,"variables")[2]),envir = data);n<-NROW(y)
    indsmooth<-attr(ft,"specials")$bsp-1;lensmooth<-length(indsmooth)
    indNsmoothI<-attr(ft,"specials")$I-1
    indNsmooth<-seq(1,nt)
    for(i in sort(c(indsmooth,indNsmoothI),decreasing=T)) indNsmooth<-indNsmooth[-(i)]
    lenNsmooth<-length(indNsmooth)
    dataMAT<-matrix(rep(1,l=n),ncol=1);P<-matrix(0,ncol=1,nrow=1)
    sm.pS<-NULL; smoothts<-NULL;paramts<-NULL;dimsP<-0
    if(length(indNsmooth)){ ## building the parameteric part
      allParsInd<-parMat<-names<-allfirstOrd<-allhigherOrd<-NULL;ss<-list()
      j<-0
      for(i in (indNsmooth)){
          j<-j+1
          firstOrd<-termsl[i]
          obj<-scale(eval(as.expression(parse(text=firstOrd)),envir = data));names<-c(names,firstOrd)
          parMat<-cbind(parMat,obj)
          allfirstOrd<-c(allfirstOrd,firstOrd)
          allParsInd<-c(allParsInd,j)
      }
      for(i in (indNsmoothI)){
          j<-j+1
#           a<-(strsplit(termsl[i],'\\^')[[1]])[1]
#           b<-substr(a,start=3,stop=nchar(a))
          a<-(strsplit(termsl[i],'\\^')[[1]])
          b<-substr(a[1],start=3,stop=nchar(a[1]))
          p<-strsplit(a[2],')')[[1]]
          allhigherOrd<-c(allhigherOrd,b)
          dd<-scale(eval(parse(text=b),envir = data))^as.numeric(p)
          obj<-matrix(dd,nrow=n);names<-c(names,b)
          parMat<-cbind(parMat,obj)
#           ;colnames(parMat[,ncol(parMat)])<-b
          allParsInd<-c(allParsInd,j)
      }
      colnames(parMat)<-names
      allPars<-c(allfirstOrd,allhigherOrd);uniquePars<-unique(c(allfirstOrd,allhigherOrd))# ;allParsInd<-c(indNsmooth,indNsmoothI)
      for(i in 1:length(uniquePars)){ss[[i]]<-allParsInd[allPars==uniquePars[i]]}
      names(ss)<-uniquePars
      dimsP<-(unlist(lapply(ss,length)))
      parMat<-as.matrix(parMat[,unlist(ss)])
     if(!length(indsmooth)) paramts[[j]]<-obj
     for(i in 1:length(ss)){
     if(!length(indsmooth))  obj<-as.matrix(parMat[,(1+sum(dimsP[0:(i-1)])):(sum(dimsP[0:i]))])
       else obj<-as.matrix((diag(n)-rep(1,l=n)%*%t(rep(1,l=n))/n)%*%parMat[,(1+sum(dimsP[0:(i-1)])):(sum(dimsP[0:i]))])
     paramts[[i]]<-obj
     dataMAT<-cbind(dataMAT,obj)
     Pj<-diag(0,ncol(obj))
     zeros<-matrix(0,ncol=ncol(P),nrow=ncol(Pj))
     P<-cbind(rbind(P,zeros),rbind(t(zeros),Pj));j<-j+1
     }
    }
    j<-1 ## building the NON parameteric part
    for(i in (indsmooth)){
        obj<-eval(parse(text=termsl[i]),envir = data)
        smoothts[[j]]<-matrix(obj$B,nrow=n)
        dataMAT<-cbind(dataMAT,smoothts[[j]])
        Pj<-obj$P
        sm.pS<-c(sm.pS,obj$sm.p)
        zeros<-matrix(0,ncol=ncol(P),nrow=ncol(Pj))
        P<-cbind(rbind(P,zeros),rbind(t(zeros),Pj));j<-j+1
    }
    if(!attr(ft,"intercept")) {
        dataMAT<-dataMAT[,-1];P<-matrix(P[-1,],ncol=ncol(P));P<-P[,-1]
        int<-0
        }
    if(length(lNAl) & is.environment(data)) { ### putting missing values back, when needed
    for(k in 1:length(all.vars(ff))){
        nn<-all.vars(ff)[k]
        assign(nn,listOR[[k]],envir = data)
        }
    }
    if(!length(lNAl)) lNAl<-NULL# FALSE
    list(dataMAT=dataMAT,params=paramts,smooths=smoothts,y=y,P=P,intercept=int,sm.p=sm.pS,seqNA=lNAl,dimsP=dimsP)
}








### used when building the penalty matrices
matNum<-function(num,dim1,dim2) matrix(num,ncol=dim2,nrow=dim1)





######### functions to generate over/under dispersed data

#' Random generator of Double Exponantial Poisson-type data (counts)
#' 
#' @param size sample size needed
#' @param lambda mean parameter/function
#' @param gamma dispersion parameter/function. \code{gamma = 1} (the default) corresponds to the Poisson distribution
#' 
#' @export
rdoublepois<-function(size=1,lambda,gamma=1){
    theta<-gamma;rho<-1/theta
    pp<-NULL;ss<-NULL;
    lambda<-rep(lambda,l=size)
    rho<-rep(1/gamma,l=size)
    for(i in 1:size){
       ll<-lambda[i];tt<-rho[i]
       for(k in 1:141){
           y<-(k-1)
           pp[k]<-(tt^(1/2))*(exp(-y)*(y^y)/factorial(y))*(exp(1)*ll/y)^(tt*y)
           probpos<-length(pp[pp!=Inf]) ### possible values
       }
    ss[i]<-sample(seq(0,(probpos-1)),size=1,prob=pp[1:(probpos)])
    }
    ss
}



#' Random generator of Double Exponantial Binomial-type data (counts)
#' 
#' @param size sample size needed
#' @param n number of trials (called \code{size} in \code{rbinom})
#' @param prob vector of probabilities
#' @param theta dispersion parameter/function. \code{theta = 1} (the default) corresponds to the Binomial distribution
#' 
#' @export
#' @importFrom splines splineDesign
#' @importFrom MASS ginv
#' @importFrom robustbase covMcd
rdoublebinom<-function(size,n,prob,theta){
prob<-rep(prob,l=size);res<-list()
theta<-rep(theta,l=size)
n<-rep(n,l=size);p<-NULL
ss<-NULL
for( i in 1:size){
y<-seq(0,n[i]);p<-NULL
for(j in 1:length(y)){
p[j]<-sqrt(1/theta[i])*choose(n[i],y[j])*((((prob[i])^y[j])*((1-prob[i])^(n[i]-y[j])))^(1/theta[i]))  *((((y[j]/n[i])^y[j])*((1-y[j]/n[i])^(n[i]-y[j])))^(1-1/theta[i]))
}
samp<-sample(y,size=1,prob=p)
ss<-c(ss,samp)
}
res$sample<-ss
res$n<-n
res
}




# library(MASS);library(splines)
# importFrom splines splineDesign
# importFrom MASS ginv
# autoImports
# import MASS
# import splines

