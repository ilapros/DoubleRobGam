## ----echo=FALSE----------------------------------------------------------
library(DoubleRobGam)
# source("..//R//DoubleRobGam.R")
# source("..//R//robFunctions.R")
# source("..//R//generalFunctions.R")

## ----ragweeed,cache=TRUE-------------------------------------------------
library(SemiPar);data(ragweed)
ragweed<-ragweed[ragweed$year==1991,]  # day.in.seas
ragweed<-ragweed[order(ragweed$day.in.seas),] # for plotting

rag1a<-DoubleGam(formulaM=ragweed ~ day.in.seas+I(day.in.seas^2),
       data=ragweed,family="poisson",selection="none")

## ----eval=TRUE,cache=TRUE------------------------------------------------
rag1b <- DoubleGam(ragweed ~ bsp(day.in.seas),
                   data = ragweed, family="poisson")

## ----plotWeed------------------------------------------------------------
par(mfrow=c(1,1), bty = "l", lwd=2)
plot(ragweed$day.in.seas, ragweed$ragweed ,pch=19, col=8)
lines(ragweed$day.in.seas, rag1a$fitted.values, lty=2, col=4)
lines(ragweed$day.in.seas, rag1b$fitted.values)
legend("topright", col=c(4,1), lty = c(2,1), bty = "n", 
      legend = c("quadratic fit", "non-parametric fit"))

## ----echo=TRUE,message=FALSE,eval=TRUE-----------------------------------
args(bsp)

## ----eval=TRUE,cache=TRUE------------------------------------------------
rag2<-DoubleGam(ragweed~bsp(day.in.seas),
                formulaG= ~bsp(day.in.seas,nknots=14), 
                data=ragweed, family="poisson")

## ----eval=TRUE,cache=TRUE------------------------------------------------
names(rag2)
names(rag2$fitM)
names(rag2$fitG)

## ----RagweedMeanDisp, fig.show='asis',fig.align='center',cache=TRUE------
par(mfrow=c(2,1),mai = c(0.3,0.3,0.3,0.3))
plot(ragweed$day.in.seas,rag2$fitM$fitted, type= "l")
### add the result obtained if no dispersion is estimated
lines(ragweed$day.in.seas,rag1b$fitted,lty=2, type= "l", col=4)
plot(ragweed$day.in.seas,rag2$fitG$fitted, type= "l")

## ----RagweedMeanDispDirect, fig.show='asis',fig.align='center'-----------
plot(rag2)

## ------------------------------------------------------------------------
class(rag2)
summary(rag2)
rag2

## ----collapse=TRUE,cache=TRUE--------------------------------------------
library(SemiPar);data(ustemp)
temp1<-DoubleGam(formulaM=min.temp~bsp(latitude,p=2,nknots=8,sm.p=0.006)+
                           bsp(longitude,p=2,nknots=7,sm.p=0.034),
                 formulaG=~bsp(latitude,p=2,nknots=6,sm.p=1400)+
                           bsp(longitude,p=2,nknots=6,sm.p=17000),
                 selection="GCV",trace=TRUE,data=ustemp)

## ----cache=TRUE----------------------------------------------------------
temp1<-DoubleGam(formulaM=min.temp~bsp(latitude,nknots=9,p=2)+
                           bsp(longitude,nknots=9,p=2,sm.p=0.2),
                 formulaG=~bsp(latitude,p=2,nknots=8)+
                           bsp(longitude,p=2,nknots=8,sm.p=0.54),
                 data=ustemp,
          control=DoubleGamControl(maxitM=25,maxitG=25,maxitOUT=5))

## ----cache=TRUE----------------------------------------------------------
temp1<-DoubleGam(formulaM=min.temp~bsp(latitude,nknots=9,p=2)+
                                   bsp(longitude,nknots=9,p=2,sm.p=0.2),
                 formulaG=~bsp(latitude,p=2,nknots=8)+
                           bsp(longitude,p=2,nknots=8,sm.p=0.54),
                 data=ustemp,
      control=DoubleGamControl(maxitM=25,maxitG=25,maxitOUT=5,acc=0.07))

## ----cache=TRUE----------------------------------------------------------
temp1 <- DoubleGam(formulaM=min.temp~bsp(latitude,nknots=9,p=2)+
                                   bsp(longitude,nknots=9,p=2,sm.p=0.2),
                 formulaG=~bsp(latitude,p=2,nknots=8)+
                           bsp(longitude,p=2,nknots=8,sm.p=0.54),
                 data=ustemp,
         control=DoubleGamControl(lambdaM=c(.05,200),
                                  lambdaG=c(.3,14),maxitOUT=4))


## ------------------------------------------------------------------------
args(DoubleGamControl)

## ----cache=TRUE----------------------------------------------------------
tempRob1<-DoubleRobGam(min.temp~bsp(latitude,nknots=9,p=2)+
                                bsp(longitude,nknots=9,p=2),
          formulaG=~bsp(latitude,p=2,nknots=8)+
          bsp(longitude,p=2,nknots=8),data=ustemp)

## ------------------------------------------------------------------------
summary(tempRob1)

## ----ustempMeanVarRob, fig.show='asis',fig.align='center'----------------
plot(tempRob1) 

## ----cache=TRUE----------------------------------------------------------
tempRob2<-DoubleRobGam(min.temp~bsp(latitude,nknots=9,p=2)
    +bsp(longitude,nknots=9,p=2),
    formulaG=~bsp(latitude,p=2,nknots=8,sm.p=500)+
    bsp(longitude,p=2,nknots=8,sm.p=500),data=ustemp,
    control=DoubleRobGamControl(tccM=1.2,tccG=1.8))

## ----USrob2, fig.show='asis',fig.align='center'--------------------------
plot(tempRob2, ci.plot = FALSE)

