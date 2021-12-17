## ----eval=TRUE----------------------------------------------------------------

tls<-function(X,y,W,q) {
  n<-nrow(X)
  Tn<-numeric(n)
  fit1<-list()
  X<-as.matrix(X)
#first stage
  for (i in 1:n) {
    #searching for gamma
    gamma<-q[i]
    D1<-matrix(0,n,n)
    diag(D1)<-as.numeric(q<=gamma)
    D2<-diag(x=1,n,n)-D1
    #instrument variable
    Z<-cbind(X,W%*%X,W%*%D1%*%W%*%X,W%*%D2%*%W%*%X)
    fit1[[i]]<-lm(W%*%y~Z+0)
  }
  Wy_hat<-matrix(0,n,n)
  for (i in 1:n) {
    Wy_hat[,i]<-fit1[[i]]$fitted.values
  }
#second stage
  fit2<-list()
  estimates<-list()
  sse<-numeric(n)
  for (i in 1:n){
     X_1<-cbind(D1%*%Wy_hat[,i],D2%*%Wy_hat[,i],X) 
     fit2[[i]]<-lm(y~X_1+0)
     estimates[[i]]<-fit2[[i]]$coefficients
     sse[i]<-t(fit2[[i]]$residuals)%*%(fit2[[i]]$residuals)
  }
  gammahat<-q[which(sse==min(sse))][1]
  return(c(gammahat,unlist(slopecoefhat<-estimates[which(sse==min(sse))][1])))
}


## ----eval=TRUE----------------------------------------------------------------

sml<-function(t,n,beta1,beta2,delta,gam){
  d<-5
  g<-matrix(0,t,d)
  W<-matrix(0,n,n)
  for (i in 1:n) {
    for (j in 1:n) {
      if(abs(j-i)<=3)
        W[i,j]<-1
    }
  }
  W <- sweep(W, 1, apply(W, 1, function(x) sqrt(sum(x^2))), "/")
  for (i in 1:t) {
    X1<-rnorm(n);X2<-rnorm(n);X<-cbind(X1,X2)
    q<-rnorm(n,2,1);u<-rnorm(n)
    D1<-matrix(0,n,n)
    diag(D1)<-as.numeric(q<=gam)
    D2<-diag(x=1,n,n)-D1
    
    Lambda<-solve(diag(x=1,n,n)-beta1*D1%*%W-beta2*D2%*%W)
    y<-Lambda%*%(X%*%delta+u)
    
    fit<-tls(X,y,W,q)
    g[i,1]<-fit[1]
    g[i,2:5]<-fit[2:5]
  }  
  g<-data.frame(g)
  for (i in 1:5) {
    result<-lapply(g,function(x) quantile(x,probs = c(0.05,0.5,0.95),na.rm=TRUE)
    )
  }
  result<-t(as.data.frame(result))
  rownames(result)<-c("gamma","beta1","beta2","delta1","delta2")
  knitr::kable(result)
}


## ----eval=TRUE----------------------------------------------------------------

t<-1000;n<-100
gamma<-2;beta1<-0.2;beta2<-0.5;delta<-c(1,1);
sml(t,n,beta1,beta2,delta,gamma)


