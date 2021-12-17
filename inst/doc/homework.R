## -----------------------------------------------------------------------------
#print('Hello World!')

## -----------------------------------------------------------------------------
dosage <- c(20, 30, 40, 45, 60) 
drugA <- c(16, 20, 27, 40, 60) 
drugB <- c(22, 30, 35, 41, 50 )
plot(dosage,drugA,type = 'b',pch=15,xlab="Drug Dosage", ylab="Drug Response")
lines(dosage, drugB, type="b",pch=17) 
legend("topleft", title="Drug Type", c("A","B"),pch=c(15, 17))

## -----------------------------------------------------------------------------
knitr::kable(head(mtcars))

## -----------------------------------------------------------------------------
rRayleigh=function(n,sigma){
  u <- runif(n)
  x <- (-2*sigma^2*log(1-u))^{1/2} 
  #F(x)=1-e^(-x^2/2*sigma^(1/2)), 0<=x<=1
hist(x, prob = TRUE,main =paste('sigma=',as.character(sigma),sep='')) 
  y <- seq(0,sigma*5, .01)
  lines(y,(y/sigma^2)*exp(-(1/2)*y^2/sigma^2))
}
for (i in 1:5) {
  rRayleigh(1e4,i)
}



## -----------------------------------------------------------------------------
Mix<-function(n,p){
  X1 <- rnorm(n)
  X2 <- rnorm(n,3,1)
  r <- sample(c(0,1),n,replace=TRUE,prob = c(1-p,p))
  Z <- r*X1+(1-r)*X2
  hist(Z,probability = TRUE,main = paste('p=',as.character(p),sep=''))
  
}
Mix(1e4,0.75)

for (i in seq(0.1,0.9,0.1)){
  Mix(1e4,i)
}

## -----------------------------------------------------------------------------

pgp=function(lambda,a,b,t){
  
  N=rpois(1,lambda*t)
  y<-numeric(length = N)
  for (i in 1:N) {
    y[i]<-rgamma(1,a,b)
  }
  Xt<-sum(y)
  return(Xt)
}

#lamda=3,alpha=2,beta=3
pgpx<-function(x){
  Z<-numeric(length = x)
  for (i in 1:x) {
    Z[i]<-pgp(3,2,3,i)
  }
  return(Z)
}

pgp50<-pgpx(50)
m<-seq(1,50,1)
plot(m,pgp50,type ="l",col="red",xlab='t',ylab='X(t)',main = 'Smulation of compound Poisson(λ)–Gamma process')


#Estimate the mean and the variance of X(10) for
#several choices of the parameters and compare with the theoretical values.
test<-function(l,a,b){
  x1<-numeric(length = 10)
   for (i in 1:10) {
      x1[i]<-pgp(l,a,b,10)
    }
   m1<-mean(x1)
   mt1<-l*10*(a/b)
   v1<-var(x1)
   vt1<-l*10*(a/(b^2)+(a/b)^2)
   print(paste("empirical mean is ",m1))
   print(paste("theoretical mean is",mt1))
   print(paste("empirical variance is ",v1))
   print(paste("theoretical variance is",vt1))
}

#lamda=3,alpha=2,gamma=3,t=10
test(3,2,3)
#lamda=4,alpha=5,gamma=6,t=10
test(4,5,6)


## -----------------------------------------------------------------------------
m <- 1e4;
mybetacdf<-function(x){
  t <- runif(m, min=0, max=x)
  theta.hat <- mean((1/beta(3,3))*t^2*(1-t)^2) * x
  return(theta.hat)
}

for (i in seq(0.1,0.9,0.1)) {
  print(c(mybetacdf(i),pbeta(i,3,3)))
}


## -----------------------------------------------------------------------------
set.seed(1234)
MC.Phi <- function(sigma, R = 10000, antithetic = TRUE) {
  u <- runif(R/2,0,1)
  if (!antithetic) v <- runif(R/2) else v <- 1 - u
  u <- c(u, v)
  x <- (-2*sigma^2*log(1-u))^{1/2} 
  x
}

MC1 <- MC.Phi(3, R = 1000, anti = FALSE)
MC2 <- MC.Phi(3, R = 1000)
print(paste("the variance without antithetic variables is ",var(MC1)))
print(paste("the variance using antithetic variables is ",var(MC2)))
print((var(MC1)-var(MC2))/(var(MC1)))


## -----------------------------------------------------------------------------

g<-function(x){
  ((x^2)/sqrt(2*pi))*exp(-x^2/2)
}
f1<-function(x){
  (1/sqrt(2*pi))*exp(-((x-2)^2)/2)
}

f2<-function(x){
  dgamma(x,2,2)
}

set.seed(1234)
m <- 10000
theta.hat <- se <- numeric(2)

x <- rnorm(m,2,1) #using f1
fg <- g(x)/f1(x)
theta.hat[1]<-mean(fg)
se[1]<-sd(fg)

x <- rgamma(m,2,2) #using f1
fg <- g(x)/f2(x)
theta.hat[2]<-mean(fg)
se[2]<-sd(fg)

rbind(theta.hat,se)


## -----------------------------------------------------------------------------

M <- 10000; k <- 10 
r <- M/k #replicates per stratum
N <- 50 #number of times to repeat the estimation
T2 <- numeric(k)
est <- matrix(0, N, 2)
g<-function(x)
{
  x^2/(sqrt(2*pi))*exp(-x^2/2)
}

for (i in 1:N) {
#odinary
est[i, 1] <- 1/2-mean(g(runif(M)))
#statified
#jth stratum
for(j in 1:k)T2[j]<-mean(g(runif(M/k,(j-1)/k,j/k)))
est[i, 2] <- 1/2-mean(T2)
}
round(apply(est,2,mean),4)
round(apply(est,2,sd),5)

## -----------------------------------------------------------------------------
n <- 20
alpha <- .05
UCL<-LCL<-numeric(1000)
set.seed(100)
for (i in 1:1000) {
  x <- rchisq(n, df = 2)
  LCL[i]<-mean(x)-qt(1-alpha/2,df = (n-1))*sqrt(var(x))/sqrt(n)
  UCL[i]<-mean(x)+qt(1-alpha/2,df = (n-1))*sqrt(var(x))/sqrt(n)

}

sum(UCL > 2 & LCL<2)
mean(UCL > 2 & LCL<2)


## -----------------------------------------------------------------------------
set.seed(123)
n <- 20
alpha <- .05
m <- 10000 #number of replicates

#chi(1)
mu1 <- 1
p1 <- numeric(m) #storage for p-values
for (j in 1:m) {
x <- rchisq(n,df=1)
ttest <- t.test(x, mu = mu1)
p1[j] <- ttest$p.value # t-test p-value
}
p1.hat <- mean(p1 < alpha)

#Uniform(0,2)
mu2<- 1
p2 <- numeric(m) #storage for p-values
for (j in 1:m) {
x <- runif(n,0,2)
ttest <- t.test(x, mu = mu2)
p2[j] <- ttest$p.value # t-test p-value
}
p2.hat <- mean(p2 < alpha)

#Exponential(1)
mu3<- 1
p3 <- numeric(m) #storage for p-values
for (j in 1:m) {
x <- rexp(n,1)
ttest <- t.test(x, mu = mu3)
p3[j] <- ttest$p.value # t-test p-value
}
p3.hat <- mean(p3 < alpha)

print(c(p1.hat,p2.hat,p3.hat))


## ----eval=FALSE---------------------------------------------------------------
#  d<-3 # set dimension=3
#  n <- c(10, 20, 30, 50, 100, 500)
#  cv<-qchisq(0.95,d*(d+1)*(d+2)/6) #crit. values for each n
#  print(cv)

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(1234)
#  sk <- function(x) {
#  #computes the sample skewness statistics.
#  Sigma<-cov(x)
#  xbar <- apply(x, 2, mean)
#  b<-0
#  for (p in 1:n[i]) {
#    for (q in i:n[i]) {
#      b<-b+((x[p,]-t(xbar))%*%solve(Sigma)%*%t(x[q,]-t(xbar)))^3
#    }
#  
#  }
#  b<-b/(n[i]^2)
#  return( n[i]*b/6 )
#  }
#  
#  #we are doing length(n) different simulations
#  p <- numeric(length(n)) #to store sim. results
#  m <- 100 #num. repl. each sim.
#  
#  for (i in 1:length(n)) {
#  sktests <- numeric(m) #test decisions
#    for (j in 1:m) {
#      #generate multivariate normal random vector
#      library(MASS)
#      Sigma0 <- matrix(c(1,0,0,0,1,0,0,0,1),3,3)
#      mean0<-c(0,0,0)
#      x<-mvrnorm(n[i],mean0,Sigma0)
#  
#      #test decision is 1 (reject) or 0
#       sktests[j] <- as.integer(sk(x) >= cv)
#    }
#  p[i] <- mean(sktests) #proportion rejected
#  }
#  
#  result<-matrix(nrow = 2,ncol = length(n))
#  result[1,]<-n
#  result[2,]<-p
#  dimnames(result)[[1]]<-c("n","estimate")
#  knitr::kable(result)

## ----eval=FALSE---------------------------------------------------------------
#  
#  set.seed(123)
#  alpha <- .1
#  n <- 30
#  m <- 100
#  epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
#  N <- length(epsilon)
#  pwr <- numeric(N)
#  
#  
#  sk2 <- function(x) {
#  #computes the sample skewness statistics.
#  Sigma<-cov(x)
#  xbar <- apply(x, 2, mean)
#  b<-0
#  
#  for (p in 1:n) {
#    for (q in 1:n) {
#      b<-b+((x[p,]-t(xbar))%*%solve(Sigma)%*%t(x[q,]-t(xbar)))^3
#    }
#  }
#  b<-b/(n^2)
#  return( n*b/6 )
#  }
#  #critical value for the skewness test
#  
#  Sigma0 <- matrix(c(1,0,0,0,1,0,0,0,1),3,3)
#  mean0=c(0,0,0)
#  library(MASS)
#  
#  for (j in 1:N) { #for each epsilon
#  e <- epsilon[j]
#  sktests <- numeric(m)
#    for (i in 1:m) { #for each replicate
#      x<-matrix(0,n,d)
#      for (z in 1:n) {
#        k <- sample(c(1, 10), replace = TRUE,
#        size = 1, prob = c(1-e, e))
#        x[z,]<-mvrnorm(1,mean0,k*Sigma0)
#      }
#  
#    sktests[i] <- as.integer(abs(sk2(x)) >= cv)
#  }
#  pwr[j] <- mean(sktests)
#  }
#  

## ----eval=FALSE---------------------------------------------------------------
#  
#  plot(epsilon, pwr, type = "b",xlab = bquote(epsilon), ylim = c(0,1))
#  abline(h = .1, lty = 3)
#  se <- sqrt(pwr * (1-pwr) / m) #add standard errors
#  lines(epsilon, pwr+se, lty = 3)
#  lines(epsilon, pwr-se, lty = 3)
#  

## -----------------------------------------------------------------------------

set.seed(123)
library(boot); library(bootstrap);set.seed(12345)

b.comp1 <- function(x,i) {
  ev <- eigen(cov(x[i,]), symmetric = TRUE)
  lambda <- ev$values
  max(lambda)/(sum(lambda))
  
}

x<-scor
n <- nrow(x)
obj <- boot(data=x,statistic=b.comp1,R=200)
round(c(original=obj$t0,bias=mean(obj$t)-obj$t0,
se=sd(obj$t)),4)



## -----------------------------------------------------------------------------
set.seed(123)
# Jackknife
theta.hat<-obj$t0
theta.jack <- numeric(n)

for(i in 1:n){
theta.jack[i] <- b.comp1(x,(1:n)[-i])
}
bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)
se.jack <-(n-1)*sqrt(var(theta.jack) / n)

round(c(original=theta.hat,bias.jack=bias.jack,se.jack=se.jack),4)

## -----------------------------------------------------------------------------
library(boot)
boot.ci(obj, conf = 0.95, type = c("perc", "bca"))

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(123)
#  #norm
#  sknorm<-0
#  m<-1000
#  ci.norm<-ci.basic<-ci.perc<-ci.bca<-matrix(NA,m,2)
#  
#  b.sk <- function(x,i) {
#    xbar <- mean(x[i])
#    m3 <- mean((x[i] - xbar)^3)
#    m2 <- mean((x[i] - xbar)^2)
#    return( m3 / m2^1.5 )
#  }
#  for (i in 1:m) {
#    x<-rnorm(1000)
#    n <- nrow(x)
#    obj <- boot(data=x,statistic=b.sk,R=200)
#  
#    ci<-boot.ci(obj, conf = 0.95, type = c("norm","basic","perc"))
#  
#    ci.norm[i,]<-ci$norm[2:3]
#    ci.basic[i,]<-ci$basic[4:5]
#    ci.perc[i,]<-ci$percent[4:5]
#  }
#  
#  cat('norm =',mean(ci.norm[,1]<=sknorm & ci.norm[,2]>=sknorm),
#  'basic =',mean(ci.basic[,1]<=sknorm & ci.basic[,2]>=sknorm),
#  'perc =',mean(ci.perc[,1]<=sknorm & ci.perc[,2]>=sknorm))
#  

## ----eval=FALSE---------------------------------------------------------------
#  #The proportion of times that the confidence intervals miss on the left, and
#  cat('norml =',mean(sknorm<=ci.norm[,1]),
#  'basicl =',mean(sknorm<=ci.basic[,1]),
#  'percl =',mean(sknorm<=ci.perc[,1]))
#  
#  #The porportion of times that the confidence intervals miss on the right.
#  cat('normr =',mean(sknorm>=ci.norm[,2]),
#  'basicr =',mean(sknorm>=ci.basic[,2]),
#  'percr =',mean(sknorm>=ci.perc[,2]))
#  

## ----eval=FALSE---------------------------------------------------------------
#  #chisq(5)=Gamma(5/2,1/2)
#  set.seed(123)
#  skchisq5<-((2*(5/2))/(1/2)^3)/((10)^1.5)
#  
#  for (i in 1:m) {
#    x<-rchisq(1000,5)
#    n <- nrow(x)
#    obj <- boot(data=x,statistic=b.sk,R=200)
#  
#    ci<-boot.ci(obj, conf = 0.95, type = c("norm","basic","perc"))
#  
#    ci.norm[i,]<-ci$norm[2:3]
#    ci.basic[i,]<-ci$basic[4:5]
#    ci.perc[i,]<-ci$percent[4:5]
#  }
#  
#  cat('norm =',mean(ci.norm[,1]<=skchisq5 & ci.norm[,2]>=skchisq5),
#  'basic =',mean(ci.basic[,1]<=skchisq5 & ci.basic[,2]>=skchisq5),
#  'perc =',mean(ci.perc[,1]<=skchisq5 & ci.perc[,2]>=skchisq5))
#  

## ----eval=FALSE---------------------------------------------------------------
#  #The proportion of times that the confidence intervals miss on the left,
#  
#  cat('norml =',mean(skchisq5<=ci.norm[,1]),
#  'basicl =',mean(skchisq5<=ci.basic[,1]),
#  'percl =',mean(skchisq5<=ci.perc[,1]))
#  # The porportion of times that the confidence intervals miss on the right.
#  cat('normr =',mean(skchisq5>=ci.norm[,2]),
#  'basicr =',mean(skchisq5>=ci.basic[,2]),
#  'percr =',mean(skchisq5>=ci.perc[,2]))

## -----------------------------------------------------------------------------

attach(chickwts) # chicken weights for various
x <- as.vector(weight[feed == "sunflower"])
y <- as.vector(weight[feed == "linseed"])
detach(chickwts)

z <- c(x, y)

R <- 999;
K <- 1:24;
n<-length(x);
set.seed(12345)
reps <- numeric(R);
t0 <- cor.test(x, y,method = "spearman")$statistic

for (i in 1:R) {
  k <- sample(K, size = n, replace = FALSE)
  x1 <- z[k]
  y1 <- z[-k] #complement of x1
  reps[i] <- cor.test(x1, y1,method = "spearman")$statistic
}

p<- mean(abs(c(t0, reps)) >= abs(t0))
round(c(p,cor.test(x,y)$p.value),3)


## ----eval=FALSE---------------------------------------------------------------
#  library(RANN)
#  library(boot)
#  library(energy)
#  library(Ball)
#  
#  #n1 = n2 = 50; mu1 = mu2 = c(0,0);
#  
#  m <- 1e3;
#  k<-3; p<-2;
#  mu <- 0.3;
#  set.seed(12345)
#  n1 <- n2 <- 50;
#  R<-999;
#  n <- n1+n2;
#  N = c(n1,n2)
#  
#  Tn <- function(z, ix, sizes,k) {
#    n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
#    if(is.vector(z)) z <- data.frame(z,0);
#    z <- z[ix, ];
#    NN <- nn2(data=z, k=k+1)
#    block1 <- NN$nn.idx[1:n1,-1]
#    block2 <- NN$nn.idx[(n1+1):n,-1]
#    i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
#    (i1 + i2) / (k * n)
#  }
#  
#  
#  eqdist.nn <- function(z,sizes,k){
#  boot.obj <- boot(data=z,statistic=Tn,R=R,
#  sim = "permutation", sizes = sizes,k=k)
#  ts <- c(boot.obj$t0,boot.obj$t)
#  p.value <- mean(ts>=ts[1])
#  list(statistic=ts[1],p.value=p.value)
#  }
#  
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#     x <- matrix(rnorm(n1*p),ncol=p);
#     y <- cbind(rnorm(n2,0,1.5),rnorm(n2,0,2));
#     z <- rbind(x,y)
#     p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#     p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#     p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999)$p.value
#  }
#  alpha <- 0.1;
#  pow <- colMeans(p.values<alpha)
#  pow
#  

## ----eval=FALSE---------------------------------------------------------------
#  
#  #mu1 = c(0,0);mu2 =c(0,0.5)
#  for(i in 1:m){
#     x <- matrix(rnorm(n1*p),ncol=p);
#     y <- cbind(rnorm(n2,0,1.5),rnorm(n2,0.5,2));
#     z <- rbind(x,y)
#     p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#     p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#     p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999)$p.value
#  }
#  
#  alpha <- 0.1;
#  pow <- colMeans(p.values<alpha)
#  pow
#  

## ----eval=FALSE---------------------------------------------------------------
#  y<-matrix(0,n2,p)
#  for(i in 1:m){
#     x <- matrix(rt(n1*p,1),ncol=p);
#  
#     for (j in 1:n2) {
#      X1 <- rnorm(2)
#      X2 <- rnorm(2,0.5,1)
#      r <- sample(c(0,1),1,replace=TRUE,prob = c(1-0.3,0.3))
#      Z <- r*X1+(1-r)*X2
#      y[j,] <- Z
#     }
#  
#     z <- rbind(x,y)
#     p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#     p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#     p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999)$p.value
#  }
#  
#  alpha <- 0.1;
#  pow <- colMeans(p.values<alpha)
#  pow

## ----eval=FALSE---------------------------------------------------------------
#  n1<-10
#  n2<-100
#  n <- n1+n2;
#  N = c(n1,n2)
#  for(i in 1:m){
#     x <- matrix(rnorm(n1*p),ncol=p);
#     y <- cbind(rnorm(n2,0,1.5),rnorm(n2,0,2));
#     z <- rbind(x,y)
#     p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#     p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#     p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999)$p.value
#  }
#  alpha <- 0.1;
#  pow <- colMeans(p.values<alpha)
#  pow

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(100)
#  rw.Metropolis <- function(n, sigma, x0, N) {
#  x <- numeric(N)
#  x[1] <- x0
#  u <- runif(N)
#  k <- 0
#  for (i in 2:N) {
#  y <- rnorm(1, x[i-1], sigma)
#  if (u[i] <= (dt(y, n) / dt(x[i-1], n)))
#  x[i] <- y else {
#  x[i] <- x[i-1]
#  k <- k + 1
#  }
#  }
#  return(list(x=x, k=k))
#  }
#  
#  
#  n <- 1 #degrees of freedom for target Student t dist.
#  N <- 2000
#  sigma <-2.5
#  x0 <- 20
#  
#  rw<- rw.Metropolis(n, sigma, x0, N)
#  
#  #number of candidate points rejected
#  print(rw$k/N)
#  
#  refline<-qt(c(0.025,0.975),df=n)
#  plot(rw$x,type = "l",xlab = bquote(sigma==.(round(sigma,3))),
#         ylab = "X",ylim = range(rw$x))
#  abline(h=refline)
#  
#  

## ----eval=FALSE---------------------------------------------------------------
#  
#  b <- 1001 #discard the burnin sample
#  y <- rw$x[b:N]
#  
#  a <-seq(.1, .9, .1)
#  Q <- qt(a, n) #theoretical quantiles
#  
#  Qrw <- quantile(y, a)
#  print(round(cbind(Q, Qrw), 3)) #not shown
#  
#  knitr::kable(round(cbind(Q, Qrw), 3))
#  

## ----eval=FALSE---------------------------------------------------------------
#  
#  #initialize constants and parameters
#  N <- 5000 #length of chain
#  burn <- 1000 #burn-in length
#  X <- matrix(0, N, 2) #the chain, a bivariate sample
#  a<-2;b<-3;n<-4#the paramater
#  
#  X[1, ] <- c(0.2, 0.3) #initialize
#  for (i in 2:N) {
#    y <- X[i-1, 2]
#    X[i, 1] <- rbinom(1,n,y)
#    x <- X[i, 1]
#    X[i, 2] <- rbeta(1,x+a,n-x+b)
#  }
#  
#  b <- burn + 1
#  f <- X[b:N, ]
#  
#  plot(f, main="", cex=.5, xlab=bquote(x),
#  ylab=bquote(y), ylim=range(X[,2]))
#  
#  

## ----eval=FALSE---------------------------------------------------------------
#  Gelman.Rubin <- function(psi) {
#  # psi[i,j] is the statistic psi(X[i,1:j])
#  # for chain in i-th row of X
#  psi <- as.matrix(psi)
#  n <- ncol(psi)
#  k <- nrow(psi)
#  psi.means <- rowMeans(psi) #row means
#  B <- n * var(psi.means) #between variance est.
#  psi.w <- apply(psi, 1, "var") #within variances
#  W <- mean(psi.w) #within est.
#  v.hat <- W*(n-1)/n + (B/n) #upper variance est.
#  r.hat <- v.hat / W #G-R statistic
#  return(r.hat)
#  }
#  
#  #9.1
#  normal.chain <- function(sigma, N, X1) {
#  #generates a Metropolis chain for t(1)
#  #with Normal(X[t], sigma^2) proposal distribution
#  #and starting value X1
#    x <- numeric(N)
#    x[1] <- X1
#    u <- runif(N)
#    k <- 0
#    for (i in 2:N) {
#    y <- rnorm(1, x[i-1], sigma)
#    if (u[i] <= (dt(y, n) / dt(x[i-1], n))){
#      x[i] <- y
#    }
#    else {
#      x[i] <- x[i-1]
#      k<- k + 1 }
#    }
#  return(x)
#  }
#  
#  sigma <- 2.5 #parameter of proposal distribution
#  k <- 4 #number of chains to generate
#  n <- 15000 #length of chains
#  b <- 1000 #burn-in length
#  #choose overdispersed initial values
#  x0 <- c(-10, -5, 5, 10)
#  #generate the chains
#  X <- matrix(0, nrow=k, ncol=n)
#  for (i in 1:k){
#    X[i, ] <- normal.chain(sigma, n, x0[i])}
#  #trace plot
#  plot(1:n,X[1,],type = "l")
#  lines(1:n,X[2,],type = "l",col=2)
#  lines(1:n,X[3,],type = "l",col=3)
#  lines(1:n,X[4,],type = "l",col=4)
#  

## ----eval=FALSE---------------------------------------------------------------
#  psi <- t(apply(X, 1, cumsum))
#  for (i in 1:nrow(psi))
#  psi[i,] <- psi[i,] / (1:ncol(psi))
#  print(Gelman.Rubin(psi))
#  
#  #plot psi for the four chains
#  par(mfrow=c(2,2))
#  for (i in 1:k)
#  plot(psi[i, (b+1):n], type="l",
#  xlab=i, ylab=bquote(psi))
#  par(mfrow=c(1,1)) #restore default
#  
#  #plot the sequence of R-hat statistics
#  rhat <- rep(0, n)
#  for (j in (b+1):n)
#      rhat[j] <- Gelman.Rubin(psi[,1:j])
#  plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
#  abline(h=1.1, lty=2)
#  
#  

## -----------------------------------------------------------------------------
f<-function(k,a,d)
{
  m1<-exp(-log(factorial(k))-k*log(2))
  m2<-exp((2*k+2)*log(sqrt(sum(a^2)))-log(2*k+1)-log(2*k+2))
  m3<-exp(lgamma((d+1)/2)+lgamma(k+3/2)-lgamma(k+d/2+1))
  ((-1)^k)*m1*m2*m3
}



## -----------------------------------------------------------------------------
g<-function(n){
  sum<-0
  for (i in 0:n) {
    sum<-sum+f(i,a,d)
  }
  sum
}




## -----------------------------------------------------------------------------
#set d=2#
a<-c(1,2)
d<-2
result<-matrix(nrow = 2,ncol = 7)
n<-round(c(1,5,10,50,100,200,300),0)
result[1,]<-n
result[2,]<-c(g(1),g(5),g(10),g(50),g(100),g(200),g(300))

dimnames(result)[[1]]<-c("n","value")

knitr::kable(result)


## -----------------------------------------------------------------------------

g <- function(a,k) {
pt(sqrt((a^2*k)/(k+1-a^2)),k,lower.tail = FALSE)
}

Z<-matrix(0,nrow = 25,ncol=3)
j<-0

for (i in c(4:25,100,500,1000)){
  h<-function(a){
  g(a,i)-g(a,i-1)
  }
  res <- uniroot(h,lower =0.01,upper = sqrt(i)-0.01)
  j<-j+1
  Z[j,1]<-res$root
  Z[j,2]<-res$f.root
  Z[j,3]<-res$iter
}


result<-matrix(nrow = 25,ncol = 4)
result[,1]<-c(4:25,100,500,1000)
result[,2:4]<-Z

dimnames(result)[[2]]<-c("k","root","f.root","iter")

knitr::kable(result)


## -----------------------------------------------------------------------------

g <- function(a,k) {
  f<-function(u){
    (1+u^2/(k-1))^(-k/2)
  }

  c<-sqrt(((a^2)*k)/(k+1-a^2))
  d<-exp(log(2)+lgamma((k+1)/2)-(1/2)*log(pi*(k))-lgamma(k/2))
  it<-integrate(f,lower = 0,upper = c)
  d*it$value
}

Z<-matrix(0,nrow = 25,ncol=3)
j<-0

for (i in c(4:25,100,500,1000)){
  h<-function(a){
  g(a,i)-g(a,i-1)
  }
  res <- uniroot(h,lower =0.01,upper = sqrt(i)-0.01)
  j<-j+1
  Z[j,1]<-res$root
  Z[j,2]<-res$f.root
  Z[j,3]<-res$iter
}


result<-matrix(nrow = 25,ncol = 4)
result[,1]<-c(4:25,100,500,1000)
result[,2:4]<-Z

dimnames(result)[[2]]<-c("k","root","f.root","iter")

knitr::kable(result)


## ----eval=FALSE---------------------------------------------------------------
#  #compute MLE of lambda
#  y<-c(0.54,0.48, 0.33,0.43,0.91,0.21,0.85,1.00,1.00,1.00)
#  mlogL <- function(lambda=1) {
#  # minus log-likelihood
#  return(-(length(y)*log(lambda)-lambda*sum(y)))
#  }
#  
#  library(stats4)
#  fit <- mle(mlogL)
#  as.numeric(c(1/mean(y),fit@coef,sqrt(fit@vcov)))
#  
#  library(imager)
#  imgjpg<-load.image("E:/研一上/统计计算/作业/11-18.jpg")
#  plot(imgjpg,axes=FALSE)

## -----------------------------------------------------------------------------
trims <- c(0, 0.1, 0.2, 0.5)
x <- rcauchy(100)
lapply(trims, function(trim) mean(x, trim = trim))
lapply(trims, mean, x = x)

## -----------------------------------------------------------------------------
#11.1.2.3
formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)
#lapply()
result<-lapply(formulas, function(x) lm(x,data = mtcars))
unlist(lapply(result,function(mod) summary(mod)$r.squared))

## -----------------------------------------------------------------------------
#11.1.2.4
set.seed(123)
bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
mtcars[rows, ]
})

result<-lapply(bootstraps, function(x) lm(mpg ~ disp,data = x))
unlist(lapply(result,function(mod) summary(mod)$r.squared))

## -----------------------------------------------------------------------------
#Use vapply() to: a) Compute the standard deviation of every column in a numeric data frame.
x<-rnorm(10,0,1)
y<-rnorm(10,0,2)
z<-rnorm(10,0,3)
d<-data.frame(x,y,z)

vapply(d, sd, FUN.VALUE = numeric(1))


## -----------------------------------------------------------------------------
# Compute the standard deviation of every numeric column in a mixed data frame. (Hint: you’ll need to use vapply() twice.)
d <- data.frame(x = 1:10,y=21:30, z = letters[1:10])

a<-vapply(d, function(x) is.numeric(x) ,FUN.VALUE = logical(1) )
vapply(d[,a], function(x) sd(x), FUN.VALUE = numeric(1))


## -----------------------------------------------------------------------------
set.seed(1)
library(Rcpp)
library(microbenchmark)
#write Rcpp function
cppFunction(' NumericMatrix gibbsC(int N,int n,int a,int b) {
      NumericMatrix X(N,2);
     X(0,0)=0.2;
     X(0,1)=0.3;
     for(int i = 1; i < N; i++) {
      double y = X(i-1, 1);
      X(i, 0) = rbinom(1,n,y)[0];
      double x = X(i, 0);
      X(i, 1) = rbeta(1,x+a,n-x+b)[0];
  }
      return X;
    }')
    
#R function
gibbsR<-function(N,n,a,b){
  X <- matrix(0, N, 2)
  X[1, ] <- c(0.2, 0.3) #initialize
  for (i in 2:N) {
   y <- X[i-1, 2]
   X[i, 1] <- rbinom(1,n,y)
   x <- X[i, 1]
   X[i, 2] <- rbeta(1,x+a,n-x+b)
  }
  X
}

N<-1000
a<-2;b<-3;n<-4

gibbsr=gibbsR(N,n,a,b)
gibbsc=gibbsC(N,n,a,b)

qqplot(gibbsr[,1],gibbsc[,1])
qqplot(gibbsr[,2],gibbsc[,2])

ts <- microbenchmark(gibbsr=gibbsR(N,n,a,b),
gibbsc=gibbsC(N,n,a,b))

summary(ts)[,c(1,3,5,6)]


