## 1. looking at bartlett rates with discrete data

## 1a. simple logistic regression using Jensen '93 formula. seems that
## even with discrete predictor the 1/n^2 rate can be achieved. rate
## seems to improve when x has more support points. does not seem to
## hold when conditioning on the covariate.

## set.seed(1)
## n <- 3e1
beta <- c(1,0)
ns <- 5:15
start <- Sys.time()
## gaps <- sapply(ns, function(n) {
gaps <- parallel::mclapply(ns, mc.cores=4, FUN=function(n) {
    cat('.')
    ## x <- rpois(n,10)
    stats <- replicate(n^4, {
        ## x <- rnorm(n)
        ## x <- rpois(n,lambda=2)
        x <- rbinom(n,5,1/3)
        p <- plogis(beta[1]+beta[2]*x)
        y <- rbinom(n,1,p)
        ## glm0 <- glm(y~1+x,family=binomial(link=logit))
        l <- function(beta,x,y)beta[2]*sum(x*y) + beta[1]*sum(y) - sum(log(1+exp(beta[1]+beta[2]*x)))
        mle.full <- optim(par=c(0,0),fn=function(beta)-l(beta,x,y))$par
        mle.redu <- optimize(function(beta1)-l(c(beta1,0),x,y),interval=c(-10,10))$min
        W <- 2*(l(mle.full,x,y) - l(c(mle.redu,0),x,y))
        s <- sapply(1:4,function(j)mean(x^j))
        bartlett <- function(beta)1/4*(4-exp(beta[1])-exp(-beta[1]))*(1/(s[2]-s[1]^2)^2*(3*s[2]^2-4*s[1]*s[3]+s[4])-1) +  5/12*(exp(beta[1])+exp(-beta[1])-2)*(1/(s[2]-s[1]^2)^3*(4*s[2]^3-3*s[1]^2*s[2]^2+4*s[1]^3*s[3]-6*s[1]*s[2]*s[3]+s[3]^2)-1)
        ## W.bartlett <- W / (1+bartlett(mle.full)/n)
        W.bartlett <- W / (1+bartlett(c(mle.redu,0))/n)
        c(W=W,W.bartlett=W.bartlett)
    })
    ## plot(ecdf(stats))
    ## curve(pchisq(x,df=1),add=TRUE,col=2)
    ## ecdf.W <- ecdf(pchisq(stats,df=1,lower=FALSE))
    ## ks.W <- max(abs(ecdf.W(knots(ecdf.W)) - knots(ecdf.W)))
    ecdfs <- apply(stats,1,function(x)ecdf(pchisq(x,df=1,lower=FALSE)))
    ks <- sapply(ecdfs,function(ecdf)max(abs(ecdf(knots(ecdf)) - knots(ecdf))))
})
print(Sys.time()-start)
gaps <- simplify2array(gaps)
matplot(ns,t(gaps))
apply(gaps,1,function(x)coef(lm(log(x)~log(ns)))[2])
## plot(ns,gaps)
## coef(lm(log(gaps)~log(ns)))[2]



## 1b. fuller sim to see effect of # of x support points affects rate
ms <- 1:5
start <- Sys.time()
by.m <- sapply(ms, function(m) {
    cat('.')
    beta <- c(1,0)
    ns <- 5:15
    ## gaps <- sapply(ns, function(n) {
    replicate(5, {
        gaps <- parallel::mclapply(ns, mc.cores=4, FUN=function(n) {
            ## x <- rpois(n,10)
            stats <- replicate(n^4, {
                ## x <- rnorm(n)
                ## x <- rpois(n,lambda=2)
                x <- rbinom(n,m,1/3)
                p <- plogis(beta[1]+beta[2]*x)
                y <- rbinom(n,1,p)
                ## glm0 <- glm(y~1+x,family=binomial(link=logit))
                l <- function(beta,x,y)beta[2]*sum(x*y) + beta[1]*sum(y) - sum(log(1+exp(beta[1]+beta[2]*x)))
                mle.full <- optim(par=c(0,0),fn=function(beta)-l(beta,x,y))$par
                mle.redu <- optimize(function(beta1)-l(c(beta1,0),x,y),interval=c(-10,10))$min
                W <- 2*(l(mle.full,x,y) - l(c(mle.redu,0),x,y))
                s <- sapply(1:4,function(j)mean(x^j))
                bartlett <- function(beta)1/4*(4-exp(beta[1])-exp(-beta[1]))*(1/(s[2]-s[1]^2)^2*(3*s[2]^2-4*s[1]*s[3]+s[4])-1) +  5/12*(exp(beta[1])+exp(-beta[1])-2)*(1/(s[2]-s[1]^2)^3*(4*s[2]^3-3*s[1]^2*s[2]^2+4*s[1]^3*s[3]-6*s[1]*s[2]*s[3]+s[3]^2)-1)
                ## W.bartlett <- W / (1+bartlett(mle.full)/n)
                W.bartlett <- W / (1+bartlett(c(mle.redu,beta[2]))/n)
                c(W=W,W.bartlett=W.bartlett)
            })
            ecdfs <- apply(stats,1,function(x)ecdf(pchisq(x,df=1,lower=FALSE)))
            ks <- sapply(ecdfs,function(ecdf)max(abs(ecdf(knots(ecdf)) - knots(ecdf))))
        })
        gaps <- simplify2array(gaps)
        apply(gaps,1,function(x)coef(lm(log(x)~log(ns)))[2])['W.bartlett']
    })
})
print(Sys.time()-start)
apply(by.m,2,median)



## mma check
x <- c(.5,1.5)
s <- sapply(1:4,function(j)mean(x^j))
bartlett <- function(beta)1/4*(4-exp(beta[1])-exp(-beta[1]))*(1/(s[2]-s[1]^2)^2*(3*s[2]^2-4*s[1]*s[3]+s[4])-1) +  5/12*(exp(beta[1])+exp(-beta[1])-2)*(1/(s[2]-s[1]^2)^3*(4*s[2]^3-3*s[1]^2*s[2]^2+4*s[1]^3*s[3]-6*s[1]*s[2]*s[3]+s[3]^2)-1)
1+bartlett(c(2,NA))/2


## 2. laplace approx

## 2a. cox & barndorff-nielson

## 2a-1. normal approx--O(1/z^5) error as expected
z <- seq(0,3,len=20)
true <- 1-pnorm(z)
try <- dnorm(z)*(1/z-1/z^3)
plot(true-try)
coef(lm(log(abs(true-try))~log(z)))[2]


## 2a-2. gamma function--relative error 1/z^-1.5, as expected
z <- seq(1,10,len=20)
## z <- c(.5,1,1.5,2,3)
true <- gamma(z+1)
## try <- z^(z+1) * exp(-z)*sqrt(2*pi)/sqrt(z)
## try <- sqrt(2*pi)*z^(z+1/2)*exp(-z*(1+1/12/z))
try <- sqrt(2*pi)*z^(z+1/2)*exp(-z)*(1+1/12/z)
## true
## try
plot(true-try)
coef(lm(log(abs((true-try)/true))~log(z)))[2]


## 2b. normal-binomial likelihood
mu <- 0
tau2 <- 2
n <- 30
pairs <- replicate(1e2, {
    theta <- rnorm(1,sd=sqrt(tau2))
    e <- rbinom(1,n,plogis(mu+theta))
    g <- function(theta)dnorm(theta/sqrt(tau2))
    d1g <- function(theta)-theta/tau2*dnorm(theta/sqrt(tau2))
    d2g <- function(theta)-1/tau2*dnorm(theta/sqrt(tau2)) + dnorm(theta/sqrt(tau2))*(theta/tau2)^2
    ## d2g <- function(theta)-dnorm(theta)+theta^2*dnorm(theta)
    r <- function(theta)-e/n*log(plogis(mu+theta))-(1-e/n)*log(1-plogis(mu+theta))
    ## curve(r,-4,4)
    theta0 <- optimize(f=r,c(-5,5))$min
    dr <- function(theta)plogis(mu+theta)-e/n
    d2r <- function(theta)plogis(mu+theta)*(1-plogis(mu+theta))
    d3r <- function(theta)plogis(mu+theta)*(1-3*plogis(mu+theta)+2*plogis(mu+theta)^2)
    d4r <- function(theta)plogis(mu+theta)-7*plogis(mu+theta)^2+12*plogis(mu+theta)^3-6*plogis(mu+theta)^4
    ## dr(theta0)
    ## d2r(theta0)
    try1 <- exp(-n*r(theta0))*g(theta0)*sqrt(2*pi)/sqrt(n)/sqrt(d2r(theta0))
    try2 <- try1/g(theta0) * (g(theta0)+1/n*(   d2g(theta0)/2/d2r(theta0) - d3r(theta0)*d1g(theta0)/2/d2r(theta0)^2 - d4r(theta0)*g(theta0)/8/d2r(theta0)^2 + 5*d3r(theta0)^2*g(theta0)/12/d2r(theta0)^3   )  )
    true <- integrate(function(theta)plogis(mu+theta)^e*(1-plogis(mu+theta))^(n-e)*dnorm(theta/sqrt(tau2)),-20,20)$val
    log(c(try1=try1,try2=try2,true=true))
})
## plot(pairs[1,],pairs[2,],pairs[3,])
plot(pairs['true',],pairs['try1',])
points(pairs['true',],pairs['try2',],col=2)
abline(0,1)


## 2c. compute MLE's
set.seed(1)
n <- 20
mu.null <- 1
tau2.null <- 2
m <- rpois(n,lambda=20)
theta <- rnorm(n,sd=sqrt(tau2.null))
e <- rbinom(n,m,plogis(mu.null+theta))
l.laplace <- function(mu,tau2,e,m) {
    g <- function(theta)dnorm(theta/sqrt(tau2))/sqrt(tau2)
    d1g <- function(theta)-theta/tau2*dnorm(theta/sqrt(tau2))/sqrt(tau2)
    d2g <- function(theta)(-1/tau2*dnorm(theta/sqrt(tau2)) + dnorm(theta/sqrt(tau2))*(theta/tau2)^2)/sqrt(tau2)
    l.units.laplace <- sapply(1:n, function(i) {
        r <- function(theta) -e[i]/m[i]*log(plogis(mu+theta))-(1-e[i]/m[i])*log(1-plogis(mu+theta))
        dr <- function(theta)plogis(mu+theta)-e[i]/m[i]
        d2r <- function(theta)plogis(mu+theta)*(1-plogis(mu+theta))
        d3r <- function(theta)plogis(mu+theta)*(1-3*plogis(mu+theta)+2*plogis(mu+theta)^2)
        d4r <- function(theta)plogis(mu+theta)-7*plogis(mu+theta)^2+12*plogis(mu+theta)^3-6*plogis(mu+theta)^4
        theta.min <- optimize(f=r,c(-5,5))$min
        ## theta.min <- qlogis(e[i]/m[i])-mu # breaks on zero cells
        try1 <- exp(-m[i]*r(theta.min))*g(theta.min)*sqrt(2*pi)/sqrt(m[i])/sqrt(d2r(theta.min))
        try2 <- try1/g(theta.min) * (g(theta.min)+1/m[i]*(   d2g(theta.min)/2/d2r(theta.min) - d3r(theta.min)*d1g(theta.min)/2/d2r(theta.min)^2 - d4r(theta.min)*g(theta.min)/8/d2r(theta.min)^2 + 5*d3r(theta.min)^2*g(theta.min)/12/d2r(theta.min)^3   )  )
        log(try2)
    })
    return(sum(l.units.laplace))
}
## l.laplace(mu,tau2,e,m)
op <- par(mfrow=c(1,2))
mus <- seq(-5,5,len=1e2)
ls <- sapply(mus,function(mu)l.laplace(mu,tau2.null,e,m))
mu.star <- optimize(f=function(mu)l.laplace(mu,tau2.null,e,m),interval=c(-5,5),maximum=TRUE)$max
plot(mus,ls)
abline(v=mu.null)
abline(v=mu.star,col=2)
tau2s <- seq(0.3,10,len=1e2)
ls <- sapply(tau2s,function(tau2)l.laplace(mu.null,tau2,e,m))
tau2.star <- optimize(f=function(tau2)l.laplace(mu.null,tau2,e,m),interval=c(.1,20),maximum=TRUE)$max
plot(tau2s,ls)
abline(v=tau2.null)
abline(v=tau2.star,col=2)
par(op)
optim(par=c(0,1),fn=function(par)-l.laplace(par[1],par[2],e,m),control=list(trace=0))$par
metafor::rma.glmm(measure='PLO',xi=e,ni=m,data=data.frame(study = 1:n, e=e, m=m))



## 2d. reparametrize in terms of standard normal theta
## normal-binomial likelihood
mu <- 0
tau2 <- 10
n <- 30
pairs <- replicate(1e2, {
    theta <- rnorm(1)
    e <- rbinom(1,n,plogis(mu+sqrt(tau2)*theta))
    g <- dnorm
    d1g <- function(theta)-theta*dnorm(theta)
    d2g <- function(theta)-dnorm(theta) +theta^2*dnorm(theta)
    r <- function(theta)-e/n*log(plogis(mu+sqrt(tau2)*theta))-(1-e/n)*log(1-plogis(mu+sqrt(tau2)*theta))
    ## curve(r,-4,4)
    theta0 <- optimize(f=r,c(-5,5))$min
    dr <- function(theta)sqrt(tau2)*( plogis(mu+sqrt(tau2)*theta)-e/n )
    d2r <- function(theta)tau2*plogis(mu+sqrt(tau2)*theta)*(1-plogis(mu+sqrt(tau2)*theta))
    d3r <- function(theta)tau2^(3/2)*plogis(mu+sqrt(tau2)*theta)*(1-3*plogis(mu+sqrt(tau2)*theta)+2*plogis(mu+sqrt(tau2)*theta)^2)
    d4r <- function(theta)tau2^2*plogis(mu+sqrt(tau2)*theta)*(1-7*plogis(mu+sqrt(tau2)*theta)+12*plogis(mu+sqrt(tau2)*theta)^2-6*plogis(mu+sqrt(tau2)*theta)^3)
    ## dr(theta0)
    ## d2r(theta0)
    try1 <- exp(-n*r(theta0))*g(theta0)*sqrt(2*pi)/sqrt(n)/sqrt(d2r(theta0))
    try2 <- try1/g(theta0) * (g(theta0)+1/n*(   d2g(theta0)/2/d2r(theta0) - d3r(theta0)*d1g(theta0)/2/d2r(theta0)^2 - d4r(theta0)*g(theta0)/8/d2r(theta0)^2 + 5*d3r(theta0)^2*g(theta0)/12/d2r(theta0)^3   )  )
    true <- integrate(function(theta)plogis(mu+theta)^e*(1-plogis(mu+theta))^(n-e)*dnorm(theta/sqrt(tau2))/sqrt(tau2),-20,20)$val
    ##     true <- true <- integrate(
    ##     function(theta)
    ##         plogis(mu + sqrt(tau2) * theta)^e *
    ##         (1 - plogis(mu + sqrt(tau2) * theta))^(n - e) *
    ##         dnorm(theta),
    ##     -20, 20
    ## )$val
    log(c(try1=try1,try2=try2,true=true))
})
## plot(pairs[1,],pairs[2,],pairs[3,])
plot(pairs['true',],pairs['try1',])
points(pairs['true',],pairs['try2',],col=2)
abline(0,1)


## compute MLE's
## set.seed(1)
n <- 20
mu.null <- 1
tau2.null <- 2
m <- rpois(n,lambda=20)
theta <- rnorm(n)
e <- rbinom(n,m,plogis(mu.null+sqrt(tau2.null)*theta))
l.laplace <- function(mu,tau2,e,m) {
    n <- length(m)
    g <- dnorm
    d1g <- function(theta)-theta*dnorm(theta)
    d2g <- function(theta)-dnorm(theta) +theta^2*dnorm(theta)
    l.units.laplace <- sapply(1:n, function(i) {
        ## browser()
        r <- function(theta)-e[i]/m[i]*log(plogis(mu+sqrt(tau2)*theta))-(1-e[i]/m[i])*log(1-plogis(mu+sqrt(tau2)*theta))
        ## curve(r,-4,4)
        theta0 <- optimize(f=r,c(-5,5))$min
        dr <- function(theta)sqrt(tau2)*( plogis(mu+sqrt(tau2)*theta)-e[i]/m[i] )
        d2r <- function(theta)tau2*plogis(mu+sqrt(tau2)*theta)*(1-plogis(mu+sqrt(tau2)*theta))
        d3r <- function(theta)tau2^(3/2)*plogis(mu+sqrt(tau2)*theta)*(1-3*plogis(mu+sqrt(tau2)*theta)+2*plogis(mu+sqrt(tau2)*theta)^2)
        d4r <- function(theta)tau2^2*plogis(mu+sqrt(tau2)*theta)*(1-7*plogis(mu+sqrt(tau2)*theta)+12*plogis(mu+sqrt(tau2)*theta)^2-6*plogis(mu+sqrt(tau2)*theta)^3)
        theta.min <- optimize(f=r,c(-5,5))$min
        ## theta.min <- (qlogis(e[i]/m[i])-mu)/sqrt(tau2) # breaks on zero cells
        try1 <- exp(-m[i]*r(theta0))*g(theta0)*sqrt(2*pi)/sqrt(m[i])/sqrt(d2r(theta0))
        try2 <- try1/g(theta0) * (g(theta0)+1/m[i]*(   d2g(theta0)/2/d2r(theta0) - d3r(theta0)*d1g(theta0)/2/d2r(theta0)^2 - d4r(theta0)*g(theta0)/8/d2r(theta0)^2 + 5*d3r(theta0)^2*g(theta0)/12/d2r(theta0)^3   )  )
        log(try2)
    })
    return(sum(l.units.laplace))
}
## l.laplace(mu,tau2,e,m)
op <- par(mfrow=c(1,2))
mus <- seq(-5,5,len=1e2)
ls <- sapply(mus,function(mu)l.laplace(mu,tau2.null,e,m))
mu.star <- optimize(f=function(mu)l.laplace(mu,tau2.null,e,m),interval=c(-5,5),maximum=TRUE)$max
plot(mus,ls)
abline(v=mu.null)
abline(v=mu.star,col=2)
tau2s <- seq(0.3,10,len=1e2)
ls <- sapply(tau2s,function(tau2)l.laplace(mu.null,tau2,e,m))
tau2.star <- optimize(f=function(tau2)l.laplace(mu.null,tau2,e,m),interval=c(.1,20),maximum=TRUE)$max
plot(tau2s,ls)
abline(v=tau2.null)
abline(v=tau2.star,col=2)
par(op)
optim(par=c(0,1),fn=function(par)-l.laplace(par[1],par[2],e,m),control=list(trace=0))$par
metafor::rma.glmm(measure='PLO',xi=e,ni=m,data=data.frame(study = 1:n, e=e, m=m))


## 2e. read in a likelihood formula exported from mma
source('utils.R')
n <- 20
mu.null <- 1
tau2.null <- 2
m <- rpois(n,lambda=20)
theta <- rnorm(n)
e <- rbinom(n,m,plogis(mu.null+sqrt(tau2.null)*theta))
rm(theta)
beta <- c(tau2.null,mu.null)
mma.str <- paste(readLines('mma-formulas/try.txt',warn=FALSE),collapse='')
mma.str
mma.str <- mma.to.R(
    mma.str,
    funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial"),
    rename = c(
        Sqrt = "sqrt",
        Log = "log",
        LogisticSigmoid = "plogis",
        Binomial = "choose"
    )
)
## j <- 2
## k <- 3
E <- exp(1)
## theta <- 1
## l.term.str <- eval(parse(text=mma.str))
## mma.str <- gsub(pattern=c('m\\[j\\]'),replacement=c('m'),mma.str)
l.term.fun <- make.fun(
    gsub(pattern=c('m\\[j\\]'),replacement=c('m'),mma.str),
    args = c('theta','m','k')
)
sapply(1:n, function(j) 
    sapply(1:(m[j]-1), function(k) {
        ## browser()
        integrate(function(theta)l.term.fun(theta,m[j],k),-Inf,Inf)$val
    }) |> sum()
    ) |> sum()



## read in multiple mma strings, build information matrix
source('utils.R')
n <- 20
mu.null <- 1
tau2.null <- 2
m <- rpois(n,lambda=20)
theta <- rnorm(n)
e <- rbinom(n,m,plogis(mu.null+sqrt(tau2.null)*theta))
rm(theta)
beta <- c(tau2.null,mu.null)
## idx.str <- c('11')
idx.strings <- apply(expand.grid(1:2,1:2),1,paste,collapse='')
k <- sapply(idx.strings, function(idx.str) {
    mma.str <- paste(readLines(paste0('mma-formulas/k-unit-',idx.str,'.txt'),warn=FALSE),collapse='')
    ## mma.str
    mma.str <- mma.to.R(mma.str,
                        funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial"),
                        rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose")
                        )
    E <- exp(1)
    l.deriv <- make.fun(
        gsub(pattern=c('m\\[j\\]'),replacement=c('m'),mma.str),
        args = c('theta','m','k')
    )
    sapply(1:n, function(j) 
        sapply(1:(m[j]-1), function(k) {
            ## browser()
            integrate(function(theta)l.deriv(theta,m[j],k)*dnorm(theta),-Inf,Inf)$val
        }) |> sum()
        ) |> sum()
})
inf.mat <- matrix(k,nrow=2)
inf.mat
se <- diag(sqrt(-solve(inf.mat)))





## 2f. encapsulate, check CIs. match metafor's pretty well.
source('utils.R')
mle.laplace <- function(e,m) {
    l.laplace <- function(mu,tau2,e,m) {
        n <- length(m)
        g <- dnorm
        d1g <- function(theta)-theta*dnorm(theta)
        d2g <- function(theta)-dnorm(theta) +theta^2*dnorm(theta)
        l.units.laplace <- sapply(1:n, function(i) {
            ## browser()
            r <- function(theta)-e[i]/m[i]*log(plogis(mu+sqrt(tau2)*theta))-(1-e[i]/m[i])*log(1-plogis(mu+sqrt(tau2)*theta))
            ## curve(r,-4,4)
            theta0 <- optimize(f=r,c(-5,5))$min
            dr <- function(theta)sqrt(tau2)*( plogis(mu+sqrt(tau2)*theta)-e[i]/m[i] )
            d2r <- function(theta)tau2*plogis(mu+sqrt(tau2)*theta)*(1-plogis(mu+sqrt(tau2)*theta))
            d3r <- function(theta)tau2^(3/2)*plogis(mu+sqrt(tau2)*theta)*(1-3*plogis(mu+sqrt(tau2)*theta)+2*plogis(mu+sqrt(tau2)*theta)^2)
            d4r <- function(theta)tau2^2*plogis(mu+sqrt(tau2)*theta)*(1-7*plogis(mu+sqrt(tau2)*theta)+12*plogis(mu+sqrt(tau2)*theta)^2-6*plogis(mu+sqrt(tau2)*theta)^3)
            theta.min <- optimize(f=r,c(-5,5))$min
            ## theta.min <- (qlogis(e[i]/m[i])-mu)/sqrt(tau2) # breaks on zero cells
            try1 <- exp(-m[i]*r(theta0))*g(theta0)*sqrt(2*pi)/sqrt(m[i])/sqrt(d2r(theta0))
            try2 <- try1/g(theta0) * (g(theta0)+1/m[i]*(   d2g(theta0)/2/d2r(theta0) - d3r(theta0)*d1g(theta0)/2/d2r(theta0)^2 - d4r(theta0)*g(theta0)/8/d2r(theta0)^2 + 5*d3r(theta0)^2*g(theta0)/12/d2r(theta0)^3   )  )
            log(try2)
        })
        return(sum(l.units.laplace))
    }
    mle <- optim(par=c(0,1),fn=function(par)-l.laplace(par[1],par[2],e,m),control=list(trace=0))$par
    return(structure(mle,names=c('mu','tau2')))
}
se.laplace <- function(e,m,mle) {
    beta <- c(mle[c('tau2','mu')])
    idx.strings <- apply(expand.grid(1:2,1:2),1,paste,collapse='')
    k <- sapply(idx.strings, function(idx.str) {
        mma.str <- paste(readLines(paste0('mma-formulas/k-unit-',idx.str,'.txt'),warn=FALSE),collapse='')
        ## mma.str
        mma.str <- mma.to.R(mma.str,
                            funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial"),
                            rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose")
                            )
        E <- exp(1)
        l.deriv <- make.fun(
            gsub(pattern=c('m\\[j\\]'),replacement=c('m'),mma.str),
            args = c('theta','m','k')
        )
        sapply(1:n, function(j) 
            sapply(1:(m[j]-1), function(k) {
                ## browser()
                integrate(function(theta)l.deriv(theta,m[j],k)*dnorm(theta),-Inf,Inf)$val
            }) |> sum()
            ) |> sum()
    })
    inf.mat <- matrix(k,nrow=2)
    inf.mat
    se <- sqrt(-diag(solve(inf.mat)))
    return(se)
}
## mle <- mle.laplace(e,m)
n <- 20
mu.null <- 1
tau2.null <- 2
for(i in 1:10) {
    m <- rpois(n,lambda=20)
    theta <- rnorm(n)
    e <- rbinom(n,m,plogis(mu.null+sqrt(tau2.null)*theta))
    rm(theta)
    rma0 <- metafor::rma.glmm(measure='PLO',xi=e,ni=m,data=data.frame(study = 1:n, e=e, m=m))
    mle <- structure(unname(c(rma0$beta,rma0$tau2)),names=c('mu','tau2'))
    mle
    print(c(rma0$se,se.laplace(e,m,mle)[2]))
}


## 2g. get all the full model k terms

source('utils.R')
lawley <- function(p,m,beta,l.derivs) {
    ## k.vals <- sapply(mma.strings, function(mma.str) {
    ##     l.deriv <- str.to.fun(
    ##         gsub(pattern=c('m\\[j\\]'),replacement=c('m'),mma.str),
    ##         args = c('theta','m','k'),
    ##         binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt,m=m,beta=beta)
    ##     )
    k.vals <- sapply(l.derivs, function(l.deriv) {
        ## browser()
        sapply(1:length(m), function(j) 
            sapply(1:(m[j]-1), function(k) {
                ## browser()
                integrate(function(theta)l.deriv(theta,m[j],x,beta)*dnorm(theta),-Inf,Inf)$val
            }) |> sum()
            ) |> sum()
    })
    idx.strings <- gsub('.txt','',names(k.vals),fixed=TRUE) |> sub(pattern='.*k-',replace='') |> strsplit(split='-')
    idx.strings <- lapply(idx.strings,function(x)if(length(x)==1)c(x,0)else x)
    k.df <- data.frame(inner=sapply(idx.strings,`[`,1),outer=sapply(idx.strings,`[`,2),k.val=k.vals)
    k <- function(inner,outer=0)k.df$k.val[k.df$inner==paste(inner,collapse='') & k.df$outer==paste(outer,collapse='')]
    inf <- -matrix(apply(expand.grid(1:p,1:p),1,k),nrow=p)
    inv.inf <- solve(inf)
    lawley.4 <- apply(do.call(expand.grid,rep(list(1:p),4)), 1, function(idx.str) {
        r <- idx.str[1]; s <- idx.str[2]; t <- idx.str[3]; u <- idx.str[4]
        inv.inf[r,s] * inv.inf[t,u] * (1/4*k(c(r,s,t,u)) - k(c(r,s,t),u) + k(c(r,t),c(s,u)))
    })
    lawley.6 <- apply(do.call(expand.grid,rep(list(1:p),6)), 1, function(idx.str) {
        r <- idx.str[1]; s <- idx.str[2]; t <- idx.str[3]; u <- idx.str[4]; v <- idx.str[5]; w <- idx.str[6]
        -inv.inf[r,s] * inv.inf[t,u] * inv.inf[v,w] * (1/6*k(c(r,t,v))*k(c(s,u,w)) + 1/4*k(c(r,t,u))*k(c(s,v,w)) - k(c(r,t,v))*k(c(s,w),u) - k(c(r,t,u))*k(c(s,w),v) + k(c(r,t),v)*k(c(s,w),u) + k(c(r,t),u)*k(c(s,w),v))
    })
    ## sum(l4)
    ## sum(l6)
    return( sum(lawley.4) - sum(lawley.6) )
}
read.mma <- function(mma.str.path) {
    filelist <- dir(mma.str.path)
    filelist <- paste0(mma.str.path,filelist[grep('\\.txt',filelist)])
    mma.strings <- sapply(filelist,function(file)paste(readLines(file,warn=FALSE),collapse=''))
    mma.strings <- sapply(mma.strings,
                          function(mma.str) mma.to.R(mma.str,
                                                     funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial"),
                                                     rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose")
                                                     )
                          )
    l.derivs <- sapply(mma.strings,function(mma.str) 
        str.to.fun(
            gsub(pattern=c('m\\[j\\]'),replacement=c('m'),mma.str),
            args = c('theta','m','x','beta'),
            binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt)
        ))
}
set.seed(1)
n <- 20
mu.null <- 1
tau2.null <- 2
## beta <- runif(2)#c(mle[c('tau2','mu')]) ## should be passed to make.fun
l.derivs.full <- read.mma('./mma-formulas/laplace/full/')
l.derivs.redu <- read.mma('./mma-formulas/laplace/redu/')
start <- Sys.time()
m <- rpois(n,lambda=20)
m <- rep(20,n)
theta <- rnorm(n)
e <- rbinom(n,m,plogis(mu.null+sqrt(tau2.null)*theta))
rm(theta)
rma0 <- metafor::rma.glmm(measure='PLO',xi=e,ni=m,data=data.frame(study = 1:n, e=e, m=m))
mle <- structure(unname(c(rma0$beta,rma0$tau2)),names=c('mu','tau2'))
beta <- rev(mle)
## beta <- c(tau2.null,mu.null)##!!!!!!!!!!!!!!!!!
eps.full <- lawley(p=2,m=m,beta=beta,l.derivs=l.derivs.full)
eps.redu <- lawley(p=1,m=m,beta=beta,l.derivs=l.derivs.redu)
correction <- 1 + eps.full - eps.redu
correction
print(Sys.time() - start)


## 2h. compute LR statistic using glmer. check LR statistic is asy chi
## squared.
library(lme4)
set.seed(1)
n <- 20
mu.null <- 1
tau2.null <- 2
start <- Sys.time()
m <- rpois(n,lambda=20)
LR.stats <- replicate(1e2, {
    theta <- rnorm(n)
    e <- rbinom(n,m,plogis(mu.null+sqrt(tau2.null)*theta))
    rm(theta)
    dat <- data.frame(e = e, non = m - e, study = factor(seq_along(e)),off = mu.null)
    fit.alt <- glmer(cbind(e, non) ~ 1 + (1 | study), data = dat,family = binomial(link = "logit"),nAGQ = 7)
    fit.null <- glmer(cbind(e, non) ~ 0 + offset(off) + (1 | study),data = dat, family = binomial(link = "logit"), nAGQ = 7)
    LR <- 2 * (as.numeric(logLik(fit.alt)) - as.numeric(logLik(fit.null)))
    ## pval <- pchisq(LR, df = 1, lower.tail = FALSE)
    ## c(LR = LR, p.value = pval)
})
print(Sys.time() - start)
plot(ecdf(LR.stats))
curve(pchisq(x,df=1),add=TRUE,col=2)






curve(x*(5*x+6)/2/(x+1)/(x+3),0,10)
curve(log(x+1),add=TRUE,lty=2)
curve(x-x^2+x^3/3,add=T,col=2)

curve(3*x/2*(x+6)/(x+3)^2,0,5)
curve(log(x+1),add=TRUE,lty=2)
curve(x-x^2/2,add=TRUE,col=2)

curve(27*(x+1),0,1)
curve((x+3)^3,add=T,col=2)



source('utils.R')
require(lme4)
ns <- 5:10
mu.null <- 0
tau2.null <- 2
## beta <- runif(2)#c(mle[c('tau2','mu')]) ## should be passed to make.fun
l.derivs.full <- read.mma('./mma-formulas/full/')
l.derivs.redu <- read.mma('./mma-formulas/redu/')
start <- Sys.time()
## m <- rpois(n,lambda=60)
by.n <- parallel::mclapply(ns, mc.cores=1,FUN=function(n) {
    m <- rep(30,n)
    stats <- replicate(2, {
        cat('.')
        theta <- rnorm(n)
        e <- rbinom(n,m,plogis(mu.null+sqrt(tau2.null)*theta))
        ## if(any(e==0)||any(e==m))return(c(NA,NA))
        ## rm(theta)
        dat <- data.frame(e = e, non = m - e, study = factor(seq_along(e)),off = mu.null)
        fit.alt <- glmer(cbind(e, non) ~ 1 + (1 | study), data = dat,family = binomial(link = "logit"),nAGQ = 7)
        fit.null <- glmer(cbind(e, non) ~ 0 + offset(off) + (1 | study),data = dat, family = binomial(link = "logit"), nAGQ = 7)
        LR <- 2 * (as.numeric(logLik(fit.alt)) - as.numeric(logLik(fit.null)))
        tau2.mle.redu <- VarCorr(fit.null)$study |> as.numeric() |> unname()
        beta <- c(tau2.mle.redu,mu.null)
        eps.full <- lawley(p=2,m=m,beta=beta,l.derivs=l.derivs.full)
        eps.redu <- lawley(p=1,m=m,beta=beta,l.derivs=l.derivs.redu)
        correction <- 1 + eps.full - eps.redu
        ## if(correction<0)break#browser()
        ## correction <- 1#!!!
        if(correction>10)browser()
        c(LR=LR,LR.bartlett=LR/correction,correction=correction)
    })
})





## 2i. template of sim. assumes read.mma and lawley loaded.
source('utils.R')
require(lme4)
n <- 5
mu.null <- 0
tau2.null <- 1
## beta <- runif(2)#c(mle[c('tau2','mu')]) ## should be passed to make.fun
l.derivs.full <- read.mma('./mma-formulas/full/')
l.derivs.redu <- read.mma('./mma-formulas/redu/')
start <- Sys.time()
## m <- rpois(n,lambda=60)
m <- rep(25,n)
stats <- replicate(n^3, {
    cat('.')
    theta <- rnorm(n)
    e <- rbinom(n,m,plogis(mu.null+sqrt(tau2.null)*theta))
    if(any(e==0)||any(e==m))return(c(NA,NA))
    ## rm(theta)
    dat <- data.frame(e = e, non = m - e, study = factor(seq_along(e)),off = mu.null)
    fit.alt <- glmer(cbind(e, non) ~ 1 + (1 | study), data = dat,family = binomial(link = "logit"),nAGQ = 7)
    fit.null <- glmer(cbind(e, non) ~ 0 + offset(off) + (1 | study),data = dat, family = binomial(link = "logit"), nAGQ = 7)
    ## LR <- 2 * (as.numeric(logLik(fit.alt)) - as.numeric(logLik(fit.null)))
    loglik.full <- as.numeric(logLik(fit.alt))
    loglik.redu <- as.numeric(logLik(fit.null))
    LR <- 2 * (loglik.full - loglik.redu)
    tau2.mle.redu <- VarCorr(fit.null)$study |> as.numeric() |> unname()
    beta <- c(tau2.mle.redu,mu.null)
    eps.full <- lawley(p=2,m=m,beta=beta,l.derivs=l.derivs.full)
    eps.redu <- lawley(p=1,m=m,beta=beta,l.derivs=l.derivs.redu)
    correction <- 1 + eps.full - eps.redu
    if(correction<0)browser()
    ## correction <- 1#!!!
    if(abs(eps.full)>7 || abs(eps.redu)>7)browser()    
    c(LR=LR,LR.bartlett=LR/correction,correction=correction,eps.full=eps.full,eps.redu=eps.redu,loglik.full=loglik.full,loglik.redu=loglik.redu)
})
## correction
print(Sys.time() - start)
plot(ecdf(stats['LR',]))
plot(ecdf(stats['LR.bartlett',]),add=TRUE,col=2)
curve(pchisq(x,df=1),add=TRUE,col=3)
apply(stats,1,mean)
apply(stats,1,summary)




## 2j. process sim results
save.path <- './sim/260706/'
## png.file <- './figs/gaussian.png'
filelist <- dir(save.path)
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]
print(length(filelist))
by.sim <- lapply(filelist, function(file) {
    ## print(file)
    if('m'%in%ls())rm(m)
    if('m.mean'%in%ls())rm(m.mean)
    load(paste0(save.path,file))
    list(stats=stats,n=n,m=    if('m'%in%ls())unique(m)else m.mean)
    ## by.n
})
## sim <- by.sim[[1]]
long <- sapply(by.sim,function(sim) {
    ## cat('.')
    stats <- sim$stats
    if(is.list(stats))    stats <- stats[!sapply(stats,function(v)any(is.na(v)))] |> simplify2array()
    ## LR <- stats['LR',]; LR.bartlett <- stats['LR.bartlett',]
    stats <- stats[c('LR','LR.bartlett'),]
    ecdfs <- apply(stats,1,function(x)ecdf(pchisq(x,df=1,lower=FALSE)))
    ks <- sapply(ecdfs,function(ecdf)max(abs(ecdf(knots(ecdf)) - knots(ecdf))))
    c(ks,n=sim$n,m=sim$m)
}) |> t()
long <- aggregate(cbind(LR,LR.bartlett) ~ ., FUN=median, data=long)
long <- subset(long,m==55)
plot(LR ~ n, type='l',data=long,ylim=range(c(long$LR,long$LR.bartlett)))
lines(LR.bartlett ~ n, data=long,col=2)
coef(lm(log(LR)~log(n),data=long))[2]
coef(lm(log(LR.bartlett)~log(n),data=long))[2]





## python check
source('utils.R')
                                        # require(lme4)
n <- 10
m <- rep(20,n)
p <- 1
mu.null <- 0
tau2.null <- 1/2
beta <- c(tau2.null,mu.null)
## beta <- runif(2)#c(mle[c('tau2','mu')]) ## should be passed to make.fun
l.derivs.full <- read.mma('./mma-formulas/full/')
l.derivs.redu <- read.mma('./mma-formulas/redu/')
k.vals <- sapply(l.derivs.redu, function(l.deriv) {
    ## browser()
    sapply(1:length(m), function(j) 
        sapply(1:(m[j]-1), function(k) {
            ## browser()
            integrate(function(theta)l.deriv(theta,m[j],k,beta)*dnorm(theta),-Inf,Inf)$val
        }) |> sum()
        ) |> sum()
})
idx.strings <- gsub('.txt','',names(k.vals),fixed=TRUE) |> sub(pattern='.*k-',replace='') |> strsplit(split='-')
idx.strings <- lapply(idx.strings,function(x)if(length(x)==1)c(x,0)else x)
k.df <- data.frame(inner=sapply(idx.strings,`[`,1),outer=sapply(idx.strings,`[`,2),k.val=k.vals)
k <- function(inner,outer=0)k.df$k.val[k.df$inner==paste(inner,collapse='') & k.df$outer==paste(outer,collapse='')]
inf <- -matrix(apply(expand.grid(1:p,1:p),1,k),nrow=p)
inv.inf <- solve(inf)


## 3. using gauss-hermite quadrature


## 3a. compute bartlett correction in single-arm binomial model
source('utils.R')
tic <- Sys.time()
p <- 1
mma.str.path <- './mma-formulas/quad/redu/'
## p <- 2
## mma.str.path <- './mma-formulas/full/'
mu.null <- 0
tau2.null <- 1/2
m <- 20
n.GH <- 80
GH <- fastGHQuad::gaussHermiteData(n.GH)
filelist <- dir(mma.str.path)
filelist <- paste0(mma.str.path,filelist[grep('\\.txt',filelist)])
mma.strings <- sapply(filelist,function(file)paste(readLines(file,warn=FALSE),collapse=''))
k.vals <- parallel::mclapply(mma.strings,mc.cores=parallel::detectCores()-2,FUN=function(mma.str) {
    cat('.')
    mma.str <- mma.to.R(mma.str,
                        funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial"),
                        rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose")
                        )
    bindings <- list(u=plogis(mu.null + sqrt(2)*sqrt(tau2.null)*GH$x),
                     m=m,
                     K=n.GH,
                     theta=GH$x,
                     w=GH$w,
                     beta=c(tau2.null,mu.null),
                     Pi=pi)
    pmf <- str.to.fun(str=mma.str,args = c('x'),bindings=bindings)
    sum(sapply(0:m,pmf))
}) |> simplify2array()
idx.strings <- gsub('.txt','',names(k.vals),fixed=TRUE) |> sub(pattern='.*k-',replace='') |> strsplit(split='-')
idx.strings <- lapply(idx.strings,function(x)if(length(x)==1)c(x,0)else x)
k.df <- data.frame(inner=sapply(idx.strings,`[`,1),outer=sapply(idx.strings,`[`,2),k.val=k.vals)
k <- function(inner,outer=0)k.df$k.val[k.df$inner==paste(sort(inner),collapse='') & k.df$outer==paste(sort(outer),collapse='')]
inf <- -matrix(apply(expand.grid(1:p,1:p),1,k),nrow=p)
inv.inf <- solve(inf)
lawley.4 <- apply(do.call(expand.grid,rep(list(1:p),4)), 1, function(idx.str) {
    r <- idx.str[1]; s <- idx.str[2]; t <- idx.str[3]; u <- idx.str[4]
    inv.inf[r,s] * inv.inf[t,u] * (1/4*k(c(r,s,t,u)) - k(c(r,s,t),u) + k(c(r,t),c(s,u)))
})
lawley.6 <- apply(do.call(expand.grid,rep(list(1:p),6)), 1, function(idx.str) {
    r <- idx.str[1]; s <- idx.str[2]; t <- idx.str[3]; u <- idx.str[4]; v <- idx.str[5]; w <- idx.str[6]
    -inv.inf[r,s] * inv.inf[t,u] * inv.inf[v,w] * (1/6*k(c(r,t,v))*k(c(s,u,w)) + 1/4*k(c(r,t,u))*k(c(s,v,w)) - k(c(r,t,v))*k(c(s,w),u) - k(c(r,t,u))*k(c(s,w),v) + k(c(r,t),v)*k(c(s,w),u) + k(c(r,t),u)*k(c(s,w),v))
})
sum(lawley.4) - sum(lawley.6) 
print(Sys.time()-tic)


## encapsulate

source('utils.R')
## p <- 1
## mma.str.path <- './mma-formulas/redu/'
## p <- 2
## mma.str.path <- './mma-formulas/full/'
lawley <- function(k.strings,p,m,mu.null,tau2.hat,n.GH=80) {
    GH <- fastGHQuad::gaussHermiteData(n.GH)
    k.vals <- parallel::mclapply(k.strings,mc.cores=parallel::detectCores()-2,FUN=function(k.str) {
        cat('.')
        bindings <- list(u=plogis(mu.null + sqrt(2)*sqrt(tau2.hat)*GH$x),
                         m=m,
                         K=n.GH,
                         theta=GH$x,
                         w=GH$w,
                         beta=c(tau2.hat,mu.null),
                         Pi=pi)
        pmf <- str.to.fun(str=k.str,args = c('x'),bindings=bindings)
        sum(sapply(0:m,pmf))
    }) |> simplify2array()
    idx.strings <- gsub('.txt','',names(k.vals),fixed=TRUE) |> sub(pattern='.*k-',replace='') |> strsplit(split='-')
    idx.strings <- lapply(idx.strings,function(x)if(length(x)==1)c(x,0)else x)
    k.df <- data.frame(inner=sapply(idx.strings,`[`,1),outer=sapply(idx.strings,`[`,2),k.val=k.vals)
    ## print(k.df)
    k <- function(inner,outer=0)k.df$k.val[k.df$inner==paste(sort(inner),collapse='') & k.df$outer==paste(sort(outer),collapse='')]
    inf <- -matrix(apply(expand.grid(1:p,1:p),1,k),nrow=p)
    inv.inf <- solve(inf)
    lawley.4 <- apply(do.call(expand.grid,rep(list(1:p),4)), 1, function(idx.str) {
        r <- idx.str[1]; s <- idx.str[2]; t <- idx.str[3]; u <- idx.str[4]
        inv.inf[r,s] * inv.inf[t,u] * (1/4*k(c(r,s,t,u)) - k(c(r,s,t),u) + k(c(r,t),c(s,u)))
    })
    lawley.6 <- apply(do.call(expand.grid,rep(list(1:p),6)), 1, function(idx.str) {
        r <- idx.str[1]; s <- idx.str[2]; t <- idx.str[3]; u <- idx.str[4]; v <- idx.str[5]; w <- idx.str[6]
        -inv.inf[r,s] * inv.inf[t,u] * inv.inf[v,w] * (1/6*k(c(r,t,v))*k(c(s,u,w)) + 1/4*k(c(r,t,u))*k(c(s,v,w)) - k(c(r,t,v))*k(c(s,w),u) - k(c(r,t,u))*k(c(s,w),v) + k(c(r,t),v)*k(c(s,w),u) + k(c(r,t),u)*k(c(s,w),v))
    })
    bias <- out <- sum(lawley.4) - sum(lawley.6)
    return(out)
}
tic <- Sys.time()
mu.null <- 0
tau2.null <- 1/2
m <- 20
k.strings <- lapply(c(full='./mma-formulas/quad/full/', redu='./mma-formulas/quad/redu/'), function(mma.str.path) {
    filelist <- dir(mma.str.path)
    filelist <- paste0(mma.str.path,filelist[grep('\\.txt',filelist)])
    mma.strings <- sapply(filelist,function(file)paste(readLines(file,warn=FALSE),collapse=''))
    k.strings <- sapply(mma.strings, function(mma.str) mma.to.R(mma.str,
                                                                funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial"),
                                                                rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose")
                                                                ))
})
eps.redu <- lawley(k.strings=k.strings[['redu']],p=1,m=m,mu.null=0,tau2.hat=1/2)
eps.full <- lawley(k.strings=k.strings[['full']],p=2,m=m,mu.null=0,tau2.hat=1/2)
print(Sys.time()-tic)
## > eps.full
## [1] 1.734166
## > eps.redu
## [1] 0.4079947


## 3b. rate sim

## set.seed(1)
tic <- Sys.time()
bartlett.C <- 1.734166 - 0.4079947 #tied to mi=20 for all i
mu.true <- mu.null <- 0
tau2.true <- tau2.null <- 1/2
ns <- 10:20
by.n <- parallel::mclapply(ns,mc.cores=4,FUN=function(n) {
    ## by.n <- sapply(ns, function(n) {
    cat('.')
    m <- rep(20,n)
    stats <- replicate(n^3, {
        theta <- rnorm(n)
        pr <- plogis(mu.true+sqrt(tau2.true)*theta)
        x <- rbinom(n,size=m,prob=pr)
        GH <- fastGHQuad::gaussHermiteData(30)
        loglik.unit <- function(mu,tau2,mi,xi,GH=NULL) {
            if(is.null(GH)) GH <- fastGHQuad::gaussHermiteData(20)
            pr <- function(theta)plogis(mu+sqrt(tau2)*theta)
            integrand <- function(theta)pr(theta)^xi*(1-pr(theta))^(mi-xi)
            integral.approx <- choose(mi,xi)*sum(GH$w*integrand(sqrt(2)*GH$x))/sqrt(pi)
            return(log(integral.approx))
        }
        ## loglik.unit(mu,tau2,m[2],x[3])
        loglik <- function(mu,tau2,m,x) sum(sapply(seq_along(x),function(i)loglik.unit(mu,tau2,m[i],x[i],GH=GH)))
        ## loglik(mu.null,tau2.null,m,x)
        logtau2.mle.redu <- optimize(function(logtau2)loglik(mu.null,exp(logtau2),m,x),interval=c(-5,5),maximum=TRUE)$max
        tau2.mle.redu <- exp(logtau2.mle.redu)
        mle.full <- optim(par=c(mu.null,log(tau2.mle.redu)),fn=function(par)-loglik(par[1],exp(par[2]),m,x),method='BFGS')$par
        mu.mle.full <- mle.full[1]; tau2.mle.full <- exp(mle.full[2])
        ## print(mle.full)
        LR <- 2 * (loglik(mu.mle.full,tau2.mle.full,m,x) - loglik(mu.null,tau2.mle.redu,m,x))
        LR.bartlett <- LR / (1+bartlett.C/n)
        c(LR=LR,LR.bartlett=LR.bartlett)
    })
    ecdfs <- apply(stats,1,function(x)ecdf(pchisq(x,df=1,lower=FALSE)))
    ks <- sapply(ecdfs,function(ecdf)max(abs(ecdf(knots(ecdf)) - knots(ecdf))))
}) |> simplify2array()
print(Sys.time()-tic)
plot(ns,by.n['LR',],ylim=range(by.n),type='l')
lines(ns,by.n['LR.bartlett',],col=2)
apply(by.n,1,function(ks)coef(lm(log(ks)~log(ns)))[2])



## process  rate sim (IID data, using R with gauss-hermite)
save.path <- './sim/260714/'
save.path <- './sim/260722/'
## png.file <- './figs/gaussian.png'
filelist <- dir(save.path)
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]
print(length(filelist))
by.sim <- sapply(filelist, function(file) {
    load(paste0(save.path,file))
    apply(by.n,1,function(ks)coef(lm(log(ks)~log(ns)))[2])
})
apply(by.sim,1,median)





## 3c. debugging--ehcking why laplace is producing such different k's.

## compare quadrature and laplace k's. similar for large m (eg 80)
## when using slow 2nd order approx. large m and 1st order still way
## off.  get quad k's
source('utils.R')
m <- 20
mu.null <- 0
tau2.null <- 1/2
n.GH <- 20
GH <- fastGHQuad::gaussHermiteData(n.GH)
mma.str.path <- './mma-formulas/quad-old/redu/'
filelist <- dir(mma.str.path)
filelist <- paste0(mma.str.path,filelist[grep('\\.txt',filelist)])
mma.strings <- sapply(filelist,function(file)paste(readLines(file,warn=FALSE),collapse=''))
k.vals <- parallel::mclapply(mma.strings,mc.cores=parallel::detectCores()-2,FUN=function(mma.str) {
    cat('.')
    mma.str <- mma.to.R(mma.str,
                        funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial"),
                        rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose")
                        )
    bindings <- list(u=plogis(mu.null + sqrt(2)*sqrt(tau2.null)*GH$x),
                     m=m,
                     K=n.GH,
                     theta=GH$x,
                     w=GH$w,
                     beta=c(tau2.null,mu.null),
                     Pi=pi)
    pmf <- str.to.fun(str=mma.str,args = c('x'),bindings=bindings)
    sum(sapply(0:m,pmf))
}) |> simplify2array()
idx.strings <- gsub('.txt','',names(k.vals),fixed=TRUE) |> sub(pattern='.*k-',replace='') |> strsplit(split='-')
idx.strings <- lapply(idx.strings,function(x)if(length(x)==1)c(x,0)else x)
k.df <- data.frame(inner=sapply(idx.strings,`[`,1),outer=sapply(idx.strings,`[`,2),k.val=k.vals)
colnames(k.df) <- NULL
k.df.quad <- k.df
## read in laplace k's
beta <- c(tau2.null,mu.null)
l.derivs <- read.mma('./mma-formulas/laplace/redu/')
k.vals <- sapply(l.derivs, function(l.deriv) {
    sapply(1:(m-1), function(x) {
        integrate(function(theta)l.deriv(theta,m,x,beta)*dnorm(theta),-Inf,Inf)$val
    }) |> sum()
})
idx.strings <- gsub('.txt','',names(k.vals),fixed=TRUE) |> sub(pattern='.*k-',replace='') |> strsplit(split='-')
idx.strings <- lapply(idx.strings,function(x)if(length(x)==1)c(x,0)else x)
k.df <- data.frame(inner=sapply(idx.strings,`[`,1),outer=sapply(idx.strings,`[`,2),k.val=k.vals)
colnames(k.df) <- NULL
k.df.laplace <- k.df
print(k.df.quad)
print(k.df.laplace)


## look at what step the big error is creeping in. seems when the expectation is taken of the inner derivative.

source('utils.R')
file <- './mma-formulas/tryLaplace.txt'
mma.str <-paste(readLines(file,warn=FALSE),collapse='')
mma.str <-  mma.to.R(mma.str,
                     funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial"),
                     rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose"))
l.deriv <- 
    str.to.fun(
        gsub(pattern=c('m\\[j\\]'),replacement=c('m'),mma.str),
        args = c('theta','m','x','beta'),
        binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt)
    )
m <- 20
## x <- 7
beta <- c(1/2,0)
pr <- function(theta)plogis(beta[2]+sqrt(beta[1])*theta)
## ## 1. check loglikelihood
## approx <- function(x){
##     r <- function(theta)-(x/m*log(pr(theta))+(1-x/m)*log(1-pr(theta)))
##     r.deriv2 <- function(theta)-beta[1]*(pr(theta)-1)*pr(theta)
##     theta.min <- 1/sqrt(beta[1])*(beta[2]+log(x/(m-x)))
##     log(exp(-m*r(theta.min))*dnorm(theta.min)*sqrt(2*pi)/sqrt(m)/sqrt(r.deriv2(theta.min)))
## }
## approx <- function(x)l.deriv(NA,m,x,beta)
## true <- function(x)log(integrate(function(theta)pr(theta)^x*(1-pr(theta))^(m-x)*dnorm(theta),-Inf,Inf)$val)
## vals <- sapply(1:(m-1),function(x)c(approx=approx(x),true=true(x)))
## plot(1:(m-1),vals['approx',],ylim=range(vals))
## points(1:(m-1),vals['true',],col=2)
## ## 2. check expected loglikelihood ie kUnit[{}]
## sapply(1:(m-1),function(x)integrate(function(theta)l.deriv(theta,m,x,beta)*dnorm(theta),-Inf,Inf)$val) |> sum()
## sapply(1:(m-1),function(x)log(integrate(function(theta)pr(theta)^x*(1-pr(theta))^(m-x)*dnorm(theta),-Inf,Inf)$val)*integrate(function(theta2)dbinom(x,m,pr(theta2))*dnorm(theta2),-Inf,Inf)$val) |> sum()
## 3. check derivative wrt beta[1] of loglik
approx <- function(x)integrate(function(theta)l.deriv(theta,m,x,beta)*dnorm(theta),-Inf,Inf)$val
I <- function(x)integrate(function(theta)pr(theta)^x*(1-pr(theta))^(m-x)*dnorm(theta),-Inf,Inf)$val
true <- function(x)1/I(x) * integrate(function(theta)dnorm(theta)*(x*pr(theta)^x*(1-pr(theta))^(m-x+1)-pr(theta)^(x+1)*(m-x)*(1-pr(theta))^(m-x))*theta/2/sqrt(beta[1]),-Inf,Inf)$val
vals <- sapply(1:(m-1),function(x)c(approx=approx(x),true=true(x)))
plot(1:(m-1),vals['approx',],ylim=range(vals))
points(1:(m-1),vals['true',],col=2)
abline(h=0)
## 4. expected derivative derivative wrt beta[1] of loglik ie
## k[{1}]. this is where the approximation breaks down.
pmf <- function(x)choose(m,x)*I(x)
pmf.vals <- sapply(1:(m-1),function(x)pmf(x))
points(1:(m-1),-log(pmf.vals),col=3)
sum(sapply(1:(m-1),function(x)approx(x)*pmf(x)))
sum(sapply(1:(m-1),function(x)true(x)*pmf(x)))


## 3d. try to speed up G-H

source('utils.R')
## m <- 20
ms <- rep(20,10)
mu.null <- 0
tau2.null <- 1/2
beta <- c(tau2.null,mu.null)
n.GH <- 20
GH <- fastGHQuad::gaussHermiteData(n.GH)
## read in g derivs, these depend on m, theta, and x
mma.str.path <- './mma-formulas/quad/g/'
filelist <- dir(mma.str.path)
filelist <- paste0(mma.str.path,filelist[grep('\\.txt',filelist)])
## file <- './mma-formulas/quad/g/g-10.txt'
g.derivs <- lapply(filelist, function(file) {
    mma.str <-paste(readLines(file,warn=FALSE),collapse='')
    mma.str <- gsub('Pi','pi',mma.str)
    mma.str <-  mma.to.R(mma.str,
                         funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial"),
                         rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose"))
    g.deriv <- 
        str.to.fun(mma.str,
                   args = c('theta','x','m'),
                   binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt,beta=beta)
                   )
})
orders.g <- paste0('gSum',substr(filelist,nchar(filelist)-5,nchar(filelist)-4))
mma.str.path <- './mma-formulas/quad/k/'
filelist <- dir(mma.str.path)
filelist <- paste0(mma.str.path,filelist[grep('\\.txt',filelist)])
## file <- './mma-formulas/quad/g/g-10.txt'
k.funs <- lapply(filelist, function(file) {
    mma.str <-paste(readLines(file,warn=FALSE),collapse='')
    ## mma.str <-  mma.to.R(mma.str,
    ##                      funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial"),
    ##                      rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose"))
    matches.idx <- gregexpr("gSum[[:digit:]][[:digit:]]", mma.str)
    variables <- regmatches(mma.str, matches.idx)[[1]] |> unique()
    ## all.vars <<- c(all.vars,variables)
    fun.str <- sprintf(
        "function(%s) { %s }", 
        paste(variables, collapse = ", "), 
        mma.str
    )
    k.fun <- l.deriv.fun <- eval(parse(text = fun.str))
    list(fun=k.fun,args=variables)
})
orders.k <- gsub('.txt','',filelist,fixed=TRUE) |> sub(pattern='.*k-',replace='') |> strsplit(split='-')
orders.k <- lapply(orders.k,function(x)if(length(x)==1)c(x,0)else x)
k.vals.by.m <- sapply(ms, function(m) {
    cat('.')
    ## for given m, sum g.derivs over theta, leaving functions of x
    ## g.deriv <- g.derivs[[1]]
    g.sums.x <- lapply(g.derivs, function(g.deriv) {
        g.deriv.mat <- outer(X=GH$x,Y=0:m,function(theta,x)g.deriv(theta,x,m))
        g.deriv.x <- colSums(GH$w * g.deriv.mat)
    })
    ## orders <- substr(filelist,nchar(filelist)-5,nchar(filelist)-4)
    ## colnames(g.sums.x) <- variable.names
    names(g.sums.x) <- orders.g
    ## for given m, compute k values by weighted sum over x
    k.vals <- sapply(seq_along(k.funs), function(i) {
        k.fun <- k.funs[[i]]$fun; args <- k.funs[[i]]$args
        k.fun.x <- do.call(k.fun, g.sums.x[args])
        k.val <- sum( choose(m,0:m) * k.fun.x )
    })
    return(k.vals)
})
k.vals <- rowSums(k.vals.by.m)
bartlett.C <- lawley(k.vals,2) - lawley(k.vals,1)


## 3d-1 convert mma strings to R strings to avoid eval(parse) slowdown

## g orders needed for input to mma:
all.vars <- unique(all.vars)
paste0(paste0('{',substr(all.vars,5,5),',',substr(all.vars,6,6),'}'),collapse=', ')
[1] "{1,0}, {0,0}, {2,0}, {3,0}, {4,0}, {0,1}, {1,1}, {2,1}, {3,1}, {0,2}, {1,2}, {2,2}, {0,3}, {1,3}, {0,4}"

## convert mma strings of k functions to R functions, save the R
## functions to disk in a single file.
mma.str.path <- './mma-formulas/quad/k/'
filelist <- dir(mma.str.path)
filelist <- paste0(mma.str.path,filelist[grep('\\.txt',filelist)])
## file <- './mma-formulas/quad/g/g-10.txt'
ks.fun <- lapply(filelist, function(file) {
    mma.str <-paste(readLines(file,warn=FALSE),collapse='')
    matches.idx <- gregexpr("gSum[[:digit:]][[:digit:]]", mma.str)
    variables <- regmatches(mma.str, matches.idx)[[1]] |> unique()
    fun.str <- sprintf(
        "function(%s) { %s }", 
        paste(variables, collapse = ", "), 
        mma.str
    )
    k.fun <- l.deriv.fun <- eval(parse(text = fun.str))
    list(fun=k.fun,args=variables)
})
orders.k <- gsub('.txt','',filelist,fixed=TRUE) |> sub(pattern='.*k-',replace='') |> strsplit(split='-')
orders.k <- lapply(orders.k,function(x)if(length(x)==1)c(x,0)else x)

ks.fun <- lapply(seq_along(ks.fun), function(i)c(ks.fun[[i]],order=list(orders.k[[i]])))
saveRDS(ks.fun,file='ks.fun.rds',compress=FALSE)



## 3d-2 convert mma strings of g functions to R functions of theta,beta,x, and m, save the R functions to disk.
mma.str.path <- './mma-formulas/quad/g/'
filelist <- dir(mma.str.path)
filelist <- paste0(mma.str.path,filelist[grep('\\.txt',filelist)])
## file <- './mma-formulas/quad/g/g-10.txt'
gs.fun <- lapply(filelist, function(file) {
    mma.str <-paste(readLines(file,warn=FALSE),collapse='')
    mma.str <- gsub('Pi','pi',mma.str)
    mma.str <-  mma.to.R(mma.str,
                         funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial"),
                         rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose"))
    g.deriv <- 
        str.to.fun(mma.str,
                   args = c('beta','theta','x','m'),
                   binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt))
})
g.orders <- paste0('gSum',substr(filelist,nchar(filelist)-5,nchar(filelist)-4))
names(gs.fun) <- g.orders
saveRDS(gs.fun,file='gs.fun.rds',compress=FALSE)




## read in pmf derivatives -- using pmf recurrences in mma
source('utils.R')
mma.str.path <- './mma-formulas/pmf-binomial-prop/'
filelist <- dir(mma.str.path)
filelist <- paste0(mma.str.path,filelist[grep('pmf-[[:digit:]][[:digit:]]\\.txt',filelist)])
pmf.derivs.fun <- lapply(filelist, function(file) {
    mma.str <-paste(readLines(file,warn=FALSE),collapse='')
    mma.str <-  mma.to.R(mma.str,
                         funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial","factorial","factMoment","pmf"),
                         rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose",factorial="factorial"))
    pmf <- function(x,params,prob) {
        m <- params['m']
        dbinom(x=x,size=m,prob=prob)
    }
    ## factMoment <- function(k,params,odds) {
    ##     ## print(params)
    ##     m.T <- params['m.T']
    ##     m.C <- params['m.C']
    ##     s <- params['s']
    ##     supp.x <- max(0,s-m.C):min(s,m.T)
    ##     sum ( apply(outer(supp.x, 0:(k-1), `-`), 1, prod) * pmf(supp.x,params,odds) )
    ## }
    pmf.deriv <- 
        str.to.fun(mma.str,
                   args = c('beta','params','theta','x'),
                   ## args = c('idx','x','theta','beta','params'),
                   binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt,pmf=pmf))
    ## })
})
pmf.deriv.orders <- paste0('pmfSum',substr(filelist,nchar(filelist)-5,nchar(filelist)-4))
names(pmf.derivs.fun) <- pmf.deriv.orders
## saveRDS(pmf.derivs.fun,file='pmf-derivs-fun-binomial-prop.rds',compress=FALSE)



## 3d-3. further speedup: reorganize sped up version to allow for parallelization, use saved
## R functions rather than mma strings

source('utils.R')
## m <- 20
mu.null <- 0
tau2.null <- 1/2
beta <- c(tau2.null,mu.null)

## 1. read in k's as functions of g derivatives. these are rational
## expressions of the g derivatives and are common to any parameters
## and data.

ks.fun <- readRDS('ks.fun.rds')
k.orders <- lapply(ks.fun,function(k.fun)k.fun$order)

## 2. read in g's depending on m, theta, x, and beta. obtain ma-level
## g(m,theta,x) by currying mle's.

gs.fun <- readRDS('gs.fun.rds')
gs.ma <- lapply(gs.fun,function(g.fun)function(theta,x,m)g.fun(beta,theta,x,m))

## 3. fixing m, for each g get g values at all support x points by
## weighted summing over theta for each support point x, leaving
## g(x1),...,g(xm)=g(0),...,g(m). apply k functions. this is study
## level and is the parallelizable part.

ks.study <- function(m,gs.ma,ks.fun,GH=NULL) {
    if(is.null(GH)) {
        GH <- fastGHQuad::gaussHermiteData(30)
        names(GH) <- c('theta.grid','w')
    }
    ## sum ma-level g's over theta, getting a vector g(x,m) with
    ## unit-level components
    gs.unit <- lapply(gs.ma, function(g) {
        g.mat <- outer(X=GH$theta.grid,Y=0:m,function(theta,x)g(theta,x,m))
        g.unit <- colSums(GH$w * g.mat)
    })
    ## apply ks.fun to unit-level g's to get unit-level k's, take sum
    ## weighted by binomial coefficient to get study-level k
    k.study <- sapply(seq_along(ks.fun), function(i) {
        k.fun <- ks.fun[[i]]$fun; args <- ks.fun[[i]]$args
        k.fun.unit <- do.call(k.fun, gs.unit[args])
        k.val <- sum( choose(m,0:m) * k.fun.unit )
    })
    return(k.study)
}

ms <- rep(20,10)
n.GH <- 20
GH <- fastGHQuad::gaussHermiteData(n.GH)
names(GH) <- c('theta.grid','w')
k.vals <- ks.study(ms[1],gs.ma,ks.fun)
k.df <- data.frame(inner=sapply(k.orders,`[`,1),outer=sapply(k.orders,`[`,2),k.val=k.vals)
bartlett.C <- lawley(k.df,2) - lawley(k.df,1)





## 3d-4. compare speed of my computation of the mle's with glmer's. glmer much slower.

set.seed(1)
mu.true <- mu.null <- 0
tau2.true <- tau2.null <- 1/2
n <- 10
m <- 20
theta <- rnorm(n)
pr <- plogis(mu.true+sqrt(tau2.true)*theta)
x <- rbinom(n,size=m,prob=pr)
times <- replicate(1e2, {
    tic <- Sys.time()
    GH <- fastGHQuad::gaussHermiteData(30)
    loglik.unit <- function(mu,tau2,mi,xi,GH=NULL) {
        if(is.null(GH)) GH <- fastGHQuad::gaussHermiteData(20)
        pr <- function(theta)plogis(mu+sqrt(tau2)*theta)
        integrand <- function(theta)pr(theta)^xi*(1-pr(theta))^(mi-xi)
        integral.approx <- choose(mi,xi)*sum(GH$w*integrand(sqrt(2)*GH$x))/sqrt(pi)
        return(log(integral.approx))
    }
    ## loglik.unit(mu,tau2,m[2],x[3])
    loglik <- function(mu,tau2,m,x) sum(sapply(seq_along(x),function(i)loglik.unit(mu,tau2,m,x[i],GH=GH)))
    ## loglik(mu.null,tau2.null,m,x)
    logtau2.mle.redu <- optimize(function(logtau2)loglik(mu.null,exp(logtau2),m,x),interval=c(-5,5),maximum=TRUE)$max
    tau2.mle.redu <- exp(logtau2.mle.redu)
    mle.full <- optim(par=c(mu.null,log(tau2.mle.redu)),fn=function(par)-loglik(par[1],exp(par[2]),m,x),method='BFGS')$par
    mu.mle.full <- mle.full[1]; tau2.mle.full <- exp(mle.full[2])
    LR <- 2 * (loglik(mu.mle.full,tau2.mle.full,m,x) - loglik(mu.null,tau2.mle.redu,m,x))
    toc <- Sys.time()
    toc - tic
})
hist(times)
summary(times)

set.seed(1)
mu.true <- mu.null <- 0
tau2.true <- tau2.null <- 1/2
n <- 10
m <- 20
theta <- rnorm(n)
pr <- plogis(mu.true+sqrt(tau2.true)*theta)
x <- rbinom(n,size=m,prob=pr)
times <- replicate(1e2, {
    tic <- Sys.time()
    dat <- data.frame(x=x, non = m - x, study = factor(seq_along(x)),off = mu.null)
    fit.alt <- glmer(cbind(x, non) ~ 1 + (1 | study), data = dat,family = binomial(link = "logit"),nAGQ = 7)
    fit.null <- glmer(cbind(x, non) ~ 0 + offset(off) + (1 | study),data = dat, family = binomial(link = "logit"), nAGQ = 7)
    LR <- 2 * (as.numeric(logLik(fit.alt)) - as.numeric(logLik(fit.null)))
    toc <- Sys.time()
    toc - tic
})
hist(times,add=TRUE,col=2)
summary(times)




## re-parameterize the sinh-ainh family

## 3e. make code generic to allow for other models



## 3e-1 computing mle's in model with normal RE. swap out loglik.unit.
set.seed(1)
mu.true <- mu.null <- 0
tau2.true <- tau2.null <- 1/2
n <- 10
m <- rep(20,n)
params=cbind(m=m)
theta <- rnorm(n)
pr <- plogis(mu.true+sqrt(tau2.true)*theta)
x <- rbinom(n,size=m,prob=pr)
## times <- replicate(1e2, {
## tic <- Sys.time()
GH <- fastGHQuad::gaussHermiteData(30)
GH$x <- sqrt(2)*GH$x; GH$w <- GH$w/sqrt(pi)
loglik.unit <- function(mu,tau2,xi,paramsi,GH=NULL) {
    mi <- paramsi['m']
    pr <- function(theta)plogis(mu+sqrt(tau2)*theta)
    integrand <- function(theta)pr(theta)^xi*(1-pr(theta))^(mi-xi)
    integral.approx <- if(is.null(GH)) integrate(function(theta)integrand(theta)*dnorm(theta),-Inf,Inf)$val
    else
        choose(mi,xi)*sum(GH$w*integrand(GH$x))
    return(log(integral.approx))
}
loglik <- function(mu,tau2,x,params,GH=NULL) sum(sapply(seq_along(x),function(i)loglik.unit(mu,tau2,x[i],params[i,],GH=GH)))
logtau2.mle.redu <- optimize(function(logtau2)loglik(mu.null,exp(logtau2),x,params),interval=c(-5,5),maximum=TRUE)$max
tau2.mle.redu <- exp(logtau2.mle.redu)
mle.full <- optim(par=c(mu.null,log(tau2.mle.redu)),fn=function(par)-loglik(par[1],exp(par[2]),x,params),method='BFGS')$par
mu.mle.full <- mle.full[1]; tau2.mle.full <- exp(mle.full[2])
LR <- 2 * (loglik(mu.mle.full,tau2.mle.full,x,params) - loglik(mu.null,tau2.mle.redu,x,params))
LR




## 3e-2 read in R-format formulas and compute correction

source('utils.R')
## m <- 20
mu.null <- 0
tau2.null <- 1/2
beta <- c(tau2.null,mu.null)

## 1. read in k's as functions of f derivatives. these are rational
## expressions of the g derivatives and are common to any model,
## parameters and data.


## ## k's: this is for all models
## mma.str.path <- './mma-formulas/k/'
## filelist <- dir(mma.str.path)
## filelist <- paste0(mma.str.path,filelist[grep('\\.txt',filelist)])
## ## file <- './mma-formulas/quad/g/g-10.txt'
## ks.fun <- lapply(filelist, function(file) {
##     mma.str <-paste(readLines(file,warn=FALSE),collapse='')
##     matches.idx <- gregexpr("pmfSum[[:digit:]][[:digit:]]", mma.str)
##     variables <- regmatches(mma.str, matches.idx)[[1]] |> unique()
##     fun.str <- sprintf(
##         "function(%s) { %s }", 
##         paste(variables, collapse = ", "), 
##         mma.str
##     )
##     k.fun <- l.deriv.fun <- eval(parse(text = fun.str))
##     list(fun=k.fun,args=variables)
## })
## k.orders <- gsub('.txt','',filelist,fixed=TRUE) |> sub(pattern='.*k-',replace='') |> strsplit(split='-')
## k.orders <- lapply(k.orders,function(x)if(length(x)==1)c(x,0)else x)
## ks.fun <- lapply(seq_along(ks.fun), function(i)c(ks.fun[[i]],order=list(k.orders[[i]])))

## ## saveRDS(ks.fun,file='ks.fun.rds',compress=FALSE)
ks.fun <- readRDS('ks.fun.rds')
k.orders <- lapply(ks.fun,function(k.fun)k.fun$order)

## 2. read in pmf derivatives depending on m, theta, x, and
## beta. obtain ma-level f(m,theta,x) derivatives by currying mle's.

## pmfs: this is specific to the model (eg single-arm logit-normal
mma.str.path <- './mma-formulas/pmf/'
filelist <- dir(mma.str.path)
filelist <- paste0(mma.str.path,filelist[grep('\\.txt',filelist)])
## file <- './mma-formulas/quad/g/g-10.txt'
pmf.derivs.fun <- lapply(filelist, function(file) {
    mma.str <-paste(readLines(file,warn=FALSE),collapse='')
    ## mma.str <- gsub('Pi','pi',mma.str)
    mma.str <-  mma.to.R(mma.str,
                         funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial"),
                         rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose"))
    g.deriv <- 
        str.to.fun(mma.str,
                   args = c('beta','theta','x','m'),
                   binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt))
})
pmf.deriv.orders <- paste0('pmfSum',substr(filelist,nchar(filelist)-5,nchar(filelist)-4))
names(pmf.derivs.fun) <- pmf.deriv.orders
pmf.derivs.fun.ma <- lapply(pmf.derivs.fun,function(pmf.deriv.fun)function(theta,x,m)pmf.deriv.fun(beta,theta,x,m))

## ## saveRDS(pmf.derivs.fun,file='pmf.derivs.fun.rds',compress=FALSE)
## pmf.derivs.fun <- readRDS('pmf.derivs.fun.rds')
## pmf.derivs.fun.ma <- lapply(pmf.derivs.fun,function(pmf.deriv.fun)function(theta,x,m)pmf.deriv.fun(beta,theta,x,m))


## 3. fixing m, for each pmf derivative get pmf deriv values at all
## support x points by weighted summing over theta for each support
## point x, leaving g(x1),...,g(xm)=g(0),...,g(m). apply k
## functions. this is study level and is the parallelizable part.

ks.study <- function(m,pmf.derivs.fun.ma,ks.fun,GH=NULL) {
    if(is.null(GH)) {
        GH <- fastGHQuad::gaussHermiteData(30)
        GH$x <- sqrt(2)*GH$x
        GH$w <- GH$w / sqrt(pi)
        names(GH) <- c('theta.grid','w')
    }
    ## sum ma-level g's over theta, getting a vector g(x,m) with
    ## unit-level components
    pmf.derivs.unit <- lapply(pmf.derivs.fun.ma, function(pmf.deriv.unit) {
        g.mat <- outer(X=GH$theta.grid,Y=0:m,function(theta,x)pmf.deriv.unit(theta,x,m))
        g.unit <- colSums(GH$w * g.mat)
    })
    ## apply ks.fun to unit-level g's to get unit-level k's, take sum
    ## weighted by binomial coefficient to get study-level k
    k.study <- sapply(seq_along(ks.fun), function(i) {
        k.fun <- ks.fun[[i]]$fun; args <- ks.fun[[i]]$args
        k.fun.unit <- do.call(k.fun, pmf.derivs.unit[args])
        k.val <- sum( choose(m,0:m) * k.fun.unit )
    })
    return(k.study)
}

ms <- rep(20,10)
k.vals <- ks.study(ms[1],pmf.derivs.fun.ma,ks.fun)
k.df <- data.frame(inner=sapply(k.orders,`[`,1),outer=sapply(k.orders,`[`,2),k.val=k.vals)
bartlett.C <- lawley(k.df,2) - lawley(k.df,1)
print(bartlett.C)




## 4. two-arm studies via hypergeometric-normal model


## 4a. code for computing the mle's
## require(BiasedUrn)
set.seed(1)
mu.true <- mu.null <- 0
tau2.true <- tau2.null <- 1/2
n <- 10
m.T <- rep(20,n)
m.C <- rep(15,n)
s <- sapply(1:n,function(i)sample(max(m.T[i],m.C[i]),1))
theta <- rnorm(n)
odds <- exp(mu.true+sqrt(tau2.true)*theta)
x <- x.T <- sapply(1:n, function(i) BiasedUrn::rFNCHypergeo(1,m1=m.T[i],m2=m.C[i],n=s[i],odds=odds[i]))
params <- cbind(m.T=m.T,m.C=m.C,s=s)
GH <- fastGHQuad::gaussHermiteData(30)
GH$x <- sqrt(2)*GH$x; GH$w <- GH$w/sqrt(pi)
loglik.unit <- function(mu,tau2,xi,paramsi,GH=NULL) {
    m.Ti <- paramsi['m.T']
    m.Ci <- paramsi['m.C']
    si <- paramsi['s']
    ## pmf <- Vectorize(function(theta)BiasedUrn::dFNCHypergeo(xi,m.Ti,m.Ci,si,exp(mu+sqrt(tau2)*theta))) # slow
    ## pmf <- function(theta) { # numerically unsafe
    ##     numer <- choose(m.Ti,xi)*choose(m.Ci,si-xi)*exp((mu+sqrt(tau2)*theta)*xi) 
    ##     denom <- rowSums(  outer(theta,0:m.Ti,function(thetaj,xj)choose(m.Ti,xj)*choose(m.Ci,si-xj)*exp((mu+sqrt(tau2)*thetaj)*xj) ) )
    ##     numer / denom
    ## }
    pmf <- function(theta) {
        ## logpmf.raw <- function(xi,theta)lgamma(m.Ti+1) - lgamma(xi+1) -lgamma(m.Ti-xi+1) + lgamma(m.Ci+1) - lgamma(si-xi+1) - lgamma(m.Ci-si+xi+1) + (mu+sqrt(tau2)*theta)*xi
        ## lognumer <- logpmf.raw(xi,theta)
        ## logpmf.mat <- outer(0:m.Ti,theta,logpmf.raw) # speed this up
        logpmf.part <- function(xi)lgamma(m.Ti+1) - lgamma(xi+1) -lgamma(m.Ti-xi+1) + lgamma(m.Ci+1) - lgamma(si-xi+1) - lgamma(m.Ci-si+xi+1)
        lognumer <- logpmf.part(xi)+ (mu+sqrt(tau2)*theta)*xi
        logpmf.mat <- logpmf.part(0:m.Ti) + outer(0:m.Ti,theta)*sqrt(tau2) + (0:m.Ti)*mu
        colmaxes <- apply(logpmf.mat,2,max)
        logdenom <- colmaxes + log(rowSums(exp(t(logpmf.mat) - colmaxes)))
        exp(lognumer - logdenom)
    }
    integral.approx <- if(is.null(GH)) integrate(function(theta)pmf(theta)*dnorm(theta),-Inf,Inf)$val
    else
        sum(GH$w*pmf(GH$x))
    return(log(integral.approx))
}
loglik <- function(mu,tau2,x,params,GH=NULL) sum(sapply(seq_along(x),function(i)loglik.unit(mu,tau2,x[i],params[i,],GH=GH)))
logtau2.mle.redu <- optimize(function(logtau2)loglik(mu.null,exp(logtau2),x,params,GH),interval=c(-5,5),maximum=TRUE)$max
tau2.mle.redu <- exp(logtau2.mle.redu)
mle.full <- optim(par=c(mu.null,log(tau2.mle.redu)),fn=function(par)-loglik(par[1],exp(par[2]),x,params,GH),method='BFGS')$par
mu.mle.full <- mle.full[1]; tau2.mle.full <- exp(mle.full[2])
LR <- 2 * (loglik(mu.mle.full,tau2.mle.full,x,params,GH) - loglik(mu.null,tau2.mle.redu,x,params,GH))
## LR
print(tau2.mle.redu)
print(c(mu.mle.full,tau2.mle.full))


## 4b. encapsulate

## LR.gaussian <- function(x,params,mu.null,loglik.unit,n.GH=30) {
##         GH <- fastGHQuad::gaussHermiteData(30)
##         GH$x <- sqrt(2)*GH$x; GH$w <- GH$w/sqrt(pi)
##         loglik <- function(mu,tau2,x,params,GH=NULL) sum(sapply(seq_along(x),function(i)loglik.unit(mu,tau2,x[i],params[i,],GH=GH)))
##         logtau2.mle.redu <- optimize(function(logtau2)loglik(mu.null,exp(logtau2),x,params,GH),interval=c(-5,5),maximum=TRUE)$max
##         tau2.mle.redu <- exp(logtau2.mle.redu)
##         mle.full <- optim(par=c(mu.null,log(tau2.mle.redu)),fn=function(par)-loglik(par[1],exp(par[2]),x,params,GH),method='BFGS')$par
##         mu.mle.full <- mle.full[1]; tau2.mle.full <- exp(mle.full[2])
##         LR <- 2 * (loglik(mu.mle.full,tau2.mle.full,x,params,GH) - loglik(mu.null,tau2.mle.redu,x,params,GH))
##         return(list(mle.full=c(mu=mu.mle.full,tau2=tau2.mle.full),mle.redu=c(mu=mu.null,tau2=tau2.mle.redu), LR=LR))
## }
cond.pmf.hyperg <- function(x,theta,mu,tau2,params) {
    m.T <- params[,'m.T']; m.C <- params[,'m.C']; s <- params[,'s']
    numer <- choose(m.T,x)*choose(m.C,s-x)*exp((mu+sqrt(tau2)*theta)*x)
    supp <- max(0,s-m.C):min(s,m.T)
    denom <- rowSums(  outer(theta,supp,function(thetaj,xj)choose(m.T,xj)*choose(m.C,s-xj)*exp((mu+sqrt(tau2)*thetaj)*xj) ) )
    numer / denom
}
cond.pmf.hyperg <- function(x,theta,mu,tau2,params) {
    ## logpmf.raw <- function(xi,theta)lgamma(m.Ti+1) - lgamma(xi+1) -lgamma(m.Ti-xi+1) + lgamma(m.Ci+1) - lgamma(si-xi+1) - lgamma(m.Ci-si+xi+1) + (mu+sqrt(tau2)*theta)*xi
    ## lognumer <- logpmf.raw(xi,theta)
    ## logpmf.mat <- outer(0:m.Ti,theta,logpmf.raw) # speed this up
    m.T <- params[,'m.T']; m.C <- params[,'m.C']; s <- params[,'s']
    logpmf.part <- function(xi)lgamma(m.T+1) - lgamma(xi+1) -lgamma(m.T-xi+1) + lgamma(m.C+1) - lgamma(s-xi+1) - lgamma(m.C-s+xi+1)
    lognumer <- logpmf.part(xi)+ (mu+sqrt(tau2)*theta)*xi
    ## logpmf.mat <- logpmf.part(0:m.T) + outer(0:m.T,theta)*sqrt(tau2) + (0:m.T)*mu
    supp <- max(0,s-m.C):min(s,m.T)
    logpmf.mat <- logpmf.part(supp) + outer(supp,theta)*sqrt(tau2) + supp*mu
    colmaxes <- apply(logpmf.mat,2,max)        
    logdenom <- colmaxes + log(rowSums(exp(t(logpmf.mat) - colmaxes)))
    exp(lognumer - logdenom)
}
loglik.unit <- function(mu,tau2,xi,paramsi,cond.pmf,GH=NULL) {
    ## browser()
    if(is.null(GH)) {
        likelihood <- integrate(function(theta)cond.pmf(xi,theta,mu,tau2,paramsi)*dnorm(theta),-Inf,Inf)$val 
    } else {
        likelihood <- sum(GH$w*cond.pmf(xi,GH$x,mu,tau2,paramsi))
    }
    return(log(likelihood))
}
set.seed(1)
mu <- runif(1); tau2 <- rexp(1); xi <- rpois(1,10); paramsi <- cbind(m.T=20,m.C=15,s=xi+5)
GH <- fastGHQuad::gaussHermiteData(30)
GH$x <- sqrt(2)*GH$x; GH$w <- GH$w/sqrt(pi)
## loglik.unit.hyperg(mu,tau2,xi,paramsi,GH=GH)
loglik.unit(mu,tau2,xi,paramsi,cond.pmf.hyperg,GH=GH)

## 4c. diagnostic--check mle's -- fine
source('utils.R')
n <- 10
mu.true <- mu.null <- 0
tau2.true <- tau2.null <- 1/2
## n <- 10
ns <- 5:15
## by.n <- sapply(ns, function(n) {
m.T <- rep(20,n)
m.C <- rep(15,n)
s <- sapply(1:n,function(i)sample(max(m.T[i],m.C[i]),1))
mles <- replicate(1e2, {
    ## theta <- rnorm(n)
    ## odds <- exp(mu.true+sqrt(tau2.true)*theta)
    ## x <- x.T <- sapply(1:n, function(i) BiasedUrn::rFNCHypergeo(1,m1=m.T[i],m2=m.C[i],n=s[i],odds=odds[i]))
    params <- cbind(m.T=m.T,m.C=m.C,s=s)
    x <- dgp.hyperg(n,mu.null,tau2.null,params)
    LR.gaussian(x,params,mu.null,cond.pmf.hyperg)$mle.full
})
## ## metafor. mle's for mu match.
## ai <- x.T
## bi <- m.T - x.T
## ci <- s - x.T
## di <- m.C - ci
## fit <- metafor::rma.glmm(
##   ai = ai, bi = bi, 
##   ci = ci, di = di, 
##   measure = "OR", 
##   model = "CM.EL", 
##   method = "ML"
## )
## summary(fit)
op <- par(mfrow=c(1,2))
hist(mles['mu',])
abline(v=mu.true,col=2)
abline(v=mean(mles['mu',]),col=3)
hist(mles['tau2',])
abline(v=tau2.true,col=2)
abline(v=mean(mles['tau2',]),col=3)
par(op)


## 4d. diagnostic--check LR rate--fine
## set.seed(1)
source('utils.R')
for(j in 1:20) {
    mu.true <- mu.null <- 0
    tau2.true <- tau2.null <- 1/2
    ## n <- 10
    ns <- 5:15
    ## by.n <- sapply(ns, function(n) {
    by.n <- parallel::mclapply(ns,mc.cores=4,FUN= function(n) {
        cat('.')
        m.T <- rep(20,n)
        m.C <- rep(15,n)
        s <- sapply(1:n,function(i)sample(max(m.T[i],m.C[i]),1))
        ## s <- rep(10,n)
        stats <- replicate(round(n^2.3), {
            theta <- rnorm(n)
            odds <- exp(mu.true+sqrt(tau2.true)*theta)
            x <- x.T <- sapply(1:n, function(i) BiasedUrn::rFNCHypergeo(1,m1=m.T[i],m2=m.C[i],n=s[i],odds=odds[i]))
            params <- cbind(m.T=m.T,m.C=m.C,s=s)
            LR <- LR.gaussian(x,params,mu.null,cond.pmf.hyperg)$LR
        })
        ## stats <- as.matrix(stats,nrow=1)
        ecdf.LR <- ecdf(pchisq(stats,df=1,lower=FALSE))
        ks <- max(abs(ecdf.LR(knots(ecdf.LR)) - knots(ecdf.LR)))
    })
    by.n <- simplify2array(by.n)
    plot(ns,by.n)
    print(coef(lm(log(by.n) ~ log(ns)))[2])
}



## 4e. diagnostic -- check z-statistics to verify k formula for information -- these
## look good. skew for tau2 is very strong there, need n >= about 40
## for match.
n <- 40
mu.true <- mu.null <- 0
tau2.true <- tau2.null <- 1/2
inf.22 <- -k.df[k.df[,'inner']=='22'  & k.df[,'outer']=='0','k.val']
inf.11 <- -k.df[k.df[,'inner']=='11'  & k.df[,'outer']=='0','k.val']
inf.12 <- -k.df[k.df[,'inner']=='12'  & k.df[,'outer']=='0','k.val']
inf <- matrix(c(inf.11,inf.12,inf.12,inf.22),nrow=2)
vcov.hat <- solve(inf)
## n <- 10
ns <- 5:15
## by.n <- sapply(ns, function(n) {
m.T <- rep(20,n)
m.C <- rep(15,n)
s <- sapply(1:n,function(i)sample(max(m.T[i],m.C[i]),1))
z.stats <- replicate(1e2, {
    theta <- rnorm(n)
    odds <- exp(mu.true+sqrt(tau2.true)*theta)
    x <- x.T <- sapply(1:n, function(i) BiasedUrn::rFNCHypergeo(1,m1=m.T[i],m2=m.C[i],n=s[i],odds=odds[i]))
    params <- cbind(m.T=m.T,m.C=m.C,s=s)
    ## mu.hat <- LR.gaussian(x,params,mu.null,loglik.unit.hyperg)$mle.full['mu']
    ## (mu.hat - mu.true) / sqrt(vcov.hat[2,2]/n)
    ## mles <- rev(LR.gaussian(x,params,mu.null,loglik.unit.hyperg)$mle.full)
    mles <- rev(LR.gaussian(x,params,mu.null,cond.pmf.hyperg)$mle.full)
    (mles-c(tau2.true,mu.true)) / sqrt(diag(vcov.hat)/n)
})
op <- par(mfrow=c(1,2))
qqnorm(z.stats['mu',]); abline(0,1)
qqnorm(z.stats['tau2',]); abline(0,1)
par(op)



## 4f. bartlett correction for the hyperg-normal model

source('utils.R')
## m <- 20
mu.null <- 0
tau2.null <- 1/2
beta <- c(tau2.null,mu.null)

ks.fun <- readRDS('ks.fun.rds')
k.orders <- lapply(ks.fun,function(k.fun)k.fun$order)

## ## read in pmf derivatives--using generic params
## source('utils.R')
## mma.str.path <- './mma-formulas/pmf-hypergOLD/'
## filelist <- dir(mma.str.path)
## filelist <- paste0(mma.str.path,filelist[grep('\\.txt',filelist)])
## ## file <- './mma-formulas/quad/g/g-10.txt'
## pmf.derivs.fun <- lapply(filelist, function(file) {
##     mma.str <-paste(readLines(file,warn=FALSE),collapse='')
##     ## mma.str <- gsub('Pi','pi',mma.str)
##     mma.str <- gsub('Inactive[Sum]','Sum',mma.str,fixed=TRUE)
##     mma.str <- gsub('mC','m.C',mma.str,fixed=TRUE)
##     mma.str <- gsub('mT','m.T',mma.str,fixed=TRUE)
##     mma.str <-  mma.to.R(mma.str,
##                          funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial","Hypergeometric2F1"),
##                          rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose",Hypergeometric2F1="gsl::hyperg_2F1"))
##     g.deriv <- 
##         str.to.fun(mma.str,
##                    args = c('beta','theta','x','params'),
##                    binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt))
## })
## pmf.deriv.orders <- paste0('pmfSum',substr(filelist,nchar(filelist)-5,nchar(filelist)-4))
## names(pmf.derivs.fun) <- pmf.deriv.orders
## pmf.derivs.fun.ma <- lapply(pmf.derivs.fun,function(pmf.deriv.fun)function(theta,x,params)pmf.deriv.fun(beta,theta,x,params))
## ## need to vectorize these generically
## pmf.derivs.fun.ma <- lapply(pmf.derivs.fun.ma,function(pmf.deriv.fun.ma) Vectorize(pmf.deriv.fun.ma,c('theta','x')))

## read in pmf derivatives -- using pmf recurrences in mma
source('utils.R')
mma.str.path <- './mma-formulas/pmf-hyperg/'
filelist <- dir(mma.str.path)
filelist <- paste0(mma.str.path,filelist[grep('pmf-[[:digit:]][[:digit:]]\\.txt',filelist)])
## file <- './mma-formulas/quad/g/g-10.txt'
## file <- filelist[1]
pmf.derivs.fun <- lapply(filelist, function(file) {
    mma.str <-paste(readLines(file,warn=FALSE),collapse='')
    mma.str <-  mma.to.R(mma.str,
                         funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial","factorial","factMoment","pmf"),
                         rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose",factorial="factorial"))
    ## only.theta.args <- regmatches(mma.str,gregexpr('thetaPart[[:digit:]][[:digit:]]',mma.str))[[1]]
    ## only.theta.args <- unique(only.theta.args)
    pmf <- function(x,params,odds) {
        m.T <- params['m.T']
        m.C <- params['m.C']
        s <- params['s']
        BiasedUrn::dFNCHypergeo(x,m1=m.T,m2=m.C,n=s,odds=odds)        
    }
    factMoment <- function(k,params,odds) {
        ## print(params)
        m.T <- params['m.T']
        m.C <- params['m.C']
        s <- params['s']
        supp.x <- max(0,s-m.C):min(s,m.T)
        sum ( apply(outer(supp.x, 0:(k-1), `-`), 1, prod) * pmf(supp.x,params,odds) )
    }
    pmf.deriv <- 
        str.to.fun(mma.str,
                   args = c('beta','params','theta','x'),
                   ## args = c('idx','x','theta','beta','params'),
                   binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt,pmf=pmf,factMoment=factMoment))
    ## })
})
pmf.deriv.orders <- paste0('pmfSum',substr(filelist,nchar(filelist)-5,nchar(filelist)-4))
names(pmf.derivs.fun) <- pmf.deriv.orders
## saveRDS(pmf.derivs.fun,file='pmf-derivs-fun-hyperg.rds',compress=FALSE)

m <- m.T <- 20; m.C <- 15; s <- 10;
support <- max(0,s-m.C):min(s,m.T)
k.vals <- ks.study(pmf.derivs.fun.ma,ks.fun,params=c(m.C=m.C,m.T=m.T,s=s),support=support)
k.df <- data.frame(inner=sapply(k.orders,`[`,1),outer=sapply(k.orders,`[`,2),k.val=k.vals)
bartlett.C <- lawley(k.df,2) - lawley(k.df,1)
print(bartlett.C)




## 4g. sim template

source('utils.R')
## check bartlett rates        
## set.seed(1)
mu.true <- mu.null <- 0
tau2.true <- tau2.null <- 1/2
## n <- 10
ns <- 5:15
bartlett.C <- 1.36433398617429
tic <- Sys.time()
## by.n <- sapply(ns, function(n) {
repeat{
    by.n <- parallel::mclapply(ns,mc.cores=4,FUN= function(n) {
        cat('.')
        m.T <- rep(20,n)
        m.C <- rep(15,n)
        s <- rep(10,n)
        ## s <- rep(10,n)
        stats <- replicate(round(n^2.3), {
            theta <- rnorm(n)
            odds <- exp(mu.true+sqrt(tau2.true)*theta)
            x <- x.T <- sapply(1:n, function(i) BiasedUrn::rFNCHypergeo(1,m1=m.T[i],m2=m.C[i],n=s[i],odds=odds[i]))
            params <- cbind(m.T=m.T,m.C=m.C,s=s)
            LR <- LR.gaussian(x,params,mu.null,cond.pmf.hyperg)$LR
            ## LR <- LR.gaussian(x,params,mu.null,loglik.unit.hyperg)$LR
            LR.bartlett <- LR / (1+bartlett.C/n)
            c(LR=LR,LR.bartlett=LR.bartlett)
        })
        ecdfs <- apply(stats,1,function(x)ecdf(pchisq(x,df=1,lower=FALSE)))
        ks <- sapply(ecdfs,function(ecdf)max(abs(ecdf(knots(ecdf)) - knots(ecdf))))
    })
    by.n <- simplify2array(by.n)
    print(Sys.time() - tic)
    plot(ns,by.n['LR',],type='l',ylim=range(by.n))
    lines(ns,by.n['LR.bartlett',],col=2)
    apply(by.n,1,function(ks)coef(lm(log(ks) ~ log(ns)))[2]) |> print()
}




## 5. poisson-normal model

## 5a. diagnostics

## check condpmf.poisson gives consistent mle's --fine
source('utils.R')
n <- 10
mu.null <- runif(1)
tau2.null <- runif(1,1/2,2)
## loglik.unit.poisson <- function(mu,tau2,xi,paramsi,GH=NULL)loglik.unit(mu,tau2,xi,paramsi,cond.pmf.poisson,GH)
mles <- replicate(1e2, {
    theta <- rnorm(n)
    T <- 4 + rexp(n)
    lambda <- T*exp(mu.null+sqrt(tau2.null)*theta)
    x <- rpois(n,lambda)
    params <- cbind(T=T)
    mles <- LR.gaussian(x,params,mu.null,cond.pmf.poisson)$mle.full
    ## LR.gaussian(x,params,mu.null,loglik.unit.poisson)$mle.full
})
op <- par(mfrow=c(1,2))
hist(mles['mu',])
abline(v=mu.null,col=2)
abline(v=mean(mles['mu',]),col=3)
hist(mles['tau2',])
abline(v=tau2.null,col=2)
abline(v=mean(mles['tau2',]),col=3)
par(op)




## 5b. bartlett correction 

source('utils.R')
## m <- 20
mu.null <- 0
tau2.null <- 1/2
beta <- c(tau2.null,mu.null)

ks.fun <- readRDS('ks.fun.rds')
k.orders <- lapply(ks.fun,function(k.fun)k.fun$order)

## ## read in pmf derivatives
## mma.str.path <- './mma-formulas/pmf-poisson/'
## filelist <- dir(mma.str.path)
## filelist <- paste0(mma.str.path,filelist[grep('\\.txt',filelist)])
## ## file <- './mma-formulas/quad/g/g-10.txt'
## pmf.derivs.fun <- lapply(filelist, function(file) {
##     mma.str <-paste(readLines(file,warn=FALSE),collapse='')
##     ## mma.str <- gsub('Pi','pi',mma.str)
##     ## mma.str <- gsub('Inactive[Sum]','Sum',mma.str,fixed=TRUE)
##     ## mma.str <- gsub('mC','m.C',mma.str,fixed=TRUE)
##     ## mma.str <- gsub('mT','m.T',mma.str,fixed=TRUE)
##     mma.str <-  mma.to.R(mma.str,
##                          funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial","factorial"),
##                          rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose",factorial="factorial"))
##     g.deriv <- 
##         str.to.fun(mma.str,
##                    args = c('beta','theta','x','params'),
##                    binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt))
## })
## pmf.deriv.orders <- paste0('pmfSum',substr(filelist,nchar(filelist)-5,nchar(filelist)-4))
## names(pmf.derivs.fun) <- pmf.deriv.orders
## saveRDS(pmf.derivs.fun,file='pmf-derivs-fun-poisson.rds',compress=FALSE)
## pmf.derivs.fun.ma <- lapply(pmf.derivs.fun,function(pmf.deriv.fun)function(theta,x,params)pmf.deriv.fun(beta,theta,x,params))

## read in pmf derivatives, after recursive formulation (using dpois)
mma.str.path <- './mma-formulas/pmf-poisson/'
filelist <- dir(mma.str.path)
filelist <- paste0(mma.str.path,filelist[grep('\\.txt',filelist)])
## file <- './mma-formulas/quad/g/g-10.txt'
pmf.derivs.fun <- lapply(filelist, function(file) {
    mma.str <-paste(readLines(file,warn=FALSE),collapse='')
    ## mma.str <- gsub('Pi','pi',mma.str)
    ## mma.str <- gsub('Inactive[Sum]','Sum',mma.str,fixed=TRUE)
    ## mma.str <- gsub('mC','m.C',mma.str,fixed=TRUE)
    ## mma.str <- gsub('mT','m.T',mma.str,fixed=TRUE)
    mma.str <-  mma.to.R(mma.str,
                         funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial","factorial","dpois"),
                         rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose",factorial="factorial",dpois="dpois"))
    g.deriv <- 
        str.to.fun(mma.str,
                   args = c('beta','theta','x','params'),
                   binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt,dpois=dpois))
})
pmf.deriv.orders <- paste0('pmfSum',substr(filelist,nchar(filelist)-5,nchar(filelist)-4))
names(pmf.derivs.fun) <- pmf.deriv.orders
## saveRDS(pmf.derivs.fun,file='pmf-derivs-fun-poisson.rds',compress=FALSE)
## pmf.derivs.fun.ma <- lapply(pmf.derivs.fun,function(pmf.deriv.fun)function(theta,x,params)pmf.deriv.fun(beta,theta,x,params))


T <- 1
k.vals <- ks.study(pmf.derivs.fun.ma,ks.fun,params=c(T=T),support=0:50)
k.df <- data.frame(inner=sapply(k.orders,`[`,1),outer=sapply(k.orders,`[`,2),k.val=k.vals)
bartlett.C <- lawley(k.df,2) - lawley(k.df,1)
print(bartlett.C)




## 5c. diagnostic -- check z-statistics -- look fine for large enough n
k.22 <- k.df[k.df[,'inner']=='22'  & k.df[,'outer']=='0','k.val']
k.11 <- k.df[k.df[,'inner']=='11'  & k.df[,'outer']=='0','k.val']
k.12 <- k.df[k.df[,'inner']=='12'  & k.df[,'outer']=='0','k.val']
n <- 1e2
mu.true <- mu.null <- 0
tau2.true <- tau2.null <- 1/2
## k.11 <-  -0.3897629 # information in tau2
## k.22 <- -0.7201051 # information in mu
## k.12 <- -0.2889417
inf <- -matrix(c(k.11,k.12,k.12,k.22),nrow=2)
vcov.hat <- solve(inf)
T <- rep(1,n)
z.stats <- replicate(1e2, {
    theta <- rnorm(n)
    ## T <- 4 + rexp(n)
    lambda <- T*exp(mu.null+sqrt(tau2.null)*theta)
    x <- rpois(n,lambda)
    params <- cbind(T=T)
    mles <- LR.gaussian(x,params,mu.null,cond.pmf.poisson)$mle.full
    (rev(mles)-c(tau2.true,mu.true)) / sqrt(diag(vcov.hat)/n)
})
op <- par(mfrow=c(1,2))
qqnorm(z.stats['mu',]); abline(0,1)
qqnorm(z.stats['tau2',]); abline(0,1)
par(op)


## 5d. sim template

source('utils.R')
## check bartlett rates        
## set.seed(1)
mu.true <- mu.null <- 0
tau2.true <- tau2.null <- 1/2
## n <- 10
ns <- 5:15
## bartlett.C <- 1.36433398617429 # hyperg
bartlett.C <- 1.17164232232645 # poisson
tic <- Sys.time()
## by.n <- sapply(ns, function(n) {
repeat{
    by.n <- parallel::mclapply(ns,mc.cores=4,FUN= function(n) {
        cat('.')
        ## m.T <- rep(20,n)
        ## m.C <- rep(15,n)
        ## s <- rep(10,n)
        ## s <- rep(10,n)
        T <- rep(1,n)
        stats <- replicate(round(n^2.3), {
            theta <- rnorm(n)
            ## T <- 4 + rexp(n)
            lambda <- T*exp(mu.null+sqrt(tau2.null)*theta)
            x <- rpois(n,lambda)
            params <- cbind(T=T)
            ## theta <- rnorm(n)
            ## odds <- exp(mu.true+sqrt(tau2.true)*theta)
            ## x <- x.T <- sapply(1:n, function(i) BiasedUrn::rFNCHypergeo(1,m1=m.T[i],m2=m.C[i],n=s[i],odds=odds[i]))
            ## params <- cbind(m.T=m.T,m.C=m.C,s=s)
            LR <- LR.gaussian(x,params,mu.null,cond.pmf.poisson)$LR
            ## LR <- LR.gaussian(x,params,mu.null,loglik.unit.hyperg)$LR
            LR.bartlett <- LR / (1+bartlett.C/n)
            c(LR=LR,LR.bartlett=LR.bartlett)
        })
        ecdfs <- apply(stats,1,function(x)ecdf(pchisq(x,df=1,lower=FALSE)))
        ks <- sapply(ecdfs,function(ecdf)max(abs(ecdf(knots(ecdf)) - knots(ecdf))))
    })
    by.n <- simplify2array(by.n)
    print(Sys.time() - tic)
    plot(ns,by.n['LR',],type='l',ylim=range(by.n))
    lines(ns,by.n['LR.bartlett',],col=2)
    apply(by.n,1,function(ks)coef(lm(log(ks) ~ log(ns)))[2]) |> print()
}


## 6. redo single-arm binomial with generic code


## 6a. diagnostics

## check condpmf.poisson gives consistent mle's --fine
source('utils.R')
n <- 10
mu.null <- runif(1)
tau2.null <- runif(1,1/2,2)
## m <- 20
## loglik.unit.poisson <- function(mu,tau2,xi,paramsi,GH=NULL)loglik.unit(mu,tau2,xi,paramsi,cond.pmf.poisson,GH)
mles <- replicate(1e2, {
    ## theta <- rnorm(n)
    ## T <- 4 + rexp(n)
    ## lambda <- T*exp(mu.null+sqrt(tau2.null)*theta)
    ## x <- rpois(n,lambda)
    ## params <- cbind(T=T)
    theta <- rnorm(n)
    pr <- plogis(mu.null+sqrt(tau2.null)*theta)
    m <- rpois(n,lambda=20)
    x <- rbinom(n,m,pr)
    params <- cbind(m=m)
    mles <- LR.gaussian(x,params,mu.null,cond.pmf.binomial)$mle.full
    ## LR.gaussian(x,params,mu.null,loglik.unit.poisson)$mle.full
})
op <- par(mfrow=c(1,2))
hist(mles['mu',])
abline(v=mu.null,col=2)
abline(v=mean(mles['mu',]),col=3)
hist(mles['tau2',])
abline(v=tau2.null,col=2)
abline(v=mean(mles['tau2',]),col=3)
par(op)




## 6b. bartlett correction 

source('utils.R')
## m <- 20
mu.null <- 0
tau2.null <- 1/2
beta <- c(tau2.null,mu.null)

ks.fun <- readRDS('ks.fun.rds')
k.orders <- lapply(ks.fun,function(k.fun)k.fun$order)

## read in pmf derivatives
mma.str.path <- './mma-formulas/pmf-binomial/'
filelist <- dir(mma.str.path)
filelist <- paste0(mma.str.path,filelist[grep('\\.txt',filelist)])
## file <- './mma-formulas/quad/g/g-10.txt'
pmf.derivs.fun <- lapply(filelist, function(file) {
    mma.str <-paste(readLines(file,warn=FALSE),collapse='')
    ## mma.str <- gsub('Pi','pi',mma.str)
    ## mma.str <- gsub('Inactive[Sum]','Sum',mma.str,fixed=TRUE)
    ## mma.str <- gsub('mC','m.C',mma.str,fixed=TRUE)
    ## mma.str <- gsub('mT','m.T',mma.str,fixed=TRUE)
    mma.str <-  mma.to.R(mma.str,
                         funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial","factorial"),
                         rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose",factorial="factorial"))
    g.deriv <- 
        str.to.fun(mma.str,
                   args = c('beta','theta','x','params'),
                   binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt))
})
pmf.deriv.orders <- paste0('pmfSum',substr(filelist,nchar(filelist)-5,nchar(filelist)-4))
names(pmf.derivs.fun) <- pmf.deriv.orders
## saveRDS(pmf.derivs.fun,file='pmf-derivs-fun-binomial.rds',compress=FALSE)
pmf.derivs.fun.ma <- lapply(pmf.derivs.fun,function(pmf.deriv.fun)function(theta,x,params)pmf.deriv.fun(beta,theta,x,params))

m <- 20
k.vals <- ks.study(pmf.derivs.fun.ma,ks.fun,params=c(m=m),support=0:m)
k.df <- data.frame(inner=sapply(k.orders,`[`,1),outer=sapply(k.orders,`[`,2),k.val=k.vals)
bartlett.C <- lawley(k.df,2) - lawley(k.df,1)
print(bartlett.C)





## 6c. diagnostic -- check z-statistics -- look fine 
k.22 <- k.df[k.df[,'inner']=='22'  & k.df[,'outer']=='0','k.val']
k.11 <- k.df[k.df[,'inner']=='11'  & k.df[,'outer']=='0','k.val']
k.12 <- k.df[k.df[,'inner']=='12'  & k.df[,'outer']=='0','k.val']
n <- 1e2
mu.true <- mu.null <- 0
tau2.true <- tau2.null <- 1/2
## k.11 <-  -0.3897629 # information in tau2
## k.22 <- -0.7201051 # information in mu
## k.12 <- -0.2889417
inf <- -matrix(c(k.11,k.12,k.12,k.22),nrow=2)
vcov.hat <- solve(inf)
T <- rep(1,n)
z.stats <- replicate(1e2, {
    theta <- rnorm(n)
    pr <- plogis(mu.null+sqrt(tau2.null)*theta)
    m <- rpois(n,lambda=20)
    x <- rbinom(n,m,pr)
    params <- cbind(m=m)
    mles <- LR.gaussian(x,params,mu.null,cond.pmf.binomial)$mle.full
    (rev(mles)-c(tau2.true,mu.true)) / sqrt(diag(vcov.hat)/n)
})
op <- par(mfrow=c(1,2))
qqnorm(z.stats['mu',]); abline(0,1)
qqnorm(z.stats['tau2',]); abline(0,1)
par(op)







## 7. binomial--incidence rate ratio model

## 7a. diagnostics
## check condpmf.poisson gives consistent mle's --fine
source('utils.R')
n <- 10
mu.null <- runif(1)
tau2.null <- runif(1,1/2,2)
## T.t <- runif(n)+.5#rep(1,n)
## T.c <- runif(n)+.5#rep(1.5,n)
## s <- rpois(n,lambda=10)
## params <- cbind(T.t=T.t,T.c=T.c,s=s)
mles <- replicate(1e2, {
    T.t <- runif(n)+.5#rep(1,n)
    T.c <- runif(n)+.5#rep(1.5,n)
    s <- rpois(n,lambda=10)
    params <- cbind(T.t=T.t,T.c=T.c,s=s)
    x <- dgp.binomial.irr(n,mu.null,tau2.null,params)
    mles <- LR.gaussian(x,params,mu.null,cond.pmf.binomial.irr)$mle.full
    ## LR.gaussian(x,params,mu.null,loglik.unit.poisson)$mle.full
})
op <- par(mfrow=c(1,2))
hist(mles['mu',])
abline(v=mu.null,col=2)
abline(v=mean(mles['mu',]),col=3)
hist(mles['tau2',])
abline(v=tau2.null,col=2)
abline(v=mean(mles['tau2',]),col=3)
par(op)



## 7b. bartlett correction 

source('utils.R')
## m <- 20
mu.null <- 0
tau2.null <- 1/2
beta <- c(tau2.null,mu.null)

ks.fun <- readRDS('ks.fun.rds')
k.orders <- lapply(ks.fun,function(k.fun)k.fun$order)

## ## read in pmf derivatives
## mma.str.path <- './mma-formulas/pmf-binomial-irr/'
## filelist <- dir(mma.str.path)
## filelist <- paste0(mma.str.path,filelist[grep('\\.txt',filelist)])
## ## file <- './mma-formulas/quad/g/g-10.txt'
## pmf.derivs.fun <- lapply(filelist, function(file) {
##     mma.str <-paste(readLines(file,warn=FALSE),collapse='')
##     mma.str <-  mma.to.R(mma.str,
##                          funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial","factorial"),
##                          rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose",factorial="factorial"))
##     g.deriv <- 
##         str.to.fun(mma.str,
##                    args = c('beta','theta','x','params'),
##                    binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt))
## })
## pmf.deriv.orders <- paste0('pmfSum',substr(filelist,nchar(filelist)-5,nchar(filelist)-4))
## names(pmf.derivs.fun) <- pmf.deriv.orders
## ## saveRDS(pmf.derivs.fun,file='pmf-derivs-fun-binomial-irr.rds',compress=FALSE)
## pmf.derivs.fun.ma <- lapply(pmf.derivs.fun,function(pmf.deriv.fun)function(theta,x,params)pmf.deriv.fun(beta,theta,x,params))


## read in pmf derivatives -- using pmf recurrences in mma
source('utils.R')
mma.str.path <- './mma-formulas/pmf-binomial-irr/'
filelist <- dir(mma.str.path)
filelist <- paste0(mma.str.path,filelist[grep('pmf-[[:digit:]][[:digit:]]\\.txt',filelist)])
pmf.derivs.fun <- lapply(filelist, function(file) {
    mma.str <-paste(readLines(file,warn=FALSE),collapse='')
    mma.str <-  mma.to.R(mma.str,
                         funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial","factorial","factMoment","pmf"),
                         rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose",factorial="factorial"))
    pmf <- function(x,params,prob) {
        s <- params['s']
        dbinom(x=x,size=s,prob=prob)
    }
    ## factMoment <- function(k,params,odds) {
    ##     ## print(params)
    ##     m.T <- params['m.T']
    ##     m.C <- params['m.C']
    ##     s <- params['s']
    ##     supp.x <- max(0,s-m.C):min(s,m.T)
    ##     sum ( apply(outer(supp.x, 0:(k-1), `-`), 1, prod) * pmf(supp.x,params,odds) )
    ## }
    pmf.deriv <- 
        str.to.fun(mma.str,
                   args = c('beta','params','theta','x'),
                   ## args = c('idx','x','theta','beta','params'),
                   binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt,pmf=pmf))
    ## })
})
pmf.deriv.orders <- paste0('pmfSum',substr(filelist,nchar(filelist)-5,nchar(filelist)-4))
names(pmf.derivs.fun) <- pmf.deriv.orders
## saveRDS(pmf.derivs.fun,file='pmf-derivs-fun-binomial-irr.rds',compress=FALSE)




s <- 30
T.t <- 1
T.c <- 1.5
params <- c(s=s,T.t=T.t,T.c=T.c)
support <- 0:s
ks.df <- ks.study(pmf.derivs.fun,ks.fun,params=params,supp.x=support,mu.null=mu.null,tau2.null=tau2.null)
## k.df <- data.frame(inner=sapply(k.orders,`[`,1),outer=sapply(k.orders,`[`,2),k.val=k.vals)
bartlett.C <- lawley(ks.df,2) - lawley(ks.df,1)
print(bartlett.C)
## 1.36297164231076





## 7c. diagnostic -- check z-statistics -- look fine 
k.22 <- k.df[k.df[,'inner']=='22'  & k.df[,'outer']=='0','k.val']
k.11 <- k.df[k.df[,'inner']=='11'  & k.df[,'outer']=='0','k.val']
k.12 <- k.df[k.df[,'inner']=='12'  & k.df[,'outer']=='0','k.val']
n <- 20
mu.true <- mu.null <- 0
tau2.true <- tau2.null <- 1/2
## k.11 <-  -0.3897629 # information in tau2
## k.22 <- -0.7201051 # information in mu
## k.12 <- -0.2889417
inf <- -matrix(c(k.11,k.12,k.12,k.22),nrow=2)
vcov.hat <- solve(inf)
s <- rep(10,n)
T.t <- rep(1,n)
T.c <- rep(1.5,n)
params <- cbind(s=s,T.t=T.t,T.c=T.c)
z.stats <- replicate(1e2, {
    x <- dgp.binomial.irr(n,mu.null,tau2.null,params)
    mles <- LR.gaussian(x,params,mu.null,cond.pmf.binomial.irr)$mle.full
    (rev(mles)-c(tau2.true,mu.true)) / sqrt(diag(vcov.hat)/n)
})
op <- par(mfrow=c(1,2))
qqnorm(z.stats['mu',]); abline(0,1)
qqnorm(z.stats['tau2',]); abline(0,1)
par(op)





## 8 speed up hypergeom calculations

## 8a. separate calculation of normalizing constant and
## derivatives. only depends on tehta, x is constant.

source('utils.R')
## read in theta-only components
mma.str.path <- './mma-formulas/pmf-hyperg/'
filelist <- dir(mma.str.path)
filelist <- paste0(mma.str.path,filelist[grep('theta-part-[[:digit:]][[:digit:]]\\.txt',filelist)])
## file <- filelist[1]
theta.parts.fun <- lapply(filelist, function(file) {
    mma.str <-paste(readLines(file,warn=FALSE),collapse='')
    mma.str <-  mma.to.R(mma.str,
                         funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial","factorial"),
                         rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose",factorial="factorial"))
    theta.part <- 
        str.to.fun(mma.str,
                   args = c('beta','theta','x','params'),
                   binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt))
})
orders.theta.parts <- substr(filelist,nchar(filelist)-5,nchar(filelist)-4)
names(theta.parts.fun) <- paste0('thetaPart',orders.theta.parts)
## saveRDS(theta.parts.fun,file='theta-parts-fun-hyperg.rds',compress=FALSE)

## read in pmf derivatives
mma.str.path <- './mma-formulas/pmf-hyperg/'
filelist <- dir(mma.str.path)
filelist <- paste0(mma.str.path,filelist[grep('pmf-[[:digit:]][[:digit:]]\\.txt',filelist)])
## file <- './mma-formulas/quad/g/g-10.txt'
## file <- filelist[1]
pmf.derivs.fun <- lapply(filelist, function(file) {
    mma.str <-paste(readLines(file,warn=FALSE),collapse='')
    mma.str <-  mma.to.R(mma.str,
                         funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial","factorial"),
                         rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose",factorial="factorial"))
    only.theta.args <- regmatches(mma.str,gregexpr('thetaPart[[:digit:]][[:digit:]]',mma.str))[[1]]
    only.theta.args <- unique(only.theta.args)
    g.deriv <- 
        str.to.fun(mma.str,
                   args = c('beta','params','theta','x',only.theta.args),
                   ## args = c('idx','x','theta','beta','params'),
                   binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt))
    ## })
})
pmf.deriv.orders <- paste0('pmfSum',substr(filelist,nchar(filelist)-5,nchar(filelist)-4))
pmf.deriv.args <- sapply(pmf.derivs.fun,function(f) {
    vars <- all.vars(body(f))
    vars <- unique(vars[grep('thetaPart',vars)])
    ## substr(vars,nchar(vars)-1,nchar(vars))
})
pmf.derivs.fun <- lapply(1:length(pmf.derivs.fun), function(i)structure(pmf.derivs.fun[[i]],args=unname(pmf.deriv.args[[i]])))
names(pmf.derivs.fun) <- pmf.deriv.orders
## saveRDS(pmf.derivs.fun,file='pmf-derivs-fun-hyperg.rds',compress=FALSE)


## pmf.derivs.fun.ma <- lapply(pmf.derivs.fun,function(pmf.deriv.fun)function(theta,x,params)pmf.deriv.fun(beta,theta,x,params))

ks.fun <- readRDS('ks-fun.rds')



beta <- c(1/2,0)
m <- m.T <- 20; m.C <- 15; s <- 10;
supp.x <- max(0,s-m.C):min(s,m.T)
params <- c(m.T=m.T,m.C=m.C,s=s)
n.GH <- 30
## with(list(theta=GH$theta.grid),
## f(1:n,1:n,GH$theta.grid,beta,params)
GH <- fastGHQuad::gaussHermiteData(n.GH)
GH$x <- sqrt(2)*GH$x
GH$w <- GH$w / sqrt(pi)
names(GH) <- c('theta.grid','w')

## theta <- GH$theta.grid[1]
## odds <- exp(beta[2]+sqrt(beta[1])*theta)
## x <- supp.x[1]
## BiasedUrn::dFNCHypergeo(x,m1=m.T,m2=m.C,n=s,odds=odds)


## f <- theta.parts[[1]]
theta.parts <- lapply(theta.parts.fun, function(f) {
    theta.part.mat <- outer(GH$theta.grid,supp.x,function(theta,x)f(beta,theta,x,params))
    rowSums(theta.part.mat)
})

## f <- pmf.derivs.fun[[1]]
pmf.derivs.unit <- lapply(pmf.derivs.fun, function(f) {
    arg.list <- c(list(
        theta = rep(GH$theta.grid,length(supp.x)),
        x = rep(supp.x,each=length(GH$theta.grid))),
        lapply(theta.parts[attr(f,'args')], function(v)rep(v,length(supp.x))),
        beta=list(beta),
        params=list(params)
        )
    theta.part.vals <- do.call(f,arg.list)
    dim(theta.part.vals) <- c(length(GH$theta.grid),length(supp.x))
    GH$w %*% theta.part.vals
})


k.orders <- lapply(ks.fun,function(k.fun)k.fun$order)
k.study <- sapply(seq_along(ks.fun), function(i) {
    k.fun <- ks.fun[[i]]$fun; args <- ks.fun[[i]]$args
    k.fun.unit <- do.call(k.fun, pmf.derivs.unit[args])
    k.val <- sum( k.fun.unit)
})
k.vals <- k.study
k.df <- data.frame(inner=sapply(k.orders,`[`,1),outer=sapply(k.orders,`[`,2),k.val=k.vals)
bartlett.C <- lawley(k.df,2) - lawley(k.df,1)
print(bartlett.C)





## 8b. encapsulate new verusion of study-level ks function, allowing for
## theta-only parts of the pmf derivative to be passed in. rework
## ks.fun to store order and args as attributes not list elements.


## k's: this is for all models
mma.str.path <- './mma-formulas/k/'
filelist <- dir(mma.str.path)
filelist <- paste0(mma.str.path,filelist[grep('\\.txt',filelist)])
## file <- './mma-formulas/quad/g/g-10.txt'
ks.fun <- lapply(filelist, function(file) {
    mma.str <-paste(readLines(file,warn=FALSE),collapse='')
    matches.idx <- gregexpr("pmfSum[[:digit:]][[:digit:]]", mma.str)
    variables <- regmatches(mma.str, matches.idx)[[1]] |> unique()
    fun.str <- sprintf(
        "function(%s) { %s }", 
        paste(variables, collapse = ", "), 
        mma.str
    )
    k.fun <- l.deriv.fun <- eval(parse(text = fun.str))
    structure(k.fun,args=variables)
    ## list(fun=k.fun,args=variables)
})
k.orders <- gsub('.txt','',filelist,fixed=TRUE) |> sub(pattern='.*k-',replace='') |> strsplit(split='-')
k.orders <- sapply(k.orders,function(x)if(length(x)==1)c(x,0)else x)
## ks.fun <- lapply(seq_along(ks.fun), function(i)c(ks.fun[[i]],order=list(k.orders[[i]])))
ks.fun <- sapply(seq_along(ks.fun), function(i)c(structure(ks.fun[[i]],order=k.orders[,i])))
saveRDS(ks.fun,file='ks-fun.rds',compress=FALSE)


ks.fun <- readRDS('ks-fun.rds')

ks.study <- function(pmf.derivs.fun,ks.fun,theta.parts.fun=list(),params,supp.x,GH=NULL) {
    ## browser()
    if(is.null(GH)) {
        GH <- fastGHQuad::gaussHermiteData(50)
        GH$x <- sqrt(2)*GH$x
        GH$w <- GH$w / sqrt(pi)
        names(GH) <- c('theta.grid','w')
    }
    theta.parts <- lapply(theta.parts.fun, function(f) {
        theta.part.mat <- outer(GH$theta.grid,supp.x,function(theta,x)f(beta,theta,x,params))
        rowSums(theta.part.mat)
    })
    pmf.derivs.unit <- lapply(pmf.derivs.fun, function(f) {
        arg.list <- c(list(
            theta = rep(GH$theta.grid,length(supp.x)),
            x = rep(supp.x,each=length(GH$theta.grid))),
            lapply(theta.parts[attr(f,'args')], function(v)rep(v,length(supp.x))),
            beta=list(beta),
            params=list(params)
            )
        pmf.deriv.unit <- do.call(f,arg.list)
        dim(pmf.deriv.unit) <- c(length(GH$theta.grid),length(supp.x))
        GH$w %*% pmf.deriv.unit
    })
    k.vals <- sapply(seq_along(ks.fun), function(i) {
        k.fun <- ks.fun[[i]]
        args <- attr(k.fun,'args')
        k.fun.unit <- do.call(k.fun, pmf.derivs.unit[args])
        k.val <- sum( k.fun.unit)
    })
    k.orders <- sapply(ks.fun,function(k.fun)attr(k.fun,'order'))
    ## k.df <- data.frame(inner=sapply(k.orders,`[`,1),outer=sapply(k.orders,`[`,2),k.val=k.vals)
    ks.df <- data.frame(inner=k.orders[1,],outer=k.orders[2,],k.val=k.vals)
    return(ks.df)
}

## 8c using recurrence relations for pmf instead

## 8c-1 compare speeds of outer and sapply over a matrix. seem
## comparable. this needed to know because biasedurn::dhyperg is only
## vectorized over the data.
set.seed(1)
start <- Sys.time()
out <- replicate(1e2, {
    xs <- rpois(500,30)
    lambdas <- rpois(500,30)
    ## uu <- outer(xs,lambdas,function(x,lambda)dpois(x,lambda))
    vv <- sapply(lambdas,function(lambda)dpois(xs,lambda))
})
print(Sys.time()-start)


## 8c-2 read in pmf derivatives, compute bartlett C, check z-stats

## source('utils.R')
## mma.str.path <- './mma-formulas/pmf-hyperg/'
## filelist <- dir(mma.str.path)
## filelist <- paste0(mma.str.path,filelist[grep('pmf-[[:digit:]][[:digit:]]\\.txt',filelist)])
## ## file <- './mma-formulas/quad/g/g-10.txt'
## ## file <- filelist[1]
## pmf.derivs.fun <- lapply(filelist, function(file) {
##     mma.str <-paste(readLines(file,warn=FALSE),collapse='')
##     mma.str <-  mma.to.R(mma.str,
##                          funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial","factorial","factMoment","pmf"),
##                          rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose",factorial="factorial"))
##     ## only.theta.args <- regmatches(mma.str,gregexpr('thetaPart[[:digit:]][[:digit:]]',mma.str))[[1]]
##     ## only.theta.args <- unique(only.theta.args)
##     pmf <- function(x,params,odds) {
##         m.T <- params['m.T']
##         m.C <- params['m.C']
##         s <- params['s']
##         BiasedUrn::dFNCHypergeo(x,m1=m.T,m2=m.C,n=s,odds=odds)        
##     }
##     factMoment <- function(k,params,odds) {
##         ## print(params)
##         m.T <- params['m.T']
##         m.C <- params['m.C']
##         s <- params['s']
##         supp.x <- max(0,s-m.C):min(s,m.T)
##         sum ( apply(outer(supp.x, 0:(k-1), `-`), 1, prod) * pmf(supp.x,params,odds) )
##     }
##     pmf.deriv <- 
##         str.to.fun(mma.str,
##                    args = c('beta','params','theta','x'),
##                    ## args = c('idx','x','theta','beta','params'),
##                    binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt,pmf=pmf,factMoment=factMoment))
##     ## })
## })
## pmf.deriv.orders <- paste0('pmfSum',substr(filelist,nchar(filelist)-5,nchar(filelist)-4))
## names(pmf.derivs.fun) <- pmf.deriv.orders


## updating pmf call so it can apply to either my or biasUrn's version
source('utils.R')
mma.str.path <- './mma-formulas/pmf-hyperg/'
filelist <- dir(mma.str.path)
filelist <- paste0(mma.str.path,filelist[grep('pmf-[[:digit:]][[:digit:]]\\.txt',filelist)])
## file <- './mma-formulas/quad/g/g-10.txt'
## file <- filelist[1]
pmf.derivs.fun <- lapply(filelist, function(file) {
    mma.str <-paste(readLines(file,warn=FALSE),collapse='')
    mma.str <-  mma.to.R(mma.str,
                         funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial","factorial","factMoment","pmf"),
                         rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose",factorial="factorial"))
    ## only.theta.args <- regmatches(mma.str,gregexpr('thetaPart[[:digit:]][[:digit:]]',mma.str))[[1]]
    ## only.theta.args <- unique(only.theta.args)
    pmf <- function(x,params,odds) {
        m.T <- params['m.T']
        m.C <- params['m.C']
        s <- params['s']
        .cond.pmf.hyperg(x=x,m1=m.T,m2=m.C,n=s,odds=odds) 
        ## BiasedUrn::dFNCHypergeo(x,m1=m.T,m2=m.C,n=s,odds=odds)      
    }
    ## factMoment <- function(k,params,odds) {
    ##     ## print(params)
    ##     m.T <- params['m.T']
    ##     m.C <- params['m.C']
    ##     s <- params['s']
    ##     supp.x <- max(0,s-m.C):min(s,m.T)
    ##     sum ( apply(outer(supp.x, 0:(k-1), `-`), 1, prod) * pmf(supp.x,params,odds) )
    ## }
    ## vectorize in odds ratio
    factMoment <- function(k,params,odds) {
        ## print(params)
        m.T <- params['m.T']
        m.C <- params['m.C']
        s <- params['s']
        supp.x <- max(0,s-m.C):min(s,m.T)
        ## sum ( apply(outer(supp.x, 0:(k-1), `-`), 1, prod) * pmf(supp.x,params,odds) )
        apply(outer(supp.x, 0:(k-1), `-`), 1, prod)%*%outer(supp.x,odds,function(xi,oddsi)pmf(xi,params,oddsi)) |> c()
    }
    pmf.deriv <- 
        str.to.fun(mma.str,
                   args = c('beta','params','theta','x'),
                   ## args = c('idx','x','theta','beta','params'),
                   binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt,pmf=pmf,factMoment=factMoment))
    ## })
})
pmf.deriv.orders <- paste0('pmfSum',substr(filelist,nchar(filelist)-5,nchar(filelist)-4))
names(pmf.derivs.fun) <- pmf.deriv.orders
## saveRDS(pmf.derivs.fun,file='pmf-derivs-fun-hyperg.rds',compress=FALSE)


## ks.fun <- readRDS('ks-fun.rds')
## pmf.derivs.fun <- readRDS('pmf-derivs-fun-hyperg.rds')
## beta <- c(1/2,0)
## m <- m.T <- 20; m.C <- 15; s <- 10;
## supp.x <- max(0,s-m.C):min(s,m.T)
## params <- c(m.T=m.T,m.C=m.C,s=s)
## ## params <- c(m.T = 174, m.C = 177, s = 4, x.min = 0, x.max = 4)
## n.GH <- 1e2
## ## with(list(theta=GH$theta.grid),
## ## f(1:n,1:n,GH$theta.grid,beta,params)
## GH <- fastGHQuad::gaussHermiteData(n.GH)
## GH$x <- sqrt(2)*GH$x
## GH$w <- GH$w / sqrt(pi)
## names(GH) <- c('theta.grid','w')
## ## pmf.deriv.fun <- pmf.derivs.fun[[5]]
## pmf.derivs.unit <- lapply(pmf.derivs.fun, function(pmf.deriv.fun) {
##     units.mat <- sapply(GH$theta,function(thetai)pmf.deriv.fun(beta,params,thetai,supp.x))
##     ## g.mat <- outer(X=GH$theta.grid,Y=supp.x,function(theta,x)pmf.deriv.fun(beta=beta,theta=theta,x=x,params=params)) 
##     colSums(GH$w * t(units.mat))
## })
## k.vals <- sapply(seq_along(ks.fun), function(i) {
##     k.fun <- ks.fun[[i]]
##     args <- attr(k.fun,'args')
##     k.fun.unit <- do.call(k.fun, pmf.derivs.unit[args])
##     k.val <- sum( k.fun.unit )
## })
## k.orders <- sapply(ks.fun,function(k.fun)attr(k.fun,'order'))
## ## k.df <- data.frame(inner=sapply(k.orders,`[`,1),outer=sapply(k.orders,`[`,2),k.val=k.vals)
## ks.df <- data.frame(inner=k.orders[1,],outer=k.orders[2,],k.val=k.vals)
## k.22 <- ks.df[ks.df[,'inner']=='22'  & ks.df[,'outer']=='0','k.val']
## k.11 <- ks.df[ks.df[,'inner']=='11'  & ks.df[,'outer']=='0','k.val']
## k.12 <- ks.df[ks.df[,'inner']=='12'  & ks.df[,'outer']=='0','k.val']
## print(k.22); print(k.11)
## bartlett.C <- lawley(ks.df,2) - lawley(ks.df,1)
## dput(bartlett.C)

## set.seed(1)
rm(list=ls())
source('utils.R')
ks.fun <- readRDS('ks-fun.rds')
pmf.derivs.fun <- readRDS('pmf-derivs-fun-hyperg.rds')
n <- 1e1
mu.null <- 0; tau2.null <- 3
paramsi <- c(m.T=20,m.C=15,s=10)
paramsi <- c(paramsi,x.min=max(0,paramsi['s']-paramsi['m.C']),x.max=min(paramsi['s'],paramsi['m.T']))
## paramsi <- c(m.T = 174, m.C = 177, s = 4, x.min = 0, x.max = 4)
n.GH <- 50
## GH <- fastGHQuad::gaussHermiteData(n.GH)
## GH$x <- sqrt(2)*GH$x
## GH$w <- GH$w / sqrt(pi)
## names(GH) <- c('theta.grid','w')
## pmf.deriv.fun <- pmf.derivs.fun[[5]]
ks.study.fit <- ks.study(pmf.derivs.fun,ks.fun,paramsi,mu.null,tau2.null,n.GH=n.GH,debug=TRUE)
ks.df <- ks.study.fit$ks.df
k.22 <- ks.df[ks.df[,'inner']=='22'  & ks.df[,'outer']=='0','k.val']
k.11 <- ks.df[ks.df[,'inner']=='11'  & ks.df[,'outer']=='0','k.val']
k.12 <- ks.df[ks.df[,'inner']=='12'  & ks.df[,'outer']=='0','k.val']
print(k.22); print(k.11)
bartlett.C <- lawley(ks.df,2) - lawley(ks.df,1)
dput(bartlett.C)
inf <- -matrix(c(k.11,k.12,k.12,k.22),nrow=2)
vcov.hat <- solve(inf)
## z.stats <- replicate(1e2, {
##     params <- cbind(m.T=rep(m.T,n),m.C=rep(m.C,n),s=rep(s,n))
##     x <- dgp.hyperg(n,mu.null,tau2.null,params)
##     mles <- LR.gaussian(x,params,mu.null,cond.pmf.hyperg,n.GH=50)$mle.full
##     (rev(mles)-c(tau2.null,mu.null)) / sqrt(diag(vcov.hat)/n)
## })
mles <- replicate(1e2, {
    params <- matrix(paramsi,nrow=n,ncol=length(paramsi),byrow=TRUE)
    colnames(params) <- names(paramsi)
    x <- dgp.hyperg(n,mu.null,tau2.null,params)
    LR.gaussian(x,params,mu.null,cond.pmf.hyperg,n.GH=n.GH)$mle.full |> rev()
})
z.stats <- (mles - c(tau2.null,mu.null)) / sqrt(diag(vcov.hat)/n)
op <- par(mfrow=c(1,2))
qqnorm(z.stats['mu',]); abline(0,1)
qqnorm(z.stats['tau2',]); abline(0,1)
par(op)


##  9. using saved formulas, run through all models

## binomial-prop
source('utils.R')
ks.fun <- readRDS('ks-fun.rds')
mu.null <- 0
tau2.null <- 1/2

pmf.derivs.fun <- readRDS('pmf-derivs-fun-binomial-prop.rds')
m <- 20
supp.x <- 0:m
params <- c(m=m)

ks.df <- ks.study(pmf.derivs.fun=pmf.derivs.fun,ks.fun=ks.fun,params=params,supp.x=supp.x,mu.null=mu.null,tau2.null=tau2.null)
bartlett.C <- lawley(ks.df,2) - lawley(ks.df,1)
dput(bartlett.C)
## 1.32617146423958


## hypergeometric
source('utils.R')
ks.fun <- readRDS('ks-fun.rds')
mu.null <- 0
tau2.null <- 1/2

pmf.derivs.fun <- readRDS('pmf-derivs-fun-hyperg.rds')
m <- m.T <- 30; m.C <- 30; s <- 15;
supp.x <- max(0,s-m.C):min(s,m.T)
params <- c(m.T=m.T,m.C=m.C,s=s)

ks.df <- ks.study(pmf.derivs.fun=pmf.derivs.fun,ks.fun=ks.fun,params=params,supp.x=supp.x,mu.null=mu.null,tau2.null=tau2.null,n.GH=1e2)
bartlett.C <- lawley(ks.df,2) - lawley(ks.df,1)
dput(bartlett.C)
## 1.37680959521052




## poisson
source('utils.R')
mu.null <- -4
tau2.null <- .2
ks.fun <- readRDS('ks-fun.rds')

pmf.derivs.fun <- readRDS('pmf-derivs-fun-poisson.rds')
T <- 1e3
supp.x <- 0:1e3
params <- c(T=T)

ks.df <- ks.study(pmf.derivs.fun=pmf.derivs.fun,ks.fun=ks.fun,params=params,supp.x=supp.x,mu.null=mu.null,tau2.null=tau2.null)
bartlett.C <- lawley(ks.df,2) - lawley(ks.df,1)
dput(bartlett.C)
## 1.47249385400558



## binomial-IRR 
source('utils.R')
ks.fun <- readRDS('ks-fun.rds')
mu.null <- 0
tau2.null <- 1/2

pmf.derivs.fun <- readRDS('pmf-derivs-fun-binomial-irr.rds')
s <- 30
T.t <- 1
T.c <- 1.5
params <- c(s=s,T.t=T.t,T.c=T.c)
supp.x <- 0:s

ks.df <- ks.study(pmf.derivs.fun=pmf.derivs.fun,ks.fun=ks.fun,params=params,supp.x=supp.x,mu.null=mu.null,tau2.null=tau2.null)
bartlett.C <- lawley(ks.df,2) - lawley(ks.df,1)
print(bartlett.C)
## 1.36297164231076






## 10 data analysis, INID data

source('utils.R')
mu.null <- -3.5
tau2.null <- .81
ks.fun <- readRDS('ks-fun.rds')
pmf.derivs.fun <- readRDS('pmf-derivs-fun-binomial-prop.rds')
niel.weise <- read.table('niel-weise.txt')

## set.seed(2)
## n <- 10
## ms <- sample(niel.weise$n0,n,replace=TRUE)
## ## ms <- rep(round(median(ms)),n)
## ks <- lapply(ms, function(m) {
##     supp.x <- 0:m
##     params <- c(m=m)
##     ks.df.study <- ks.study(pmf.derivs.fun=pmf.derivs.fun,ks.fun=ks.fun,params=params,supp.x=supp.x,mu.null=mu.null,tau2.null=tau2.null,n.GH=1e3)    
## })
## k.vals.ma <- sapply(ks,function(df)df$k.val) |> rowSums()
## ks.df.ma <- data.frame(inner=ks[[1]]$inner, outer=ks[[1]]$outer,k.val=k.vals.ma)
## bartlett.C <- lawley(ks.df.ma,2) - lawley(ks.df.ma,1)
## dput(bartlett.C)





## 10a binomial-IRR
source('utils.R')
ks.fun <- readRDS('ks-fun.rds')
pmf.derivs.fun <- readRDS('pmf-derivs-fun-binomial-irr.rds')
stijnen.6 <- read.table('stijnen-table-VI.txt')
x <- x.T <- stijnen.6$Y1
x.C <- stijnen.6$Y0
T.t <- stijnen.6$T1/1000
T.c <- stijnen.6$T0/1000
s <- x.C+x.T
params <- cbind(T.t=T.t,T.c=T.c,s=s,x.min=0,x.max=s)
LR.fit <- LR.gaussian(x=x,params=params,mu.null=0,cond.pmf=cond.pmf.binomial.irr,logtau2.bounds=c(-10,10),par.start=par.start.binomial.irr(x,params),n.GH=50)
## print(fit$mle.full) #matches stijnen and metafor
## require(metafor) 
## fit <- rma.glmm(measure="IRLN", xi=Y0, ti=T0/1000, data=stijnen.vi, method="ML")
## summary(fit)
start <- Sys.time()
profile.fit <- profile.LR(LR.fit,grid=list(offset=3.5,size=20),logtau2.bounds=c(-5,5))
print(Sys.time() - start)
CIs <- confint.LR(profile.fit,alpha=.05)
print(CIs)
## > > $LR
##        CI.left     CI.right
## [1,] -1.053636 -0.005980577

## $LR.bartlett
##        CI.left   CI.right
## [1,] -1.119055 0.04583006



## diagnostic for n.GH
n.GHs <- seq(50,250,by=50)
cond.pmf <- cond.pmf.binomial.irr
## cond.pmf <- cond.pmf.poisson
by.n.GH <- lapply(n.GHs,function(n.GH) {
    LR.fit <- LR.gaussian(x=x,params=params,mu.null=mu.null,cond.pmf=cond.pmf,par.start=par.start.binomial.irr(x,params),logtau2.bounds=c(-10,10),n.GH=n.GH)
    profile.fit <- profile.LR(LR.fit,grid=list(offset=1.5,size=10,mc.cores=5))
    profile.fit
})

stats <- lapply(by.n.GH,FUN=function(m)m$profile['LR.bartlett',])    
stats <- simplify2array(stats)
col <- cm.colors(ncol(stats))
matplot(by.n.GH[[1]]$mu.grid,stats,type='l',lty=1,col=col)
legend('topleft',lty=1,col=col,legend=n.GHs)
abline(h=qchisq(1-.05,df=1))


## 10b binomial-prop.  !!ugly curve for large n.GH [update: needed to
## increase n.GH. for small n.GH the curve looks fine but the
## information is unstable. once k terms stabilize the bartlett curve and CIs are essentially the same as the uncorrected LR's].
source('utils.R')
ks.fun <- readRDS('ks-fun.rds')
stijnen.1 <- read.table('stijnen-table-I.txt')
stijnen.1 <- stijnen.1[-15,] # get rid of uninformative study
## niel.weise <- read.table('niel-weise.txt')
pmf.derivs.fun <- readRDS('pmf-derivs-fun-binomial-prop.rds')
x <- stijnen.1$Y1
ms <- stijnen.1$n1
params <- cbind(m=ms,x.min=0,x.max=ms)
LR.fit <- LR.gaussian(x=x,params=params,mu.null=0,cond.pmf=cond.pmf.binomial.prop,logtau2.bounds=c(-5,5),par.start=par.start.binomial.prop(x,params),n.GH=400)
start <- Sys.time()
profile.fit <- profile.LR(LR.fit,grid=list(offset=1.2,size=20),logtau2.bounds=c(-10,10),mc.cores=5)
print(Sys.time() - start)
CIs <- confint.LR(profile.fit,alpha=.05)
print(CIs)
## control group
## $LR
##        CI.left  CI.right
## [1,] -3.937548 -2.916536
## $LR.bartlett
##       CI.left CI.right
## [1,] -3.96277 -2.89521

## treatment group
## > $LR
##        CI.left  CI.right
## [1,] -5.572486 -4.132205
## $LR.bartlett
##        CI.left  CI.right
## [1,] -5.606172 -4.109588


## diagnostic for n.GH
n.GHs <- seq(100,400,by=100)
cond.pmf <- cond.pmf.binomial.prop
## cond.pmf <- cond.pmf.poisson
by.n.GH <- lapply(n.GHs,function(n.GH) {
    LR.fit <- LR.gaussian(x=x,params=params,mu.null=mu.null,cond.pmf=cond.pmf,par.start=par.start.binomial.prop(x,params),logtau2.bounds=c(-10,10),n.GH=n.GH)
    profile.fit <- profile.LR(LR.fit,grid=list(offset=1.5,size=10,mc.cores=5))
    profile.fit
})

stats <- lapply(by.n.GH,FUN=function(m)m$profile['LR.bartlett',])    
stats <- simplify2array(stats)
col <- cm.colors(ncol(stats))
matplot(by.n.GH[[1]]$mu.grid,stats,type='l',lty=1,col=col)
legend('topleft',lty=1,col=col,legend=n.GHs)
abline(h=qchisq(1-.05,df=1))
## abline(v=CIs[,'LR.bartlett'])

lapply(by.n.GH,confint.LR)

## patwardhan data
source('utils.R')
ks.fun <- readRDS('ks-fun.rds')
stijnen.v <- read.table('stijnen-table-V.txt')
pmf.derivs.fun <- readRDS('pmf-derivs-fun-binomial-prop.rds')
x <- stijnen.v$TP
ms <- stijnen.v$n1
params <- cbind(m=ms,x.min=0,x.max=ms)
LR.fit <- LR.gaussian(x=x,params=params,mu.null=0,cond.pmf=cond.pmf.binomial.prop,logtau2.bounds=c(-5,5),n.GH=1e2)
start <- Sys.time()
profile.fit <- profile.LR(LR.fit,grid=list(offset=1.5,size=20),logtau2.bounds=c(-5,5))
print(Sys.time() - start)
CIs <- confint.LR(profile.fit,alpha=.05)
print(CIs)
plogis(CIs)



## 10c poisson data analysis
source('utils.R')
ks.fun <- readRDS('ks-fun.rds')
pmf.derivs.fun <- readRDS('pmf-derivs-fun-poisson.rds')
stijnen.vi <- read.table('stijnen-table-VI.txt')
n.GH <- 3e2
x <- stijnen.vi$Y1
T <- stijnen.vi$T1/1000
mu.null <- 0
params <- cbind(T=T,x.min=0,x.max=2.5e3)
LR.fit <- LR.gaussian(x=x,params=params,mu.null=mu.null,cond.pmf=cond.pmf.poisson,logtau2.bounds=c(-5,5),n.GH=n.GH)
mu.null <- LR.fit$mle.full['mu']; tau2.null <- LR.fit$mle.full['tau2']
## print(fit$mle.full) #matches stijnen and metafor
## require(metafor) 
## fit <- rma.glmm(measure="IRLN", xi=Y0, ti=T0/1000, data=stijnen.vi, method="ML")
## summary(fit)
start <- Sys.time()
profile.fit <- profile.LR(LR.fit,grid=list(offset=3.5,size=20),logtau2.bounds=c(-5,5),mc.cores=5)
print(Sys.time() - start)
CIs <- confint.LR(profile.fit,alpha=.05)
print(CIs)


## poisson diagnostic
source('utils.R')
ks.fun <- readRDS('ks-fun.rds')
pmf.derivs.fun <- readRDS('pmf-derivs-fun-poisson.rds')
stijnen.vi <- read.table('stijnen-table-VI.txt')
x <- stijnen.vi$Y0
T <- stijnen.vi$T0/1000
mu.null <- 0
n.GH <- seq(100,800,by=100)
x.max <- seq(1e3,1.5e4,by=1e3)
pars <- expand.grid(n.GH=n.GH,x.max=x.max)
profile.fits <- apply(pars, 1, function(par) {
    n.GH <- par['n.GH']; x.max <- par['x.max']
    # browser()
    params <- cbind(T=T,x.min=0,x.max=x.max)
    par.start <- par.start.poisson(x,params)
    LR.fit <- LR.gaussian(x=x,params=params,mu.null=mu.null,cond.pmf=cond.pmf.poisson,logtau2.bounds=c(-5,5),par.start=par.start,n.GH=n.GH)
    mu.null <- LR.fit$mle.full['mu']; tau2.null <- LR.fit$mle.full['tau2']
    start <- Sys.time()
    profile.fit <- tryCatch({
        fit <- profile.LR(LR.fit,grid=list(offset=3.5,size=10),logtau2.bounds=c(-5,5),mc.cores=parallel::detectCores()-2)
        print(confint.LR(fit,plot=FALSE))
        fit},
        error=function(e)NA)
    print(Sys.time() - start)
    profile.fit
})

stats <- lapply(profile.fits,FUN=function(m)m$profile['LR.bartlett',]) 
stats <- simplify2array(stats)
col <- cm.colors(ncol(stats))
matplot(profile.fits[[1]]$mu.grid,stats,type='l',lty=1,col=col)
legend('topleft',lty=1,col=col,legend=pars$x.max)
abline(h=(qchisq(1-.05,df=1)))
lapply(profile.fits,confint.LR)



# n.GH <- seq(100,500,by=100)
# x.max <- seq(1e3,5e3,by=1e3)
# pars <- expand.grid(n.GH=n.GH,x.max=x.max)
load('../poisson diagnostic.RData')
source('utils.R')
CIs <- lapply(profile.fits,function(profile.fit) {
    tryCatch(confint.LR(profile.fit,plot=FALSE)$LR.bartlett,error=function(e)NA)
})
## unlist(sapply(CIs, function(mat)dim(mat)[1]) ) 
## bad.idx <- unlist(sapply(CIs, function(mat)dim(mat)[1]) )!=1
bad.idx <- sapply(CIs, function(mat)dim(mat)[1])
bad.idx[sapply(bad.idx,is.null)] <- NA
bad.idx <- unlist(bad.idx) !=1
CIs[bad.idx] <- NA
CIs.str <- sapply(CIs, function(CI) paste0(round(CI,4),collapse=', '))
CI.mat <- matrix(CIs.str,ncol=5,byrow=FALSE)
rownames(CI.mat) <- unique(pars$n.GH)
colnames(CI.mat) <- unique(pars$x.max)


idx <- which(pars$n.GH==700)
idx <- which(pars$x.max>=6e3)
stats <- lapply(profile.fits[idx],FUN=function(m)m$profile['LR.bartlett',])    
stats <- simplify2array(stats)
col <- cm.colors(ncol(stats))
matplot(profile.fits[[1]]$mu.grid,stats,type='l',lty=1,col=col)
legend('topleft',lty=1,col=col,legend=pars$x.max[idx])
abline(h=(qchisq(1-.05,df=1)))




## 10d hyperg. n.GH>55 or so gets NaNs. biasedUrn hyperg pmf returns them
## when the odds is very large. oddly this is only when passing a
## vector of x's, otherwise, biasedUrn just returns the expected 0
## therm (example below). but non-vectorized it is too slow.
source('utils.R')
ks.fun <- readRDS('ks-fun.rds')
pmf.derivs.fun <- readRDS('pmf-derivs-fun-hyperg.rds')
stijnen.1 <- read.table('stijnen-table-I.txt')
stijnen.1 <- stijnen.1[-15,] # get rid of uninformative study
x <- stijnen.1$Y1
s <- with(stijnen.1, Y0+Y1)
m.T <- stijnen.1$n1
m.C <- stijnen.1$n0
n.GH <- 50
## supp.x <- max(0,s-m.C):min(s,m.T)
params <- cbind(m.T=m.T,m.C=m.C,s=s,x.min=pmax(0,s-m.C),x.max=pmin(s,m.T))
LR.fit <- LR.gaussian(x=x,params=params,mu.null=0,cond.pmf=cond.pmf.hyperg,par.start=par.start.hyperg(x,params),n.GH=n.GH,control = list(maxit = 5e3))
start <- Sys.time()
profile.fit <- profile.LR(LR.fit,grid=list(offset=1.5,size=20),logtau2.bounds=c(-10,10),mc.cores=2)
print(Sys.time() - start)
CIs <- confint.LR(profile.fit)
print(CIs)
## $LR
##        CI.left   CI.right
## [1,] -2.216696 -0.7077608

## $LR.bartlett
##        CI.left   CI.right
## [1,] -2.253979 -0.6763749


## biasedUrn stange behavior.

beta <- c(6.46818227043636, 0.728696470009073)
paramsi <- c(m.T = 151, m.C = 157, s = 11, x.min = 0, x.max = 11)
thetai <- 14.3671500775778
m.T <- paramsi['m.T']
m.C <- paramsi['m.C']
s <- paramsi['s']
odds <- exp(1)^(thetai * sqrt(beta[1]) + beta[2])
BiasedUrn::dFNCHypergeo(5,m1=m.T,m2=m.C,n=s,odds=odds)        



## diagnostic for n.GH
n.GHs <- c(20,30,40,50)
cond.pmf <- cond.pmf.hyperg
## cond.pmf <- cond.pmf.poisson
by.n.GH <- lapply(n.GHs,function(n.GH) {
    LR.fit <- LR.gaussian(x=x,params=params,mu.null=0,cond.pmf=cond.pmf,par.start=par.start.hyperg(x,params),logtau2.bounds=c(-10,10),n.GH=n.GH)
    profile.fit <- profile.LR(LR.fit,grid=list(offset=1.5,size=10,mc.cores=5))
    profile.fit
})

stats <- lapply(by.n.GH,FUN=function(m)m$profile['LR.bartlett',])    
stats <- simplify2array(stats)
col <- cm.colors(ncol(stats))
matplot(by.n.GH[[1]]$mu.grid,stats,type='l',lty=1,col=col)
legend('topleft',lty=1,col=col,legend=n.GHs)
abline(h=qchisq(1-.05,df=1))
## abline(v=CIs[,'LR.bartlett'])






## refactoring code
devtools::load_all('./meta.glmm')
## ## OLD WAY:
## source('utils.R')
## ks.fun <- readRDS('ks-fun.rds')
## LR.fit <- LR.gaussian(x=x,params=params,mu.null=0,cond.pmf=cond.pmf.binomial.prop,logtau2.bounds=c(-10,10),par.start=par.start.binomial.prop(x,params),n.GH=400)
## start <- Sys.time()
## profile.fit <- profile.LR(LR.fit,grid=list(offset=1,size=10),logtau2.bounds=c(-10,10),mc.cores=5)
## print(Sys.time() - start)
## CIs <- confint.LR(profile.fit,alpha=.05)
## source('utils.R')
stijnen.1 <- read.table('stijnen-table-I.txt')
stijnen.1 <- stijnen.1[-15,] # get rid of uninformative study
## pmf.derivs.fun <- readRDS('pmf-derivs-fun-binomial-prop.rds')
x <- stijnen.1$Y1
ms <- stijnen.1$n1
params <- cbind(m=ms,x.min=0,x.max=ms)
n.GH <- 50
meta.glmm.fit <- meta.glmm(x=x,params=params,n.GH=n.GH,model.str='bn-prop',control = list(maxit = 5e3))
CI <- confint.meta.glmm(meta.glmm.fit,grid=list(offset=1,size=10),mc.cores=5,bartlett=TRUE)



## 11. adaptive GH

source('utils.R')
ks.fun <- readRDS('ks-fun.rds')
pmf.derivs.fun <- readRDS('pmf-derivs-fun-poisson.rds')
stijnen.vi <- read.table('stijnen-table-VI.txt')
n.GH <- 3e2
x <- stijnen.vi$Y1
T <- stijnen.vi$T1/1000
mu.null <- 0
params <- cbind(T=T,x.min=0,x.max=2.5e3)
LR.fit <- LR.gaussian(x=x,params=params,mu.null=mu.null,cond.pmf=cond.pmf.poisson,logtau2.bounds=c(-5,5),n.GH=n.GH)
mu.null <- LR.fit$mle.full['mu']; tau2.null <- LR.fit$mle.full['tau2']

GH.LR <- function(xi, paramsi, mu, tau2, n.GH=20, tol=1e-10, adaptive=TRUE) {
    ## browser()
    rule <- fastGHQuad::gaussHermiteData(n.GH)
    if(!adaptive)return(list(wts=rule$w/sqrt(pi), nodes=sqrt(2)*rule$x))
    integrand.mode <- function(xi, paramsi, mu, tau2) {
        ## browser()
        z <- 0
        for (j in 1:20) {
            lambda <- paramsi[,"T"] *
                exp(mu + sqrt(tau2) * z)
            f <- sqrt(tau2) * (xi - lambda) - z
            fp <- -1 - tau2 * lambda
            z.new <- z - f / fp
            if (abs(z.new - z) < tol) {
                break
            }
            z <- z.new
        }
        mode <- z.new
        return(mode)
    }
    var.map <- function(zi,paramsi,mu,tau2)paramsi[,'T']*exp(mu+sqrt(tau2)*zi)
    mode <- integrand.mode(xi,paramsi,mu,tau2)
    scale <- (1 + tau2 * var.map(mode,paramsi, mu, tau2))^(-1 / 2)
    theta <- mode + sqrt(2) * scale * rule$x
    w <- sqrt(2) * scale *
        rule$w *
        exp(rule$x^2) *
        dnorm(theta)
    return(list(wts=w,nodes=theta))
}
LR.gaussian <- function(x,params,mu.null,cond.pmf,logtau2.bounds=c(-3,3),par.start=c(mu=0,tau2=1),n.GH=1e2,...) {
    ## browser()
    GH <- fastGHQuad::gaussHermiteData(n.GH)
    GH$x <- sqrt(2)*GH$x; GH$w <- GH$w/sqrt(pi)
    loglik.unit <- function(mu,tau2,xi,paramsi,GH) {
        ## browser()
        ## print(c(mu,tau2,xi,paramsi))
        ## if(tau2>1e3)browser()
        GH <- GH.LR(xi=xi,paramsi=paramsi,mu=mu,tau2=tau2,n.GH=n.GH,adaptive=TRUE)
        if(log(tau2)<logtau2.bounds[1] || log(tau2)>logtau2.bounds[2]) {
            ## warning('tau2 boundary  hit, returning minimum negative')
            return(-1e10)#     return(-.Machine$double.xmax)
        }
        likelihood <- sum(GH$wts*cond.pmf(xi,GH$nodes,mu,tau2,paramsi))
        ## likelihood <- sum(GH$w*cond.pmf(xi,GH$x,mu,tau2,paramsi))
        ## print(c(mu,tau2,log(likelihood),paramsi))
        ## print(tau2)
        return(log(likelihood))
    }
    loglik <- function(mu,tau2,x,params,GH) {
        out <- sum(sapply(seq_along(x),function(i)loglik.unit(mu,tau2,x[i],params[i,,drop=FALSE],GH)))
        ## print(out)        
        ## if(is.nan(out))browser()
        if(is.nan(out))out <- -1e10#-.Machine$double.xmax
        return(out)
        }
    logtau2.mle.redu <- optimize(function(logtau2)loglik(mu.null,exp(logtau2),x,params,GH),interval=logtau2.bounds,maximum=TRUE)$max
    tau2.mle.redu <- exp(logtau2.mle.redu)
    start.val <- loglik(par.start[1],par.start[2],x,params,GH)
    if(is.nan(start.val) || !is.finite(start.val))browser()
    fit.full <- optim(par=c(par.start['mu'],log(par.start['tau2'])),fn=function(par)-loglik(par[1],exp(par[2]),x,params,GH),method='BFGS',...)
    if(fit.full$convergence != 0L)warning('full fit convergence failure')
    mle.full <- fit.full$par
    mu.mle.full <- mle.full[1]; tau2.mle.full <- exp(mle.full[2])
    LR <- 2 * (loglik(mu.mle.full,tau2.mle.full,x,params,GH) - loglik(mu.null,tau2.mle.redu,x,params,GH))
    return(list(mle.full=c(mu=unname(mu.mle.full),tau2=unname(tau2.mle.full)),mle.redu=c(mu=unname(mu.null),tau2=unname(tau2.mle.redu)), LR=LR, loglik.fun=loglik, x=x, params=params, cond.pmf=cond.pmf,logtau2.bounds=logtau2.bounds, n.GH=n.GH))
}
ks.fun <- readRDS('ks-fun.rds')
pmf.derivs.fun <- readRDS('pmf-derivs-fun-poisson.rds')
stijnen.vi <- read.table('stijnen-table-VI.txt')
n.GH <- 30
x <- stijnen.vi$Y1
T <- stijnen.vi$T1/1000
mu.null <- 0
params <- cbind(T=T,x.min=0,x.max=2.5e3)
LR.fit <- LR.gaussian(x=x,params=params,mu.null=mu.null,cond.pmf=cond.pmf.poisson,logtau2.bounds=c(-5,5),n.GH=n.GH)
LR.fit$LR



devtools::load_all('./meta.glmm')
## source('utils.R')
## ks.fun <- readRDS('ks-fun.rds')
## pmf.derivs.fun <- readRDS('pmf-derivs-fun-poisson.rds')
stijnen.vi <- read.table('stijnen-table-VI.txt')
n.GH <- 1e2
x <- stijnen.vi$Y1
T <- stijnen.vi$T1/1000
mu.null <- 0
params <- cbind(T=T,x.min=0,x.max=5e2)
meta.glmm.fit <- meta.glmm(x=x,params=params,n.GH=n.GH,model.str='pn',control = list(maxit = 5e3),adaptive=FALSE)

CI <- confint.meta.glmm(meta.glmm.fit,grid=list(offset=3,size=10),mc.cores=5,bartlett=TRUE)








## 12, second re-factoring

source('models.R')
source('./meta.glmm/R/meta.glmm.R')
stijnen.vi <- read.table('stijnen-table-VI.txt')
n.GH <- 30
y <- stijnen.vi$Y1
T <- stijnen.vi$T1/1000
params <- cbind(T=T,y.min=0,y.max=1e2)
mu.null <- 0
tau2.null <- 1
mod <- model('pn')
quad <- make.quad(integrand.mode,mod$var.map,n.GH=n.GH,adaptive=TRUE)
## quad(x[1],params[1,],mu.null,tau2.null)
loglik <- make.loglik(mod$cond.pmf,params,y,quad)
## loglik(mu.null,tau2.null)
nn.data <- mod$nn.transform(y,params)
## par.start <- mod$par.start(y,params)
mle0 <- mle(loglik,par.start)$mle
profile.LR <- make.profile.LR(loglik,mle0,par.start,logtau2.bounds=c(-20,3))
profile.LR(0)
grid <- list(offset=2,size=10)
mu.grid <- seq(mle0['mu']-grid$offset*sqrt(mle0['tau2']),mle0['mu']+grid$offset*sqrt(mle0['tau2']),len=grid$size)
## profile <- profile(mu.grid,profile.LR, mle,bartlett=TRUE)
bartlett.C <- make.bartlett.C(params,quad,mod$pmf.derivs.fun)
## bartlett.C(mu.null,tau2.null)
profiled <- profile(mu.grid,profile.LR,mle0,bartlett.C,mc.cores=5)
profiled <- rbind(LR=profiled$LR,LR.bartlett=profiled$LR.bartlett)
confint(mu.grid=mu.grid,profile=profiled)
## using constructor/convenience functions:

## 12a redoing data analysis
    
## 12a-1 binomial-prop:
ks.fun <- readRDS('ks-fun.rds')
pmf.derivs.fun.pn <- readRDS('pmf-derivs-fun-poisson.rds')
pmf.derivs.fun.bn.prop <- readRDS('pmf-derivs-fun-binomial-prop.rds')
pmf.derivs.fun.bn.irr <- readRDS('pmf-derivs-fun-binomial-irr.rds')
pmf.derivs.fun.hn <- readRDS('pmf-derivs-fun-hyperg.rds')
source('models.R')
source('./meta.glmm/R/meta.glmm.R')

stijnen.1 <- read.table('stijnen-table-I.txt')
y <- stijnen.1$Y1
ms <- stijnen.1$n1
params <- cbind(m=ms,y.min=0,y.max=ms)
n.GH <- 350
start <- Sys.time()
fit <- meta.glmm(y,params,'bn-prop',n.GH=n.GH,adaptive=FALSE,logtau2.bounds=c(-20,3))#,control=list(trace=4,ndeps = c(1e-7, 1e-7)))
CIs <- confint.meta.glmm(fit,mc.cores=5,bartlett=TRUE,plot=TRUE,mu.grid=list(offset=1,size=20),col=1)
print(CIs)
print(Sys.time() - start)


## varying n.GH
stijnen.1 <- read.table('stijnen-table-I.txt')
## y <- stijnen.1$Y0
## ms <- stijnen.1$n0
y <- stijnen.1$Y1
ms <- stijnen.1$n1
params <- cbind(m=ms,y.min=0,y.max=ms)
n.GHs <- seq(1e2,6e2,by=100)
## n.GHs <- seq(10,1e2,by=10)
by.n.GH <- lapply(n.GHs, function(n.GH) {
    print(n.GH)
    start <- Sys.time()
    fit <- meta.glmm(y,params,'bn-prop',n.GH=n.GH,adaptive=TRUE,logtau2.bounds=c(-20,3))#,control=list(trace=4,ndeps = c(1e-7, 1e-7)))
    CIs <- confint.meta.glmm(fit,mc.cores=parallel::detectCores()-2,,plot=FALSE)
    print(CIs)
    print(Sys.time() - start)
    return(CIs)
})
## save.image('bnprop diagnostic ctrl nonadaptive.RData')
## save.image('bnprop diagnostic ctrl adaptive.RData')
lapply(lapply(by.n.GH,`[[`,'CIs'),`[[`,'LR.bartlett')

idx <- n.GHs %in% c(20,40,60,80)
n.GHs[idx]
lapply(lapply(by.n.GH,`[[`,'CIs'),`[[`,'LR.bartlett')[idx]

idx <- n.GHs %in% c(50,150,250,350)
n.GHs[idx]
lapply(lapply(by.n.GH,`[[`,'CIs'),`[[`,'LR.bartlett')[idx]

with(by.n.GH[[8]],plot(mu.grid,profile['LR',],type='l'))
with(by.n.GH[[8]],lines(mu.grid,profile['LR.bartlett',],type='l',lty=2))
alpha <- .05
abline(h=qchisq(1-alpha,df=1))

dd


## converged CI for control group:
##                              [[8]]
##        CI.left  CI.right
## [1,] -4.089781 -2.980016


## 12a-2 hyperg
source('models.R')
source('./meta.glmm/R/meta.glmm.R')
stijnen.1 <- read.table('stijnen-table-I.txt')
stijnen.1 <- stijnen.1[-15,] # get rid of uninformative study
y <- stijnen.1$Y1
s <- with(stijnen.1, Y0+Y1)
m.T <- stijnen.1$n1
m.C <- stijnen.1$n0
n.GH <- 20
## supp.x <- max(0,s-m.C):min(s,m.T)
params <- cbind(m.T=m.T,m.C=m.C,s=s,y.min=pmax(0,s-m.C),y.max=pmin(s,m.T))
fit <- meta.glmm(y,params,'hn',n.GH=n.GH,adaptive=TRUE,logtau2.bounds=c(-20,3))#,control=list(trace=4,ndeps = c(1e-7, 1e-7)))
CIs <- confint.meta.glmm(fit,mc.cores=5,bartlett=TRUE,mu.grid=list(offset=1.5,size=10))



## 12a-3 binomial-IRR
source('models.R')
source('./meta.glmm/R/meta.glmm.R')
stijnen.6 <- read.table('stijnen-table-VI.txt')
y <- x.T <- stijnen.6$Y1
x.C <- stijnen.6$Y0
T.T <- stijnen.6$T1/1000
T.C <- stijnen.6$T0/1000
s <- x.C+x.T
params <- cbind(T.T=T.T,T.C=T.C,s=s,y.min=0,y.max=s)
n.GH <- 15
fit <- meta.glmm(y,params,'bn-irr',n.GH=n.GH,adaptive=TRUE,logtau2.bounds=c(-20,3))#,control=list(trace=4,ndeps = c(1e-7, 1e-7)))
CIs <- confint.meta.glmm(fit,mc.cores=5,bartlett=TRUE,mu.grid=list(offset=2,size=20))
print(CIs)




## 12a-4 poisson
ks.fun <- readRDS('ks-fun.rds')
pmf.derivs.fun.pn <- readRDS('pmf-derivs-fun-poisson.rds')
source('models.R')
source('./meta.glmm/R/meta.glmm.R')
stijnen.vi <- read.table('stijnen-table-VI.txt')
n.GH <- 40
y <- stijnen.vi$Y1
T <- stijnen.vi$T1/1000
params <- cbind(T=T,y.min=0,y.max=850)
fit <- meta.glmm(y,params,'pn',n.GH=n.GH,adaptive=TRUE,logtau2.bounds=c(-20,3))
confint.meta.glmm(fit,mc.cores=5,mu.grid=list(offset=2.,size=20))

## varyng n.GH and y.max
n.GH <- c(seq(10,20,by=5),1e2)
y.max <- c(seq(5e2,1.5e3,by=5e2),5e3)#seq(1e2,3e2,by=50)
n.GH <- c(seq(10,20,by=5),1e2)
y.max <- c(seq(5e2,1.5e3,by=5e2),5e3)#seq(1e2,3e2,by=50)
pars <- expand.grid(n.GH=n.GH,y.max=y.max)
CIs <- apply(pars, 1, function(par) {
    print(par)
    n.GH <- par['n.GH']; y.max <- par['y.max']
    params <- cbind(T=T,y.min=0,y.max=y.max)
    fit <- meta.glmm(y,params,'pn',n.GH=n.GH,adaptive=TRUE,logtau2.bounds=c(-20,3))
    CI <- tryCatch({
        fit <- confint.meta.glmm(fit,mc.cores=parallel::detectCores()-1,mu.grid=list(offset=1.5,size=10),plot=FALSE)
        print(fit)
        fit
    },    error=function(e)NA)
    return(CI)
})

## save.image('pn diagnostic--adaptive--trt 4.RData')
## save.image('pn diagnostic--adaptive--trt 3.RData')
## save.image('pn diagnostic--adaptive--trt 2.RData')
## save.image('pn diagnostic--adaptive--trt.RData')
## save.image('pn diagnostic--adaptive.RData').


lapply(CIs[!is.na(CIs)],`[[`,'CIs') |> sapply(`[[`,'LR.bartlett')  |> t() |> round(4) |> apply(1,paste,collapse=',') |> matrix(nrow=length(unique(pars$n.GH)))

lapply(CIs[!is.na(CIs)],`[[`,'CIs') |> sapply(`[[`,'LR.bartlett')  |> sapply(function(elt)elt[1,]) |> t() |> round(4) |> apply(1,paste,collapse=',') |> matrix(nrow=length(unique(pars$n.GH)))



## 12b R package creation

usethis::create_package("meta.glmm", open = FALSE)
usethis::proj_set("meta.glmm")

setwd('meta.glmm')
usethis::proj_set(".")
ks.fun <- readRDS('../ks-fun.rds')
## pmf.derivs.fun.binomial.prop <- readRDS('../pmf-derivs-fun-binomial-prop.rds')
## pmf.derivs.fun.hyperg <- readRDS('../pmf-derivs-fun-hyperg.rds')
## pmf.derivs.fun.poisson <- readRDS('../pmf-derivs-fun-poisson.rds')
## pmf.derivs.fun.binomial.irr <- readRDS('../pmf-derivs-fun-binomial-irr.rds')
pmf.derivs.fun.bn.prop <- readRDS('../pmf-derivs-fun-binomial-prop.rds')
pmf.derivs.fun.hn <- readRDS('../pmf-derivs-fun-hyperg.rds')
pmf.derivs.fun.pn <- readRDS('../pmf-derivs-fun-poisson.rds')
pmf.derivs.fun.bn.irr <- readRDS('../pmf-derivs-fun-binomial-irr.rds')
## usethis::use_data(ks.fun, pmf.derivs.fun.binomial.prop, pmf.derivs.fun.hyperg, pmf.derivs.fun.poisson, pmf.derivs.fun.binomial.irr, internal=TRUE,base_path=".")
save(
  ks.fun, 
  pmf.derivs.fun.bn.prop, 
  pmf.derivs.fun.hn, 
  pmf.derivs.fun.pn, 
  pmf.derivs.fun.bn.irr, 
  file = "R/sysdata.rda",
  compress = "xz"
)


dir.create("meta.glmm/data", showWarnings = FALSE)

nw08 <- stijnen.6 <- read.table('stijnen-table-VI.txt')
save(nw08, file='meta.glmm/data/nw008.rda', compress='xz')
nw07 <- stijnen.1 <- read.table('stijnen-table-I.txt')
save(nw07, file='meta.glmm/data/nw007.rda', compress='xz')


devtools::load_all('./meta.glmm')




## 12c sim using R package

devtools::load_all('meta.glmm')
source('utils.R') ## for sinh-asinh reparameterization


sas.params <- sas.invert.params(skew,target.kurt=0)
re.distr <- function(n)rsas(n,eps=sas.params['eps'],delta=sas.params['delta'])
## RNGkind("L'Ecuyer-CMRG")#!!!!!!
## set.seed(1)#!!!!!!!!!
## re.distr <- rnorm #!!!!
## tau2.null <- 1/2


ns <- 10:20
mu.null <- -1.353
tau2.null <- .83^2
model.str <- 'hn'
n.GH <- 2e1
## stijnen.1 <- read.table('../../stijnen-table-I.txt')
## stijnen.1 <- stijnen.1[-15,] # get rid of uninformative study
stijnen.params <- with(nw07, cbind(m.T=n1,m.C=n0,s=Y0+Y1))
stijnen.params <- with(as.data.frame(stijnen.params), cbind(stijnen.params,y.min=pmax(0,s-m.C),y.max=pmin(s,m.T)) )                  
mc.idx <- sample(nrow(stijnen.params),max(ns),replace=TRUE)
mc.params <- stijnen.params[mc.idx,]
model0 <- model(model.str)
quad <- make.quad(model0$moment.map,model0$cond.pmf,n.GH=n.GH,adaptive=TRUE,adaptive.tol=1e-5)
start <- Sys.time()
## by.n <- parallel::mclapply(seq_along(ns),mc.cores=6,FUN= function(i) {
by.n <- sapply(seq_along(ns), function(i) {
    ## repeat{
    cat('.')
    n <- ns[i]
    params <- mc.params[1:n,]
    bartlett.C <- make.bartlett.C(params,quad,model0$pmf.derivs.fun)
    bartlett.C0 <- bartlett.C(mu.null,tau2.null)$bartlett.C
    stats <- replicate(round(n^4), {        
        y <- model0$dgp(n,mu.null,tau2.null,params,re.distr=re.distr)
        fit <- meta.glmm(y,params,model.str,n.GH,adaptive=FALSE,logtau2.bounds=c(-20,3))        
        LR <- fit$profile.LR(mu.null)$LR
        ## loglik <- make.loglik(model0$cond.pmf,params,y,quad)
        ## logtau2.mle.null <- optimize(function(logtau2)loglik(mu.null,exp(logtau2)),interval=c(-20,3),maximum=TRUE)$max
        ## tau2.mle.null <- exp(logtau2.mle.null)
        ## par.start0 <- c(mu=mu.null,tau2=tau2.mle.null)
        ## mle0 <- mle(loglik,par.start0,method='BFGS')$mle
        ## profile.LR <- make.profile.LR(loglik,mle0,par.start0,logtau2.bounds=c(-20,3))
        ## LR <- profile.LR(mu.null)$LR
        LR.bartlett <- LR / (1+bartlett.C0)
        c(LR=LR,LR.bartlett=LR.bartlett)
    })
    ecdfs <- apply(stats,1,function(x)ecdf(pchisq(x,df=1,lower=FALSE)))
    ks <- sapply(ecdfs,function(ecdf)max(abs(ecdf(knots(ecdf)) - knots(ecdf))))
    means <- rowMeans(stats-1)
    c(ks=ks,means=means)
})
by.n <- simplify2array(by.n)
print(Sys.time() - start)
## plot(ns,by.n['ks.LR',],type='l',ylim=range(by.n[1:2,]))
## lines(ns,by.n['ks.LR.bartlett',],col=2)
## apply(by.n[1:2,],1,function(ks)coef(lm(log(ks) ~ log(ns)))[2]) |> print()




## ---investigate why par.start with IV/DL doing worse than reduced
## mle. culprit seems to be using IV/D-L estimates. was preivously
## using (mu.null,tau2.mle.redu), which does a lot better.
source('utils.R')
try.mle <- function(x,params,mu.null,cond.pmf,logtau2.bounds=c(-20,3),n.GH=1e2,...) {
    ## browser()
    GH <- fastGHQuad::gaussHermiteData(n.GH)
    GH$x <- sqrt(2)*GH$x; GH$w <- GH$w/sqrt(pi)
    loglik.unit <- function(mu,tau2,xi,paramsi,GH) {
        ## browser()
        ## print(c(mu,tau2,xi,paramsi))
        ## if(tau2>1e3)browser()
        if(log(tau2)<logtau2.bounds[1] || log(tau2)>logtau2.bounds[2]) {
            ## warning('tau2 boundary  hit, returning minimum negative')
            return(-1e10)#     return(-.Machine$double.xmax)
        }
        likelihood <- sum(GH$w*cond.pmf(xi,GH$x,mu,tau2,paramsi))
        ## print(c(mu,tau2,log(likelihood),paramsi))
        ## print(tau2)
        return(log(likelihood))
    }
    loglik <- function(mu,tau2,x,params,GH) {
        out <- sum(sapply(seq_along(x),function(i)loglik.unit(mu,tau2,x[i],params[i,],GH)))
        ## print(out)        
        ## if(is.nan(out))browser()
        if(is.nan(out))out <- -1e10#-.Machine$double.xmax
        return(out)
        }
    logtau2.mle.redu <- optimize(function(logtau2)loglik(mu.null,exp(logtau2),x,params,GH),interval=logtau2.bounds,maximum=TRUE)$max
    tau2.mle.redu <- exp(logtau2.mle.redu)
    par.start <- c(mu=mu.null,tau2=tau2.mle.redu)
    ## par.start <- par.start.binomial.prop(x,params)
    ## print(par.start)
    start.val <- loglik(par.start[1],par.start[2],x,params,GH)
    if(is.nan(start.val) || !is.finite(start.val))browser()
    fit.full <- optim(par=c(par.start['mu'],log(par.start['tau2'])),fn=function(par)-loglik(par[1],exp(par[2]),x,params,GH),method='BFGS',...)
    if(fit.full$convergence != 0L)warning('full fit convergence failure')
    mle.full <- fit.full$par
    mle.full[2] <- exp(mle.full[2])
    return(list(mle.full=mle.full,tau2.mle.redu=tau2.mle.redu))
}
ks.fun <- readRDS('ks-fun.rds')
pmf.derivs.fun <- pmf.derivs.fun.bn.prop <- readRDS('pmf-derivs-fun-binomial-prop.rds')
source('models.R')
source('./meta.glmm/R/meta.glmm.R')
stijnen.1 <- read.table('stijnen-table-I.txt')
y <- stijnen.1$Y1
ms <- stijnen.1$n1
## ms <- rep(20,n)
mu.null <- -3.5
tau2.null <- .81
n <- 10
params <- cbind(m=ms,y.min=0,y.max=ms)
re.distr <- rnorm
ms <- ms[1:n]
params <- params[1:n,]
n.GH <- 3e2
source('utils.R')
set.seed(10)
pairs <- replicate(1e3, {
    x <- dgp.binomial.prop(n,mu.null,tau2.null,params,re.distr=re.distr)
    mle <- try.mle(x,params,mu.null,cond.pmf.binomial.prop,n.GH=1e2,logtau2.bounds=c(-20,3))$tau2.mle.redu
    IV <- par.start.binomial.prop(x,params)['tau2']
    c(mle,IV)
})
plot(pairs[1,],pairs[2,])
abline(0,1)
plot((pairs[1,]-pairs[2,])/pairs[2,])
op <- par(mfrow=c(1,2))
hist(pairs[1,]); abline(v=tau2.null,col=2); abline(v=mean(pairs[1,]),col=3)
hist(pairs[2,]); abline(v=tau2.null,col=2); abline(v=mean(pairs[2,]),col=3)
par(op)


ks <- lapply(ms, function(m) {
    supp.x <- 0:m
    paramsi <- c(m=m,x.min=0,x.max=m)
    ks.df.study <- ks.study(pmf.derivs.fun=pmf.derivs.fun,ks.fun=ks.fun,params=paramsi,mu=mu.null,tau2=tau2.null,n.GH=1e2)    
})
k.vals.ma <- sapply(ks,function(df)df$ks.df$k.val) |> rowSums()
ks.df <- data.frame(inner=ks[[1]]$ks.df$inner, outer=ks[[1]]$ks.df$outer,k.val=k.vals.ma)
k.22 <- ks.df[ks.df[,'inner']=='22'  & ks.df[,'outer']=='0','k.val']
k.11 <- ks.df[ks.df[,'inner']=='11'  & ks.df[,'outer']=='0','k.val']
k.12 <- ks.df[ks.df[,'inner']=='12'  & ks.df[,'outer']=='0','k.val']
inf <- -matrix(c(k.11,k.12,k.12,k.22),nrow=2)
vcov.hat <- solve(inf)
z.stats <- replicate(1e2, {
    x <- dgp.binomial.prop(n,mu.null,tau2.null,params,re.distr=re.distr)
    mles <- try.mle(x,params,mu.null,cond.pmf.binomial.prop,n.GH=1e2,logtau2.bounds=c(-20,3))$mle.full
    (rev(mles)-c(tau2.null,mu.null)) / sqrt(diag(vcov.hat))
})
op <- par(mfrow=c(1,2))
qqnorm(z.stats['mu',]); abline(0,1)
qqnorm(z.stats['tau2',]); abline(0,1)
par(op)



## 13 custom mle using newton-raphson instead of optim

## 13a gradient and hessian of loglik


n.GH <- 1e2
GH <- fastGHQuad::gaussHermiteData(n.GH)
m <- 15
yi <- 4
paramsi <- c(m=m,y.min=0,y.max=m)
eta <- function(z,mu,tau2)mu+z*sqrt(tau2)
deta.mu <- function(z,mu,tau2)1
d2eta.mu <- function(z,mu,tau2)0
deta.tau2 <- function(z,mu,tau2)z/2/sqrt(tau2)
d2eta.tau2 <- function(z,mu,tau2)-z/4/sqrt(tau2)^3
d2eta.mu.tau2 <- function(z,mu,tau2)0
f.eta <- function(eta)dbinom(yi,m,plogis(eta))
moment.1 <- function(paramsi,eta) paramsi['m']*plogis(eta)
moment.2 <- function(paramsi,eta)    paramsi['m']*plogis(eta)*plogis(eta,lower=FALSE)
df.eta <- function(eta)f.eta(eta)*(yi-moment.1(paramsi,eta))
d2f.eta <- function(eta)f.eta(eta)*((yi-moment.1(paramsi,eta))^2 - moment.2(paramsi,eta))
## moment.2(paramsi,z,mu,tau2)
## df.eta(4)
## I <- function(f.eta,mu,tau2){
##     rule <- list(z=sqrt(2)*GH$x,wt=GH$w/sqrt(pi),log.wt=log(GH$w/sqrt(pi)))
##     sum(rule$wt * f.eta(eta(rule$z,mu,tau2)) )
##     }
I <- function(f.z){
    rule <- list(z=sqrt(2)*GH$x,wt=GH$w/sqrt(pi),log.wt=log(GH$w/sqrt(pi)))
    sum(rule$wt * f.z(rule$z) )
}
mus <- seq(-10,10,len=1e2)
tau2s <- seq(0,10,len=1e2)
mu0 <- 0; tau20 <- .7
l <- function(mu,tau2)log(I(function(z)f.eta(eta(z,mu,tau2))))
## l(mu0,tau20)
plot(mus, sapply(mus,function(mu)l(mu,tau20)),type='l')
dl.mu <- function(mu,tau2)I(function(z)df.eta(eta(z,mu,tau2))*deta.mu(z,mu,tau2))/ I(function(z)f.eta(eta(z,mu,tau2)))
abline(a=l(mu0,tau20)-mu0*dl.mu(mu0,tau20),b= dl.mu(mu0,tau20))
plot(tau2s, sapply(tau2s,function(tau2)l(mu0,tau2)),type='l')
dl.tau2 <- function(mu,tau2)I(function(z)df.eta(eta(z,mu,tau2))*deta.tau2(z,mu,tau2)) / I(function(z)f.eta(eta(z,mu,tau2)))
abline(a=l(mu0,tau20)-tau20*dl.tau2(mu0,tau20),b= dl.tau2(mu0,tau20))
d2l.mu <- function(mu,tau2) (I(function(z)d2f.eta(eta(z,mu,tau2))*deta.mu(z,mu,tau2)^2 + df.eta(eta(z,mu,tau2))*d2eta.mu(z,mu,tau2))*I(function(z)f.eta(eta(z,mu,tau2))) - I(function(z)df.eta(eta(z,mu,tau2))*deta.mu(z,mu,tau2))^2 ) / I(function(z)f.eta(eta(z,mu,tau2)))^2
plot(mus, sapply(mus,function(mu)l(mu,tau20)),type='l')
approx <- function(mu)l(mu0,tau20) + dl.mu(mu0,tau20)*(mu-mu0) + 1/2*d2l.mu(mu0,tau20)*(mu-mu0)^2 
lines(mus, sapply(mus,approx),col=2)
d2l.tau2 <- function(mu,tau2) (I(function(z)d2f.eta(eta(z,mu,tau2))*deta.tau2(z,mu,tau2)^2 + df.eta(eta(z,mu,tau2))*d2eta.tau2(z,mu,tau2))*I(function(z)f.eta(eta(z,mu,tau2))) - I(function(z)df.eta(eta(z,mu,tau2))*deta.tau2(z,mu,tau2))^2 ) / I(function(z)f.eta(eta(z,mu,tau2)))^2
plot(tau2s, sapply(tau2s,function(tau2)l(mu0,tau2)),type='l')
approx <- function(tau2)l(mu0,tau20) + dl.tau2(mu0,tau20)*(tau2-tau20) + 1/2*d2l.tau2(mu0,tau20)*(tau2-tau20)^2 
lines(tau2s, sapply(tau2s,approx),col=2)
d2l.mu.tau2 <- function(mu,tau2) (I(function(z)d2f.eta(eta(z,mu,tau2))*deta.tau2(z,mu,tau2)*deta.mu(z,mu,tau2) + df.eta(eta(z,mu,tau2))*d2eta.mu.tau2(z,mu,tau2))*I(function(z)f.eta(eta(z,mu,tau2))) - I(function(z)df.eta(eta(z,mu,tau2))*deta.tau2(z,mu,tau2))^2 ) / I(function(z)f.eta(eta(z,mu,tau2)))^2
## d2l.mu.tau2(mu0,tau20)


## 13b reparametrize eta -> mu,tau2 to be compatible with current model
## routines

n.GH <- 1e2
GH <- fastGHQuad::gaussHermiteData(n.GH)
m <- 15
yi <- 4
paramsi <- c(m=m,y.min=0,y.max=m)
eta <- function(z,mu,tau2)mu+z*sqrt(tau2)
deta.mu <- function(z,mu,tau2)1
d2eta.mu <- function(z,mu,tau2)0
deta.tau2 <- function(z,mu,tau2)z/2/sqrt(tau2)
d2eta.tau2 <- function(z,mu,tau2)-z/4/sqrt(tau2)^3
d2eta.mu.tau2 <- function(z,mu,tau2)0
## f.eta <- function(eta)dbinom(yi,m,plogis(eta))
f.eta <- function(z,mu,tau2)dbinom(yi,m,plogis(eta(z,mu,tau2)))
## moment.1 <- function(paramsi,eta) paramsi['m']*plogis(eta)
## moment.2 <- function(paramsi,eta)    paramsi['m']*plogis(eta)*plogis(eta,lower=FALSE)
moment.1 <- function(z,mu,tau2,paramsi) paramsi['m']*plogis(eta(z,mu,tau2))
moment.2 <- function(z,mu,tau2,paramsi) {eta0 <- eta(z,mu,tau2);   paramsi['m']*plogis(eta0)*plogis(eta0,lower=FALSE)}
df.eta <- function(z,mu,tau2)f.eta(z,mu,tau2)*(yi-moment.1(z,mu,tau2,paramsi))
d2f.eta <- function(z,mu,tau2)f.eta(z,mu,tau2)*((yi-moment.1(z,mu,tau2,paramsi))^2 - moment.2(z,mu,tau2,paramsi))
I <- function(f.z){
    rule <- list(z=sqrt(2)*GH$x,wt=GH$w/sqrt(pi),log.wt=log(GH$w/sqrt(pi)))
    sum(rule$wt * f.z(rule$z) )
}
mus <- seq(-10,10,len=1e2)
tau2s <- seq(0,10,len=1e2)
mu0 <- 0; tau20 <- .7
l <- function(mu,tau2)log(I(function(z)f.eta(z,mu,tau2)))
## l(mu0,tau20)
plot(mus, sapply(mus,function(mu)l(mu,tau20)),type='l')
dl.mu <- function(mu,tau2)I(function(z)df.eta(z,mu,tau2)*deta.mu(z,mu,tau2))/ I(function(z)f.eta(z,mu,tau2))
abline(a=l(mu0,tau20)-mu0*dl.mu(mu0,tau20),b= dl.mu(mu0,tau20))
plot(tau2s, sapply(tau2s,function(tau2)l(mu0,tau2)),type='l')
dl.tau2 <- function(mu,tau2)I(function(z)df.eta(z,mu,tau2)*deta.tau2(z,mu,tau2)) / I(function(z)f.eta(z,mu,tau2))
abline(a=l(mu0,tau20)-tau20*dl.tau2(mu0,tau20),b= dl.tau2(mu0,tau20))
d2l.mu <- function(mu,tau2) (I(function(z)d2f.eta(z,mu,tau2)*deta.mu(z,mu,tau2)^2 + df.eta(z,mu,tau2)*d2eta.mu(z,mu,tau2))*I(function(z)f.eta(z,mu,tau2)) - I(function(z)df.eta(z,mu,tau2)*deta.mu(z,mu,tau2))^2 ) / I(function(z)f.eta(z,mu,tau2))^2
plot(mus, sapply(mus,function(mu)l(mu,tau20)),type='l')
approx <- function(mu)l(mu0,tau20) + dl.mu(mu0,tau20)*(mu-mu0) + 1/2*d2l.mu(mu0,tau20)*(mu-mu0)^2 
lines(mus, sapply(mus,approx),col=2)
d2l.tau2 <- function(mu,tau2) (I(function(z)d2f.eta(z,mu,tau2)*deta.tau2(z,mu,tau2)^2 + df.eta(z,mu,tau2)*d2eta.tau2(z,mu,tau2))*I(function(z)f.eta(z,mu,tau2)) - I(function(z)df.eta(z,mu,tau2)*deta.tau2(z,mu,tau2))^2 ) / I(function(z)f.eta(z,mu,tau2))^2
plot(tau2s, sapply(tau2s,function(tau2)l(mu0,tau2)),type='l')
approx <- function(tau2)l(mu0,tau20) + dl.tau2(mu0,tau20)*(tau2-tau20) + 1/2*d2l.tau2(mu0,tau20)*(tau2-tau20)^2 
lines(tau2s, sapply(tau2s,approx),col=2)
d2l.mu.tau2 <- function(mu,tau2) (I(function(z)d2f.eta(z,mu,tau2)*deta.tau2(z,mu,tau2)*deta.mu(z,mu,tau2) + df.eta(z,mu,tau2)*d2eta.mu.tau2(z,mu,tau2))*I(function(z)f.eta(z,mu,tau2)) - I(function(z)df.eta(z,mu,tau2)*deta.tau2(z,mu,tau2))^2 ) / I(function(z)f.eta(z,mu,tau2))^2
d2l.mu.tau2(mu0,tau20)



## 13c encapsulate
source('models.R')
source('./meta.glmm/R/meta.glmm.R')
n.GH <- 1e2
GH <- fastGHQuad::gaussHermiteData(n.GH)
eta <- function(z,mu,tau2)mu+z*sqrt(tau2)
deta.mu <- function(z,mu,tau2)1
d2eta.mu <- function(z,mu,tau2)0
deta.tau2 <- function(z,mu,tau2)z/2/sqrt(tau2)
d2eta.tau2 <- function(z,mu,tau2)-z/4/sqrt(tau2)^3
d2eta.mu.tau2 <- function(z,mu,tau2)0
moment.1 <- function(z,mu,tau2,paramsi) paramsi['m']*plogis(eta(z,mu,tau2))
moment.2 <- function(z,mu,tau2,paramsi) {eta0 <- eta(z,mu,tau2);   paramsi['m']*plogis(eta0)*plogis(eta0,lower=FALSE)}
moment.map <- moment.map.bn.prop
moment.map <- function(z, paramsi, mu, tau2) {
    logp <- plogis(mu+sqrt(tau2)*z,log.p=TRUE) |> unname()
    logp.bar <- plogis(mu+sqrt(tau2)*z,log.p=TRUE,lower.tail=FALSE) |> unname()
    out <- exp(structure(unname(c(log(paramsi['m'])+logp,  log(paramsi['m'])+logp+logp.bar)), names=c('first','second')))
    return(out)
}
mean.map <- function(z, paramsi, mu, tau2) {
    logp <- plogis(mu+sqrt(tau2)*z,log.p=TRUE) 
    out <- exp(log(paramsi['m'])+logp)
    out <- unname(out)
    return(out)
    }
var.map <- function(z, paramsi, mu, tau2) {
    logp <- plogis(mu+sqrt(tau2)*z,log.p=TRUE) 
    logp.bar <- plogis(mu+sqrt(tau2)*z,log.p=TRUE,lower.tail=FALSE) 
    out <- exp(log(paramsi['m'])+logp+logp.bar)
    out <- unname(out)
    return(out)
    }
cond.pmf <- cond.pmf.binomial.prop
## above to be supplied by model, applies to any study in the model
derivs.study <- function(yi,paramsi,mu,tau2) {
    ## browser()
    ## f.eta <- function(z,mu,tau2)dbinom(yi,paramsi['m'],plogis(eta(z,mu,tau2)))
    f.eta <- function(z,mu,tau2)cond.pmf(x=yi,theta=z,mu=mu,tau2=tau2,params=paramsi,log=FALSE)
    ## df.eta <- function(z,mu,tau2)f.eta(z,mu,tau2)*(yi-moment.1(z,mu,tau2,paramsi))
    df.eta <- function(z,mu,tau2)f.eta(z,mu,tau2)*(yi-mean.map(z,paramsi,mu,tau2))
    ## d2f.eta <- function(z,mu,tau2)f.eta(z,mu,tau2)*((yi-moment.1(z,mu,tau2,paramsi))^2 - moment.2(z,mu,tau2,paramsi))
    d2f.eta <- function(z,mu,tau2)f.eta(z,mu,tau2)*((yi-moment.1(z,mu,tau2,paramsi))^2 -var.map(z,paramsi,mu,tau2))
    I <- function(f.z){
        rule <- list(z=sqrt(2)*GH$x,wt=GH$w/sqrt(pi),log.wt=log(GH$w/sqrt(pi)))
        sum(rule$wt * f.z(rule$z) )
    }
    l <- function(mu,tau2)log(I(function(z)f.eta(z,mu,tau2)))
    dl.mu <- function(mu,tau2)I(function(z)df.eta(z,mu,tau2)*deta.mu(z,mu,tau2))/ I(function(z)f.eta(z,mu,tau2))
    dl.tau2 <- function(mu,tau2)I(function(z)df.eta(z,mu,tau2)*deta.tau2(z,mu,tau2)) / I(function(z)f.eta(z,mu,tau2))
    d2l.mu <- function(mu,tau2) (I(function(z)d2f.eta(z,mu,tau2)*deta.mu(z,mu,tau2)^2 + df.eta(z,mu,tau2)*d2eta.mu(z,mu,tau2))*I(function(z)f.eta(z,mu,tau2)) - I(function(z)df.eta(z,mu,tau2)*deta.mu(z,mu,tau2))^2 ) / I(function(z)f.eta(z,mu,tau2))^2
    d2l.tau2 <- function(mu,tau2) (I(function(z)d2f.eta(z,mu,tau2)*deta.tau2(z,mu,tau2)^2 + df.eta(z,mu,tau2)*d2eta.tau2(z,mu,tau2))*I(function(z)f.eta(z,mu,tau2)) - I(function(z)df.eta(z,mu,tau2)*deta.tau2(z,mu,tau2))^2 ) / I(function(z)f.eta(z,mu,tau2))^2
    d2l.mu.tau2 <- function(mu,tau2) (I(function(z)d2f.eta(z,mu,tau2)*deta.tau2(z,mu,tau2)*deta.mu(z,mu,tau2) + df.eta(z,mu,tau2)*d2eta.mu.tau2(z,mu,tau2))*I(function(z)f.eta(z,mu,tau2)) - I(function(z)df.eta(z,mu,tau2)*deta.tau2(z,mu,tau2))*I(function(z)df.eta(z,mu,tau2)*deta.mu(z,mu,tau2)) ) / I(function(z)f.eta(z,mu,tau2))^2
    objective <- l(mu,tau2)
    gradient <- c(dl.mu(mu,tau2),dl.tau2(mu,tau2))
    hessian <- {hess.11 <- d2l.mu(mu,tau2); hess.12 <- d2l.mu.tau2(mu,tau2); hess.22 <- d2l.tau2(mu,tau2); matrix(c(hess.11,hess.12,hess.12,hess.22),nrow=2)}
    return(list(objective=objective,gradient=gradient,hessian=hessian))
}
## ks.fun <- readRDS('ks-fun.rds')
## pmf.derivs.fun.bn.prop <- readRDS('pmf-derivs-fun-binomial-prop.rds')
stijnen.1 <- read.table('stijnen-table-I.txt')
y <- stijnen.1$Y0
ms <- stijnen.1$n0
params <- cbind(m=ms,y.min=0,y.max=ms)
## mu <- 0; tau2 <- 1
theta <- c(mu=0,tau2=1)
## theta <- c(mu=-4.8,tau2=.83)
## theta <- c(mu=-3.97,tau2=.24)
theta <- par.start(y,params,nn.transform.bn.prop)
## theta <- c(mu = -7.40722632338218, tau2 = 24.2754804611624)
tol <- 1e-7
for(i in 1:20) {
    derivs <- lapply(seq_along(y), function(i)derivs.study(y[i],params[i,],theta['mu'],theta['tau2']))
    grad <- sapply(derivs,`[[`,'gradient') |> rowSums()
    hess <- Reduce(`+`,lapply(derivs,`[[`,'hessian'))
    step <-  - solve(hess)%*%grad
    theta.new <-  theta + c(step)
    if(max(abs(theta-theta.new))<tol)break
    theta <- theta.new
    dput(theta)
}

theta.start <- par.start(y,params,nn.transform.bn.prop)
theta.start <- c(mu=-4.8,tau2=.83)
obj <- function(theta)sum(sapply(seq_along(y), function(i)derivs.study(y[i],params[i,],theta['mu'],theta['tau2'])$obj))
optim(theta.start,obj,method='Nelder')

require(nloptr)
obj <- function(theta)-sum(sapply(seq_along(y), function(i)derivs.study(y[i],params[i,],theta[1],theta[2])$obj))
fit <- nloptr(
    x0     = c(0, 1),
    eval_f = obj,
    lb     = c(-Inf, 0),
    ub     = c( Inf, Inf),
    opts   = list(
        algorithm = "NLOPT_LN_BOBYQA",
        xtol_abs  = 1e-8,
        ftol_abs  = 1e-8,
        maxeval   = 1e5
    )
)
mu.hat   <- fit$solution[1]
tau.hat  <- fit$solution[2]
tau2.hat <- tau.hat^2
c(mu = mu.hat, tau2 = tau2.hat, logLik = -fit$objective)
res_mle <- metafor::rma.glmm(
  measure = "PLO", 
  xi = Y1, 
  ni = n1, 
  data = stijnen.1,
  model = "UM.FS", # Unconditional model with fixed study effects / random intercept
  method = "ML"    # Maximum Likelihood
)
# Extract mu (beta) and tau2
mle_estimates <- c(
  mu   = as.numeric(res_mle$beta),
  tau2 = as.numeric(res_mle$tau2)
)
print(mle_estimates)



## 13d speed up code, apply to bn-prop rate sim. update: newton keeps overshooting. cant get through a sim. moving to nloptr.

source('models.R')
source('./meta.glmm/R/meta.glmm.R')
## stijnen.1 <- read.table('stijnen-table-I.txt')
## y <- stijnen.1$Y0
## ms <- stijnen.1$n0
## params <- cbind(m=ms,y.min=0,y.max=ms)
n.GH <- 5e2
GH <- fastGHQuad::gaussHermiteData(n.GH)
eta.map <- function(z, paramsi, mu, tau2) {
    eta <- mu+z*sqrt(tau2)
    grad <- cbind(mu=1,tau2=z/2/sqrt(tau2))
    hessian <- cbind(mu2=0, mu.tau2=0, tau22=-z/4/sqrt(tau2)^3)
    return(list(eta=eta,grad=grad,hessian=hessian))
}
cond.pmf <- cond.pmf.binomial.prop
moment.map <- function(z, paramsi, mu, tau2) {
    logp <- plogis(mu+sqrt(tau2)*z,log.p=TRUE) |> unname()
    logp.bar <- plogis(mu+sqrt(tau2)*z,log.p=TRUE,lower.tail=FALSE) |> unname()
    cbind(log(paramsi['m'])+logp, +logp.bar)
    out <- log(paramsi['m'])+logp
    out <- unname(out)
    out <- exp( cbind(mean=out,var=out+logp.bar) )
    ## out <- exp(structure(unname(c(log(paramsi['m'])+logp,  log(paramsi['m'])+logp+logp.bar)), names=c('first','second')))
    return(out)
}
## max(abs(moment.map(rule$z,paramsi,mu,tau2)[,'mean']-mean.map(rule$z,paramsi,mu,tau2)))
## yi <- y[1]
## paramsi <- params[1,]
## mu <- 0
## tau2 <- 1
derivs.study <- function(yi,paramsi,mu,tau2) {
    rule <- list(z=sqrt(2)*GH$x,wt=GH$w/sqrt(pi),log.wt=log(GH$w/sqrt(pi)))
    moments <- moment.map(rule$z,paramsi,mu,tau2)
    eta.derivs <- eta.map(rule$z,paramsi,mu,tau2)
    f.eta.z <- cond.pmf(x=yi,theta=rule$z,mu=mu,tau2=tau2,params=paramsi,log=FALSE)
                                        #f.eta(rule$z,mu,tau2)
    ## df.eta.z <- df.eta(rule$z,mu,tau2)
    df.eta.z <- f.eta.z*(yi-moments[,'mean'])
    deta.mu.z <- eta.derivs$grad[,'mu']#deta.mu(rule$z,mu,tau2)
    deta.tau2.z <- eta.derivs$grad[,'tau2']#deta.tau2(rule$z,mu,tau2)
    ## d2f.eta.z <- d2f.eta(rule$z,mu,tau2)
    d2f.eta.z <- f.eta.z*( (yi-moments[,'mean'])^2 - moments[,'var'] )
    d2eta.mu.z <- eta.derivs$hessian[,'mu2']#d2eta.mu(rule$z,mu,tau2)
    d2eta.tau2.z <- eta.derivs$hessian[,'tau22']#d2eta.tau2(rule$z,mu,tau2)
    d2eta.mu.tau2.z <- eta.derivs$hessian[,'mu.tau2']#d2eta.mu.tau2(rule$z,mu,tau2)
    I.f <- f.eta.z%*%rule$wt#I(function(z)f.eta(z,mu,tau2))
    I.f.mu <- rule$wt%*%(df.eta.z*deta.mu.z)#I(function(z)df.eta(z,mu,tau2)*deta.mu(z,mu,tau2))
    I.f.tau2 <- rule$wt%*%(df.eta.z*deta.tau2.z)#I(function(z)df.eta(z,mu,tau2)*deta.tau2(z,mu,tau2))
    dl.mu <- I.f.mu / I.f
    dl.tau2 <- I.f.tau2 / I.f
    d2l.mu <- (rule$wt%*%(d2f.eta.z*deta.mu.z^2 + df.eta.z*d2eta.mu.z) * I.f - I.f.mu^2 ) / I.f^2
    d2l.tau2 <- ( rule$wt%*%(d2f.eta.z*deta.tau2.z^2 + df.eta.z*d2eta.tau2.z) *I.f - I.f.tau2^2 ) /I.f^2
    d2l.mu.tau2 <- (rule$wt%*%(d2f.eta.z*deta.tau2.z*deta.mu.z + df.eta.z*d2eta.mu.tau2.z)*I.f - I.f.mu*I.f.tau2 ) / I.f^2
    ## objective <- l(mu,tau2)
    gradient <- c(dl.mu,dl.tau2)
    hessian <- matrix(c(d2l.mu,d2l.mu.tau2,d2l.mu.tau2,d2l.tau2),nrow=2)
    list(gradient=gradient,hessian=hessian)
}
data <- structure(c(3, 0, 6, 1, 2, 6, 4, 5, 136, 177, 117, 35, 64, 180,69, 117), dim = c(8L, 2L), dimnames = list(NULL, c("y", "m")))
y <- data[,'y']; m <- data[,'m']
params <- cbind(m=m,y.min=0,y.max=m)
tol <- 1e-7
theta <- par.start(y,params,nn.transform.bn.prop)
for(i in 1:20) {
    derivs <- lapply(seq_along(y), function(i)derivs.study(y[i],params[i,],theta['mu'],theta['tau2']))
    grad <- sapply(derivs,`[[`,'gradient') |> rowSums()
    hess <- Reduce(`+`,lapply(derivs,`[[`,'hessian'))
    print(hess)
    step <-  - solve(hess)%*%grad
    theta.new <-  theta + c(step)
    ## if(theta.new['tau2']<0)browser()
        print(theta.new)
    while(theta.new['tau2']<0) {
        print(theta.new)
        step <- step/2
        theta.new  <- theta+c(step)
        print(theta.new)
    }
    if(max(abs(theta-theta.new))<tol)break
    theta <- theta.new
    ## dput(theta)
}

try.mle <- function(y,params,par.start, tol=1e-7){
    ## browser()
    theta <- par.start
    for(i in 1:20) {
        derivs <- lapply(seq_along(y), function(i)derivs.study(y[i],params[i,],theta['mu'],theta['tau2']))
        grad <- sapply(derivs,`[[`,'gradient') |> rowSums()
        hess <- Reduce(`+`,lapply(derivs,`[[`,'hessian'))
        step <-  - solve(hess)%*%grad
        theta.new <-  theta + c(step)
        ## if(theta.new['tau2']<0)browser()
        while(theta.new['tau2']<0) {
            step <- step/2
            theta.new  <- theta+c(step)
            ## print(theta.new)
            }
        if(max(abs(theta-theta.new))<tol)break
        theta <- theta.new
        ## dput(theta)
    }
    return(theta)
}

source('models.R')
source('utils.R')
source('./meta.glmm/R/meta.glmm.R')
## stijnen.1 <- read.table('stijnen-table-I.txt')
## y <- stijnen.1$Y1
## ms <- stijnen.1$n1
## params <- cbind(m=ms,y.min=0,y.max=ms)
## ## mu <- 0; tau2 <- 1
## theta <- c(mu=0,tau2=1)
## ## theta <- c(mu=-4.8,tau2=.83)
## ## theta <- c(mu=-3.97,tau2=.24)
## par.start0 <- par.start(y,params,nn.transform.bn.prop)
## mle(y,params,par.start0)
## set.seed(1)
skew <- kurt <- 0
ks.fun <- readRDS('ks-fun.rds')
pmf.derivs.fun.bn.prop <- readRDS('pmf-derivs-fun-binomial-prop.rds')
sas.params <- sas.invert.params(skew,kurt)
re.distr <- function(n)rsas(n,eps=sas.params['eps'],delta=sas.params['delta'])
## RNGkind("L'Ecuyer-CMRG")#!!!!!!
## set.seed(1)#!!!!!!!!!
## re.distr <- rnorm #!!!!
## tau2.null <- 1/2
ns <- 5:15
mu.null <- -3.5
tau2.null <- .81
model.str <- 'bn-prop'
n.GH <- 3e2
stijnen.1 <- read.table('stijnen-table-I.txt')
y <- stijnen.1$Y0
ms <- stijnen.1$n0
stijnen.params <- cbind(m=ms,y.min=0,y.max=ms)
mc.idx <- sample(nrow(stijnen.params),max(ns),replace=TRUE)
mc.params <- stijnen.params[mc.idx,]
model0 <- model(model.str)
quad <- make.quad(model0$moment.map,model0$cond.pmf,n.GH=n.GH,adaptive=FALSE,adaptive.tol=1e-5)
start <- Sys.time()
by.n <- parallel::mclapply(seq_along(ns),mc.cores=6,FUN= function(i) {
## by.n <- sapply(seq_along(ns), function(i) {
    ## repeat{
    cat('.')
    n <- ns[i]
    params <- mc.params[1:n,]
    bartlett.C <- make.bartlett.C(params,quad,model0$pmf.derivs.fun)
    bartlett.C0 <- bartlett.C(mu.null,tau2.null)$bartlett.C
    stats <- replicate(round(n^3), {
        ## browser()
        y <- model0$dgp(n,mu.null,tau2.null,params,re.distr=re.distr)
        ## dput(cbind(y,params[,'m']))
        loglik <- make.loglik(model0$cond.pmf,params,y,quad)
        ## logtau2.mle.null <- optimize(function(logtau2)loglik(mu.null,exp(logtau2)),interval=c(-20,3),maximum=TRUE)$max
        ## tau2.mle.null <- exp(logtau2.mle.null)
        ## par.start0 <- c(mu=mu.null+rnorm(1)/sqrt(n),tau2=tau2.mle.null)
        ## par.start0 <- par.start(y,params,nn.transform.bn.prop)
        ## mle0 <- try.mle(y,params,par.start0)
        ## print(mle0)
        ## mle0 <- mle(loglik,par.start0,method='BFGS')$mle
        obj <- function(par)-loglik(par[1],par[2])
        fit <- nloptr::nloptr(
                           x0     = c(0, 1),
                           eval_f = obj,
                           lb     = c(-Inf, 0),
                           ub     = c( Inf, Inf),
                           opts   = list(
                               algorithm = "NLOPT_LN_BOBYQA",
                               xtol_abs  = 1e-8,
                               ftol_abs  = 1e-8,
                               maxeval   = 1e5
                           )
                       )
        mle0 <- c(mu=fit$solution[1],tau2=fit$solution[2])
        profile.LR <- make.profile.LR(loglik,mle0,par.start0,logtau2.bounds=c(-20,3))
        LR <- profile.LR(mu.null)$LR
        LR.bartlett <- LR / (1+bartlett.C0)
        c(LR=LR,LR.bartlett=LR.bartlett)
    })
    ecdfs <- apply(stats,1,function(x)ecdf(pchisq(x,df=1,lower=FALSE)))
    ks <- sapply(ecdfs,function(ecdf)max(abs(ecdf(knots(ecdf)) - knots(ecdf))))
    means <- rowMeans(stats-1)
    c(ks=ks,means=means)
})
by.n <- simplify2array(by.n)
print(Sys.time() - start)
plot(ns,by.n['ks.LR',],type='l',ylim=range(by.n[1:2,]))
lines(ns,by.n['ks.LR.bartlett',],col=2)
apply(by.n[1:2,],1,function(ks)coef(lm(log(ks) ~ log(ns)))[2]) |> print()



## y <- c(0L, 0L, 3L, 1L, 0L)
## params <- structure(c(174, 345, 113, 345, 44, 177, 362, 122, 362, 35, 4, 
## 4, 3, 4, 4, 0, 0, 0, 0, 0, 4, 4, 3, 4, 4), dim = c(5L, 5L), dimnames = list(
##     NULL, c("m.T", "m.C", "s", "y.min", "y.max")))





## ## 14. diagnostics


## ## nonadaptive quadrature
## ks.fun <- readRDS('ks-fun.rds')
## pmf.derivs.fun.pn <- readRDS('pmf-derivs-fun-poisson.rds')
## source('models.R')
## source('./meta.glmm/R/meta.glmm2.R')
## stijnen.vi <- read.table('stijnen-table-VI.txt')
## n.GH <- 1e2
## y <- stijnen.vi$Y1
## T <- stijnen.vi$T1/1000
## params <- cbind(T=T,y.min=0,y.max=1e4)
## n.GH <- seq(100,500,by=200)
## y.max <- seq(1e3,1.3e4,by=3e3)
## pars <- expand.grid(n.GH=n.GH,y.max=y.max)
## ## CIs <- apply(pars, 1, function(par) {
## par <- unlist(pars[13,])
## print(par)
## n.GH <- par['n.GH']; y.max <- par['y.max']
## n.GH <- 20;y.max <- 1e3
## params <- cbind(T=T,y.min=0,y.max=y.max)


## ## bn-prop
## ks.fun <- readRDS('ks-fun.rds')
## pmf.derivs.fun.bn.prop <- readRDS('pmf-derivs-fun-binomial-prop.rds')
## source('models.R')
## source('./meta.glmm/R/meta.glmm2.R')
## stijnen.1 <- read.table('stijnen-table-I.txt')
## y <- stijnen.1$Y1
## ms <- stijnen.1$n1
## params <- cbind(m=ms,y.min=0,y.max=ms)
## fit <- meta.glmm(y,params,'bn-prop',n.GH=n.GH,adaptive=FALSE,logtau2.bounds=c(-20,3),debug=TRUE)
## CI <- confint.meta.glmm(fit,mc.cores=parallel::detectCores()-1,mu.grid=list(offset=2.5,size=10),plot=TRUE)
## mle.null <- fit$profile.LR( -6)$mle.null
## bartlett.fit <- fit$bartlett.C(mle.null['mu'],mle.null['tau2'])

## ## mle.null <- fit$profile.LR( -1.1726610)$mle.null
## ## fit <- fit$bartlett.C(mle.null['mu'],mle.null['tau2'])



## ## pn CIs--adaptive
## ## ks.fun <- readRDS('ks-fun.rds')
## ## pmf.derivs.fun.pn <- readRDS('pmf-derivs-fun-poisson.rds')
## ## source('models.R')
## ## source('./meta.glmm/R/meta.glmm.R')
## devtools::load_all('meta.glmm')
## stijnen.vi <- read.table('stijnen-table-VI.txt')
## n.GH <- 1e2
## y <- stijnen.vi$Y1
## T <- stijnen.vi$T1/1000
## params <- cbind(T=T,y.min=0,y.max=1e4)
## ## n.GH <- c(seq(10,20,by=5),1e2)
## ## y.max <- c(seq(5e2,1.5e3,by=5e2),5e3)#seq(1e2,3e2,by=50)
## n.GH <- 20
## y.max <- 1e2
## params <- cbind(T=T,y.min=0,y.max=y.max)
## fit <- meta.glmm(y,params,'pn',n.GH=n.GH,adaptive=TRUE,logtau2.bounds=c(-20,3),debug=TRUE)
## CI <- confint.meta.glmm(fit,mc.cores=parallel::detectCores()-1,mu.grid=list(offset=1.5,size=10),plot=TRUE)
## mle.null <- fit$profile.LR( -1.1726610)$mle.null
## bartlett.fit <- fit$bartlett.C(mle.null['mu'],mle.null['tau2'])

## diagnostic.pmf(bartlett.fit,3)
## diagnostic.k(bartlett.fit,3)











## look at whats going wrong with bn-prop sim
RNGkind("L'Ecuyer-CMRG")
set.seed(1)
devtools::load_all('meta.glmm')
source('utils.R') ## for sinh-asinh reparameterization
sas.params <- sas.invert.params(target.skew=0,target.kurt=0)
re.distr <- function(n)rsas(n,eps=sas.params['eps'],delta=sas.params['delta'])
ns <- 5:15
mu.null <- -3.5
tau2.null <- .81
model.str <- 'bn-prop'
n.GH <- 20
stijnen.params <- with(nw07, cbind(m=n0,y.min=0,y.max=n0))
mc.idx <- sample(nrow(stijnen.params),max(ns),replace=TRUE)
mc.params <- stijnen.params[mc.idx,]
model0 <- model(model.str)
quad <- make.quad(model0$moment.map,model0$cond.pmf,n.GH=n.GH,adaptive=TRUE,adaptive.tol=1e-6)
start <- Sys.time()
## by.n <- parallel::mclapply(seq_along(ns),mc.cores=6,FUN= function(i) {
by.n <- sapply(seq_along(ns), function(i) {
    ## repeat{
    cat('.')
    n <- ns[i]
    params <- mc.params[1:n,]
    bartlett.C <- make.bartlett.C(params,quad,model0$pmf.derivs.fun)
    bartlett.C0 <- bartlett.C(mu.null,tau2.null)$bartlett.C
    stats <- replicate(round(n^2), {        
        y <- model0$dgp(n,mu.null,tau2.null,params,re.distr=re.distr)
        fit <- meta.glmm(y,params,model.str,n.GH,adaptive=FALSE,logtau2.bounds=c(-20,3))        
        LR <- fit$profile.LR(mu.null)$LR
        LR.bartlett <- LR / (1+bartlett.C0)
        c(LR=LR,LR.bartlett=LR.bartlett)
    })
    ecdfs <- apply(stats,1,function(x)ecdf(pchisq(x,df=1,lower=FALSE)))
    ks <- sapply(ecdfs,function(ecdf)max(abs(ecdf(knots(ecdf)) - knots(ecdf))))
    means <- rowMeans(stats-1)
    c(ks=ks,means=means)
})
by.n <- simplify2array(by.n)
print(Sys.time() - start)
plot(ns,by.n['ks.LR',],type='l',ylim=range(by.n[1:2,]))
lines(ns,by.n['ks.LR.bartlett',],col=2)
apply(by.n[1:2,],1,function(ks)coef(lm(log(ks) ~ log(ns)))[2]) |> print()

