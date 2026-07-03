## simple logistic regression using Jensen '93 formula. seems that
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



## fuller sim to see effect of # of x support points affects rate
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


## laplace approx

## cox & barndorff-nielson

## normal approx--O(1/z^5) error as expected
z <- seq(0,3,len=20)
true <- 1-pnorm(z)
try <- dnorm(z)*(1/z-1/z^3)
plot(true-try)
coef(lm(log(abs(true-try))~log(z)))[2]


## gamma function--relative error 1/z^-1.5, as expected
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


## normal-binomial likelihood
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


## compute MLE's
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



## reparametrize in terms of standard normal theta
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


## read in a likelihood formula exported from mma
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





## encapsulate, check CIs. match metafor's pretty well.
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


## get all the full model k terms

source('utils.R')
lawley <- function(p,m,beta,mma.strings) {
    k.vals <- sapply(mma.strings, function(mma.str) {
        l.deriv <- str.to.fun(
            gsub(pattern=c('m\\[j\\]'),replacement=c('m'),mma.str),
            args = c('theta','m','k'),
            binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt,m=m,beta=beta)
        )
        ## browser()
        sapply(1:length(m), function(j) 
            sapply(1:(m[j]-1), function(k) {
                ## browser()
                integrate(function(theta)l.deriv(theta,m[j],k)*dnorm(theta),-Inf,Inf)$val
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
        inv.inf[r,s] * inv.inf[t,u] * inv.inf[v,w] * (1/6*k(c(r,t,v))*k(c(s,u,w)) + 1/4*k(c(r,t,u))*k(c(s,u,w)) - k(c(r,t,v))*k(c(s,w),u) - k(c(r,t,u))*k(c(s,w),v) + k(c(r,t),v)*k(c(s,w),u) + k(c(r,t),u)*k(c(s,w),v))
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
}
set.seed(1)
n <- 20
mu.null <- 1
tau2.null <- 2
m <- rpois(n,lambda=20)
## theta <- rnorm(n)
## e <- rbinom(n,m,plogis(mu.null+sqrt(tau2.null)*theta))
## rm(theta)
beta <- runif(2)#c(mle[c('tau2','mu')]) ## should be passed to make.fun
eps.full <- lawley(p=2,m=m,beta=beta,mma.strings=read.mma('./mma-formulas/full/'))
eps.redu <- lawley(p=1,m=m,beta=beta,mma.strings=read.mma('./mma-formulas/redu/'))
correction <- 1 + eps.full - eps.redu
correction
