## integrand.mode.poisson <- function(yi, paramsi, mu, tau2, tol=1e-10) {
##     ## browser()
##     z <- 0
##     for (j in 1:20) {
##         lambda <- paramsi["T"] *
##             exp(mu + sqrt(tau2) * z)
##         f <- sqrt(tau2) * (yi - lambda) - z
##         fp <- -1 - tau2 * lambda
##         z.new <- z - f / fp
##         if (abs(z.new - z) < tol) {
##             break
##         }
##         z <- z.new
##     }
##     mode <- z.new
##     return(mode)
## }
## par.start.poisson <- function(x,params) {
##     x.cc <- x; x.cc[x.cc=0L] <-  x.cc[x.cc=0L]+.5
##     y.nn <- log(x.cc)/params[,'T']
##     sigma2.nn <- 1/x.cc
##     w.fe <- 1/sigma2.nn
##     mu.fe <- y.nn %*% w.fe / sum(w.fe)    
##     Q <- w.fe%*%(y.nn-c(mu.fe))^2
##     tau2.DL <- max(0L, (Q - (length(y.nn)-1)) / (sum(w.fe)-sum(w.fe^2)/sum(w.fe)))
##     tau2.start <- if(tau2.DL>0L)tau2.DL else var(y.nn)/1e4
##     return(c(mu=mu.fe,tau2=tau2.start))
## }


## .cond.pmf.hyperg <- function(x,theta,mu,tau2,params) {
##     ## x <- 0; theta <-  -18.95964; mu <- 0; tau2 <- 0.3071743;params <- c(m.T=116,m.C=117,s=3,x.min=0,x.max=3)
##     m.T <- params[,'m.T']
##     m.C <- params[,'m.C']
##     s <- params[,'s']
##     odds <- exp(mu+sqrt(tau2)*theta)
##     ## print(x);print(s);print(params)
##     ## browser()
##     ## print(c(m.T,m.C,s,odds))
##         BiasedUrn::dFNCHypergeo(x=x,m1=m.T,m2=m.C,n=s,odds=odds)            
## }
## cond.pmf.hyperg <-  Vectorize(.cond.pmf.hyperg,vectorize.args=c("theta"))

## cond.pmf.hyperg <- function(x,theta,mu,tau2,params)  <- sapply(theta, function(thetai) _cond.pmf.hyperg(x,thetai,mu,tau2,params)

## my version. dotted version matches biasedUrn's call syntax.
## check liao rosen. rewrite using logsumexp.
.cond.pmf.hyperg <- function(x,m1,m2,n,odds,log=FALSE) {
    ## print(c(length(x),length(m1),length(m2),length(n),length(odds)))
    ## if(length(x)==6&&length(odds)==10)browser()
    ## m.T <- params[,'m.T']
    ## m.C <- params[,'m.C']
    ## s <- params[,'s']
    ## print(c(length(odds),length(x)))
    numer <- choose(m1,x)*choose(m2,n-x)*odds^x
    supp.x <- max(0,n-m2):min(n,m1)
    denom <- rowSums(  outer(odds,supp.x,function(oddsj,xj)choose(m1,xj)*choose(m2,n-xj)*oddsj^xj ) )
    ## numer / denom
    out <- if(log) log(numer)-log(denom) else numer / denom #!!!!!!
}
cond.pmf.hn <- cond.pmf.hyperg <- function(x,theta,mu,tau2,params,log=FALSE) {
    if(is.null(dim(params))) {
        m.T <- params['m.T'];m.C <- params['m.C'];   s <- params['s']
    } else {
        m.T <- params[,'m.T'];m.C <- params[,'m.C']; s <- params[,'s']
        }
    ## m.T <- params[,'m.T']
    ## m.C <- params[,'m.C']
    ## s <- params[,'s']
    odds <- exp(mu+sqrt(tau2)*theta)
    ## browser()
    out <- .cond.pmf.hyperg(x=x,m1=m.T,m2=m.C,n=s,odds=odds,log=log)
    ## print(out)
    return(out)
}
## ## safer but slower:
## cond.pmf.hyperg <- function(x,theta,mu,tau2,params) {
##     m.T <- params[,'m.T']; m.C <- params[,'m.C']; s <- params[,'s']
##     logpmf.part <- function(xi)lgamma(m.T+1) - lgamma(xi+1) -lgamma(m.T-xi+1) + lgamma(m.C+1) - lgamma(s-xi+1) - lgamma(m.C-s+xi+1)
##     lognumer <- logpmf.part(x)+ (mu+sqrt(tau2)*theta)*x
##     ## logpmf.mat <- logpmf.part(0:m.T) + outer(0:m.T,theta)*sqrt(tau2) + (0:m.T)*mu
##     supp <- max(0,s-m.C):min(s,m.T)
##     logpmf.mat <- logpmf.part(supp) + outer(supp,theta)*sqrt(tau2) + supp*mu
##     colmaxes <- apply(logpmf.mat,2,max)        
##     logdenom <- colmaxes + log(rowSums(exp(t(logpmf.mat) - colmaxes)))
##     exp(lognumer - logdenom)
## }
## n <- sample(2:10,1)
## x <- sample(1:n,1)
## m1 <- 8; m2 <- 10
## odds <- rchisq(1,df=5)
## BiasedUrn::dFNCHypergeo(x=x,m1=m1,m2=m2,n=n,odds=odds)
## .cond.pmf.hyperg(x,m1,m2,n,odds)
## cond.pmf.hyperg <- function(x,theta,mu,tau2,params) {
##     m.T <- params[,'m.T']
##     m.C <- params[,'m.C']
##     s <- params[,'s']
##     numer <- choose(m.T,x)*choose(m.C,s-x)*exp((mu+sqrt(tau2)*theta)*x)
##     supp <- max(0,s-m.C):min(s,m.T)
##     denom <- rowSums(  outer(theta,supp,function(thetaj,xj)choose(m.T,xj)*choose(m.C,s-xj)*exp((mu+sqrt(tau2)*thetaj)*xj) ) )
##     numer / denom
## }
dgp.hn <- dgp.hyperg <- function(n,mu,tau2,params,re.distr=rnorm) {
    m.T <- params[,'m.T']
    m.C <- params[,'m.C']
    s <- params[,'s']
    theta <- re.distr(n)
    odds <- exp(mu+sqrt(tau2)*theta)
    x <- x.T <- sapply(1:n, function(i) BiasedUrn::rFNCHypergeo(1,m1=m.T[i],m2=m.C[i],n=s[i],odds=odds[i]))
}
nn.transform.hn <- function(y,params) {
    cells <- data.frame(y.1.T=y, y.1.C =params[,'s']-y)
    cells <- within(cells, {
        y.0.T <- params[,'m.T']-y.1.T;
        y.0.C <- params[,'m.C']-y.1.C
    })
    idx.cc <- apply(cells,1,function(v)any(v==0L))
    cells[idx.cc,] <- cells[idx.cc,]+.5
    y.nn <- with(cells, log(y.1.T/y.0.T) - log(y.1.C/y.0.C))
    sigma2.nn <- with(cells, 1/y.1.T + 1/y.0.T + 1/y.1.C + 1/y.0.C)
    return(data.frame(y=y.nn,sigma2=sigma2.nn))
}
moment.map.hn <- function(zi, paramsi, mu, tau2) {
    eta <- mu+sqrt(tau2)*zi
    supp.y <- paramsi['y.min']:paramsi['y.max']
    log.pr <- lchoose(paramsi['m.T'],supp.y) + lchoose(paramsi['m.C'],paramsi['s']-supp.y) + supp.y*eta
    log.pr <- log.pr - max(log.pr)
    pr <- exp(log.pr); pr <- pr/sum(pr)
    fnchy.mean <- sum(pr*supp.y)
    fnchy.var <- sum(pr*(supp.y-fnchy.mean)^2)
    out <- c(first=fnchy.mean,second=fnchy.var)
    }
## newton.step.hn <- function(zi, yi, paramsi, mu, tau2) {
##     ## browser()
##     ## eta <- mu+sqrt(tau2)*zi
##     ## supp.y <- paramsi['y.min']:paramsi['y.max']
##     ## log.pr <- lchoose(paramsi['m.T'],supp.y) + lchoose(paramsi['m.C'],paramsi['s']-supp.y) + supp.y*eta
##     ## log.pr <- log.pr - max(log.pr)
##     ## pr <- exp(log.pr); pr <- pr/sum(pr)
##     ## fnchy.mean <- sum(pr*supp.y)
##     ## fnchy.var <- sum(pr*(supp.y-fnchy.mean)^2)
##     moments <- moment.map.hn(zi,paramsi,mu,tau2)
##     deriv <- sqrt(tau2)*(yi-moments['first'])-zi
##     hess  <- -tau2*moments['second']-1
##     out <-  - deriv / hess    
##     ## omega <- exp(mu+sqrt(tau2)*zi)
##     ## supp.y <- paramsi['y.min']:paramsi['y.max']
##     ## P <- sapply(0:2, function(m) sum(choose(paramsi['m.T'],supp.y) * choose(paramsi['m.C'],paramsi['s']-supp.y)*omega^supp.y*supp.y^m))
##     ## deriv <- sqrt(tau2)*(yi-P[1+1]/P[0+1])-zi
##     ## hess  <- -tau2*(P[2+1]/P[0+1]-(P[1+1]/P[0+1])^2)-1
##     ## out <-  - deriv / hess
##     return(out)
## }
## var.map.hn <- function(zi,paramsi,mu,tau2) {
##     omega <- exp(mu+sqrt(tau2)*zi)
##     supp.y <- paramsi['y.min']:paramsi['y.max']
##     P <- sapply(0:2, function(m) sum(choose(paramsi['m.T'],supp.y) * choose(paramsi['m.C'],paramsi['s']-supp.y)*omega^supp.y*supp.y^m))
##     hess  <- P[2+1]/P[0+1]-(P[1+1]/P[0+1])^2
## }

cond.pmf.pn <- cond.pmf.poisson <- function(x,theta,mu,tau2,params,log=FALSE) {
    T <- if(is.null(dim(params))) params['T'] else params[,'T']
    lambda <- T*exp(mu+sqrt(tau2)*theta)
    dpois(x,lambda=lambda,log=log)
}
dgp.pn <- dgp.poisson <- function(n,mu,tau2,params,re.distr=rnorm) {
    T <- params[,'T']
    theta <- re.distr(n)
    lambda <- T*exp(mu+sqrt(tau2)*theta)
    rpois(n,lambda)
    ## pr <- function(theta)plogis( mu+sqrt(tau2)*theta)
    ## rbinom(n,size=m,prob=pr(theta))
}
nn.transform.pn <- function(y,params) {
    y.cc <- y; y.cc[y.cc=0L] <-  y.cc[y.cc=0L]+.5
    y.nn <- log(y.cc)/params[,'T']
    sigma2.nn <- 1/y.cc
    return(data.frame(y=y.nn,sigma2=sigma2.nn))
}
## cond.pmf.derivs.poisson <- function(zi, yi, paramsi, mu, tau2) {
##     lambda <- paramsi["T"] * exp(mu + sqrt(tau2) * zi)
##     deriv <- sqrt(tau2) * (yi - lambda) - zi
##     hess <- -1 - tau2 * lambda
##     return(structure(unname(c(deriv,hess)), names=c('first','second')))
## }
## newton.step.pn <- function(zi, yi, paramsi, mu, tau2) {
##     lambda <- paramsi["T"] * exp(mu + sqrt(tau2) * zi)
##     out <- (sqrt(tau2)*(yi-lambda) - zi) / (tau2*lambda+1)
##     return(out)
## }
## var.map.pn <- function(zi,paramsi,mu,tau2)paramsi['T']*exp(mu+sqrt(tau2)*zi)
moment.map.pn <- function(zi, paramsi, mu, tau2) {
    lambda <- paramsi["T"] * exp(mu + sqrt(tau2) * zi)
    out <- structure(unname(c(lambda,lambda)), names=c('first','second'))
    return(out)
    }

cond.pmf.bn.prop <- cond.pmf.binomial.prop <- function(x,theta,mu,tau2,params,log=FALSE) {
    ## m <- params[,'m']
    m <- if(is.null(dim(params))) params['m'] else params[,'m']
    pr <- function(theta)plogis(mu+sqrt(tau2)*theta)
    dbinom(x,size=m,prob=pr(theta),log=log)
    ## integrand <- function(theta)pr(theta)^xi*(1-pr(theta))^(mi-xi)
}
dgp.bn.prop <- dgp.binomial.prop <- function(n,mu,tau2,params,re.distr=rnorm) {
    m <- params[,'m']
    theta <- re.distr(n)
    pr <- function(theta)plogis( mu+sqrt(tau2)*theta)
    rbinom(n,size=m,prob=pr(theta))
}
nn.transform.bn.prop <- function(y,params) {
    y.1 <- y
    y.0 <- params[,'m'] - y
    idy.cc <- y.0==0 | y.1==0
    y.0[idy.cc] <- y.0[idy.cc]+.5
    y.1[idy.cc] <- y.1[idy.cc]+.5
    y.nn <- log(y.1/y.0)
    sigma2.nn <- 1/y.0 + 1/y.1
    return(data.frame(y=y.nn,sigma2=sigma2.nn))
}
## newton.step.bn.prop <- function(zi, yi, paramsi, mu, tau2) {
##     ## browser()
##     p <- plogis(mu+sqrt(tau2)*zi)
##     out <- (sqrt(tau2)*(yi - paramsi['m']*p)-zi) / (paramsi['m']*tau2*p*(1-p)+1)
##     ## deriv <- sqrt(tau2)*(yi - paramsi['m']*p)
##     ## hess <- -paramsi['m']*tau2*p*(1-p)
##     ## if(is.nan(deriv)||is.nan(hess))browser()
##     ## return(structure(unname(c(deriv,hess)), names=c('first','second')))
##     return(out)
## }
## var.map.bn.prop <- function(zi,paramsi,mu,tau2) {
##     p <- plogis(mu+sqrt(tau2)*zi)
##     paramsi['m']*p*(1-p)
##     }
moment.map.bn.prop <- function(zi, paramsi, mu, tau2) {
    pr <- plogis(mu+sqrt(tau2)*zi)
    out <- structure(unname(c(paramsi['m']*pr,  paramsi['m']*pr*(1-pr))), names=c('first','second'))
    return(out)
    }

cond.pmf.bn.irr <- cond.pmf.binomial.irr <- function(x,theta,mu,tau2,params,log=FALSE) {
    ## browser()
    if(is.null(dim(params))) {
        T.t <- params['T.t']; T.c <- params['T.c'];   s <- params['s']
    } else {
        T.t <- params[,'T.t']; T.c <- params[,'T.c'];   s <- params[,'s']
    }
    ## T.t <- params[,'T.t']
    ## T.c <- params[,'T.c']
    ## s <- params[,'s']
    pr <- function(theta)plogis(log(T.t/T.c) + mu+sqrt(tau2)*theta)
    dbinom(x,size=s,prob=pr(theta),log=log)
}
dgp.bn.irr <- dgp.binomial.irr <- function(n,mu,tau2,params,re.distr=rnorm) {
    T.t <- params[,'T.t']
    T.c <- params[,'T.c']
    s <- params[,'s']
    theta <- re.distr(n)
    pr <- function(theta)plogis(log(T.t/T.c) + mu+sqrt(tau2)*theta)
    rbinom(n,size=s,prob=pr(theta))
}
nn.transform.bn.irr <- function(y,params) {
    idy.cc <- y==0L | y==params[,'s']        
    y.T.nn <- y; y.C.nn <- params[,'s']-y;
    y.T.nn[idy.cc] <- y.T.nn[idy.cc]+.5; y.C.nn[idy.cc] <- y.C.nn[idy.cc]+.5
    y.nn <- log(y.T.nn/params[,'T.T']) - log(y.C.nn/params[,'T.C'])
    sigma2.nn <- 1/y.T.nn + 1/y.C.nn
    return(data.frame(y=unname(y.nn),sigma2=unname(sigma2.nn)))
}
moment.map.bn.irr <- function(zi, paramsi, mu, tau2) {
    pr <- plogis(log(paramsi['T.T'])-log(paramsi['T.C'])+mu+sqrt(tau2)*zi)
    out <- structure(unname(c(paramsi['s']*pr,  paramsi['s']*pr*(1-pr))), names=c('first','second'))
    return(out)
    }

model <- function(model.str) {
    out <- switch(model.str,
                  'bn-prop' = list(
                      cond.pmf = cond.pmf.bn.prop,
                      dgp = dgp.bn.prop,
                      moment.map = moment.map.bn.prop,
                      nn.transform = nn.transform.bn.prop,
                      pmf.derivs.fun = pmf.derivs.fun.bn.prop
                  ),
                  'pn' = list(
                      cond.pmf = cond.pmf.pn,
                      dgp = dgp.pn,
                      moment.map = moment.map.pn,
                      nn.transform = nn.transform.pn,
                      pmf.derivs.fun = pmf.derivs.fun.pn
                  ),
                  'hn' = list(
                      cond.pmf = cond.pmf.hn,
                      ## newton.step = newton.step.hn,
                      dgp = dgp.hn,
                      moment.map = moment.map.hn,
                      nn.transform = nn.transform.hn,
                      ## var.map = var.map.hn,
                      pmf.derivs.fun = pmf.derivs.fun.hn
                  ),
                  'bn-irr' = list(
                      cond.pmf = cond.pmf.bn.irr,
                      dgp = dgp.bn.irr,
                      moment.map = moment.map.bn.irr,
                      nn.transform = nn.transform.bn.irr,
                      pmf.derivs.fun = pmf.derivs.fun.bn.irr
                  ),
                  )
    return(out)
}
