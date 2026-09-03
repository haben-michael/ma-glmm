## log scale version
## my version. dotted version matches biasedUrn's call syntax.
.cond.pmf.hyperg <- function(x,m1,m2,n,odds,log=FALSE) {
    lnumer <- lchoose(m1,x) + lchoose(m2,n-x) + x*log(odds)
    supp.x <- max(0,n-m2):min(n,m1)
    ldenom.mat <-  outer(odds,supp.x,function(oddsj,xj)lchoose(m1,xj) + lchoose(m2,n-xj) + xj*log(oddsj) )
    lrowmaxes <- apply(ldenom.mat,1,max)
    ldenom  <- lrowmaxes + log(rowSums( exp(ldenom.mat - lrowmaxes)))
    out <- if(log) lnumer - ldenom else exp(lnumer-ldenom)
    return(out)
}

cond.pmf.hn <- cond.pmf.hyperg <- function(x,theta,mu,tau2,params,log=FALSE) {
    if(is.null(dim(params))) {
        m.T <- params['m.T'];m.C <- params['m.C'];   s <- params['s']
    } else {
        m.T <- params[,'m.T'];m.C <- params[,'m.C']; s <- params[,'s']
        }
    odds <- exp(mu+sqrt(tau2)*theta)
    out <- .cond.pmf.hyperg(x=x,m1=m.T,m2=m.C,n=s,odds=odds,log=log)
    return(out)
}

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
moment.map.pn <- function(zi, paramsi, mu, tau2) {
    lambda <- paramsi["T"] * exp(mu + sqrt(tau2) * zi)
    out <- structure(unname(c(lambda,lambda)), names=c('first','second'))
    return(out)
    }

cond.pmf.bn.prop <- cond.pmf.binomial.prop <- function(x,theta,mu,tau2,params,log=FALSE) {
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
moment.map.bn.prop <- function(zi, paramsi, mu, tau2) {
    ## pr <- plogis(mu+sqrt(tau2)*zi)    
    ## out <- structure(unname(c(paramsi['m']*pr,  paramsi['m']*pr*(1-pr))), names=c('first','second'))
    logp <- plogis(mu+sqrt(tau2)*zi,log.p=TRUE)
    logp.bar <- plogis(mu+sqrt(tau2)*zi,log.p=TRUE,lower.tail=FALSE)
    out <- exp(structure(unname(c(log(paramsi['m'])+logp,  log(paramsi['m'])+logp+logp.bar)), names=c('first','second')))
    return(out)
    }

cond.pmf.bn.irr <- cond.pmf.binomial.irr <- function(x,theta,mu,tau2,params,log=FALSE) {
    if(is.null(dim(params))) {
        T.t <- params['T.T']; T.c <- params['T.C'];   s <- params['s']
    } else {
        T.t <- params[,'T.T']; T.c <- params[,'T.C'];   s <- params[,'s']
    }
    pr <- function(theta)plogis(log(T.t/T.c) + mu+sqrt(tau2)*theta)
    dbinom(x,size=s,prob=pr(theta),log=log)
}
dgp.bn.irr <- dgp.binomial.irr <- function(n,mu,tau2,params,re.distr=rnorm) {
    T.t <- params[,'T.T']
    T.c <- params[,'T.C']
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
    out <- c(out, model.str=model.str)
    return(out)
}
