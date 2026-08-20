

## re-parameterize the sinh-arcsinh family

rsas <- function(n,eps,delta) {
    z <- rnorm(n)
    sinh((asinh(z)+eps)/delta)
}

dsas <- function(x,eps,delta) {
    S2 <- sinh(delta*asinh(x)-eps)^2
    delta * sqrt( (1+S2) / (2*pi*(1+x^2)) ) * exp(-S2/2)
}

sas.moment <- function(k,eps,delta) {
    integrate(function(x) x^k * dsas(x,eps,delta), -Inf, Inf)$val
}

sas.skew <- function(m) {
    var <- m[2] - m[1]^2
    skew.unstand <- m[3] - 3*m[2]*m[1] + 3*m[1]*m[1]^2 - m[1]^3
    skew.unstand / var^(3/2)
}

sas.kurt <- function(m) {
    var <- m[2] - m[1]^2
    kurt.unstand <- m[4] - 4*m[3]*m[1] + 6*m[2]*m[1]^2 - 4*m[1]*m[1]^3 + m[1]^4
    
    kurt.unstand / (var^2) - 3 
}

sas.invert.params <- function(target.skew, target.kurt) {
    obj <- function(par) {
        eps <- par[1] 
        delta <- exp(par[2])
        m <- sapply(1:4, function(k) sas.moment(k, eps, delta))
        (sas.skew(m) - target.skew)^2 + (sas.kurt(m) - target.kurt)^2
    }
    
    argmin <- optim(par=c(0,0), fn=obj)$par
    return( c(eps=argmin[1], delta=exp(argmin[2])) )
}



ks.study <- function(paramsi,mu,tau2,quad,pmf.derivs.fun,debug=FALSE)  {
    beta <- c(tau2,mu)
    ## browser()
    pmf.derivs.unit <- lapply(pmf.derivs.fun, function(pmf.deriv.fun) {
        supp.y <-  paramsi['y.min']:paramsi['y.max']
        ## units.mat <- sapply(quad.rule$z,function(zi)pmf.deriv.fun(beta=beta,params=params,theta=zi,x=supp.y))
        ## colSums(quad.rule$wt * t(units.mat))
        quad.rules <- lapply(supp.y, function(yi)quad(yi,paramsi,mu,tau2))
        nodes.mat <- sapply(quad.rules,`[[`,'wt')
        ## browser()
        units.mat <- sapply(seq_along(supp.y), function(i)pmf.deriv.fun(beta=beta,params=paramsi,theta=quad.rules[[i]]$z,x=supp.y[i]))
        colSums(nodes.mat*units.mat)
    })
    ks.unit <- sapply(seq_along(ks.fun), function(i) {
        k.fun <- ks.fun[[i]]
        args <- attr(k.fun,'args')
        k.unit <- do.call(k.fun, pmf.derivs.unit[args])
        ## k.val <- sum( k.unit )
    })
    k.vals <- colSums(ks.unit,na.rm=TRUE)#!!!!!!!!!!!!!!1 
    k.orders <- sapply(ks.fun,function(k.fun)attr(k.fun,'order'))
    ks.df <- data.frame(inner=k.orders[1,],outer=k.orders[2,],k.val=k.vals)
    ## out <- list(ks.df=ks.df, ks.unit=ks.unit, pmf.derivs.unit=pmf.derivs.unit)
    out <- list(ks.df=ks.df)
    if(debug) out <- c(out,debug=list(list(ks.unit=ks.unit, pmf.derivs.unit=pmf.derivs.unit)))
    return(out)
}
## ks.study(params[1,],mu.null,tau2.null,quad,pmf.derivs.fun)
ks.ma <- function(params,mu,tau2,quad,pmf.derivs.fun,debug=FALSE) {
    ## browser()
    ks.study.fits <- lapply(1:nrow(params), function(i) {
        ks.study.fit <- ks.study(pmf.derivs.fun=pmf.derivs.fun,paramsi=params[i,],mu=mu,tau2=tau2,quad=quad,debug=debug)
    })
    ## k.vals.ma <- sapply(ks,function(df)df$k.val) |> rowSums()
    k.vals.ma <- sapply(ks.study.fits,function(fit)fit$ks.df$k.val) |> rowSums()
    ks.df.ma <- data.frame(inner=ks.study.fits[[1]]$ks.df$inner, outer=ks.study.fits[[1]]$ks.df$outer,k.val=k.vals.ma)
    ks.study.fits <- sapply(ks.study.fits,`[[<-`,'ks.df',NULL)
    out <- list(ks.df.ma=ks.df.ma)
    if(debug) out <- c(out,debug=list(ks.study.fits))
    return(out)
    ## return(list(ks.df.ma=ks.df.ma,ks.study.fits=ks.study.fits))
}
## ks.ma(params,mu.null,tau2.null,quad,pmf.derivs.fun)
lawley <- function(k.df,p) {
    ## k.df <- data.frame(inner=sapply(orders.k,`[`,1),outer=sapply(orders.k,`[`,2),k.val=k.vals)
    ## print(k.df)
    k <- function(inner,outer=0)k.df$k.val[k.df$inner==paste(sort(inner),collapse='') & k.df$outer==paste(sort(outer),collapse='')]
    information <- -matrix(apply(expand.grid(1:p,1:p),1,k),nrow=p)
    inv.inf <- solve(information)
    lawley.4 <- apply(do.call(expand.grid,rep(list(1:p),4)), 1, function(idx.str) {
        r <- idx.str[1]; s <- idx.str[2]; t <- idx.str[3]; u <- idx.str[4]
        inv.inf[r,s] * inv.inf[t,u] * (1/4*k(c(r,s,t,u)) - k(c(r,s,t),u) + k(c(r,t),c(s,u)))
    })
    lawley.6 <- apply(do.call(expand.grid,rep(list(1:p),6)), 1, function(idx.str) {
        r <- idx.str[1]; s <- idx.str[2]; t <- idx.str[3]; u <- idx.str[4]; v <- idx.str[5]; w <- idx.str[6]
        -inv.inf[r,s] * inv.inf[t,u] * inv.inf[v,w] * (1/6*k(c(r,t,v))*k(c(s,u,w)) + 1/4*k(c(r,t,u))*k(c(s,v,w)) - k(c(r,t,v))*k(c(s,w),u) - k(c(r,t,u))*k(c(s,w),v) + k(c(r,t),v)*k(c(s,w),u) + k(c(r,t),u)*k(c(s,w),v))
    })
    out <- bias <- sum(lawley.4) - sum(lawley.6)
    return(out)
}
make.bartlett.C <- function(params,quad,pmf.derivs.fun,debug=FALSE){
    bartlett.C <- function(mu,tau2) {
        ks.ma.fit <- ks.ma(params=params,mu=mu,tau2=tau2,quad=quad,pmf.derivs.fun=pmf.derivs.fun,debug=debug)
        ks.df <- ks.ma.fit$ks.df
        bartlett.C <- lawley(ks.df,2) - lawley(ks.df,1)
        return(list(bartlett.C=bartlett.C))
    }
}
## make.quad <- function(integrand.mode,var.map,n.GH,adaptive=TRUE) {
##     GH <- fastGHQuad::gaussHermiteData(n.GH)
##     if(!adaptive) {
##         out <- function(...) {
##             list(z=sqrt(2)*GH$x,wt=GH$w/sqrt(pi))
##         }
##     } else {
##         out <- function(yi,paramsi,mu,tau2) {
##         ## browser()
##             mode <- integrand.mode(yi,paramsi,mu,tau2)
##             scale <- (1 + tau2 * var.map(mode,paramsi, mu, tau2))^(-1 / 2)
##             z <- mode + sqrt(2) * scale * GH$x
##             wt <- sqrt(2) * scale *
##                 GH$w *
##                 exp(GH$x^2) *
##                 dnorm(z)
##             return(list(z=z,wt=wt))
##         }
##     }
##     return(out)
## }
## using moment map
make.quad <- function(moment.map,n.GH,adaptive=TRUE,adaptive.tol=NULL) {
    GH <- fastGHQuad::gaussHermiteData(n.GH)
    if(!adaptive) {
        out <- function(...) {
            list(z=sqrt(2)*GH$x,wt=GH$w/sqrt(pi),log.wt=log(GH$w/sqrt(pi)))
        }
    } else {
        newton.step <- function(zi, yi, paramsi, mu, tau2) {
            moments <- moment.map(zi,paramsi,mu,tau2)
            deriv <- sqrt(tau2)*(yi-moments['first'])-zi
            hess  <- -tau2*moments['second']-1
            out <-  - deriv / hess    
            return(out)
        }
        integrand.mode <- function(yi, paramsi, mu, tau2) {
            zi <- 0 #!!!!!!!!
            ## browser()
            for (j in 1:20) { #!!!!
                zi.new <- zi + newton.step(zi, yi, paramsi, mu, tau2)
                ## if(!is.finite(zi.new))browser()
                ## 1+1
                ## print(zi.new)
                ## print(c(zi.new,zi))
                if (abs(zi.new - zi) < adaptive.tol) {
                    break
                }
                zi <- zi.new
            }
            mode <- zi.new
            return(mode)
        }
        out <- function(yi,paramsi,mu,tau2) {
            ## browser()
            mode <- integrand.mode(yi,paramsi,mu,tau2)
            ## scale <- (1 + tau2 * var.map(mode,paramsi, mu, tau2))^(-1 / 2)
            var.star <- moment.map(mode,paramsi,mu,tau2)['second']
            scale <- (1 + tau2 * var.star)^(-1 / 2)
            z <- mode + sqrt(2) * scale * GH$x
            ## wt <- sqrt(2) * scale * 
            ##     GH$w *
            ##     exp(GH$x^2) *
            ##     dnorm(z)
            wt <- 1/sqrt(pi) * scale * 
                GH$w *
                exp(GH$x^2 - z^2/2)                       
            return(list(z=z,wt=wt,log.wt=log(wt)))
        }
    }
    return(out)
}
make.loglik <- function(cond.pmf,params,y,quad) {
    loglik.unit <- function(mu,tau2,yi,paramsi) {
        rule <- quad(yi,paramsi,mu,tau2)
        ## likelihood <- sum(rule$wt*cond.pmf(yi,rule$z,mu,tau2,paramsi))
        ## out <- log(likelihood)
        lsummands <- rule$log.wt + cond.pmf(yi,rule$z,mu,tau2,paramsi,log=TRUE)
        max.lsummands <- max(lsummands)
        out <- max.lsummands + log(sum(exp(lsummands - max.lsummands)))
        ## if(tau2>1e2)browser()
        ## print(c(mu,tau2,out))
        return(out)
    }
    loglik <- function(mu,tau2) {
        ## browser()
        ## print(c(mu,tau2))
        out <- sum(sapply(seq_along(y),function(i)loglik.unit(mu,tau2,y[i],params[i,])))
        ## if(!is.finite(out))browser()
        ## print(c(mu,tau2,out))
        ## ## if(is.nan(out))browser()
        ## if(is.nan(out))out <- -1e10#-.Machine$double.xmax
        return(out)
    }
    return(loglik)
}
make.profile.LR <- function(loglik,mle,par.start,logtau2.bounds) {
    LR.profile <- function(mu.null) {
        logtau2.mle.null <- optimize(function(logtau2)loglik(mu.null,exp(logtau2)),interval=logtau2.bounds,maximum=TRUE)$max
        tau2.mle.null <- exp(logtau2.mle.null)
        LR <- 2 * (loglik(mle['mu'],mle['tau2']) - loglik(mu.null,tau2.mle.null))
        out <- list(mle.null=c(mu=unname(mu.null),tau2=unname(tau2.mle.null)), LR=LR)
        return(out)
    }
    return(LR.profile)
}
## make.integrand.mode <- function(newton.step, tol=1e-10) {
##     integrand.mode <- function(yi, paramsi, mu, tau2) {
##         zi <- 0 #!!!!!!!!
##         ## browser()
##         for (j in 1:20) { #!!!!
##             zi.new <- zi + newton.step(zi, yi, paramsi, mu, tau2)
##             if(!is.finite(zi.new))browser()
##             ## 1+1
##             ## print(zi.new)
##             ## print(c(zi.new,zi))
##             if (abs(zi.new - zi) < tol) {
##                 break
##             }
##             zi <- zi.new
##         }
##         mode <- zi.new
##         return(mode)
##     }
##     return(integrand.mode)
## }
## z.deriv <- function(zi, yi, paramsi, mu, tau2) {
##     ## browser()
##     p <- plogis(mu+sqrt(tau2)*zi)
##     sqrt(tau2)*( yi-paramsi['m']*p )
##     }
## curve(cond.pmf.hyperg(yi,x,mu,tau2,paramsi)*dnorm(x),-5,5)
mle <- function(loglik,par.start,logtau2.bounds=c(-20,3),...) {
    ## browser()
    ## print(par.start)
    fit <- optim(par=c(par.start['mu'],log(par.start['tau2'])),fn=function(par)-loglik(par[1],exp(par[2])),...)
    ## fit <- optim(
    ##     par = c(par.start['mu'], log(par.start['tau2'])),
    ##     fn = function(par) -loglik(par[1], exp(par[2])),
    ##     method = 'L-BFGS-B',
    ##     lower = c(-Inf, logtau2.bounds[1]),
    ##     upper = c( Inf, logtau2.bounds[2]),
    ##     ...
    ## )
    mu.mle <- unname(fit$par[1]); tau2.mle <- unname(exp(fit$par[2]))
    out <- list(mle=c(mu=mu.mle,tau2=tau2.mle), status=fit$convergence)
}
profile <- function(mu.grid,profile.LR,mle,bartlett.C=NULL,mc.cores=1) {
    ## browser()
    profile.fits <- lapply(mu.grid, function(mu.null)profile.LR(mu.null))
    LR.profile <- sapply(profile.fits,`[[`,'LR')
    tau2.grid <- sapply(profile.fits,`[[`,'mle.null')['tau2',]
    out <- list(LR=LR.profile,mu.grid=mu.grid,tau2.grid=tau2.grid)
    if(!is.null(bartlett.C)) {
        bartlett.C.fits <- parallel::mclapply(seq_along(mu.grid),mc.cores=mc.cores,FUN=function(i)    bartlett.C(mu.grid[i],tau2.grid[i]))
        bartlett.C.profile <- sapply(bartlett.C.fits,`[[`,'bartlett.C')
        LR.bartlett.profile <- LR.profile / (1+bartlett.C.profile)
        out <- c(out, list(LR.bartlett=LR.bartlett.profile,bartlett.C=bartlett.C.profile))
        ## if(debug) {
        ##     debug.profile <- lapply(bartlett.C.fits,`[[`,'debug')
        ##     out <- c(out,debug=list(debug.profile))
        ##     ## browser()
        ## }
    }
    return(out)
}
confint <- function(mu.grid, profile, alpha=.05, plot=TRUE,...) {
    q <- qchisq(1-alpha,df=1)
    ## browser()
    if(plot) {
        ## plot(mu.grid,profile['LR',],ylim=range(profile),type='o')
        ## if(!is.null(LR.bartlett.profile)) lines(mu.grid,profile['LR.bartlett',],lty=2,type='o')
        matplot(mu.grid,t(profile),type='l',...)
        abline(h=q)
    }
    ## interpolate CIs
    CIs <- apply(profile - q, 1, simplify=FALSE, FUN=function(vals) {
        left.idx <- which(diff(sign(vals))!=0)
        stopifnot(length(left.idx)>0)
        right.idx <- left.idx+1 
        left.wt <- abs(vals[right.idx]/(vals[right.idx]-vals[left.idx]))
        right.wt <- abs(vals[left.idx]/(vals[right.idx]-vals[left.idx]))
        zeroes <- left.wt*mu.grid[left.idx] + right.wt*mu.grid[right.idx]
        if(vals[1]<0)zeroes <- c(-Inf,zeroes)
        if(vals[length(vals)]<0)zeroes <- c(zeroes,Inf)
        n.zeroes <- length(zeroes)
        CI.left <- zeroes[(1:n.zeroes) %% 2 ==1]
        CI.right <- zeroes[(1:n.zeroes) %% 2 ==0]
        cbind(CI.left,CI.right)
    } )
    out <- list(CIs=CIs)
    ## if(debug) out <- c(out,profile.fit$debug)
    return(out)
}
par.start <- function(y,params,nn.transform) {
    ## browser()
    ## 1+1
    with(nn.transform(y,params), {
        w.fe <- 1/sigma2
        mu.fe <- y %*% w.fe / sum(w.fe)    
        Q <- w.fe%*%(y-c(mu.fe))^2
        tau2.DL <- max(0L, (Q - (length(y)-1)) / (sum(w.fe)-sum(w.fe^2)/sum(w.fe)))
        tau2.start <- if(tau2.DL==0L) var(y) else tau2.DL #!!!!!
        c(mu=mu.fe,tau2=tau2.start)
    })
}
meta.glmm <- function(y,params,model.str,n.GH,adaptive,logtau2.bounds,...) {
    model <- model(model.str)
    ## integrand.mode <- make.integrand.mode(model$newton.step)
    quad <- make.quad(model$moment.map,n.GH=n.GH,adaptive=adaptive,adaptive.tol=1e-10)#!!!!!!
    loglik <- make.loglik(model$cond.pmf,params,y,quad)
    ## par.start0 <- model$par.start(y,params)
    par.start0 <- par.start(y,params,model$nn.transform)
    ## browser()
    mle0 <- mle(loglik,par.start0,logtau2.bounds,...)$mle
    profile.LR <- make.profile.LR(loglik,mle0,par.start0,logtau2.bounds=logtau2.bounds)
    bartlett.C <- make.bartlett.C(params,quad,model$pmf.derivs.fun)
    out <- list(model.str=model.str, quad=quad, loglik=loglik, mle=mle0, profile.LR=profile.LR, bartlett.C=bartlett.C, par.start=par.start0)
    return(out)
}
## LR.profile.meta.glmm <- function(y,params,modelmu.null
confint.meta.glmm <- function(meta.glmm.fit,alpha=.05,plot=TRUE,debug=FALSE,bartlett=TRUE,mu.grid=list(offset=1.5,size=20),mc.cores=1,...) {
    ## browser()
    mle0 <- meta.glmm.fit$mle
    mu.grid <- seq(mle0['mu']-mu.grid$offset*sqrt(mle0['tau2']),mle0['mu']+mu.grid$offset*sqrt(mle0['tau2']),len=mu.grid$size)    
    profile0 <- profile(mu.grid,meta.glmm.fit$profile.LR,mle0,if(bartlett)meta.glmm.fit$bartlett.C else NULL)
    confint(mu.grid=mu.grid, profile=rbind(LR=profile0$LR,LR.bartlett=profile0$LR.bartlett), alpha=alpha, plot=plot,...)
}



## ## .cond.pmf.hyperg <- function(x,theta,mu,tau2,params) {
## ##     ## x <- 0; theta <-  -18.95964; mu <- 0; tau2 <- 0.3071743;params <- c(m.T=116,m.C=117,s=3,x.min=0,x.max=3)
## ##     m.T <- params[,'m.T']
## ##     m.C <- params[,'m.C']
## ##     s <- params[,'s']
## ##     odds <- exp(mu+sqrt(tau2)*theta)
## ##     ## print(x);print(s);print(params)
## ##     ## browser()
## ##     ## print(c(m.T,m.C,s,odds))
## ##         BiasedUrn::dFNCHypergeo(x=x,m1=m.T,m2=m.C,n=s,odds=odds)            
## ## }
## ## cond.pmf.hyperg <-  Vectorize(.cond.pmf.hyperg,vectorize.args=c("theta"))

## ## cond.pmf.hyperg <- function(x,theta,mu,tau2,params)  <- sapply(theta, function(thetai) _cond.pmf.hyperg(x,thetai,mu,tau2,params)


## ## my version. dotted version matches biasedUrn's call syntax.
## .cond.pmf.hyperg <- function(x,m1,m2,n,odds) {
##     ## m.T <- params[,'m.T']
##     ## m.C <- params[,'m.C']
##     ## s <- params[,'s']
##     numer <- choose(m1,x)*choose(m2,n-x)*odds^x
##     supp.x <- max(0,n-m2):min(n,m1)
##     denom <- rowSums(  outer(odds,supp.x,function(oddsj,xj)choose(m1,xj)*choose(m2,n-xj)*oddsj^xj ) )
##     numer / denom
## }
## cond.pmf.hyperg <- function(x,theta,mu,tau2,params) {
##     m.T <- params[,'m.T']
##     m.C <- params[,'m.C']
##     s <- params[,'s']
##     odds <- exp(mu+sqrt(tau2)*theta)
##     out <- .cond.pmf.hyperg(x=x,m1=m.T,m2=m.C,n=s,odds=odds)
##     ## print(out)
##     return(out)
## }

## ## ## safer but slower:
## ## cond.pmf.hyperg <- function(x,theta,mu,tau2,params) {
## ##     m.T <- params[,'m.T']; m.C <- params[,'m.C']; s <- params[,'s']
## ##     logpmf.part <- function(xi)lgamma(m.T+1) - lgamma(xi+1) -lgamma(m.T-xi+1) + lgamma(m.C+1) - lgamma(s-xi+1) - lgamma(m.C-s+xi+1)
## ##     lognumer <- logpmf.part(x)+ (mu+sqrt(tau2)*theta)*x
## ##     ## logpmf.mat <- logpmf.part(0:m.T) + outer(0:m.T,theta)*sqrt(tau2) + (0:m.T)*mu
## ##     supp <- max(0,s-m.C):min(s,m.T)
## ##     logpmf.mat <- logpmf.part(supp) + outer(supp,theta)*sqrt(tau2) + supp*mu
## ##     colmaxes <- apply(logpmf.mat,2,max)        
## ##     logdenom <- colmaxes + log(rowSums(exp(t(logpmf.mat) - colmaxes)))
## ##     exp(lognumer - logdenom)
## ## }

## ## n <- sample(2:10,1)
## ## x <- sample(1:n,1)
## ## m1 <- 8; m2 <- 10
## ## odds <- rchisq(1,df=5)
## ## BiasedUrn::dFNCHypergeo(x=x,m1=m1,m2=m2,n=n,odds=odds)
## ## .cond.pmf.hyperg(x,m1,m2,n,odds)


## ## cond.pmf.hyperg <- function(x,theta,mu,tau2,params) {
## ##     m.T <- params[,'m.T']
## ##     m.C <- params[,'m.C']
## ##     s <- params[,'s']
## ##     numer <- choose(m.T,x)*choose(m.C,s-x)*exp((mu+sqrt(tau2)*theta)*x)
## ##     supp <- max(0,s-m.C):min(s,m.T)
## ##     denom <- rowSums(  outer(theta,supp,function(thetaj,xj)choose(m.T,xj)*choose(m.C,s-xj)*exp((mu+sqrt(tau2)*thetaj)*xj) ) )
## ##     numer / denom
## ## }
## dgp.hyperg <- function(n,mu,tau2,params,re.distr=rnorm) {
##     m.T <- params[,'m.T']
##     m.C <- params[,'m.C']
##     s <- params[,'s']
##     theta <- re.distr(n)
##     odds <- exp(mu+sqrt(tau2)*theta)
##     x <- x.T <- sapply(1:n, function(i) BiasedUrn::rFNCHypergeo(1,m1=m.T[i],m2=m.C[i],n=s[i],odds=odds[i]))
## }


## par.start.hyperg <- function(x,params) {
##     cells <- data.frame(y.1.T=x, y.1.C =params[,'s']-x)
##     cells <- within(cells, {
##         y.0.T <- params[,'m.T']-y.1.T;
##         y.0.C <- params[,'m.C']-y.1.C
##     })
##     ## y.1.T <- x; y.0.T <- params[,'m.T']-y.1.T
##     ## y.1.C <- params[,'s']-y.1.T; y.0.C <- params[,'m.C']-y.1.C
##     ## idx.cc <- y.1.T==0L | y.0.T==0L | y.1.C==0L | y.0.C==0L
##     idx.cc <- apply(cells,1,function(v)any(v==0L))
##     cells[idx.cc,] <- cells[idx.cc,]+.5
##     y.nn <- with(cells, log(y.1.T/y.0.T) - log(y.1.C/y.0.C))
##     sigma2.nn <- with(cells, 1/y.1.T + 1/y.0.T + 1/y.1.C + 1/y.0.C)
##     w.fe <- 1/sigma2.nn
##     mu.fe <- y.nn %*% w.fe / sum(w.fe)    
##     Q <- w.fe%*%(y.nn-c(mu.fe))^2
##     tau2.DL <- max(0L, (Q - (length(y.nn)-1)) / (sum(w.fe)-sum(w.fe^2)/sum(w.fe)))
##     tau2.start <- if(tau2.DL>0L)tau2.DL else var(y.nn)/1e4
##     return(c(mu=mu.fe,tau2=tau2.start))
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

## cond.pmf.poisson <- function(x,theta,mu,tau2,params) {
##     T <- params[,'T']
##     lambda <- T*exp(mu+sqrt(tau2)*theta)
##     dpois(x,lambda=lambda)
## }
## dgp.poisson <- function(n,mu,tau2,params,re.distr=rnorm) {
##     T <- params[,'T']
##     theta <- re.distr(n)
##     lambda <- T*exp(mu+sqrt(tau2)*theta)
##     rpois(n,lambda)
##     ## pr <- function(theta)plogis( mu+sqrt(tau2)*theta)
##     ## rbinom(n,size=m,prob=pr(theta))
## }



## cond.pmf.binomial.prop <- function(x,theta,mu,tau2,params) {
##     m <- params[,'m']
##     pr <- function(theta)plogis(mu+sqrt(tau2)*theta)
##     dbinom(x,size=m,prob=pr(theta))
##     ## integrand <- function(theta)pr(theta)^xi*(1-pr(theta))^(mi-xi)
## }
## dgp.binomial.prop <- function(n,mu,tau2,params,re.distr=rnorm) {
##     m <- params[,'m']
##     theta <- re.distr(n)
##     pr <- function(theta)plogis( mu+sqrt(tau2)*theta)
##     rbinom(n,size=m,prob=pr(theta))
## }

## par.start.binomial.prop <- function(x,params) {
##     x.1 <- x
##     x.0 <- params[,'m'] - x
##     idx.cc <- x.0==0 | x.1==0
##     x.0[idx.cc] <- x.0[idx.cc]+.5
##     x.1[idx.cc] <- x.1[idx.cc]+.5
##     y.nn <- log(x.1/x.0)
##     sigma2.nn <- 1/x.0 + 1/x.1
##     w.fe <- 1/sigma2.nn
##     mu.fe <- y.nn %*% w.fe / sum(w.fe)    
##     Q <- w.fe%*%(y.nn-c(mu.fe))^2
##     tau2.DL <- max(0L, (Q - (length(y.nn)-1)) / (sum(w.fe)-sum(w.fe^2)/sum(w.fe)))
##     tau2.start <- if(tau2.DL>0L)tau2.DL else var(y.nn)/1e4
##     return(c(mu=mu.fe,tau2=tau2.start))
## }


## par.start.binomial.irr <- function(x,params) {
##     idx.cc <- x==0L | x==params[,'s']        
##     x.T.nn <- x; x.C.nn <- params[,'s']-x;
##     x.T.nn[idx.cc] <- x.T.nn[idx.cc]+.5; x.C.nn[idx.cc] <- x.C.nn[idx.cc]+.5
##     y.nn <- log(x.T.nn/params[,'T.t']) - log(x.C.nn/params[,'T.c'])
##     sigma2.nn <- 1/x.T.nn + 1/x.C.nn
##     w.fe <- 1/sigma2.nn
##     mu.fe <- y.nn %*% w.fe / sum(w.fe)    
##     Q <- w.fe%*%(y.nn-c(mu.fe))^2
##     tau2.DL <- max(0L, (Q - (length(y.nn)-1)) / (sum(w.fe)-sum(w.fe^2)/sum(w.fe)))
##     tau2.start <- if(tau2.DL>0L)tau2.DL else max(1,var(y.nn))/1e4
##     return(c(mu=mu.fe,tau2=tau2.start))
## }
## cond.pmf.binomial.irr <- function(x,theta,mu,tau2,params) {
##     T.t <- params[,'T.t']
##     T.c <- params[,'T.c']
##     s <- params[,'s']
##     pr <- function(theta)plogis(log(T.t/T.c) + mu+sqrt(tau2)*theta)
##     dbinom(x,size=s,prob=pr(theta))
## }
## dgp.binomial.irr <- function(n,mu,tau2,params,re.distr=rnorm) {
##     T.t <- params[,'T.t']
##     T.c <- params[,'T.c']
##     s <- params[,'s']
##     theta <- re.distr(n)
##     pr <- function(theta)plogis(log(T.t/T.c) + mu+sqrt(tau2)*theta)
##     rbinom(n,size=s,prob=pr(theta))
## }




## ## ##  rework ks.fun to store order and args as attributes not list
## ## ## elements.
## ## ks.study <- function(pmf.derivs.fun,ks.fun,params,supp.x,mu.null,tau2.null,n.GH=50) {
## ##     ## browser()
## ##     beta <- c(tau2.null,mu.null)
## ##     ## if(is.null(GH)) {
## ##     GH <- fastGHQuad::gaussHermiteData(n.GH)
## ##     GH$x <- sqrt(2)*GH$x
## ##     GH$w <- GH$w / sqrt(pi)
## ##     names(GH) <- c('theta.grid','w')
## ##     ## }
## ##     pmf.derivs.unit <- lapply(pmf.derivs.fun, function(pmf.deriv.fun) {
## ##     ## browser()
## ##         ## supp <- max(0,s-m.C):min(s,m.T)
## ##         ## g.mat <- outer(X=GH$theta.grid,Y=supp,function(theta,x)pmf.deriv.fun.ma(theta,x,m,m.C,s)) # these are globals!!!
## ##         ## g.mat <- outer(X=GH$theta.grid,Y=supp.x,function(theta,x)pmf.deriv.fun(beta=beta,theta=theta,x=x,params=params)) 
## ##         ## g.unit <- colSums(GH$w * g.mat)
## ##         units.mat <- sapply(GH$theta,function(thetai)pmf.deriv.fun(beta=beta,params=params,theta=thetai,x=supp.x))
## ##         colSums(GH$w * t(units.mat))
## ##     })
## ##     k.vals <- sapply(seq_along(ks.fun), function(i) {
## ##         k.fun <- ks.fun[[i]]
## ##         args <- attr(k.fun,'args')
## ##         k.fun.unit <- do.call(k.fun, pmf.derivs.unit[args])
## ##         k.val <- sum( k.fun.unit )
## ##     })
## ##     k.orders <- sapply(ks.fun,function(k.fun)attr(k.fun,'order'))
## ##     ## k.df <- data.frame(inner=sapply(k.orders,`[`,1),outer=sapply(k.orders,`[`,2),k.val=k.vals)
## ##     ks.df <- data.frame(inner=k.orders[1,],outer=k.orders[2,],k.val=k.vals)
## ##     return(ks.df)
## ## }


## ks.study <- function(pmf.derivs.fun,ks.fun,params,mu,tau2,n.GH=1e2,debug=FALSE) {
##     ## browser()
##     beta <- c(tau2,mu)
##     ## if(is.null(GH)) {
##     GH <- fastGHQuad::gaussHermiteData(n.GH)
##     GH$x <- sqrt(2)*GH$x
##     GH$w <- GH$w / sqrt(pi)
##     names(GH) <- c('theta.grid','w')
##     ## }
##     pmf.derivs.unit <- lapply(pmf.derivs.fun, function(pmf.deriv.fun) {
##         ## browser()
##         supp.x <-  params['x.min']:params['x.max']
##         units.mat <- sapply(GH$theta.grid,function(thetai)pmf.deriv.fun(beta=beta,params=params,theta=thetai,x=supp.x))
##         colSums(GH$w * t(units.mat))
##     })
##     ## k.vals <- sapply(seq_along(ks.fun), function(i) {
##     ##     k.fun <- ks.fun[[i]]
##     ##     args <- attr(k.fun,'args')
##     ##     k.unit <- do.call(k.fun, pmf.derivs.unit[args])
##     ##     k.val <- sum( k.unit )
##     ## })
##     ks.unit <- sapply(seq_along(ks.fun), function(i) {
##         k.fun <- ks.fun[[i]]
##         args <- attr(k.fun,'args')
##         k.unit <- do.call(k.fun, pmf.derivs.unit[args])
##         ## k.val <- sum( k.unit )
##     })
##     k.vals <- colSums(ks.unit) 
##     k.orders <- sapply(ks.fun,function(k.fun)attr(k.fun,'order'))
##     ks.df <- data.frame(inner=k.orders[1,],outer=k.orders[2,],k.val=k.vals)
##     ## out <- list(ks.df=ks.df, ks.unit=ks.unit, pmf.derivs.unit=pmf.derivs.unit)
##     out <- list(ks.df=ks.df)
##     if(debug) out <- c(out,debug=list(list(ks.unit=ks.unit, pmf.derivs.unit=pmf.derivs.unit)))
##     return(out)
## }


## ## ks.ma <- function(params,ks.fun,pmf.derivs.fun,mu,tau2,n.GH=50) {
## ##     ## browser()
## ##     ks <- lapply(1:nrow(params), function(i) {
## ##         ## if(i==16)browser()
## ##         supp.x <- params[i,'x.min']:params[i,'x.max']
## ##         ks.df.study <- ks.study(pmf.derivs.fun=pmf.derivs.fun,ks.fun=ks.fun,params=params[i,],mu,tau2,n.GH=n.GH)$ks.df
## ##     })
## ##     k.vals.ma <- sapply(ks,function(df)df$k.val) |> rowSums()
## ##     ks.df.ma <- data.frame(inner=ks[[1]]$inner, outer=ks[[1]]$outer,k.val=k.vals.ma)
## ##     return(ks.df.ma)
## ## }

## ## with debugging info
## ks.ma <- function(params,ks.fun,pmf.derivs.fun,mu,tau2,n.GH=1e2,debug=FALSE) {
##     ## browser()
##     ks.study.fits <- lapply(1:nrow(params), function(i) {
##         ## if(i==16)browser()
##         ## supp.x <- params[i,'x.min']:params[i,'x.max']
##         ks.study.fit <- ks.study(pmf.derivs.fun=pmf.derivs.fun,ks.fun=ks.fun,params=params[i,],mu,tau2,n.GH=n.GH,debug=debug)
##     })
##     ## k.vals.ma <- sapply(ks,function(df)df$k.val) |> rowSums()
##     k.vals.ma <- sapply(ks.study.fits,function(fit)fit$ks.df$k.val) |> rowSums()
##     ks.df.ma <- data.frame(inner=ks.study.fits[[1]]$ks.df$inner, outer=ks.study.fits[[1]]$ks.df$outer,k.val=k.vals.ma)
##     ks.study.fits <- sapply(ks.study.fits,`[[<-`,'ks.df',NULL)
##     out <- list(ks.df.ma=ks.df.ma)
##     if(debug) out <- c(out,debug=list(ks.study.fits))
##     return(out)
##     ## return(list(ks.df.ma=ks.df.ma,ks.study.fits=ks.study.fits))
## }


## lawley <- function(k.df,p) {
##     ## k.df <- data.frame(inner=sapply(orders.k,`[`,1),outer=sapply(orders.k,`[`,2),k.val=k.vals)
##     ## print(k.df)
##     k <- function(inner,outer=0)k.df$k.val[k.df$inner==paste(sort(inner),collapse='') & k.df$outer==paste(sort(outer),collapse='')]
##     information <- -matrix(apply(expand.grid(1:p,1:p),1,k),nrow=p)
##     inv.inf <- solve(information)
##     lawley.4 <- apply(do.call(expand.grid,rep(list(1:p),4)), 1, function(idx.str) {
##         r <- idx.str[1]; s <- idx.str[2]; t <- idx.str[3]; u <- idx.str[4]
##         inv.inf[r,s] * inv.inf[t,u] * (1/4*k(c(r,s,t,u)) - k(c(r,s,t),u) + k(c(r,t),c(s,u)))
##     })
##     lawley.6 <- apply(do.call(expand.grid,rep(list(1:p),6)), 1, function(idx.str) {
##         r <- idx.str[1]; s <- idx.str[2]; t <- idx.str[3]; u <- idx.str[4]; v <- idx.str[5]; w <- idx.str[6]
##         -inv.inf[r,s] * inv.inf[t,u] * inv.inf[v,w] * (1/6*k(c(r,t,v))*k(c(s,u,w)) + 1/4*k(c(r,t,u))*k(c(s,v,w)) - k(c(r,t,v))*k(c(s,w),u) - k(c(r,t,u))*k(c(s,w),v) + k(c(r,t),v)*k(c(s,w),u) + k(c(r,t),u)*k(c(s,w),v))
##     })
##     out <- bias <- sum(lawley.4) - sum(lawley.6)
##     return(out)
## }

## ## bartlett.LR <- function(params,ks.fun,pmf.derivs.fun,mu,tau2,n.GH=50,debug=FALSE) {
## ##     ks.ma.fit <- ks.ma(params=params,ks.fun=ks.fun,pmf.derivs.fun=pmf.derivs.fun,mu=mu,tau2=tau2,n.GH=n.GH)
## ##     ks.df.ma <- ks.ma.fit$ks.df.ma
## ##     bartlett.C <- lawley(ks.df.ma,2) - lawley(ks.df.ma,1)
## ##     return(bartlett.C)
## ## }



## ## pass in starting values for optim, stop starting with tau2.mle.redu
## LR.gaussian <- function(x,params,mu.null,cond.pmf,logtau2.bounds=c(-3,3),par.start=c(mu=0,tau2=1),n.GH=1e2,...) {
##     ## browser()
##     GH <- fastGHQuad::gaussHermiteData(n.GH)
##     GH$x <- sqrt(2)*GH$x; GH$w <- GH$w/sqrt(pi)
##     loglik.unit <- function(mu,tau2,xi,paramsi,GH) {
##         ## browser()
##         ## print(c(mu,tau2,xi,paramsi))
##         ## if(tau2>1e3)browser()
##         if(log(tau2)<logtau2.bounds[1] || log(tau2)>logtau2.bounds[2]) {
##             ## warning('tau2 boundary  hit, returning minimum negative')
##             return(-1e10)#     return(-.Machine$double.xmax)
##         }
##         likelihood <- sum(GH$w*cond.pmf(xi,GH$x,mu,tau2,paramsi))
##         ## print(c(mu,tau2,log(likelihood),paramsi))
##         ## print(tau2)
##         return(log(likelihood))
##     }
##     loglik <- function(mu,tau2,x,params,GH) {
##         out <- sum(sapply(seq_along(x),function(i)loglik.unit(mu,tau2,x[i],params[i,,drop=FALSE],GH)))
##         ## print(out)        
##         ## if(is.nan(out))browser()
##         if(is.nan(out))out <- -1e10#-.Machine$double.xmax
##         return(out)
##         }
##     logtau2.mle.redu <- optimize(function(logtau2)loglik(mu.null,exp(logtau2),x,params,GH),interval=logtau2.bounds,maximum=TRUE)$max
##     tau2.mle.redu <- exp(logtau2.mle.redu)
##     start.val <- loglik(par.start[1],par.start[2],x,params,GH)
##     if(is.nan(start.val) || !is.finite(start.val))browser()
##     fit.full <- optim(par=c(par.start['mu'],log(par.start['tau2'])),fn=function(par)-loglik(par[1],exp(par[2]),x,params,GH),method='BFGS',...)
##     if(fit.full$convergence != 0L)warning('full fit convergence failure')
##     mle.full <- fit.full$par
##     mu.mle.full <- mle.full[1]; tau2.mle.full <- exp(mle.full[2])
##     LR <- 2 * (loglik(mu.mle.full,tau2.mle.full,x,params,GH) - loglik(mu.null,tau2.mle.redu,x,params,GH))
##     return(list(mle.full=c(mu=unname(mu.mle.full),tau2=unname(tau2.mle.full)),mle.redu=c(mu=unname(mu.null),tau2=unname(tau2.mle.redu)), LR=LR, loglik.fun=loglik, x=x, params=params, cond.pmf=cond.pmf,logtau2.bounds=logtau2.bounds, n.GH=n.GH))
## }

## profile.LR <- function(LR.fit,grid=list(offset=1,size=20),logtau2.bounds=LR.fit$logtau2.bounds,n.GH=LR.fit$n.GH,mc.cores=1) {
##     ## browser()
##     GH <- fastGHQuad::gaussHermiteData(n.GH)
##     GH$x <- sqrt(2)*GH$x; GH$w <- GH$w/sqrt(pi)
##     ## names(GH)[which(names(GH)=='x')] <- 'theta.grid'
##     x <- LR.fit$x
##     params <- LR.fit$params
##     cond.pmf <- LR.fit$cond.pmf
##     mu.mle.full <- LR.fit$mle.full['mu']; tau2.mle.full <- LR.fit$mle.full['tau2'];
##     loglik <- LR.fit$loglik.fun
##     loglik.profile <- function(mu) {
##         logtau2.mle.profile <- optimize(function(logtau2)loglik(mu,exp(logtau2),x,params,GH),interval=logtau2.bounds,maximum=TRUE)$max
##         tau2.mle.profile <- exp(logtau2.mle.profile)
##         return(tau2.mle.profile)
##     }
##     mu.profile.grid <- seq(mu.mle.full-grid$offset*sqrt(tau2.mle.full),mu.mle.full+grid$offset*sqrt(tau2.mle.full),len=grid$size)
##     tau2.profile.grid <- sapply(mu.profile.grid,function(mu)loglik.profile(mu))
##     ## add support for parallelization here
##     ## profiled.stats <- lapply(seq_along(mu.profile.grid),function(i) {
##     ##     cat('.')
##     ##     LR.profile <- 2 * (loglik(mu.mle.full,tau2.mle.full,x,params,GH) - loglik(mu.profile.grid[i],tau2.profile.grid[i],x,params,GH))
##     ##     ## ks.profile <- ks.ma(params,ks.fun,pmf.derivs.fun,mu=mu.profile.grid[i],tau2=tau2.profile.grid[i],n.GH=n.GH)
##     ##     ks.fit.profile <- ks.ma(params,ks.fun,pmf.derivs.fun,mu=mu.profile.grid[i],tau2=tau2.profile.grid[i],n.GH=n.GH)
##     ##     ks.profile <- ks.fit.profile$ks.df.ma
##     ##     bartlett.C.profile <- lawley(ks.profile,2) - lawley(ks.profile,1)
##     ##     return(list(LR=LR.profile,LR.bartlett=LR.profile / (1+bartlett.C.profile),bartlett.C=bartlett.C.profile,ks=ks.profile,ks.fit=ks.fit.profile))
##     ## })
##     ## ks.profile <- lapply(profiled.stats,`[[`,'ks')
##     ## profiled.stats <- lapply(profiled.stats,`[[<-`,'ks',NULL) |> sapply(unlist)
##     ## return(list(profile=profiled.stats,mu.grid=mu.profile.grid,tau2.grid=tau2.profile.grid,ks.profile=ks.profile))
##     ## profiled.stats <- lapply(seq_along(mu.profile.grid),function(i) {
##     profiled.stats <- parallel::mclapply(seq_along(mu.profile.grid),mc.cores=mc.cores,FUN=function(i) {
##         cat('.')
##         LR.profile <- 2 * (loglik(mu.mle.full,tau2.mle.full,x,params,GH) - loglik(mu.profile.grid[i],tau2.profile.grid[i],x,params,GH))
##         ## ks.profile <- ks.ma(params,ks.fun,pmf.derivs.fun,mu=mu.profile.grid[i],tau2=tau2.profile.grid[i],n.GH=n.GH)
##         ks.fit.profile <- ks.ma(params,ks.fun,pmf.derivs.fun,mu=mu.profile.grid[i],tau2=tau2.profile.grid[i],n.GH=n.GH)
##         ks.df.ma.profile <- ks.fit.profile$ks.df.ma
##         ## browser()
##         bartlett.C.profile <- lawley(ks.df.ma.profile,2) - lawley(ks.df.ma.profile,1)
##         return(list(LR=LR.profile,LR.bartlett=LR.profile / (1+bartlett.C.profile),bartlett.C=bartlett.C.profile,ks.fit=ks.fit.profile))
##     })
##     ks.fit.profile <- lapply(profiled.stats,function(lst)lst$ks.fit)
##     ks.df.ma.profile <- lapply(ks.fit.profile,function(lst)lst$ks.df.ma)
##     profiled.stats <- lapply(profiled.stats,`[[<-`,'ks.fit',NULL) |> sapply(unlist)
##     return(list(profile=profiled.stats,mu.grid=mu.profile.grid,tau2.grid=tau2.profile.grid,ks.fit.profile=ks.fit.profile,ks.df.ma=ks.df.ma.profile))
## }

## ## confint.LR <- function(fit,logtau2.bounds=c(-5,5),grid=list(offset=1,size=20),alpha=.05,n.GH=fit$n.GH) {
## confint.LR <- function(profile.fit,alpha=.05,plot=TRUE,...) {
##     ## browser()
##     ## profiled.stats <- profile.fit$profile
##     ## LR.profile <- sapply(profiled.stats,function(x)x$LR)
##     ## LR.bartlett.profile <- sapply(profiled.stats,function(x)x$LR.bartlett)
##     profile <- profile.fit$profile
##     mu.grid <- profile.fit$mu.grid
##     tau2.grid <- profile.fit$tau2.grid
##     q <- qchisq(1-alpha,df=1)
##     if(plot) {
##         plot(mu.grid,profile['LR',],ylim=range(profile[c('LR','LR.bartlett'),]),type='o')
##         lines(mu.grid,profile['LR.bartlett',],lty=2,type='o')
##         abline(h=q)
##     }
##     CIs <- apply(profile[c('LR','LR.bartlett'),] - q, 1, simplify=FALSE, FUN=function(vals) {
##         ## left.idx <- which(diff(sign(vals))!=0)
##         ## stopifnot(length(left.idx)==2)
##         ## stopifnot(left.idx[2]<length(vals))
##         ## stopifnot(vals[left.idx[2]]<0)
##         ## right.idx <- left.idx+1 
##         ## left.wt <- abs(vals[right.idx]/(vals[right.idx]-vals[left.idx]))
##         ## right.wt <- abs(vals[left.idx]/(vals[right.idx]-vals[left.idx]))
##         ## left.wt*mu.grid[left.idx] + right.wt*mu.grid[right.idx]
##         left.idx <- which(diff(sign(vals))!=0)
##         stopifnot(length(left.idx)>0)
##         right.idx <- left.idx+1 
##         left.wt <- abs(vals[right.idx]/(vals[right.idx]-vals[left.idx]))
##         right.wt <- abs(vals[left.idx]/(vals[right.idx]-vals[left.idx]))
##         zeroes <- left.wt*mu.grid[left.idx] + right.wt*mu.grid[right.idx]
##         if(vals[1]<0)zeroes <- c(-Inf,zeroes)
##         if(vals[length(vals)]<0)zeroes <- c(zeroes,Inf)
##         n.zeroes <- length(zeroes)
##         CI.left <- zeroes[(1:n.zeroes) %% 2 ==1]
##         CI.right <- zeroes[(1:n.zeroes) %% 2 ==0]
##         cbind(CI.left,CI.right)
##     } )
##     return(CIs)
## }

## ## vals <- profile['LR.bartlett',]-q
## ## ## stopifnot(vals[1]>0 && vals[length(vals)]>0)
## ## ## signs <- sign(vals)
## ## left.idx <- which(diff(sign(vals))!=0)
## ## right.idx <- left.idx+1 
## ## left.wt <- abs(vals[right.idx]/(vals[right.idx]-vals[left.idx]))
## ## right.wt <- abs(vals[left.idx]/(vals[right.idx]-vals[left.idx]))
## ## zeroes <- left.wt*mu.grid[left.idx] + right.wt*mu.grid[right.idx]
## ## if(vals[1]<0)zeroes <- c(-Inf,zeroes)
## ## if(vals[length(vals)]<0)zeroes <- c(zeroes,Inf)
## ## n.zeroes <- length(zeroes)
## ## CI.left <- zeroes[(1:n.zeroes) %% 2 ==1]
## ## CI.right <- zeroes[(1:n.zeroes) %% 2 ==0]
## ## cbind(CI.left,CI.right)



## ## REFACTORED


## meta.glmm <- function(x,params,model.str,GH=NULL,n.GH=50,logtau2.bounds=c(-20,3),pmf.derivs.fun=NULL,...) {
##     if(is.null(GH)) {
##         GH <- fastGHQuad::gaussHermiteData(n.GH)
##         GH$x <- sqrt(2)*GH$x; GH$w <- GH$w/sqrt(pi)
##         names(GH)[which(names(GH)=='x')] <- 'theta.grid'
##     }
##     loglik.unit <- function(mu,tau2,xi,paramsi,GH) {
##         ## browser()
##         ## print(xi); print(paramsi); print(mu);print(tau2)
##         if(log(tau2)<logtau2.bounds[1] || log(tau2)>logtau2.bounds[2]) {
##             ## warning('tau2 boundary  hit, returning minimum negative')
##             return(-.Machine$double.xmax)
##         }
##         likelihood <- sum(GH$w*cond.pmf(xi,GH$theta.grid,mu,tau2,paramsi))
##         ## print(c(mu,tau2,log(likelihood),paramsi))
##         ## print(tau2)
##         return(log(likelihood))
##     }
##     loglik <- function(mu,tau2) {
##         out <- sum(sapply(seq_along(x),function(i)loglik.unit(mu,tau2,x[i],params[i,,drop=FALSE],GH)))
##         ## print(out)        
##         ## if(is.nan(out))browser()
##         if(is.nan(out))out <- -1e10#-.Machine$double.xmax
##         return(out)
##     }
##     ## load appropriate k formulas if necessary
##     ## if(is.null(pmf.derivs.fun)) {
##     ##     pmf.derivs.filename <- switch(model.str,
##     ##                                   'bn-prop' = 'pmf-derivs-fun-binomial-prop.rds', 'hn'='pmf-derivs-fun-hyperg.rds', 'pn'='pmf-derivs-fun-poisson.rds', 'bn-irr'='pmf-derivs-fun-binomial-irr.rds')
##     ##     pmf.derivs.fun <- readRDS(pmf.derivs.filename)
##     ## }
##     ## assign cond.pmf,dgp,par.start
##     cond.pmf <- switch(model.str,
##                        'bn-prop' = cond.pmf.binomial.prop, 'hn'=cond.pmf.hyperg, 'pn'=cond.pmf.poisson, 'bn-irr'=cond.pmf.binomial.irr)
##     integrand.mode <- switch(model.str,
##                        'bn-prop' = integrand.mode.binomial.prop, 'hn'=integrand.mode.hyperg, 'pn'=integrand.mode.poisson, 'bn-irr'=integrand.mode.binomial.irr)
##     var.map <- switch(model.str,
##                        'bn-prop' = var.map.binomial.prop, 'hn'=var.map.hyperg, 'pn'=var.map.poisson, 'bn-irr'=var.map.binomial.irr)
##     dgp <- switch(model.str,
##                   'bn-prop' = dgp.binomial.prop, 'hn'=dgp.hyperg, 'pn'=dgp.poisson, 'bn-irr'=dgp.binomial.irr)
##     par.start <- switch(model.str,
##                         'bn-prop' = par.start.binomial.prop, 'hn'=par.start.hyperg, 'pn'=par.start.poisson, 'bn-irr'=par.start.binomial.irr)
##     par.start <- par.start(x,params) # clean up, bring par.start here
##     ## compute full mle
##     start.val <- loglik(par.start[1],par.start[2])
##     if(is.nan(start.val) || !is.finite(start.val))browser()
##     fit.full <- optim(par=c(par.start['mu'],log(par.start['tau2'])),fn=function(par)-loglik(par[1],exp(par[2])),method='BFGS',...)
##     if(fit.full$convergence != 0L)warning('full fit convergence failure')
##     mle.full <- fit.full$par
##     mu.mle.full <- mle.full[1]; tau2.mle.full <- exp(mle.full[2])
##     out <- list(
##         model.str=model.str,
##         mle.full=c(mu=unname(mu.mle.full),tau2=unname(tau2.mle.full)),
##         GH=GH,n.GH=n.GH,
##         loglik=loglik,
##         params=params,
##         par.start=par.start,
##         logtau2.bounds=logtau2.bounds,
##         integrand.mode=integrand.mode,
##         var.map=var.map)
##     return(out)
## }

## ## support for adaptive GH
## meta.glmm <- function(x,params,model.str,GH=NULL,n.GH=50,adaptive=TRUE,logtau2.bounds=c(-20,3),pmf.derivs.fun=NULL,...) {
##     if(is.null(GH)&&!adaptive) {
##         GH <- fastGHQuad::gaussHermiteData(n.GH)
##         GH$x <- sqrt(2)*GH$x; GH$w <- GH$w/sqrt(pi)
##         names(GH)[which(names(GH)=='x')] <- 'theta.grid'
##     }
##     ## assign cond.pmf,dgp,etc
##     cond.pmf <- switch(model.str,
##                        'bn-prop' = cond.pmf.binomial.prop, 'hn'=cond.pmf.hyperg, 'pn'=cond.pmf.poisson, 'bn-irr'=cond.pmf.binomial.irr)
##     integrand.mode <- switch(model.str,
##                        'bn-prop' = integrand.mode.binomial.prop, 'hn'=integrand.mode.hyperg, 'pn'=integrand.mode.poisson, 'bn-irr'=integrand.mode.binomial.irr)
##     var.map <- switch(model.str,
##                        'bn-prop' = var.map.binomial.prop, 'hn'=var.map.hyperg, 'pn'=var.map.poisson, 'bn-irr'=var.map.binomial.irr)
##     dgp <- switch(model.str,
##                   'bn-prop' = dgp.binomial.prop, 'hn'=dgp.hyperg, 'pn'=dgp.poisson, 'bn-irr'=dgp.binomial.irr)
##     par.start <- switch(model.str,
##                         'bn-prop' = par.start.binomial.prop, 'hn'=par.start.hyperg, 'pn'=par.start.poisson, 'bn-irr'=par.start.binomial.irr)
##     par.start <- par.start(x,params) # clean up, bring par.start here
##     GH.LR <- function(xi, paramsi, mu, tau2, n.GH=20) {
##         ## browser()
##         rule <- fastGHQuad::gaussHermiteData(n.GH)
##         mode <- integrand.mode(xi,paramsi,mu,tau2)
##         scale <- (1 + tau2 * var.map(mode,paramsi, mu, tau2))^(-1 / 2)
##         theta <- mode + sqrt(2) * scale * rule$x
##         w <- sqrt(2) * scale *
##             rule$w *
##             exp(rule$x^2) *
##             dnorm(theta)
##         return(list(w=w,theta.grid=theta))
## }
##     loglik.unit <- function(mu,tau2,xi,paramsi) {
##         ## browser()
##         ## print(xi); print(paramsi); print(mu);print(tau2)
##         if(adaptive)
##             GH <- GH.LR(xi,paramsi,mu,tau2,n.GH)
##         if(log(tau2)<logtau2.bounds[1] || log(tau2)>logtau2.bounds[2]) {
##             ## warning('tau2 boundary  hit, returning minimum negative')
##             return(-.Machine$double.xmax)
##         }
##         likelihood <- sum(GH$w*cond.pmf(xi,GH$theta.grid,mu,tau2,paramsi))
##         ## print(c(mu,tau2,log(likelihood),paramsi))
##         ## print(tau2)
##         return(log(likelihood))
##     }
##     loglik <- function(mu,tau2) {
##         out <- sum(sapply(seq_along(x),function(i)loglik.unit(mu,tau2,x[i],params[i,,drop=FALSE])))
##         ## print(out)        
##         ## if(is.nan(out))browser()
##         if(is.nan(out))out <- -1e10#-.Machine$double.xmax
##         return(out)
##     }
##     ## compute full mle
##     ## browser()
##     start.val <- loglik(par.start[1],par.start[2])
##     if(is.nan(start.val) || !is.finite(start.val))browser()
##     fit.full <- optim(par=c(par.start['mu'],log(par.start['tau2'])),fn=function(par)-loglik(par[1],exp(par[2])),method='BFGS',...)
##     if(fit.full$convergence != 0L)warning('full fit convergence failure')
##     mle.full <- fit.full$par
##     mu.mle.full <- mle.full[1]; tau2.mle.full <- exp(mle.full[2])
##     out <- list(
##         model.str=model.str,
##         mle.full=c(mu=unname(mu.mle.full),tau2=unname(tau2.mle.full)),
##         GH=GH,n.GH=n.GH,
##         loglik=loglik,
##         params=params,
##         par.start=par.start,
##         logtau2.bounds=logtau2.bounds,
##         integrand.mode=integrand.mode,
##         var.map=var.map)
##     return(out)
## }


## LR.meta.glmm <- function(meta.glmm.fit,mu.null) {
##     ## browser()
##     loglik <- meta.glmm.fit$loglik
##     x <- meta.glmm.fit$x
##     mu.mle.full <- meta.glmm.fit$mle.full['mu']
##     tau2.mle.full <- meta.glmm.fit$mle.full['tau2']
##     par.start <- meta.glmm.fit$par.start
##     logtau2.bounds <- meta.glmm.fit$logtau2.bounds
##     logtau2.mle.redu <- optimize(function(logtau2)loglik(mu.null,exp(logtau2)),interval=logtau2.bounds,maximum=TRUE)$max
##     tau2.mle.redu <- exp(logtau2.mle.redu)
##     LR <- 2 * (loglik(mu.mle.full,tau2.mle.full) - loglik(mu.null,tau2.mle.redu))
##     out <- list(mle.redu=c(mu=unname(mu.null),tau2=unname(tau2.mle.redu)), LR=LR)
##     return(out)
##     ## return(list(mle.full=c(mu=unname(mu.mle.full),tau2=unname(tau2.mle.full)),mle.redu=c(mu=unname(mu.null),tau2=unname(tau2.mle.redu)), LR=LR, loglik.fun=loglik, x=x, params=params, cond.pmf=cond.pmf,logtau2.bounds=logtau2.bounds, n.GH=n.GH))
## }
## bartlett.meta.glmm <- function(meta.glmm.fit,params,mu,tau2,n.GH=50,debug=FALSE) {
##     cat('.')
##     ## ks.fun <- readRDS('ks-fun.rds')
##     pmf.derivs.fun <- switch(meta.glmm.fit$model.str,
##                              'bn-prop' = pmf.derivs.fun.binomial.prop, 'hn'=pmf.derivs.fun.hyperg, 'pn'=pmf.derivs.fun.poisson, 'bn-irr'=pmf.derivs.fun.binomial.irr)
##     ## browser()
##     ks.ma.fit <- ks.ma(params=meta.glmm.fit$params,ks.fun=ks.fun,pmf.derivs.fun=pmf.derivs.fun,mu=mu,tau2=tau2,n.GH=n.GH,debug=debug)
##     ks.df.ma <- ks.ma.fit$ks.df.ma
##     bartlett.C <- lawley(ks.df.ma,2) - lawley(ks.df.ma,1)
##     out <- list(bartlett.C=bartlett.C)
##     if(debug) out <- c(out, debug=list(ks.ma.fit$debug))
##     ## browser()
##     return(out)
## }
## profile.meta.glmm <- function(meta.glmm.fit,grid=list(offset=1,size=20),mc.cores=1,bartlett=TRUE,debug=FALSE) {
##     ## browser()
##     x <- meta.glmm.fit$x
##     params <- meta.glmm.fit$params
##     cond.pmf <- meta.glmm.fit$cond.pmf
##     logtau2.bounds <- meta.glmm.fit$logtau2.bounds
##     GH <- meta.glmm.fit$GH
##     mu.mle.full <- meta.glmm.fit$mle.full['mu']; tau2.mle.full <- meta.glmm.fit$mle.full['tau2'];
##     loglik <- meta.glmm.fit$loglik.fun
##     pmf.derivs.fun <- meta.glmm.fit$pmf.derivs.fun
##     mu.profile.grid <- seq(mu.mle.full-grid$offset*sqrt(tau2.mle.full),mu.mle.full+grid$offset*sqrt(tau2.mle.full),len=grid$size)
##     profile.fits <- lapply(mu.profile.grid, function(mu.null)LR.meta.glmm(meta.glmm.fit,mu.null))
##     ## LR.profile <- sapply(profile.fits,function(fit)fit$LR)
##     LR.profile <- sapply(profile.fits,`[[`,'LR')
##     out <- list(LR=LR.profile,mu.grid=mu.profile.grid)
##     if(bartlett) {
##         tau2.profile.grid <- sapply(profile.fits,function(fit)fit$mle.redu['tau2'])
##         bartlett.C.fits <- parallel::mclapply(seq_along(mu.profile.grid),mc.cores=mc.cores,FUN=function(i)
##         ## bartlett.C.fits <- lapply(seq_along(mu.profile.grid),FUN=function(i)
##             ## bartlett.LR(params=params,pmf.derivs.fun=pmf.derivs.fun,mu=mu.profile.grid[i],tau2=tau2.profile.grid[i],n.GH=n.GH,debug=debug))
##             bartlett.meta.glmm(meta.glmm.fit,params=params,mu=mu.profile.grid[i],tau2=tau2.profile.grid[i],n.GH=n.GH,debug=debug))
##         bartlett.C.profile <- sapply(bartlett.C.fits,`[[`,'bartlett.C')
##         LR.bartlett.profile <- LR.profile / (1+bartlett.C.profile)
##         out <- c(out, list(LR.bartlett=LR.bartlett.profile,tau2.grid=tau2.profile.grid))
##         if(debug) {
##             debug.profile <- lapply(bartlett.C.fits,`[[`,'debug')
##             out <- c(out,debug=list(debug.profile))
##             ## browser()
##         }
##     }
##     return(out)
## }
## confint.meta.glmm <- function(meta.glmm.fit,alpha=.05,plot=TRUE,debug=FALSE,bartlett=TRUE,grid=list(offset=1,size=20),mc.cores=1,...) {
##     ## browser()
##     profile.fit <- profile.meta.glmm(meta.glmm.fit,grid=grid,mc.cores=mc.cores,bartlett=bartlett,debug=debug)
##     LR.profile <- profile.fit$LR
##     LR.bartlett.profile <- profile.fit$LR.bartlett
##     ## profile <- profile.fit$profile
##     mu.grid <- profile.fit$mu.grid
##     ## tau2.grid <- profile.fit$tau2.grid
##     profile <- rbind(LR=LR.profile,LR.bartlett=LR.bartlett.profile)
##     q <- qchisq(1-alpha,df=1)
##     ## browser()
##     if(plot) {
##         plot(mu.grid,profile['LR',],ylim=range(profile),type='o')
##         if(bartlett) lines(mu.grid,profile['LR.bartlett',],lty=2,type='o')
##         abline(h=q)
##     }
##     ## interpolate CIs
##     CIs <- apply(profile - q, 1, simplify=FALSE, FUN=function(vals) {
##         left.idx <- which(diff(sign(vals))!=0)
##         stopifnot(length(left.idx)>0)
##         right.idx <- left.idx+1 
##         left.wt <- abs(vals[right.idx]/(vals[right.idx]-vals[left.idx]))
##         right.wt <- abs(vals[left.idx]/(vals[right.idx]-vals[left.idx]))
##         zeroes <- left.wt*mu.grid[left.idx] + right.wt*mu.grid[right.idx]
##         if(vals[1]<0)zeroes <- c(-Inf,zeroes)
##         if(vals[length(vals)]<0)zeroes <- c(zeroes,Inf)
##         n.zeroes <- length(zeroes)
##         CI.left <- zeroes[(1:n.zeroes) %% 2 ==1]
##         CI.right <- zeroes[(1:n.zeroes) %% 2 ==0]
##         cbind(CI.left,CI.right)
##     } )
##     out <- list(CIs=CIs)
##     if(debug) out <- c(out,profile.fit$debug)
##     return(out)
## }


