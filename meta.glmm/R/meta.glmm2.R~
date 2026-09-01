

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
        ## cat('.')
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
make.quad <- function(moment.map,cond.pmf,n.GH,adaptive=TRUE,adaptive.tol=NULL) {
    GH <- fastGHQuad::gaussHermiteData(n.GH)
    if(!adaptive) {
     quad <- function(...) {
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
                ## zi.new <- zi + newton.step(zi, yi, paramsi, mu, tau2)
                ## if(!is.finite(zi.new))browser()
                ## 1+1
                ## print(zi.new)
                ## print(c(zi.new,zi))
                step <- newton.step(zi, yi, paramsi, mu, tau2)
                zi.new <- zi + step
                ## if(!is.finite(cond.pmf(yi,zi.new,mu,tau2,paramsi,log=TRUE)))browser()
                while(cond.pmf(yi,zi.new,mu,tau2,paramsi,log=TRUE) < cond.pmf(yi,zi,mu,tau2,paramsi,log=TRUE)) {
                    step <- step / 2
                    zi.new <- zi + step
                }
                if (abs(zi.new - zi) < adaptive.tol) {
                    break
                }
                zi <- zi.new
            }
            mode <- zi.new
            return(mode)
        }
        quad <- function(yi,paramsi,mu,tau2) {
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
    return(quad)
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
    ## fit <- optim(par=c(par.start['mu'],log(par.start['tau2'])),fn=function(par)-loglik(par[1],exp(par[2])),...)
    ## mu.mle <- unname(fit$par[1]); tau2.mle <- unname(exp(fit$par[2]))
    ## out <- list(mle=c(mu=mu.mle,tau2=tau2.mle), status=fit$convergence)
    obj <- function(par)-loglik(par[1],par[2])
    fit <- nloptr::nloptr(
                       x0     = par.start,
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
    out <-  list(mle=c(mu=fit$solution[1],tau2=fit$solution[2]), status=fit$status)
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
    out <- list(CIs=CIs,mu.grid=mu.grid,profile=profile)
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
    quad <- make.quad(model$moment.map,model$cond.pmf,n.GH=n.GH,adaptive=adaptive,adaptive.tol=1e-4)#!!!!!!
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
    mle0 <- meta.glmm.fit$mle
    mu.grid <- seq(mle0['mu']-mu.grid$offset*sqrt(mle0['tau2']),mle0['mu']+mu.grid$offset*sqrt(mle0['tau2']),len=mu.grid$size)    
    profile0 <- profile(mu.grid,meta.glmm.fit$profile.LR,mle0,if(bartlett)meta.glmm.fit$bartlett.C else NULL)
    out <- confint(mu.grid=mu.grid, profile=rbind(LR=profile0$LR,LR.bartlett=profile0$LR.bartlett), alpha=alpha, plot=plot,...)
    return(out)
}


