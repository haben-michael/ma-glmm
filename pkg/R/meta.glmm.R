

ks.study <- function(paramsi,mu,tau2,quad,pmf.derivs.fun,debug=FALSE)  {
    beta <- c(tau2,mu)
    supp.y <-  paramsi['y.min']:paramsi['y.max']
    quad.rules <- lapply(supp.y, function(yi)quad(yi,paramsi,mu,tau2))
    weight.mat <- sapply(quad.rules,`[[`,'wt')
    units.mats <- lapply(pmf.derivs.fun, function(pmf.deriv.fun) {
        units.mat <- sapply(seq_along(supp.y), function(i)pmf.deriv.fun(beta=beta,params=paramsi,theta=quad.rules[[i]]$z,x=supp.y[i]))
        weight.mat*units.mat
    })
    units.mats <- simplify2array(units.mats)
    pmf.derivs.unit <- apply(units.mats,2:3,sum)
    ks.unit <- sapply(seq_along(ks.fun), function(i) {
        k.fun <- ks.fun[[i]]
        args <- attr(k.fun,'args')
        k.unit <- do.call(k.fun, as.data.frame(pmf.derivs.unit[,args]))
    })
    k.vals <- colSums(ks.unit)
    k.orders <- sapply(ks.fun,function(k.fun)attr(k.fun,'order'))
    ks.df <- data.frame(inner=k.orders[1,],outer=k.orders[2,],k.val=k.vals)
    out <- list(ks.df=ks.df)
    if(debug) out <- c(out,debug=list(list(ks.unit=ks.unit, units.mats=units.mats, nodes=sapply(quad.rules,`[[`,'z'))))
    return(out)
}

ks.ma <- function(params,mu,tau2,quad,pmf.derivs.fun,debug=FALSE) {
    ks.study.fits <- lapply(1:nrow(params), function(i) {
        ks.study.fit <- ks.study(pmf.derivs.fun=pmf.derivs.fun,paramsi=params[i,],mu=mu,tau2=tau2,quad=quad,debug=debug)
    })
    k.vals.ma <- sapply(ks.study.fits,function(fit)fit$ks.df$k.val) |> rowSums()
    ks.df.ma <- data.frame(inner=ks.study.fits[[1]]$ks.df$inner, outer=ks.study.fits[[1]]$ks.df$outer,k.val=k.vals.ma)
    ks.study.fits <- sapply(ks.study.fits,`[[<-`,'ks.df',NULL)
    out <- list(ks.df=ks.df.ma)
    if(debug) out <- c(out,debug=list(ks.study.fits))
    return(out)
}

lawley <- function(k.df,p) {
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
        out <- list(bartlett.C=bartlett.C)
        if(debug) out <- c(out,debug=list(ks.ma.fit$debug))
        return(out)
    }
}


## 1. loop iterations for N-R is hard-coded
## 2. make bisection in N-R optional
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
            zi <- 0 
            for (j in 1:20) { #!!!!
                step <- newton.step(zi, yi, paramsi, mu, tau2)
                zi.new <- zi + step
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
            mode <- integrand.mode(yi,paramsi,mu,tau2)
            var.star <- moment.map(mode,paramsi,mu,tau2)['second']
            scale <- (1 + tau2 * var.star)^(-1 / 2)
            z <- mode + sqrt(2) * scale * GH$x
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
        lsummands <- rule$log.wt + cond.pmf(yi,rule$z,mu,tau2,paramsi,log=TRUE)
        max.lsummands <- max(lsummands)
        out <- max.lsummands + log(sum(exp(lsummands - max.lsummands)))
        return(out)
    }
    loglik <- function(mu,tau2) {
        out <- sum(sapply(seq_along(y),function(i)loglik.unit(mu,tau2,y[i],params[i,])))
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


## 1. logtau2.bounds not currently used--from old optim code
## 2. parameters passed to nloptr are hard-coded
mle <- function(loglik,par.start,logtau2.bounds=c(-20,3),...) {
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

## 1. update to return an S3 object that works with native confint
profile <- function(mu.grid,profile.LR,mle,bartlett.C=NULL,mc.cores=1) {
    profile.fits <- lapply(mu.grid, function(mu.null)profile.LR(mu.null))
    LR.profile <- sapply(profile.fits,`[[`,'LR')
    tau2.grid <- sapply(profile.fits,`[[`,'mle.null')['tau2',]
    out <- list(LR=LR.profile,mu.grid=mu.grid,tau2.grid=tau2.grid)
    if(!is.null(bartlett.C)) {
        bartlett.C.fits <- parallel::mclapply(seq_along(mu.grid),mc.cores=mc.cores,FUN=function(i)    bartlett.C(mu.grid[i],tau2.grid[i]))
        bartlett.C.profile <- sapply(bartlett.C.fits,`[[`,'bartlett.C')
        LR.bartlett.profile <- LR.profile / (1+bartlett.C.profile)
        out <- c(out, list(LR.bartlett=LR.bartlett.profile,bartlett.C=bartlett.C.profile))
    }
    return(out)
}

## 1. just use R's confint which has better interpolation, need to update profile() first
confint <- function(mu.grid, profile, alpha=.05, plot=TRUE,...) {
    q <- qchisq(1-alpha,df=1)
    if(plot) {
        matplot(mu.grid,t(profile),type='l',col=1,lty=2:1,xlab=bquote(mu),ylab='profile LR',...)
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

## 1. choice of starting tau2 when tau2.DL==0 is currently var(y)
par.start <- function(y,params,nn.transform) {
    with(nn.transform(y,params), {
        w.fe <- 1/sigma2
        mu.fe <- y %*% w.fe / sum(w.fe)    
        Q <- w.fe%*%(y-c(mu.fe))^2
        tau2.DL <- max(0L, (Q - (length(y)-1)) / (sum(w.fe)-sum(w.fe^2)/sum(w.fe)))
        tau2.start <- if(tau2.DL==0L) var(y) else tau2.DL #!!!!!
        c(mu=mu.fe,tau2=tau2.start)
    })
}

## 1. adaptive.tol is hard-coded. change parameters so adaptive.tol is
## passed, defaults to NULL which is fixed rule .
meta.glmm <- function(y,params,model.str,n.GH,adaptive,logtau2.bounds,debug=FALSE,...) {
    model <- model(model.str)
    quad <- make.quad(model$moment.map,model$cond.pmf,n.GH=n.GH,adaptive=adaptive,adaptive.tol=1e-6)#!!!!!!
    loglik <- make.loglik(model$cond.pmf,params,y,quad)
    par.start0 <- par.start(y,params,model$nn.transform)
    mle0 <- mle(loglik,par.start0,logtau2.bounds,...)$mle
    profile.LR <- make.profile.LR(loglik,mle0,par.start0,logtau2.bounds=logtau2.bounds)
    bartlett.C <- make.bartlett.C(params,quad,model$pmf.derivs.fun,debug=debug)
    out <- list(model.str=model.str, quad=quad, loglik=loglik, mle=mle0, profile.LR=profile.LR, bartlett.C=bartlett.C, par.start=par.start0)
    return(out)
}

confint.meta.glmm <- function(meta.glmm.fit,alpha=.05,plot=TRUE,debug=FALSE,bartlett=TRUE,mu.grid=list(offset=1.5,size=20),mc.cores=1,...) {
    mle0 <- meta.glmm.fit$mle
    mu.grid <- seq(mle0['mu']-mu.grid$offset*sqrt(mle0['tau2']),mle0['mu']+mu.grid$offset*sqrt(mle0['tau2']),len=mu.grid$size)    
    profile0 <- profile(mu.grid,meta.glmm.fit$profile.LR,mle0,if(bartlett)meta.glmm.fit$bartlett.C else NULL)
    out <- confint(mu.grid=mu.grid, profile=rbind(LR=profile0$LR,LR.bartlett=profile0$LR.bartlett), alpha=alpha, plot=plot,...)
    return(out)
}


diagnostic.pmf <- function(bartlett.fit, study.no,...) {
    op <- par(mfrow=c(3,5))
    for(i in 1:15){
        deriv <- attr(bartlett.fit$debug[[study.no]]$units.mat,'dimnames')[[3]][i]
        label <- paste0('(',paste(strsplit(gsub('[^0-9]','',deriv),'')[[1]],collapse=','),')')
        matplot(bartlett.fit$debug[[study.no]]$units.mats[,,i],pch=1,type='l',lty=1,col=1,ylab='',xlab=label,...)
    }
    par(op)
}

## 1. currently uses indices on x axis, can pass support of y as ... . may need to reorganize so bartlett.fit gets the params.
## 2. maybe add asinh transformation option
diagnostic.k <- function(bartlett.fit,study.no,...)
    matplot(bartlett.fit$debug[[study.no]]$ks.unit[,],type='l',pch=1,col=1,xlab='y',ylab='k',lty=1,...)
