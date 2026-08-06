
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






mma.to.R <- function(
    x,
    funs = c(
      "Sqrt", "Log", "Exp",
      "Sin", "Cos", "Tan",
      "ArcSin", "ArcCos", "ArcTan",
      "Power",
      "LogisticSigmoid",
      "Binomial"
    ),
    rename = c(
      Sqrt = "sqrt",
      Log = "log",
      Exp = "exp",
      Sin = "sin",
      Cos = "cos",
      Tan = "tan",
      ArcSin = "asin",
      ArcCos = "acos",
      ArcTan = "atan",
      LogisticSigmoid = "plogis",
      Binomial = "choose"
    )
) {
  stopifnot(length(x) == 1L)

  chars <- strsplit(x, "", fixed = TRUE)[[1]]
  n <- length(chars)

  is_ident_start <- function(z) grepl("[A-Za-z.$`]", z)
  is_ident_char  <- function(z) grepl("[A-Za-z0-9.$`]", z)
  is_space       <- function(z) grepl("\\s", z)

  split_top_level <- function(text) {
    text <- trimws(text)

    if (!nzchar(text)) {
      return(character())
    }

    z <- strsplit(text, "", fixed = TRUE)[[1]]

    paren_depth   <- 0L
    bracket_depth <- 0L
    brace_depth   <- 0L
    start <- 1L
    pieces <- character()

    for (i in seq_along(z)) {
      ch <- z[i]

      if (ch == "(") {
        paren_depth <- paren_depth + 1L
      } else if (ch == ")") {
        paren_depth <- paren_depth - 1L
      } else if (ch == "[") {
        bracket_depth <- bracket_depth + 1L
      } else if (ch == "]") {
        bracket_depth <- bracket_depth - 1L
      } else if (ch == "{") {
        brace_depth <- brace_depth + 1L
      } else if (ch == "}") {
        brace_depth <- brace_depth - 1L
      } else if (
        ch == "," &&
        paren_depth == 0L &&
        bracket_depth == 0L &&
        brace_depth == 0L
      ) {
        pieces <- c(
          pieces,
          trimws(substr(text, start, i - 1L))
        )
        start <- i + 1L
      }
    }

    c(
      pieces,
      trimws(substr(text, start, nchar(text)))
    )
  }
  parse_iterator <- function(text) {
    text <- trimws(text)

    if (
      nchar(text) < 2L ||
      substr(text, 1L, 1L) != "{" ||
      substr(text, nchar(text), nchar(text)) != "}"
    ) {
      stop("Unsupported Sum iterator: ", text)
    }

    inside <- substr(text, 2L, nchar(text) - 1L)
    parts <- split_top_level(inside)

    if (length(parts) == 2L) {
      ans <- list(
        variable = trimws(parts[1L]),
        lower = "1",
        upper = trimws(parts[2L]),
        step = NULL
      )
    } else if (length(parts) == 3L) {
      ans <- list(
        variable = trimws(parts[1L]),
        lower = trimws(parts[2L]),
        upper = trimws(parts[3L]),
        step = NULL
      )
    } else if (length(parts) == 4L) {
      ans <- list(
        variable = trimws(parts[1L]),
        lower = trimws(parts[2L]),
        upper = trimws(parts[3L]),
        step = trimws(parts[4L])
      )
    } else {
      stop("Unsupported Sum iterator: ", text)
    }

    if (!grepl("^[A-Za-z.$`][A-Za-z0-9.$`]*$", ans$variable)) {
      stop("Invalid Sum index variable: ", ans$variable)
    }

    ans
  }

  make_sequence <- function(iterator) {
    if (is.null(iterator$step)) {
      if (iterator$lower == "1") {
        paste0("seq_len(", iterator$upper, ")")
      } else {
        paste0(
          "seq.int(",
          iterator$lower, ",",
          iterator$upper,
          ")"
        )
      }
    } else {
      paste0(
        "seq.int(",
        iterator$lower, ",",
        iterator$upper,
        ",by=", iterator$step,
        ")"
      )
    }
  }

  convert_one_sum <- function(expr, iterator) {
    expr <- trimws(expr)
    sequence <- make_sequence(iterator)

    match <- regexec(
      "^([A-Za-z.$`][A-Za-z0-9.$`]*)\\[([A-Za-z.$`][A-Za-z0-9.$`]*)\\]$",
      expr,
      perl = TRUE
    )

    fields <- regmatches(expr, match)[[1]]

    if (
      length(fields) == 3L &&
      fields[3L] == iterator$variable
    ) {
      vector_name <- fields[2L]

      if (
        iterator$lower == "1" &&
        iterator$upper == "n" &&
        is.null(iterator$step)
      ) {
        return(paste0("sum(", vector_name, ")"))
      }

      return(
        paste0(
          "sum(",
          vector_name, "[", sequence, "]",
          ")"
        )
      )
    }

    paste0(
      "sum(vapply(",
      sequence,
      ",function(", iterator$variable, ") ",
      expr,
      ",numeric(1)))"
    )
  }

  convert_sum <- function(inner_text) {
    arguments <- split_top_level(inner_text)

    if (length(arguments) < 2L) {
      stop("Sum requires a summand and at least one iterator")
    }

    expr <- arguments[1L]
    iterators <- lapply(arguments[-1L], parse_iterator)

    for (iterator in rev(iterators)) {
      expr <- convert_one_sum(expr, iterator)
    }

    expr
  }

  parse_until <- function(i, close = NULL) {
    out <- character()

    while (i <= n) {
      ch <- chars[i]

      if (!is.null(close) && ch == close) {
        return(
          list(
            text = paste0(out, collapse = ""),
            i = i + 1L
          )
        )
      }

      if (is_ident_start(ch)) {
        start <- i
        i <- i + 1L

        while (i <= n && is_ident_char(chars[i])) {
          i <- i + 1L
        }

        head <- paste0(chars[start:(i - 1L)], collapse = "")

        # Allow whitespace between symbol and bracket.
        ws_start <- i

        while (i <= n && is_space(chars[i])) {
          i <- i + 1L
        }

        ws <- if (i > ws_start) {
          paste0(chars[ws_start:(i - 1L)], collapse = "")
        } else {
          ""
        }

        if (i <= n && chars[i] == "[") {
          inner <- parse_until(i + 1L, close = "]")

          if (head == "Sum") {
            out <- c(out, convert_sum(inner$text))
          } else if (head %in% funs) {
            new_head <- if (head %in% names(rename)) {
              unname(rename[[head]])
            } else {
              head
            }

            out <- c(out, new_head, "(", inner$text, ")")
          } else {
            # Preserve indexing/application brackets.
            out <- c(out, head, ws, "[", inner$text, "]")
          }

          i <- inner$i
        } else {
          out <- c(out, head, ws)
        }

      } else if (ch == "[") {
        inner <- parse_until(i + 1L, close = "]")
        out <- c(out, "[", inner$text, "]")
        i <- inner$i

      } else {
        out <- c(out, ch)
        i <- i + 1L
      }
    }

    if (!is.null(close)) {
      stop("Unmatched bracket: expected closing ", close)
    }

    list(
      text = paste0(out, collapse = ""),
      i = i
    )
  }

  parse_until(1L)$text
}
## make.fun <- function(str, args, envir = parent.frame()) {
##   expr <- parse(text = str)[[1L]]

##   f <- eval(
##     call(
##       "function",
##       as.pairlist(setNames(rep(list(quote(expr = )), length(args)), args)),
##       expr
##     ),
##     envir = envir
##   )

##   f
## }

str.to.fun <- function(str, args, bindings = list(), parent = baseenv()) {
  expr <- parse(text = str)[[1L]]

  env <- list2env(bindings, parent = parent)

  f <- eval(
    call(
      "function",
      as.pairlist(setNames(rep(list(quote(expr = )), length(args)), args)),
      expr
    ),
    envir = env
  )

  environment(f) <- env
  f
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
.cond.pmf.hyperg <- function(x,theta,mu,tau2,params) {
    ## x <- 0; theta <-  -18.95964; mu <- 0; tau2 <- 0.3071743;params <- c(m.T=116,m.C=117,s=3,x.min=0,x.max=3)
        m.T <- params[,'m.T']
        m.C <- params[,'m.C']
        s <- params[,'s']
    odds <- exp(mu+sqrt(tau2)*theta)
    ## print(x);print(s);print(params)
    ## browser()
        BiasedUrn::dFNCHypergeo(x=x,m1=m.T,m2=m.C,n=s,odds=odds)            
}
cond.pmf.hyperg <-  Vectorize(.cond.pmf.hyperg,vectorize.args=c("x","theta"))
## cond.pmf.hyperg <- function(x,theta,mu,tau2,params)  <- sapply(theta, function(thetai) _cond.pmf.hyperg(x,thetai,mu,tau2,params)

## cond.pmf.hyperg <- function(x,theta,mu,tau2,params) {
##     m.T <- params[,'m.T']
##     m.C <- params[,'m.C']
##     s <- params[,'s']
##     numer <- choose(m.T,x)*choose(m.C,s-x)*exp((mu+sqrt(tau2)*theta)*x)
##     supp <- max(0,s-m.C):min(s,m.T)
##     denom <- rowSums(  outer(theta,supp,function(thetaj,xj)choose(m.T,xj)*choose(m.C,s-xj)*exp((mu+sqrt(tau2)*thetaj)*xj) ) )
##     numer / denom
## }
dgp.hyperg <- function(n,mu,tau2,params) {
    m.T <- params[,'m.T']
    m.C <- params[,'m.C']
    s <- params[,'s']
    theta <- rnorm(n)
    odds <- exp(mu+sqrt(tau2)*theta)
    x <- x.T <- sapply(1:n, function(i) BiasedUrn::rFNCHypergeo(1,m1=m.T[i],m2=m.C[i],n=s[i],odds=odds[i]))
}




cond.pmf.poisson <- function(x,theta,mu,tau2,params) {
    T <- params[,'T']
    lambda <- T*exp(mu+sqrt(tau2)*theta)
    dpois(x,lambda=lambda)
}
dgp.poisson <- function(n,mu,tau2,params,re.distr=rnorm) {
    T <- params[,'T']
    theta <- re.distr(n)
    lambda <- T*exp(mu+sqrt(tau2)*theta)
    rpois(n,lambda)
    ## pr <- function(theta)plogis( mu+sqrt(tau2)*theta)
    ## rbinom(n,size=m,prob=pr(theta))
}



cond.pmf.binomial.prop <- function(x,theta,mu,tau2,params) {
    m <- params[,'m']
    pr <- function(theta)plogis(mu+sqrt(tau2)*theta)
    dbinom(x,size=m,prob=pr(theta))
    ## integrand <- function(theta)pr(theta)^xi*(1-pr(theta))^(mi-xi)
}
dgp.binomial.prop <- function(n,mu,tau2,params,re.distr=rnorm) {
    m <- params[,'m']
    theta <- re.distr(n)
    pr <- function(theta)plogis( mu+sqrt(tau2)*theta)
    rbinom(n,size=m,prob=pr(theta))
}



cond.pmf.binomial.irr <- function(x,theta,mu,tau2,params) {
    T.t <- params[,'T.t']
    T.c <- params[,'T.c']
    s <- params[,'s']
    pr <- function(theta)plogis(log(T.t/T.c) + mu+sqrt(tau2)*theta)
    dbinom(x,size=s,prob=pr(theta))
}
dgp.binomial.irr <- function(n,mu,tau2,params,re.distr=rnorm) {
    T.t <- params[,'T.t']
    T.c <- params[,'T.c']
    s <- params[,'s']
    theta <- re.distr(n)
    pr <- function(theta)plogis(log(T.t/T.c) + mu+sqrt(tau2)*theta)
    rbinom(n,size=s,prob=pr(theta))
}




## ##  rework ks.fun to store order and args as attributes not list
## ## elements.
## ks.study <- function(pmf.derivs.fun,ks.fun,params,supp.x,mu.null,tau2.null,n.GH=50) {
##     ## browser()
##     beta <- c(tau2.null,mu.null)
##     ## if(is.null(GH)) {
##     GH <- fastGHQuad::gaussHermiteData(n.GH)
##     GH$x <- sqrt(2)*GH$x
##     GH$w <- GH$w / sqrt(pi)
##     names(GH) <- c('theta.grid','w')
##     ## }
##     pmf.derivs.unit <- lapply(pmf.derivs.fun, function(pmf.deriv.fun) {
##     ## browser()
##         ## supp <- max(0,s-m.C):min(s,m.T)
##         ## g.mat <- outer(X=GH$theta.grid,Y=supp,function(theta,x)pmf.deriv.fun.ma(theta,x,m,m.C,s)) # these are globals!!!
##         ## g.mat <- outer(X=GH$theta.grid,Y=supp.x,function(theta,x)pmf.deriv.fun(beta=beta,theta=theta,x=x,params=params)) 
##         ## g.unit <- colSums(GH$w * g.mat)
##         units.mat <- sapply(GH$theta,function(thetai)pmf.deriv.fun(beta=beta,params=params,theta=thetai,x=supp.x))
##         colSums(GH$w * t(units.mat))
##     })
##     k.vals <- sapply(seq_along(ks.fun), function(i) {
##         k.fun <- ks.fun[[i]]
##         args <- attr(k.fun,'args')
##         k.fun.unit <- do.call(k.fun, pmf.derivs.unit[args])
##         k.val <- sum( k.fun.unit )
##     })
##     k.orders <- sapply(ks.fun,function(k.fun)attr(k.fun,'order'))
##     ## k.df <- data.frame(inner=sapply(k.orders,`[`,1),outer=sapply(k.orders,`[`,2),k.val=k.vals)
##     ks.df <- data.frame(inner=k.orders[1,],outer=k.orders[2,],k.val=k.vals)
##     return(ks.df)
## }


ks.study <- function(pmf.derivs.fun,ks.fun,params,mu,tau2,n.GH=50) {
    ## browser()
    beta <- c(tau2,mu)
    ## if(is.null(GH)) {
    GH <- fastGHQuad::gaussHermiteData(n.GH)
    GH$x <- sqrt(2)*GH$x
    GH$w <- GH$w / sqrt(pi)
    names(GH) <- c('theta.grid','w')
    ## }
    pmf.derivs.unit <- lapply(pmf.derivs.fun, function(pmf.deriv.fun) {
        ## browser()
        supp.x <-  params['x.min']:params['x.max']
        units.mat <- sapply(GH$theta.grid,function(thetai)pmf.deriv.fun(beta=beta,params=params,theta=thetai,x=supp.x))
        colSums(GH$w * t(units.mat))
    })
    k.vals <- sapply(seq_along(ks.fun), function(i) {
        k.fun <- ks.fun[[i]]
        args <- attr(k.fun,'args')
        k.unit <- do.call(k.fun, pmf.derivs.unit[args])
        k.val <- sum( k.unit )
    })
    ks.unit <- sapply(seq_along(ks.fun), function(i) {
        k.fun <- ks.fun[[i]]
        args <- attr(k.fun,'args')
        k.unit <- do.call(k.fun, pmf.derivs.unit[args])
        ## k.val <- sum( k.unit )
    })
    k.vals <- colSums(ks.unit) 
    k.orders <- sapply(ks.fun,function(k.fun)attr(k.fun,'order'))
    ks.df <- data.frame(inner=k.orders[1,],outer=k.orders[2,],k.val=k.vals)
    out <- list(ks.df=ks.df, ks.unit=ks.unit, pmf.derivs.unit=pmf.derivs.unit)
    return(out)
}


ks.ma <- function(params,ks.fun,pmf.derivs.fun,mu,tau2,n.GH=50) {
    ## browser()
    ks <- lapply(1:nrow(params), function(i) {
        ## if(i==16)browser()
        supp.x <- params[i,'x.min']:params[i,'x.max']
        ks.df.study <- ks.study(pmf.derivs.fun=pmf.derivs.fun,ks.fun=ks.fun,params=params[i,],mu,tau2,n.GH=n.GH)$ks.df
    })
    k.vals.ma <- sapply(ks,function(df)df$k.val) |> rowSums()
    ks.df.ma <- data.frame(inner=ks[[1]]$inner, outer=ks[[1]]$outer,k.val=k.vals.ma)
    return(ks.df.ma)
}

## with debugging info
ks.ma <- function(params,ks.fun,pmf.derivs.fun,mu,tau2,n.GH=50) {
    ## browser()
    ks.study.fits <- lapply(1:nrow(params), function(i) {
        ## if(i==16)browser()
        supp.x <- params[i,'x.min']:params[i,'x.max']
        ks.study.fit <- ks.study(pmf.derivs.fun=pmf.derivs.fun,ks.fun=ks.fun,params=params[i,],mu,tau2,n.GH=n.GH)
    })
    ## k.vals.ma <- sapply(ks,function(df)df$k.val) |> rowSums()
    k.vals.ma <- sapply(ks.study.fits,function(fit)fit$ks.df$k.val) |> rowSums()
    ks.df.ma <- data.frame(inner=ks.study.fits[[1]]$ks.df$inner, outer=ks.study.fits[[1]]$ks.df$outer,k.val=k.vals.ma)
    return(list(ks.df.ma=ks.df.ma,ks.study.fits=ks.study.fits))
}


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
    bias <- out <- sum(lawley.4) - sum(lawley.6)
    return(bias)
}


## work with cond pmf, compute unit loglik inside function
LR.gaussian <- function(x,params,mu.null,cond.pmf,logtau2.bounds=c(-5,5),n.GH=50) {
    ## browser()
        GH <- fastGHQuad::gaussHermiteData(n.GH)
        GH$x <- sqrt(2)*GH$x; GH$w <- GH$w/sqrt(pi)
        ## loglik.unit <- function(mu,tau2,xi,paramsi,GH=NULL)loglik.unit(mu,tau2,xi,paramsi,cond.pmf,GH)
        ## loglik.unit <- function(mu,tau2,xi,paramsi,GH=GH) {
        ##     ## browser()
        ##     if(is.null(GH)) {
        ##         likelihood <- integrate(function(theta)cond.pmf(xi,theta,mu,tau2,paramsi)*dnorm(theta),-Inf,Inf)$val 
        ##     } else {
        ##         likelihood <- sum(GH$w*cond.pmf(xi,GH$x,mu,tau2,paramsi))
        ##     }
        ##     return(log(likelihood))
        ## }
        loglik.unit <- function(mu,tau2,xi,paramsi,GH) {
            ## browser()
            ## print(xi); print(paramsi); print(mu);print(tau2)
            likelihood <- sum(GH$w*cond.pmf(xi,GH$x,mu,tau2,paramsi))
            return(log(likelihood))
        }
        loglik <- function(mu,tau2,x,params,GH) sum(sapply(seq_along(x),function(i)loglik.unit(mu,tau2,x[i],params[i,,drop=FALSE],GH)))
        logtau2.mle.redu <- optimize(function(logtau2)loglik(mu.null,exp(logtau2),x,params,GH),interval=logtau2.bounds,maximum=TRUE)$max
        tau2.mle.redu <- exp(logtau2.mle.redu)
        mle.full <- optim(par=c(mu.null,log(tau2.mle.redu)),fn=function(par)-loglik(par[1],exp(par[2]),x,params,GH),method='BFGS')$par
        mu.mle.full <- mle.full[1]; tau2.mle.full <- exp(mle.full[2])
        LR <- 2 * (loglik(mu.mle.full,tau2.mle.full,x,params,GH) - loglik(mu.null,tau2.mle.redu,x,params,GH))
        return(list(mle.full=c(mu=mu.mle.full,tau2=tau2.mle.full),mle.redu=c(mu=mu.null,tau2=tau2.mle.redu), LR=LR, loglik.fun=loglik, x=x, params=params, cond.pmf=cond.pmf,logtau2.bounds=logtau2.bounds, n.GH=n.GH))
}

profile.LR <- function(LR.fit,grid=list(offset=1,size=20),logtau2.bounds=fit$logtau2.bounds,n.GH=LR.fit$n.GH) {
    ## browser()
    GH <- fastGHQuad::gaussHermiteData(n.GH)
    GH$x <- sqrt(2)*GH$x; GH$w <- GH$w/sqrt(pi)
    ## names(GH)[which(names(GH)=='x')] <- 'theta.grid'
    x <- LR.fit$x
    params <- LR.fit$params
    cond.pmf <- LR.fit$cond.pmf
    mu.mle.full <- LR.fit$mle.full['mu']; tau2.mle.full <- LR.fit$mle.full['tau2'];
    loglik <- LR.fit$loglik.fun
    loglik.profile <- function(mu) {
        logtau2.mle.profile <- optimize(function(logtau2)loglik(mu,exp(logtau2),x,params,GH),interval=logtau2.bounds,maximum=TRUE)$max
        tau2.mle.profile <- exp(logtau2.mle.profile)
        return(tau2.mle.profile)
    }
    mu.profile.grid <- seq(mu.mle.full-grid$offset*sqrt(tau2.mle.full),mu.mle.full+grid$offset*sqrt(tau2.mle.full),len=grid$size)
    tau2.profile.grid <- sapply(mu.profile.grid,function(mu)loglik.profile(mu))
    ## add support for parallelization here
    ## profiled.stats <- lapply(seq_along(mu.profile.grid),function(i) {
    ##     cat('.')
    ##     LR.profile <- 2 * (loglik(mu.mle.full,tau2.mle.full,x,params,GH) - loglik(mu.profile.grid[i],tau2.profile.grid[i],x,params,GH))
    ##     ## ks.profile <- ks.ma(params,ks.fun,pmf.derivs.fun,mu=mu.profile.grid[i],tau2=tau2.profile.grid[i],n.GH=n.GH)
    ##     ks.fit.profile <- ks.ma(params,ks.fun,pmf.derivs.fun,mu=mu.profile.grid[i],tau2=tau2.profile.grid[i],n.GH=n.GH)
    ##     ks.profile <- ks.fit.profile$ks.df.ma
    ##     bartlett.C.profile <- lawley(ks.profile,2) - lawley(ks.profile,1)
    ##     return(list(LR=LR.profile,LR.bartlett=LR.profile / (1+bartlett.C.profile),bartlett.C=bartlett.C.profile,ks=ks.profile,ks.fit=ks.fit.profile))
    ## })
    ## ks.profile <- lapply(profiled.stats,`[[`,'ks')
    ## profiled.stats <- lapply(profiled.stats,`[[<-`,'ks',NULL) |> sapply(unlist)
    ## return(list(profile=profiled.stats,mu.grid=mu.profile.grid,tau2.grid=tau2.profile.grid,ks.profile=ks.profile))
    profiled.stats <- lapply(seq_along(mu.profile.grid),function(i) {
        cat('.')
        LR.profile <- 2 * (loglik(mu.mle.full,tau2.mle.full,x,params,GH) - loglik(mu.profile.grid[i],tau2.profile.grid[i],x,params,GH))
        ## ks.profile <- ks.ma(params,ks.fun,pmf.derivs.fun,mu=mu.profile.grid[i],tau2=tau2.profile.grid[i],n.GH=n.GH)
        ks.fit.profile <- ks.ma(params,ks.fun,pmf.derivs.fun,mu=mu.profile.grid[i],tau2=tau2.profile.grid[i],n.GH=n.GH)
        ks.df.ma.profile <- ks.fit.profile$ks.df.ma
        ## browser()
        bartlett.C.profile <- lawley(ks.df.ma.profile,2) - lawley(ks.df.ma.profile,1)
        return(list(LR=LR.profile,LR.bartlett=LR.profile / (1+bartlett.C.profile),bartlett.C=bartlett.C.profile,ks.fit=ks.fit.profile))
    })
    ks.fit.profile <- lapply(profiled.stats,function(lst)lst$ks.fit.profile)
    ks.df.ma.profile <- lapply(ks.fit.profile,function(lst)lst$ks.df.ma)
    profiled.stats <- lapply(profiled.stats,`[[<-`,'ks.fit',NULL) |> sapply(unlist)
    return(list(profile=profiled.stats,mu.grid=mu.profile.grid,tau2.grid=tau2.profile.grid,ks.fit.profile=ks.fit.profile,ks.df.ma=ks.df.ma.profile))
}
## confint.LR <- function(fit,logtau2.bounds=c(-5,5),grid=list(offset=1,size=20),alpha=.05,n.GH=fit$n.GH) {
confint.LR <- function(profile.fit,alpha=.05,plot=TRUE) {
    ## browser()
    ## profiled.stats <- profile.fit$profile
    ## LR.profile <- sapply(profiled.stats,function(x)x$LR)
    ## LR.bartlett.profile <- sapply(profiled.stats,function(x)x$LR.bartlett)
    profile <- profile.fit$profile
    mu.grid <- profile.fit$mu.grid
    tau2.grid <- profile.fit$tau2.grid
    q <- qchisq(1-alpha,df=1)
    if(plot) {
        plot(mu.grid,profile['LR',],ylim=range(profile[c('LR','LR.bartlett'),]),type='o')
        lines(mu.grid,profile['LR.bartlett',],lty=2,type='o')
        abline(h=q)
    }
    CIs <- apply(profile[c('LR','LR.bartlett'),] - q, 1, function(vals) {        
        left.idx <- which(diff(sign(vals))!=0)
        ## print(vals)
        ## print(left.idx)
        stopifnot(length(left.idx)==2)
        stopifnot(left.idx[2]<length(vals))
        stopifnot(vals[left.idx[2]]<0)
        right.idx <- left.idx+1 
        left.wt <- abs(vals[right.idx]/(vals[right.idx]-vals[left.idx]))
        right.wt <- abs(vals[left.idx]/(vals[right.idx]-vals[left.idx]))
        left.wt*mu.grid[left.idx] + right.wt*mu.grid[right.idx]
    } )
    return(CIs)
}


