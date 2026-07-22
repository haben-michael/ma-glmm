# mma.to.R <- function(x,
#                      funs = c(
#                        "Sqrt", "Log", "Exp",
#                        "Sin", "Cos", "Tan",
#                        "ArcSin", "ArcCos", "ArcTan",
#                        "Power",
#                        "LogisticSigmoid",
#                        "Binomial"
#                      ),
#                      rename = c(
#                        Sqrt = "sqrt",
#                        Log = "log",
#                        Exp = "exp",
#                        Sin = "sin",
#                        Cos = "cos",
#                        Tan = "tan",
#                        ArcSin = "asin",
#                        ArcCos = "acos",
#                        ArcTan = "atan",
#                        LogisticSigmoid = "plogis",
#                        Binomial = "choose"
#                      )) {
#   stopifnot(length(x) == 1)

#   chars <- strsplit(x, "", fixed = TRUE)[[1]]
#   n <- length(chars)

#   is_ident_start <- function(z) grepl("[A-Za-z.$`]", z)
#   is_ident_char  <- function(z) grepl("[A-Za-z0-9.$`]", z)
#   is_space       <- function(z) grepl("\\s", z)

#   parse_until <- function(i, close = NULL) {
#     out <- character()

#     while (i <= n) {
#       ch <- chars[i]

#       if (!is.null(close) && ch == close) {
#         return(list(text = paste0(out, collapse = ""), i = i + 1))
#       }

#       if (is_ident_start(ch)) {
#         start <- i
#         i <- i + 1

#         while (i <= n && is_ident_char(chars[i])) {
#           i <- i + 1
#         }

#         head <- paste0(chars[start:(i - 1)], collapse = "")

#         # Allow whitespace between symbol and bracket, though InputForm usually omits it.
#         ws_start <- i
#         while (i <= n && is_space(chars[i])) {
#           i <- i + 1
#         }
#         ws <- if (i > ws_start) paste0(chars[ws_start:(i - 1)], collapse = "") else ""

#         if (i <= n && chars[i] == "[") {
#           if (head %in% funs) {
#             new_head <- if (head %in% names(rename)) unname(rename[[head]]) else head
#             inner <- parse_until(i + 1, close = "]")
#             out <- c(out, new_head, "(", inner$text, ")")
#             i <- inner$i
#           } else {
#             # Preserve indexing/application brackets, but still convert known functions inside.
#             inner <- parse_until(i + 1, close = "]")
#             out <- c(out, head, ws, "[", inner$text, "]")
#             i <- inner$i
#           }
#         } else {
#           out <- c(out, head, ws)
#         }

#       } else if (ch == "[") {
#         # Unheaded bracket; preserve, but recurse inside.
#         inner <- parse_until(i + 1, close = "]")
#         out <- c(out, "[", inner$text, "]")
#         i <- inner$i

#       } else {
#         out <- c(out, ch)
#         i <- i + 1
#       }
#     }

#     if (!is.null(close)) {
#       stop("Unmatched bracket: expected closing ", close)
#     }

#     list(text = paste0(out, collapse = ""), i = i)
#   }

#   parse_until(1)$text
# }


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



## ## WRONG SIGN ON INFORMATION MATRIX
## lawley <- function(p,m,beta,l.derivs) {
##     n <- length(m)
##     k.vals <- sapply(l.derivs, function(l.deriv) {
##         ## browser()
##         sapply(1:length(m), function(j) 
##             sapply(1:(m[j]-1), function(k) {
##                 ## browser()
##                 integrate(function(theta)l.deriv(theta,m[j],k,beta)*dnorm(theta),-Inf,Inf)$val
##             }) |> sum()
##             ) |> sum()
##     })
##     idx.strings <- gsub('.txt','',names(k.vals),fixed=TRUE) |> sub(pattern='.*k-',replace='') |> strsplit(split='-')
##     idx.strings <- lapply(idx.strings,function(x)if(length(x)==1)c(x,0)else x)
##     k.df <- data.frame(inner=sapply(idx.strings,`[`,1),outer=sapply(idx.strings,`[`,2),k.val=k.vals)
##     k <- function(inner,outer=0)k.df$k.val[k.df$inner==paste(inner,collapse='') & k.df$outer==paste(outer,collapse='')]
##     inf <- -matrix(apply(expand.grid(1:p,1:p),1,k),nrow=p)
##     inv.inf <- solve(inf)
##     lawley.4 <- apply(do.call(expand.grid,rep(list(1:p),4)), 1, function(idx.str) {
##         r <- idx.str[1]; s <- idx.str[2]; t <- idx.str[3]; u <- idx.str[4]
##         inv.inf[r,s] * inv.inf[t,u] * (1/4*k(c(r,s,t,u)) - k(c(r,s,t),u) + k(c(r,t),c(s,u)))
##     })
##     lawley.6 <- apply(do.call(expand.grid,rep(list(1:p),6)), 1, function(idx.str) {
##         r <- idx.str[1]; s <- idx.str[2]; t <- idx.str[3]; u <- idx.str[4]; v <- idx.str[5]; w <- idx.str[6]
##         inv.inf[r,s] * inv.inf[t,u] * inv.inf[v,w] * (1/6*k(c(r,t,v))*k(c(s,u,w)) + 1/4*k(c(r,t,u))*k(c(s,v,w)) - k(c(r,t,v))*k(c(s,w),u) - k(c(r,t,u))*k(c(s,w),v) + k(c(r,t),v)*k(c(s,w),u) + k(c(r,t),u)*k(c(s,w),v))
##     })
##     ## sum(l4)
##     ## sum(l6)
##     return( sum(lawley.4) - sum(lawley.6) )
## }
## lawley <- function(k.vals,k.p) {
##     k.df <- data.frame(inner=sapply(orders.k,`[`,1),outer=sapply(orders.k,`[`,2),k.val=k.vals)
##     ## print(k.df)
##     k <- function(inner,outer=0)k.df$k.val[k.df$inner==paste(sort(inner),collapse='') & k.df$outer==paste(sort(outer),collapse='')]
##     inf <- -matrix(apply(expand.grid(1:p,1:p),1,k),nrow=p)
##     inv.inf <- solve(inf)
##     lawley.4 <- apply(do.call(expand.grid,rep(list(1:p),4)), 1, function(idx.str) {
##         r <- idx.str[1]; s <- idx.str[2]; t <- idx.str[3]; u <- idx.str[4]
##         inv.inf[r,s] * inv.inf[t,u] * (1/4*k(c(r,s,t,u)) - k(c(r,s,t),u) + k(c(r,t),c(s,u)))
##     })
##     lawley.6 <- apply(do.call(expand.grid,rep(list(1:p),6)), 1, function(idx.str) {
##         r <- idx.str[1]; s <- idx.str[2]; t <- idx.str[3]; u <- idx.str[4]; v <- idx.str[5]; w <- idx.str[6]
##         -inv.inf[r,s] * inv.inf[t,u] * inv.inf[v,w] * (1/6*k(c(r,t,v))*k(c(s,u,w)) + 1/4*k(c(r,t,u))*k(c(s,v,w)) - k(c(r,t,v))*k(c(s,w),u) - k(c(r,t,u))*k(c(s,w),v) + k(c(r,t),v)*k(c(s,w),u) + k(c(r,t),u)*k(c(s,w),v))
##     })
##     bias <- out <- sum(lawley.4) - sum(lawley.6)
##     return(bias)
## }



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


LR.gaussian <- function(x,params,mu.null,loglik.unit,n.GH=30) {
        GH <- fastGHQuad::gaussHermiteData(30)
        GH$x <- sqrt(2)*GH$x; GH$w <- GH$w/sqrt(pi)
        loglik <- function(mu,tau2,x,params,GH=NULL) sum(sapply(seq_along(x),function(i)loglik.unit(mu,tau2,x[i],params[i,,drop=FALSE],GH=GH)))
        logtau2.mle.redu <- optimize(function(logtau2)loglik(mu.null,exp(logtau2),x,params,GH),interval=c(-5,5),maximum=TRUE)$max
        tau2.mle.redu <- exp(logtau2.mle.redu)
        mle.full <- optim(par=c(mu.null,log(tau2.mle.redu)),fn=function(par)-loglik(par[1],exp(par[2]),x,params,GH),method='BFGS')$par
        mu.mle.full <- mle.full[1]; tau2.mle.full <- exp(mle.full[2])
        LR <- 2 * (loglik(mu.mle.full,tau2.mle.full,x,params,GH) - loglik(mu.null,tau2.mle.redu,x,params,GH))
        return(list(mle.full=c(mu=mu.mle.full,tau2=tau2.mle.full),mle.redu=c(mu=mu.null,tau2=tau2.mle.redu), LR=LR))
}

loglik.unit.hyperg <- function(mu,tau2,xi,paramsi,GH=NULL) {
    m.Ti <- paramsi['m.T']
    m.Ci <- paramsi['m.C']
    si <- paramsi['s']
    pmf <- function(theta) { # numerically unsafe
        numer <- choose(m.Ti,xi)*choose(m.Ci,si-xi)*exp((mu+sqrt(tau2)*theta)*xi)
        supp <- max(0,si-m.Ci):min(si,m.Ti)
        denom <- rowSums(  outer(theta,supp,function(thetaj,xj)choose(m.Ti,xj)*choose(m.Ci,si-xj)*exp((mu+sqrt(tau2)*thetaj)*xj) ) )
        ## denom <- rowSums(  outer(theta,0:m.Ti,function(thetaj,xj)choose(m.Ti,xj)*choose(m.Ci,si-xj)*exp((mu+sqrt(tau2)*thetaj)*xj) ) )
        numer / denom
    }
    ## pmf <- function(theta) {
    ##     ## logpmf.raw <- function(xi,theta)lgamma(m.Ti+1) - lgamma(xi+1) -lgamma(m.Ti-xi+1) + lgamma(m.Ci+1) - lgamma(si-xi+1) - lgamma(m.Ci-si+xi+1) + (mu+sqrt(tau2)*theta)*xi
    ##     ## lognumer <- logpmf.raw(xi,theta)
    ##     ## logpmf.mat <- outer(0:m.Ti,theta,logpmf.raw) # speed this up
    ##     logpmf.part <- function(xi)lgamma(m.Ti+1) - lgamma(xi+1) -lgamma(m.Ti-xi+1) + lgamma(m.Ci+1) - lgamma(si-xi+1) - lgamma(m.Ci-si+xi+1)
    ##     lognumer <- logpmf.part(xi)+ (mu+sqrt(tau2)*theta)*xi
    ##     ## logpmf.mat <- logpmf.part(0:m.Ti) + outer(0:m.Ti,theta)*sqrt(tau2) + (0:m.Ti)*mu
    ##     supp <- max(0,si-m.Ci):min(si,m.Ti)
    ##     logpmf.mat <- logpmf.part(supp) + outer(supp,theta)*sqrt(tau2) + supp*mu
    ##     colmaxes <- apply(logpmf.mat,2,max)        
    ##     logdenom <- colmaxes + log(rowSums(exp(t(logpmf.mat) - colmaxes)))
    ##     exp(lognumer - logdenom)
    ## }
    integral.approx <- if(is.null(GH)) integrate(function(theta)pmf(theta)*dnorm(theta),-Inf,Inf)$val
    else
        sum(GH$w*pmf(GH$x))
    return(log(integral.approx))
}
