mma.to.R <- function(x,
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
                     )) {
  stopifnot(length(x) == 1)

  chars <- strsplit(x, "", fixed = TRUE)[[1]]
  n <- length(chars)

  is_ident_start <- function(z) grepl("[A-Za-z.$`]", z)
  is_ident_char  <- function(z) grepl("[A-Za-z0-9.$`]", z)
  is_space       <- function(z) grepl("\\s", z)

  parse_until <- function(i, close = NULL) {
    out <- character()

    while (i <= n) {
      ch <- chars[i]

      if (!is.null(close) && ch == close) {
        return(list(text = paste0(out, collapse = ""), i = i + 1))
      }

      if (is_ident_start(ch)) {
        start <- i
        i <- i + 1

        while (i <= n && is_ident_char(chars[i])) {
          i <- i + 1
        }

        head <- paste0(chars[start:(i - 1)], collapse = "")

        # Allow whitespace between symbol and bracket, though InputForm usually omits it.
        ws_start <- i
        while (i <= n && is_space(chars[i])) {
          i <- i + 1
        }
        ws <- if (i > ws_start) paste0(chars[ws_start:(i - 1)], collapse = "") else ""

        if (i <= n && chars[i] == "[") {
          if (head %in% funs) {
            new_head <- if (head %in% names(rename)) unname(rename[[head]]) else head
            inner <- parse_until(i + 1, close = "]")
            out <- c(out, new_head, "(", inner$text, ")")
            i <- inner$i
          } else {
            # Preserve indexing/application brackets, but still convert known functions inside.
            inner <- parse_until(i + 1, close = "]")
            out <- c(out, head, ws, "[", inner$text, "]")
            i <- inner$i
          }
        } else {
          out <- c(out, head, ws)
        }

      } else if (ch == "[") {
        # Unheaded bracket; preserve, but recurse inside.
        inner <- parse_until(i + 1, close = "]")
        out <- c(out, "[", inner$text, "]")
        i <- inner$i

      } else {
        out <- c(out, ch)
        i <- i + 1
      }
    }

    if (!is.null(close)) {
      stop("Unmatched bracket: expected closing ", close)
    }

    list(text = paste0(out, collapse = ""), i = i)
  }

  parse_until(1)$text
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




lawley <- function(p,m,beta,l.derivs) {
    n <- length(m)
    k.vals <- sapply(l.derivs, function(l.deriv) {
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
    lawley.4 <- apply(do.call(expand.grid,rep(list(1:p),4)), 1, function(idx.str) {
        r <- idx.str[1]; s <- idx.str[2]; t <- idx.str[3]; u <- idx.str[4]
        inv.inf[r,s] * inv.inf[t,u] * (1/4*k(c(r,s,t,u)) - k(c(r,s,t),u) + k(c(r,t),c(s,u)))
    })
    lawley.6 <- apply(do.call(expand.grid,rep(list(1:p),6)), 1, function(idx.str) {
        r <- idx.str[1]; s <- idx.str[2]; t <- idx.str[3]; u <- idx.str[4]; v <- idx.str[5]; w <- idx.str[6]
        inv.inf[r,s] * inv.inf[t,u] * inv.inf[v,w] * (1/6*k(c(r,t,v))*k(c(s,u,w)) + 1/4*k(c(r,t,u))*k(c(s,v,w)) - k(c(r,t,v))*k(c(s,w),u) - k(c(r,t,u))*k(c(s,w),v) + k(c(r,t),v)*k(c(s,w),u) + k(c(r,t),u)*k(c(s,w),v))
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
            args = c('theta','m','k','beta'),
            binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt)
        ))
}
