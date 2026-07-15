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



## WRONG SIGN ON INFORMATION MATRIX
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
