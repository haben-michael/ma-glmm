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
