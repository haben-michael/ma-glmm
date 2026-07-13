log1pexp <- function(z) {
    ifelse(z > 0, z + log1p(exp(-z)), log1p(exp(z)))
}

safe_bern_loglik <- function(s, n) {
    ## Supremum of Bernoulli log-likelihood with s successes in n trials.
    if (n == 0) return(0)
    if (s == 0 || s == n) return(0)
    ph <- s / n
    s * log(ph) + (n - s) * log1p(-ph)
}

logistic_full_mle_grouped <- function(x, N, S, start = c(0, 0), maxit = 50, tol = 1e-10) {
    ## Fast Newton for grouped logistic likelihood.
    ## This assumes no separation. For small n, separation can occur;
    ## see notes below.
    b <- start

    for (iter in seq_len(maxit)) {
        eta <- b[1] + b[2] * x
        p <- plogis(eta)

        score <- c(
            sum(S - N * p),
            sum(x * (S - N * p))
        )

        w <- N * p * (1 - p)

        info <- matrix(c(
            sum(w),
            sum(w * x),
            sum(w * x),
            sum(w * x^2)
        ), 2, 2)

        if (rcond(info) < 1e-12) break

        step <- solve(info, score)
        b.new <- b + step

        if (max(abs(step)) < tol) {
            b <- b.new
            break
        }

        b <- b.new
    }

    b
}

group_loglik <- function(beta, x, N, S) {
    eta <- beta[1] + beta[2] * x
    sum(S * eta - N * log1pexp(eta))
}

bartlett_jensen <- function(beta, x, N) {
    n <- sum(N)
    s <- sapply(1:4, function(j) sum(N * x^j) / n)

    v <- s[2] - s[1]^2

    term1 <- (
        (3 * s[2]^2 - 4 * s[1] * s[3] + s[4]) / v^2 - 1
    )

    term2 <- (
        (4 * s[2]^3 - 3 * s[1]^2 * s[2]^2 +
             4 * s[1]^3 * s[3] - 6 * s[1] * s[2] * s[3] + s[3]^2) / v^3 - 1
    )

    0.25 * (4 - exp(beta[1]) - exp(-beta[1])) * term1 +
        (5 / 12) * (exp(beta[1]) + exp(-beta[1]) - 2) * term2
}

compositions <- function(n, k) {
    ## All k-part compositions of n, returned as matrix.
    if (k == 1) return(matrix(n, nrow = 1))
    out <- list()
    rec <- function(prefix, remaining, slots) {
        if (slots == 1) {
            out[[length(out) + 1L]] <<- c(prefix, remaining)
        } else {
            for (z in 0:remaining) {
                rec(c(prefix, z), remaining - z, slots - 1L)
            }
        }
    }
    rec(integer(0), n, k)
    do.call(rbind, out)
}

all_S_given_N <- function(N) {
    grids <- lapply(N, function(Nj) 0:Nj)
    as.matrix(expand.grid(grids))
}

exact_logistic_cdf <- function(n, m, beta1 = 1, prob_x = 1/3) {
    x <- 0:m
    K <- length(x)

    pi.x <- dbinom(x, size = m, prob = prob_x)
    p <- plogis(beta1)

    Ns <- compositions(n, K)

    vals <- numeric(0)
    vals.b <- numeric(0)
    probs <- numeric(0)

    for (r in seq_len(nrow(Ns))) {
        N <- Ns[r, ]

        logp.N <- dmultinom(N, prob = pi.x, log = TRUE)

        Ss <- all_S_given_N(N)

        for (a in seq_len(nrow(Ss))) {
            S <- Ss[a, ]

            logp.S <- sum(dbinom(S, size = N, prob = p, log = TRUE))
            pr <- exp(logp.N + logp.S)

            if (pr == 0) next

            T <- sum(S)

            ## Reduced model under beta2 = 0.
            l.redu <- safe_bern_loglik(T, n)
            beta.redu <- c(qlogis((T + 0.5) / (n + 1)), 0)

            ## Better exact restricted beta:
            ## if T is 0 or n, the MLE is infinite and l.redu is 0.
            if (T > 0 && T < n) {
                beta.redu[1] <- qlogis(T / n)
            }
v.x <- sum(N * x^2) / n - (sum(N * x) / n)^2

if (v.x < 1e-12) {
    W <- 0
    Wb <- 0
} else {
    beta.full <- logistic_full_mle_grouped(x, N, S, start = beta.redu)
    l.full <- group_loglik(beta.full, x, N, S)
    W <- 2 * (l.full - l.redu)

    bc <- bartlett_jensen(beta.redu, x, N)
    Wb <- W / (1 + bc / n)
}
            ## Full grouped logistic MLE.
            # ## Start at restricted value.
            # beta.full <- logistic_full_mle_grouped(x, N, S, start = beta.redu)

            # l.full <- group_loglik(beta.full, x, N, S)

            # W <- 2 * (l.full - l.redu)

            # ## Jensen evaluates the factor at the restricted nuisance estimate.
            # ## You can switch to beta.full to match your simulation exactly.
            # bc <- bartlett_jensen(beta.redu, x, N)

            # Wb <- W / (1 + bc / n)

            vals <- c(vals, W)
            vals.b <- c(vals.b, Wb)
            probs <- c(probs, pr)
        }
    }

    probs <- probs / sum(probs)

    list(W = vals, Wb = vals.b, prob = probs)
}

exact_ks <- function(w, prob) {
    keep <- is.finite(w) & is.finite(prob) & prob > 0
    w <- w[keep]
    prob <- prob[keep]
    prob <- prob / sum(prob)

    u <- pchisq(w, df = 1, lower.tail = FALSE)

    o <- order(u)
    u <- u[o]
    prob <- prob[o]

    key <- signif(u, 14)
    split.prob <- rowsum(prob, key, reorder = FALSE)

    u.unique <- as.numeric(rownames(split.prob))
    p.unique <- as.numeric(split.prob[, 1])

    o2 <- order(u.unique)
    u.unique <- u.unique[o2]
    p.unique <- p.unique[o2]

    F.right <- cumsum(p.unique)
    F.left <- F.right - p.unique

    max(abs(F.right - u.unique), abs(F.left - u.unique))
}
res <- exact_logistic_cdf(n = 6, m = 3)
exact_ks(res$W, res$prob)
exact_ks(res$Wb, res$prob)
max(rowsum(res$prob, signif(pchisq(res$Wb, 1, lower.tail = FALSE), 14)))

by.m <- sapply(1:5,function(m) {
    cat('.')
ns <- 5:15
stats <- parallel::mclapply(ns,mc.cores=6,FUN=function(n) {
res <- exact_logistic_cdf(n = n, m = m)
c(LR=exact_ks(res$W, res$prob),LR.bartlett=exact_ks(res$Wb, res$prob))
})
stats <- simplify2array(stats)
apply(stats,1,function(ks)coef(lm(log(ks)~log(ns)))[2])
})

# > by.m
#                   [,1]      [,2]      [,3]      [,4]      [,5]
# LR          -0.8744659 -1.353304 -1.538647 -1.634168 -1.660474
# LR.bartlett -1.4806439 -1.835398 -1.575617 -1.895621 -1.946495