


import numpy as np
import pandas as pd
from scipy.special import roots_hermite, logsumexp, gammaln, expit
from scipy.optimize import minimize_scalar
from scipy.stats import chi2, kstest
import matplotlib.pyplot as plt


# ---------- Quadrature setup ----------

def gh_nodes(q=80):
    """
    For Z ~ N(0,1),
        E[f(Z)] ≈ sum_j w_std[j] * f(z[j])
    using Gauss-Hermite quadrature.
    """
    nodes, weights = roots_hermite(q)
    z = np.sqrt(2.0) * nodes
    w_std = weights / np.sqrt(np.pi)
    log_w_std = np.log(w_std)
    return z, log_w_std


Z_GH, LOGW_GH = gh_nodes(q=80)


# ---------- Logistic-normal binomial marginal likelihood ----------

def log_binom_coef(m, x):
    return gammaln(m + 1) - gammaln(x + 1) - gammaln(m - x + 1)


def log_unit_prob(x, m, mu, tau2, z_gh=Z_GH, logw_gh=LOGW_GH):
    """
    log P_beta(X=x), where

        theta ~ N(0,1)
        X | theta ~ Binomial(m, expit(mu + sqrt(tau2)*theta))

    computed by Gauss-Hermite quadrature.

    Includes the binomial coefficient. For LR maximization it cancels,
    but including it is useful for diagnostics.
    """
    if tau2 <= 0 or not np.isfinite(tau2):
        return -np.inf

    eta = mu + np.sqrt(tau2) * z_gh

    # Stable logs:
    # log p = -log(1 + exp(-eta))
    # log(1-p) = -log(1 + exp(eta))
    logp = -np.logaddexp(0.0, -eta)
    log1mp = -np.logaddexp(0.0, eta)

    log_terms = (
        logw_gh
        + log_binom_coef(m, x)
        + x * logp
        + (m - x) * log1mp
    )

    return logsumexp(log_terms)


def loglik_reduced_tau2(log_tau2, x, m, mu0):
    """
    Reduced-model log likelihood as a function of log(tau2).
    """
    tau2 = np.exp(log_tau2)
    return sum(log_unit_prob(xi, m, mu0, tau2) for xi in x)


def fit_reduced_tau2(x, m, mu0, bounds=(-10.0, 5.0)):
    """
    Maximize reduced log likelihood over tau2 > 0 via log(tau2).

    Bounds are for log(tau2). Adjust if your true tau2 is outside this range.
    """
    obj = lambda log_tau2: -loglik_reduced_tau2(log_tau2, x, m, mu0)

    opt = minimize_scalar(
        obj,
        bounds=bounds,
        method="bounded",
        options={"xatol": 1e-8}
    )

    if not opt.success:
        return {
            "success": False,
            "tau2_hat": np.nan,
            "loglik_hat": np.nan,
            "message": opt.message,
        }

    tau2_hat = np.exp(opt.x)
    loglik_hat = -opt.fun

    return {
        "success": True,
        "tau2_hat": tau2_hat,
        "loglik_hat": loglik_hat,
        "message": opt.message,
    }


def lr_reduced_point_null(x, m, mu0, tau2_0):
    """
    W = 2{ ell(tau2_hat; mu0) - ell(tau2_0; mu0) }.
    """
    fit = fit_reduced_tau2(x, m, mu0)

    if not fit["success"]:
        return np.nan, np.nan, np.nan

    ll_hat = fit["loglik_hat"]
    ll_null = loglik_reduced_tau2(np.log(tau2_0), x, m, mu0)

    W = 2.0 * (ll_hat - ll_null)

    # Numerical guard: tiny negative values can arise from optimizer tolerance.
    if W < 0 and W > -1e-8:
        W = 0.0

    return W, fit["tau2_hat"], ll_null


# ---------- Simulation ----------

def simulate_one_lr(n, m, mu0, tau2_0, rng):
    theta = rng.normal(size=n)
    p = expit(mu0 + np.sqrt(tau2_0) * theta)
    x = rng.binomial(m, p, size=n)

    W, tau2_hat, ll_null = lr_reduced_point_null(x, m, mu0, tau2_0)
    return W, tau2_hat


def ks_to_chisq1(values):
    values = np.asarray(values)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan
    return kstest(values, chi2(df=1).cdf).statistic


def run_reduced_lr_experiment(
    n_values=range(5, 16),
    m=20,
    mu0=0.0,
    tau2_0=0.5,
    reps_fun=lambda n: n**3,
    bias_const=None,
    seed=1,
):
    """
    bias_const:
        If your Lawley calculation at sample size n has bias_n ≈ C/n,
        pass C here.

        Example:
            if bias = 0.6024822126848147 was computed at n_ref=10,
            then bias_const = 10 * bias.

        Then Bartlett-adjusted statistic is
            W_B = W / (1 + bias_const/n).
    """
    rng = np.random.default_rng(seed)
    rows = []

    for n in n_values:
        reps = int(reps_fun(n))
        W = np.empty(reps)
        tau2_hat = np.empty(reps)

        for b in range(reps):
            W[b], tau2_hat[b] = simulate_one_lr(n, m, mu0, tau2_0, rng)

        W = W[np.isfinite(W)]

        if bias_const is None:
            W_bart = np.full_like(W, np.nan)
            correction = np.nan
            ks_bart = np.nan
            mean_bart = np.nan
        else:
            bias_n = bias_const / n
            correction = 1.0 + bias_n
            W_bart = W / correction
            ks_bart = ks_to_chisq1(W_bart)
            mean_bart = np.mean(W_bart)

        rows.append({
            "n": n,
            "m": m,
            "reps": len(W),
            "mean_LR": np.mean(W),
            "var_LR": np.var(W, ddof=1),
            "ks_LR": ks_to_chisq1(W),
            "correction": correction,
            "mean_LR_bartlett": mean_bart,
            "ks_LR_bartlett": ks_bart,
            "median_tau2_hat": np.nanmedian(tau2_hat),
        })

        print(
            f"n={n:2d}, reps={len(W):5d}, "
            f"mean LR={np.mean(W):.4f}, KS={ks_to_chisq1(W):.4f}, "
            f"corr={correction if np.isfinite(correction) else np.nan:.4f}, "
            f"KS bart={ks_bart if np.isfinite(ks_bart) else np.nan:.4f}"
        )

    return pd.DataFrame(rows)



# bias_at_n_ref = np.float64(0.040799473719434154)
# n_ref = 10
# bias_const = n_ref * bias_at_n_ref
bias_const = np.float64(0.40799473719433044);

res = run_reduced_lr_experiment(
    n_values=range(5, 16),
    m=20,
    mu0=0.0,
    tau2_0=0.5,
    reps_fun=lambda n: n**4,
    bias_const=bias_const,
    seed=40,
)

print(res)


plt.plot(res["n"], res["ks_LR"], marker="o", label="LR")
plt.plot(res["n"], res["ks_LR_bartlett"], marker="o", label="Bartlett-adjusted LR")
plt.xlabel("n")
plt.ylabel("KS distance to chi-square(1)")
plt.legend()
plt.show()


def loglog_slope(y, n):
    y = np.asarray(y)
    n = np.asarray(n)
    ok = np.isfinite(y) & (y > 0)
    return np.polyfit(np.log(n[ok]), np.log(y[ok]), deg=1)[0]

print("LR slope:", loglog_slope(res["ks_LR"], res["n"]))
print("Bartlett slope:", loglog_slope(res["ks_LR_bartlett"], res["n"]))