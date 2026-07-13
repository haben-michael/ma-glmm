    # rng = np.random.default_rng(seed=2)
    # ks = np.empty(len(ns))
    # ks_bartlett = np.empty(len(ns))
    # for i,n in enumerate(ns):
    #     print(n)
    #     reps = n**1
        
    # W = np.empty(reps)
    # W_bartlett = np.empty(reps)
    # for j in range(reps):
    #     theta = rng.normal(size=n)
    #     pr = expit(mu_null + theta*np.sqrt(tau2_null))
    #     x = rng.binomial(n=m,p=pr,size=n)

    #     loglik = lambda mu, tau2 : sum(log_unit_prob(xi,m,mu,tau2) for xi in x)
    #     def obj(params):
    #         mu, logtau2 = params 
    #         return -loglik(mu,np.exp(logtau2))
    #     opt = minimize(obj, x0=[mu_null,np.log(tau2_null)],method='BFGS')
    #     if not opt.success:
    #         W[j] = W_bartlett[j] = np.nan
    #     else:
    #         mu_mle, logtau2_mle = opt['x']
    #         tau2_mle = np.exp(logtau2_mle)

    #         W[j] = 2 * (loglik(mu_mle,tau2_mle) - loglik(mu_null,tau2_null))        
    #         W_bartlett[j] = W[j] / (1+bartlett_C/p/n)

#     ks[i] = kstest(W, chi2(df=2).cdf).statistic
#     ks_bartlett[i] = kstest(W_bartlett, chi2(df=2).cdf).statistic

# np.polyfit(np.log(list(ns)),np.log(ks),deg=1)[0]
# np.polyfit(np.log(list(ns)),np.log(ks_bartlett),deg=1)[0]

import numpy as np
import pandas as pd
from scipy.special import roots_hermite, logsumexp, gammaln, expit
from scipy.optimize import minimize,minimize_scalar
from scipy.stats import chi2, kstest
import matplotlib.pyplot as plt
from multiprocessing import Pool

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

# p=1
# bartlett_C = np.float64(0.40799473719433044)
p=2
bartlett_C =np.float64(1.7341662014725463)
mu_null = 0
tau2_null = 1/2
m = 20
ns = range(5,15+1) # change to K

# def sim(n,m,mu_null,tau2_null,bartlett_C):
def sim(args):
    seed, n = args
    rng = np.random.default_rng(seed)
    theta = rng.normal(size=n)
    pr = expit(mu_null + theta*np.sqrt(tau2_null))
    x = rng.binomial(n=m,p=pr,size=n)
    loglik = lambda mu, tau2 : sum(log_unit_prob(xi,m,mu,tau2) for xi in x)

    # full mle of mu and tau2
    def obj(params):
        mu, logtau2 = params 
        return -loglik(mu,np.exp(logtau2))
    opt = minimize(obj, x0=[mu_null,np.log(tau2_null)],method='BFGS')
    if not opt.success:
        out = [np.nan]*p
    else:
        mu_mle, logtau2_mle = opt.x
        tau2_mle = np.exp(logtau2_mle)
        W = 2 * (loglik(mu_mle,tau2_mle) - loglik(mu_null,tau2_null))
        W_bartlett = W / (1+bartlett_C/p/n)
        out = [W, W_bartlett]
    return out
    # def obj(params):
    #     logtau2 = params 
    #     return -loglik(mu_null,np.exp(logtau2))
    # opt = minimize_scalar(obj,bounds=[-10,10],method="bounded",options={"xatol": 1e-8})
    # if not opt.success:
    #     out = [np.nan,np.nan]
    # else:
    #     logtau2_mle = opt.x
    #     tau2_mle = np.exp(logtau2_mle)
    #     W = 2 * (loglik(mu_null,tau2_mle) - loglik(mu_null,tau2_null))        
    #     W_bartlett = W / (1+bartlett_C/p/n)
    #     out = [W, W_bartlett]
    # return out


if __name__ == "__main__":

    ks = np.empty(len(ns))
    ks_bartlett = np.empty(len(ns))
    for i,n in enumerate(ns):
        print(n)
        reps = n**4
        with Pool() as pool:
            stats = pool.map(sim, [(seed+700,n) for seed in range(reps)])
        W = np.array([pair[0] for pair in stats])
        W_bartlett = np.array([pair[1] for pair in stats])
        # print(W)
        ks[i] = kstest(W[np.isfinite(W)], chi2(df=p).cdf).statistic
        ks_bartlett[i] = kstest(W_bartlett[np.isfinite(W_bartlett)], chi2(df=p).cdf).statistic
        print( {"n":n, "reps":reps, "mean LR":np.mean(W), "KS":ks[i], "corr":(1+bartlett_C/p/n), "KS bart":ks_bartlett[i]} )

    rate = np.polyfit(np.log(list(ns)),np.log(ks),deg=1)[0]
    rate_bartlett = np.polyfit(np.log(list(ns)),np.log(ks_bartlett),deg=1)[0]
    print(rate)
    print(rate_bartlett)
