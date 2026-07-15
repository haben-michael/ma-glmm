import numpy as np
from scipy.special import roots_hermite, logsumexp, gammaln, expit
from scipy.optimize import minimize,minimize_scalar
from scipy.stats import chi2, kstest
import matplotlib.pyplot as plt
from multiprocessing import Pool
import pandas as pd
import math, time, pickle, random

GH_NODES, GH_WEIGHTS = roots_hermite(75)
BARTLETT_C = {'redu':np.float64(0.40799473719433044), 'full':np.float64(1.7341662014725463)}
MU_NULL = 0
TAU2_NULL = 1/2
M = 20
N_EXP = 4

def l_unit(x, m, mu, tau2, gh_nodes=GH_NODES, gh_weights=GH_WEIGHTS):
    pr = lambda theta : expit(mu + np.sqrt(tau2)*theta)
    cond_pmf = lambda theta : math.comb(m,x) * pr(theta)**x * (1-pr(theta))**(m-x)
    out = np.log(sum(cond_pmf(np.sqrt(2)*gh_nodes)/np.sqrt(np.pi)*gh_weights))
    return(out)



# def sim(n,m,mu_null,tau2_null,bartlett_C):
def sim(n):
# def sim(args):
    # seed, n = args
    # rng = np.random.default_rng(seed)
    rng = np.random.default_rng()
    theta = rng.normal(size=n)
    pr = expit(MU_NULL + theta*np.sqrt(TAU2_NULL))
    x = rng.binomial(n=M,p=pr,size=n)
    loglik = lambda mu, tau2 : sum(l_unit(xi,M,mu,tau2) for xi in x)

    # full mle of mu and tau2
    def obj(params):
        mu, logtau2 = params 
        return -loglik(mu,np.exp(logtau2))
    opt = minimize(obj, x0=[MU_NULL,np.log(TAU2_NULL)],method='BFGS')
    if not opt.success:
        mu_mle_full, tau2_mle_full = [np.nan,np.nan]
    else:
        mu_mle_full, logtau2_mle_full = opt.x
        tau2_mle_full = np.exp(logtau2_mle_full)
    def obj(params):
        logtau2 = params 
        return -loglik(MU_NULL,np.exp(logtau2))
    opt = minimize_scalar(obj,bounds=[-10,10],method="bounded",options={"xatol": 1e-8})
    if not opt.success:
        tau2_mle_redu = np.nan
    else:
        logtau2_mle_redu = opt.x
        tau2_mle_redu = np.exp(logtau2_mle_redu)
    W = 2 * (loglik(mu_mle_full,tau2_mle_full) - loglik(MU_NULL,tau2_mle_redu))
    W_bartlett = W / (1+(BARTLETT_C['full']-BARTLETT_C['redu'])/n)   
    return [W, W_bartlett]     


if __name__ == "__main__":
    tic = time.time()

    #ns = range(5,15+1) # change to K
    ns = range(15,25+1)    
    ks = np.empty(len(ns))
    ks_bartlett = np.empty(len(ns))
    for i,n in enumerate(ns):
        print(n)
        reps = n**N_EXP
        W = np.empty(reps)
        W_bartlett = np.empty(reps)
        for j in range(reps):
            stats = sim(n)
            W[j] = stats[0]
            W_bartlett[j] = stats[1]
        # with Pool() as pool:    
        #     stats = pool.map(sim, [n]*reps)
        # W = np.array([pair[0] for pair in stats])
        # W_bartlett = np.array([pair[1] for pair in stats])
        # print(W)
        ks[i] = kstest(W[np.isfinite(W)], chi2(df=1).cdf).statistic
        ks_bartlett[i] = kstest(W_bartlett[np.isfinite(W_bartlett)], chi2(df=1).cdf).statistic
        # print( {"n":n, "reps":reps, "mean LR":np.mean(W[np.isfinite(W)]), "mean LR Bartlett":np.mean(W_bartlett[np.isfinite(W_bartlett)]), "KS":ks[i], "corr":(W/W_bartlett-1)[0], "KS bart":ks_bartlett[i]} )

    with open('save' + str(random.randint(1,int(1e30))) + '.pkl', 'wb') as f:
        pickle.dump({'n_exp':N_EXP,'ns':list(ns),'ks':list(ks),'ks_bartlett':list(ks_bartlett)}, f)
    # rate = np.polyfit(np.log(list(ns)),np.log(ks),deg=1)[0]
    # rate_bartlett = np.polyfit(np.log(list(ns)),np.log(ks_bartlett),deg=1)[0]
    # print(rate)
    # print(rate_bartlett)

    toc = time.time()
    print(toc - tic)
