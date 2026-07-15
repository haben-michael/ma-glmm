import itertools
import time
from multiprocessing import Pool

import sympy as sy
import numpy as np

from scipy.integrate import quad
from scipy.special import expit


def _bound_to_float(bound, values):
    """
    Convert a SymPy integration bound to a Python float, handling +/- oo.
    """
    bound = bound.subs(values)

    if bound == sy.oo:
        return np.inf
    if bound == -sy.oo:
        return -np.inf

    return float(bound.evalf())


def _numeric_integral(integral_obj, values, modules, quad_opts):
    """
    Numerically evaluate a single SymPy Integral after substituting parameter values.
    Currently assumes a definite one-dimensional integral.
    """
    if len(integral_obj.limits) != 1:
        raise NotImplementedError(
            "This helper currently supports only one-dimensional integrals."
        )

    limit = integral_obj.limits[0]

    if len(limit) != 3:
        raise NotImplementedError(
            "This helper currently supports only definite integrals."
        )

    var, lower, upper = limit

    lower_num = _bound_to_float(lower, values)
    upper_num = _bound_to_float(upper, values)

    integrand = integral_obj.function.subs(values)

    # After substituting parameter values, only the integration variable should remain.
    remaining = integrand.free_symbols - {var}

    if remaining:
        raise ValueError(
            f"Integral still contains unsubstituted symbols: {remaining}"
        )

    f = sy.lambdify(var, integrand, modules=modules)

    val, err = quad(
        lambda t: f(t),
        lower_num,
        upper_num,
        **quad_opts,
    )

    return val


def eval_with_numeric_integrals(
    expr,
    values,
    modules=None,
    quad_opts=None,
    verbose=False,
):
    if modules is None:
        modules = [
            {"LogisticSigmoid": expit},
            "numpy",
        ]

    if quad_opts is None:
        quad_opts = {
            "epsabs": 1e-10,
            "epsrel": 1e-10,
            "limit": 200,
        }

    expr_num = expr.subs(values)

    integral_cache = {}

    while True:
        integrals = sorted(
            expr_num.atoms(sy.Integral),
            key=lambda I: sy.count_ops(I),
        )

        if not integrals:
            break

        # Take the smallest remaining integral, usually the innermost one.
        I = integrals[0]

        if I not in integral_cache:
            val = _numeric_integral(
                I,
                values={},
                modules=modules,
                quad_opts=quad_opts,
            )
            integral_cache[I] = val

            if verbose:
                print("Integral:")
                print(I)
                print("Value:", val)
                print()

        expr_num = expr_num.xreplace({I: sy.Float(integral_cache[I])})

    outer_f = sy.lambdify((), expr_num, modules=modules)
    out = outer_f()

    if np.ndim(out) == 0:
        return float(out)

    return out

class LogisticSigmoid(sy.Function):
    def fdiff(self, argindex=1):
        return self * (1 - self)
    
    # def _eval_evalf(self, prec):
    #     # 1. Get the argument as a SymPy float
    #     arg = self.args[0]._eval_evalf(prec)

m, x, theta, theta2 = sy.symbols('m x theta theta2')
beta = sy.IndexedBase('beta')
# n = 10

pdf_normal = 1 / sy.sqrt(2 * sy.pi) * sy.exp(-theta**2 / 2)
z = beta[2] + sy.sqrt(beta[1]) * theta
pr = LogisticSigmoid(z)

# integrand = sy.binomial(m, x) * (pr**x) * ((1 - pr)**(m - x)) * pdf_normal
integrand =  (pr**x) * ((1 - pr)**(m - x)) * pdf_normal
L_unit = sy.Integral(integrand, (theta,-sy.oo,sy.oo))
l_unit = sy.log(L_unit)

def k_term(args):
    order, mu_null, tau2_null, m_obs = args
    print(order)
    inner_order, outer_order = order
    inner_deriv = l_unit.diff(*[beta[i] for i in inner_order])
    inner_deriv = inner_deriv.replace(sy.Integral, lambda expr, var: sy.Integral(sy.factor(expr), var))
    binomial_pmf =  sy.binomial(m_obs,x) * (pr**x * (1-pr)**(m_obs-x) * pdf_normal).subs(theta,theta2)
    expected = sy.Integral(
        sy.summation(inner_deriv * binomial_pmf, (x,0,m_obs)),
                        (theta2,-sy.oo,sy.oo)
                        )
    if len(outer_order)==0:
        out = eval_with_numeric_integrals(expected,values={beta[1]:tau2_null,beta[2]:mu_null,m:m_obs})
    else:
        outer_deriv = expected.diff(*[beta[i] for i in outer_order])
        out = eval_with_numeric_integrals(outer_deriv,values={beta[1]:tau2_null,beta[2]:mu_null,m:m_obs})
    print(out)
    return(out)




if __name__ == "__main__":


    def lawley(p,mu_mle,tau2_mle,m_obs):
        idx = list(range(1,p+1))
        simple_orders = [[(tup,()) for tup in itertools.product(idx,repeat=rep)] for rep in [2,3,4]]
        mixed_orders = [itertools.product(itertools.product(idx,repeat=rep[0]), itertools.product(idx,repeat=rep[1]))  for rep in [[2,1],[2,2],[3,1]]]
        orders_unsorted = [item for sublist in (simple_orders+mixed_orders) for item in sublist]
        orders = [ (tuple(sorted(order[0])),tuple(sorted(order[1]))) for order in orders_unsorted]
        orders = list(set(orders))

        with Pool() as pool:    
            k_vals = pool.map(k_term, [(order,mu_mle,tau2_mle,m_obs) for order in orders])
        k_vals = dict(zip(orders,k_vals))

        def k(inner_order, outer_order=()):
            inner_order = tuple(sorted(inner_order)); outer_order = tuple(sorted(outer_order))
            return( k_vals[(inner_order,outer_order)])
        inf_mat = np.matrix( [k(inner_order) for inner_order in itertools.product(idx,repeat=2)] ).reshape(p,p) # rename--this is the negation
        inv_inf = np.linalg.inv(inf_mat) # ditto

        lawley4 = {}
        for tup in itertools.product(idx,repeat=4):
            r, s, t, u = tup
            lawley4[tup] = inv_inf[r-1,s-1] * inv_inf[t-1,u-1] * (1/4*k((r,s,t,u))
                                                                - k((r,s,t),(u,))
                                                                + k((r,t),(s,u)))

        lawley6 = {}
        for tup in itertools.product(idx,repeat=6):
            r,s,t,u,v,w = tup
            lawley6[tup] = inv_inf[r-1,s-1] * inv_inf[t-1,u-1] * inv_inf[v-1,w-1] * (1/6*k((r,t,v))*k((s,u,w)) 
                                                                    + 1/4*k((r,t,u))*k((s,v,w))
                                                                    - k((r,t,v))*k((s,w),(u,))
                                                                    - k((r,t,u))*k((s,w),(v,)) 
                                                                    + k((r,t),(v,))*k((s,w),(u,))
                                                                    + k((r,t),(u,))*k((s,w),(v,)))
            
        bias = sum(lawley4.values()) - sum(lawley6.values())
        print(k_vals)
        return bias
    # print(bias)

    mu_mle = 0
    tau2_mle = 1/2
    m_obs = 20
    p = 1
    tic = time.time()
    bias = lawley(p,mu_mle,tau2_mle,m_obs)
    toc = time.time()
    print(toc - tic)
    print(bias)

